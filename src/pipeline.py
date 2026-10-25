"""Phase 1 pipeline: GO-defined oxidative-stress seeds -> STRING networks -> QC -> weighted graphs.

Run from the project root:  python -m extremonet --config config/phase1_oxidative_stress.yaml
"""

from __future__ import annotations

import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd
import yaml

from .annotation.annotate import add_network_only_nodes, annotation_coverage, build_node_table
from .config import config_fingerprint
from .data_collection.http_cache import CachedHttpClient
from .data_collection.string_client import StringClient
from .data_collection.uniprot_client import UniProtClient, build_seed_query
from .logging_utils import configure_logging, get_logger
from .network_construction.graph_builder import (
    annotate_edges, build_graph, degree_table, summary_statistics,
)
from .preprocessing.parsers import parse_string_mapping, parse_string_network, parse_uniprot_seeds
from .provenance import collect_outputs, relative_posix, software_versions, write_checksums, write_manifest
from .quality_control.reports import run_markdown, species_qc_markdown, write_json
from .quality_control.validators import reason_counts, validate_edges, validate_mapping, validate_seeds

log = get_logger("pipeline")

STANDING_LIMITATIONS = [
    "STRING edges are functional associations (combined evidence channels), not physical "
    "interactions. `in_physical_subnetwork` marks pairs in STRING's physical-complex subnetwork, "
    "which is itself a scored prediction; `has_experiments_channel_score` reflects STRING's experiments "
    "channel, which can include evidence transferred from other organisms. Neither is equivalent to "
    "an experimentally confirmed physical interaction (e.g. from BioGRID/IntAct).",
    "Seeds are defined by GO annotation in UniProtKB, so the network covers only proteins already "
    "annotated to the configured terms; poorly annotated proteomes yield smaller seed sets.",
    "With string.add_nodes = 0 only interactions among seeds are retrieved; associations to "
    "non-seed proteins are outside the graph by design.",
]


def _run_id(config: dict[str, Any], results_dir: Path) -> str:
    """UTC timestamp + config fingerprint; suffixed if a run with that id already exists."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base = f"{stamp}_{config_fingerprint(config)}"
    run_id, n = base, 1
    while (results_dir / run_id).exists():
        n += 1
        run_id = f"{base}-{n}"
    return run_id


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def _excluded_table(stage: str, frame: pd.DataFrame, id_cols: list[str]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=["stage", "record", "exclusion_reason", "detail"])
    detail_cols = [c for c in ("combined_score", "string_taxon_id", "uniprot_organism_id", "network_type")
                   if c in frame.columns]
    return pd.DataFrame({
        "stage": stage,
        "record": frame[id_cols].astype(str).agg(" | ".join, axis=1),
        "exclusion_reason": frame["exclusion_reason"],
        "detail": frame[detail_cols].astype(str).agg("; ".join, axis=1) if detail_cols else "",
    })


def run_species(species_key: str, species: dict[str, Any], config: dict[str, Any], run_id: str,
                string: StringClient, uniprot: UniProtClient, string_version: str,
                processed_dir: Path, results_dir: Path) -> dict[str, Any]:
    ctx = {"species": species_key}
    log.info("species started", extra={"ctx": ctx})
    go_terms = config["seeds"]["go_terms"]
    reviewed_only = bool(species.get("reviewed_only", False))
    expected_taxon = species["string_taxon_id"]
    threshold = float(config["network"]["edge_score_threshold"])
    retrieval = int(config["string"]["retrieval_required_score"])
    add_nodes = int(config["string"].get("add_nodes", 0))

    # 1. seeds ---------------------------------------------------------------------------------
    seed_raw = uniprot.fetch_seed_proteins(go_terms, species["uniprot_taxon_query"], reviewed_only, species_key)
    seeds_all = parse_uniprot_seeds(seed_raw.text)
    organism_entries, _ = uniprot.fetch_entry_count(species["uniprot_taxon_query"], reviewed_only, species_key)
    allowed = {species["ncbi_taxon_id"], species["string_taxon_id"]}
    seeds, seeds_excluded = validate_seeds(seeds_all, allowed)
    log.info("seeds validated", extra={"ctx": {**ctx, "retrieved": len(seeds_all), "valid": len(seeds)}})
    if seeds.empty:
        raise RuntimeError(f"{species_key}: no valid seed proteins; see the UniProt query and QC output")

    # 2. identifier mapping (queried with the configured NCBI taxon; STRING resolves it) -------
    mapping_rows, _ = string.map_identifiers(list(seeds["uniprot_accession"]), species["ncbi_taxon_id"], species_key)
    mapping = parse_string_mapping(mapping_rows).assign(mapping_pass="batch")
    # A batch answer holds one row per STRING protein, so accessions that share a protein can be
    # missing from it. Re-query those one at a time so collapses are not reported as gaps.
    missing = sorted(set(seeds["uniprot_accession"]) - set(mapping["uniprot_accession"]))
    max_single = int(config["string"].get("max_individual_remap", 50))
    if len(missing) > max_single:
        log.warning("too many unmapped accessions to re-query individually; reported as unmapped",
                    extra={"ctx": {**ctx, "unmapped": len(missing), "max_individual_remap": max_single}})
        missing = []
    for accession in missing:
        rows, _ = string.map_identifiers([accession], species["ncbi_taxon_id"], species_key, endpoint_suffix="_single")
        if rows:
            mapping = pd.concat([mapping, parse_string_mapping(rows).assign(mapping_pass="single")],
                                ignore_index=True)
    if missing:
        log.info("individual re-query of batch-unmapped accessions",
                 extra={"ctx": {**ctx, "queried": len(missing),
                                "resolved": int((mapping["mapping_pass"] == "single").sum())}})
    mapped, mapping_excluded = validate_mapping(seeds, mapping, expected_taxon)
    returned_taxa = sorted({int(t) for t in pd.to_numeric(mapping["string_taxon_id"], errors="coerce").dropna()})
    log.info("identifiers mapped", extra={"ctx": {**ctx, "mapped": len(mapped), "string_taxa": returned_taxa}})
    if mapped.empty:
        raise RuntimeError(f"{species_key}: no seed mapped to STRING taxon {expected_taxon}")

    nodes = build_node_table(mapped, species_key, species, go_terms)
    node_ids = set(nodes["string_id"])

    # 3. networks ------------------------------------------------------------------------------
    retained: dict[str, pd.DataFrame] = {}
    excluded: dict[str, pd.DataFrame] = {}
    retrieved_counts: dict[str, int] = {}
    for network_type in config["string"]["network_types"]:
        rows, _ = string.network(sorted(node_ids), expected_taxon, retrieval, network_type, add_nodes, species_key)
        edges = parse_string_network(rows, network_type)
        retrieved_counts[network_type] = len(edges)
        keep, drop = validate_edges(edges, node_ids, expected_taxon, threshold, allow_new_nodes=add_nodes > 0)
        retained[network_type], excluded[network_type] = keep, drop
        _write_tsv(edges, processed_dir / f"string_{network_type}_records_parsed.tsv")
        log.info("edges validated", extra={"ctx": {**ctx, "network_type": network_type,
                                                   "retrieved": len(edges), "retained": len(keep)}})

    physical = retained.get("physical")
    source_label = f"STRING v{string_version} functional association"
    edges, physical_orphans = annotate_edges(
        retained["functional"], physical, float(config["network"].get("experimental_score_threshold", 0.0)),
        source_label)
    if add_nodes > 0:
        nodes = add_network_only_nodes(nodes, edges.rename(columns={"node_u": "string_id_a", "node_v": "string_id_b"}),
                                       species_key, species)
    if len(physical_orphans):
        log.warning("physical edges without a functional counterpart (not added to graph)",
                    extra={"ctx": {**ctx, "count": len(physical_orphans)}})

    # 4. graph ---------------------------------------------------------------------------------
    graph = build_graph(nodes, edges, keep_isolated=bool(config["network"].get("keep_isolated_nodes", True)),
                        graph_attrs={
                            "species": species["label"], "ncbi_taxon_id": species["ncbi_taxon_id"],
                            "string_taxon_id": expected_taxon, "string_version": string_version,
                            "edge_score_threshold": threshold, "edge_weight": "STRING combined score (0-1)",
                            "run_id": run_id,
                        })
    stats = summary_statistics(graph)
    connected = set(edges["node_u"]) | set(edges["node_v"])
    nodes["has_retained_edge"] = nodes["string_id"].isin(connected)

    coverage = {
        "uniprot_entries_for_organism": organism_entries,
        "seed_records_retrieved": len(seeds_all),
        "seed_fraction_of_organism_entries": round(len(seeds_all) / organism_entries, 4) if organism_entries else None,
        "seeds_valid": len(seeds),
        "seeds_mapped_to_string": len(mapped),
        "mapping_rate": round(len(mapped) / len(seeds), 4),
        "nodes_in_graph": graph.number_of_nodes(),
        "nodes_with_edges": len(connected),
        "isolated_nodes": stats["isolated_nodes"],
        "connected_fraction_of_mapped_seeds": round(len(connected & node_ids) / len(node_ids), 4),
    }
    ann = annotation_coverage(nodes)

    limitations = list(STANDING_LIMITATIONS)
    if expected_taxon != species["ncbi_taxon_id"]:
        limitations.append(
            f"STRING represents {species['label']} by taxon {expected_taxon} (a single strain), not the "
            f"configured species taxon {species['ncbi_taxon_id']}; results describe that strain's proteome.")
    if stats["isolated_nodes"]:
        limitations.append(f"{stats['isolated_nodes']} of {graph.number_of_nodes()} seed proteins have no "
                           f"association at combined score >= {threshold} with another seed.")
    if len(mapping_excluded):
        limitations.append(f"{len(mapping_excluded)} valid seed records were not carried into the graph at "
                           "the identifier-mapping stage (see excluded_records.tsv).")
    if ann["uniprot_xref_disagrees_with_string_mapping"]["count"]:
        limitations.append(f"{ann['uniprot_xref_disagrees_with_string_mapping']['count']} proteins have a UniProt "
                           "STRING cross-reference that differs from the STRING API mapping.")

    qc = {
        "run_id": run_id, "species_key": species_key, "species_label": species["label"],
        "string_version": string_version, "edge_score_threshold": threshold,
        "retrieval_required_score": retrieval,
        "seed_query": build_seed_query(go_terms, species["uniprot_taxon_query"], reviewed_only),
        "taxonomy": {
            "ncbi_taxon_id": species["ncbi_taxon_id"], "expected_string_taxon_id": expected_taxon,
            "string_taxa_returned": returned_taxa,
            "uniprot_organism_ids_in_seeds": sorted(set(seeds_all["uniprot_organism_id"])),
            "all_retained_match": bool(
                all(t == expected_taxon for t in nodes["string_id"].map(lambda s: int(s.split(".")[0])))),
        },
        "excluded": {
            "seeds": reason_counts(seeds_excluded), "mapping": reason_counts(mapping_excluded),
            "functional_edges": reason_counts(excluded["functional"]),
            "physical_edges": reason_counts(excluded["physical"]) if "physical" in excluded else {},
        },
        "edges": {
            "functional_records_retrieved": retrieved_counts["functional"],
            "functional_retained": len(retained["functional"]),
            "physical_records_retrieved": retrieved_counts.get("physical"),
            "physical_retained": len(physical) if physical is not None else None,
            "physical_without_functional": len(physical_orphans),
        },
        "network_coverage": coverage, "annotation_coverage": ann, "summary_statistics": stats,
        "limitations": limitations,
    }

    # 5. outputs -------------------------------------------------------------------------------
    excluded_all = pd.concat([
        _excluded_table("seed", seeds_excluded, ["uniprot_accession"]),
        _excluded_table("mapping", mapping_excluded, ["uniprot_accession", "string_id"]),
        *[_excluded_table(f"edge_{t}", excluded[t], ["string_id_a", "string_id_b"]) for t in excluded],
    ], ignore_index=True)
    _write_tsv(seeds_all, processed_dir / "uniprot_seed_records_parsed.tsv")
    _write_tsv(mapping, processed_dir / "string_id_mapping_parsed.tsv")
    _write_tsv(nodes, processed_dir / "nodes.tsv")
    _write_tsv(edges, processed_dir / "edges_functional_retained.tsv")
    if physical is not None:
        _write_tsv(physical, processed_dir / "edges_physical_retained.tsv")
    _write_tsv(excluded_all, processed_dir / "excluded_records.tsv")

    results_dir.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(graph, results_dir / "network.graphml")
    _write_tsv(edges, results_dir / "weighted_edge_list.tsv")
    _write_tsv(nodes, results_dir / "protein_table.tsv")
    _write_tsv(degree_table(graph), results_dir / "degree_table.tsv")
    write_json(results_dir / "summary_statistics.json", stats)
    write_json(results_dir / "qc_report.json", qc)
    (results_dir / "qc_report.md").write_text(species_qc_markdown(qc), encoding="utf-8")
    log.info("species finished", extra={"ctx": {**ctx, "nodes": stats["nodes"], "edges": stats["edges"]}})
    return qc


def run_pipeline(config: dict[str, Any], session: Any = None, sleep: Any = None) -> dict[str, Any]:
    paths = {k: Path(v) for k, v in config["paths"].items()}
    run_id = _run_id(config, paths["results_dir"])
    results_root = paths["results_dir"] / run_id
    processed_root = paths["processed_dir"] / run_id
    configure_logging(config.get("logging", {}).get("level", "INFO"), results_root / "logs" / "pipeline.jsonl")
    log.info("run started", extra={"ctx": {"run_id": run_id, "mode": config["mode"]}})

    client_kwargs: dict[str, Any] = {}
    if session is not None:
        client_kwargs["session"] = session
    if sleep is not None:
        client_kwargs["sleep"] = sleep
    http = CachedHttpClient.from_config(config, **client_kwargs)
    string = StringClient(http, config["string"]["base_url"], config["string"]["caller_identity"])
    uniprot = UniProtClient(http, config["seeds"].get("uniprot_base_url", "https://rest.uniprot.org/uniprotkb"))

    version_info, _ = string.version()
    string_version = str(version_info.get("string_version"))

    results_root.mkdir(parents=True, exist_ok=True)
    public_config = {k: v for k, v in config.items() if not k.startswith("_")}
    (results_root / "config_used.yaml").write_text(yaml.safe_dump(public_config, sort_keys=False), encoding="utf-8")

    species_qc = {}
    for key, spec in config["species"].items():
        species_qc[key] = run_species(key, spec, config, run_id, string, uniprot, string_version,
                                      processed_root / key, results_root / key)

    comparison = [{
        "species": qc["species_label"], "ncbi_taxon": qc["taxonomy"]["ncbi_taxon_id"],
        "string_taxon": qc["taxonomy"]["expected_string_taxon_id"],
        "seeds_retrieved": qc["network_coverage"]["seed_records_retrieved"],
        "seeds_mapped": qc["network_coverage"]["seeds_mapped_to_string"],
        **{k: qc["summary_statistics"][k] for k in (
            "nodes", "edges", "isolated_nodes", "density", "connected_components", "largest_component_nodes",
            "degree_mean", "average_clustering_unweighted", "edges_in_physical_subnetwork",
            "edges_with_experiments_channel_score")},
    } for qc in species_qc.values()]
    _write_tsv(pd.DataFrame(comparison), results_root / "species_comparison.tsv")

    served = dict(Counter(r["served_from"] for r in http.requests_log))
    run_info = {
        "run_id": run_id, "mode": config["mode"], "string_version": string_version,
        "edge_score_threshold": config["network"]["edge_score_threshold"],
        "config_path": config.get("_source_path", "<in-memory>"),
        "config_fingerprint": config_fingerprint(config), "served_from_counts": served,
    }
    (results_root / "run_report.md").write_text(run_markdown(run_info, comparison), encoding="utf-8")
    log.info("run outputs written", extra={"ctx": {"results": str(results_root)}})

    # Manifest last, so it can checksum every output. Log handlers are flushed first.
    for handler in get_logger().handlers:
        handler.flush()
    root = Path.cwd()
    log_dir = results_root / "logs"
    outputs = collect_outputs([results_root, processed_root], root, exclude=(log_dir,))
    manifest = {
        **run_info,
        "project": config.get("project", {}),
        "string_release": version_info,
        "software": software_versions(),
        "config": public_config,
        "taxonomy_resolution": {k: qc["taxonomy"] for k, qc in species_qc.items()},
        "data_sources": {
            "string": {"base_url": config["string"]["base_url"], "documentation": "https://string-db.org/help/api/",
                       "license": "CC BY 4.0 (https://string-db.org/cgi/access)"},
            "uniprot": {"base_url": config["seeds"].get("uniprot_base_url"),
                        "documentation": "https://www.uniprot.org/help/api", "license": "CC BY 4.0"},
        },
        "raw_inputs": http.requests_log,
        "outputs": outputs,
        "outputs_not_checksummed": {
            "path": relative_posix(log_dir, root),
            "reason": "run log is still being appended to when the manifest is written",
        },
        "limitations": {k: qc["limitations"] for k, qc in species_qc.items()},
    }
    manifest_path = paths["manifests_dir"] / f"{run_id}.json"
    write_manifest(manifest_path, manifest)
    shutil.copyfile(manifest_path, results_root / "manifest.json")
    write_checksums(results_root / "SHA256SUMS.txt", outputs)
    (paths["results_dir"] / "LATEST_RUN.txt").write_text(run_id + "\n", encoding="utf-8")
    log.info("run finished", extra={"ctx": {"run_id": run_id, "manifest": str(manifest_path), "served_from": served}})
    return {"run_id": run_id, "results_dir": results_root, "processed_dir": processed_root,
            "manifest": manifest_path, "species": species_qc}
