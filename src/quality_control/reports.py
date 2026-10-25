"""Human-readable (Markdown) and machine-readable (JSON) data-quality reports."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _table(rows: list[tuple[str, Any]]) -> str:
    out = ["| item | value |", "|---|---|"]
    out += [f"| {k} | {v} |" for k, v in rows]
    return "\n".join(out)


def _reasons(counts: dict[str, int]) -> str:
    if not counts:
        return "_none_"
    return "\n".join(f"- `{reason}`: {count}" for reason, count in counts.items())


def species_qc_markdown(qc: dict[str, Any]) -> str:
    cov, ann, stats = qc["network_coverage"], qc["annotation_coverage"], qc["summary_statistics"]
    tax = qc["taxonomy"]
    lines = [
        f"# Data-quality report: {qc['species_label']}",
        "",
        f"Run `{qc['run_id']}` · STRING {qc['string_version']} · edge threshold "
        f"{qc['edge_score_threshold']} (retrieved at STRING required_score {qc['retrieval_required_score']})",
        "",
        "## Species assignment",
        _table([
            ("NCBI taxon (configured)", tax["ncbi_taxon_id"]),
            ("STRING taxon (expected)", tax["expected_string_taxon_id"]),
            ("STRING taxa returned by identifier mapping", tax["string_taxa_returned"]),
            ("UniProt organism ids among seeds", tax["uniprot_organism_ids_in_seeds"]),
            ("all retained records match expected STRING taxon", tax["all_retained_match"]),
        ]),
        "",
        "## Seeds (UniProtKB, GO-defined)",
        f"Query: `{qc['seed_query']}`",
        "",
        _table([
            ("UniProtKB entries for organism (context)", cov["uniprot_entries_for_organism"]),
            ("seed records retrieved", cov["seed_records_retrieved"]),
            ("seed fraction of organism entries", cov["seed_fraction_of_organism_entries"]),
            ("seeds passing validation", cov["seeds_valid"]),
            ("seeds mapped to STRING", cov["seeds_mapped_to_string"]),
            ("mapping rate", cov["mapping_rate"]),
        ]),
        "",
        "Excluded seed records:",
        _reasons(qc["excluded"]["seeds"]),
        "",
        "Excluded mapping records:",
        _reasons(qc["excluded"]["mapping"]),
        "",
        "## Edges",
        _table([
            ("functional records retrieved", qc["edges"]["functional_records_retrieved"]),
            ("functional edges retained", qc["edges"]["functional_retained"]),
            ("physical-subnetwork records retrieved", qc["edges"]["physical_records_retrieved"]),
            ("physical-subnetwork edges retained", qc["edges"]["physical_retained"]),
            ("physical edges without functional counterpart", qc["edges"]["physical_without_functional"]),
        ]),
        "",
        "Excluded functional records:",
        _reasons(qc["excluded"]["functional_edges"]),
        "",
        "Excluded physical records:",
        _reasons(qc["excluded"]["physical_edges"]),
        "",
        "## Network coverage",
        _table([
            ("nodes in graph", cov["nodes_in_graph"]),
            ("nodes with >=1 retained edge", cov["nodes_with_edges"]),
            ("isolated nodes (kept, flagged)", cov["isolated_nodes"]),
            ("fraction of mapped seeds connected", cov["connected_fraction_of_mapped_seeds"]),
        ]),
        "",
        "## Annotation coverage",
        _table([(k, v if not isinstance(v, dict) else f"{v['count']} ({v['fraction']})") for k, v in ann.items()]),
        "",
        "## Graph summary statistics",
        _table(list(stats.items())),
        "",
        "## Limitations",
        *[f"- {item}" for item in qc["limitations"]],
        "",
    ]
    return "\n".join(lines)


def run_markdown(run: dict[str, Any], comparison_rows: list[dict[str, Any]]) -> str:
    keys = list(comparison_rows[0].keys()) if comparison_rows else []
    header = "| " + " | ".join(keys) + " |"
    sep = "|" + "---|" * len(keys)
    body = ["| " + " | ".join(str(r[k]) for k in keys) + " |" for r in comparison_rows]
    lines = [
        f"# Run report `{run['run_id']}`",
        "",
        f"- Mode: `{run['mode']}`",
        f"- STRING version: {run['string_version']}",
        f"- Edge threshold (combined score): {run['edge_score_threshold']}",
        f"- Config: `{run['config_path']}` (fingerprint `{run['config_fingerprint']}`)",
        f"- Requests served: {run['served_from_counts']}",
        "",
        "## Side-by-side network description",
        "Descriptive statistics only. Seed sets are defined independently per species from GO",
        "annotation, so differences reflect annotation depth as well as biology; no cross-species",
        "orthology mapping is performed in Phase 1.",
        "",
        header, sep, *body,
        "",
        "Per-species details: `<species>/qc_report.md`.",
        "",
    ]
    return "\n".join(lines)
