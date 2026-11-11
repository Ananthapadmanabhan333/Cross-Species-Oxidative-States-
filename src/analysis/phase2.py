"""Phase 2 orchestration: verify Phase 1 -> retrieve evidence -> analyse -> compare -> rank -> evaluate.

Run from the project root:  python -m extremonet.analysis --config config/phase2_analysis.yaml
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
import yaml

from ..config import config_fingerprint, load_config
from ..data_collection import evidence_sources as ev
from ..data_collection.http_cache import CachedHttpClient
from ..data_collection.string_client import StringClient
from ..logging_utils import configure_logging, get_logger
from ..provenance import collect_outputs, relative_posix, software_versions, write_checksums, write_manifest
from . import alignment, communities, coverage, embedding, enrichment, ranking, similarity, topology
from .inputs import graph_at_threshold, load_species, verify_phase1_run

log = get_logger("phase2")
DR, HS = "deinococcus_radiodurans", "homo_sapiens"


def derived_seed(master: int, name: str) -> int:
    return int(hashlib.sha256(f"{master}:{name}".encode()).hexdigest(), 16) % (2 ** 31)


def _tsv(frame: pd.DataFrame, path: Path) -> None:
    """Floats are written with 12 significant digits: networkx sums over string sets whose order
    depends on per-process hash randomisation, which perturbs results at ~1e-16. Rounding keeps
    reruns byte-identical without affecting any reported precision."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, float_format="%.12g")


def _round_floats(value: Any) -> Any:
    """Round floats to 12 significant digits, as `_tsv` does, so JSON summaries are also
    byte-identical across hash seeds (e.g. Louvain modularity differs at ~1e-16 otherwise)."""
    if isinstance(value, float):
        return float(f"{value:.12g}")
    if isinstance(value, dict):
        return {k: _round_floats(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round_floats(v) for v in value]
    return value


def _json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _round_floats(json.loads(json.dumps(payload, default=_default)))
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if np.isnan(obj) else float(obj)
    if isinstance(obj, (set, tuple)):
        return sorted(obj) if isinstance(obj, set) else list(obj)
    if isinstance(obj, (pd.DataFrame, np.ndarray)):
        return None
    return str(obj)


def load_phase2_config(path: str | Path, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    for dotted, value in (overrides or {}).items():
        node = cfg
        *parents, leaf = dotted.split(".")
        for key in parents:
            node = node.setdefault(key, {})
        node[leaf] = value
    weights = cfg["ranking"]["weights"]
    if set(weights) != set(ranking.COMPONENTS):
        raise ValueError(f"ranking.weights must define exactly {ranking.COMPONENTS}")
    if any(w < 0 for w in weights.values()) or abs(sum(weights.values()) - 1.0) > 1e-6:
        raise ValueError("ranking.weights must be non-negative and sum to 1")
    phase1 = load_config(cfg["phase1"]["config"], {"mode": cfg["mode"]})
    if min(cfg["network"]["sensitivity_thresholds"]) < phase1["string"]["retrieval_required_score"] / 1000:
        raise ValueError("a sensitivity threshold is below the Phase 1 retrieval score")
    cfg["_phase1"] = phase1
    cfg["_source_path"] = str(path)
    return cfg


# ------------------------------------------------------------------------------ per-species analysis
def analyse_species(sp, threshold: float, cfg: dict, string_version: str, master_seed: int) -> dict[str, Any]:
    graph = graph_at_threshold(sp, threshold, string_version)
    metrics = topology.node_metrics(graph, cfg["topology"]["hub_percentile"], cfg["topology"]["bottleneck_percentile"])
    com = cfg["communities"]
    seeds = [derived_seed(master_seed, f"louvain:{sp.key}:{i}") for i in range(com["n_seeds"])]
    stab = {res: communities.stability(graph, res, seeds) for res in com["resolutions_for_stability"]}
    if com["resolution"] not in stab:
        stab[com["resolution"]] = communities.stability(graph, com["resolution"], seeds)
    primary = stab[com["resolution"]]
    roles = topology.module_roles(graph, primary["partition"])
    node_stab = communities.node_stability(primary)
    metrics = metrics.merge(roles, on="string_id").merge(node_stab.drop(columns="module"), on="string_id")
    return {"graph": graph, "metrics": metrics, "stability": stab, "partition": primary["partition"],
            "stats": topology.graph_statistics(graph)}


# ------------------------------------------------------------------------------------------- main
def run_phase2(cfg: dict[str, Any], session: Any = None, sleep: Any = None) -> dict[str, Any]:
    p1 = cfg["_phase1"]
    root = Path.cwd()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    public = {k: v for k, v in cfg.items() if not k.startswith("_")}
    run_id = f"{stamp}_{config_fingerprint(public)}"
    out = Path(cfg["paths"]["results_dir"]) / run_id
    n = 1
    while out.exists():
        n += 1
        out = Path(cfg["paths"]["results_dir"]) / f"{run_id}-{n}"
    run_id = out.name
    configure_logging(p1.get("logging", {}).get("level", "INFO"), out / "logs" / "phase2.jsonl")
    master = int(cfg["random_seed"])
    log.info("phase 2 started", extra={"ctx": {"run_id": run_id, "mode": cfg["mode"]}})

    # 1. verify Phase 1 inputs ------------------------------------------------------------------
    p1_run = cfg["phase1"]["run_id"]
    verification = verify_phase1_run(Path(p1["paths"]["manifests_dir"]) / f"{p1_run}.json", root)
    string_version = verification["string_version"]
    log.info("phase 1 inputs verified", extra={"ctx": verification})
    processed = Path(p1["paths"]["processed_dir"]) / p1_run
    species = {k: load_species(processed, k, p1["species"][k]) for k in (DR, HS)}

    # 2. evidence retrieval ---------------------------------------------------------------------
    kwargs = {k: v for k, v in (("session", session), ("sleep", sleep)) if v is not None}
    http = CachedHttpClient.from_config(p1, **kwargs)
    string = StringClient(http, p1["string"]["base_url"], p1["string"]["caller_identity"])
    go_terms = p1["seeds"]["go_terms"]
    evidence: dict[str, dict[str, Any]] = {}
    for key, sp in species.items():
        spec, ids = sp.spec, list(sp.nodes["string_id"])
        evidence[key] = {
            "annotations": ev.string_functional_annotation(string, ids, spec["string_taxon_id"], key),
            "uniprot": ev.uniprot_evidence(http, p1["seeds"]["uniprot_base_url"], go_terms, spec["uniprot_taxon_query"],
                                           bool(spec.get("reviewed_only")), key),
            "quickgo": ev.quickgo_annotations(http, go_terms, spec["ncbi_taxon_id"], key),
            "ppi_enrichment": [ev.string_ppi_enrichment(string, ids, spec["string_taxon_id"], int(round(t * 1000)), key)
                               for t in cfg["network"]["sensitivity_thresholds"]],
        }
    dr_to_hs = ev.string_homology_best(string, list(species[DR].nodes["string_id"]), species[DR].spec["string_taxon_id"],
                                       species[HS].spec["string_taxon_id"], DR)
    hs_to_dr = ev.string_homology_best(string, list(species[HS].nodes["string_id"]), species[HS].spec["string_taxon_id"],
                                       species[DR].spec["string_taxon_id"], HS)
    kegg_cfg = cfg.get("kegg", {})
    ko_links = ({k: ev.kegg_ko_links(http, kegg_cfg["org_codes"][k]) for k in (DR, HS)}
                if kegg_cfg.get("enabled") else {DR: {}, HS: {}})
    curated = _load_curated(cfg.get("curated_evidence"), species)
    log.info("evidence retrieved", extra={"ctx": {"requests": len(http.requests_log),
                                                  "curated_rows": len(curated)}})

    sim = cfg["similarity"]
    profiles = {k: similarity.build_profiles(sp.nodes, evidence[k]["annotations"], evidence[k]["uniprot"],
                                             evidence[k]["quickgo"], k, sim["go_categories"], sim["pathway_category"],
                                             set(sim["kegg_excluded_maps"]), ko_links[k],
                                             curated[curated["species_key"] == k])
                for k, sp in species.items()}
    signature_pruning = similarity.separate_frequent_signatures(profiles, float(sim.get("max_signature_fraction", 1.0)))

    # 3. baseline analysis at the primary threshold ---------------------------------------------
    thr = float(cfg["network"]["edge_score_threshold"])
    base = {k: analyse_species(sp, thr, cfg, string_version, master) for k, sp in species.items()}
    for k, res in base.items():
        d = out / k
        _tsv(res["metrics"], d / "node_metrics.tsv")
        _tsv(topology.degree_distribution(res["graph"]), d / "degree_distribution.tsv")
        _tsv(pd.DataFrame([s["summary"] for s in res["stability"].values()]), d / "louvain_stability.tsv")
        _json(res["stats"], d / "graph_statistics.json")
        nx.write_graphml(_graphml_safe(res["graph"], res["metrics"]), d / "network_with_metrics.graphml")

    # 3b. experimentally confirmed physical interactions (IntAct) --------------------------------
    intact_summary: dict[str, Any] = {}
    if cfg.get("intact", {}).get("enabled"):
        for k, sp in species.items():
            records = ev.intact_within_set(http, _accessions(sp.nodes), k, int(cfg["intact"].get("chunk", 60)))
            table, summ = _experimental_confirmation(base[k]["graph"], sp, records)
            intact_summary[k] = summ
            _tsv(records, out / "evidence" / f"intact_records_{k}.tsv")
            _tsv(table, out / k / "edges_with_experimental_confirmation.tsv")

    # 4. module enrichment ------------------------------------------------------------------------
    enr_cfg = cfg["enrichment"]
    enrichment_summary = {}
    module_tables = {}
    for k, res in base.items():
        table, summary = enrichment.module_enrichment(res["partition"], evidence[k]["annotations"], enr_cfg["categories"],
                                                      enr_cfg["min_module_size"], enr_cfg["min_term_size"])
        table["significant"] = table["fdr_bh"] <= enr_cfg["fdr"]
        summary["significant_tests"] = int(table["significant"].sum())
        summary["modules_with_significant_term"] = int(table.loc[table["significant"], "module"].nunique())
        enrichment_summary[k] = summary
        _tsv(table, out / k / "module_enrichment.tsv")
        module_tables[k] = _module_table(res, table, profiles[k])
        _tsv(module_tables[k], out / k / "modules.tsv")

    # 5. size-matched context for the smaller network -------------------------------------------
    # the larger network is subsampled to the size of the smaller one
    small, large = sorted((DR, HS), key=lambda k: base[k]["graph"].number_of_nodes())
    null = topology.size_matched_null(base[large]["graph"], base[small]["graph"].number_of_nodes(),
                                      cfg["topology"]["size_matched_samples"], derived_seed(master, "size_matched"))
    _tsv(null, out / "cross_species" / f"size_matched_null_{large}.tsv")
    dr_stats = base[small]["stats"]
    size_matched = {metric: topology.empirical_position(null[metric], float(dr_stats[col]))
                    for metric, col in (("edges", "edges"), ("density", "density"),
                                        ("connected_nodes", "connected_nodes"),
                                        ("largest_component_nodes", "largest_component_nodes"))}
    size_matched["reference_network"] = large
    size_matched["observed_network"] = small

    # 6. cross-species: annotation overlap, topology, modules -----------------------------------
    overlap = _annotation_overlap(evidence, profiles)
    _tsv(overlap["per_category"], out / "cross_species" / "annotation_overlap_by_category.tsv")
    _tsv(overlap["shared_terms"], out / "cross_species" / "shared_terms.tsv")
    topo_cmp = _topology_comparison(base)
    _tsv(topo_cmp, out / "cross_species" / "topology_comparison.tsv")
    cores = {k: similarity.module_core_terms(base[k]["partition"], profiles[k]) for k in (DR, HS)}
    module_corr = _module_correspondence(cores, module_tables, sim["measure"])
    _tsv(module_corr, out / "cross_species" / "module_correspondence.tsv")

    # 7. alignment -------------------------------------------------------------------------------
    tables = similarity.HomologyTables.from_frames(dr_to_hs, hs_to_dr)
    dr_nodes, hs_nodes = set(species[DR].nodes["string_id"]), set(species[HS].nodes["string_id"])
    rbh = {(a, b) for a, (b, s) in tables.dr_best.items()
           if b in hs_nodes and tables.hs_best.get(b, (None,))[0] == a and s >= sim["min_homology_bitscore"]}
    al_cfg = cfg["alignment"]
    align = alignment.run_alignment(base[DR]["graph"], base[HS]["graph"], dr_to_hs, hs_to_dr, al_cfg["alphas"],
                                    al_cfg["max_iterations"], al_cfg["tolerance"], al_cfg["null_permutations"],
                                    derived_seed(master, "alignment_null"), rbh, al_cfg["rewired_null_samples"])
    _tsv(align["pairs"], out / "cross_species" / "alignment_pairs.tsv")
    _tsv(pd.DataFrame(align["summary"]), out / "cross_species" / "alignment_summary.tsv")
    homology_table = _homology_table(dr_to_hs, hs_to_dr, species, dr_nodes, hs_nodes, rbh)
    _tsv(homology_table, out / "cross_species" / "homology_best_hits.tsv")

    # 7b. coverage variant: KO-projected D. radiodurans network ---------------------------------
    cov_cfg, coverage_summary = cfg.get("coverage_variant", {}), {"enabled": False}
    if cov_cfg.get("enabled") and ko_links[DR]:
        human_kos = set().union(*(p.ko for p in profiles[HS].values()))
        variant = coverage.ko_projected_network(
            string, list(species[DR].nodes["string_id"]), species[DR].spec["string_taxon_id"], human_kos,
            ko_links[DR], int(p1["string"]["retrieval_required_score"]), thr)
        vg = variant["graph"]
        ko_by_dr = {n: set(filter(None, str(d.get("shared_kos", "")).split(";"))) for n, d in vg.nodes(data=True)}
        prior_frame = coverage.ko_prior_frame(vg, profiles[HS], ko_by_dr)
        ko_pairs = set(zip(prior_frame["query_id"], prior_frame["hit_id"]))
        v_align = alignment.run_alignment(vg, base[HS]["graph"], prior_frame, prior_frame.iloc[0:0], cov_cfg["alphas"],
                                          cov_cfg["max_iterations"], al_cfg["tolerance"], al_cfg["null_permutations"],
                                          derived_seed(master, "coverage_alignment"), ko_pairs,
                                          cov_cfg["rewired_null_samples"])
        v_stats = topology.graph_statistics(vg)
        v_metrics = topology.node_metrics(vg)
        coverage_summary = {"enabled": True, "method": "ko_projection", **variant["summary"],
                            "graph_statistics": v_stats, "alignment": v_align["summary"],
                            "note": "orthology-defined; used for coverage sensitivity only, not for ranking"}
        d = out / "coverage_variant"
        _tsv(variant["nodes"], d / "nodes.tsv")
        _tsv(v_metrics, d / "node_metrics.tsv")
        _tsv(pd.DataFrame(v_align["summary"]), d / "alignment_summary.tsv")
        _tsv(v_align["pairs"], d / "alignment_pairs.tsv")
        _json(coverage_summary, d / "coverage_summary.json")
        nx.write_graphml(vg, d / "network.graphml")

    # 8. candidate ranking -----------------------------------------------------------------------
    rk = cfg["ranking"]

    def rank_at(results: dict[str, dict], measure: str = sim["measure"]) -> pd.DataFrame:
        cores_t = {k: similarity.module_core_terms(results[k]["partition"], profiles[k]) for k in (DR, HS)}
        pairs = ranking.build_pairs(profiles[DR], profiles[HS], results[DR]["metrics"], results[HS]["metrics"],
                                    results[DR]["partition"], results[HS]["partition"], cores_t[DR], cores_t[HS],
                                    tables, sim["min_homology_bitscore"], measure)
        return pairs, ranking.rank_pairs(pairs, rk["weights"], rk["min_evidence_coverage"],
                                         rk["go_similarity_for_shared_function"], rk["pathway_similarity_for_shared_role"])

    pairs, ranked = rank_at(base)
    sens = ranking.weight_sensitivity(pairs, ranked, rk["weights"], rk["sensitivity"]["dirichlet_samples"],
                                      rk["sensitivity"]["dirichlet_concentration"], rk["sensitivity"]["top_k"],
                                      derived_seed(master, "weights"))
    ranked = ranked.merge(sens["rank_intervals"], on=["dr_string_id", "hs_string_id"], how="left")
    _tsv(ranked, out / "candidates" / "candidate_pairs_all.tsv")
    best_per_dr = (ranked[ranked["eligible"]].sort_values("score", ascending=False)
                   .groupby("dr_string_id").head(3).sort_values(["dr_string_id", "score"], ascending=[True, False]))
    _tsv(best_per_dr, out / "candidates" / "top3_human_per_dr_protein.tsv")
    _tsv(pd.DataFrame(sens["leave_one_out"]), out / "candidates" / "weights_leave_one_out.tsv")
    measure_rows = []
    for measure in sim["measures_for_sensitivity"]:
        _, ranked_m = (None, ranked) if measure == sim["measure"] else rank_at(base, measure)
        measure_rows.append({"measure": measure,
                             f"top{rk['sensitivity']['top_k']}_overlap_with_primary":
                                 ranking.top_k_overlap(ranked, ranked_m, rk["sensitivity"]["top_k"]),
                             "candidate_analogues": int((ranked_m["category"] == "candidate_analogue").sum()),
                             **{f"rank_of_{a}_{b}": _pair_rank(ranked_m, a, b)
                                for a, b in rk.get("reference_pairs", [])}})
    _tsv(pd.DataFrame(measure_rows), out / "evaluation" / "similarity_measure_sensitivity.tsv")

    # 9. threshold sensitivity (network, modules, ranking) ---------------------------------------
    k_top = rk["sensitivity"]["top_k"]
    thr_rows, ranking_thr = [], []
    for t in cfg["network"]["sensitivity_thresholds"]:
        res_t = base if t == thr else {k: analyse_species(sp, t, cfg, string_version, master) for k, sp in species.items()}
        for k in (DR, HS):
            st = res_t[k]["stats"]
            summ = res_t[k]["stability"][cfg["communities"]["resolution"]]["summary"]
            ppi = next(r for r in evidence[k]["ppi_enrichment"] if r["required_score"] == int(round(t * 1000)))
            thr_rows.append({"species": k, "threshold": t, **{c: st[c] for c in (
                "edges", "connected_nodes", "isolated_nodes", "density", "connected_components",
                "largest_component_nodes", "average_clustering")},
                "louvain_ari_mean": summ["ari_mean"], "louvain_modules_median": summ["modules_median"],
                "modularity": summ["modularity_medoid"],
                "string_ppi_enrichment_edges": ppi["number_of_edges"],
                "string_ppi_expected_edges": ppi["expected_number_of_edges"],
                "string_ppi_enrichment_p": ppi["p_value"]})
        _, ranked_t = (None, ranked) if t == thr else rank_at(res_t)
        ranking_thr.append({"threshold": t, f"top{k_top}_overlap_with_primary": ranking.top_k_overlap(ranked, ranked_t, k_top),
                            "eligible_pairs": int(ranked_t["eligible"].sum()),
                            "candidate_analogues": int((ranked_t["category"] == "candidate_analogue").sum())})
    threshold_table = pd.DataFrame(thr_rows)
    _tsv(threshold_table, out / "evaluation" / "threshold_sensitivity_networks.tsv")
    _tsv(pd.DataFrame(ranking_thr), out / "evaluation" / "threshold_sensitivity_ranking.tsv")

    # 10. embeddings -----------------------------------------------------------------------------
    emb_cfg, emb_summary = cfg["embedding"], {}
    if emb_cfg.get("enabled"):
        for k in (DR, HS):
            g = base[k]["graph"]
            if g.number_of_edges() < emb_cfg["min_edges_for_evaluation"]:
                emb_summary[k] = {"evaluated": False, "reason": f"{g.number_of_edges()} edges < "
                                  f"{emb_cfg['min_edges_for_evaluation']} required for a meaningful held-out split"}
                continue
            lp = embedding.link_prediction(g, emb_cfg["dimensions"], emb_cfg["window"], emb_cfg["negative"],
                                           emb_cfg["test_fraction"], emb_cfg["n_splits"], derived_seed(master, f"embed:{k}"))
            _tsv(lp, out / "evaluation" / f"link_prediction_{k}.tsv")
            emb_summary[k] = {"evaluated": True, "summary": embedding.summarise_link_prediction(lp)}

    # 11. proxy evaluation against RBH orthologs -------------------------------------------------
    proxy = _proxy_evaluation(ranked)
    _json(proxy, out / "evaluation" / "proxy_evaluation.json")

    # 12. evidence report tables, figures, report, manifest --------------------------------------
    protein_evidence = _protein_evidence_table(profiles)
    _tsv(protein_evidence, out / "candidates" / "protein_evidence.tsv")
    category_counts = ranked.groupby(["category", "eligible"]).size().unstack(fill_value=0).reset_index()
    _tsv(category_counts, out / "candidates" / "category_counts.tsv")
    for k in (DR, HS):
        _tsv(evidence[k]["quickgo"], out / "evidence" / f"quickgo_annotations_{k}.tsv")
        _tsv(evidence[k]["uniprot"], out / "evidence" / f"uniprot_evidence_{k}.tsv")
        _tsv(evidence[k]["annotations"], out / "evidence" / f"string_functional_annotation_{k}.tsv")

    summary = {
        "run_id": run_id, "random_seed": master, "phase1": verification, "string_version": string_version,
        "primary_threshold": thr, "graph_statistics": {k: base[k]["stats"] for k in base},
        "louvain": {k: [s["summary"] for s in base[k]["stability"].values()] for k in base},
        "enrichment": enrichment_summary, "size_matched_context": size_matched,
        "annotation_overlap": overlap["summary"], "alignment": align["summary"],
        "coverage_variant": coverage_summary,
        "experimental_interactions": intact_summary,
        "signature_pruning": signature_pruning,
        "curated_evidence": {"rows": len(curated),
                             "proteins": int(curated["string_id"].nunique()) if len(curated) else 0,
                             "source": cfg.get("curated_evidence")},
        "rbh_pairs_within_networks": len(rbh), "weight_sensitivity": {k: v for k, v in sens.items() if k != "rank_intervals"},
        "ranking_threshold_sensitivity": ranking_thr, "similarity_measure_sensitivity": measure_rows,
        "embedding": emb_summary, "proxy_evaluation": proxy,
        "module_correspondence_top": module_corr.head(10).to_dict("records"),
        "ppi_enrichment": {k: evidence[k]["ppi_enrichment"] for k in evidence},
        "category_counts": ranked.groupby("category")["eligible"].agg(["size", "sum"]).rename(
            columns={"size": "pairs", "sum": "eligible"}).reset_index().to_dict("records"),
        "evidence_availability": _evidence_availability(profiles),
    }
    _json(summary, out / "analysis_summary.json")
    (out / "config_used.yaml").write_text(yaml.safe_dump(public, sort_keys=False), encoding="utf-8")

    from . import figures, report  # matplotlib is imported only here
    figures.make_all(out, base, null, dr_stats, threshold_table, ranked, sens, emb_summary, out / "figures")
    (out / "technical_report.md").write_text(
        report.build(summary, ranked, module_tables, base, align, homology_table, cfg), encoding="utf-8")

    for handler in get_logger().handlers:
        handler.flush()
    log_dir = out / "logs"
    outputs = collect_outputs([out], root, exclude=(log_dir,))
    manifest = {
        "run_id": run_id, "phase": 2, "mode": cfg["mode"], "config": public, "phase1_verification": verification,
        "software": software_versions(), "random_seed": master,
        "raw_inputs": http.requests_log, "outputs": outputs,
        "outputs_not_checksummed": {"path": relative_posix(log_dir, root), "reason": "log still being written"},
        "served_from_counts": dict(Counter(r["served_from"] for r in http.requests_log)),
        "data_sources": {
            "string": "https://string-db.org/help/api/ (CC BY 4.0); KEGG maps as redistributed by STRING",
            "uniprot": "https://www.uniprot.org/help/api (CC BY 4.0)",
            "quickgo": "https://www.ebi.ac.uk/QuickGO/api (GO annotations, CC BY 4.0)",
        },
    }
    manifest_path = Path(cfg["paths"]["manifests_dir"]) / f"phase2_{run_id}.json"
    write_manifest(manifest_path, manifest)
    shutil.copyfile(manifest_path, out / "manifest.json")
    write_checksums(out / "SHA256SUMS.txt", outputs)
    (Path(cfg["paths"]["results_dir"]) / "LATEST_RUN.txt").write_text(run_id + "\n", encoding="utf-8")
    log.info("phase 2 finished", extra={"ctx": {"run_id": run_id, "served_from": manifest["served_from_counts"]}})
    return {"run_id": run_id, "results_dir": out, "manifest": manifest_path, "summary": summary, "ranked": ranked}


# --------------------------------------------------------------------------------------- helpers
def _accessions(nodes: pd.DataFrame) -> list[str]:
    accs: set[str] = set()
    for one, many in zip(nodes["uniprot_accession"], nodes.get("uniprot_accessions", nodes["uniprot_accession"])):
        for value in (one, many):
            if isinstance(value, str):
                accs.update(a for a in value.split(";") if a)
    return sorted(accs)


def _experimental_confirmation(graph: nx.Graph, sp, records: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Flag each STRING edge (at the analysis threshold) that IntAct supports with a physical
    interaction type, and summarise the overlap in both directions."""
    acc_to_sid: dict[str, str] = {}
    for sid, one, many in zip(sp.nodes["string_id"], sp.nodes["uniprot_accession"],
                              sp.nodes.get("uniprot_accessions", sp.nodes["uniprot_accession"])):
        for value in (one, many):
            if isinstance(value, str):
                for acc in value.split(";"):
                    acc_to_sid.setdefault(acc, sid)
    phys = records[records["physical"]].copy()
    phys["u"] = phys["uniprot_a"].map(acc_to_sid)
    phys["v"] = phys["uniprot_b"].map(acc_to_sid)
    phys = phys.dropna(subset=["u", "v"])
    self_interactions = int((phys["u"] == phys["v"]).sum())
    phys = phys[phys["u"] != phys["v"]]
    pairs: dict[tuple[str, str], dict] = {}
    for r in phys.itertuples():
        key = tuple(sorted((r.u, r.v)))
        entry = pairs.setdefault(key, {"records": 0, "methods": set(), "pubs": set(), "miscore": None,
                                       "binary": False})
        entry["records"] += 1
        entry["binary"] = entry["binary"] or bool(r.physical_binary)
        entry["methods"].add(r.detection_method.split("(")[-1].rstrip(")"))
        entry["pubs"].update(p for p in r.publications.split(";") if p)
        if r.intact_miscore is not None and not pd.isna(r.intact_miscore):
            entry["miscore"] = max(entry["miscore"] or 0.0, float(r.intact_miscore))
    rows = []
    for u, v, d in sorted(graph.edges(data=True)):
        key = tuple(sorted((u, v)))
        hit = pairs.get(key)
        rows.append({"node_u": key[0], "node_v": key[1], "combined_score": d.get("combined_score"),
                     "in_physical_subnetwork": bool(d.get("in_physical_subnetwork")),
                     "intact_confirmed": hit is not None,
                     "intact_confirmed_binary": bool(hit and hit["binary"]),
                     "intact_records": hit["records"] if hit else 0,
                     "intact_max_miscore": hit["miscore"] if hit else None,
                     "intact_methods": ";".join(sorted(hit["methods"])) if hit else "",
                     "intact_publications": ";".join(sorted(hit["pubs"])) if hit else ""})
    table = pd.DataFrame(rows, columns=["node_u", "node_v", "combined_score", "in_physical_subnetwork",
                                        "intact_confirmed", "intact_confirmed_binary", "intact_records",
                                        "intact_max_miscore", "intact_methods", "intact_publications"])
    string_edges = {tuple(sorted(e)) for e in graph.edges}
    retrieved = {tuple(sorted(p)) for p in zip(sp.functional_records["string_id_a"], sp.functional_records["string_id_b"])}
    n_edges = len(table)
    confirmed = int(table["intact_confirmed"].sum()) if n_edges else 0
    confirmed_binary = int(table["intact_confirmed_binary"].sum()) if n_edges else 0
    phys_sub = table[table["in_physical_subnetwork"]] if n_edges else table
    binary_pairs = {k for k, e in pairs.items() if e["binary"]}
    summary = {
        "accessions_queried": len(acc_to_sid), "intact_records": len(records),
        "physical_records": int(records["physical"].sum()) if len(records) else 0,
        "physical_pairs_between_seeds": len(pairs), "binary_physical_pairs_between_seeds": len(binary_pairs),
        "self_interactions": self_interactions,
        "string_edges": n_edges, "string_edges_confirmed": confirmed,
        "fraction_string_edges_confirmed": round(confirmed / n_edges, 4) if n_edges else None,
        "string_edges_confirmed_binary": confirmed_binary,
        "fraction_string_edges_confirmed_binary": round(confirmed_binary / n_edges, 4) if n_edges else None,
        "string_physical_subnetwork_edges": len(phys_sub),
        "string_physical_subnetwork_edges_confirmed": int(phys_sub["intact_confirmed"].sum()) if len(phys_sub) else 0,
        "string_physical_subnetwork_edges_confirmed_binary":
            int(phys_sub["intact_confirmed_binary"].sum()) if len(phys_sub) else 0,
        "intact_pairs_not_in_string_graph": len(set(pairs) - string_edges),
        "intact_pairs_absent_from_string_records_ge_0_4": len(set(pairs) - retrieved),
    }
    return table, summary


CURATED_COLUMNS = ["species_key", "string_id", "uniprot", "gene", "evidence_type", "effect_on_role", "assay",
                   "finding", "pmid", "doi", "curation_basis", "curated_on"]


def _load_curated(path: str | None, species: dict) -> pd.DataFrame:
    """Literature-curated findings; every row must name a network protein of its species."""
    if not path:
        return pd.DataFrame(columns=CURATED_COLUMNS)
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(CURATED_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"{path} lacks columns {sorted(missing)}")
    bad_effect = ~frame["effect_on_role"].isin(["supports", "contradicts"])
    bad_type = ~frame["evidence_type"].isin(["perturbation", "biochemical"])
    bad_pmid = ~frame["pmid"].str.fullmatch(r"\d+")
    unknown, wrong_accession = [], []
    for r in frame.itertuples():
        if r.species_key not in species or r.string_id not in set(species[r.species_key].nodes["string_id"]):
            unknown.append(r.string_id)
            continue
        node = species[r.species_key].nodes.set_index("string_id").loc[r.string_id]
        accessions = {a for v in (node.get("uniprot_accession"), node.get("uniprot_accessions"))
                      if isinstance(v, str) for a in v.split(";")}
        if r.uniprot not in accessions:
            wrong_accession.append(f"{r.string_id}:{r.uniprot}")
    if bad_effect.any() or bad_type.any() or bad_pmid.any() or unknown or wrong_accession:
        raise ValueError(f"invalid curated evidence rows in {path}: effect={int(bad_effect.sum())} "
                         f"type={int(bad_type.sum())} pmid={int(bad_pmid.sum())} unknown_proteins={unknown} "
                         f"accession_mismatch={wrong_accession}")
    return frame[CURATED_COLUMNS]


def _graphml_safe(graph: nx.Graph, metrics: pd.DataFrame) -> nx.Graph:
    g = graph.copy()
    cols = ["degree_pct", "betweenness", "pagerank", "module", "hub_bottleneck_class", "participation_coefficient"]
    for row in metrics[["string_id", *cols]].to_dict("records"):
        for c in cols:
            v = row[c]
            if v is not None and not (isinstance(v, float) and np.isnan(v)):
                g.nodes[row["string_id"]][c] = v.item() if hasattr(v, "item") else v
    return g


def _module_table(res: dict, enr: pd.DataFrame, profiles: dict) -> pd.DataFrame:
    rows = []
    sig = enr[enr["significant"]]
    by_module: dict[int, list[str]] = {}
    for n, m in res["partition"].items():
        by_module.setdefault(m, []).append(n)
    metrics = res["metrics"].set_index("string_id")
    for m, members in sorted(by_module.items()):
        top = sig[sig["module"] == m].head(5)
        rows.append({
            "module": m, "size": len(members),
            "members": ";".join(sorted(f"{profiles[n].name}" for n in members)),
            "mean_coassignment": float(np.nanmean(metrics.loc[members, "coassignment_with_module"]))
            if len(members) > 1 else np.nan,
            "significant_terms": int((sig["module"] == m).sum()),
            "top_terms": " | ".join(f"{r.description} ({r.term}, FDR={r.fdr_bh:.2g})" for r in top.itertuples()),
        })
    return pd.DataFrame(rows)


def _annotation_overlap(evidence: dict, profiles: dict) -> dict[str, Any]:
    rows, shared_rows = [], []
    per_species_terms = {}
    for k in (DR, HS):
        ann = evidence[k]["annotations"]
        per_species_terms[k] = {cat: set(grp["term"]) for cat, grp in ann.groupby("category")}
    cats = sorted(set(per_species_terms[DR]) | set(per_species_terms[HS]))
    for cat in cats:
        a = per_species_terms[DR].get(cat, set())
        b = per_species_terms[HS].get(cat, set())
        if cat == "KEGG":
            a = {similarity.kegg_map(t) for t in a} - {None}
            b = {similarity.kegg_map(t) for t in b} - {None}
        shared = a & b
        rows.append({"category": cat, "dr_terms": len(a), "hs_terms": len(b), "shared_terms": len(shared),
                     "jaccard": len(shared) / len(a | b) if a | b else None,
                     "fraction_of_dr_terms_shared": len(shared) / len(a) if a else None,
                     "comparable_across_species": cat not in {"NetworkNeighborAL", "COMPARTMENTS", "TISSUES", "RCTM"}})
        if cat in {"Process", "Function", "KEGG", "Pfam", "InterPro"}:
            desc = {}
            for k in (DR, HS):
                ann = evidence[k]["annotations"]
                for t, d in zip(ann["term"], ann["description"]):
                    desc[similarity.kegg_map(t) if cat == "KEGG" else t] = d
            for t in sorted(shared):
                shared_rows.append({"category": cat, "term": t, "description": desc.get(t, "")})
    per_cat = pd.DataFrame(rows)
    comparable = per_cat[per_cat["comparable_across_species"]]
    return {"per_category": per_cat, "shared_terms": pd.DataFrame(shared_rows),
            "summary": {"note": "Descriptive overlap of annotation vocabularies; no statistical test is implied.",
                        "categories": comparable.to_dict("records")}}


def _topology_comparison(base: dict) -> pd.DataFrame:
    rows = []
    for k in (DR, HS):
        m = base[k]["metrics"]
        conn = m[~m["is_isolated"]]
        rows.append({"species": k, **base[k]["stats"],
                     "hubs": int(m["is_hub"].sum()), "bottlenecks": int(m["is_bottleneck"].sum()),
                     "median_betweenness_connected": float(conn["betweenness"].median()) if len(conn) else None,
                     "median_clustering_connected": float(conn["clustering"].median()) if len(conn) else None})
    return pd.DataFrame(rows)


def _pair_rank(ranked: pd.DataFrame, dr_name: str, hs_name: str) -> float | None:
    hit = ranked[(ranked["dr_name"] == dr_name) & (ranked["hs_name"] == hs_name)]
    return float(hit["rank"].iloc[0]) if len(hit) and pd.notna(hit["rank"].iloc[0]) else None


def _module_correspondence(cores: dict, module_tables: dict, measure: str) -> pd.DataFrame:
    ic = similarity.information_content([*cores[DR].values(), *cores[HS].values()])
    rows = []
    for ma, ta in cores[DR].items():
        for mb, tb in cores[HS].items():
            s = similarity.simgic(ta, tb, ic, measure)
            if s is not None:
                rows.append({"dr_module": ma, "hs_module": mb, "core_term_similarity": s,
                             "shared_core_terms": len(ta & tb), "dr_core_terms": len(ta), "hs_core_terms": len(tb)})
    frame = pd.DataFrame(rows, columns=["dr_module", "hs_module", "core_term_similarity", "shared_core_terms",
                                        "dr_core_terms", "hs_core_terms"])
    if frame.empty:
        return frame
    names = {k: module_tables[k].set_index("module")["members"] for k in (DR, HS)}
    frame["dr_members"] = frame["dr_module"].map(names[DR])
    frame["hs_module_size"] = frame["hs_module"].map(module_tables[HS].set_index("module")["size"])
    return frame.sort_values("core_term_similarity", ascending=False).reset_index(drop=True)


def _homology_table(dr_to_hs, hs_to_dr, species, dr_nodes, hs_nodes, rbh) -> pd.DataFrame:
    names = {k: dict(zip(sp.nodes["string_id"], sp.nodes["string_preferred_name"])) for k, sp in species.items()}
    frame = dr_to_hs.copy()
    frame["dr_name"] = frame["query_id"].map(names[DR])
    frame["hit_in_human_network"] = frame["hit_id"].isin(hs_nodes)
    frame["hit_name"] = frame["hit_id"].map(names[HS])
    frame["reciprocal_best_hit"] = [(a, b) in rbh for a, b in zip(frame["query_id"], frame["hit_id"])]
    return frame.sort_values("bitscore", ascending=False).reset_index(drop=True)


MIN_PROXY_POSITIVES = 10


def _proxy_evaluation(ranked: pd.DataFrame) -> dict[str, Any]:
    """No ground truth for cross-species functional analogy exists. As a limited proxy, ask whether
    individual components rank RBH ortholog pairs above other pairs. Annotation components are
    partly transferred by homology (IEA/IBA), so their AUC is inflated by construction (leakage);
    the topology-only 'role' component has no such dependency."""
    labels = ranked["category"].isin(["ortholog_rbh", "ortholog_ko"])
    out: dict[str, Any] = {"ground_truth": "none available; orthologs (reciprocal best hits or shared KEGG KO) "
                                           "used as a proxy positive set",
                           "positives": int(labels.sum()), "components": {}}
    for comp in ("role", "module", "pathway", "go", "ec"):
        vals = ranked[comp].astype(float)
        mask = vals.notna()
        pos, neg = vals[mask & labels].to_numpy(), vals[mask & ~labels].to_numpy()
        if len(pos) < 3 or len(neg) < 3:
            out["components"][comp] = {"evaluable": False, "positives_with_value": int(len(pos))}
            continue
        out["components"][comp] = {"evaluable": len(pos) >= MIN_PROXY_POSITIVES,
                                   "auc_descriptive": embedding._auc(pos, neg),
                                   "positives_with_value": int(len(pos)), "negatives_with_value": int(len(neg)),
                                   # module similarity is computed from GO/KEGG core terms, so it shares
                                   # the annotation-transfer leakage of the annotation components
                                   "leakage_risk": comp in ("go", "ec", "pathway", "module")}
    out["note"] = (f"An AUC is called evaluable only with >= {MIN_PROXY_POSITIVES} positives; below that it is "
                   "listed as descriptive and must not be read as a performance estimate.")
    return out


def _protein_evidence_table(profiles: dict) -> pd.DataFrame:
    rows = []
    for k, prof in profiles.items():
        for p in prof.values():
            rows.append({"species": k, "string_id": p.string_id, "name": p.name, "uniprot": p.accession,
                         "ec": ";".join(sorted(p.ec)), "panther": ";".join(sorted(p.panther)),
                         "pfam": ";".join(sorted(p.pfam)), "interpro": len(p.interpro), "go_terms": len(p.go_terms),
                         "kegg_ko": ";".join(sorted(p.ko)), "perturbation_score": p.perturbation_score,
                         "curated_findings": " | ".join(f"{c['effect_on_role']}: {c['finding']} (PMID:{c['pmid']})"
                                                        for c in p.curated),
                         "kegg_maps": ";".join(sorted(p.pathways)), "literature_level": p.literature_level,
                         "evidence_codes": ";".join(sorted(p.evidence_codes)),
                         "perturbation_evidence": p.perturbation,
                         "experimental_references": ";".join(sorted(p.experimental_references)),
                         "negated_go": ";".join(sorted(p.negated_go))})
    return pd.DataFrame(rows)


def _evidence_availability(profiles: dict) -> dict[str, Any]:
    out = {}
    for k, prof in profiles.items():
        vals = list(prof.values())
        n = len(vals)
        out[k] = {"proteins": n,
                  "with_ec": sum(bool(p.ec) for p in vals), "with_go": sum(bool(p.go_terms) for p in vals),
                  "with_kegg_map": sum(bool(p.pathways) for p in vals),
                  "with_panther": sum(bool(p.panther) for p in vals), "with_pfam": sum(bool(p.pfam) for p in vals),
                  "with_quickgo_annotation": sum(p.literature_level is not None for p in vals),
                  "with_experimental_go": sum(p.literature_level == 1.0 for p in vals),
                  "with_perturbation_evidence": sum(p.perturbation for p in vals),
                  "with_kegg_ko": sum(bool(p.ko) for p in vals),
                  "with_curated_literature": sum(bool(p.curated) for p in vals)}
    return out
