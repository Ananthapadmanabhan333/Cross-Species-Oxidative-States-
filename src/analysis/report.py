"""Generated technical report. Every number is read from the current run's computed results."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .ranking import CATEGORY_DESCRIPTIONS, CATEGORY_ORDER

DR, HS = "deinococcus_radiodurans", "homo_sapiens"
LABEL = {DR: "*D. radiodurans*", HS: "*H. sapiens*"}


def _fmt(v: Any, nd: int = 3) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    if isinstance(v, (float, np.floating)):
        if v == 0:
            return "0"
        if abs(v) < 1e-3:
            return f"{v:.2e}"
        return f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return str(v)


def _table(frame: pd.DataFrame, cols: list[str] | None = None) -> str:
    if frame is None or frame.empty:
        return "_none_"
    frame = frame[cols] if cols else frame
    head = "| " + " | ".join(frame.columns) + " |"
    sep = "|" + "---|" * len(frame.columns)
    body = ["| " + " | ".join(_fmt(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def build(summary: dict, ranked: pd.DataFrame, module_tables: dict, base: dict, align: dict,
          homology_table: pd.DataFrame, cfg: dict) -> str:
    s = summary
    gs = s["graph_statistics"]
    stats_rows = pd.DataFrame([{"metric": m, LABEL[DR]: gs[DR].get(m), LABEL[HS]: gs[HS].get(m)} for m in (
        "nodes", "edges", "density", "connected_nodes", "isolated_nodes", "connected_components",
        "non_trivial_components", "largest_component_nodes", "degree_mean", "degree_max", "transitivity",
        "average_clustering", "lcc_average_shortest_path_hops", "lcc_diameter_hops", "degree_assortativity")])
    ppi = pd.DataFrame([{"species": LABEL[k], "threshold": r["required_score"] / 1000, "edges": r["number_of_edges"],
                         "expected_edges": r["expected_number_of_edges"], "p_value": r["p_value"]}
                        for k in (DR, HS) for r in s["ppi_enrichment"][k]])
    sm = pd.DataFrame([{"statistic": k, **{kk: round(vv, 4) for kk, vv in v.items()}}
                       for k, v in s["size_matched_context"].items() if isinstance(v, dict)])
    louv = pd.DataFrame([{"species": LABEL[k], **row} for k in (DR, HS) for row in s["louvain"][k]])
    enr = pd.DataFrame([{"species": LABEL[k], **{c: v[c] for c in ("background_size", "modules_total", "modules_tested",
                                                                   "tests", "significant_tests",
                                                                   "modules_with_significant_term")}}
                        for k, v in s["enrichment"].items()])
    ov = pd.DataFrame(s["annotation_overlap"]["categories"])
    al = pd.DataFrame(s["alignment"])
    avail = pd.DataFrame([{"species": LABEL[k], **v} for k, v in s["evidence_availability"].items()])
    cats = pd.DataFrame(s["category_counts"])
    if not cats.empty:
        cats["order"] = cats["category"].map({c: i for i, c in enumerate(CATEGORY_ORDER)})
        cats = cats.sort_values("order").drop(columns="order")
        cats["definition"] = cats["category"].map(CATEGORY_DESCRIPTIONS)

    top_cols = ["rank", "dr_name", "hs_name", "category", "score", "evidence_coverage", "rank_q05", "rank_q95",
                "supporting_evidence", "contradictory_evidence"]
    eligible = ranked[ranked["eligible"]]
    top = eligible.head(cfg["ranking"]["report_top_n"])
    analogues = eligible[eligible["category"] == "candidate_analogue"].head(15)
    role_analogues = eligible[eligible["analogous_network_role"] & eligible["category"].isin(
        ["candidate_analogue", "shared_function_annotation", "similar_pathway_role"])].head(10)
    best_dr = (eligible.sort_values("score", ascending=False).groupby("dr_string_id").head(1)
               .sort_values("score", ascending=False))

    ws = s["weight_sensitivity"]
    top_key = f"top{cfg['ranking']['sensitivity']['top_k']}_overlap_median"
    loo = pd.DataFrame(ws["leave_one_out"])
    thr_rank = pd.DataFrame(s["ranking_threshold_sensitivity"])
    emb_lines = []
    for k, v in s["embedding"].items():
        if v.get("evaluated"):
            emb_lines.append(f"**{LABEL[k]}** (held-out edges, {cfg['embedding']['n_splits']} splits):\n\n"
                             + _table(pd.DataFrame(v["summary"])))
        else:
            emb_lines.append(f"**{LABEL[k]}**: not evaluated: {v['reason']}.")
    proxy = s["proxy_evaluation"]
    proxy_rows = pd.DataFrame([{"component": c, **v} for c, v in proxy["components"].items()])

    hom = homology_table.head(16)[["dr_name", "hit_id", "hit_name", "bitscore", "hit_in_human_network",
                                   "reciprocal_best_hit"]] if not homology_table.empty else pd.DataFrame()

    return "\n".join([
        f"# Phase 2 technical report: run `{s['run_id']}`",
        "",
        "Generated by `extremonet.analysis` from the run's own tables. Computational evidence only: nothing",
        "here demonstrates biological equivalence or experimental validity.",
        "",
        f"- Phase 1 input: run `{s['phase1']['run_id']}` (STRING {s['string_version']}); "
        f"{s['phase1']['outputs_verified']} outputs and {s['phase1']['raw_inputs_verified']} raw responses re-verified by SHA-256.",
        f"- Primary edge threshold: {s['primary_threshold']} (combined score). Master random seed: {cfg['random_seed']}.",
        "",
        "## 1. Baseline topology",
        _table(stats_rows),
        "",
        "STRING PPI-enrichment test (are the seeds more associated than a random protein set of the same size",
        "from the same genome?). P-values are as returned by STRING; 1e-64 appears to be its reporting floor.",
        "Because the seeds share GO annotations, which STRING's text-mining and database channels partly",
        "reflect, this enrichment is expected and is not evidence of a novel pathway:",
        "",
        _table(ppi),
        "",
        f"Size-matched context: {cfg['topology']['size_matched_samples']} random subsets of "
        f"{gs[DR]['nodes']} human seed proteins (descriptive context, not a hypothesis test):",
        "",
        _table(sm),
        "",
        "## 2. Communities (Louvain)",
        _table(louv, ["species", "resolution", "n_runs", "ari_mean", "ari_min", "modules_median",
                      "non_singleton_modules_medoid", "modularity_medoid", "medoid_seed"]),
        "",
        "Module enrichment (hypergeometric vs the species' seed-network background, BH-FDR within species):",
        "",
        _table(enr),
        "",
        "KEGG map membership is assigned through orthology (KO), so a bacterial protein can sit in a map",
        "named after a eukaryotic organelle (e.g. dra04146 'Peroxisome'); such terms are not organelle evidence.",
        "",
        f"{LABEL[DR]} modules:",
        "",
        _table(module_tables[DR], ["module", "size", "members", "significant_terms", "top_terms"]),
        "",
        f"{LABEL[HS]}: largest modules (full table in `homo_sapiens/modules.tsv`):",
        "",
        _table(module_tables[HS].sort_values("size", ascending=False).head(10),
               ["module", "size", "mean_coassignment", "significant_terms", "top_terms"]),
        "",
        "## 3. Evidence availability",
        _table(avail),
        "",
        f"Curated literature evidence: {s.get('curated_evidence', {}).get('rows', 0)} PubMed-verified findings for "
        f"{s.get('curated_evidence', {}).get('proteins', 0)} proteins (`{s.get('curated_evidence', {}).get('source')}`), "
        "including findings that contradict a simple protective role; see `candidates/protein_evidence.tsv`.",
        "",
        "## 4. Cross-species comparison",
        "### 4.1 Annotation vocabulary overlap (descriptive)",
        _table(ov, ["category", "dr_terms", "hs_terms", "shared_terms", "jaccard", "fraction_of_dr_terms_shared"]),
        "",
        _pruning_line(s.get("signature_pruning")),
        "",
        "### 4.2 Sequence homology (STRING best hits, D. radiodurans -> human)",
        _table(hom),
        "",
        f"Reciprocal best hits with both proteins inside the two networks: {s['rbh_pairs_within_networks']}.",
        "",
        "### 4.3 Network alignment (IsoRank)",
        _table(al),
        "",
        "`p_vs_random_mapping` only shows that the optimiser beats a random assignment. `p_vs_rewired_null`",
        "compares against degree-preserving rewired D. radiodurans graphs aligned the same way",
        f"({gs[DR]['edges']} source edges; `rewired_graphs_distinct` counts rewirings that differ from the",
        "observed graph). With a homology prior, conservation is partly built in: STRING transfers",
        "association evidence between orthologs of different organisms.",
        "",
        "Module correspondence (core-term similarity of non-singleton modules; top 10):",
        "",
        _table(pd.DataFrame(s["module_correspondence_top"])),
        "",
        "### 4.4 Experimentally confirmed physical interactions (IntAct)",
        *_intact_lines(s.get("experimental_interactions", {})),
        "",
        "### 4.5 Coverage variant: KO-projected D. radiodurans network",
        *_coverage_lines(s.get("coverage_variant", {})),
        "",
        "## 5. Candidate ranking",
        "Weights (config `ranking.weights`): " + ", ".join(f"{k}={v}" for k, v in cfg["ranking"]["weights"].items())
        + f". Missing components are excluded from each pair's denominator; pairs need evidence coverage >= "
          f"{cfg['ranking']['min_evidence_coverage']} to be ranked.",
        "",
        _table(cats),
        "",
        f"### 5.1 Top {cfg['ranking']['report_top_n']} pairs overall",
        _table(top, top_cols),
        "",
        "### 5.2 Candidate analogues (no detected homology, shared EC plus independent evidence)",
        _table(analogues, top_cols),
        "",
        "### 5.3 Non-homologous pairs with similar network roles (role similarity >= 0.75) and functional evidence",
        _table(role_analogues, top_cols),
        "",
        "### 5.4 Best-scoring human partner for each D. radiodurans protein",
        _table(best_dr, ["dr_name", "dr_uniprot", "hs_name", "hs_uniprot", "category", "score", "evidence_coverage",
                         "uncertainty"]),
        "",
        "## 6. Robustness",
        f"Dirichlet weight perturbation ({ws['dirichlet']['samples']} samples, concentration "
        f"{ws['dirichlet']['concentration']}): median Spearman vs base = {_fmt(ws['dirichlet']['spearman_median'])}, "
        f"5th percentile = {_fmt(ws['dirichlet']['spearman_q05'])}; median top-{cfg['ranking']['sensitivity']['top_k']} "
        f"overlap = {_fmt(ws['dirichlet'][top_key])}.",
        "",
        "Leave-one-component-out:",
        "",
        _table(loo),
        "",
        "Similarity-measure sensitivity (GO, pathway and module components recomputed with each measure):",
        "",
        _table(pd.DataFrame(s["similarity_measure_sensitivity"])),
        "",
        "STRING threshold sensitivity of the ranking (modules and roles recomputed at each threshold):",
        "",
        _table(thr_rank),
        "",
        "Network-level threshold sensitivity: `evaluation/threshold_sensitivity_networks.tsv`; figure "
        "`figures/threshold_sensitivity.png`.",
        "",
        "## 7. Embeddings (NetMF, within species only)",
        *emb_lines,
        "",
        "## 8. Evaluation against a proxy",
        f"Ground truth: {proxy['ground_truth']} ({proxy['positives']} positive pairs).",
        "",
        _table(proxy_rows),
        "",
        proxy.get("note", ""),
        "",
        "AUCs of annotation-derived components (`leakage_risk = True`: EC, GO, pathway, module) are inflated by",
        "construction; only `role` (topology) is free of annotation leakage. For annotation components, GO/EC",
        "annotations are often transferred between homologs. No accuracy for analogue prediction is claimed.",
        "",
        "## 9. Limitations",
        "- The *D. radiodurans* network has very few edges; topological, module and alignment statistics for it",
        "  have little power, and role/module components are not assessable for its isolated proteins.",
        "- Seed sets are GO-defined; annotation depth differs strongly between the organisms.",
        "- STRING edges are functional associations, not physical interactions; homology is best-hit based and",
        "  cannot exclude remote homology.",
        "- GO similarity uses propagated STRING annotations, an IC corpus dominated by human proteins and the",
        f"  `{cfg['similarity']['measure']}` normalisation; see the similarity-measure sensitivity table.",
        "- Perturbation evidence comes only from GO mutant-phenotype evidence codes; no phenotype database was used.",
        "- Literature evidence is limited to references attached to experimental GO annotations (QuickGO).",
        "",
    ])


def _coverage_lines(cov: dict) -> list[str]:
    if not cov.get("enabled"):
        return ["_not run_"]
    g = cov["graph_statistics"]
    lines = [
        f"D. radiodurans genes (original seeds included) in a KEGG KO group carried by a human seed: "
        f"{cov['projected_dr_genes']} (human seed KOs: {cov['human_seed_kos']}). Union with the "
        f"{cov['original_seeds']} original seeds: {cov['mapped_to_string']} proteins mapped to STRING "
        f"({cov['mapped_to_string'] - cov['original_seeds_mapped']} new; {cov['unmapped']} unmapped).",
        "Caveat: STRING transfers association evidence between orthologs of different organisms, so",
        "edge conservation between KO-orthologous pairs is partly built into the input scores.",
        f"At the primary threshold: {g['nodes']} nodes, {g['edges']} edges, {g['isolated_nodes']} isolated, "
        f"largest component {g['largest_component_nodes']} nodes, average clustering {_fmt(g['average_clustering'])}.",
        "",
        "Alignment of the variant to the human network (prior for alpha < 1: shared KEGG KO):",
        "",
        _table(pd.DataFrame(cov["alignment"])),
        "",
        f"{cov['note']}.",
    ]
    return lines


def _intact_lines(summary: dict) -> list[str]:
    if not summary:
        return ["_not run_"]
    rows = pd.DataFrame([{"species": LABEL.get(k, k), **v} for k, v in summary.items()])
    return [
        "STRING edges are functional associations. IntAct records with a physical interaction type and both",
        "partners among the seeds mark which of them are experimentally confirmed physical interactions.",
        "The `_binary` columns exclude co-complex association (MI:0914), which places two proteins in one",
        "complex without showing that they bind each other.",
        "",
        _table(rows, ["species", "accessions_queried", "physical_pairs_between_seeds", "string_edges",
                      "string_edges_confirmed", "fraction_string_edges_confirmed",
                      "string_edges_confirmed_binary", "fraction_string_edges_confirmed_binary",
                      "string_physical_subnetwork_edges", "string_physical_subnetwork_edges_confirmed",
                      "string_physical_subnetwork_edges_confirmed_binary",
                      "intact_pairs_not_in_string_graph", "self_interactions"]),
        "",
        "Zero records for a species means none are deposited in IntAct, not that the proteins do not interact.",
    ]


def _pruning_line(pruning: dict | None) -> str:
    if not pruning:
        return ""
    folds = pruning["common_fold_signatures"]
    items = [f"{k}: " + (", ".join(f"{s} ({v:.0%})" for s, v in d.items()) or "none") for k, d in folds.items()]
    return (f"Common-fold signatures (carried by > {pruning['max_fraction']:.0%} of the {pruning['proteins']} network "
            "proteins) count only as weak `distant_homolog` evidence: " + "; ".join(items) + ".")
