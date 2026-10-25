"""Transparent, configurable ranking of cross-species candidate pairs.

Each (D. radiodurans protein, human protein) pair gets eight evidence components in [0, 1] or
missing. The score is the weighted mean over *available* components only:

    score = sum_i w_i s_i / sum_i w_i      (i over components that are not missing)

so missing evidence never lowers a score. ``evidence_coverage`` = sum of available weights is
reported with every score; pairs below ``min_evidence_coverage`` are kept in the table but not
ranked. Categories are assigned by explicit rules (``assign_category``), independently of the score.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .similarity import (
    HomologyTables, ProteinEvidence, curated_contradictions, ec_similarity, homology, information_content,
    role_similarity, simgic,
)

COMPONENTS = ["ec", "go", "pathway", "literature", "perturbation", "module", "role", "homology"]
ROLE_FEATURES = ["degree_pct", "betweenness_pct", "closeness_pct", "clustering", "participation_coefficient"]

CATEGORY_ORDER = ["ortholog_rbh", "ortholog_ko", "homolog", "candidate_analogue", "shared_function_annotation",
                  "similar_pathway_role", "speculative", "unsupported"]
CATEGORY_DESCRIPTIONS = {
    "ortholog_rbh": "reciprocal best STRING sequence hits (orthology candidate; not proof of equal function)",
    "ortholog_ko": "same KEGG Orthology group, not a reciprocal best hit (curated orthology; paralogs possible)",
    "homolog": "detectable homology (shared PANTHER family, one-way best hit, shared Pfam domain or shared "
               "InterPro entry)",
    "candidate_analogue": "no detected homology (no best hit, PANTHER, Pfam or InterPro signal); shared catalytic "
                          "activity (EC) plus independent pathway, GO or experimental evidence",
    "shared_function_annotation": "no detected homology; shared EC or GO function annotation only",
    "similar_pathway_role": "no detected homology; shared KEGG pathway maps without shared function annotation",
    "speculative": "no detected homology; only network-role or module similarity",
    "unsupported": "no component above the category thresholds",
}


def assign_category(row: dict[str, Any], go_thr: float, path_thr: float) -> str:
    tier = row["homology_tier"]
    if tier == "reciprocal_best_hit":
        return "ortholog_rbh"
    if tier == "shared_kegg_ko":
        return "ortholog_ko"
    if tier != "none_detected":
        return "homolog"
    ec = row["ec"] if row["ec"] is not None else 0.0
    go = row["go"] if row["go"] is not None else 0.0
    path = row["pathway"] if row["pathway"] is not None else 0.0
    both_experimental = row["literature"] == 1.0
    if ec >= 0.5 and (path >= path_thr or go >= go_thr or both_experimental):
        return "candidate_analogue"
    if ec >= 0.5 or go >= go_thr:
        return "shared_function_annotation"
    if path >= path_thr:
        return "similar_pathway_role"
    if (row["role"] or 0) >= 0.5 or (row["module"] or 0) >= 0.5:
        return "speculative"
    return "unsupported"


VALIDATION = {
    "ortholog_rbh": "Confirm orthology with a gene-tree analysis (e.g. eggNOG/OrthoDB families); test whether the "
                    "human gene complements a D. radiodurans {dr} deletion under H2O2, paraquat or ionising radiation.",
    "ortholog_ko": "Confirm the ortholog assignment with a gene tree; test whether the human gene complements a "
                   "D. radiodurans {dr} deletion under H2O2, paraquat or ionising radiation.",
    "homolog": "Resolve the homology relationship (ortholog vs paralog) phylogenetically, then test complementation "
               "of a D. radiodurans {dr} deletion by the human gene under oxidative challenge.",
    "candidate_analogue": "Assay the shared activity ({ec}) of both purified proteins under matched conditions; express "
                          "the human gene in a D. radiodurans {dr} deletion strain and measure survival under H2O2/paraquat; "
                          "compare folds structurally to confirm that the proteins are not remote homologs.",
    "shared_function_annotation": "Check whether the shared annotations rest on experimental or inferred evidence; "
                                  "measure the shared activity directly for both proteins.",
    "similar_pathway_role": "Measure co-regulation of both genes with their pathway partners under oxidative stress "
                            "(transcriptomics/proteomics) in each organism.",
    "speculative": "Low priority: obtain direct functional evidence (activity assay or mutant phenotype) before any "
                   "comparison.",
    "unsupported": "Not recommended for follow-up on current evidence.",
}


def score_matrix(values: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """values: pairs x components with NaN for missing. Returns (score, coverage)."""
    available = ~np.isnan(values)
    w = np.where(available, weights[None, :], 0.0)
    coverage = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        score = np.where(coverage > 0, np.nansum(np.where(available, values, 0.0) * w, axis=1) / coverage, np.nan)
    return score, coverage


def build_pairs(dr: dict[str, ProteinEvidence], hs: dict[str, ProteinEvidence],
                dr_metrics: pd.DataFrame, hs_metrics: pd.DataFrame,
                dr_modules: dict[str, int], hs_modules: dict[str, int],
                dr_cores: dict[int, set[str]], hs_cores: dict[int, set[str]],
                homology_tables: HomologyTables, min_bitscore: float, measure: str = "ochiai") -> pd.DataFrame:
    go_ic = information_content([p.go_terms for p in (*dr.values(), *hs.values())])
    path_ic = information_content([p.pathways for p in (*dr.values(), *hs.values())])
    module_ic = information_content([*dr_cores.values(), *hs_cores.values()])
    dr_m = dr_metrics.set_index("string_id")
    hs_m = hs_metrics.set_index("string_id")
    rows = []
    for a in sorted(dr):
        pa, ma = dr[a], dr_m.loc[a] if a in dr_m.index else None
        core_a = dr_cores.get(dr_modules.get(a))
        for b in sorted(hs):
            pb, mb = hs[b], hs_m.loc[b] if b in hs_m.index else None
            ec, shared_ec = ec_similarity(pa.ec, pb.ec)
            hom, tier, hom_detail = homology(pa, pb, homology_tables, min_bitscore)
            literature = ((pa.literature_level + pb.literature_level) / 2
                          if pa.literature_level is not None and pb.literature_level is not None else None)
            # mean over the proteins that HAVE perturbation evidence: absence is never scored
            pert = [s for s in (pa.perturbation_score, pb.perturbation_score) if s is not None]
            core_b = hs_cores.get(hs_modules.get(b))
            negated_conflict = sorted((pa.negated_go & pb.go_terms) | (pb.negated_go & pa.go_terms))
            rows.append({
                "dr_string_id": a, "dr_name": pa.name, "dr_uniprot": pa.accession,
                "hs_string_id": b, "hs_name": pb.name, "hs_uniprot": pb.accession,
                "ec": ec, "go": simgic(pa.go_terms, pb.go_terms, go_ic, measure),
                "pathway": simgic(pa.pathways, pb.pathways, path_ic, measure),
                "literature": literature,
                "perturbation": sum(pert) / len(pert) if pert else None,
                "perturbation_evidence_for": len(pert),
                "module": simgic(core_a, core_b, module_ic, measure) if core_a and core_b else None,
                "role": role_similarity(ma, mb, ROLE_FEATURES),
                "homology": hom, "homology_tier": tier,
                "shared_ec": ";".join(shared_ec), "dr_ec": ";".join(sorted(pa.ec)), "hs_ec": ";".join(sorted(pb.ec)),
                "shared_panther": ";".join(hom_detail["shared_panther"]),
                "shared_pfam": ";".join(hom_detail["shared_pfam"]),
                "shared_interpro": ";".join(hom_detail["shared_interpro"]),
                "shared_ko": ";".join(hom_detail["shared_ko"]),
                "dr_ko": ";".join(sorted(pa.ko)), "hs_ko": ";".join(sorted(pb.ko)),
                "ko_conflict": hom_detail["ko_conflict"],
                "curated_contradictions": " | ".join(curated_contradictions(pa) + curated_contradictions(pb)),
                "shared_pathways": ";".join(sorted(pa.pathways & pb.pathways)),
                "dr_best_human_hit": hom_detail["best_hit_dr_to_hs"],
                "dr_best_hit_bitscore": hom_detail["bitscore_dr_to_hs"],
                "dr_isolated": bool(ma["is_isolated"]) if ma is not None else None,
                "hs_isolated": bool(mb["is_isolated"]) if mb is not None else None,
                "dr_evidence_codes": ";".join(sorted(pa.evidence_codes)),
                "hs_evidence_codes": ";".join(sorted(pb.evidence_codes)),
                "references": ";".join(sorted(pa.experimental_references | pb.experimental_references)),
                "negated_go_conflicts": ";".join(negated_conflict),
            })
    return pd.DataFrame(rows)


def contradictions(row: pd.Series) -> str:
    items = []
    if row["ec"] == 0.0:
        items.append(f"different EC numbers ({row['dr_ec']} vs {row['hs_ec']})")
    if row["go"] == 0.0:
        items.append("no shared GO BP/MF term")
    if row["role"] is not None and not pd.isna(row["role"]) and row["role"] < 0.3:
        items.append(f"dissimilar network roles (role similarity {row['role']:.2f})")
    if row["negated_go_conflicts"]:
        items.append(f"NOT-qualified GO annotation conflicts: {row['negated_go_conflicts']}")
    if row.get("ko_conflict") is True and row["homology_tier"] == "reciprocal_best_hit":
        items.append(f"reciprocal best hit but different KEGG Orthology groups ({row['dr_ko']} vs {row['hs_ko']})")
    if isinstance(row.get("curated_contradictions"), str) and row["curated_contradictions"]:
        items.append(f"curated perturbation evidence against a simple protective role: {row['curated_contradictions']}")
    return "; ".join(items)


def supporting(row: pd.Series) -> str:
    items = []
    if row["shared_ec"]:
        items.append(f"shared EC {row['shared_ec']}")
    if row["go"] is not None and not pd.isna(row["go"]) and row["go"] > 0:
        items.append(f"GO similarity {row['go']:.2f}")
    if row["shared_pathways"]:
        items.append(f"shared KEGG maps {row['shared_pathways']}")
    if row["homology_tier"] != "none_detected":
        items.append(f"homology: {row['homology_tier']}")
    if row["literature"] == 1.0:
        items.append("experimental GO evidence for both proteins")
    if row["perturbation"] is not None and not pd.isna(row["perturbation"]):
        who = "both proteins" if row.get("perturbation_evidence_for") == 2 else "one protein"
        items.append(f"perturbation evidence for {who} (score {row['perturbation']:.2f})")
    if row["module"] is not None and not pd.isna(row["module"]) and row["module"] > 0:
        items.append(f"module profile similarity {row['module']:.2f}")
    if row["role"] is not None and not pd.isna(row["role"]):
        items.append(f"role similarity {row['role']:.2f}")
    return "; ".join(items)


def uncertainty_notes(row: pd.Series) -> str:
    notes = []
    if row["dr_isolated"]:
        notes.append("D. radiodurans protein has no edge at this threshold: role/module not assessable")
    if row["dr_evidence_codes"] and set(row["dr_evidence_codes"].split(";")) <= {"IEA"}:
        notes.append("D. radiodurans GO evidence is electronic only (IEA)")
    if row["homology_tier"] == "none_detected":
        notes.append("no homology detected (best hit, PANTHER, Pfam, InterPro); very remote homology is not excluded")
    if row["evidence_coverage"] < 0.6:
        notes.append(f"low evidence coverage ({row['evidence_coverage']:.2f})")
    return "; ".join(notes)


def rank_pairs(pairs: pd.DataFrame, weights: dict[str, float], min_coverage: float,
               go_thr: float, path_thr: float) -> pd.DataFrame:
    frame = pairs.copy()
    w = np.array([weights[c] for c in COMPONENTS], dtype=float)
    values = frame[COMPONENTS].astype(float).to_numpy()
    frame["score"], frame["evidence_coverage"] = score_matrix(values, w)
    frame["components_available"] = (~np.isnan(values)).sum(axis=1)
    frame["eligible"] = frame["evidence_coverage"] >= min_coverage
    frame["category"] = [assign_category(r, go_thr, path_thr) for r in frame.to_dict("records")]
    frame["analogous_network_role"] = (frame["homology_tier"] == "none_detected") & (frame["role"].astype(float) >= 0.75)
    # round before sorting and break ties by identifiers, so the order is reproducible bit-for-bit
    frame["score"] = frame["score"].round(12)
    frame = frame.sort_values(["eligible", "score", "evidence_coverage", "dr_string_id", "hs_string_id"],
                              ascending=[False, False, False, True, True], kind="mergesort")
    frame["rank"] = np.where(frame["eligible"], np.arange(1, len(frame) + 1), np.nan)
    frame["supporting_evidence"] = frame.apply(supporting, axis=1)
    frame["contradictory_evidence"] = frame.apply(contradictions, axis=1)
    frame["uncertainty"] = frame.apply(uncertainty_notes, axis=1)
    frame["proposed_validation"] = [
        VALIDATION[c].format(dr=r["dr_name"], ec=r["shared_ec"] or "shared activity")
        for c, r in zip(frame["category"], frame.to_dict("records"))]
    return frame.reset_index(drop=True)


def weight_sensitivity(pairs: pd.DataFrame, ranked: pd.DataFrame, weights: dict[str, float],
                       samples: int, concentration: float, top_k: int, seed: int) -> dict[str, Any]:
    """Dirichlet perturbation of the weights and leave-one-component-out, on the eligible pairs."""
    eligible = ranked[ranked["eligible"]]
    key = list(zip(eligible["dr_string_id"], eligible["hs_string_id"]))
    lookup = pairs.set_index(["dr_string_id", "hs_string_id"]).loc[key]
    values = lookup[COMPONENTS].astype(float).to_numpy()
    base_w = np.array([weights[c] for c in COMPONENTS], dtype=float)
    base_score, _ = score_matrix(values, base_w)
    base_top = set(np.argsort(-base_score, kind="stable")[:top_k])
    rng = np.random.default_rng(seed)
    rhos, overlaps, ranks = [], [], []
    for _ in range(samples):
        w = rng.dirichlet(concentration * base_w / base_w.sum())
        s, _ = score_matrix(values, w)
        rhos.append(spearmanr(base_score, s).statistic)
        overlaps.append(len(base_top & set(np.argsort(-s, kind="stable")[:top_k])) / top_k)
        ranks.append(pd.Series(-s).rank(method="min").to_numpy())
    ranks_arr = np.array(ranks)
    rank_q = pd.DataFrame({"dr_string_id": eligible["dr_string_id"].to_numpy(),
                           "hs_string_id": eligible["hs_string_id"].to_numpy(),
                           "rank_q05": np.quantile(ranks_arr, 0.05, axis=0),
                           "rank_median": np.median(ranks_arr, axis=0),
                           "rank_q95": np.quantile(ranks_arr, 0.95, axis=0)})
    loo = []
    for i, comp in enumerate(COMPONENTS):
        w = base_w.copy()
        w[i] = 0.0
        s, _ = score_matrix(values, w)
        loo.append({"dropped_component": comp,
                    "spearman_vs_base": float(spearmanr(base_score, s, nan_policy="omit").statistic),
                    f"top{top_k}_overlap": len(base_top & set(np.argsort(-np.nan_to_num(s, nan=-1), kind="stable")[:top_k])) / top_k})
    return {
        "dirichlet": {"samples": samples, "concentration": concentration, "seed": seed,
                      "spearman_median": float(np.median(rhos)), "spearman_q05": float(np.quantile(rhos, 0.05)),
                      f"top{top_k}_overlap_median": float(np.median(overlaps)),
                      f"top{top_k}_overlap_q05": float(np.quantile(overlaps, 0.05))},
        "leave_one_out": loo, "rank_intervals": rank_q, "eligible_pairs": len(eligible),
    }


def top_k_overlap(a: pd.DataFrame, b: pd.DataFrame, k: int) -> float:
    ka = set(zip(*[a[a["eligible"]].head(k)[c] for c in ("dr_string_id", "hs_string_id")]))
    kb = set(zip(*[b[b["eligible"]].head(k)[c] for c in ("dr_string_id", "hs_string_id")]))
    return len(ka & kb) / k if k else float("nan")
