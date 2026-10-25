"""Module over-representation analysis.

For each module (size >= min_module_size) and each term annotating >= min_term_size background
proteins, a one-sided hypergeometric test asks whether the term is over-represented in the module
relative to the background. The background is the species' whole network node set: because the
seeds were themselves selected by oxidative-stress GO terms, a genome background would trivially
"enrich" those terms. P-values are corrected with Benjamini-Hochberg across *all* tests of a species.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def module_enrichment(membership: dict[str, int], annotations: pd.DataFrame, categories: list[str],
                      min_module_size: int = 3, min_term_size: int = 2) -> tuple[pd.DataFrame, dict]:
    background = set(membership)
    ann = annotations[annotations["category"].isin(categories) & annotations["string_id"].isin(background)]
    term_members = ann.groupby(["category", "term"])["string_id"].apply(set)
    descriptions = ann.drop_duplicates(["category", "term"]).set_index(["category", "term"])["description"]
    modules: dict[int, set[str]] = {}
    for node, m in membership.items():
        modules.setdefault(m, set()).add(node)
    N = len(background)
    rows = []
    tested_modules = [m for m, members in modules.items() if len(members) >= min_module_size]
    for m in tested_modules:
        members = modules[m]
        n = len(members)
        for (category, term), annotated in term_members.items():
            K = len(annotated)
            if K < min_term_size:
                continue
            hits = members & annotated
            k = len(hits)
            if k == 0:
                continue
            rows.append({
                "module": m, "module_size": n, "category": category, "term": term,
                "description": descriptions.get((category, term), ""),
                "k_in_module": k, "K_in_background": K, "N_background": N,
                "fold_enrichment": (k / n) / (K / N),
                "p_value": float(hypergeom.sf(k - 1, N, K, n)),
                "members": ";".join(sorted(hits)),
            })
    frame = pd.DataFrame(rows, columns=["module", "module_size", "category", "term", "description", "k_in_module",
                                        "K_in_background", "N_background", "fold_enrichment", "p_value", "members"])
    frame["fdr_bh"] = benjamini_hochberg(frame["p_value"].to_numpy()) if len(frame) else []
    summary = {"background_size": N, "modules_total": len(modules), "modules_tested": len(tested_modules),
               "tests": len(frame), "categories": categories,
               "annotated_background_nodes": int(ann["string_id"].nunique())}
    return frame.sort_values(["fdr_bh", "p_value", "module"]).reset_index(drop=True), summary
