"""IsoRank global network alignment (Singh, Xu & Berger 2008) with an optional homology prior.

R = alpha * A1n R A2n^T + (1 - alpha) * E, iterated to convergence (R normalised to sum 1), where
Ain is the column-normalised weighted adjacency and E is the prior (STRING homology bitscores,
normalised). alpha = 1 uses topology alone. A one-to-one alignment is read off greedily from R.

Edge conservation (fraction of D. radiodurans edges mapped onto human edges) is compared with two
nulls:
  * random one-to-one mappings - a weak baseline: it only shows that the optimiser beats chance,
    because a dense target network lets almost any small graph be embedded well;
  * degree-preserving rewired copies of the D. radiodurans graph aligned with the *same* procedure -
    asks whether the observed topology aligns better than graphs with the same degree sequence.
With very few D. radiodurans edges both tests have little power; results are reported as computed.
"""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
import pandas as pd


def _column_normalised(graph: nx.Graph, nodes: list[str]) -> np.ndarray:
    a = nx.to_numpy_array(graph, nodelist=nodes, weight="weight")
    col = a.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(col > 0, a / col, 0.0)


def isorank(g1: nx.Graph, g2: nx.Graph, prior: np.ndarray | None, alpha: float,
            max_iter: int = 100, tol: float = 1e-9) -> tuple[np.ndarray, list[str], list[str], dict]:
    n1, n2 = sorted(g1.nodes), sorted(g2.nodes)
    a1, a2 = _column_normalised(g1, n1), _column_normalised(g2, n2)
    if prior is None or prior.sum() == 0:
        if alpha < 1.0:
            raise ValueError("alpha < 1 needs a non-empty prior")
        e = np.full((len(n1), len(n2)), 1.0 / (len(n1) * len(n2)))
    else:
        e = prior / prior.sum()
    r = e.copy()
    delta, it = np.inf, 0
    for it in range(1, max_iter + 1):
        nxt = alpha * (a1 @ r @ a2.T) + (1 - alpha) * e
        total = nxt.sum()
        nxt = nxt / total if total > 0 else e.copy()
        delta = float(np.abs(nxt - r).sum())
        r = nxt
        if delta < tol:
            break
    return r, n1, n2, {"iterations": it, "final_l1_change": delta, "converged": delta < tol}


def greedy_matching(r: np.ndarray, n1: list[str], n2: list[str]) -> list[tuple[str, str, float]]:
    order = np.dstack(np.unravel_index(np.argsort(-r, axis=None, kind="stable"), r.shape))[0]
    used1, used2, pairs = set(), set(), []
    for i, j in order:
        if i in used1 or j in used2:
            continue
        used1.add(i), used2.add(j)
        pairs.append((n1[i], n2[j], float(r[i, j])))
        if len(used1) == len(n1) or len(used2) == len(n2):
            break
    return pairs


def edge_conservation(g1: nx.Graph, g2: nx.Graph, mapping: dict[str, str]) -> tuple[int, int]:
    conserved = sum(1 for u, v in g1.edges if g2.has_edge(mapping[u], mapping[v]))
    return conserved, g1.number_of_edges()


def conservation_null(g1: nx.Graph, g2: nx.Graph, n_perm: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n1, n2 = sorted(g1.nodes), np.array(sorted(g2.nodes))
    out = np.empty(n_perm)
    for i in range(n_perm):
        targets = rng.choice(n2, size=len(n1), replace=False)
        out[i] = edge_conservation(g1, g2, dict(zip(n1, targets)))[0]
    return out


def rewired_null(g1: nx.Graph, g2: nx.Graph, prior: np.ndarray | None, alpha: float, n_samples: int, seed: int,
                 max_iter: int, tol: float) -> tuple[np.ndarray, int]:
    """Edge conservation of IsoRank alignments of degree-preserving rewired copies of g1."""
    rng = np.random.default_rng(seed)
    observed_edges = {tuple(sorted(e)) for e in g1.edges}
    values, distinct = [], 0
    for _ in range(n_samples):
        g = g1.copy()
        m = g.number_of_edges()
        try:
            nx.double_edge_swap(g, nswap=10 * m, max_tries=1000 * m, seed=int(rng.integers(2 ** 31)))
        except (nx.NetworkXAlgorithmError, nx.NetworkXError):
            pass  # small graphs admit few swaps; keep whatever rewiring succeeded
        if {tuple(sorted(e)) for e in g.edges} != observed_edges:
            distinct += 1
        r, n1, n2, _ = isorank(g, g2, prior, alpha, max_iter, tol)
        mapping = {a: b for a, b, _ in greedy_matching(r, n1, n2)}
        values.append(edge_conservation(g, g2, mapping)[0])
    return np.array(values, dtype=float), distinct


def homology_prior(n1: list[str], n2: list[str], dr_to_hs: pd.DataFrame, hs_to_dr: pd.DataFrame) -> np.ndarray:
    idx1, idx2 = {n: i for i, n in enumerate(n1)}, {n: j for j, n in enumerate(n2)}
    prior = np.zeros((len(n1), len(n2)))
    for frame, flip in ((dr_to_hs, False), (hs_to_dr, True)):
        for r in frame.itertuples():
            a, b = (r.hit_id, r.query_id) if flip else (r.query_id, r.hit_id)
            if a in idx1 and b in idx2:
                prior[idx1[a], idx2[b]] = max(prior[idx1[a], idx2[b]], float(r.bitscore))
    return prior


def run_alignment(g1: nx.Graph, g2: nx.Graph, dr_to_hs: pd.DataFrame, hs_to_dr: pd.DataFrame, alphas: list[float],
                  max_iter: int, tol: float, n_perm: int, seed: int, rbh_pairs: set[tuple[str, str]],
                  n_rewired: int = 100) -> dict[str, Any]:
    n1, n2 = sorted(g1.nodes), sorted(g2.nodes)
    if len(n1) > len(n2):
        return {"summary": [{"skipped": f"source network ({len(n1)} nodes) is larger than target ({len(n2)}); "
                                        "a one-to-one alignment of every source node is impossible"}],
                "pairs": pd.DataFrame(), "rbh_pairs_within_networks": len(rbh_pairs)}
    prior = homology_prior(n1, n2, dr_to_hs, hs_to_dr)
    null = conservation_null(g1, g2, n_perm, seed)
    results, tables = [], []
    for alpha in alphas:
        if alpha < 1.0 and prior.sum() == 0:
            results.append({"alpha": alpha, "skipped": "no homology prior between network members"})
            continue
        r, _, _, info = isorank(g1, g2, prior if alpha < 1.0 else None, alpha, max_iter, tol)
        pairs = greedy_matching(r, n1, n2)
        mapping = {a: b for a, b, _ in pairs}
        conserved, total = edge_conservation(g1, g2, mapping)
        p_emp = (1 + int((null >= conserved).sum())) / (1 + n_perm)
        rewired, distinct = rewired_null(g1, g2, prior if alpha < 1.0 else None, alpha, n_rewired, seed + 1,
                                         max_iter, tol)
        p_rewired = (1 + int((rewired >= conserved).sum())) / (1 + n_rewired)
        results.append({
            "alpha": alpha, **info, "aligned_pairs": len(pairs),
            "edges_conserved": conserved, "edges_in_dr": total,
            "edge_conservation": conserved / total if total else None,
            "null_mean_conserved": float(null.mean()), "null_q95_conserved": float(np.quantile(null, 0.95)),
            "p_vs_random_mapping": p_emp,
            "rewired_null_samples": n_rewired, "rewired_graphs_distinct": distinct,
            "rewired_null_mean_conserved": float(rewired.mean()),
            "p_vs_rewired_null": p_rewired,
            "aligned_pairs_that_are_rbh": sum(1 for a, b, _ in pairs if (a, b) in rbh_pairs),
        })
        tables.append(pd.DataFrame([{"alpha": alpha, "dr_string_id": a, "hs_string_id": b, "isorank_score": s,
                                     "is_rbh": (a, b) in rbh_pairs,
                                     "dr_edge_count": g1.degree(a)} for a, b, s in pairs]))
    return {"summary": results, "pairs": pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(),
            "rbh_pairs_within_networks": len(rbh_pairs)}
