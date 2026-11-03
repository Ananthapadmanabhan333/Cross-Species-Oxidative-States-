"""NetMF node embeddings (Qiu et al. 2018) and leakage-free link-prediction evaluation.

NetMF is the closed-form matrix that DeepWalk implicitly factorises:
    M = vol(G) / (b * T) * (sum_{r=1..T} P^r) D^-1,  embedding = SVD_d(log(max(M, 1)))
It is deterministic given the graph, needs only numpy, and avoids random-walk sampling noise.

Embeddings are computed and evaluated *within one species only*. Embeddings of two independently
factorised graphs live in unrelated coordinate systems, so no cross-species embedding comparison is
made. For evaluation, held-out edges are removed before the embedding is computed (no leakage), and
the embedding is compared with classical neighbourhood baselines on the same split.
"""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
import pandas as pd


def netmf_embedding(graph: nx.Graph, dim: int, window: int, negative: float) -> tuple[np.ndarray, list[str]]:
    nodes = sorted(graph.nodes)
    a = nx.to_numpy_array(graph, nodelist=nodes, weight="weight")
    deg = a.sum(axis=1)
    vol = deg.sum()
    inv = np.where(deg > 0, 1.0 / np.where(deg > 0, deg, 1), 0.0)
    p = inv[:, None] * a
    acc, power = np.zeros_like(a), np.eye(len(nodes))
    for _ in range(window):
        power = power @ p
        acc += power
    m = (vol / (negative * window)) * acc * inv[None, :]
    log_m = np.log(np.maximum(m, 1.0))
    u, s, _ = np.linalg.svd(log_m)
    k = min(dim, len(nodes))
    return u[:, :k] * np.sqrt(s[:k]), nodes


def _unit_rows(emb: np.ndarray, rel_tol: float = 1e-9) -> np.ndarray:
    """Normalise embedding rows to unit length. Rows of nodes without signal (e.g. isolated in the
    training graph) are zero only up to SVD rounding (~1e-15); normalising that noise would give
    arbitrary, hardware-dependent cosines, so such rows are set exactly to zero (cosine 0)."""
    norm = np.linalg.norm(emb, axis=1)
    keep = norm > rel_tol * norm.max() if norm.size and norm.max() > 0 else np.zeros(len(norm), bool)
    unit = np.zeros_like(emb)
    unit[keep] = emb[keep] / norm[keep, None]
    return unit


def _auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC with ties counted as 1/2."""
    scores = np.concatenate([pos, neg])
    ranks = pd.Series(scores).rank(method="average").to_numpy()
    r_pos = ranks[: len(pos)].sum()
    return float((r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def link_prediction(graph: nx.Graph, dim: int, window: int, negative: float, test_fraction: float,
                    n_splits: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    edges = sorted(tuple(sorted(e)) for e in graph.edges)
    nodes = sorted(graph.nodes)
    rows = []
    for split in range(n_splits):
        perm = rng.permutation(len(edges))
        n_test = max(1, int(round(test_fraction * len(edges))))
        test = [edges[i] for i in perm[:n_test]]
        train = graph.copy()
        train.remove_edges_from(test)
        negatives: set[tuple[str, str]] = set()
        while len(negatives) < n_test:
            u, v = rng.choice(nodes, size=2, replace=False)
            pair = tuple(sorted((u, v)))
            if not graph.has_edge(*pair):
                negatives.add(pair)
        neg = sorted(negatives)
        emb, order = netmf_embedding(train, dim, window, negative)
        index = {n: i for i, n in enumerate(order)}
        unit = _unit_rows(emb)

        def cosine(pairs):
            # Rounded so that cosines equal up to floating-point noise tie exactly in the AUC.
            return np.round(np.array([unit[index[u]] @ unit[index[v]] for u, v in pairs]), 10)

        def classical(fn, pairs):
            return np.array([s for _, _, s in fn(train, pairs)], dtype=float)

        methods = {
            "netmf_cosine": cosine,
            "common_neighbours": lambda ps: np.array([len(list(nx.common_neighbors(train, u, v))) for u, v in ps], float),
            "jaccard": lambda ps: classical(nx.jaccard_coefficient, ps),
            "adamic_adar": lambda ps: classical(_safe_adamic_adar, ps),
            "preferential_attachment": lambda ps: classical(nx.preferential_attachment, ps),
        }
        for name, fn in methods.items():
            rows.append({"split": split, "method": name, "auc": _auc(fn(test), fn(neg)),
                         "n_test_edges": n_test, "n_negatives": len(neg)})
    return pd.DataFrame(rows)


def _safe_adamic_adar(graph: nx.Graph, pairs):
    """Adamic-Adar that skips degree-1 common neighbours (log 1 = 0 would divide by zero)."""
    for u, v in pairs:
        score = sum(1.0 / np.log(graph.degree(w)) for w in nx.common_neighbors(graph, u, v) if graph.degree(w) > 1)
        yield u, v, score


def summarise_link_prediction(frame: pd.DataFrame) -> list[dict[str, Any]]:
    grouped = frame.groupby("method")["auc"].agg(["mean", "std", "min", "max", "count"]).reset_index()
    return grouped.sort_values("mean", ascending=False).round(4).to_dict("records")
