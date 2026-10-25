"""Louvain community detection with recorded seeds and stability assessment.

Louvain (Blondel et al. 2008) as implemented in ``networkx.community.louvain_communities``,
weighted by STRING combined score. Isolated nodes form singleton modules. Stability across random
seeds is summarised by the pairwise adjusted Rand index (ARI) and a node co-assignment matrix; the
reported partition is the medoid run (highest mean ARI to all other runs at that resolution).
"""

from __future__ import annotations

from itertools import combinations
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd


def louvain_partition(graph: nx.Graph, resolution: float, seed: int, weight: str = "weight") -> dict[str, int]:
    communities = nx.community.louvain_communities(graph, weight=weight, resolution=resolution, seed=seed)
    # deterministic labels: by decreasing size, then smallest member id
    ordered = sorted(communities, key=lambda c: (-len(c), min(c)))
    return {node: i for i, comm in enumerate(ordered) for node in comm}


def adjusted_rand_index(a: list[int], b: list[int]) -> float:
    """Hubert & Arabie (1985) ARI from the contingency table."""
    a_arr, b_arr = np.asarray(a), np.asarray(b)
    _, a_idx = np.unique(a_arr, return_inverse=True)
    _, b_idx = np.unique(b_arr, return_inverse=True)
    table = np.zeros((a_idx.max() + 1, b_idx.max() + 1))
    np.add.at(table, (a_idx, b_idx), 1)

    def comb2(x: np.ndarray) -> float:
        return float((x * (x - 1) / 2).sum())

    n = len(a_arr)
    sum_cells, sum_a, sum_b = comb2(table), comb2(table.sum(1)), comb2(table.sum(0))
    total = n * (n - 1) / 2
    expected = sum_a * sum_b / total if total else 0.0
    maximum = (sum_a + sum_b) / 2
    if maximum == expected:
        return 1.0
    return (sum_cells - expected) / (maximum - expected)


def modularity(graph: nx.Graph, membership: dict[str, int], resolution: float = 1.0) -> float | None:
    if graph.number_of_edges() == 0:
        return None
    groups: dict[int, set[str]] = {}
    for n, m in membership.items():
        groups.setdefault(m, set()).add(n)
    return float(nx.community.modularity(graph, groups.values(), weight="weight", resolution=resolution))


def stability(graph: nx.Graph, resolution: float, seeds: list[int]) -> dict[str, Any]:
    nodes = sorted(graph.nodes)
    runs = [louvain_partition(graph, resolution, s) for s in seeds]
    labels = [[r[n] for n in nodes] for r in runs]
    aris = [adjusted_rand_index(labels[i], labels[j]) for i, j in combinations(range(len(runs)), 2)]
    mean_ari_per_run = []
    for i in range(len(runs)):
        others = [adjusted_rand_index(labels[i], labels[j]) for j in range(len(runs)) if j != i]
        mean_ari_per_run.append(float(np.mean(others)) if others else 1.0)
    medoid = int(np.argmax(mean_ari_per_run))
    lab = np.array(labels)
    coassign = (lab[:, :, None] == lab[:, None, :]).mean(axis=0)
    n_modules = [len(set(r.values())) for r in runs]
    return {
        "resolution": resolution, "seeds": seeds, "runs": runs, "medoid_index": medoid,
        "medoid_seed": seeds[medoid], "partition": runs[medoid], "nodes": nodes,
        "coassignment": coassign,
        "summary": {
            "resolution": resolution, "n_runs": len(runs),
            "ari_mean": float(np.mean(aris)) if aris else 1.0,
            "ari_min": float(np.min(aris)) if aris else 1.0,
            "modules_median": float(np.median(n_modules)),
            "modules_min": int(min(n_modules)), "modules_max": int(max(n_modules)),
            "non_singleton_modules_medoid": sum(1 for m in set(runs[medoid].values())
                                                if sum(1 for v in runs[medoid].values() if v == m) > 1),
            "modularity_medoid": modularity(graph, runs[medoid], resolution),
            "medoid_seed": seeds[medoid],
        },
    }


def node_stability(result: dict[str, Any]) -> pd.DataFrame:
    """Mean co-assignment frequency of each node with the other members of its medoid module
    (1 = always together; NaN for singletons)."""
    nodes, part, co = result["nodes"], result["partition"], result["coassignment"]
    index = {n: i for i, n in enumerate(nodes)}
    rows = []
    for n in nodes:
        mates = [index[v] for v in nodes if part[v] == part[n] and v != n]
        rows.append({"string_id": n, "module": part[n],
                     "coassignment_with_module": float(co[index[n], mates].mean()) if mates else np.nan})
    return pd.DataFrame(rows)
