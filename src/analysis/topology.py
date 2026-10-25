"""Node- and graph-level topology, with explicit handling of disconnected graphs.

* Shortest-path measures use edge length = 1 / combined score (stronger association = shorter).
* Betweenness and PageRank are defined on disconnected graphs; closeness is computed with the
  Wasserman-Faust correction (networkx ``wf_improved``), i.e. scaled by component size.
* Cross-species comparisons use within-species percentile ranks over *connected* nodes, because
  raw centrality values depend on network size.
"""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
import pandas as pd


def add_distance(graph: nx.Graph) -> nx.Graph:
    g = graph.copy()
    for _, _, data in g.edges(data=True):
        data["distance"] = 1.0 / data["weight"]
    return g


def _percentile(values: pd.Series, mask: pd.Series) -> pd.Series:
    """Percentile rank (0-1, average ties) among masked rows; NaN elsewhere."""
    out = pd.Series(np.nan, index=values.index)
    if mask.sum() == 0:
        return out
    if mask.sum() == 1:
        out[mask] = 1.0
        return out
    out[mask] = values[mask].rank(pct=True, method="average")
    return out


def node_metrics(graph: nx.Graph, hub_pct: float = 0.8, bottleneck_pct: float = 0.8) -> pd.DataFrame:
    g = add_distance(graph)
    nodes = list(g.nodes)
    components = {n: i for i, comp in enumerate(sorted(nx.connected_components(g), key=len, reverse=True))
                  for n in comp}
    comp_size = {n: len(c) for c in nx.connected_components(g) for n in c}
    betweenness = nx.betweenness_centrality(g, weight="distance", normalized=True) if len(g) > 2 else {}
    closeness = nx.closeness_centrality(g, distance="distance", wf_improved=True) if len(g) > 1 else {}
    pagerank = nx.pagerank(g, weight="weight") if g.number_of_edges() else {n: 1 / len(g) for n in g}
    cores = nx.core_number(g) if g.number_of_edges() else {n: 0 for n in g}
    frame = pd.DataFrame({
        "string_id": nodes,
        "string_preferred_name": [g.nodes[n].get("string_preferred_name") for n in nodes],
        "degree": [g.degree(n) for n in nodes],
        "weighted_degree": [round(g.degree(n, weight="weight"), 6) for n in nodes],
        "clustering": [nx.clustering(g, n) for n in nodes],
        "weighted_clustering": [nx.clustering(g, n, weight="weight") for n in nodes],
        "betweenness": [betweenness.get(n, 0.0) for n in nodes],
        "closeness": [closeness.get(n, 0.0) for n in nodes],
        "pagerank": [pagerank[n] for n in nodes],
        "core_number": [cores[n] for n in nodes],
        "component": [components[n] for n in nodes],
        "component_size": [comp_size[n] for n in nodes],
    })
    connected = frame["degree"] > 0
    frame["is_isolated"] = ~connected
    for col in ("degree", "weighted_degree", "betweenness", "closeness", "pagerank", "clustering"):
        frame[f"{col}_pct"] = _percentile(frame[col], connected)
    is_hub = connected & (frame["degree_pct"] >= hub_pct)
    is_bottleneck = connected & (frame["betweenness"] > 0) & (frame["betweenness_pct"] >= bottleneck_pct)
    frame["is_hub"], frame["is_bottleneck"] = is_hub, is_bottleneck
    frame["hub_bottleneck_class"] = np.select(
        [~connected, is_hub & is_bottleneck, is_hub, is_bottleneck],
        ["isolated", "hub-bottleneck", "hub-nonbottleneck", "nonhub-bottleneck"], "nonhub-nonbottleneck")
    return frame.sort_values("string_id").reset_index(drop=True)


def module_roles(graph: nx.Graph, membership: dict[str, int]) -> pd.DataFrame:
    """Guimera & Amaral (2005) within-module degree z-score and participation coefficient."""
    rows = []
    module_nodes: dict[int, list[str]] = {}
    for n, m in membership.items():
        module_nodes.setdefault(m, []).append(n)
    within = {n: sum(1 for v in graph[n] if membership[v] == membership[n]) for n in graph}
    for m, members in module_nodes.items():
        k = np.array([within[n] for n in members], dtype=float)
        mu, sd = k.mean(), k.std()
        for n, kn in zip(members, k):
            degree = graph.degree(n)
            if degree == 0:
                participation = np.nan
            else:
                counts: dict[int, int] = {}
                for v in graph[n]:
                    counts[membership[v]] = counts.get(membership[v], 0) + 1
                participation = 1.0 - sum((c / degree) ** 2 for c in counts.values())
            rows.append({"string_id": n, "module": m, "module_size": len(members),
                         "within_module_degree": int(kn),
                         "within_module_z": (kn - mu) / sd if sd > 0 else 0.0,
                         "participation_coefficient": participation})
    return pd.DataFrame(rows).sort_values("string_id").reset_index(drop=True)


def graph_statistics(graph: nx.Graph) -> dict[str, Any]:
    n, m = graph.number_of_nodes(), graph.number_of_edges()
    degrees = np.array([d for _, d in graph.degree()], dtype=float)
    connected = [c for c in nx.connected_components(graph)]
    lcc = graph.subgraph(max(connected, key=len)).copy() if connected else nx.Graph()
    stats: dict[str, Any] = {
        "nodes": n, "edges": m,
        "density": nx.density(graph) if n > 1 else 0.0,
        "isolated_nodes": int((degrees == 0).sum()),
        "connected_nodes": int((degrees > 0).sum()),
        "connected_components": len(connected),
        "non_trivial_components": sum(1 for c in connected if len(c) > 1),
        "largest_component_nodes": lcc.number_of_nodes(),
        "largest_component_edges": lcc.number_of_edges(),
        "degree_mean": float(degrees.mean()) if n else None,
        "degree_median": float(np.median(degrees)) if n else None,
        "degree_max": int(degrees.max()) if n else None,
        "transitivity": nx.transitivity(graph),
        "average_clustering": nx.average_clustering(graph) if n else None,
    }
    stats["lcc_average_shortest_path_hops"] = (nx.average_shortest_path_length(lcc)
                                               if lcc.number_of_nodes() > 1 else None)
    stats["lcc_diameter_hops"] = nx.diameter(lcc) if lcc.number_of_nodes() > 1 else None
    try:
        stats["degree_assortativity"] = (float(nx.degree_assortativity_coefficient(graph))
                                         if m >= 3 else None)
    except (ValueError, ZeroDivisionError):
        stats["degree_assortativity"] = None
    if stats["degree_assortativity"] is not None and not np.isfinite(stats["degree_assortativity"]):
        stats["degree_assortativity"] = None
    return {k: (round(v, 6) if isinstance(v, float) else v) for k, v in stats.items()}


def degree_distribution(graph: nx.Graph) -> pd.DataFrame:
    degrees = pd.Series([d for _, d in graph.degree()])
    counts = degrees.value_counts().sort_index()
    frame = pd.DataFrame({"degree": counts.index, "count": counts.values})
    frame["fraction"] = frame["count"] / frame["count"].sum()
    frame["ccdf"] = frame["fraction"][::-1].cumsum()[::-1]
    return frame


def size_matched_null(reference: nx.Graph, k: int, samples: int, seed: int) -> pd.DataFrame:
    """Statistics of induced subgraphs on k uniformly sampled nodes of `reference` (the larger
    network), to put a smaller network's statistics in a size-matched context."""
    rng = np.random.default_rng(seed)
    nodes = np.array(sorted(reference.nodes))
    rows = []
    for i in range(samples):
        sub = reference.subgraph(rng.choice(nodes, size=k, replace=False))
        degrees = np.array([d for _, d in sub.degree()])
        comps = list(nx.connected_components(sub))
        rows.append({"sample": i, "edges": sub.number_of_edges(), "density": nx.density(sub),
                     "connected_nodes": int((degrees > 0).sum()),
                     "largest_component_nodes": max(len(c) for c in comps)})
    return pd.DataFrame(rows)


def empirical_position(null: pd.Series, observed: float) -> dict[str, float]:
    """Share of null samples below / at-or-above the observed value (descriptive context, not a
    test of a biological hypothesis)."""
    return {"observed": observed, "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "null_q05": float(null.quantile(0.05)), "null_q95": float(null.quantile(0.95)),
            "fraction_null_below": float((null < observed).mean()),
            "fraction_null_at_or_above": float((null >= observed).mean())}
