"""Construct documented, weighted NetworkX graphs from QC-retained STRING records.

Graph model (one undirected graph per species):
  * node   = STRING protein id (e.g. ``243230.DR_1998``); attributes from the node table
  * edge   = STRING functional association retained by QC
  * weight = STRING combined score (0-1), also stored as ``combined_score``
  * channel scores (``experimental_score`` ...) are kept on every edge
  * ``in_physical_subnetwork`` / ``physical_score``: the pair also appears in STRING's physical
    subnetwork at the same threshold. ``has_experiments_channel_score``: experiments-channel score above
    the configured threshold. Neither flag means "experimentally confirmed physical interaction".
"""

from __future__ import annotations

import math
from typing import Any

import networkx as nx
import pandas as pd

from ..preprocessing.parsers import STRING_CHANNELS

CHANNEL_COLUMNS = [f"{name}_score" for name in STRING_CHANNELS.values()]
EDGE_COLUMNS = ["node_u", "node_v", "preferred_name_a", "preferred_name_b", "combined_score",
                *CHANNEL_COLUMNS, "in_physical_subnetwork", "physical_score", "has_experiments_channel_score",
                "evidence_source", "source_record"]


def annotate_edges(functional: pd.DataFrame, physical: pd.DataFrame | None,
                   experimental_threshold: float, source_label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach physical-subnetwork membership and experimental-evidence flags to functional edges.

    Returns the annotated functional edges and the physical edges that have *no* functional
    counterpart (expected to be empty; reported rather than silently added or dropped).
    """
    edges = functional.copy()
    if physical is not None and not physical.empty:
        phys = physical[["node_u", "node_v", "combined_score"]].rename(columns={"combined_score": "physical_score"})
        edges = edges.merge(phys, on=["node_u", "node_v"], how="left")
        keys = set(zip(edges["node_u"], edges["node_v"]))
        orphan = physical[[k not in keys for k in zip(physical["node_u"], physical["node_v"])]].copy()
    else:
        edges["physical_score"] = float("nan")
        orphan = pd.DataFrame(columns=physical.columns if physical is not None else [])
    edges["in_physical_subnetwork"] = edges["physical_score"].notna()
    edges["has_experiments_channel_score"] = edges["experimental_score"].fillna(0) > experimental_threshold
    edges["evidence_source"] = source_label
    return edges[EDGE_COLUMNS].sort_values(["node_u", "node_v"]).reset_index(drop=True), orphan


def _clean(value: Any) -> Any:
    """GraphML cannot store None/NaN; drop such attributes instead of writing placeholders."""
    if value is None or value is pd.NA:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if hasattr(value, "item"):  # numpy scalar -> python scalar
        return value.item()
    return value


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame, keep_isolated: bool = True,
                graph_attrs: dict[str, Any] | None = None) -> nx.Graph:
    graph = nx.Graph(**(graph_attrs or {}))
    node_ids = set(nodes["string_id"])
    connected = set(edges["node_u"]) | set(edges["node_v"])
    missing = connected - node_ids
    if missing:
        raise ValueError(f"{len(missing)} edge endpoints are not in the node table, e.g. {sorted(missing)[:3]}")

    for record in nodes.to_dict("records"):
        sid = record["string_id"]
        if not keep_isolated and sid not in connected:
            continue
        attrs = {k: v for k, v in ((k, _clean(v)) for k, v in record.items()) if v is not None}
        attrs["is_isolated"] = sid not in connected
        graph.add_node(sid, **attrs)

    for record in edges.to_dict("records"):
        attrs = {k: v for k, v in ((k, _clean(v)) for k, v in record.items())
                 if v is not None and k not in ("node_u", "node_v")}
        attrs["weight"] = attrs["combined_score"]
        graph.add_edge(record["node_u"], record["node_v"], **attrs)
    return graph


def summary_statistics(graph: nx.Graph) -> dict[str, Any]:
    """Descriptive topology only. No biological interpretation is attached to these numbers."""
    n, m = graph.number_of_nodes(), graph.number_of_edges()
    degrees = [d for _, d in graph.degree()]
    strengths = [s for _, s in graph.degree(weight="weight")]
    components = sorted((len(c) for c in nx.connected_components(graph)), reverse=True) if n else []
    weights = [d["weight"] for _, _, d in graph.edges(data=True)]
    stats: dict[str, Any] = {
        "nodes": n,
        "edges": m,
        "isolated_nodes": sum(1 for d in degrees if d == 0),
        "density": round(nx.density(graph), 6) if n > 1 else 0.0,
        "connected_components": len(components),
        "largest_component_nodes": components[0] if components else 0,
        "largest_component_fraction": round(components[0] / n, 4) if n else None,
        "degree_mean": round(sum(degrees) / n, 4) if n else None,
        "degree_max": max(degrees) if degrees else None,
        "weighted_degree_mean": round(sum(strengths) / n, 4) if n else None,
        "edge_weight_min": min(weights) if weights else None,
        "edge_weight_median": float(pd.Series(weights).median()) if weights else None,
        "edge_weight_max": max(weights) if weights else None,
        "average_clustering_unweighted": round(nx.average_clustering(graph), 4) if n else None,
        "edges_in_physical_subnetwork": sum(1 for *_, d in graph.edges(data=True) if d.get("in_physical_subnetwork")),
        "edges_with_experiments_channel_score": sum(1 for *_, d in graph.edges(data=True)
                                                if d.get("has_experiments_channel_score")),
    }
    return stats


def degree_table(graph: nx.Graph) -> pd.DataFrame:
    rows = [{
        "string_id": node,
        "string_preferred_name": data.get("string_preferred_name"),
        "degree": graph.degree(node),
        "weighted_degree": round(graph.degree(node, weight="weight"), 4),
        "is_seed": data.get("is_seed"),
    } for node, data in graph.nodes(data=True)]
    frame = pd.DataFrame(rows, columns=["string_id", "string_preferred_name", "degree", "weighted_degree", "is_seed"])
    return frame.sort_values(["degree", "weighted_degree", "string_id"], ascending=[False, False, True]).reset_index(drop=True)
