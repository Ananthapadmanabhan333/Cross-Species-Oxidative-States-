import networkx as nx
import pandas as pd
import pytest

from extremonet.annotation.annotate import annotation_coverage, build_node_table, matched_terms
from extremonet.network_construction.graph_builder import (
    annotate_edges, build_graph, degree_table, summary_statistics,
)
from extremonet.preprocessing.parsers import parse_string_network
from extremonet.quality_control.validators import validate_edges

from .conftest import NETWORKS

SPECIES = {"label": "Deinococcus radiodurans", "ncbi_taxon_id": 1299, "string_taxon_id": 243230}
NODE_IDS = [f"243230.DR_T00{i}" for i in range(1, 5)]


def nodes_frame():
    mapped = pd.DataFrame({
        "string_id": NODE_IDS,
        "uniprot_accession": ["P00001", "P00002", "Q00003", "Q00004"],
        "string_preferred_name": ["tstA", "tstB", "DR_T003", "tstD"],
        "string_annotation": ["a", "b", "annotation not available", "d"],
        "gene_primary": ["tstA", "tstB", "", "tstD"],
        "protein_name": ["A", "B", "C", "D"],
        "uniprot_reviewed": [True, True, False, False],
        "go_ids": ["GO:0006979; GO:0005737", "GO:0016209", "GO:0034599", "GO:0006979"],
        "xref_agrees": pd.array([True, True, True, pd.NA], dtype="boolean"),
    })
    return build_node_table(mapped, "deinococcus_radiodurans", SPECIES, ["GO:0006979", "GO:0016209"])


def edges_frames():
    out = {}
    for kind in ("functional", "physical"):
        parsed = parse_string_network(NETWORKS[(243230, kind)], kind)
        out[kind], _ = validate_edges(parsed, set(NODE_IDS), 243230, 0.7)
    return annotate_edges(out["functional"], out["physical"], 0.0, "STRING v12.5 functional association")


def test_matched_terms_lists_only_verbatim_seed_terms():
    assert matched_terms("GO:0006979; GO:0005737", ["GO:0006979", "GO:0016209"]) == "GO:0006979"
    assert matched_terms("GO:0034599", ["GO:0006979"]) == ""
    assert matched_terms(None, ["GO:0006979"]) == ""


def test_edges_carry_physical_and_experimental_flags():
    edges, orphans = edges_frames()
    assert orphans.empty
    flags = edges.set_index(["node_u", "node_v"])
    first = flags.loc[("243230.DR_T001", "243230.DR_T002")]
    assert bool(first["in_physical_subnetwork"]) and first["physical_score"] == 0.80
    assert bool(first["has_experiments_channel_score"])
    second = flags.loc[("243230.DR_T002", "243230.DR_T003")]
    assert not bool(second["in_physical_subnetwork"]) and not bool(second["has_experiments_channel_score"])


def test_physical_edge_without_functional_counterpart_is_reported():
    functional = parse_string_network(NETWORKS[(243230, "functional")][:1], "functional")
    physical = parse_string_network(NETWORKS[(243230, "functional")][1:2], "physical")
    f, _ = validate_edges(functional, set(NODE_IDS), 243230, 0.7)
    p, _ = validate_edges(physical, set(NODE_IDS), 243230, 0.7)
    _, orphans = annotate_edges(f, p, 0.0, "x")
    assert len(orphans) == 1


def test_graph_is_weighted_and_keeps_isolated_seeds():
    edges, _ = edges_frames()
    graph = build_graph(nodes_frame(), edges, keep_isolated=True, graph_attrs={"species": "D. radiodurans"})
    assert graph.number_of_nodes() == 4 and graph.number_of_edges() == 2
    assert graph["243230.DR_T001"]["243230.DR_T002"]["weight"] == 0.95
    assert graph.nodes["243230.DR_T004"]["is_isolated"] is True
    assert graph.nodes["243230.DR_T001"]["uniprot_accession"] == "P00001"
    assert graph.graph["species"] == "D. radiodurans"
    without = build_graph(nodes_frame(), edges, keep_isolated=False)
    assert "243230.DR_T004" not in without


def test_graph_rejects_edges_to_unknown_nodes():
    edges, _ = edges_frames()
    with pytest.raises(ValueError, match="not in the node table"):
        build_graph(nodes_frame().iloc[:1], edges)


def test_graphml_round_trip(tmp_path):
    edges, _ = edges_frames()
    graph = build_graph(nodes_frame(), edges)
    path = tmp_path / "g.graphml"
    nx.write_graphml(graph, path)
    loaded = nx.read_graphml(path)
    assert loaded.number_of_edges() == 2
    assert float(loaded["243230.DR_T001"]["243230.DR_T002"]["weight"]) == 0.95


def test_summary_statistics_and_degree_table():
    edges, _ = edges_frames()
    graph = build_graph(nodes_frame(), edges)
    stats = summary_statistics(graph)
    assert stats["nodes"] == 4 and stats["edges"] == 2 and stats["isolated_nodes"] == 1
    assert stats["connected_components"] == 2 and stats["largest_component_nodes"] == 3
    assert stats["edges_in_physical_subnetwork"] == 1 and stats["edges_with_experiments_channel_score"] == 1
    assert stats["edge_weight_max"] == 0.95
    table = degree_table(graph)
    assert table.iloc[0]["string_id"] == "243230.DR_T002" and table.iloc[0]["degree"] == 2


def test_annotation_coverage_counts_missing_values():
    cov = annotation_coverage(nodes_frame())
    assert cov["nodes"] == 4
    assert cov["with_gene_name"] == {"count": 3, "fraction": 0.75}
    assert cov["with_string_annotation"]["count"] == 3        # "annotation not available" excluded
    assert cov["with_seed_go_term_listed_verbatim"]["count"] == 3
    assert cov["uniprot_reviewed"]["count"] == 2
