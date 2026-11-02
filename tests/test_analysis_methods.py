"""Unit tests for Phase 2 methods on small, hand-checkable SYNTHETIC graphs and annotation sets."""

import math

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

from extremonet.analysis import alignment, communities, embedding, enrichment, ranking, similarity, topology
from extremonet.data_collection.evidence_sources import evidence_class


def weighted(edges):
    g = nx.Graph()
    for u, v, w in edges:
        g.add_edge(u, v, weight=w)
    return g


def two_cliques():
    g = nx.Graph()
    for block in (["a1", "a2", "a3", "a4"], ["b1", "b2", "b3", "b4"]):
        for i, u in enumerate(block):
            for v in block[i + 1:]:
                g.add_edge(u, v, weight=0.9)
    g.add_edge("a1", "b1", weight=0.7)
    g.add_node("iso")
    return g


# ------------------------------------------------------------------------------------ topology
def test_node_metrics_star_hub_and_isolate():
    g = weighted([("hub", f"leaf{i}", 0.8) for i in range(5)])
    g.add_node("alone")
    m = topology.node_metrics(g).set_index("string_id")
    assert m.loc["hub", "degree"] == 5 and bool(m.loc["hub", "is_hub"])
    assert m.loc["hub", "hub_bottleneck_class"] == "hub-bottleneck"
    assert m.loc["alone", "hub_bottleneck_class"] == "isolated"
    assert math.isnan(m.loc["alone", "degree_pct"])         # percentiles only over connected nodes
    assert m.loc["leaf0", "betweenness"] == 0.0


def test_graph_statistics_handles_disconnected_graph():
    s = topology.graph_statistics(two_cliques())
    assert s["connected_components"] == 2 and s["isolated_nodes"] == 1
    assert s["largest_component_nodes"] == 8 and s["lcc_diameter_hops"] == 3


def test_module_roles_participation():
    g = two_cliques()
    membership = {n: (0 if n.startswith("a") else 1 if n.startswith("b") else 2) for n in g}
    roles = topology.module_roles(g, membership).set_index("string_id")
    assert roles.loc["a2", "participation_coefficient"] == 0.0     # all links inside its module
    assert roles.loc["a1", "participation_coefficient"] == pytest.approx(1 - (3 / 4) ** 2 - (1 / 4) ** 2)
    assert math.isnan(roles.loc["iso", "participation_coefficient"])


def test_size_matched_null_is_seeded():
    g = two_cliques()
    a = topology.size_matched_null(g, 4, 20, seed=7)
    b = topology.size_matched_null(g, 4, 20, seed=7)
    pd.testing.assert_frame_equal(a, b)
    pos = topology.empirical_position(a["edges"], 6.0)
    assert 0 <= pos["fraction_null_below"] <= 1


# --------------------------------------------------------------------------------- communities
def test_ari_known_values():
    assert communities.adjusted_rand_index([0, 0, 1, 1], [5, 5, 9, 9]) == 1.0
    # Hubert & Arabie example: two 3-element groups vs one merged split
    assert communities.adjusted_rand_index([0, 0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2]) == pytest.approx(0.2424, abs=1e-3)


def test_louvain_recovers_cliques_and_is_stable():
    result = communities.stability(two_cliques(), 1.0, seeds=[1, 2, 3, 4])
    part = result["partition"]
    assert part["a1"] == part["a4"] != part["b1"] == part["b4"]
    assert part["iso"] not in {part["a1"], part["b1"]}
    assert result["summary"]["ari_mean"] == 1.0
    ns = communities.node_stability(result).set_index("string_id")
    assert ns.loc["a2", "coassignment_with_module"] == 1.0 and math.isnan(ns.loc["iso", "coassignment_with_module"])


# ---------------------------------------------------------------------------------- enrichment
def test_benjamini_hochberg_matches_definition():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    expected = np.array([0.04, 0.05333333, 0.05333333, 0.2])
    np.testing.assert_allclose(enrichment.benjamini_hochberg(p), expected, rtol=1e-6)


def test_module_enrichment_hypergeometric():
    membership = {f"n{i}": (0 if i < 4 else 1) for i in range(10)}
    ann = pd.DataFrame({"string_id": ["n0", "n1", "n2", "n5"], "category": "Process",
                        "term": "GO:1", "description": "test term"})
    table, summary = enrichment.module_enrichment(membership, ann, ["Process"], min_module_size=3, min_term_size=2)
    row = table[table["module"] == 0].iloc[0]
    assert (row.k_in_module, row.K_in_background, row.N_background) == (3, 4, 10)
    assert row.p_value == pytest.approx(hypergeom.sf(2, 10, 4, 4))
    assert summary["modules_tested"] == 2 and "fdr_bh" in table


# ---------------------------------------------------------------------------------- similarity
def test_similarity_measures():
    ic = {"x": 1.0, "y": 2.0, "z": 3.0}
    a, b = {"x", "y"}, {"y", "z"}
    assert similarity.simgic(a, b, ic, "simgic") == pytest.approx(2 / 6)
    assert similarity.simgic(a, b, ic, "ochiai") == pytest.approx(2 / math.sqrt(3 * 5))
    assert similarity.simgic(a, b, ic, "min") == pytest.approx(2 / 3)
    assert similarity.simgic(set(), b, ic) is None            # missing, not dissimilar
    with pytest.raises(ValueError):
        similarity.simgic(a, b, ic, "bogus")


def test_information_content():
    ic = similarity.information_content([{"a"}, {"a", "b"}, set()])
    assert ic["a"] == 0.0 and ic["b"] == pytest.approx(math.log(2))


def test_ec_similarity_levels():
    assert similarity.ec_similarity({"1.15.1.1"}, {"1.15.1.1", "1.8.-.-"}) == (1.0, ["1.15.1.1"])
    assert similarity.ec_similarity({"1.11.1.6"}, {"1.11.1.24"}) == (0.5, ["1.11.1.-"])
    assert similarity.ec_similarity({"1.11.1.6"}, {"2.7.1.1"})[0] == 0.0
    assert similarity.ec_similarity(set(), {"1.1.1.1"})[0] is None
    assert similarity.ec_similarity({"1.-.-.-"}, {"1.-.-.-"})[0] == 0.0   # incomplete ECs never "shared"


def test_kegg_map_and_evidence_class():
    assert similarity.kegg_map("dra00380") == similarity.kegg_map("hsa00380") == "00380"
    assert similarity.kegg_map("GO:0001") is None
    assert evidence_class("IMP") == "experimental" and evidence_class("IEA") == "electronic"


def prot(sid, **kw):
    return similarity.ProteinEvidence(sid, "x", None, sid, **kw)


def test_homology_tiers():
    tables = similarity.HomologyTables({"d1": ("h1", 200.0), "d2": ("h9", 80.0)}, {"h1": ("d1", 190.0)})
    assert similarity.homology(prot("d1"), prot("h1"), tables, 50)[1] == "reciprocal_best_hit"
    assert similarity.homology(prot("d2", panther={"P1"}), prot("h2", panther={"P1"}), tables, 50)[1] == "shared_panther_family"
    assert similarity.homology(prot("d2"), prot("h9"), tables, 50)[1] == "one_way_best_hit"
    assert similarity.homology(prot("d3", pfam={"PF1"}), prot("h3", pfam={"PF1"}), tables, 50)[1] == "shared_pfam_domain"
    assert similarity.homology(prot("d3", interpro={"IPR1"}), prot("h3", interpro={"IPR1"}), tables, 50)[1] \
        == "shared_interpro_entry"
    assert similarity.homology(prot("d3", ko={"K1"}, pfam={"PF1"}), prot("h3", ko={"K1"}, pfam={"PF1"}),
                               tables, 50)[1] == "shared_kegg_ko"
    _, _, detail = similarity.homology(prot("d3", ko={"K1"}), prot("h3", ko={"K2"}), tables, 50)
    assert detail["ko_conflict"] is True
    score, tier, _ = similarity.homology(prot("d3"), prot("h3"), tables, 50)
    assert (score, tier) == (0.0, "none_detected")
    assert similarity.homology(prot("d2"), prot("h9"), tables, 100)[1] == "none_detected"   # below min bitscore


def test_frequent_signatures_become_weak_common_fold_tier():
    common = {"IPR_PLOOP"}
    dr = {f"d{i}": prot(f"d{i}", interpro=common | ({"IPR_RARE"} if i == 0 else set()), pfam={"PF_X"}) for i in range(5)}
    hs = {f"h{i}": prot(f"h{i}", interpro=common | ({"IPR_RARE"} if i == 0 else set())) for i in range(15)}
    result = similarity.separate_frequent_signatures({"dr": dr, "hs": hs}, max_fraction=0.2)
    assert set(result["common_fold_signatures"]["interpro"]) == {"IPR_PLOOP"}   # 20/20 carry it
    assert set(result["common_fold_signatures"]["pfam"]) == {"PF_X"}           # 5/20 = 25 % > 20 %
    assert dr["d0"].interpro == {"IPR_RARE"} and dr["d0"].common_folds == {"IPR_PLOOP", "PF_X"}
    tables = similarity.HomologyTables({}, {})
    assert similarity.homology(dr["d0"], hs["h0"], tables, 50)[1] == "shared_interpro_entry"
    score, tier, _ = similarity.homology(dr["d1"], hs["h1"], tables, 50)
    assert (score, tier) == (0.2, "shared_common_fold")                         # still homology, weakly
    assert ranking.assign_category(base_row(homology_tier="shared_common_fold", ec=1.0, pathway=0.9), 0.3, 0.3) \
        == "distant_homolog"                                                     # never an analogue


def test_role_similarity_not_assessable_for_isolates():
    a = pd.Series({"is_isolated": False, "degree_pct": 0.9, "betweenness_pct": 0.8})
    b = pd.Series({"is_isolated": False, "degree_pct": 0.7, "betweenness_pct": 0.8})
    assert similarity.role_similarity(a, b, ["degree_pct", "betweenness_pct"]) == pytest.approx(0.9)
    iso = pd.Series({"is_isolated": True, "degree_pct": np.nan, "betweenness_pct": np.nan})
    assert similarity.role_similarity(a, iso, ["degree_pct"]) is None


# ------------------------------------------------------------------------------------- ranking
def test_missing_evidence_is_not_negative():
    w = np.array([0.5, 0.5])
    values = np.array([[0.8, np.nan], [0.8, 0.0], [np.nan, np.nan]])
    score, coverage = ranking.score_matrix(values, w)
    assert score[0] == pytest.approx(0.8) and coverage[0] == 0.5     # missing component ignored
    assert score[1] == pytest.approx(0.4)                             # observed 0 does count
    assert np.isnan(score[2]) and coverage[2] == 0


def base_row(**kw):
    row = {"homology_tier": "none_detected", "ec": None, "go": None, "pathway": None, "literature": None,
           "role": None, "module": None}
    row.update(kw)
    return row


def test_category_rules():
    assert ranking.assign_category(base_row(homology_tier="reciprocal_best_hit"), 0.3, 0.3) == "ortholog_rbh"
    assert ranking.assign_category(base_row(homology_tier="shared_interpro_entry", ec=1.0), 0.3, 0.3) == "homolog"
    assert ranking.assign_category(base_row(homology_tier="shared_kegg_ko"), 0.3, 0.3) == "ortholog_ko"
    assert ranking.assign_category(base_row(ec=1.0, pathway=0.5), 0.3, 0.3) == "candidate_analogue"
    assert ranking.assign_category(base_row(ec=1.0, literature=1.0), 0.3, 0.3) == "candidate_analogue"
    assert ranking.assign_category(base_row(ec=1.0), 0.3, 0.3) == "shared_function_annotation"
    assert ranking.assign_category(base_row(pathway=0.4), 0.3, 0.3) == "similar_pathway_role"
    assert ranking.assign_category(base_row(role=0.9), 0.3, 0.3) == "speculative"
    assert ranking.assign_category(base_row(go=float("nan"), role=float("nan")), 0.3, 0.3) == "unsupported"


def synthetic_pairs():
    rng = np.random.default_rng(0)
    rows = []
    for i in range(30):
        vals = {c: (float(rng.random()) if rng.random() > 0.3 else None) for c in ranking.COMPONENTS}
        vals["homology"] = 0.0
        rows.append({"dr_string_id": f"d{i % 5}", "dr_name": f"d{i % 5}", "dr_uniprot": None,
                     "hs_string_id": f"h{i}", "hs_name": f"h{i}", "hs_uniprot": None, **vals,
                     "homology_tier": "none_detected", "shared_ec": "", "dr_ec": "", "hs_ec": "",
                     "shared_panther": "", "shared_pfam": "", "shared_interpro": "", "shared_pathways": "",
                     "dr_best_human_hit": None, "dr_best_hit_bitscore": None, "dr_isolated": False,
                     "hs_isolated": False, "dr_evidence_codes": "IEA", "hs_evidence_codes": "IDA",
                     "references": "", "negated_go_conflicts": ""})
    return pd.DataFrame(rows)


def test_tied_scores_are_ordered_by_identifiers():
    pairs = synthetic_pairs().head(4).copy()
    for c in ranking.COMPONENTS:
        pairs[c] = 0.5                                   # every pair scores exactly the same
    pairs = pairs.iloc[::-1].reset_index(drop=True)       # input order must not matter
    weights = {c: 1 / len(ranking.COMPONENTS) for c in ranking.COMPONENTS}
    ranked = ranking.rank_pairs(pairs, weights, 0.3, 0.3, 0.3)
    keys = list(zip(ranked["dr_string_id"], ranked["hs_string_id"]))
    assert keys == sorted(keys)


def test_set_similarity_is_independent_of_set_iteration_order():
    ic = {f"t{i}": 0.1 * (i + 1) ** 1.7 for i in range(40)}
    a = {f"t{i}" for i in range(0, 40, 2)}
    b = {f"t{i}" for i in range(0, 40, 3)}
    reference = similarity.simgic(a, b, ic, "ochiai")
    for perm in (sorted(a, reverse=True), sorted(a, key=hash)):
        assert similarity.simgic(set(perm), set(b), ic, "ochiai") == reference


def test_rank_pairs_and_weight_sensitivity_are_deterministic():
    pairs = synthetic_pairs()
    weights = {c: 1 / len(ranking.COMPONENTS) for c in ranking.COMPONENTS}
    ranked = ranking.rank_pairs(pairs, weights, 0.3, 0.3, 0.3)
    assert ranked["rank"].dropna().is_monotonic_increasing
    assert (ranked.loc[ranked["eligible"], "score"].diff().dropna() <= 1e-12).all()
    assert {"supporting_evidence", "contradictory_evidence", "uncertainty", "proposed_validation"} <= set(ranked)
    a = ranking.weight_sensitivity(pairs, ranked, weights, 20, 50, 5, seed=3)
    b = ranking.weight_sensitivity(pairs, ranked, weights, 20, 50, 5, seed=3)
    assert a["dirichlet"] == b["dirichlet"]
    assert len(a["leave_one_out"]) == len(ranking.COMPONENTS)


# ----------------------------------------------------------------------------------- alignment
def test_isorank_with_prior_recovers_identity_on_isomorphic_graphs():
    g1 = weighted([("a", "b", 0.9), ("b", "c", 0.8), ("c", "d", 0.7)])
    g2 = nx.relabel_nodes(g1, {"a": "A", "b": "B", "c": "C", "d": "D"})
    n1, n2 = sorted(g1), sorted(g2)
    prior = np.eye(4) + 0.01
    r, _, _, info = alignment.isorank(g1, g2, prior, alpha=0.5, max_iter=200, tol=1e-12)
    mapping = {a: b for a, b, _ in alignment.greedy_matching(r, n1, n2)}
    assert mapping == {"a": "A", "b": "B", "c": "C", "d": "D"} and info["converged"]
    assert alignment.edge_conservation(g1, g2, mapping) == (3, 3)


def test_isorank_alpha_below_one_requires_prior():
    g = weighted([("a", "b", 0.9)])
    with pytest.raises(ValueError):
        alignment.isorank(g, g, None, alpha=0.5)


def test_rewired_null_keeps_degree_sequence():
    g1 = two_cliques()
    g2 = nx.relabel_nodes(two_cliques(), lambda n: n.upper())
    values, distinct = alignment.rewired_null(g1, g2, None, 1.0, 5, seed=1, max_iter=50, tol=1e-9)
    assert len(values) == 5 and (values <= g1.number_of_edges()).all() and 0 <= distinct <= 5


# ----------------------------------------------------------------------------------- embedding
def test_netmf_shape_and_determinism():
    g = two_cliques()
    e1, nodes = embedding.netmf_embedding(g, 4, 3, 1.0)
    e2, _ = embedding.netmf_embedding(g, 4, 3, 1.0)
    assert e1.shape == (g.number_of_nodes(), 4) and nodes == sorted(g)
    np.testing.assert_allclose(e1, e2)


def test_link_prediction_never_sees_held_out_edges(monkeypatch):
    g = nx.relabel_nodes(nx.connected_caveman_graph(6, 6), str)
    nx.set_edge_attributes(g, 0.9, "weight")
    seen = []
    original = embedding.netmf_embedding

    def spy(graph, *args):
        seen.append({tuple(sorted(e)) for e in graph.edges})
        return original(graph, *args)
    monkeypatch.setattr(embedding, "netmf_embedding", spy)
    result = embedding.link_prediction(g, 8, 3, 1.0, 0.2, 2, seed=5)
    full = {tuple(sorted(e)) for e in g.edges}
    assert all(len(train) < len(full) for train in seen)     # held-out edges removed before embedding
    assert set(result["method"]) >= {"netmf_cosine", "adamic_adar"}
    assert (result["auc"].between(0, 1)).all()


def test_auc():
    assert embedding._auc(np.array([3.0, 4.0]), np.array([1.0, 2.0])) == 1.0
    assert embedding._auc(np.array([1.0]), np.array([1.0])) == 0.5
