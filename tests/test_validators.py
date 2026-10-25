import pandas as pd

from extremonet.preprocessing.parsers import parse_string_mapping, parse_string_network, parse_uniprot_seeds
from extremonet.quality_control.validators import (
    is_uniprot_accession, reason_counts, string_id_taxon, validate_edges, validate_mapping, validate_seeds,
)

from .conftest import DR_SEEDS, HS_SEEDS, NETWORKS, UNIPROT_HEADER, edge


def seeds_frame(rows):
    return parse_uniprot_seeds("\n".join([UNIPROT_HEADER, *rows]) + "\n")


def test_identifier_helpers():
    assert is_uniprot_accession("P04040") and is_uniprot_accession("A0A023GPI8")
    assert not is_uniprot_accession("not-an-acc") and not is_uniprot_accession(None)
    assert string_id_taxon("243230.DR_1998") == 243230
    assert string_id_taxon("ENSP000001") is None


def test_seed_validation_reasons():
    kept, dropped = validate_seeds(seeds_frame(DR_SEEDS), {1299, 243230})
    assert reason_counts(dropped) == {"invalid_uniprot_accession": 1, "unexpected_organism": 1}
    assert len(kept) + len(dropped) == len(DR_SEEDS)
    kept, dropped = validate_seeds(seeds_frame(HS_SEEDS), {9606})
    assert reason_counts(dropped) == {"duplicate_accession": 1}


def mapping_rows(pairs):
    return parse_string_mapping([
        {"queryItem": acc, "queryIndex": i, "stringId": sid, "ncbiTaxonId": tax,
         "preferredName": acc, "annotation": "x"} for i, (acc, sid, tax) in enumerate(pairs)])


def test_mapping_validation_flags_unmapped_wrong_taxon_and_collapse():
    seeds, _ = validate_seeds(seeds_frame(HS_SEEDS[:5]), {9606})
    mapping = mapping_rows([
        ("P10001", "9606.ENSPTEST01", 9606),
        ("P10002", "9606.ENSPTEST01", 9606),   # same STRING protein as P10001
        ("P10003", "10090.ENSMUSP1", 10090),   # wrong species
    ])                                          # P10004 unmapped
    kept, dropped = validate_mapping(seeds, mapping, expected_string_taxon=9606)
    assert reason_counts(dropped) == {"collapsed_to_same_string_id": 1, "string_taxon_mismatch": 1,
                                      "unmapped_in_string": 1}
    assert list(kept["string_id"]) == ["9606.ENSPTEST01"]
    assert kept.iloc[0]["uniprot_accessions"] == "P10001;P10002"


def test_mapping_records_xref_disagreement_without_excluding():
    seeds, _ = validate_seeds(seeds_frame(HS_SEEDS[:3]), {9606})
    mapping = mapping_rows([("P10001", "9606.ENSPTEST01", 9606), ("P10002", "9606.ENSPTEST02", 9606),
                            ("P10003", "9606.ENSPTEST03", 9606)])
    kept, dropped = validate_mapping(seeds, mapping, 9606)
    assert dropped.empty
    assert kept.set_index("uniprot_accession")["xref_agrees"].to_dict() == {
        "P10001": True, "P10002": True, "P10003": False}


def test_edge_validation_reasons_and_duplicate_keeps_highest_score():
    edges = parse_string_network(NETWORKS[(243230, "functional")], "functional")
    nodes = {f"243230.DR_T00{i}" for i in range(1, 5)}
    kept, dropped = validate_edges(edges, nodes, 243230, score_threshold=0.7)
    assert reason_counts(dropped) == {"below_threshold": 1, "duplicate_edge": 1, "self_loop": 1}
    pair = kept[(kept.node_u == "243230.DR_T001") & (kept.node_v == "243230.DR_T002")]
    assert pair["combined_score"].tolist() == [0.95]
    assert len(kept) + len(dropped) == len(edges)


def test_edge_validation_species_and_score_checks():
    rows = [
        edge("9606.A", "9606.B", None, 9606),            # missing_score
        edge("9606.A", "10090.C", 0.9, 9606),            # species_mismatch (id prefix)
        edge("9606.A", "9606.B", 0.9, 10090),            # species_mismatch (reported taxon)
        edge("9606.A", "9606.B", 1.7, 9606),             # score_out_of_range
        edge("9606.A", "9606.Z", 0.9, 9606),             # endpoint_not_in_node_set
        edge("bad", "9606.B", 0.9, 9606),                # invalid_string_id
        edge("9606.A", "9606.B", 0.8, 9606),             # retained
    ]
    edges = parse_string_network(rows, "functional")
    kept, dropped = validate_edges(edges, {"9606.A", "9606.B"}, 9606, 0.7)
    assert reason_counts(dropped) == {"endpoint_not_in_node_set": 1, "invalid_string_id": 1, "missing_score": 1,
                                      "score_out_of_range": 1, "species_mismatch": 2}
    assert len(kept) == 1
    kept, dropped = validate_edges(edges, {"9606.A", "9606.B"}, 9606, 0.7, allow_new_nodes=True)
    assert "endpoint_not_in_node_set" not in reason_counts(dropped)


def test_threshold_is_inclusive():
    edges = parse_string_network([edge("9606.A", "9606.B", 0.7, 9606)], "functional")
    kept, _ = validate_edges(edges, {"9606.A", "9606.B"}, 9606, 0.7)
    assert len(kept) == 1


def test_empty_edge_table():
    edges = parse_string_network([], "physical")
    kept, dropped = validate_edges(edges, {"9606.A"}, 9606, 0.7)
    assert kept.empty and dropped.empty and isinstance(kept, pd.DataFrame)
