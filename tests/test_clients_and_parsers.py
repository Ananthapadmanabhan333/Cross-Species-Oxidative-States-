import json

import pytest

from extremonet.data_collection.http_cache import CachedHttpClient, DataSourceError
from extremonet.data_collection.string_client import StringClient
from extremonet.data_collection.uniprot_client import UniProtClient, build_seed_query
from extremonet.preprocessing.parsers import (
    parse_string_mapping, parse_string_network, parse_uniprot_seeds,
)

from .conftest import DR_SEEDS, UNIPROT_HEADER, FakeResponse, FakeSession


def http(tmp_path, session):
    return CachedHttpClient(tmp_path / "raw", min_interval_seconds=0, backoff_seconds=0, max_retries=0,
                            session=session, sleep=lambda s: None)


def test_seed_query_combines_go_terms_taxon_and_review_status():
    q = build_seed_query(["GO:0006979", "GO:0016209"], "organism_id:9606", reviewed_only=True)
    assert q == "((go:0006979) OR (go:0016209)) AND (organism_id:9606) AND (reviewed:true)"
    assert "reviewed" not in build_seed_query(["GO:0006979"], "taxonomy_id:1299", reviewed_only=False)


def test_uniprot_seed_parsing_keeps_every_row(fake_session, tmp_path):
    raw = UniProtClient(http(tmp_path, fake_session)).fetch_seed_proteins(
        ["GO:0006979"], "taxonomy_id:1299", False, "dr")
    seeds = parse_uniprot_seeds(raw.text)
    assert len(seeds) == len(DR_SEEDS)
    assert seeds.loc[0, "uniprot_xref_string"] == "243230.DR_T001"
    assert seeds["uniprot_reviewed"].dtype == bool
    assert "query" in raw.meta["params"]


def test_uniprot_parser_rejects_unexpected_columns():
    with pytest.raises(ValueError, match="columns"):
        parse_uniprot_seeds("Entry\tGene\nP00001\tx\n")


def test_uniprot_entry_count_reads_header(fake_session, tmp_path):
    count, raw = UniProtClient(http(tmp_path, fake_session)).fetch_entry_count("taxonomy_id:1299", False, "dr")
    assert count == 3104
    assert raw.meta["response_headers"]["x-total-results"] == "3104"


def test_uniprot_entry_count_requires_header(tmp_path):
    session = FakeSession(lambda m, u, p: FakeResponse(body="Entry\n"))
    with pytest.raises(ValueError, match="X-Total-Results"):
        UniProtClient(http(tmp_path, session)).fetch_entry_count("organism_id:9606", True, "hs")


def test_string_mapping_sends_identity_and_sorted_identifiers(fake_session, tmp_path):
    client = StringClient(http(tmp_path, fake_session), caller_identity="unit_test")
    rows, raw = client.map_identifiers(["P00002", "P00001", "P00001"], 1299, "dr")
    method, url, params = fake_session.calls[0]
    assert method == "POST" and url.endswith("/json/get_string_ids")
    assert params["identifiers"] == "P00001\rP00002"
    assert params["caller_identity"] == "unit_test"
    assert params["species"] == 1299
    mapping = parse_string_mapping(rows)
    assert set(mapping["string_taxon_id"]) == {243230}


def test_string_network_parsing_preserves_channel_scores(fake_session, tmp_path):
    client = StringClient(http(tmp_path, fake_session))
    rows, _ = client.network(["243230.DR_T001", "243230.DR_T002"], 243230, 400, "functional", 0, "dr")
    edges = parse_string_network(rows, "functional")
    assert {"experimental_score", "textmining_score", "combined_score"} <= set(edges.columns)
    assert edges["experimental_score"].dtype == float
    assert list(edges["source_record"]) == list(range(len(edges)))
    assert (edges["network_type"] == "functional").all()


def test_string_error_object_is_raised(tmp_path):
    session = FakeSession(lambda m, u, p: FakeResponse(body=json.dumps({"Error": "species not found"})))
    with pytest.raises(DataSourceError, match="species not found"):
        StringClient(http(tmp_path, session)).version()


def test_string_non_json_is_raised(tmp_path):
    session = FakeSession(lambda m, u, p: FakeResponse(body="<html>maintenance</html>"))
    with pytest.raises(DataSourceError, match="non-JSON"):
        StringClient(http(tmp_path, session)).version()


def test_empty_uniprot_result_gives_empty_frame():
    assert parse_uniprot_seeds("").empty
    assert parse_uniprot_seeds(UNIPROT_HEADER + "\n").empty
