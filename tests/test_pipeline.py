"""End-to-end pipeline runs against the fake STRING/UniProt APIs."""

import json

import networkx as nx
import pandas as pd
import pytest

from extremonet.data_collection.http_cache import OfflineCacheMiss, sha256_file
from extremonet.pipeline import run_pipeline


@pytest.fixture
def run(config, fake_session, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return run_pipeline(config, session=fake_session, sleep=lambda s: None)


def test_outputs_exist_for_both_species(run):
    results = run["results_dir"]
    for species in ("deinococcus_radiodurans", "homo_sapiens"):
        for name in ("network.graphml", "weighted_edge_list.tsv", "protein_table.tsv", "degree_table.tsv",
                     "summary_statistics.json", "qc_report.json", "qc_report.md"):
            assert (results / species / name).is_file(), name
        for name in ("nodes.tsv", "edges_functional_retained.tsv", "excluded_records.tsv",
                     "uniprot_seed_records_parsed.tsv", "string_functional_records_parsed.tsv"):
            assert (run["processed_dir"] / species / name).is_file(), name
    for name in ("run_report.md", "species_comparison.tsv", "config_used.yaml", "manifest.json",
                 "SHA256SUMS.txt", "logs/pipeline.jsonl"):
        assert (results / name).is_file(), name


def test_species_are_kept_apart(run):
    for species, taxon in (("deinococcus_radiodurans", 243230), ("homo_sapiens", 9606)):
        graph = nx.read_graphml(run["results_dir"] / species / "network.graphml")
        assert all(node.startswith(f"{taxon}.") for node in graph.nodes)
        assert int(graph.graph["string_taxon_id"]) == taxon
        assert run["species"][species]["taxonomy"]["all_retained_match"]
    dr = run["species"]["deinococcus_radiodurans"]["taxonomy"]
    assert dr["ncbi_taxon_id"] == 1299 and dr["string_taxa_returned"] == [243230]


def test_every_record_is_retained_or_reported(run):
    dr = run["species"]["deinococcus_radiodurans"]
    assert dr["excluded"]["seeds"] == {"invalid_uniprot_accession": 1, "unexpected_organism": 1}
    assert dr["excluded"]["functional_edges"] == {"below_threshold": 1, "duplicate_edge": 1, "self_loop": 1}
    assert dr["edges"]["functional_retained"] + sum(dr["excluded"]["functional_edges"].values()) \
        == dr["edges"]["functional_records_retrieved"]
    hs = run["species"]["homo_sapiens"]
    assert hs["excluded"]["seeds"] == {"duplicate_accession": 1}
    assert hs["excluded"]["mapping"] == {"collapsed_to_same_string_id": 1, "unmapped_in_string": 1}
    assert hs["excluded"]["functional_edges"] == {"missing_score": 1}
    excluded = pd.read_csv(run["processed_dir"] / "homo_sapiens" / "excluded_records.tsv", sep="\t")
    assert len(excluded) == 4


def test_batch_unmapped_accessions_are_requeried_individually(run, fake_session):
    singles = [p["identifiers"] for m, u, p in fake_session.calls
               if u.endswith("/get_string_ids") and "\r" not in p["identifiers"] and p["species"] == 9606]
    assert sorted(singles) == ["P10004", "P10005"]   # P10005 shares a protein with P10002 in the batch
    nodes = pd.read_csv(run["processed_dir"] / "homo_sapiens" / "nodes.tsv", sep="\t")
    accessions = nodes.set_index("string_id")["uniprot_accessions"]
    assert accessions["9606.ENSPTEST02"] == "P10002;P10005"


def test_coverage_and_limitations_are_reported(run):
    dr = run["species"]["deinococcus_radiodurans"]
    assert dr["network_coverage"]["isolated_nodes"] == 1
    assert dr["network_coverage"]["uniprot_entries_for_organism"] == 3104
    assert any("243230" in item and "strain" in item for item in dr["limitations"])
    assert any("Neither is equivalent to an experimentally confirmed" in item for item in dr["limitations"])
    hs = run["species"]["homo_sapiens"]
    assert hs["annotation_coverage"]["uniprot_xref_disagrees_with_string_mapping"]["count"] == 1


def test_manifest_records_provenance_and_valid_checksums(run, tmp_path):
    manifest = json.loads(run["manifest"].read_text())
    assert manifest["string_release"]["string_version"] == "12.5"
    assert manifest["config"]["network"]["edge_score_threshold"] == 0.7
    assert {r["source"] for r in manifest["raw_inputs"]} == {"string", "uniprot"}
    for raw in manifest["raw_inputs"]:
        assert raw["retrieved_at_utc"] and raw["params"]
    for output in manifest["outputs"]:
        assert "\\" not in output["path"] and not output["path"].startswith("/")
        assert sha256_file(tmp_path / output["path"]) == output["sha256"]
    sums = (run["results_dir"] / "SHA256SUMS.txt").read_bytes()
    assert b"\r" not in sums and len(sums.splitlines()) == len(manifest["outputs"])


def test_relative_config_paths_give_portable_manifest_paths(config, fake_session, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config["paths"] = {"raw_dir": "data/raw", "processed_dir": "data/processed",
                       "manifests_dir": "data/manifests", "results_dir": "results"}
    result = run_pipeline(config, session=fake_session, sleep=lambda s: None)
    manifest = json.loads(result["manifest"].read_text())
    assert all(o["path"].startswith(("results/", "data/processed/")) for o in manifest["outputs"])


def test_offline_rerun_reproduces_results_without_network(run, config, offline_session, tmp_path):
    config["mode"] = "offline"
    second = run_pipeline(config, session=offline_session, sleep=lambda s: None)
    assert offline_session.calls == []
    for species in ("deinococcus_radiodurans", "homo_sapiens"):
        a = pd.read_csv(run["results_dir"] / species / "weighted_edge_list.tsv", sep="\t")
        b = pd.read_csv(second["results_dir"] / species / "weighted_edge_list.tsv", sep="\t")
        pd.testing.assert_frame_equal(a, b)
    manifest = json.loads(second["manifest"].read_text())
    assert {r["served_from"] for r in manifest["raw_inputs"]} == {"cache"}


def test_offline_threshold_change_reuses_downloaded_records(run, config, offline_session):
    config["mode"] = "offline"
    config["network"]["edge_score_threshold"] = 0.9
    stricter = run_pipeline(config, session=offline_session, sleep=lambda s: None)
    assert stricter["species"]["deinococcus_radiodurans"]["summary_statistics"]["edges"] == 1


def test_offline_without_cache_fails(config, offline_session, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config["mode"] = "offline"
    with pytest.raises(OfflineCacheMiss):
        run_pipeline(config, session=offline_session, sleep=lambda s: None)
