import pytest

from extremonet.config import ConfigError, config_fingerprint, load_config, validate_config

from .conftest import REPO_ROOT

CONFIG = REPO_ROOT / "config" / "phase1_oxidative_stress.yaml"


def test_repository_config_is_valid():
    cfg = load_config(CONFIG)
    assert cfg["network"]["edge_score_threshold"] == 0.7
    assert cfg["species"]["deinococcus_radiodurans"]["ncbi_taxon_id"] == 1299
    assert cfg["species"]["homo_sapiens"]["ncbi_taxon_id"] == 9606


def test_overrides_are_applied_and_validated():
    cfg = load_config(CONFIG, {"mode": "offline", "network.edge_score_threshold": 0.9})
    assert cfg["mode"] == "offline"
    assert cfg["network"]["edge_score_threshold"] == 0.9
    with pytest.raises(ConfigError, match="mode"):
        load_config(CONFIG, {"mode": "sometimes"})


def test_threshold_below_retrieval_score_is_rejected():
    with pytest.raises(ConfigError, match="below the retrieval threshold"):
        load_config(CONFIG, {"network.edge_score_threshold": 0.3})


def test_threshold_out_of_range_is_rejected():
    with pytest.raises(ConfigError):
        load_config(CONFIG, {"network.edge_score_threshold": 7})


def test_species_must_have_distinct_string_taxa():
    cfg = load_config(CONFIG)
    cfg["species"]["homo_sapiens"]["string_taxon_id"] = 243230
    with pytest.raises(ConfigError, match="used twice"):
        validate_config(cfg)


def test_go_terms_must_be_go_identifiers():
    with pytest.raises(ConfigError, match="go_terms"):
        load_config(CONFIG, {"seeds.go_terms": ["oxidative stress"]})


def test_fingerprint_changes_with_analysis_parameters():
    a = load_config(CONFIG)
    b = load_config(CONFIG, {"network.edge_score_threshold": 0.8})
    assert config_fingerprint(a) != config_fingerprint(b)
    assert config_fingerprint(a) == config_fingerprint(load_config(CONFIG))
