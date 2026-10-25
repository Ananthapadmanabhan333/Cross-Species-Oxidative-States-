"""Loading and validating the YAML pipeline configuration."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml

VALID_MODES = {"online", "refresh", "offline"}
VALID_NETWORK_TYPES = {"functional", "physical"}


class ConfigError(ValueError):
    """Raised when the configuration is missing a value or holds an inconsistent one."""


def load_config(path: str | Path, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    """Read a YAML config, apply dotted-key overrides (e.g. {"mode": "offline"}) and validate it."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    if not isinstance(config, dict):
        raise ConfigError(f"{path} does not contain a YAML mapping")
    config = copy.deepcopy(config)
    for dotted, value in (overrides or {}).items():
        _set_dotted(config, dotted, value)
    validate_config(config)
    config["_source_path"] = str(path)
    return config


def _set_dotted(config: dict[str, Any], dotted: str, value: Any) -> None:
    node = config
    *parents, leaf = dotted.split(".")
    for key in parents:
        node = node.setdefault(key, {})
    node[leaf] = value


def _require(config: dict[str, Any], dotted: str) -> Any:
    node: Any = config
    for key in dotted.split("."):
        if not isinstance(node, dict) or key not in node:
            raise ConfigError(f"missing required config value: {dotted}")
        node = node[key]
    return node


def validate_config(config: dict[str, Any]) -> None:
    for key in (
        "paths.raw_dir", "paths.processed_dir", "paths.manifests_dir", "paths.results_dir",
        "mode", "http.timeout_seconds", "http.max_retries", "string.base_url",
        "string.caller_identity", "string.retrieval_required_score", "string.network_types",
        "network.edge_score_threshold", "seeds.go_terms", "species",
    ):
        _require(config, key)

    if config["mode"] not in VALID_MODES:
        raise ConfigError(f"mode must be one of {sorted(VALID_MODES)}, got {config['mode']!r}")

    retrieval = config["string"]["retrieval_required_score"]
    if not isinstance(retrieval, int) or not 0 <= retrieval <= 1000:
        raise ConfigError("string.retrieval_required_score must be an integer in [0, 1000]")

    threshold = config["network"]["edge_score_threshold"]
    if not isinstance(threshold, (int, float)) or not 0.0 <= threshold <= 1.0:
        raise ConfigError("network.edge_score_threshold must be a number in [0, 1]")
    if threshold < retrieval / 1000:
        raise ConfigError(
            f"network.edge_score_threshold ({threshold}) is below the retrieval threshold "
            f"({retrieval / 1000}); edges between the two were never downloaded. Lower "
            "string.retrieval_required_score and re-fetch instead."
        )

    types = set(config["string"]["network_types"])
    if not types or not types <= VALID_NETWORK_TYPES:
        raise ConfigError(f"string.network_types must be a non-empty subset of {sorted(VALID_NETWORK_TYPES)}")
    if "functional" not in types:
        raise ConfigError("string.network_types must include 'functional' (it defines the graph edges)")

    go_terms = config["seeds"]["go_terms"]
    if not go_terms or not all(isinstance(t, str) and t.startswith("GO:") and t[3:].isdigit() for t in go_terms):
        raise ConfigError("seeds.go_terms must be a non-empty list of GO:nnnnnnn identifiers")

    species = config["species"]
    if not isinstance(species, dict) or not species:
        raise ConfigError("species must be a non-empty mapping")
    string_taxa = set()
    for key, spec in species.items():
        for field in ("label", "ncbi_taxon_id", "string_taxon_id", "uniprot_taxon_query"):
            if field not in spec:
                raise ConfigError(f"species.{key}.{field} is required")
        for field in ("ncbi_taxon_id", "string_taxon_id"):
            if not isinstance(spec[field], int) or spec[field] <= 0:
                raise ConfigError(f"species.{key}.{field} must be a positive integer")
        if spec["string_taxon_id"] in string_taxa:
            raise ConfigError(f"species.{key}: string_taxon_id {spec['string_taxon_id']} is used twice")
        string_taxa.add(spec["string_taxon_id"])


def config_fingerprint(config: dict[str, Any]) -> str:
    """Stable short hash of the whole config except private _ keys. `mode` is included, so online
    and offline runs of otherwise identical settings get different run-id suffixes."""
    public = {k: v for k, v in config.items() if not k.startswith("_")}
    payload = json.dumps(public, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:12]
