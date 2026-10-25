"""Client for the documented STRING REST API (https://string-db.org/help/api/).

Endpoints used:
  * ``version``        - STRING release, recorded in every manifest
  * ``get_string_ids`` - map UniProt accessions to STRING identifiers
  * ``network``        - functional (default) or physical association subnetwork among identifiers

STRING edges are *functional associations* (co-expression, text mining, genomic context, database
and experimental channels combined). ``network_type=physical`` restricts to evidence that the
proteins are part of a physical complex; this is still a scored prediction, not a list of
experimentally confirmed binary interactions.
"""

from __future__ import annotations

import json
from typing import Any

from .http_cache import CachedHttpClient, DataSourceError, RawResponse


class StringClient:
    def __init__(self, http: CachedHttpClient, base_url: str = "https://string-db.org/api",
                 caller_identity: str = "extremonet") -> None:
        self.http = http
        self.base_url = base_url.rstrip("/")
        self.caller_identity = caller_identity

    def _url(self, endpoint: str) -> str:
        return f"{self.base_url}/json/{endpoint}"

    @staticmethod
    def parse_json(raw: RawResponse) -> list[dict[str, Any]]:
        try:
            payload = json.loads(raw.text) if raw.text.strip() else []
        except json.JSONDecodeError as exc:
            raise DataSourceError(f"STRING returned non-JSON content ({raw.data_path}): {exc}") from exc
        if isinstance(payload, dict):
            # STRING reports problems as a JSON object, e.g. {"Error": "...", "ErrorMessage": "..."}
            raise DataSourceError(f"STRING returned an error object ({raw.data_path}): {payload}")
        if not isinstance(payload, list):
            raise DataSourceError(f"unexpected STRING payload type {type(payload).__name__}")
        return payload

    def version(self) -> tuple[dict[str, Any], RawResponse]:
        raw = self.http.fetch("string", "version", self._url("version"),
                              {"caller_identity": self.caller_identity}, method="GET", ext="json",
                              description="STRING database version")
        rows = self.parse_json(raw)
        if not rows:
            raise DataSourceError("STRING version endpoint returned no rows")
        return rows[0], raw

    def map_identifiers(self, identifiers: list[str], species: int, species_key: str,
                        endpoint_suffix: str = "") -> tuple[list[dict], RawResponse]:
        """Note: for a batch, STRING returns at most one row per STRING protein, so when several
        accessions resolve to the same protein only one of them is echoed back. Callers re-query
        the missing accessions individually to tell such collapses apart from true gaps."""
        params = {
            "identifiers": "\r".join(sorted(set(identifiers))),
            "species": species,
            "limit": 1,
            "echo_query": 1,
            "caller_identity": self.caller_identity,
        }
        raw = self.http.fetch("string", f"get_string_ids_{species_key}{endpoint_suffix}", self._url("get_string_ids"),
                              params, method="POST", ext="json",
                              description=f"map {len(set(identifiers))} identifiers to STRING (species={species})")
        return self.parse_json(raw), raw

    def network(self, string_ids: list[str], species: int, required_score: int, network_type: str,
                add_nodes: int, species_key: str) -> tuple[list[dict], RawResponse]:
        params = {
            "identifiers": "\r".join(sorted(set(string_ids))),
            "species": species,
            "required_score": required_score,
            "network_type": network_type,
            "add_nodes": add_nodes,
            "caller_identity": self.caller_identity,
        }
        raw = self.http.fetch("string", f"network_{network_type}_{species_key}", self._url("network"),
                              params, method="POST", ext="json",
                              description=(f"STRING {network_type} network among {len(set(string_ids))} "
                                           f"proteins, required_score={required_score}"))
        return self.parse_json(raw), raw
