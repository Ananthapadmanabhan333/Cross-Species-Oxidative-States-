"""Shared fixtures. All HTTP traffic is answered by FakeSession; no test needs the network.

The fixture records below are SYNTHETIC (made-up accessions, gene names and scores) and exist only
to exercise the code paths. They are not biological data.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import pytest
import requests

from extremonet.config import load_config

REPO_ROOT = Path(__file__).resolve().parents[1]
UNIPROT_HEADER = ("Entry\tEntry Name\tReviewed\tGene Names (primary)\tGene Names\tProtein names\t"
                  "Organism (ID)\tOrganism\tSTRING\tGene Ontology IDs\tAnnotation")


class FakeResponse:
    def __init__(self, status_code: int = 200, body: str | bytes = "", headers: dict[str, str] | None = None):
        self.status_code = status_code
        self.content = body if isinstance(body, bytes) else body.encode("utf-8")
        self.text = self.content.decode("utf-8")
        self.headers = {k.lower(): v for k, v in (headers or {}).items()}
        # requests' CaseInsensitiveDict semantics for .get
        self.headers = _CIDict(self.headers)


class _CIDict(dict):
    def get(self, key, default=None):  # type: ignore[override]
        return super().get(key.lower(), default)


class FakeSession:
    """Routes requests to a handler; records every call. A handler may return a FakeResponse or
    raise (e.g. requests.ConnectionError)."""

    def __init__(self, handler: Callable[[str, str, dict[str, Any]], FakeResponse]):
        self.handler = handler
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def request(self, method: str, url: str, **kwargs: Any) -> FakeResponse:
        params = kwargs.get("data") or kwargs.get("params") or {}
        self.calls.append((method, url, dict(params)))
        return self.handler(method, url, params)


def uniprot_row(acc, entry, reviewed, gene, protein, organism_id, organism, xref, go_ids, score="3.0"):
    return "\t".join([acc, entry, reviewed, gene, gene, protein, str(organism_id), organism, xref, go_ids, score])


# --------------------------------------------------------------------------- synthetic "databases"
DR = {"ncbi": 1299, "string": 243230}
HS = {"ncbi": 9606, "string": 9606}

DR_SEEDS = [
    uniprot_row("P00001", "TSTA_DEIRA", "reviewed", "tstA", "Test protein A", 243230, "D. radiodurans R1",
                "243230.DR_T001;", "GO:0006979; GO:0005737"),
    uniprot_row("P00002", "TSTB_DEIRA", "reviewed", "tstB", "Test protein B", 243230, "D. radiodurans R1",
                "243230.DR_T002;", "GO:0016209"),
    uniprot_row("Q00003", "Q00003_DEIRA", "unreviewed", "", "Test protein C", 243230, "D. radiodurans R1",
                "243230.DR_T003;", "GO:0034599"),          # matched via a descendant term
    uniprot_row("Q00004", "Q00004_DEIRA", "unreviewed", "tstD", "Test protein D", 243230, "D. radiodurans R1",
                "243230.DR_T004;", "GO:0006979"),          # will have no edge -> isolated
    uniprot_row("P00005", "TSTE_DEIRX", "reviewed", "tstE", "Other strain protein", 999999, "Deinococcus sp.",
                "", "GO:0006979"),                         # unexpected organism -> excluded
    uniprot_row("not-an-acc", "BAD", "unreviewed", "bad", "Malformed", 243230, "D. radiodurans R1",
                "", "GO:0006979"),                         # invalid accession -> excluded
]
HS_SEEDS = [
    uniprot_row("P10001", "TSTA_HUMAN", "reviewed", "HTSTA", "Human test A", 9606, "Homo sapiens",
                "9606.ENSPTEST01;", "GO:0006979"),
    uniprot_row("P10002", "TSTB_HUMAN", "reviewed", "HTSTB", "Human test B", 9606, "Homo sapiens",
                "9606.ENSPTEST02;", "GO:0098869"),
    uniprot_row("P10003", "TSTC_HUMAN", "reviewed", "HTSTC", "Human test C", 9606, "Homo sapiens",
                "9606.ENSPTEST99;", "GO:0016209"),        # xref disagrees with STRING mapping
    uniprot_row("P10004", "TSTD_HUMAN", "reviewed", "HTSTD", "Human test D", 9606, "Homo sapiens",
                "", "GO:0006979"),                         # STRING cannot map -> excluded
    uniprot_row("P10001", "TSTA_HUMAN", "reviewed", "HTSTA", "Human test A", 9606, "Homo sapiens",
                "9606.ENSPTEST01;", "GO:0006979"),         # duplicate accession -> excluded
    uniprot_row("P10005", "TSTE_HUMAN", "reviewed", "HTSTE", "Human test E", 9606, "Homo sapiens",
                "9606.ENSPTEST02;", "GO:0006979"),         # same STRING protein as P10002 -> collapsed
    uniprot_row("P10006", "TSTF_HUMAN", "reviewed", "HTSTF", "Human test F", 9606, "Homo sapiens",
                "9606.ENSPTEST06;", "GO:0006979"),         # mapped, no edges (isolated)
    uniprot_row("P10007", "TSTG_HUMAN", "reviewed", "HTSTG", "Human test G", 9606, "Homo sapiens",
                "9606.ENSPTEST07;", "GO:0006979"),         # mapped, no edges (isolated)
]
MAPPING = {
    "P00001": (243230, "243230.DR_T001", "tstA"), "P00002": (243230, "243230.DR_T002", "tstB"),
    "Q00003": (243230, "243230.DR_T003", "DR_T003"), "Q00004": (243230, "243230.DR_T004", "tstD"),
    "P10001": (9606, "9606.ENSPTEST01", "HTSTA"), "P10002": (9606, "9606.ENSPTEST02", "HTSTB"),
    "P10003": (9606, "9606.ENSPTEST03", "HTSTC"), "P10005": (9606, "9606.ENSPTEST02", "HTSTB"),
    "P10006": (9606, "9606.ENSPTEST06", "HTSTF"), "P10007": (9606, "9606.ENSPTEST07", "HTSTG"),
}


def edge(a, b, score, taxon, escore=0.0, names=("x", "y")):
    return {"stringId_A": a, "stringId_B": b, "preferredName_A": names[0], "preferredName_B": names[1],
            "ncbiTaxonId": taxon, "score": score, "nscore": 0, "fscore": 0, "pscore": 0, "ascore": 0.1,
            "escore": escore, "dscore": 0, "tscore": 0.5}


NETWORKS = {
    (243230, "functional"): [
        edge("243230.DR_T001", "243230.DR_T002", 0.95, 243230, escore=0.6),
        edge("243230.DR_T002", "243230.DR_T003", 0.75, 243230),
        edge("243230.DR_T001", "243230.DR_T003", 0.55, 243230),   # below 0.7 -> excluded
        edge("243230.DR_T002", "243230.DR_T001", 0.90, 243230),   # duplicate (reverse) -> excluded
        edge("243230.DR_T003", "243230.DR_T003", 0.99, 243230),   # self loop -> excluded
    ],
    (243230, "physical"): [edge("243230.DR_T001", "243230.DR_T002", 0.80, 243230, escore=0.6)],
    (9606, "functional"): [
        edge("9606.ENSPTEST01", "9606.ENSPTEST02", 0.99, 9606, escore=0.9),
        edge("9606.ENSPTEST02", "9606.ENSPTEST03", 0.72, 9606),
        edge("9606.ENSPTEST01", "9606.ENSPTEST03", None, 9606),     # missing score -> excluded
    ],
    (9606, "physical"): [],
}


def string_api_handler(method: str, url: str, params: dict[str, Any]) -> FakeResponse:
    if url.endswith("/json/version"):
        return FakeResponse(body=json.dumps([{"string_version": "12.5", "stable_address": "https://example"}]))
    if url.endswith("/json/get_string_ids"):
        rows, seen = [], set()
        for i, item in enumerate(params["identifiers"].split("\r")):
            # like STRING: a batch answer holds at most one row per STRING protein
            if item in MAPPING and MAPPING[item][1] not in seen:
                taxon, sid, name = MAPPING[item]
                seen.add(sid)
                rows.append({"queryIndex": i, "queryItem": item, "stringId": sid, "ncbiTaxonId": taxon,
                             "taxonName": "synthetic", "preferredName": name,
                             "annotation": "annotation not available" if name.startswith("DR_") else f"{name} synthetic"})
        return FakeResponse(body=json.dumps(rows))
    if url.endswith("/json/network"):
        rows = NETWORKS[(int(params["species"]), params["network_type"])]
        required = int(params["required_score"]) / 1000
        return FakeResponse(body=json.dumps([r for r in rows if r["score"] is None or r["score"] >= required]))
    if url.endswith("/uniprotkb/stream"):
        rows = DR_SEEDS if "1299" in params["query"] else HS_SEEDS
        return FakeResponse(body="\n".join([UNIPROT_HEADER, *rows]) + "\n")
    if url.endswith("/uniprotkb/search"):
        total = "3104" if "1299" in params["query"] else "20400"
        return FakeResponse(body="Entry\nP00001\n", headers={"X-Total-Results": total})
    return FakeResponse(status_code=404, body="not found")


@pytest.fixture
def fake_session() -> FakeSession:
    return FakeSession(string_api_handler)


@pytest.fixture
def offline_session() -> FakeSession:
    def refuse(method, url, params):
        raise requests.ConnectionError("network disabled in test")
    return FakeSession(refuse)


@pytest.fixture
def config(tmp_path: Path) -> dict[str, Any]:
    cfg = load_config(REPO_ROOT / "config" / "phase1_oxidative_stress.yaml")
    cfg["paths"] = {f"{name}_dir": str(tmp_path / "data" / name) for name in ("raw", "processed", "manifests")}
    cfg["paths"]["results_dir"] = str(tmp_path / "results")
    cfg["http"]["min_interval_seconds"] = 0
    cfg["http"]["backoff_seconds"] = 0
    return cfg
