import json

import pytest
import requests

from extremonet.data_collection.http_cache import (
    CacheIntegrityError, CachedHttpClient, DataSourceError, DataSourceUnavailable, OfflineCacheMiss,
    sha256_file,
)

from .conftest import FakeResponse, FakeSession

URL = "https://api.example/json/thing"


def client(tmp_path, session, mode="online", **kw):
    return CachedHttpClient(tmp_path / "raw", mode=mode, max_retries=kw.pop("max_retries", 2),
                            backoff_seconds=0, min_interval_seconds=0, session=session,
                            sleep=kw.pop("sleep", lambda s: None), **kw)


def ok_session(body='[{"a": 1}]'):
    return FakeSession(lambda m, u, p: FakeResponse(body=body, headers={"Content-Type": "application/json"}))


def test_fetch_stores_raw_body_and_provenance(tmp_path):
    session = ok_session()
    raw = client(tmp_path, session).fetch("src", "ep", URL, {"x": 1}, ext="json", description="d")
    assert raw.served_from == "network"
    assert raw.data_path.read_text() == '[{"a": 1}]'
    meta = json.loads(raw.meta_path.read_text())
    assert meta["sha256"] == sha256_file(raw.data_path)
    assert meta["params"] == {"x": "1"}
    assert meta["http_status"] == 200
    assert meta["retrieved_at_utc"]
    assert meta["response_headers"]["content-type"] == "application/json"


def test_online_mode_reuses_cache_without_network(tmp_path):
    session = ok_session()
    first = client(tmp_path, session).fetch("src", "ep", URL, {"x": 1})
    second = client(tmp_path, session).fetch("src", "ep", URL, {"x": 1})
    assert len(session.calls) == 1
    assert second.served_from == "cache"
    assert second.text == first.text


def test_different_parameters_are_cached_separately(tmp_path):
    session = ok_session()
    c = client(tmp_path, session)
    c.fetch("src", "ep", URL, {"x": 1})
    c.fetch("src", "ep", URL, {"x": 2})
    assert len(session.calls) == 2


def test_offline_mode_without_cache_fails_loudly(tmp_path, offline_session):
    with pytest.raises(OfflineCacheMiss, match="mode=online"):
        client(tmp_path, offline_session, mode="offline").fetch("src", "ep", URL, {})
    assert offline_session.calls == []


def test_offline_mode_serves_cache_and_never_calls_network(tmp_path, offline_session):
    client(tmp_path, ok_session()).fetch("src", "ep", URL, {})
    raw = client(tmp_path, offline_session, mode="offline").fetch("src", "ep", URL, {})
    assert raw.served_from == "cache"
    assert offline_session.calls == []


def test_refresh_keeps_previous_raw_files(tmp_path):
    client(tmp_path, ok_session("old")).fetch("src", "ep", URL, {})
    raw = client(tmp_path, ok_session("new"), mode="refresh").fetch("src", "ep", URL, {})
    assert raw.text == "new"
    assert len(list((tmp_path / "raw" / "src" / "ep").glob("*.meta.json"))) == 2
    # later online reads use the newest copy
    assert client(tmp_path, ok_session("unused")).fetch("src", "ep", URL, {}).text == "new"


def test_retries_transient_errors_then_succeeds(tmp_path):
    replies = iter([FakeResponse(503, "busy"), FakeResponse(429, "slow down"), FakeResponse(200, "ok")])
    waits = []
    session = FakeSession(lambda m, u, p: next(replies))
    raw = client(tmp_path, session, sleep=waits.append).fetch("src", "ep", URL, {})
    assert raw.text == "ok"
    assert len(session.calls) == 3
    assert len(waits) == 2


def test_connection_errors_and_timeouts_are_retried(tmp_path):
    errors = iter([requests.ConnectionError("down"), requests.Timeout("slow")])

    def handler(m, u, p):
        try:
            raise next(errors)
        except StopIteration:
            return FakeResponse(200, "ok")
    assert client(tmp_path, FakeSession(handler)).fetch("src", "ep", URL, {}).text == "ok"


def test_retry_after_header_is_respected(tmp_path):
    replies = iter([FakeResponse(429, "", {"Retry-After": "7"}), FakeResponse(200, "ok")])
    waits = []
    client(tmp_path, FakeSession(lambda m, u, p: next(replies)), sleep=waits.append).fetch("s", "e", URL, {})
    assert waits == [7.0]


def test_exhausted_retries_raise_unavailable(tmp_path):
    session = FakeSession(lambda m, u, p: FakeResponse(502, "bad gateway"))
    with pytest.raises(DataSourceUnavailable, match="mode=offline"):
        client(tmp_path, session, max_retries=2).fetch("src", "ep", URL, {})
    assert len(session.calls) == 3
    assert not (tmp_path / "raw" / "src" / "ep").exists()  # nothing cached from a failure


def test_client_errors_are_not_retried(tmp_path):
    session = FakeSession(lambda m, u, p: FakeResponse(400, "bad species"))
    with pytest.raises(DataSourceError, match="400"):
        client(tmp_path, session).fetch("src", "ep", URL, {})
    assert len(session.calls) == 1


def test_stale_cache_only_used_when_explicitly_allowed(tmp_path):
    client(tmp_path, ok_session("cached")).fetch("src", "ep", URL, {})
    failing = FakeSession(lambda m, u, p: FakeResponse(503, "down"))
    with pytest.raises(DataSourceUnavailable):
        client(tmp_path, failing, mode="refresh").fetch("src", "ep", URL, {})
    c = client(tmp_path, failing, mode="refresh", allow_stale_cache_on_failure=True)
    raw = c.fetch("src", "ep", URL, {})
    assert raw.text == "cached" and raw.served_from == "stale_cache"
    assert c.requests_log[-1]["served_from"] == "stale_cache"


def test_tampered_cache_is_detected(tmp_path):
    raw = client(tmp_path, ok_session("genuine")).fetch("src", "ep", URL, {})
    raw.data_path.write_text("edited by hand")
    with pytest.raises(CacheIntegrityError, match="checksum mismatch"):
        client(tmp_path, ok_session()).fetch("src", "ep", URL, {})


def test_post_sends_form_data(tmp_path):
    session = ok_session()
    client(tmp_path, session).fetch("src", "ep", URL, {"identifiers": "a\rb"}, method="POST")
    method, _, params = session.calls[0]
    assert method == "POST" and params == {"identifiers": "a\rb"}
