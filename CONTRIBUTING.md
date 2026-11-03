# Contributing

Thank you for helping. This project values reproducibility and honest reporting over new
features, so a few rules apply to every change.

## Ground rules

- **No invented data.** Never fill gaps with guessed values or silently substitute an
  identifier. If a source is unavailable, report the limitation and fail with a clear message.
- **No overclaiming.** Results are computational candidates. Do not describe them as proof of
  biological equivalence, safety, therapeutic effect or experimental validation.
- **Keep provenance.** New data must go through the cached HTTP client
  (`src/data_collection/http_cache.py`), so every response is stored with its request,
  timestamp and SHA-256.
- **Keep outputs deterministic.** Seed every random step from the config, sort before summing,
  and write floats with 12 significant digits.
- **Curated evidence** (`config/curated_evidence*.tsv`) needs a PMID and DOI for every row, a
  paraphrased finding checked against the source, and the effect on the role (`supports` or
  `contradicts`). Contradicting findings are kept.

## Development setup

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip install -e . --no-deps
```

On Windows, use `py -3.11 -m venv .venv` and `.venv/Scripts/python`.

## Before opening a pull request

```bash
.venv/bin/python -m pytest
.venv/bin/python scripts/verify_checksums.py --all
```

If your change can affect results, also re-run Phase 2 offline and compare it with the reference
run (this is what CI does):

```bash
.venv/bin/python -m extremonet.analysis --config config/phase2_analysis.yaml --mode offline
.venv/bin/python scripts/check_reproduction.py results/phase2/20261102T165336Z_e64fca252fa7 results/phase2/<new_run_id>
```

A change that alters results on purpose needs new reference runs and updated findings documents
in the same pull request, with the reason stated.

Tests must not contact real APIs. Mock HTTP with the `FakeSession` in `tests/conftest.py`. A test
that genuinely needs the network must be marked `@pytest.mark.live`, so it can be deselected with
`-m "not live"` (CI does this).

## Licensing of contributions

Code contributions are accepted under the MIT licence (`LICENSE`). Documentation, figures and
result tables are accepted under CC BY 4.0 (`LICENSE-docs`).
