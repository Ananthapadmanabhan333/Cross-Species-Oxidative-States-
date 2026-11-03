# ExtremoNet: cross-species network analysis of *Deinococcus radiodurans* and human stress defences

[![CI](https://github.com/Ananthapadmanabhan333/Cross-Species-Oxidative-States-/actions/workflows/ci.yml/badge.svg)](https://github.com/Ananthapadmanabhan333/Cross-Species-Oxidative-States-/actions/workflows/ci.yml)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
[![Code: MIT](https://img.shields.io/badge/code-MIT-green)](LICENSE)
[![Docs: CC BY 4.0](https://img.shields.io/badge/docs%20%26%20results-CC%20BY%204.0-lightgrey)](LICENSE-docs)

*Deinococcus radiodurans* survives radiation and desiccation that kill human cells, yet many of
its protective proteins have no human ortholog. This project asks whether protein-network
analysis can point to human proteins that perform **analogous functions without shared
ancestry**, for two processes: **oxidative-stress defence** and **DNA repair**.

The pipeline builds STRING functional-association networks for both species, compares them by
topology, modules and network alignment against degree-preserving null models, and ranks every
cross-species protein pair with a transparent, evidence-graded score. Every input is cached and
checksummed, and the reference analyses are reproduced offline in CI.

> **All results are computational candidates.** None demonstrates biological equivalence, safety,
> therapeutic effect or experimental validation.

## Key findings

| | Oxidative stress | DNA repair |
|---|---|---|
| Networks (*D. radiodurans* / human proteins, STRING score ≥ 0.7) | 16 / 382 | 52 / 426 |
| Topology conserved beyond degree-preserving nulls? | no (p = 0.10) | borderline (p = 0.05); with homology prior p = 0.01–0.02 |
| Clearest candidate analogues | Mn/Fe vs Cu/Zn superoxide dismutases | Fpg/MutM vs HhH-GPD DNA glycosylases (e.g. MutM–OGG1) |
| Experimentally confirmed bacterial interactions (IntAct, BioGRID, MINT) | none | none |

- **Orthologs are recovered** as a sanity check (e.g. *katA*–CAT, *sodA*–SOD2, *nth*–NTHL1,
  *ung*–UNG), not as a discovery.
- **Topology alone does not identify equivalent proteins.** Proposed analogues rest on shared
  catalytic activity with independent support, not on network position.
- **Curated mutant studies change interpretations.** A ΔMnSOD mutant survives gamma radiation
  *better*, pointing to non-protein Mn²⁺ antioxidants; *lexA* deletion does not affect survival
  despite its "DNA repair" annotation; PprA acts in DNA segregation and has no supported human
  counterpart.

Details: [oxidative stress](docs/PHASE2_FINDINGS.md) · [DNA repair](docs/DNA_REPAIR_FINDINGS.md) ·
[methods](docs/PHASE2_METHODS.md) · [paper draft (IEEE format)](docs/overleaf/paper.tex) ·
[combined report (HTML)](docs/report/deinococcus_human_networks.html)

## Status against the research plan

The project follows the research plan by Prof. Pietro Hiram Guzzi (*2026-Graph-AnanthapadmanabhanDS*).

| Phase | Done | Not done or done differently |
|---|---|---|
| 1. Network construction | STRING ≥ 0.7 networks for oxidative stress and DNA repair; GO/KEGG annotation; orthology as a feature | broader stress-response subsystem; OrthoFinder (graded homology from reciprocal best hits, KEGG KO, PANTHER, Pfam and InterPro instead) |
| 2. Topology | degree, betweenness, closeness, clustering, hub/bottleneck classes; Louvain with stability | graphlet degree vectors; Leiden; transcriptomic overlay |
| 3. Cross-species comparison | module correspondence; IsoRank alignment with degree-preserving nulls; NetMF embeddings (within species) | module-preservation statistics; MAGNA++ or graphlet alignment; cross-species embeddings |
| 4. Interpretability | transparent weighted score, rule-based categories, sensitivity analyses | — |
| 5. Validation prioritisation | shortlist with proposed experiments; curated mutant evidence | cross-check against human genetic constraint and disease associations |

## Quick start

Requires Python 3.11.

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip install -e . --no-deps
```

On Windows, use `py -3.11 -m venv .venv` and `.venv/Scripts/python` in the commands below.

Run the tests (all HTTP is mocked, no network needed):

```bash
.venv/bin/python -m pytest
```

Verify every reference run and every cached API response against its recorded SHA-256:

```bash
.venv/bin/python scripts/verify_checksums.py --all
```

## Reproducing the reference results

All raw API responses are committed under `data/raw/`, so the analyses re-run offline:

| Process | Phase 1 (networks) | Phase 2 (analysis) |
|---|---|---|
| Oxidative stress | `results/20261025T174736Z_5bbc7e723df4` | `results/phase2/20261102T165336Z_e64fca252fa7` |
| DNA repair | `results/dna_repair/20261102T161839Z_7128c0f5ca4e` | `results/dna_repair/phase2/20261102T170150Z_5a2a85f66091` |

```bash
.venv/bin/python -m extremonet.analysis --config config/phase2_analysis.yaml --mode offline
```

```bash
.venv/bin/python scripts/check_reproduction.py results/phase2/20261102T165336Z_e64fca252fa7 results/phase2/<new_run_id>
```

`check_reproduction.py` requires every result table to be byte-identical and every JSON summary
to agree to 12 significant digits. CI runs this for both processes on every push. Each Phase 2
run takes about two minutes.

Random-seed sensitivity (two alternative seeds per process are committed):

```bash
.venv/bin/python -m extremonet.analysis --mode offline --seed 1
```

```bash
.venv/bin/python scripts/compare_phase2_runs.py <reference_run_id> <seed_run_id>
```

## Running new analyses

```bash
.venv/bin/python -m extremonet --config config/phase1_oxidative_stress.yaml
```

```bash
.venv/bin/python -m extremonet.analysis --config config/phase2_analysis.yaml
```

The DNA-repair equivalents use `config/phase1_dna_repair.yaml` and `config/phase2_dna_repair.yaml`.

- `--mode online|refresh|offline` sets data acquisition: use the cache when possible, always
  re-fetch, or never touch the network.
- `--threshold 0.9` (Phase 1) overrides the STRING combined-score threshold.
- `--seed N` (Phase 2) overrides the master random seed.
- Exit codes: `0` success, `2` configuration error, `3` data acquisition failed (the message
  states the recovery path), `4` Phase 1 inputs failed checksum verification.

All parameters live in `config/`, and each run copies its configuration into its manifest.

## How it works

### Phase 1: networks

1. **Seed proteins (UniProtKB).** Entries annotated with the process GO terms, including
   descendants: GO:0006979, GO:0098869 and GO:0016209 for oxidative stress; GO:0006281 for DNA
   repair. Human uses Swiss-Prot (reviewed) entries; *D. radiodurans* uses all entries, because
   only about 526 of its 3,104 entries are reviewed.
2. **Identifier mapping (STRING `get_string_ids`).** The STRING taxon returned is checked, never
   assumed. Accessions missing from a batch answer are re-queried singly, which separates true
   gaps from accessions that share a STRING protein.
3. **Networks (STRING `network`).** Retrieved at score ≥ 0.4 and thresholded at 0.7 locally, so
   other thresholds can be applied offline.
4. **Quality control.** Ordered checks on identifiers, species, duplicates, self-loops and score
   ranges. Every excluded record is written out with its reason.
5. **Graphs and provenance.** One weighted NetworkX graph per species (GraphML and TSV), and a
   manifest with requests, timestamps, software versions and SHA-256 checksums.

Rules and data model: [docs/METHODS.md](docs/METHODS.md).

### Phase 2: cross-species analysis

- **Evidence:** STRING annotation, homology and PPI enrichment; UniProt EC, PANTHER and Pfam;
  QuickGO evidence codes; KEGG Orthology; IntAct physical interactions; literature-curated
  mutant studies (`config/curated_evidence*.tsv`, each row checked against PubMed).
- **Topology and modules:** centralities, hub/bottleneck classes, Louvain over 50 seeds with
  adjusted Rand stability, hypergeometric enrichment with Benjamini–Hochberg FDR.
- **Graded homology:** reciprocal best hit (1.0) > KEGG KO (0.9) > PANTHER (0.8) > one-way hit
  (0.6) > Pfam (0.5) > InterPro (0.4) > common fold only (0.2).
- **Alignment:** IsoRank with and without a homology prior, tested against random mappings and
  degree-preserving rewired networks.
- **Ranking:** weighted mean over the evidence components available for each pair, so missing
  evidence never counts as negative. Rule-based categories: ortholog, homolog, distant homolog,
  candidate analogue, shared function, similar pathway role, speculative.
- **Robustness:** weight perturbation, leave-one-component-out, similarity measure, STRING
  threshold, random and hash seeds; NetMF embeddings as a baseline.

Full description: [docs/PHASE2_METHODS.md](docs/PHASE2_METHODS.md).

### Functional associations vs physical interactions

STRING edges are **functional associations**, combining co-expression, text mining, genomic
context, databases and experiments, and partly transferred between organisms through orthology.
They are not physical interactions. Edges also in STRING's physical subnetwork are flagged, and
experimentally determined physical interactions are checked separately in IntAct. IntAct, BioGRID
and MINT hold none for the *D. radiodurans* seed proteins.

### Species handling

STRING v12.5 represents *D. radiodurans* by strain R1 (taxon 243230). Queried with the species
taxon 1299, it resolves to 243230. The configuration records both values, and QC confirms that
every node and edge carries the expected taxon. Human is 9606 in both systems.

### When an API is unavailable

HTTP 429/5xx errors, timeouts and connection errors are retried with exponential backoff,
honouring `Retry-After`. If retries are exhausted, the run stops with exit code 3 and names the
failing request; nothing is substituted silently. Recover by re-running later, or by restoring a
`data/raw/` directory and running `--mode offline`. Optionally,
`http.allow_stale_cache_on_failure: true` lets a failed refresh fall back to an older cached copy,
logged as a warning and recorded as `served_from: stale_cache`.

## Outputs

Each run gets an id `<UTC timestamp>_<config fingerprint>`; `LATEST_RUN.txt` names the most recent.

```
data/raw/<source>/<endpoint>/<key>__<timestamp>.{json,tsv}   raw API body, exactly as received
data/raw/<source>/<endpoint>/<key>__<timestamp>.meta.json    request, time, status, headers, sha256
data/processed/<run_id>/<species>/                           parsed tables, retained and excluded records
results/<run_id>/                                            Phase 1: graphs, QC reports, manifest, SHA256SUMS
results/phase2/<run_id>/
    technical_report.md      generated report (every number read from this run's tables)
    <species>/               node metrics, modules, enrichment, graph statistics
    cross_species/           annotation overlap, module correspondence, homology, alignment
    candidates/              every pair: 8 evidence components, score, category, evidence,
                             references, uncertainty, proposed validation, rank interval
    evaluation/              sensitivity analyses, link prediction, proxy evaluation
    coverage_variant/        KO-projected bacterial network (coverage check only)
    figures/                 degree distributions, sensitivity, stability, module network
results/dna_repair/...       the same layout for DNA repair
```

## Repository layout

```
config/                  YAML configuration and curated evidence tables
src/                     package `extremonet` (installed from src/)
  data_collection/       cached HTTP client; STRING, UniProt, QuickGO, KEGG, IntAct clients
  preprocessing/         raw responses -> typed tables (no rows dropped)
  quality_control/       validators (retain or exclude, with reasons), reports
  annotation/            protein table, annotation coverage
  network_construction/  weighted graph build, summary statistics
  analysis/              Phase 2: topology, communities, enrichment, similarity, alignment,
                         ranking, embedding, coverage, figures, report
tests/                   pytest suite; HTTP answered by a fake session
scripts/                 run_pipeline.py, verify_checksums.py, check_reproduction.py,
                         compare_phase2_runs.py
docs/                    methods, findings, combined HTML report, paper draft (docs/overleaf/)
data/, results/          cached API responses and reference runs (checksummed)
.github/workflows/       CI: tests and checksums on Linux and Windows, offline reproduction
```

## Citing

See [CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button uses it.

## Licence

- **Code** (`src/`, `tests/`, `scripts/`, configuration): [MIT](LICENSE).
- **Documentation, figures, curated evidence tables and generated results** (`docs/`, `results/`,
  `config/curated_evidence*.tsv`): [CC BY 4.0](LICENSE-docs).
- **Cached third-party data** (`data/raw/`) remains under the terms of its providers, listed below.

## Data sources and their terms

- STRING (https://string-db.org), CC BY 4.0. Requests identify the client with `caller_identity`
  and wait at least 1 s between calls, as the API documentation asks.
- UniProtKB (https://www.uniprot.org/help/api), CC BY 4.0.
- QuickGO (https://www.ebi.ac.uk/QuickGO/api), GO annotations, CC BY 4.0.
- IntAct (https://www.ebi.ac.uk/intact) via PSICQUIC, CC BY 4.0.
- KEGG REST (https://www.kegg.jp/kegg/rest/keggapi.html), KO assignments: free for academic use;
  commercial use needs a licence (https://www.kegg.jp/kegg/legal.html). Pathway maps come from
  STRING's redistribution of KEGG.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). In short: no invented data, no overclaiming, keep
provenance and determinism, and run the tests and checks before opening a pull request.
