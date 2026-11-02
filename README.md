# Cross-Species Biological Network Discovery — Phase 1

A reproducible pipeline that builds **oxidative-stress defence** protein networks for
*Deinococcus radiodurans* (NCBI taxon 1299) and *Homo sapiens* (NCBI taxon 9606) from the
**STRING** functional-association database, with full provenance, record-level quality control
and weighted NetworkX graphs.

- **Phase 1:** per-species network construction, QC and provenance (sections below).
- **Phase 2:** cross-species analysis — topology, Louvain modules with enrichment, homology and
  annotation evidence, IsoRank alignment, transparent candidate ranking, NetMF embeddings and
  robustness checks. See [Phase 2](#phase-2-cross-species-analysis),
  [docs/PHASE2_METHODS.md](docs/PHASE2_METHODS.md) and
  [docs/PHASE2_FINDINGS.md](docs/PHASE2_FINDINGS.md).

- **DNA repair:** the second process, run through the same Phase 1 and Phase 2 pipeline with its
  own configuration ([DNA repair](#dna-repair-second-process),
  [docs/DNA_REPAIR_FINDINGS.md](docs/DNA_REPAIR_FINDINGS.md)).

> **Research plan.** The guiding document (Prof. Pietro Hiram Guzzi, *2026-Graph-AnanthapadmanabhanDS*)
> and its Overleaf source were not available in this workspace when the pipeline was built, so the
> design follows the written project brief. Check the plan against `docs/METHODS.md` once the
> document is available; anything that differs is a configuration change, not a code change.

## Quick start

Requires Python 3.11. From the project root:

```bash
py -3.11 -m venv .venv                     # Windows (macOS/Linux: python3.11 -m venv .venv)
.venv/Scripts/python -m pip install -r requirements.txt -e .   # macOS/Linux: .venv/bin/python
```

**Run the pipeline (single command):**

```bash
.venv/Scripts/python -m extremonet --config config/phase1_oxidative_stress.yaml
```

Options: `--mode online|refresh|offline` overrides the config's acquisition mode, and
`--threshold 0.9` overrides the combined-score edge threshold. Exit codes: `0` success,
`2` configuration error, `3` data acquisition failed (the message states the recovery path).

**Run fully offline** from previously acquired raw data. Nothing touches the network, and the
run fails loudly if any required response is missing:

```bash
.venv/Scripts/python -m extremonet --mode offline
```

**Run the tests** (all HTTP is mocked, so no network is needed):

```bash
.venv/Scripts/python -m pytest
```

Test temporary files are written to `.pytest_tmp/` inside the project (gitignored), not to the
system temp directory.

**Verify checksums** of the latest run and of every cached raw response:

```bash
.venv/Scripts/python scripts/verify_checksums.py
```

## What the pipeline does

1. **Seed proteins (UniProtKB REST).** Oxidative-stress defence is defined operationally as
   UniProtKB entries annotated with GO:0006979 (response to oxidative stress), GO:0098869
   (cellular oxidant detoxification) or GO:0016209 (antioxidant activity), via UniProt's `go:`
   query field. Human uses Swiss-Prot (reviewed) entries only. *D. radiodurans* uses all entries,
   because only about 526 of its 3,104 entries are reviewed.
2. **Identifier mapping (STRING `get_string_ids`).** Accessions are mapped using the configured
   NCBI taxon. The STRING taxon that comes back is **checked against the expected value**, never
   assumed (see *Species handling* below).
3. **Networks (STRING `network`).** Functional and physical subnetworks among the seeds are
   downloaded at STRING `required_score=400`. The analysis threshold of **0.7** is then applied
   locally, so it can be changed offline and every below-threshold record is reported.
4. **QC.** This step checks identifier formats, species assignment, duplicate edges (unordered
   pairs), self-loops, missing or out-of-range scores, unmapped and collapsed identifiers,
   isolated nodes, annotation coverage and network coverage. Every excluded record is written
   out with a reason.
5. **Graphs.** The pipeline builds one undirected weighted graph per species (weight = STRING
   combined score) with all channel scores, evidence flags and protein annotations attached.
   Each graph is written as GraphML plus TSV tables.
6. **Provenance.** Each run writes a manifest with the STRING release, the exact requests, the
   retrieval timestamps, SHA-256 checksums of all raw and generated files, software versions and
   the full config.

## Functional associations vs physical interactions

STRING edges are **functional associations**: evidence from several channels (co-expression,
text mining, genomic context, databases, experiments) combined into one score. They are *not*
physical interactions. Every edge carries:

| attribute | meaning |
|---|---|
| `combined_score` / `weight` | STRING integrated score (0–1) |
| `experimental_score`, `database_score`, `textmining_score`, … | individual STRING channel scores |
| `in_physical_subnetwork`, `physical_score` | the pair also appears in STRING's *physical* subnetwork (evidence of being in a physical complex) at the same threshold |
| `has_experiments_channel_score` | STRING's experiments channel is > `network.experimental_score_threshold` (default 0) |

Neither flag means "experimentally confirmed physical interaction". STRING's physical
subnetwork is still a scored prediction, and the experiments channel can include evidence
transferred from orthologs in other organisms. Curated experimental interaction sets (e.g.
BioGRID or IntAct) are a candidate Phase 2 data source.

## Species handling

STRING v12.5 represents *D. radiodurans* by the sequenced strain **R1, taxon 243230**, not the
species taxon 1299. When queried with 1299, STRING resolves to 243230. The config records both
values (`ncbi_taxon_id: 1299`, `string_taxon_id: 243230`), and QC confirms that every retained
node and edge carries the expected STRING taxon prefix. The resolution is reported in each QC
report and in the manifest. Human is 9606 in both systems.

## Outputs

Each run gets an id `<UTC timestamp>_<config fingerprint>`. `results/LATEST_RUN.txt` names the most recent one.

```
data/raw/<source>/<endpoint>/<request-key>__<timestamp>.{json,tsv}  raw API body, exactly as received
data/raw/<source>/<endpoint>/<request-key>__<timestamp>.meta.json   request, params, time, status, headers, sha256
data/processed/<run_id>/<species>/
    uniprot_seed_records_parsed.tsv    every seed record retrieved
    string_id_mapping_parsed.tsv       every mapping row (pass = batch | single)
    string_{functional,physical}_records_parsed.tsv   every STRING record retrieved
    nodes.tsv                          retained proteins with annotations
    edges_functional_retained.tsv      retained edges with evidence flags
    edges_physical_retained.tsv
    excluded_records.tsv               every excluded seed / mapping / edge record + reason
data/manifests/<run_id>.json           provenance manifest (also copied into results)
results/<run_id>/
    config_used.yaml  manifest.json  SHA256SUMS.txt  run_report.md  species_comparison.tsv
    logs/pipeline.jsonl                structured JSON-lines log
    <species>/network.graphml  weighted_edge_list.tsv  protein_table.tsv  degree_table.tsv
              summary_statistics.json  qc_report.json  qc_report.md
```

## Reference run (STRING 12.5, UniProt release 2026_03, retrieved 2026-10-25)

These are the observed outputs of the pipeline, not biological conclusions.

| | *D. radiodurans* (STRING 243230) | *H. sapiens* (9606) |
|---|---|---|
| UniProtKB entries for organism | 3,104 (all) | 20,431 (reviewed) |
| GO-defined seed records | 16 | 385 |
| seeds mapped to STRING | 16 | 383 accessions → 382 proteins |
| STRING functional records retrieved (≥ 0.4) | 38 | 5,550 |
| edges retained (combined score ≥ 0.7) | 11 | 2,210 |
| isolated nodes (kept, flagged) | 9 of 16 | 70 of 382 |
| edges also in STRING physical subnetwork | 4 | 1,152 |

### Known limitations of this run

- **Seed-set asymmetry.** The *D. radiodurans* seed set (16 proteins, 0.5 % of its entries) is far
  smaller than the human one (385 proteins), and only 4 of the 16 are Swiss-Prot reviewed. The
  bacterial network is therefore very sparse: 11 edges, with 9 of 16 seeds isolated at 0.7. Much
  of this reflects GO annotation depth, not only biology, so the two graphs should not be compared
  directly as if they had equal coverage.
- **True STRING gaps.** STRING 12.5 has no entry for human P07203 (GPX1) or P30711 (GSTT1); both
  are reported as `unmapped_in_string`. P0DMV8/P0DMV9 (HSPA1A/HSPA1B) share one STRING protein
  and are collapsed into one node, with both accessions kept.
- **No neighbours added.** With `add_nodes: 0`, only seed–seed associations are retrieved.
- **No orthology.** No cross-species orthology mapping is done. `species_comparison.tsv` is a
  descriptive side-by-side table only.

## Recovery path when an API is unavailable

- HTTP 429/5xx errors, timeouts and connection errors are retried with exponential backoff,
  honouring `Retry-After`. Other 4xx errors fail immediately with the server's message.
- If retries are exhausted, the run stops with exit code 3 and names the failing request.
  Nothing is substituted silently.
- To recover, re-run later, or copy a previously acquired `data/raw/` directory into place and run
  `--mode offline`. If the cached bodies were altered, their checksums no longer match and the
  run stops.
- Optional: setting `http.allow_stale_cache_on_failure: true` lets a failed *refresh* fall back to
  an older cached copy. That fallback is logged as a WARNING and recorded as
  `served_from: stale_cache` in the manifest.

## Phase 2: cross-species analysis

Phase 2 reads one Phase 1 run (`phase1.run_id` in `config/phase2_analysis.yaml`) and re-verifies
its checksums before doing anything else. It then retrieves the additional evidence it needs:

- STRING functional annotation, best-hit homology and PPI-enrichment results;
- UniProt EC numbers, PANTHER families and Pfam domains;
- QuickGO evidence codes;
- KEGG Orthology (KO) assignments (academic use);
- experimentally determined interactions between seed proteins (IntAct).

It also reads a literature-curated evidence table, `config/curated_evidence.tsv`, in which every
row is verified against its PubMed record.

These responses are cached with provenance, exactly as in Phase 1.

**Run Phase 2 (single command):**

```bash
.venv/Scripts/python -m extremonet.analysis --config config/phase2_analysis.yaml
```

Options are `--mode online|refresh|offline` and `--seed N` (overrides the master random seed).
Exit codes are as in Phase 1, plus `4` when the Phase 1 inputs fail verification. A run takes about
2 minutes once the evidence is cached.

Reproduce the reference analysis offline and check random-seed sensitivity:

```bash
.venv/Scripts/python -m extremonet.analysis --mode offline
```

```bash
.venv/Scripts/python -m extremonet.analysis --mode offline --seed 1
```

```bash
.venv/Scripts/python scripts/compare_phase2_runs.py <reference_run_id> <seed_run_id>
```

Outputs, under `results/phase2/<run_id>/`:

```
technical_report.md        generated report (every number read from this run's tables)
analysis_summary.json      machine-readable summary   config_used.yaml  manifest.json  SHA256SUMS.txt
<species>/                 node_metrics.tsv (centralities, hub/bottleneck class, module role, stability),
                           graph_statistics.json, degree_distribution.tsv, louvain_stability.tsv,
                           modules.tsv, module_enrichment.tsv (hypergeometric + BH-FDR), network_with_metrics.graphml
cross_species/             annotation overlap, shared terms, topology comparison, size-matched null,
                           module correspondence, homology best hits, IsoRank alignment pairs and summary
candidates/                candidate_pairs_all.tsv (every pair: 8 evidence components, score, coverage,
                           category, supporting / contradictory evidence, references, uncertainty,
                           proposed validation, rank interval), top3_human_per_dr_protein.tsv,
                           protein_evidence.tsv, category_counts.tsv, weights_leave_one_out.tsv
evaluation/                threshold, similarity-measure and ranking sensitivity; link prediction;
                           proxy evaluation
evidence/                  QuickGO, UniProt and STRING annotation tables used
coverage_variant/          KO-projected D. radiodurans network, metrics and alignment (coverage check only)
figures/                   degree CCDF, threshold sensitivity, size-matched context, Louvain stability,
                           ranking weight sensitivity, D. radiodurans module network, link-prediction AUC
```

**Reference Phase 2 run.** `20261102T165336Z_e64fca252fa7` is the reference run. Its
interpretation, including all coverage limitations, is in
[docs/PHASE2_FINDINGS.md](docs/PHASE2_FINDINGS.md). Briefly:

- the five reciprocal-best-hit orthologs and three further KEGG-KO orthologs are recovered;
- the best-supported bacterial protein is peroxiredoxin Q (DR0846);
- the strongest *non-homologous* candidates are Mn/Fe vs Cu/Zn superoxide dismutases, with
  conflicting sodA mutant evidence flagged;
- network alignment shows no topology conservation beyond a degree-preserving null, even with a
  KO-projected (enlarged) bacterial network;
- NetMF embeddings do not beat neighbourhood baselines and are not used in the ranking.

No result here is experimental validation.

## DNA repair (second process)

The DNA-repair networks use the same code with separate configuration files:

- `config/phase1_dna_repair.yaml`: seeds from GO:0006281 (DNA repair), with results in
  `results/dna_repair/`;
- `config/phase2_dna_repair.yaml`: identical analysis parameters to oxidative stress, with results
  in `results/dna_repair/phase2/`;
- `config/curated_evidence_dna_repair.tsv`: PubMed-verified mutant findings for *D. radiodurans*.

The raw-data cache in `data/raw/` is shared, because identical requests are fetched only once.

```bash
.venv/Scripts/python -m extremonet --config config/phase1_dna_repair.yaml
```

```bash
.venv/Scripts/python -m extremonet.analysis --config config/phase2_dna_repair.yaml
```

```bash
.venv/Scripts/python scripts/compare_phase2_runs.py --root results/dna_repair/phase2 <reference_run_id> <seed_run_id>
```

Findings: [docs/DNA_REPAIR_FINDINGS.md](docs/DNA_REPAIR_FINDINGS.md).

## Repository layout

```
config/                       YAML configuration (all parameters)
src/                          package `extremonet` (installed from src/ via pyproject.toml)
  data_collection/            cached HTTP client, STRING and UniProt clients
  preprocessing/              raw response -> typed tables (no rows dropped)
  quality_control/            validators (retain/exclude with reasons), reports
  annotation/                 protein table, annotation coverage
  network_construction/       weighted graph build, summary statistics
  analysis/                   Phase 2: inputs, topology, communities, enrichment, similarity,
                              alignment, ranking, embedding, figures, report, phase2 (orchestrator)
  pipeline.py  cli.py  config.py  provenance.py  logging_utils.py
tests/                        pytest suite; HTTP mocked via a fake session
scripts/                      run_pipeline.py, verify_checksums.py, compare_phase2_runs.py
docs/METHODS.md               Phase 1 methods, data model and QC rules
docs/PHASE2_METHODS.md        Phase 2 methods
docs/PHASE2_FINDINGS.md       Phase 2 interpretation and limitations (oxidative stress)
docs/DNA_REPAIR_FINDINGS.md   Phase 2 interpretation and limitations (DNA repair)
```

## Data sources and licences

- STRING (https://string-db.org, API docs: https://string-db.org/help/api/), CC BY 4.0.
  Requests identify the client with `caller_identity` and wait at least 1 s between calls, as
  the documentation asks.
- UniProtKB REST API (https://www.uniprot.org/help/api), CC BY 4.0.
- QuickGO (https://www.ebi.ac.uk/QuickGO/api), GO annotations, CC BY 4.0.
- KEGG REST (https://www.kegg.jp/kegg/rest/keggapi.html) for KO assignments: free for academic
  use; commercial use needs a licence. Pathway maps are taken from STRING's redistribution.
- IntAct (https://www.ebi.ac.uk/intact), via PSICQUIC; CC BY 4.0.
