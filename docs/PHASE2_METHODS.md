# Phase 2 methods: cross-species analysis

All parameters are in `config/phase2_analysis.yaml` and are copied into every run's
`config_used.yaml` and manifest. Every stochastic step derives its seed from `random_seed` with
SHA-256 (`derived_seed(master, name)`), so a run is fully determined by its config.

## 0. Inputs and their verification

Phase 2 reads the processed tables of one Phase 1 run (`phase1.run_id`). Before any analysis, every
Phase 1 output and raw API response listed in that run's manifest is re-checksummed; a mismatch
stops the run (exit code 4). Graphs at other thresholds are rebuilt offline from the stored STRING
records (retrieved at score ≥ 0.4) using the same Phase 1 QC and graph code.

New evidence is retrieved via the Phase 1 cached HTTP client (same modes, retries, provenance):

| source | endpoint | used for |
|---|---|---|
| STRING | `functional_annotation` | per-protein GO (Process, Function, Component), KEGG maps, Pfam, InterPro, Reactome terms |
| STRING | `homology_best` (both directions) | best Smith–Waterman hits; reciprocal best hits (RBH) |
| STRING | `ppi_enrichment` (each threshold) | test of seed-set connectivity vs random genome sets |
| UniProtKB | `stream` (same seed query, extra fields) | EC numbers, PANTHER families, Pfam, KEGG genes, PubMed ids |
| QuickGO | `annotation/search` (paged) | GO evidence codes and references for the seed terms and descendants |
| KEGG REST | `link/ko/<org>` (once per organism) | KEGG Orthology (KO) group of every gene |
| curated table | `config/curated_evidence.tsv` | literature findings that GO evidence codes miss |

KEGG REST is called once per organism (`dra`, `hsa`) and cached. KEGG's terms allow free academic
use; commercial use needs a licence (https://www.kegg.jp/kegg/legal.html). Pathway maps come from
STRING's redistribution of KEGG.

**Curated evidence.** Each row of `config/curated_evidence.tsv` holds:

- the protein (which must be a node of its species' network; the run fails otherwise);
- the evidence type (`perturbation` or `biochemical`);
- the effect on the proposed role (`supports` or `contradicts`);
- the assay and a paraphrased finding;
- the PMID and DOI, and the curation basis (currently the PubMed abstract) and date.

Contradicting findings are kept, never filtered out.

## 1. Baseline topology (`analysis/topology.py`)

Node metrics: degree, weighted degree, unweighted and weighted clustering, betweenness and
closeness (edge length = 1 / combined score; closeness with the Wasserman–Faust component
correction), PageRank, k-core, component id and size. Percentile ranks are computed **within
species over connected nodes only**; isolated nodes get no percentiles and the class `isolated`.
Hubs are nodes at or above the `hub_percentile` degree percentile; bottlenecks are nodes at or
above the `bottleneck_percentile` betweenness percentile (Yu et al. 2007 classes). Guimerà–Amaral
within-module z-scores and participation coefficients use the primary Louvain partition.

Graph statistics are computed on the whole (disconnected) graph; path statistics on the largest
component. **Unequal coverage:** the larger network is subsampled `size_matched_samples` times to
the smaller network's node count, and the smaller network's statistics are located in that
distribution. This is descriptive context, not a hypothesis test.

## 2. Communities (`analysis/communities.py`, `analysis/enrichment.py`)

Weighted Louvain (`networkx.community.louvain_communities`) at each resolution in
`resolutions_for_stability`, with `n_seeds` runs per resolution. Stability = pairwise adjusted
Rand index (Hubert & Arabie 1985; checked against scikit-learn). The reported partition is the
medoid run, and node stability is the mean co-assignment frequency with the other members of that
node's module.

Enrichment: one-sided hypergeometric test of each module (size ≥ `min_module_size`) for each term
annotating ≥ `min_term_size` background proteins. Background = the species' whole seed network,
because a genome background would trivially "enrich" the GO terms used to select the seeds.
Benjamini–Hochberg FDR is applied across **all** tests of a species. Only FDR ≤ `fdr` is called
significant; overlap statistics elsewhere are explicitly descriptive.

## 3. Evidence profiles and similarity (`analysis/similarity.py`)

Per protein: GO BP+MF terms and KEGG map numbers (STRING; global/overview maps removed), InterPro
entries, EC numbers, PANTHER families and Pfam domains (UniProt), and QuickGO evidence. QuickGO
evidence is summarised as a literature level (1 = experimental code, 0.5 = phylogenetic,
computational or curator, 0.25 = electronic only), experimental PubMed references, a perturbation
flag (IMP/IGI/HMP/HGI) and NOT-qualified annotations, which are kept as contradictory evidence.
Curated findings add their PubMed reference and set the literature level to 1 when they support
the role. Each protein gets a **perturbation score** from all its perturbation evidence (GO codes
and curated rows): 1 if all of it supports the role, 0.5 if it is mixed, 0 if it only contradicts,
and missing if there is none. A pair's perturbation component is the mean over the proteins that
*have* perturbation evidence. An earlier rule gave 0.5 when only one protein had evidence, which
implicitly scored the other protein's missing evidence as negative.

Set similarity is IC-weighted, with IC(t) = −ln p(t) from both networks' annotation corpus and
shared information S = Σ IC(A∩B):

- `ochiai` (default): S / √(Σ IC(A) · Σ IC(B))
- `simgic`: S / Σ IC(A∪B) (Pesquita et al. 2008)
- `min`: S / min(Σ IC(A), Σ IC(B))

Ochiai is the default because annotation depth is highly unequal. Union normalisation gave the
*D. radiodurans* sodA – human SOD2 ortholog pair a GO similarity of 0.07, because SOD2 carries
hundreds of propagated terms. All three measures are run as a sensitivity analysis.

EC similarity: 1.0 for a shared complete EC number, 0.5 for a shared sub-subclass (three levels),
0.0 if both proteins have EC numbers but share no sub-subclass, and missing otherwise.

Homology tiers (first match wins):

| tier | score |
|---|---|
| reciprocal best hit (bitscore ≥ `min_homology_bitscore`) | 1.0 |
| shared KEGG Orthology (KO) group | 0.9 |
| shared PANTHER family | 0.8 |
| one-way best hit | 0.6 |
| shared Pfam domain | 0.5 |
| shared InterPro entry | 0.4 |
| none detected | 0.0 |

Every protein was searched, so "none detected" is an observation. A reciprocal best hit whose two
proteins carry different KO groups is reported as contradictory evidence. InterPro is included so that
"analogue" requires the absence of *any* shared signature (e.g. thioredoxin-fold peroxiredoxin
subfamilies with different Pfam models).

Role similarity is 1 − mean |Δ| over within-species percentiles of degree, betweenness and
closeness, plus clustering and participation coefficient. It is missing when either protein is
isolated, because isolation at a threshold partly reflects coverage. Module similarity compares the
"core terms" (carried by ≥ 50 % of members) of the two proteins' non-singleton modules.

## 4. Cross-species comparison

- **Annotation overlap:** term vocabularies per category (KEGG compared by map number); descriptive only.
- **Topology:** side-by-side statistics plus the size-matched context.
- **Module correspondence:** core-term similarity of every pair of non-singleton modules.
- **Alignment** (`analysis/alignment.py`): IsoRank (Singh, Xu & Berger 2008) with α ∈ `alphas`
  (α = 1 is topology only; α < 1 adds a homology prior from STRING bitscores), followed by greedy
  one-to-one matching. Edge conservation is compared with random one-to-one mappings (a weak
  baseline) and with **degree-preserving rewired copies of the source graph aligned the same way**,
  which tests whether the observed topology aligns better than chance graphs with the same degrees.

## 5. Candidate ranking (`analysis/ranking.py`)

All D. radiodurans × human node pairs are scored on eight components (`ec`, `go`, `pathway`,
`literature`, `perturbation`, `module`, `role`, `homology`).

```
score = Σ w_i s_i / Σ w_i      over components that are not missing
evidence_coverage = Σ w_i      over the same components
```

Missing evidence therefore never lowers a score. Pairs with coverage below
`min_evidence_coverage` are listed but not ranked. The weights and their rationale are in the config.

Categories are assigned by fixed rules, independently of the score:

| category | rule |
|---|---|
| `ortholog_rbh` | reciprocal best hit |
| `ortholog_ko` | shared KEGG KO group, not a reciprocal best hit |
| `homolog` | any other homology tier |
| `candidate_analogue` | no homology, EC ≥ 0.5, and at least one of: pathway ≥ threshold, GO ≥ threshold, or experimental evidence for both proteins |
| `shared_function_annotation` | no homology, EC ≥ 0.5 or GO ≥ threshold |
| `similar_pathway_role` | no homology, pathway ≥ threshold only |
| `speculative` | no homology, only role or module ≥ 0.5 |
| `unsupported` | none of the above |

Every pair carries its supporting evidence, contradictory evidence, references (experimental GO
annotations), uncertainty notes and a category-specific proposed validation.

Sensitivity analyses:

- Dirichlet perturbation of the weights: Spearman correlation, top-k overlap and 5–95 % rank interval per pair.
- Leave-one-component-out.
- Similarity measure (ochiai, simgic, min).
- STRING threshold: modules and roles recomputed at every sensitivity threshold.
- Random seed: separate runs with `--seed`, compared with `scripts/compare_phase2_runs.py`.

### Determinism

Python randomises string hashing per process, so set iteration order varies between runs. All
IC sums therefore run over sorted terms. Scores are rounded to 12 significant digits before
ranking, and ties are broken by protein identifiers with a stable sort. Tables are written with
12 significant digits. Outputs are byte-identical across `PYTHONHASHSEED` values (checked for
0–3 on the reference run).

## 5b. Coverage variant (`analysis/coverage.py`)

The KO-projected *D. radiodurans* network is built as follows:

1. Take every gene that KEGG assigns to a KO group carried by any human seed protein, together
   with the original seeds.
2. Map these genes with STRING `get_string_ids` (species 243230).
3. Retrieve their functional network at the Phase 1 retrieval score, then apply the Phase 1 edge
   QC at the primary threshold.

The variant is then aligned to the human network, with no prior (α = 1) and with a shared-KO prior
(α = 0.8), and both are tested against the same nulls as the main alignment.

The variant is defined by orthology, so it is used only to test whether the topology and alignment
conclusions depend on coverage. It is never used for candidate ranking. STRING transfers
association evidence between orthologs of different organisms, so edge conservation between
KO-orthologous pairs is partly built into the input scores.

## 6. Embeddings (`analysis/embedding.py`)

NetMF (Qiu et al. 2018) is the closed-form matrix DeepWalk factorises. It is deterministic and
needs only numpy. It is evaluated by held-out-edge link prediction within one species. The test
edges are removed *before* the embedding is computed (enforced by a test), and NetMF is compared
with common neighbours, Jaccard, Adamic–Adar and preferential attachment on the same splits.
Embeddings are never compared across species, because two independently factorised graphs do not
share a coordinate system. Graphs with fewer than `min_edges_for_evaluation` edges are not evaluated.

## 7. Evaluation

No ground truth for cross-species functional analogy exists, and none is invented. As a limited
proxy, the pipeline reports how well each component separates RBH ortholog pairs from other pairs
(AUC). An AUC is called evaluable only with ≥ 10 positives. Annotation components carry a leakage
flag, because GO and EC annotations are often transferred between homologs.
