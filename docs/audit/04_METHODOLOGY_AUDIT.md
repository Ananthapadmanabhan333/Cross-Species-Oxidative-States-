# 4. Methodology audit

For each step: what is done, how it was verified in this audit, and the weaknesses a reviewer is
likely to raise.

| Step | What is done | Verified | Weaknesses |
|---|---|---|---|
| Seed proteins | UniProtKB entries with process GO terms incl. descendants; human reviewed only, *D. radiodurans* all entries | counts rebuilt from raw UniProt responses | seed sets inherit GO annotation depth, which differs strongly between species; reviewed-only vs all-entries is asymmetric by design |
| Mapping | STRING `get_string_ids`, taxon checked (1299 → 243230) | rebuilt from raw responses; taxa confirmed | unmapped seeds (4 bacterial and 6 human for DNA repair) are reported, not recovered |
| Networks | STRING functional edges, combined score ≥ 0.7 | rebuilt from raw responses; all counts match | STRING transfers evidence between orthologs, so cross-species similarity is partly built into the input; threshold sensitivity is reported but 0.7 is a convention |
| Physical layer | STRING physical subnetwork flag; IntAct via PSICQUIC | bacterial results match; human counts were wrong (C1, fixed in v0.4.0) | v0.4.0 reports both definitions: all physical types, and binary types only (co-complex association excluded) |
| Topology | degree, betweenness, closeness, clustering, PageRank, k-core; hubs and bottlenecks | not re-derived node by node | several metrics are computed but not used in the Results |
| Modules | weighted Louvain, resolution 1.0, 50 seeds, medoid, ARI | modularity of the reported partition matches; ARI reproduced with fresh seeds | a single resolution; human ARI 0.64–0.74 means module boundaries are only moderately stable |
| Enrichment | hypergeometric vs the seed network, BH across all tests | the reported *D. radiodurans* Peroxisome result checked in the table | background choice is conservative; terms used to select the seeds cannot be enriched |
| Alignment | IsoRank (α = 1, 0.8, 0.5), greedy one-to-one matching, random and degree-preserving nulls | IsoRank re-implemented; all conserved-edge counts match; p-values are consistent integer counts | 100 rewired samples limit p to ≥ 0.0099; nine reported tests without multiplicity adjustment; conserved-edge count is not normalised for the density of the human network |
| Ranking | weighted mean over available components, coverage ≥ 0.45, rule-based categories | scores recomputed for all 28,264 pairs; categories match; eligibility bug found (C2, fixed in v0.4.0) | weights are a priori, not estimated; EC similarity depends on annotation quality; categories use fixed thresholds (0.30, 0.5) |
| Curated evidence | 21 findings for 12 bacterial proteins from PubMed abstracts | rows carry PMIDs and DOIs | abstract-level only; targeted search, not systematic; no human mutant curation |
| Embeddings | NetMF, 16 dimensions, window 5, link prediction vs neighbourhood baselines | re-implemented; AUCs match | one dimension and window setting; not used for ranking, so it is a negative baseline only |
| Evaluation | ortholog proxy, AUC evaluable with ≥ 10 positives | values checked in the output table | proxy leakage for annotation components; no ground truth for analogy |
| Robustness | Dirichlet weights, leave-one-out, similarity measure, threshold, seeds | seed-sensitivity tables checked | top-20 lists change noticeably with weights (median overlap 0.80 and 0.70) |
| Tests | 93 unit and integration tests, HTTP mocked (two regression tests for C1 and C2 added in v0.4.0) | all pass on Linux and Windows | they check that code does what it was written to do on synthetic inputs; they do not validate biology, and bugs C1 and C2 were not covered until v0.4.0 |
