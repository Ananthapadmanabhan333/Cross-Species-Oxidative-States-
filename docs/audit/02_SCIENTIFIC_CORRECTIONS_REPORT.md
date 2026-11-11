# 2. Scientific corrections report

Status key: **open** = found, not yet corrected in code or manuscript; **done** = corrected.

## A. Errors in the results

### C1. IntAct confirmation under-counted (open)

- **Evidence.** IntAct writes a multi-protein complex as several binary rows that share one
  interaction ID. `intact_within_set` (`src/data_collection/evidence_sources.py`) de-duplicates by
  interaction ID alone, which keeps one protein pair per complex. In the oxidative-stress human
  data, 222 complexes carry 1,499 further pairs that were discarded.
- **Effect.** Table I ("edges in IntAct") and the Physical interactions paragraph are wrong under
  either definition of a physical interaction:

  | Human STRING edges confirmed | Manuscript | All physical types incl. co-complex association (MI:0914) | Binary physical types only |
  |---|---|---|---|
  | Oxidative stress | 275 (12.4 %) | 522 (23.6 %) | 217 (9.8 %) |
  | DNA repair | 431 (8.2 %) | 744 (14.1 %) | 272 (5.1 %) |
  | Physical subnetwork, oxidative | 21.1 % | 40.5 % | 16.7 % |
  | Physical subnetwork, DNA repair | 16.5 % | 28.0 % | 10.4 % |

- **Conclusion affected?** No. Most STRING edges still lack a deposited physical interaction, and
  the bacterial result (no records in IntAct) is unchanged.
- **Fix.** De-duplicate on (interaction ID, partner A, partner B); decide whether co-complex
  association counts as a physical interaction; regenerate the runs; add a regression test.

### C2. Pairs at exactly 0.45 evidence coverage excluded (open)

- **Evidence.** 0.15 + 0.10 + 0.10 + 0.10 sums to 0.44999999999999996 in floating point, so pairs
  with coverage of exactly 0.45 fail the `>= 0.45` test. The method excludes only pairs *below* 0.45.
- **Effect.** 4 oxidative and 64 DNA-repair pairs are wrongly unranked. All are `speculative`, with
  scores ≤ 0.31, and would enter at about rank 2,336 and 12,687, below every pair discussed in the
  paper. The counts "4,329 and 20,020 had enough evidence to be ranked" become 4,333 and 20,084.
- **Fix.** Compare coverage with a tolerance; regenerate the runs; add a boundary test.

### C3. Role-similarity proxy: number of positives (open, manuscript only)

The manuscript says AUC 0.46 "on 16 DNA-repair ortholog pairs". The proxy set has 16 pairs, but
only 14 have a role value (`evaluation/proxy_evaluation.json`), so the AUC rests on 14.

### C4. Number of automated tests (done)

The project brief refers to 88 tests. The suite has 91 (three regression tests were added in
v0.3.0); the manuscript already says 91.

### C5. CPU-dependent results (done, v0.3.0)

NetMF normalised rounding-noise rows into random vectors, and IsoRank broke near-ties by rounding
noise. Fixing both changed the DNA-repair topology-only alignment p from 0.05 to 0.08 and the NetMF
AUCs (0.933 → 0.935, 0.877 → 0.879, 0.737 → 0.728). The manuscript reports the corrected values.

### C6. Outdated example in the methods document (open, documentation only)

`docs/PHASE2_METHODS.md` says that discarding common folds turned RecQ DR_1289–WRN into a false
analogue. In the reference data that pair shares a PANTHER family and is a homolog either way. The
manuscript uses the verified example instead (peroxiredoxin Q–GPX2).

## B. Interpretation and writing issues in the manuscript

These are for the authors to resolve; the text has not been changed.

1. **Ortholog vs homolog for DdrA.** The Introduction says DdrA, DdrB and PprA "have no human
   ortholog" (uncited), while Results report DdrA–RAD52 as a homolog via a shared Pfam domain, with
   a citation describing DdrA as a Rad52 homolog. Both can be true, but the text should define the
   difference and cite the Introduction claim.
2. **"Physical interaction" includes co-complex membership.** MI:0914 (association) means two
   proteins were found in the same complex, not that they bind directly. The terminology should
   follow whichever definition C1 adopts.
3. **Uncited background claims.** Unrelated Cu/Zn and Mn/Fe SOD families; Fpg vs HhH-GPD glycosylase
   folds and their shared 8-oxoguanine activity; UvrABC vs XPC/TFIIH; shared ligase domains;
   "bacterial interaction data are sparse and biased". Verified candidate sources are in report 3.
4. **Analyses described but never reported.** Hubs and bottlenecks, PageRank and k-core are listed
   in Methods but not used in Results. Report them or remove them.
5. **Multiple alignment tests.** Nine rewired-null tests are reported (eight in Table II, one in
   the text) without adjustment for
   multiplicity, and the two smallest p-values sit at the resolution limit of 100 samples. The
   Abstract states p = 0.01–0.02 without these qualifications.
6. **Feasibility of the proposed SOD1 experiment.** Copper delivery to human SOD1 depends on its
   chaperone CCS (Culotta et al. 1997, doi:10.1074/jbc.272.38.23469); expressing SOD1 alone in
   *D. radiodurans* may yield inactive enzyme. The proposal should address this.
7. **Weights are chosen, not estimated.** The component weights (0.15 / 0.10) are a priori choices.
   The text says so implicitly; it should say so explicitly and point to the weight-sensitivity
   result (top-20 overlap 0.80 and 0.70).
