# 1. Reproducibility report

Audit date: 2026-11-11. Manuscript: `docs/overleaf/paper.tex` (commit `271f168`).
Code and data: release `v0.3.0` plus later documentation commits.

## What was checked, and how

| Check | Method | Result |
|---|---|---|
| Integrity of inputs and outputs | `scripts/verify_checksums.py --all`: SHA-256 of every file in 9 run directories and all 225 cached API responses | 0 failures |
| Re-execution on a different CPU | Both Phase 2 analyses re-run with OpenBLAS forced to Haswell and SkylakeX kernels | 49/49 and 50/50 result files identical |
| Re-execution on a fresh machine | GitHub Actions run `41631754742` (Windows runner) re-runs both analyses offline and compares them with the reference runs | identical |
| Unit and integration tests | `pytest`, Windows and Linux in CI | 91 passed |
| Independent recomputation | `scripts/independent_audit.py` (see below) | 70 match, 9 close, 6 differ |

## Independent recomputation

`scripts/independent_audit.py` does not import the pipeline. It reads the raw cached responses
named in each run's manifest, rebuilds the networks, and re-implements Louvain stability, IsoRank,
NetMF link prediction, the score formula and the category rules from their written definitions.
Full output: `independent_audit_output.txt` and `independent_audit_results.json`.

| Area | Claims checked | Outcome |
|---|---|---|
| STRING version (12.5), seed counts, mapping, edges, isolated nodes, largest component, physical-subnetwork edges | 26 | all match |
| Random-subset connectivity (3.5 / 77 edges; 97.7 % / 95.3 %) | 4 | close (fresh random draws give 3.6 / 77.4 and 97.6 % / 95.2 %) |
| Louvain ARI and modularity | 7 | modularity matches; ARI close (0.75, 0.80, 0.67 with fresh seeds vs 0.74, 0.76, 0.64) |
| IsoRank conserved edges and rewired-null p-values | 12 | all match; every p-value corresponds to an integer count of null samples |
| NetMF and Adamic–Adar link-prediction AUC | 6 | all match to three decimals |
| Pair scores, eligibility, categories, ranks of the 16 pairs in Table III | 22 | scores and categories match; **eligibility differs for 68 pairs** (correction C2) |
| IntAct confirmation | 8 | bacterial values match; **human values differ** (correction C1) |

"Close" means the value was recomputed with different random draws (fresh Louvain seeds or fresh
random subsets), so exact equality is not expected; the agreement is within the stated tolerance.

## Limits of this reproducibility

- **The data are a snapshot.** STRING 12.5 and UniProt release 2026_03, retrieved 2026-10-25 and
  2026-11-02. An online re-run on a later date may return different records; the committed cache is
  what makes the results reproducible.
- **The independent checker was also written by Claude.** It is an independent *implementation*,
  not an independent *person*. A human re-derivation of the key numbers (Part 3) is still needed.
- **Exact reproduction in CI is checked on Windows only.** Linux runs the tests and checksums but
  not the full re-run.
- **KEGG terms.** Two cached KEGG responses are redistributed in a public repository; KEGG's terms
  restrict redistribution. This needs a decision by the authors.
