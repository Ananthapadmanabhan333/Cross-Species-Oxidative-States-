# 1. Reproducibility report

Audit date: 2026-11-11. Manuscript: `docs/overleaf/paper.tex`. First audited at commit `271f168`
(release v0.3.0); updated for release v0.4.0, which fixes corrections C1 and C2 and regenerates
the reference runs (`20261111T155419Z_860ea5376e14`, `20261111T155841Z_73788e2787df`).

## What was checked, and how

| Check | Method | Result |
|---|---|---|
| Integrity of inputs and outputs | `scripts/verify_checksums.py --all`: SHA-256 of every file in 9 run directories and all 225 cached API responses | 0 failures |
| Re-execution on a different CPU | Both Phase 2 analyses re-run with OpenBLAS forced to Haswell and SkylakeX kernels | 49/49 and 50/50 result files identical |
| Re-execution on a fresh machine | GitHub Actions (Windows runner) re-runs both analyses offline on every push and compares them with the reference runs | identical (v0.3.0 run `41631754742`; v0.4.0 checked on push) |
| Unit and integration tests | `pytest`, Windows and Linux in CI | 93 passed |
| Independent recomputation | `scripts/independent_audit.py` (see below) | v0.4.0: 76 match, 9 close, 0 differ (v0.3.0: 70 / 9 / 6) |

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
| Pair scores, eligibility, categories, ranks of the 16 pairs in Table III | 22 | all match in v0.4.0 (v0.3.0: eligibility differed for 68 pairs, correction C2) |
| IntAct confirmation | 8 | all match in v0.4.0 (v0.3.0: human values differed, correction C1) |

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
- **KEGG data are not redistributed (since v0.4.1).** KEGG's terms restrict redistribution, so the
  two KEGG responses were removed from the repository; their provenance records remain, and
  `scripts/fetch_kegg.py` restores them from KEGG if the SHA-256 still matches. KEGG had already
  changed by 2026-11-11 (dra: 2 lines removed, 1 added; hsa: 8 added), so outside the authors'
  copy the reference runs cannot be reproduced byte for byte; CI then runs the analyses on the
  current KEGG release without the comparison. The files remain in the git history before
  v0.4.1.
