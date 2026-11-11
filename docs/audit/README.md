# Audit of the manuscript (2026-11-11)

An audit of `docs/overleaf/paper.tex` against the code, the cached raw data and the run outputs.

1. [Reproducibility report](01_REPRODUCIBILITY_REPORT.md)
2. [Scientific corrections report](02_SCIENTIFIC_CORRECTIONS_REPORT.md)
3. [Reference verification report](03_REFERENCE_VERIFICATION_REPORT.md)
4. [Methodology audit](04_METHODOLOGY_AUDIT.md)
5. [AI-use disclosure](05_AI_USE_DISCLOSURE.md)
6. [Unresolved scientific limitations](06_UNRESOLVED_LIMITATIONS.md)

Supporting files: `independent_audit_output.txt` and `independent_audit_results.json`, produced by
`python scripts/independent_audit.py --json docs/audit/independent_audit_results.json`.

**Summary.** 85 numerical claims were recomputed independently from the raw data: 70 match, 9 agree
within the expected range for fresh random draws, and 6 differ. All 6 trace to two pipeline bugs,
which affect the IntAct confirmation counts and the number of ranked pairs but no reported rank,
p-value, network statistic or candidate. All 43 references exist and are cited correctly; seven
background claims need citations. These reports and the audit script were prepared by Claude; the
researcher's own verification is the next step.
