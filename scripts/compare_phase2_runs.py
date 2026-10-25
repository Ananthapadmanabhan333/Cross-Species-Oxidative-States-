"""Compare Phase 2 runs (e.g. made with different --seed values) for robustness reporting.

Usage (from the project root):
    python scripts/compare_phase2_runs.py <reference_run_id> <other_run_id> [<other_run_id> ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("results/phase2")


def load(run_id: str) -> tuple[pd.DataFrame, dict]:
    ranked = pd.read_csv(ROOT / run_id / "candidates" / "candidate_pairs_all.tsv", sep="\t")
    summary = json.loads((ROOT / run_id / "analysis_summary.json").read_text(encoding="utf-8"))
    return ranked, summary


def top(ranked: pd.DataFrame, k: int) -> set[tuple[str, str]]:
    eligible = ranked[ranked["eligible"]].sort_values("rank").head(k)
    return set(zip(eligible["dr_string_id"], eligible["hs_string_id"]))


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    ref_id, others = argv[0], argv[1:]
    ref, ref_sum = load(ref_id)
    rows = []
    for run_id in others:
        other, summ = load(run_id)
        merged = ref.merge(other, on=["dr_string_id", "hs_string_id"], suffixes=("_ref", "_other"))
        row = {"run": run_id, "seed": summ.get("random_seed")}
        for k in (10, 20, 50):
            row[f"top{k}_overlap"] = len(top(ref, k) & top(other, k)) / k
        row["spearman_score"] = merged["score_ref"].corr(merged["score_other"], method="spearman")
        row["category_agreement"] = float((merged["category_ref"] == merged["category_other"]).mean())
        for sp in ("deinococcus_radiodurans", "homo_sapiens"):
            ref_l = {r["resolution"]: r for r in ref_sum["louvain"][sp]}
            oth_l = {r["resolution"]: r for r in summ["louvain"][sp]}
            row[f"{sp}_modules_ref_vs_other"] = f"{ref_l[1.0]['modules_median']} / {oth_l[1.0]['modules_median']}"
            row[f"{sp}_ari_ref_vs_other"] = f"{ref_l[1.0]['ari_mean']:.3f} / {oth_l[1.0]['ari_mean']:.3f}"
        rows.append(row)
    table = pd.DataFrame(rows)
    print(table.to_string(index=False))
    out = ROOT / f"seed_sensitivity_{ref_id}.tsv"  # outside the checksummed run directories
    table.to_csv(out, sep="\t", index=False, float_format="%.12g")
    print(f"\nwritten: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
