"""Figures, all drawn from the computed tables of the current run (matplotlib, Agg backend)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .topology import degree_distribution  # noqa: E402

LABELS = {"deinococcus_radiodurans": "D. radiodurans", "homo_sapiens": "H. sapiens"}
COLORS = {"deinococcus_radiodurans": "#c2410c", "homo_sapiens": "#1d4ed8"}


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, metadata={"Software": "extremonet phase 2"})
    plt.close(fig)


def degree_ccdf(base: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    for k, res in base.items():
        dd = degree_distribution(res["graph"])
        dd = dd[dd["degree"] > 0]
        ax.step(dd["degree"], dd["ccdf"], where="post", label=f"{LABELS[k]} (n={res['graph'].number_of_nodes()})",
                color=COLORS[k])
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlabel("degree k (connected nodes)"), ax.set_ylabel("P(K >= k)")
    ax.set_title("Degree distribution (CCDF), threshold 0.7")
    ax.legend(frameon=False)
    _save(fig, path)


def threshold_curves(table: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    for k, grp in table.groupby("species"):
        nodes = grp["connected_nodes"] + grp["isolated_nodes"]
        axes[0].plot(grp["threshold"], grp["edges"], "o-", color=COLORS[k], label=LABELS[k])
        axes[1].plot(grp["threshold"], grp["connected_nodes"] / nodes, "o-", color=COLORS[k])
        axes[2].plot(grp["threshold"], grp["louvain_ari_mean"], "o-", color=COLORS[k])
    axes[0].set_yscale("log"), axes[0].set_ylabel("edges"), axes[0].legend(frameon=False)
    axes[1].set_ylabel("fraction of seeds connected"), axes[1].set_ylim(0, 1.02)
    axes[2].set_ylabel("Louvain mean pairwise ARI"), axes[2].set_ylim(0, 1.02)
    for ax in axes:
        ax.set_xlabel("STRING combined-score threshold")
    _save(fig, path)


def size_matched(null: pd.DataFrame, dr_stats: dict, path: Path) -> None:
    metrics = [("edges", "edges"), ("connected_nodes", "connected_nodes"), ("largest_component_nodes", "largest_component_nodes")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    for ax, (col, key) in zip(axes, metrics):
        vals = null[col]
        bins = np.arange(vals.min() - 0.5, vals.max() + 1.5) if vals.max() - vals.min() < 60 else 30
        ax.hist(vals, bins=bins, color=COLORS["homo_sapiens"], alpha=0.6, label="random human subsets")
        ax.axvline(dr_stats[key], color=COLORS["deinococcus_radiodurans"], lw=2, label="D. radiodurans")
        ax.set_xlabel(col.replace("_", " ")), ax.set_ylabel("samples")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Size-matched context: {len(null)} random subsets of the larger network, size {dr_stats['nodes']}",
                 fontsize=10)
    _save(fig, path)


def louvain_stability(base: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 3.6))
    for k, res in base.items():
        s = pd.DataFrame([v["summary"] for v in res["stability"].values()]).sort_values("resolution")
        ax.plot(s["resolution"], s["ari_mean"], "o-", color=COLORS[k], label=LABELS[k])
        ax.fill_between(s["resolution"], s["ari_min"], s["ari_mean"], color=COLORS[k], alpha=0.15)
    ax.set_xlabel("Louvain resolution"), ax.set_ylabel("pairwise ARI across seeds (mean, min)")
    ax.set_ylim(0, 1.02), ax.legend(frameon=False)
    _save(fig, path)


def ranking_intervals(ranked: pd.DataFrame, path: Path, top: int = 20) -> None:
    data = ranked[ranked["eligible"]].head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(data) + 1.2))
    y = np.arange(len(data))
    ax.hlines(y, data["rank_q05"], data["rank_q95"], color="#6b7280")
    ax.plot(data["rank"], y, "o", color="#111827", label="rank (base weights)")
    ax.set_yticks(y, [f"{a} - {b} [{c}]" for a, b, c in zip(data["dr_name"], data["hs_name"], data["category"])],
                  fontsize=7)
    ax.set_xlabel("rank (line: 5-95% under Dirichlet weight perturbation)")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, path)


def link_prediction(emb_summary: dict, path: Path) -> bool:
    evaluated = {k: v for k, v in emb_summary.items() if v.get("evaluated")}
    if not evaluated:
        return False
    k, v = next(iter(evaluated.items()))
    s = pd.DataFrame(v["summary"])
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    ax.barh(s["method"], s["mean"], xerr=s["std"], color=["#7c3aed" if m == "netmf_cosine" else "#9ca3af" for m in s["method"]])
    ax.axvline(0.5, color="black", lw=0.8, ls="--")
    ax.set_xlim(0, 1), ax.set_xlabel(f"held-out edge AUC ({LABELS[k]}, mean +/- sd over splits)")
    _save(fig, path)
    return True


def dr_network(base: dict, path: Path) -> None:
    g = base["deinococcus_radiodurans"]["graph"]
    part = base["deinococcus_radiodurans"]["partition"]
    fig, ax = plt.subplots(figsize=(6, 5))
    connected = [n for n in g if g.degree(n) > 0]
    isolated = sorted(n for n in g if g.degree(n) == 0)
    pos = nx.spring_layout(g.subgraph(connected), seed=1, weight="weight", k=1.2) if connected else {}
    for i, n in enumerate(isolated):  # isolates in a row below the connected part
        pos[n] = np.array([-1.0 + 2.0 * i / max(len(isolated) - 1, 1), -1.6])
    cmap = plt.get_cmap("tab20")
    nx.draw_networkx_edges(g, pos, ax=ax, width=[3 * d["weight"] for *_, d in g.edges(data=True)], alpha=0.6)
    nx.draw_networkx_nodes(g, pos, ax=ax, node_color=[cmap(part[n] % 20) for n in g], node_size=260)
    nx.draw_networkx_labels(g, pos, ax=ax, labels={n: g.nodes[n].get("string_preferred_name", n) for n in g}, font_size=7)
    ax.set_title("D. radiodurans seed network at 0.7 (colour = Louvain module; bottom row = isolated)", fontsize=9)
    ax.axis("off")
    _save(fig, path)


def make_all(out: Path, base: dict, null: pd.DataFrame, dr_stats: dict, threshold_table: pd.DataFrame,
             ranked: pd.DataFrame, sens: dict, emb_summary: dict, fig_dir: Path) -> list[str]:
    made = []
    for name, fn in (("degree_ccdf.png", lambda p: degree_ccdf(base, p)),
                     ("threshold_sensitivity.png", lambda p: threshold_curves(threshold_table, p)),
                     ("size_matched_context.png", lambda p: size_matched(null, dr_stats, p)),
                     ("louvain_stability.png", lambda p: louvain_stability(base, p)),
                     ("ranking_weight_sensitivity.png", lambda p: ranking_intervals(ranked, p)),
                     ("dr_network_modules.png", lambda p: dr_network(base, p))):
        fn(fig_dir / name)
        made.append(name)
    if link_prediction(emb_summary, fig_dir / "link_prediction_auc.png"):
        made.append("link_prediction_auc.png")
    return made
