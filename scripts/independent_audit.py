"""Independent audit of the manuscript's numbers, rebuilt from the raw cached API responses.

This script deliberately does NOT import the `extremonet` pipeline. It re-parses the raw STRING,
UniProt and IntAct responses named in each run's manifest, rebuilds the networks, and re-implements
the published definitions (Louvain stability, IsoRank, NetMF link prediction, the pair score and the
category rules) so that a pipeline bug cannot confirm itself. Each manuscript claim is compared with
the recomputed value and reported as MATCH, CLOSE (agrees within the stated tolerance; for stochastic
re-implementations with different random draws) or DIFFERS.

Usage (from the project root):  python scripts/independent_audit.py [--json out.json]
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path.cwd()
RESULTS = []

PROCESSES = {
    "oxidative": {"phase1": "20261025T174736Z_5bbc7e723df4", "phase2": "results/phase2/20261111T155419Z_860ea5376e14"},
    "dna_repair": {"phase1": "20261102T161839Z_7128c0f5ca4e",
                   "phase2": "results/dna_repair/phase2/20261111T155841Z_73788e2787df"},
}
SPECIES = {"deinococcus_radiodurans": 243230, "homo_sapiens": 9606}
THRESHOLD = 0.7
MASTER_SEED = 20261026
PHYSICAL_MI = {"MI:0915", "MI:0407", "MI:0914", "MI:0195", "MI:0408"}
WEIGHTS = {"ec": .15, "go": .15, "pathway": .15, "literature": .15, "perturbation": .10, "module": .10,
           "role": .10, "homology": .10}


def record(process, claim, manuscript, recomputed, status=None, tol=0.0, note=""):
    if status is None:
        if isinstance(manuscript, (int, float)) and isinstance(recomputed, (int, float)):
            diff = abs(manuscript - recomputed)
            status = "MATCH" if diff == 0 else ("CLOSE" if diff <= tol else "DIFFERS")
        else:
            status = "MATCH" if manuscript == recomputed else "DIFFERS"
    RESULTS.append({"process": process, "claim": claim, "manuscript": manuscript, "recomputed": recomputed,
                    "status": status, "note": note})


def raw_files(manifest_path, endpoint):
    m = json.loads((ROOT / manifest_path).read_text(encoding="utf-8"))
    return [ROOT / r["data_file"].replace("\\", "/") for r in m["raw_inputs"] if r["endpoint"] == endpoint]


def verify_checksum(path):
    meta = json.loads(path.with_name(path.name.rsplit(".", 1)[0] + ".meta.json").read_text(encoding="utf-8"))
    ok = hashlib.sha256(path.read_bytes()).hexdigest() == meta["sha256"]
    return ok, meta


# ---------------------------------------------------------------- network construction (Phase 1)
def build_species(phase1, species):
    manifest = f"data/manifests/{phase1}.json"
    taxon = SPECIES[species]
    (seed_file,) = raw_files(manifest, f"seeds_{species}")
    seeds = list(csv.DictReader(seed_file.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    accessions = [s["Entry"] for s in seeds]
    mapping = {}
    for f in raw_files(manifest, f"get_string_ids_{species}") + raw_files(manifest, f"get_string_ids_{species}_single"):
        for row in json.loads(f.read_text(encoding="utf-8")):
            if row["ncbiTaxonId"] == taxon and row["queryItem"] in accessions:
                mapping[row["queryItem"]] = (row["stringId"], row["preferredName"])
    nodes = {sid for sid, _ in mapping.values()}
    names = {sid: name for sid, name in mapping.values()}
    g = nx.Graph()
    g.add_nodes_from(sorted(nodes))
    (net_file,) = raw_files(manifest, f"network_functional_{species}")
    best = {}
    for e in json.loads(net_file.read_text(encoding="utf-8")):
        a, b = e["stringId_A"], e["stringId_B"]
        if a == b or a not in nodes or b not in nodes:
            continue
        key = tuple(sorted((a, b)))
        best[key] = max(best.get(key, 0.0), float(e["score"]))
    g.add_weighted_edges_from((a, b, s) for (a, b), s in best.items() if s >= THRESHOLD)
    (phys_file,) = raw_files(manifest, f"network_physical_{species}")
    phys = {tuple(sorted((e["stringId_A"], e["stringId_B"]))) for e in json.loads(phys_file.read_text(encoding="utf-8"))
            if e["stringId_A"] in nodes and e["stringId_B"] in nodes and float(e["score"]) >= THRESHOLD}
    phys &= {tuple(sorted(e)) for e in g.edges}
    acc_to_sid = {acc: sid for acc, (sid, _) in mapping.items()}
    return {"graph": g, "seeds": len(seeds), "mapped_accessions": len(mapping), "names": names,
            "acc_to_sid": acc_to_sid, "physical_edges": phys, "version_ok": None}


def network_checks(process, nets, claims):
    for sp, net in nets.items():
        g = net["graph"]
        c = claims[sp]
        tag = "D.r." if sp.startswith("dein") else "H.s."
        record(process, f"{tag} seed entries", c["seeds"], net["seeds"])
        record(process, f"{tag} mapped proteins", c["mapped"], g.number_of_nodes())
        record(process, f"{tag} edges (score >= 0.7)", c["edges"], g.number_of_edges())
        record(process, f"{tag} isolated proteins", c["isolated"], sum(1 for n in g if g.degree(n) == 0))
        record(process, f"{tag} largest component", c["largest"], len(max(nx.connected_components(g), key=len)))
        record(process, f"{tag} edges also in STRING physical subnetwork", c["physical"], len(net["physical_edges"]))


def random_subset_check(process, nets, claims, samples=5000):
    hs, dr = nets["homo_sapiens"]["graph"], nets["deinococcus_radiodurans"]["graph"]
    rng = np.random.default_rng(1)
    nodes = np.array(sorted(hs.nodes))
    k = dr.number_of_nodes()
    counts = np.array([hs.subgraph(rng.choice(nodes, k, replace=False)).number_of_edges() for _ in range(samples)])
    obs = dr.number_of_edges()
    record(process, f"mean edges of random {k}-protein human subsets", claims["mean"], round(float(counts.mean()), 1),
           tol=claims["mean_tol"], note=f"{samples} fresh random subsets (different draws from the pipeline)")
    record(process, "share of random subsets with fewer edges than D.r. (%)", claims["below_pct"],
           round(100 * float((counts < obs).mean()), 1), tol=1.0)


# ---------------------------------------------------------------- IntAct confirmation
def intact_check(process, nets, phase2, claims):
    manifest = f"{phase2}/manifest.json"
    for sp in SPECIES:
        files = raw_files(manifest, f"psicquic_{sp}")
        net = nets[sp]
        confirmed_pairs = set()
        statuses = set()
        for f in files:
            ok, meta = verify_checksum(f)
            statuses.add((ok, meta.get("http_status")))
            for line in f.read_text(encoding="utf-8").splitlines():
                cols = line.split("\t")
                if len(cols) < 12 or not any(mi in cols[11] for mi in PHYSICAL_MI):
                    continue
                ids = []
                for col in (cols[0], cols[1]):
                    acc = col.split("uniprotkb:")[-1].split("|")[0].split("-")[0] if "uniprotkb:" in col else None
                    ids.append(net["acc_to_sid"].get(acc))
                if None not in ids and ids[0] != ids[1]:
                    confirmed_pairs.add(tuple(sorted(ids)))
        g = net["graph"]
        on_edges = sum(1 for e in g.edges if tuple(sorted(e)) in confirmed_pairs)
        tag = "D.r." if sp.startswith("dein") else "H.s."
        record(process, f"{tag} STRING edges confirmed in IntAct", claims[sp]["edges"], on_edges,
               note=f"{len(files)} IntAct responses; (checksum ok, HTTP status): {sorted(statuses)}")
        phys = net["physical_edges"]
        if phys:
            pct = round(100 * sum(1 for e in phys if e in confirmed_pairs) / len(phys), 1)
            record(process, f"{tag} physical-subnetwork edges confirmed (%)", claims[sp]["phys_pct"], pct)


# ---------------------------------------------------------------- Louvain
def ari(a, b):
    table = Counter(zip(a, b))
    comb = lambda n: n * (n - 1) / 2
    sum_ij = sum(comb(v) for v in table.values())
    sa = sum(comb(v) for v in Counter(a).values())
    sb = sum(comb(v) for v in Counter(b).values())
    expected = sa * sb / comb(len(a))
    maximum = (sa + sb) / 2
    return 1.0 if maximum == expected else (sum_ij - expected) / (maximum - expected)


def louvain_check(process, nets, claims, runs=50):
    for sp, c in claims.items():
        g = nets[sp]["graph"]
        nodes = sorted(g.nodes)
        labels = []
        partitions = []
        for s in range(runs):
            comms = nx.community.louvain_communities(g, weight="weight", resolution=1.0, seed=1000 + s)
            label = {n: i for i, com in enumerate(comms) for n in com}
            labels.append([label[n] for n in nodes])
            partitions.append(comms)
        aris = [ari(labels[i], labels[j]) for i, j in combinations(range(runs), 2)]
        tag = "D.r." if sp.startswith("dein") else "H.s."
        record(process, f"{tag} Louvain mean pairwise ARI", c["ari"], round(float(np.mean(aris)), 2), tol=0.05,
               note="50 fresh seeds (pipeline uses different derived seeds)")
        if "modularity" in c:
            mods = [nx.community.modularity(g, p, weight="weight") for p in partitions]
            record(process, f"{tag} Louvain modularity (median of fresh runs)", c["modularity"],
                   round(float(np.median(mods)), 2), tol=0.02)


def pipeline_partition_modularity(process, nets, phase2, claims):
    for sp, value in claims.items():
        g = nets[sp]["graph"]
        names = {v: k for k, v in nets[sp]["names"].items()}
        rows = list(csv.DictReader((ROOT / phase2 / sp / "modules.tsv").read_text(encoding="utf-8").splitlines(),
                                   delimiter="\t"))
        comms = [{names[m] for m in r["members"].split(";") if m} for r in rows]
        covered = set().union(*comms)
        comms += [{n} for n in g if n not in covered]
        tag = "D.r." if sp.startswith("dein") else "H.s."
        record(process, f"{tag} modularity of the pipeline's reported partition", value,
               round(nx.community.modularity(g, comms, weight="weight"), 2),
               note="networkx modularity on the independently rebuilt graph")


# ---------------------------------------------------------------- IsoRank
def column_normalised(g, nodes):
    a = nx.to_numpy_array(g, nodelist=nodes, weight="weight")
    col = a.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(col > 0, a / col, 0.0)


def isorank_conserved(g1, g2, prior, alpha):
    n1, n2 = sorted(g1), sorted(g2)
    a1, a2 = column_normalised(g1, n1), column_normalised(g2, n2)
    e = np.full((len(n1), len(n2)), 1.0 / (len(n1) * len(n2))) if prior is None else prior / prior.sum()
    r = e.copy()
    for _ in range(500):
        nxt = alpha * (a1 @ r @ a2.T) + (1 - alpha) * e
        nxt /= nxt.sum()
        done = np.abs(nxt - r).sum() < 1e-9
        r = nxt
        if done:
            break
    key = np.round(r / r.max(), 10)
    order = np.dstack(np.unravel_index(np.argsort(-key, axis=None, kind="stable"), r.shape))[0]
    used1, used2, mapping = set(), set(), {}
    for i, j in order:
        if i in used1 or j in used2:
            continue
        used1.add(i), used2.add(j)
        mapping[n1[i]] = n2[j]
        if len(used1) == len(n1):
            break
    return sum(1 for u, v in g1.edges if g2.has_edge(mapping[u], mapping[v]))


def homology_prior(phase2, g1, g2):
    manifest = f"{phase2}/manifest.json"
    n1, n2 = sorted(g1), sorted(g2)
    i1, i2 = {n: i for i, n in enumerate(n1)}, {n: i for i, n in enumerate(n2)}
    prior = np.zeros((len(n1), len(n2)))
    for endpoint, flip in (("homology_best_deinococcus_radiodurans_to_9606", False),
                           ("homology_best_homo_sapiens_to_243230", True)):
        for f in raw_files(manifest, endpoint):
            for h in json.loads(f.read_text(encoding="utf-8")):
                a, b = (h["stringId_B"], h["stringId_A"]) if flip else (h["stringId_A"], h["stringId_B"])
                if a in i1 and b in i2:
                    prior[i1[a], i2[b]] = max(prior[i1[a], i2[b]], float(h["bitscore"]))
    return prior


def alignment_check(process, nets, phase2, claims):
    g1, g2 = nets["deinococcus_radiodurans"]["graph"], nets["homo_sapiens"]["graph"]
    prior = homology_prior(phase2, g1, g2)
    summary = {float(r["alpha"]): r for r in csv.DictReader(
        (ROOT / phase2 / "cross_species" / "alignment_summary.tsv").read_text(encoding="utf-8").splitlines(),
        delimiter="\t")}
    for alpha, (edges, p) in claims.items():
        conserved = isorank_conserved(g1, g2, None if alpha == 1.0 else prior, alpha)
        record(process, f"IsoRank conserved edges, alpha={alpha}", edges, conserved,
               note="independent IsoRank + greedy matching on the rebuilt graphs")
        row = summary[alpha]
        n = int(row["rewired_null_samples"])
        p_table = float(row["p_vs_rewired_null"])
        k = round(p_table * (1 + n)) - 1
        consistent = abs((1 + k) / (1 + n) - p_table) < 1e-9
        record(process, f"rewired-null p, alpha={alpha} (as reported)", p, round(p_table, 2),
               tol=0.005, note=f"p = (1 + {k}) / (1 + {n}); integer count consistent: {consistent}")


# ---------------------------------------------------------------- NetMF link prediction
def derived_seed(name):
    return int(hashlib.sha256(f"{MASTER_SEED}:{name}".encode()).hexdigest(), 16) % (2 ** 31)


def netmf(g, dim=16, window=5):
    nodes = sorted(g)
    a = nx.to_numpy_array(g, nodelist=nodes, weight="weight")
    deg = a.sum(axis=1)
    inv = np.where(deg > 0, 1.0 / np.where(deg > 0, deg, 1), 0.0)
    p = inv[:, None] * a
    acc, power = np.zeros_like(a), np.eye(len(nodes))
    for _ in range(window):
        power = power @ p
        acc += power
    m = (deg.sum() / window) * acc * inv[None, :]
    u, s, _ = np.linalg.svd(np.log(np.maximum(m, 1.0)))
    emb = u[:, :dim] * np.sqrt(s[:dim])
    norm = np.linalg.norm(emb, axis=1)
    keep = norm > 1e-9 * norm.max()
    unit = np.zeros_like(emb)
    unit[keep] = emb[keep] / norm[keep, None]
    return unit, {n: i for i, n in enumerate(nodes)}


def auc(pos, neg):
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv))
    sorted_v = allv[order]
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def link_prediction_check(process, nets, claims):
    for sp, c in claims.items():
        g = nets[sp]["graph"]
        rng = np.random.default_rng(derived_seed(f"embed:{sp}"))
        edges = sorted(tuple(sorted(e)) for e in g.edges)
        nodes = sorted(g)
        res = defaultdict(list)
        for _ in range(5):
            perm = rng.permutation(len(edges))
            n_test = max(1, int(round(0.1 * len(edges))))
            test = [edges[i] for i in perm[:n_test]]
            train = g.copy()
            train.remove_edges_from(test)
            neg = set()
            while len(neg) < n_test:
                u, v = rng.choice(nodes, size=2, replace=False)
                pair = tuple(sorted((u, v)))
                if not g.has_edge(*pair):
                    neg.add(pair)
            neg = sorted(neg)
            unit, idx = netmf(train)
            cos = lambda ps: np.round(np.array([unit[idx[u]] @ unit[idx[v]] for u, v in ps]), 10)
            aa = lambda ps: np.array([sum(1 / math.log(train.degree(w)) for w in nx.common_neighbors(train, u, v)
                                          if train.degree(w) > 1) for u, v in ps])
            res["netmf"].append(auc(cos(test), cos(neg)))
            res["adamic_adar"].append(auc(aa(test), aa(neg)))
        tag = "D.r." if sp.startswith("dein") else "H.s."
        record(process, f"{tag} NetMF link-prediction AUC (mean of 5 splits)", c["netmf"],
               round(float(np.mean(res["netmf"])), 3), tol=0.0015,
               note="independent NetMF and AUC; splits drawn with the pipeline's documented seed rule")
        record(process, f"{tag} Adamic-Adar AUC (mean of 5 splits)", c["aa"],
               round(float(np.mean(res["adamic_adar"])), 3), tol=0.0015)


# ---------------------------------------------------------------- ranking arithmetic and categories
def ranking_check(process, phase2, examples):
    rows = list(csv.DictReader((ROOT / phase2 / "candidates" / "candidate_pairs_all.tsv")
                               .read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    worst, eligible_mismatch, category_mismatch = 0.0, 0, 0
    f = lambda x: None if x in ("", "nan", None) else float(x)
    for r in rows:
        comps = {k: f(r[k]) for k in WEIGHTS}
        avail = {k: v for k, v in comps.items() if v is not None and not math.isnan(v)}
        cov = sum(WEIGHTS[k] for k in avail)
        score = sum(WEIGHTS[k] * v for k, v in avail.items()) / cov if cov else float("nan")
        if r["score"] not in ("", "nan"):
            worst = max(worst, abs(score - float(r["score"])))
        if (cov >= 0.45 - 1e-12) != (r["eligible"] == "True"):
            eligible_mismatch += 1
        category_mismatch += category(r, comps) != r["category"]
    record(process, "pair scores recomputed from components (max abs. difference)", 0.0, float(f"{worst:.1e}"),
           tol=1e-9, note=f"{len(rows)} pairs")
    record(process, "eligibility (coverage >= 0.45) disagreements", 0, eligible_mismatch)
    record(process, "category rule disagreements", 0, category_mismatch)
    by_name = {(r["dr_name"], r["hs_name"]): r for r in rows}
    for (d, h), (rank, cat) in examples.items():
        r = by_name[(d, h)]
        record(process, f"rank and category of {d}-{h}", f"{rank} {cat}", f"{r['rank']} {r['category']}")


def category(r, c):
    tier = r["homology_tier"]
    if tier == "reciprocal_best_hit":
        return "ortholog_rbh"
    if tier == "shared_kegg_ko":
        return "ortholog_ko"
    if tier in ("shared_panther_family", "one_way_best_hit", "shared_pfam_domain", "shared_interpro_entry"):
        return "homolog"
    if tier == "shared_common_fold":
        return "distant_homolog"
    ec, go, path = (c["ec"] or 0), (c["go"] or 0), (c["pathway"] or 0)
    if ec >= 0.5 and (path >= 0.30 or go >= 0.30 or c["literature"] == 1.0):
        return "candidate_analogue"
    if ec >= 0.5 or go >= 0.30:
        return "shared_function_annotation"
    if path >= 0.30:
        return "similar_pathway_role"
    if (c["role"] or 0) >= 0.5 or (c["module"] or 0) >= 0.5:
        return "speculative"
    return "unsupported"


# ---------------------------------------------------------------- claims from the manuscript
CLAIMS = {
    "oxidative": {
        "net": {"deinococcus_radiodurans": dict(seeds=16, mapped=16, edges=11, isolated=9, largest=5, physical=4),
                "homo_sapiens": dict(seeds=385, mapped=382, edges=2210, isolated=70, largest=295, physical=1152)},
        "subset": dict(mean=3.5, mean_tol=0.15, below_pct=97.7),
        "intact": {"deinococcus_radiodurans": dict(edges=0, phys_pct=0.0), "homo_sapiens": dict(edges=522, phys_pct=40.5)},
        "louvain": {"homo_sapiens": dict(ari=0.74, modularity=0.54)},
        "partition_modularity": {"homo_sapiens": 0.54},
        "alignment": {1.0: (9, 0.10), 0.8: (8, 0.63), 0.5: (7, 0.88)},
        "lp": {"homo_sapiens": dict(netmf=0.935, aa=0.949)},
        "examples": {("katA", "CAT"): ("7", "ortholog_rbh"), ("sodA", "SOD2"): ("16", "ortholog_rbh"),
                     ("DR_0846", "PRDX5"): ("1", "homolog"), ("sodA", "SOD3"): ("19", "candidate_analogue"),
                     ("sodC", "SOD2"): ("76", "candidate_analogue"), ("katA", "PRDX2"): ("96", "candidate_analogue"),
                     ("DR_1538", "GPX8"): ("11", "shared_function_annotation"),
                     ("DR_0846", "GPX2"): ("30", "distant_homolog")},
    },
    "dna_repair": {
        "net": {"deinococcus_radiodurans": dict(seeds=56, mapped=52, edges=119, isolated=11, largest=41, physical=14),
                "homo_sapiens": dict(seeds=432, mapped=426, edges=5285, isolated=23, largest=400, physical=2393)},
        "subset": dict(mean=77, mean_tol=1.5, below_pct=95.3),
        "intact": {"deinococcus_radiodurans": dict(edges=0, phys_pct=0.0), "homo_sapiens": dict(edges=744, phys_pct=28.0)},
        "louvain": {"deinococcus_radiodurans": dict(ari=0.76, modularity=0.39), "homo_sapiens": dict(ari=0.64)},
        "partition_modularity": {"deinococcus_radiodurans": 0.39},
        "alignment": {1.0: (114, 0.08), 0.8: (87, 0.02), 0.5: (89, 0.01)},
        "lp": {"homo_sapiens": dict(netmf=0.879, aa=0.952), "deinococcus_radiodurans": dict(netmf=0.728, aa=0.931)},
        "examples": {("nth", "NTHL1"): ("2", "ortholog_rbh"), ("DR_0354", "APEX1"): ("4", "ortholog_rbh"),
                     ("mutM", "NTHL1"): ("26", "candidate_analogue"), ("nth", "NEIL3"): ("31", "candidate_analogue"),
                     ("mutM", "OGG1"): ("47", "candidate_analogue"), ("DR_2584", "MPG"): ("62", "candidate_analogue"),
                     ("recO", "RAD52"): ("176", "similar_pathway_role"), ("ligA", "LIG1"): ("514", "distant_homolog")},
    },
}


def main(argv):
    for process, cfg in PROCESSES.items():
        claims = CLAIMS[process]
        version_file = raw_files(f"data/manifests/{cfg['phase1']}.json", "version")[0]
        record(process, "STRING version", "12.5", json.loads(version_file.read_text(encoding="utf-8"))[0]["string_version"])
        nets = {sp: build_species(cfg["phase1"], sp) for sp in SPECIES}
        network_checks(process, nets, claims["net"])
        random_subset_check(process, nets, claims["subset"])
        intact_check(process, nets, cfg["phase2"], claims["intact"])
        louvain_check(process, nets, claims["louvain"])
        pipeline_partition_modularity(process, nets, cfg["phase2"], claims["partition_modularity"])
        alignment_check(process, nets, cfg["phase2"], claims["alignment"])
        link_prediction_check(process, nets, claims["lp"])
        ranking_check(process, cfg["phase2"], claims["examples"])
    width = max(len(r["claim"]) for r in RESULTS)
    for r in RESULTS:
        print(f"{r['status']:8} {r['process']:10} {r['claim']:<{width}}  manuscript={r['manuscript']}  "
              f"recomputed={r['recomputed']}" + (f"  [{r['note']}]" if r["note"] else ""))
    counts = Counter(r["status"] for r in RESULTS)
    print(dict(counts))
    if "--json" in argv:
        Path(argv[argv.index("--json") + 1]).write_text(json.dumps(RESULTS, indent=2, default=str), encoding="utf-8")
    return 1 if counts.get("DIFFERS") else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
