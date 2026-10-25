"""Per-protein evidence profiles and pairwise cross-species similarity components.

Every component returns a value in [0, 1] or ``None`` when the evidence needed to judge it is
absent for either protein. ``None`` means "not assessable", never "dissimilar".
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..data_collection.evidence_sources import split_list

FULL_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


@dataclass
class ProteinEvidence:
    string_id: str
    species_key: str
    accession: str | None
    name: str | None
    go_terms: set[str] = field(default_factory=set)
    pathways: set[str] = field(default_factory=set)
    ec: set[str] = field(default_factory=set)
    panther: set[str] = field(default_factory=set)
    pfam: set[str] = field(default_factory=set)
    interpro: set[str] = field(default_factory=set)
    literature_level: float | None = None      # 1 experimental, 0.5 curated/phylo/computational, 0.25 IEA only
    evidence_codes: set[str] = field(default_factory=set)
    experimental_references: set[str] = field(default_factory=set)
    perturbation: bool = False
    negated_go: set[str] = field(default_factory=set)   # "NOT" annotations: contradictory evidence


def kegg_map(term: str) -> str | None:
    """'dra00380' / 'hsa00380' -> '00380' (species-independent KEGG map number)."""
    match = re.match(r"^[a-z]{2,4}(\d{5})$", str(term))
    return match.group(1) if match else None


def build_profiles(nodes: pd.DataFrame, annotations: pd.DataFrame, uniprot: pd.DataFrame,
                   quickgo: pd.DataFrame, species_key: str, go_categories: list[str],
                   pathway_category: str, excluded_maps: set[str]) -> dict[str, ProteinEvidence]:
    up = uniprot.set_index("uniprot_accession") if not uniprot.empty else pd.DataFrame()
    profiles: dict[str, ProteinEvidence] = {}
    ann_by_node = {sid: grp for sid, grp in annotations.groupby("string_id")}
    for row in nodes.itertuples(index=False):
        sid = row.string_id
        accessions = split_list(getattr(row, "uniprot_accessions", None)) or split_list(row.uniprot_accession)
        prof = ProteinEvidence(sid, species_key, accessions[0] if accessions else None, row.string_preferred_name)
        grp = ann_by_node.get(sid)
        if grp is not None:
            prof.go_terms = set(grp.loc[grp["category"].isin(go_categories), "term"])
            maps = {kegg_map(t) for t in grp.loc[grp["category"] == pathway_category, "term"]}
            prof.pathways = {m for m in maps if m and m not in excluded_maps}
            prof.interpro = set(grp.loc[grp["category"] == "InterPro", "term"])
        for acc in accessions:
            if acc in up.index:
                rec = up.loc[acc]
                prof.ec |= set(split_list(rec["ec_numbers"]))
                prof.panther |= {p.split(":")[0] for p in split_list(rec["panther_families"])}
                prof.pfam |= set(split_list(rec["pfam"]))
        qg = quickgo[quickgo["uniprot_accession"].isin(accessions)] if len(quickgo) else quickgo
        if len(qg):
            negated = qg["qualifier"].fillna("").str.startswith("NOT")
            prof.negated_go = set(qg.loc[negated, "go_id"])
            positive = qg[~negated]
            prof.evidence_codes = set(positive["go_evidence"])
            classes = set(positive["evidence_class"])
            if "experimental" in classes:
                prof.literature_level = 1.0
            elif classes & {"phylogenetic", "computational", "author_curator"}:
                prof.literature_level = 0.5
            elif "electronic" in classes:
                prof.literature_level = 0.25
            exp = positive[positive["evidence_class"] == "experimental"]
            prof.experimental_references = {r for r in exp["reference"].dropna() if str(r).startswith("PMID:")}
            prof.perturbation = bool(positive["is_perturbation"].any())
        profiles[sid] = prof
    return profiles


def information_content(term_sets: list[set[str]]) -> dict[str, float]:
    """IC(t) = -ln(fraction of annotated proteins carrying t)."""
    annotated = [s for s in term_sets if s]
    counts = Counter(t for s in annotated for t in s)
    total = len(annotated)
    return {t: -math.log(c / total) for t, c in counts.items()} if total else {}


SIMILARITY_MEASURES = ("ochiai", "simgic", "min")


def simgic(a: set[str] | None, b: set[str] | None, ic: dict[str, float], measure: str = "simgic") -> float | None:
    """IC-weighted set similarity. Shared information S = sum IC(a & b), normalised by
    * simgic: sum IC(a | b)            (Pesquita et al. 2008; penalises richly annotated proteins)
    * ochiai: sqrt(sum IC(a) * sum IC(b))  (IC-weighted cosine; default, see docs/PHASE2_METHODS.md)
    * min:    min(sum IC(a), sum IC(b))    (overlap coefficient; permissive for sparse annotation)
    """
    if not a or not b:
        return None
    shared = sum(ic.get(t, 0.0) for t in a & b)
    if measure == "simgic":
        denom = sum(ic.get(t, 0.0) for t in a | b)
    else:
        ia, ib = sum(ic.get(t, 0.0) for t in a), sum(ic.get(t, 0.0) for t in b)
        if measure == "ochiai":
            denom = math.sqrt(ia * ib)
        elif measure == "min":
            denom = min(ia, ib)
        else:
            raise ValueError(f"unknown similarity measure {measure!r}; use one of {SIMILARITY_MEASURES}")
    return shared / denom if denom > 0 else None


def ec_similarity(a: set[str], b: set[str]) -> tuple[float | None, list[str]]:
    """1.0 shared complete EC; 0.5 shared sub-subclass (first three levels); 0.0 both annotated
    but no shared sub-subclass; None if either has no EC number."""
    if not a or not b:
        return None, []
    full = sorted({e for e in a & b if FULL_EC.match(e)})
    if full:
        return 1.0, full

    def level3(ecs: set[str]) -> set[str]:
        return {".".join(e.split(".")[:3]) for e in ecs if len(e.split(".")) >= 3 and e.split(".")[2] != "-"}
    shared3 = sorted(level3(a) & level3(b))
    if shared3:
        return 0.5, [f"{s}.-" for s in shared3]
    return 0.0, []


@dataclass
class HomologyTables:
    dr_best: dict[str, tuple[str, float]]   # D. radiodurans protein -> (best human hit, bitscore)
    hs_best: dict[str, tuple[str, float]]   # human protein -> (best D. radiodurans hit, bitscore)

    @classmethod
    def from_frames(cls, dr_to_hs: pd.DataFrame, hs_to_dr: pd.DataFrame) -> "HomologyTables":
        def best(frame: pd.DataFrame) -> dict[str, tuple[str, float]]:
            if frame.empty:
                return {}
            top = frame.sort_values("bitscore", ascending=False).drop_duplicates("query_id")
            return {r.query_id: (r.hit_id, float(r.bitscore)) for r in top.itertuples()}
        return cls(best(dr_to_hs), best(hs_to_dr))


def homology(a: ProteinEvidence, b: ProteinEvidence, tables: HomologyTables,
             min_bitscore: float) -> tuple[float, str, dict]:
    """Tiered homology evidence. Every protein was searched (STRING homology_best), so absence of a
    hit is an observation and scores 0 - but 0 here only means 'no detectable homology', which is
    precisely what an analogue would show."""
    hit_ab = tables.dr_best.get(a.string_id)
    hit_ba = tables.hs_best.get(b.string_id)
    ab = hit_ab is not None and hit_ab[0] == b.string_id and hit_ab[1] >= min_bitscore
    ba = hit_ba is not None and hit_ba[0] == a.string_id and hit_ba[1] >= min_bitscore
    shared_panther = sorted(a.panther & b.panther)
    shared_pfam = sorted(a.pfam & b.pfam)
    shared_interpro = sorted(a.interpro & b.interpro)
    detail = {"best_hit_dr_to_hs": hit_ab[0] if hit_ab else None,
              "bitscore_dr_to_hs": hit_ab[1] if hit_ab else None,
              "best_hit_hs_to_dr": hit_ba[0] if hit_ba else None,
              "shared_panther": shared_panther, "shared_pfam": shared_pfam, "shared_interpro": shared_interpro}
    if ab and ba:
        return 1.0, "reciprocal_best_hit", detail
    if shared_panther:
        return 0.8, "shared_panther_family", detail
    if ab or ba:
        return 0.6, "one_way_best_hit", detail
    if shared_pfam:
        return 0.5, "shared_pfam_domain", detail
    if shared_interpro:
        # InterPro entries (domains, homologous superfamilies) integrate several signature databases
        # and can reveal remote homology that Pfam best-hit models miss (e.g. thioredoxin-fold
        # peroxiredoxin subfamilies). Treated as weak homology so that "analogue" claims require
        # the absence of any shared signature.
        return 0.4, "shared_interpro_entry", detail
    return 0.0, "none_detected", detail


def role_similarity(a: pd.Series | None, b: pd.Series | None, features: list[str]) -> float | None:
    if a is None or b is None or bool(a.get("is_isolated")) or bool(b.get("is_isolated")):
        return None
    diffs = [abs(a[f] - b[f]) for f in features if pd.notna(a.get(f)) and pd.notna(b.get(f))]
    return float(1.0 - np.mean(diffs)) if diffs else None


def module_core_terms(membership: dict[str, int], profiles: dict[str, ProteinEvidence],
                      min_share: float = 0.5) -> dict[int, set[str]]:
    """Terms (GO + KEGG map) carried by >= min_share of a module's members; singletons excluded."""
    modules: dict[int, list[str]] = {}
    for n, m in membership.items():
        modules.setdefault(m, []).append(n)
    cores = {}
    for m, members in modules.items():
        if len(members) < 2:
            continue
        counts = Counter(t for n in members for t in (profiles[n].go_terms | {f"KEGG:{p}" for p in profiles[n].pathways}))
        cores[m] = {t for t, c in counts.items() if c / len(members) >= min_share}
    return cores
