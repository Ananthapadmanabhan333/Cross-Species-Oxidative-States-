"""Validation of seeds, identifier mappings and edges.

Each function returns ``(retained, excluded)``. ``excluded`` holds the original record plus an
``exclusion_reason``; a record is assigned the first reason it fails, in the order documented in
each function, so counts per reason add up to the number of excluded records.
"""

from __future__ import annotations

import re

import pandas as pd

UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
STRING_ID_RE = re.compile(r"^(\d+)\.(\S+)$")


def is_uniprot_accession(value: object) -> bool:
    return isinstance(value, str) and bool(UNIPROT_ACCESSION_RE.match(value))


def string_id_taxon(value: object) -> int | None:
    """Taxon prefix of a STRING protein id ('9606.ENSP00000270142' -> 9606), or None if malformed."""
    if not isinstance(value, str):
        return None
    match = STRING_ID_RE.match(value)
    return int(match.group(1)) if match else None


def _split(frame: pd.DataFrame, reasons: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    excluded = frame[reasons.notna()].copy()
    excluded["exclusion_reason"] = reasons[reasons.notna()]
    return frame[reasons.isna()].copy(), excluded


def _assign(reasons: pd.Series, mask: pd.Series, reason: str) -> None:
    reasons[mask & reasons.isna()] = reason


# ---------------------------------------------------------------------------------------- seeds
def validate_seeds(seeds: pd.DataFrame, allowed_organism_ids: set[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Checks, in order: missing_accession, invalid_uniprot_accession, unexpected_organism,
    duplicate_accession (first occurrence kept)."""
    reasons = pd.Series(pd.NA, index=seeds.index, dtype="object")
    accession = seeds["uniprot_accession"].fillna("")
    _assign(reasons, accession.eq(""), "missing_accession")
    _assign(reasons, ~accession.map(is_uniprot_accession), "invalid_uniprot_accession")
    organism = pd.to_numeric(seeds["uniprot_organism_id"], errors="coerce")
    _assign(reasons, ~organism.isin(list(allowed_organism_ids)), "unexpected_organism")
    _assign(reasons, accession.duplicated(keep="first"), "duplicate_accession")
    return _split(seeds, reasons)


# -------------------------------------------------------------------------------------- mapping
def validate_mapping(seeds: pd.DataFrame, mapping: pd.DataFrame,
                     expected_string_taxon: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join seeds to their STRING mapping.

    Checks, in order: unmapped_in_string, invalid_string_id, string_taxon_mismatch,
    collapsed_to_same_string_id (several accessions map to one STRING protein; the first is kept
    and the others are listed in ``uniprot_accessions`` of the retained row).

    Adds ``xref_agrees``: whether UniProt's own STRING cross-reference contains the STRING id the
    STRING API returned (NA when UniProt has no cross-reference). Disagreement is reported, not
    used to exclude.
    """
    first_hit = mapping.drop_duplicates("uniprot_accession", keep="first")
    merged = seeds.merge(first_hit, on="uniprot_accession", how="left")

    reasons = pd.Series(pd.NA, index=merged.index, dtype="object")
    _assign(reasons, merged["string_id"].isna(), "unmapped_in_string")
    prefix = merged["string_id"].map(string_id_taxon)
    _assign(reasons, prefix.isna(), "invalid_string_id")
    reported_taxon = pd.to_numeric(merged["string_taxon_id"], errors="coerce")
    _assign(reasons, (prefix != expected_string_taxon) | (reported_taxon != expected_string_taxon),
            "string_taxon_mismatch")

    candidate = reasons.isna()
    groups = (merged[candidate].groupby("string_id")["uniprot_accession"]
              .agg(lambda s: ";".join(s)).rename("uniprot_accessions"))
    dup = candidate & merged["string_id"].duplicated(keep="first")
    _assign(reasons, dup, "collapsed_to_same_string_id")

    def agrees(row: pd.Series) -> object:
        xref = row.get("uniprot_xref_string")
        if not isinstance(xref, str) or not xref or pd.isna(row["string_id"]):
            return pd.NA
        return row["string_id"] in xref.split(";")

    merged["xref_agrees"] = merged.apply(agrees, axis=1).astype("boolean")
    retained, excluded = _split(merged, reasons)
    retained = retained.merge(groups, left_on="string_id", right_index=True, how="left")
    return retained, excluded


# ---------------------------------------------------------------------------------------- edges
def validate_edges(edges: pd.DataFrame, node_ids: set[str], expected_string_taxon: int,
                   score_threshold: float, allow_new_nodes: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Checks, in order:
    missing_identifier, invalid_string_id, species_mismatch, self_loop, missing_score,
    score_out_of_range, endpoint_not_in_node_set (only when ``allow_new_nodes`` is False),
    duplicate_edge (same unordered pair; the highest combined score is kept), below_threshold.
    """
    edges = edges.copy()
    reasons = pd.Series(pd.NA, index=edges.index, dtype="object")
    a, b = edges["string_id_a"], edges["string_id_b"]

    _assign(reasons, a.isna() | b.isna() | a.eq("") | b.eq(""), "missing_identifier")
    ta, tb = a.map(string_id_taxon), b.map(string_id_taxon)
    _assign(reasons, ta.isna() | tb.isna(), "invalid_string_id")
    reported = pd.to_numeric(edges["string_taxon_id"], errors="coerce")
    _assign(reasons, (ta != expected_string_taxon) | (tb != expected_string_taxon)
            | (reported.notna() & (reported != expected_string_taxon)), "species_mismatch")
    _assign(reasons, a.eq(b), "self_loop")
    score = edges["combined_score"]
    _assign(reasons, score.isna(), "missing_score")
    _assign(reasons, (score < 0) | (score > 1), "score_out_of_range")
    if not allow_new_nodes:
        _assign(reasons, ~(a.isin(node_ids) & b.isin(node_ids)), "endpoint_not_in_node_set")

    pair = [tuple(sorted(p)) if isinstance(p[0], str) and isinstance(p[1], str) else (None, None)
            for p in zip(a, b)]
    edges["node_u"] = [p[0] for p in pair]
    edges["node_v"] = [p[1] for p in pair]
    still_ok = reasons.isna()
    ranked = edges[still_ok].sort_values(["combined_score", "source_record"], ascending=[False, True])
    dup_index = ranked[ranked.duplicated(["node_u", "node_v"], keep="first")].index
    reasons.loc[dup_index] = "duplicate_edge"
    _assign(reasons, score < score_threshold, "below_threshold")
    return _split(edges, reasons)


def reason_counts(excluded: pd.DataFrame) -> dict[str, int]:
    if excluded.empty:
        return {}
    return {str(k): int(v) for k, v in excluded["exclusion_reason"].value_counts().sort_index().items()}
