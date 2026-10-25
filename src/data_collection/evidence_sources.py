"""Phase 2 evidence retrieval: annotations, homology and GO evidence codes.

Sources (all documented public APIs, all responses cached with provenance by CachedHttpClient):
  * STRING ``functional_annotation`` - per-protein GO, KEGG, Reactome, Pfam, InterPro ... terms
  * STRING ``homology_best``         - best Smith-Waterman hit of each protein in the other species
  * STRING ``ppi_enrichment``        - test of whether a protein set has more associations than
                                       expected for a random set of the same size from the genome
  * UniProtKB REST                   - EC numbers, PANTHER families, Pfam, KEGG genes, Reactome,
                                       PubMed references cited by the entry
  * QuickGO (EBI)                    - GO annotations with evidence codes and references

KEGG content is used only as redistributed inside STRING (pathway maps) and UniProt (gene
cross-references); the KEGG API itself is not called, in line with KEGG's access conditions.
"""

from __future__ import annotations

import io
import json
from typing import Any

import pandas as pd

from .http_cache import CachedHttpClient, DataSourceError
from .string_client import StringClient
from .uniprot_client import build_seed_query

UNIPROT_EVIDENCE_FIELDS = ("accession", "ec", "xref_panther", "xref_pfam", "xref_kegg", "xref_reactome",
                           "lit_pubmed_id")
UNIPROT_EVIDENCE_COLUMNS = ["uniprot_accession", "ec_numbers", "panther_families", "pfam", "kegg_genes",
                            "reactome", "uniprot_pubmed_ids"]

QUICKGO_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
QUICKGO_PAGE_SIZE = 200  # QuickGO maximum

# GO evidence codes grouped as in the GO Consortium guide (geneontology.org/docs/guide-go-evidence-codes)
EVIDENCE_CLASSES = {
    "experimental": {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"},
    "phylogenetic": {"IBA", "IBD", "IKR", "IRD"},
    "computational": {"ISS", "ISO", "ISA", "ISM", "IGC", "RCA"},
    "author_curator": {"TAS", "NAS", "IC", "ND"},
    "electronic": {"IEA"},
}
PERTURBATION_CODES = {"IMP", "IGI", "HMP", "HGI"}


def evidence_class(code: str) -> str:
    for name, codes in EVIDENCE_CLASSES.items():
        if code in codes:
            return name
    return "unknown"


# ----------------------------------------------------------------------------------- STRING
def string_functional_annotation(string: StringClient, string_ids: list[str], species: int,
                                 species_key: str) -> pd.DataFrame:
    """Long table: one row per (protein, term). All categories STRING returns are kept."""
    params = {"identifiers": "\r".join(sorted(set(string_ids))), "species": species, "allow_pubmed": 0,
              "caller_identity": string.caller_identity}
    raw = string.http.fetch("string", f"functional_annotation_{species_key}", string._url("functional_annotation"),
                            params, method="POST", ext="json",
                            description=f"STRING functional annotation of {len(set(string_ids))} proteins")
    rows = []
    for term in string.parse_json(raw):
        for sid in term.get("inputGenes", []):
            rows.append({"string_id": sid, "category": term["category"], "term": term["term"],
                         "description": term.get("description", "")})
    return pd.DataFrame(rows, columns=["string_id", "category", "term", "description"])


def string_homology_best(string: StringClient, string_ids: list[str], species: int, species_b: int,
                         species_key: str) -> pd.DataFrame:
    params = {"identifiers": "\r".join(sorted(set(string_ids))), "species": species, "species_b": species_b,
              "caller_identity": string.caller_identity}
    raw = string.http.fetch("string", f"homology_best_{species_key}_to_{species_b}", string._url("homology_best"),
                            params, method="POST", ext="json",
                            description=f"best homologs in taxon {species_b} for {len(set(string_ids))} proteins")
    frame = pd.DataFrame(string.parse_json(raw))
    cols = ["stringId_A", "stringId_B", "bitscore"]
    if frame.empty:
        return pd.DataFrame(columns=["query_id", "hit_id", "bitscore"])
    return frame[cols].rename(columns={"stringId_A": "query_id", "stringId_B": "hit_id"})


def string_ppi_enrichment(string: StringClient, string_ids: list[str], species: int, required_score: int,
                          species_key: str) -> dict[str, Any]:
    params = {"identifiers": "\r".join(sorted(set(string_ids))), "species": species,
              "required_score": required_score, "caller_identity": string.caller_identity}
    raw = string.http.fetch("string", f"ppi_enrichment_{species_key}", string._url("ppi_enrichment"),
                            params, method="POST", ext="json",
                            description=f"PPI enrichment test, required_score={required_score}")
    rows = string.parse_json(raw)
    if not rows:
        raise DataSourceError("STRING ppi_enrichment returned no rows")
    return {**rows[0], "required_score": required_score}


# ---------------------------------------------------------------------------------- UniProt
def uniprot_evidence(http: CachedHttpClient, base_url: str, go_terms: list[str], taxon_query: str,
                     reviewed_only: bool, species_key: str) -> pd.DataFrame:
    """Same entry set as the Phase 1 seed query, with orthology/function cross-reference fields."""
    query = build_seed_query(go_terms, taxon_query, reviewed_only)
    raw = http.fetch("uniprot", f"evidence_{species_key}", f"{base_url.rstrip('/')}/stream",
                     {"query": query, "format": "tsv", "fields": ",".join(UNIPROT_EVIDENCE_FIELDS)},
                     method="GET", ext="tsv", description="UniProt EC, PANTHER, Pfam, KEGG, Reactome, PubMed")
    if not raw.text.strip():
        return pd.DataFrame(columns=UNIPROT_EVIDENCE_COLUMNS)
    frame = pd.read_csv(io.StringIO(raw.text), sep="\t", dtype=str, keep_default_na=False)
    if frame.shape[1] != len(UNIPROT_EVIDENCE_FIELDS):
        raise ValueError(f"UniProt evidence TSV has {frame.shape[1]} columns, expected {len(UNIPROT_EVIDENCE_FIELDS)}")
    frame.columns = UNIPROT_EVIDENCE_COLUMNS
    return frame.apply(lambda col: col.str.strip())


def split_list(value: object, sep: str = ";") -> list[str]:
    if not isinstance(value, str) or not value.strip():
        return []
    return [v.strip() for v in value.split(sep) if v.strip()]


# ---------------------------------------------------------------------------------- QuickGO
def quickgo_annotations(http: CachedHttpClient, go_terms: list[str], ncbi_taxon: int,
                        species_key: str) -> pd.DataFrame:
    """All GO annotations to the seed terms or their descendants for the taxon (and its descendant
    taxa), with evidence codes and references. Pages are cached individually."""
    rows: list[dict[str, Any]] = []
    page, total_pages = 1, 1
    while page <= total_pages:
        params = {"taxonId": ncbi_taxon, "taxonUsage": "descendants", "goId": ",".join(go_terms),
                  "goUsage": "descendants", "goUsageRelationships": "is_a,part_of,occurs_in",
                  "limit": QUICKGO_PAGE_SIZE, "page": page}
        raw = http.fetch("quickgo", f"annotations_{species_key}", QUICKGO_URL, params, method="GET", ext="json",
                         headers={"Accept": "application/json"},
                         description=f"QuickGO annotations to {','.join(go_terms)} (+descendants), page {page}")
        payload = json.loads(raw.text)
        total_pages = int(payload.get("pageInfo", {}).get("total", 1) or 1)
        for item in payload.get("results", []):
            product = item.get("geneProductId", "")
            rows.append({
                "uniprot_accession": product.split(":", 1)[1] if product.startswith("UniProtKB:") else product,
                "go_id": item.get("goId"), "qualifier": item.get("qualifier"),
                "go_evidence": item.get("goEvidence"), "eco_code": item.get("evidenceCode"),
                "reference": item.get("reference"), "assigned_by": item.get("assignedBy"),
                "taxon_id": item.get("taxonId"), "date": item.get("date"),
            })
        page += 1
    frame = pd.DataFrame(rows, columns=["uniprot_accession", "go_id", "qualifier", "go_evidence", "eco_code",
                                        "reference", "assigned_by", "taxon_id", "date"])
    frame["evidence_class"] = frame["go_evidence"].map(lambda c: evidence_class(str(c)))
    frame["is_perturbation"] = frame["go_evidence"].isin(PERTURBATION_CODES)
    return frame
