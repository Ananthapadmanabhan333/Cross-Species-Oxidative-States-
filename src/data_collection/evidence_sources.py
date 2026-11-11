"""Phase 2 evidence retrieval: annotations, homology and GO evidence codes.

Sources (all documented public APIs, all responses cached with provenance by CachedHttpClient):
  * STRING ``functional_annotation`` - per-protein GO, KEGG, Reactome, Pfam, InterPro ... terms
  * STRING ``homology_best``         - best Smith-Waterman hit of each protein in the other species
  * STRING ``ppi_enrichment``        - test of whether a protein set has more associations than
                                       expected for a random set of the same size from the genome
  * UniProtKB REST                   - EC numbers, PANTHER families, Pfam, KEGG genes, Reactome,
                                       PubMed references cited by the entry
  * QuickGO (EBI)                    - GO annotations with evidence codes and references

  * KEGG REST ``link/ko/<org>``      - KEGG Orthology (KO) assignment of every gene of an organism
  * IntAct PSICQUIC (EBI)            - experimentally determined molecular interactions (MITAB)

KEGG REST is used once per organism (cached) for KO assignments. KEGG's terms allow free academic
use; commercial use requires a licence (https://www.kegg.jp/kegg/legal.html). Pathway maps are
taken from STRING's redistribution of KEGG.
"""

from __future__ import annotations

import io
import json
import re
from typing import Any
from urllib.parse import quote

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


# ------------------------------------------------------------------------------------- KEGG
KEGG_REST = "https://rest.kegg.jp"


def kegg_ko_links(http: CachedHttpClient, org_code: str) -> dict[str, set[str]]:
    """KEGG gene id ('dra:DR_1998') -> set of KO ids ('K03781') for one organism."""
    raw = http.fetch("kegg", f"link_ko_{org_code}", f"{KEGG_REST}/link/ko/{org_code}", {}, method="GET",
                     ext="tsv", description=f"KEGG Orthology assignments of all {org_code} genes (academic use)")
    links: dict[str, set[str]] = {}
    for line in raw.text.splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[1].startswith("ko:"):
            links.setdefault(parts[0], set()).add(parts[1][3:])
    if not links:
        raise DataSourceError(f"KEGG link/ko/{org_code} returned no KO assignments")
    return links


# ------------------------------------------------------------------------------------ IntAct
INTACT_PSICQUIC = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query"
# PSI-MI interaction types accepted as evidence of physical binding. Proximity labelling (MI:2364)
# and colocalization (MI:0403) are recorded but do not count as confirmed physical interactions.
PHYSICAL_INTERACTION_TYPES = {"MI:0915": "physical association", "MI:0407": "direct interaction",
                              "MI:0914": "association", "MI:0195": "covalent binding",
                              "MI:0408": "disulfide bond"}
INTACT_PAGE = 2500


def _mi_id(field: str) -> str | None:
    match = re.search(r"MI:\d{4}", field or "")
    return match.group(0) if match else None


def intact_within_set(http: CachedHttpClient, accessions: list[str], species_key: str,
                      chunk: int = 60) -> pd.DataFrame:
    """IntAct interactions whose BOTH partners are in `accessions` (MITAB 2.5, all pages cached).

    Queries are split into chunk x chunk blocks so that URLs stay short; duplicate rows from
    overlapping blocks are removed by IntAct interaction id."""
    accs = sorted({a for a in accessions if a})
    blocks = [accs[i:i + chunk] for i in range(0, len(accs), chunk)]
    rows: list[list[str]] = []
    for i, a in enumerate(blocks):
        for b in blocks[i:]:
            query = f"idA:({' OR '.join(a)}) AND idB:({' OR '.join(b)})"
            for swap in ((query,) if a is b else (query, f"idA:({' OR '.join(b)}) AND idB:({' OR '.join(a)})")):
                first = 0
                while True:
                    url = f"{INTACT_PSICQUIC}/{quote(swap, safe='')}"
                    raw = http.fetch("intact", f"psicquic_{species_key}", url,
                                     {"format": "tab25", "firstResult": first, "maxResults": INTACT_PAGE},
                                     method="GET", ext="tsv",
                                     description=f"IntAct interactions within {species_key} seed proteins")
                    page = [line.split("\t") for line in raw.text.splitlines() if line.strip()]
                    rows.extend(r for r in page if len(r) >= 15)
                    if len(page) < INTACT_PAGE:
                        break
                    first += INTACT_PAGE
    records = []
    for r in rows:
        a = r[0].split(":", 1)[-1].split("-")[0]  # drop isoform suffixes (P04040-2 -> P04040)
        b = r[1].split(":", 1)[-1].split("-")[0]
        interaction_id = next((x.split(":", 1)[1] for x in r[13].split("|") if x.startswith("intact:")), r[13])
        miscore = next((x.split(":", 1)[1] for x in r[14].split("|") if x.startswith("intact-miscore:")), "")
        records.append({
            "uniprot_a": a, "uniprot_b": b, "interaction_id": interaction_id,
            "detection_method": r[6], "interaction_type": r[11], "interaction_type_mi": _mi_id(r[11]),
            "publications": ";".join(x for x in r[8].split("|") if x.startswith("pubmed:")),
            "intact_miscore": float(miscore) if miscore else None,
        })
    frame = pd.DataFrame(records, columns=["uniprot_a", "uniprot_b", "interaction_id", "detection_method",
                                           "interaction_type", "interaction_type_mi", "publications",
                                           "intact_miscore"])
    # A complex is written as several binary rows that share one interaction id, so a row is a
    # duplicate (from overlapping query blocks) only if its id AND both partners repeat.
    frame = frame.drop_duplicates(["interaction_id", "uniprot_a", "uniprot_b"]).reset_index(drop=True)
    frame["physical"] = frame["interaction_type_mi"].isin(PHYSICAL_INTERACTION_TYPES)
    # Co-complex association (MI:0914) places two proteins in one complex without showing that
    # they bind each other; the binary subset excludes it.
    frame["physical_binary"] = frame["physical"] & (frame["interaction_type_mi"] != "MI:0914")
    return frame
