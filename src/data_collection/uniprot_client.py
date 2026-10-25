"""Seed-protein retrieval from the UniProtKB REST API (https://www.uniprot.org/help/api).

Oxidative-stress defence proteins are defined operationally: UniProtKB entries of the target
organism annotated with any of the configured GO terms, exactly as resolved by UniProt's ``go:``
query field. The verbatim query string is stored with every raw response.
"""

from __future__ import annotations

from typing import Any

from .http_cache import CachedHttpClient, RawResponse

SEED_FIELDS = (
    "accession", "id", "reviewed", "gene_primary", "gene_names", "protein_name",
    "organism_id", "organism_name", "xref_string", "go_id", "annotation_score",
)


def build_seed_query(go_terms: list[str], taxon_query: str, reviewed_only: bool) -> str:
    go_clause = " OR ".join(f"(go:{term.split(':', 1)[1]})" for term in go_terms)
    query = f"({go_clause}) AND ({taxon_query})"
    if reviewed_only:
        query += " AND (reviewed:true)"
    return query


class UniProtClient:
    def __init__(self, http: CachedHttpClient, base_url: str = "https://rest.uniprot.org/uniprotkb") -> None:
        self.http = http
        self.base_url = base_url.rstrip("/")

    def fetch_seed_proteins(self, go_terms: list[str], taxon_query: str, reviewed_only: bool,
                            species_key: str) -> RawResponse:
        query = build_seed_query(go_terms, taxon_query, reviewed_only)
        params: dict[str, Any] = {"query": query, "format": "tsv", "fields": ",".join(SEED_FIELDS)}
        return self.http.fetch(
            source="uniprot", endpoint=f"seeds_{species_key}", url=f"{self.base_url}/stream",
            params=params, method="GET", ext="tsv",
            description=f"UniProtKB entries annotated with {', '.join(go_terms)} for {taxon_query}",
        )

    def fetch_entry_count(self, taxon_query: str, reviewed_only: bool, species_key: str) -> tuple[int, RawResponse]:
        """Number of UniProtKB entries for the organism (context for annotation coverage)."""
        query = f"({taxon_query})" + (" AND (reviewed:true)" if reviewed_only else "")
        params = {"query": query, "format": "tsv", "fields": "accession", "size": 1}
        raw = self.http.fetch(
            source="uniprot", endpoint=f"entry_count_{species_key}", url=f"{self.base_url}/search",
            params=params, method="GET", ext="tsv",
            description=f"UniProtKB entry count for {query} (read from the X-Total-Results header)",
        )
        total = raw.headers.get("x-total-results")
        if total is None or not str(total).isdigit():
            raise ValueError(f"UniProt response for {query!r} has no usable X-Total-Results header")
        return int(total), raw
