"""Build the per-species protein (node) table and measure how well it is annotated.

Annotations come only from the retrieved records (UniProt seed rows and STRING mapping rows);
nothing is inferred. Missing values stay empty and are counted in the coverage metrics.
"""

from __future__ import annotations

import pandas as pd

STRING_NO_ANNOTATION = "annotation not available"

NODE_COLUMNS = [
    "string_id", "species_key", "species_label", "ncbi_taxon_id", "string_taxon_id",
    "string_preferred_name", "string_annotation", "uniprot_accession", "uniprot_accessions",
    "uniprot_entry_name", "uniprot_reviewed", "gene_primary", "gene_names", "protein_name",
    "uniprot_organism_id", "go_ids", "matched_seed_go_terms", "xref_agrees", "is_seed",
]


def matched_terms(go_ids: object, seed_terms: list[str]) -> str:
    """Seed GO terms listed verbatim on the entry. Empty means UniProt matched the entry through a
    descendant term of a seed term."""
    if not isinstance(go_ids, str) or not go_ids:
        return ""
    present = {t.strip() for t in go_ids.split(";")}
    return ";".join(t for t in seed_terms if t in present)


def build_node_table(mapped_seeds: pd.DataFrame, species_key: str, species: dict,
                     seed_terms: list[str]) -> pd.DataFrame:
    nodes = mapped_seeds.copy()
    nodes["species_key"] = species_key
    nodes["species_label"] = species["label"]
    nodes["ncbi_taxon_id"] = species["ncbi_taxon_id"]
    nodes["string_taxon_id"] = species["string_taxon_id"]
    nodes["matched_seed_go_terms"] = nodes["go_ids"].map(lambda g: matched_terms(g, seed_terms))
    nodes["is_seed"] = True
    for col in NODE_COLUMNS:
        if col not in nodes.columns:
            nodes[col] = pd.NA
    return nodes[NODE_COLUMNS].sort_values("string_id").reset_index(drop=True)


def add_network_only_nodes(nodes: pd.DataFrame, edges: pd.DataFrame, species_key: str,
                           species: dict) -> pd.DataFrame:
    """Nodes STRING added around the seeds (only when string.add_nodes > 0). They carry the
    STRING preferred name from the edge records and no UniProt annotation."""
    known = set(nodes["string_id"])
    extra: dict[str, str] = {}
    for col, name_col in (("string_id_a", "preferred_name_a"), ("string_id_b", "preferred_name_b")):
        for sid, name in zip(edges[col], edges[name_col]):
            if sid not in known:
                extra.setdefault(sid, name)
    if not extra:
        return nodes
    added = pd.DataFrame({
        "string_id": list(extra), "string_preferred_name": list(extra.values()),
        "species_key": species_key, "species_label": species["label"],
        "ncbi_taxon_id": species["ncbi_taxon_id"], "string_taxon_id": species["string_taxon_id"],
        "is_seed": False,
    })
    return (pd.concat([nodes, added.reindex(columns=NODE_COLUMNS)], ignore_index=True)
            .sort_values("string_id").reset_index(drop=True))


def _has_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def annotation_coverage(nodes: pd.DataFrame) -> dict[str, object]:
    total = len(nodes)

    def share(mask: pd.Series) -> dict[str, object]:
        count = int(mask.sum())
        return {"count": count, "fraction": round(count / total, 4) if total else None}

    string_annotation = nodes["string_annotation"].fillna("").astype(str).str.strip()
    return {
        "nodes": total,
        "with_uniprot_accession": share(_has_text(nodes["uniprot_accession"])),
        "uniprot_reviewed": share(nodes["uniprot_reviewed"].fillna(False).astype(bool)),
        "with_gene_name": share(_has_text(nodes["gene_primary"])),
        "with_protein_name": share(_has_text(nodes["protein_name"])),
        "with_string_annotation": share(string_annotation.ne("")
                                        & string_annotation.str.lower().ne(STRING_NO_ANNOTATION)),
        "with_seed_go_term_listed_verbatim": share(_has_text(nodes["matched_seed_go_terms"])),
        "uniprot_xref_disagrees_with_string_mapping": share(nodes["xref_agrees"].eq(False).fillna(False)),
    }
