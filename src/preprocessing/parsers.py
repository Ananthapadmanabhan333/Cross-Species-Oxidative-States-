"""Convert raw UniProt TSV and STRING JSON payloads into pandas DataFrames."""

from __future__ import annotations

import io
from typing import Any

import pandas as pd

from ..data_collection.uniprot_client import SEED_FIELDS

# Column names used throughout the pipeline for the UniProt fields requested in SEED_FIELDS.
SEED_COLUMNS = {
    "accession": "uniprot_accession",
    "id": "uniprot_entry_name",
    "reviewed": "uniprot_reviewed",
    "gene_primary": "gene_primary",
    "gene_names": "gene_names",
    "protein_name": "protein_name",
    "organism_id": "uniprot_organism_id",
    "organism_name": "uniprot_organism_name",
    "xref_string": "uniprot_xref_string",
    "go_id": "go_ids",
    "annotation_score": "uniprot_annotation_score",
}

# STRING channel scores, as documented for the `network` endpoint.
STRING_CHANNELS = {
    "nscore": "neighborhood",
    "fscore": "fusion",
    "pscore": "cooccurrence",
    "ascore": "coexpression",
    "escore": "experimental",
    "dscore": "database",
    "tscore": "textmining",
}
NETWORK_COLUMNS = ["stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "ncbiTaxonId",
                   "score", *STRING_CHANNELS]


def parse_uniprot_seeds(text: str) -> pd.DataFrame:
    """Parse a UniProt TSV whose columns were requested in SEED_FIELDS order."""
    expected = [SEED_COLUMNS[f] for f in SEED_FIELDS]
    if not text.strip():
        return pd.DataFrame(columns=expected)
    frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    if frame.shape[1] != len(SEED_FIELDS):
        raise ValueError(f"UniProt TSV has {frame.shape[1]} columns, expected {len(SEED_FIELDS)} "
                         f"({', '.join(SEED_FIELDS)})")
    frame.columns = expected
    frame = frame.apply(lambda col: col.str.strip())
    # "243230.DR_1998;" -> "243230.DR_1998" (several ids are ';'-separated)
    frame["uniprot_xref_string"] = frame["uniprot_xref_string"].str.rstrip(";")
    frame["uniprot_reviewed"] = frame["uniprot_reviewed"].str.lower().eq("reviewed")
    return frame


def parse_string_mapping(rows: list[dict[str, Any]]) -> pd.DataFrame:
    columns = ["queryItem", "queryIndex", "stringId", "ncbiTaxonId", "taxonName", "preferredName", "annotation"]
    frame = pd.DataFrame(rows)
    for col in columns:
        if col not in frame.columns:
            frame[col] = pd.NA
    return frame[columns].rename(columns={
        "queryItem": "uniprot_accession", "queryIndex": "query_index", "stringId": "string_id",
        "ncbiTaxonId": "string_taxon_id", "taxonName": "string_taxon_name",
        "preferredName": "string_preferred_name", "annotation": "string_annotation",
    })


def parse_string_network(rows: list[dict[str, Any]], network_type: str) -> pd.DataFrame:
    """One row per STRING record, unchanged except for renaming; missing fields become NA."""
    frame = pd.DataFrame(rows)
    for col in NETWORK_COLUMNS:
        if col not in frame.columns:
            frame[col] = pd.NA
    frame = frame[NETWORK_COLUMNS].copy()
    for col in ["score", *STRING_CHANNELS]:
        # float even when every value is an integer (STRING writes 0 rather than 0.0)
        frame[col] = pd.to_numeric(frame[col], errors="coerce").astype(float)
    frame = frame.rename(columns={
        "stringId_A": "string_id_a", "stringId_B": "string_id_b",
        "preferredName_A": "preferred_name_a", "preferredName_B": "preferred_name_b",
        "ncbiTaxonId": "string_taxon_id", "score": "combined_score",
        **{k: f"{v}_score" for k, v in STRING_CHANNELS.items()},
    })
    frame.insert(0, "source_record", range(len(frame)))
    frame["network_type"] = network_type
    return frame
