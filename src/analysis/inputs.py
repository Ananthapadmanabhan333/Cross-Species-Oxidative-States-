"""Load and verify Phase 1 outputs, and rebuild species graphs at any threshold offline.

Before anything is analysed, every Phase 1 output and raw API response listed in the run manifest
is re-checksummed. A mismatch stops the analysis (the input would no longer be traceable).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd

from ..data_collection.http_cache import sha256_file
from ..network_construction.graph_builder import annotate_edges, build_graph
from ..quality_control.validators import validate_edges


class InputIntegrityError(RuntimeError):
    """A Phase 1 file is missing or no longer matches its recorded checksum."""


@dataclass
class SpeciesInput:
    key: str
    spec: dict[str, Any]
    nodes: pd.DataFrame                 # Phase 1 nodes.tsv
    functional_records: pd.DataFrame    # every STRING functional record retrieved (>= 0.4)
    physical_records: pd.DataFrame | None


def verify_phase1_run(manifest_path: Path, root: Path) -> dict[str, Any]:
    if not manifest_path.is_file():
        raise InputIntegrityError(f"Phase 1 manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    problems = []
    for item in manifest["outputs"]:
        path = root / item["path"]
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            problems.append(item["path"])
    for item in manifest["raw_inputs"]:
        # Phase 1 manifests written on Windows store native separators; normalise for portability.
        path = Path(str(item["data_file"]).replace("\\", "/"))
        path = path if path.is_absolute() else root / path
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            problems.append(str(item["data_file"]))
    if problems:
        raise InputIntegrityError(f"{len(problems)} Phase 1 files fail verification, e.g. {problems[:3]}")
    return {"run_id": manifest["run_id"], "string_version": manifest["string_version"],
            "outputs_verified": len(manifest["outputs"]), "raw_inputs_verified": len(manifest["raw_inputs"]),
            "manifest": str(manifest_path)}


def load_species(processed_run_dir: Path, key: str, spec: dict[str, Any]) -> SpeciesInput:
    base = processed_run_dir / key
    nodes = pd.read_csv(base / "nodes.tsv", sep="\t", dtype={"uniprot_organism_id": str})
    functional = pd.read_csv(base / "string_functional_records_parsed.tsv", sep="\t")
    physical_path = base / "string_physical_records_parsed.tsv"
    physical = pd.read_csv(physical_path, sep="\t") if physical_path.is_file() else None
    return SpeciesInput(key, spec, nodes, functional, physical)


def graph_at_threshold(species: SpeciesInput, threshold: float, string_version: str) -> nx.Graph:
    """Re-apply Phase 1 QC and graph construction to the stored records at `threshold`."""
    node_ids = set(species.nodes["string_id"])
    taxon = species.spec["string_taxon_id"]
    functional, _ = validate_edges(species.functional_records, node_ids, taxon, threshold)
    physical = None
    if species.physical_records is not None:
        physical, _ = validate_edges(species.physical_records, node_ids, taxon, threshold)
    edges, _ = annotate_edges(functional, physical, 0.0, f"STRING v{string_version} functional association")
    nodes = species.nodes.drop(columns=[c for c in ("has_retained_edge",) if c in species.nodes.columns])
    return build_graph(nodes, edges, keep_isolated=True,
                       graph_attrs={"species": species.spec["label"], "threshold": threshold,
                                    "string_taxon_id": taxon})
