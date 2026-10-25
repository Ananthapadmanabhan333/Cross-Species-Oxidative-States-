"""Coverage variant: a KO-projected D. radiodurans network.

The GO-defined D. radiodurans seed set is small (16 proteins) because the bacterium is sparsely
annotated, not necessarily because its defence network is small. As a *coverage sensitivity
analysis*, this module adds every D. radiodurans gene that KEGG assigns to a KEGG Orthology (KO)
group carried by a human seed protein, then retrieves and QCs the STRING network of the enlarged
set exactly as in Phase 1.

The variant is defined by orthology, so it is used only to test whether topology and alignment
conclusions depend on coverage. It is never used for candidate ranking, where it would be circular.
"""

from __future__ import annotations

from typing import Any

import networkx as nx
import pandas as pd

from ..data_collection.string_client import StringClient
from ..preprocessing.parsers import parse_string_mapping, parse_string_network
from ..quality_control.validators import reason_counts, string_id_taxon, validate_edges


def ko_projected_network(string: StringClient, dr_seed_ids: list[str], dr_taxon: int, human_kos: set[str],
                         dr_ko_links: dict[str, set[str]], retrieval_score: int, threshold: float,
                         key: str = "deinococcus_radiodurans_ko_projected") -> dict[str, Any]:
    projected_genes = sorted(g for g, kos in dr_ko_links.items() if kos & human_kos)
    loci = [g.split(":", 1)[1] for g in projected_genes]
    seed_loci = {sid.split(".", 1)[1] for sid in dr_seed_ids}
    query = sorted(set(loci) | seed_loci)

    rows, _ = string.map_identifiers(query, dr_taxon, key)
    mapping = parse_string_mapping(rows)
    mapping = mapping[mapping["string_id"].map(string_id_taxon) == dr_taxon]
    mapped_ids = sorted(set(mapping["string_id"]))
    unmapped = sorted(set(query) - set(mapping["uniprot_accession"]))  # queryItem column holds the locus

    raw, _ = string.network(mapped_ids, dr_taxon, retrieval_score, "functional", 0, key)
    edges = parse_string_network(raw, "functional")
    kept, dropped = validate_edges(edges, set(mapped_ids), dr_taxon, threshold)

    names = dict(zip(mapping["string_id"], mapping["string_preferred_name"]))
    locus_to_ko = {g.split(":", 1)[1]: dr_ko_links[g] & human_kos for g in projected_genes}
    graph = nx.Graph(variant="ko_projected", string_taxon_id=dr_taxon, threshold=threshold)
    for sid in mapped_ids:
        locus = sid.split(".", 1)[1]
        graph.add_node(sid, string_preferred_name=names.get(sid, locus), is_original_seed=locus in seed_loci,
                       shared_kos=";".join(sorted(locus_to_ko.get(locus, set()))))
    for r in kept.itertuples():
        graph.add_edge(r.node_u, r.node_v, weight=float(r.combined_score), combined_score=float(r.combined_score),
                       experimental_score=float(r.experimental_score))
    nodes = pd.DataFrame([{"string_id": n, **d} for n, d in graph.nodes(data=True)])
    return {
        "graph": graph, "nodes": nodes,
        "summary": {
            "human_seed_kos": len(human_kos), "projected_dr_genes": len(projected_genes),
            "original_seeds": len(seed_loci), "queried": len(query), "mapped_to_string": len(mapped_ids),
            "unmapped": len(unmapped),
            "original_seeds_mapped": int(nodes["is_original_seed"].sum()) if len(nodes) else 0,
            "string_records_retrieved": len(edges), "edges_retained": len(kept),
            "excluded_edges": reason_counts(dropped),
        },
    }


def ko_prior_frame(graph: nx.Graph, hs_profiles: dict, ko_by_dr: dict[str, set[str]]) -> pd.DataFrame:
    """Pseudo 'homology' table (query, hit, score=1) for D. radiodurans/human pairs sharing a KO,
    used as the IsoRank prior for the variant."""
    rows = []
    for dr_id in graph.nodes:
        kos = ko_by_dr.get(dr_id, set())
        if not kos:
            continue
        for hs_id, prof in hs_profiles.items():
            if kos & prof.ko:
                rows.append({"query_id": dr_id, "hit_id": hs_id, "bitscore": 1.0})
    return pd.DataFrame(rows, columns=["query_id", "hit_id", "bitscore"])
