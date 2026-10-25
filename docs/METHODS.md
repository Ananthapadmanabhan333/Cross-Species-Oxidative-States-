# Methods, data model and QC rules (Phase 1)

## 1. Seed definition

Seeds are UniProtKB entries of the organism annotated with any configured GO term
(`seeds.go_terms`), queried as:

```
((go:0006979) OR (go:0098869) OR (go:0016209)) AND (<uniprot_taxon_query>) [AND (reviewed:true)]
```

UniProt's `go:` field may also match entries annotated to descendant terms. The column
`matched_seed_go_terms` lists the configured terms that appear verbatim on each entry. An empty
value means the entry matched through a descendant term.

`uniprot_taxon_query` is `taxonomy_id:1299` for *D. radiodurans*, which includes strains below
the species node, and `organism_id:9606` for human. An additional request reads UniProt's
`X-Total-Results` header for the organism, which gives the denominator for seed coverage.

## 2. Identifier mapping

The STRING `get_string_ids` endpoint is called with all seed accessions, `species=<ncbi_taxon_id>`,
`limit=1` and `echo_query=1`. A batch answer contains **at most one row per STRING protein**, so
when two accessions map to the same protein, only one appears. Each accession missing from the
batch answer is therefore re-queried on its own (`string.max_individual_remap`, default 50). That
separates true gaps (`unmapped_in_string`) from accessions that share a protein
(`collapsed_to_same_string_id`). Mapping rows record `mapping_pass = batch | single`.

The returned STRING taxon must equal `string_taxon_id`. For *D. radiodurans*, STRING resolves 1299
to the R1 strain, 243230. The pipeline also compares UniProt's own STRING cross-reference
(`xref_agrees`) with the API result; disagreements are counted but do not cause exclusion.

## 3. Network retrieval

For each `network_type` in `string.network_types` (`functional`, `physical`), the pipeline calls
STRING `network` with the mapped seeds, `species=<string_taxon_id>`,
`required_score=<retrieval_required_score>` (default 400) and `add_nodes` (default 0). The
analysis threshold `network.edge_score_threshold` (default 0.7) is applied afterwards. The config
validator refuses a threshold below the retrieval score, because those records were never
downloaded.

## 4. QC rules

A record receives the **first** reason it fails, in the order listed, so the per-reason counts sum
to the number of excluded records. Every retrieved record is either retained or listed in
`excluded_records.tsv`.

| stage | reasons (in order) |
|---|---|
| seed | `missing_accession`, `invalid_uniprot_accession`, `unexpected_organism` (UniProt organism not in {ncbi_taxon_id, string_taxon_id}), `duplicate_accession` |
| mapping | `unmapped_in_string`, `invalid_string_id`, `string_taxon_mismatch`, `collapsed_to_same_string_id` |
| edge | `missing_identifier`, `invalid_string_id`, `species_mismatch` (id prefix or reported taxon), `self_loop`, `missing_score`, `score_out_of_range`, `endpoint_not_in_node_set` (only when `add_nodes = 0`), `duplicate_edge` (same unordered pair; highest score kept), `below_threshold` (score < threshold; the threshold itself is inclusive) |

Physical-subnetwork edges go through the same checks. A retained physical edge that has no
retained functional counterpart is counted (`physical_without_functional`) and never added to the
graph silently.

## 5. Graph model

- One `networkx.Graph` (undirected) per species. The graph attributes record the species, both
  taxa, the STRING version, the threshold and the run id.
- Nodes are STRING ids. Their attributes are UniProt accession(s), gene and protein names, the
  reviewed status, GO ids, the STRING preferred name and annotation, `is_seed`, `is_isolated` and
  `xref_agrees`. Missing values are omitted rather than filled with placeholders.
- Edges carry `weight = combined_score`, the seven channel scores, `in_physical_subnetwork`,
  `physical_score`, `has_experiments_channel_score`, `evidence_source` and `source_record`, which
  is the row index in the parsed STRING table.
- Isolated seeds are kept by default (`network.keep_isolated_nodes`) and flagged.

## 6. Coverage metrics

- **Network coverage:** seeds retrieved / organism entries, mapping rate, nodes with at least one
  retained edge, isolated nodes.
- **Annotation coverage:** fractions of nodes that have a UniProt accession, Swiss-Prot reviewed
  status, a gene name, a protein name, a STRING annotation other than "annotation not available",
  and a seed GO term listed verbatim.

## 7. Caching and provenance

The cache key is the SHA-256 of (method, URL, sorted parameters). Raw bodies are never
overwritten; `refresh` adds a new timestamped copy. On every cache read, the stored SHA-256 is
re-checked. The manifest lists every response served, with `served_from` set to `network`,
`cache` or `stale_cache`. It also lists all generated files with checksums; the live log directory
is the only exception and is declared as such.

## 8. Phase 2 candidates (not implemented)

- Orthology mapping between the two species (e.g. eggNOG/OrthoDB) for cross-species alignment.
- Experimentally confirmed physical interactions (BioGRID/IntAct) as a separate edge layer.
- DNA-repair networks, to be started only after Phase 1 is accepted.
