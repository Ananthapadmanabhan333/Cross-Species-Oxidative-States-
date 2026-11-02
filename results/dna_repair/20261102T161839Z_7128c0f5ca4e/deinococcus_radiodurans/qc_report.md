# Data-quality report: Deinococcus radiodurans

Run `20261102T161839Z_7128c0f5ca4e` · STRING 12.5 · edge threshold 0.7 (retrieved at STRING required_score 400)

## Species assignment
| item | value |
|---|---|
| NCBI taxon (configured) | 1299 |
| STRING taxon (expected) | 243230 |
| STRING taxa returned by identifier mapping | [243230] |
| UniProt organism ids among seeds | ['1299', '243230'] |
| all retained records match expected STRING taxon | True |

## Seeds (UniProtKB, GO-defined)
Query: `((go:0006281)) AND (taxonomy_id:1299)`

| item | value |
|---|---|
| UniProtKB entries for organism (context) | 3104 |
| seed records retrieved | 56 |
| seed fraction of organism entries | 0.018 |
| seeds passing validation | 56 |
| seeds mapped to STRING | 52 |
| mapping rate | 0.9286 |

Excluded seed records:
_none_

Excluded mapping records:
- `unmapped_in_string`: 4

## Edges
| item | value |
|---|---|
| functional records retrieved | 346 |
| functional edges retained | 119 |
| physical-subnetwork records retrieved | 38 |
| physical-subnetwork edges retained | 14 |
| physical edges without functional counterpart | 0 |

Excluded functional records:
- `below_threshold`: 227

Excluded physical records:
- `below_threshold`: 24

## Network coverage
| item | value |
|---|---|
| nodes in graph | 52 |
| nodes with >=1 retained edge | 41 |
| isolated nodes (kept, flagged) | 11 |
| fraction of mapped seeds connected | 0.7885 |

## Annotation coverage
| item | value |
|---|---|
| nodes | 52 |
| with_uniprot_accession | 52 (1.0) |
| uniprot_reviewed | 23 (0.4423) |
| with_gene_name | 28 (0.5385) |
| with_protein_name | 52 (1.0) |
| with_string_annotation | 52 (1.0) |
| with_seed_go_term_listed_verbatim | 20 (0.3846) |
| uniprot_xref_disagrees_with_string_mapping | 0 (0.0) |

## Graph summary statistics
| item | value |
|---|---|
| nodes | 52 |
| edges | 119 |
| isolated_nodes | 11 |
| density | 0.089744 |
| connected_components | 12 |
| largest_component_nodes | 41 |
| largest_component_fraction | 0.7885 |
| degree_mean | 4.5769 |
| degree_max | 21 |
| weighted_degree_mean | 3.7662 |
| edge_weight_min | 0.7 |
| edge_weight_median | 0.801 |
| edge_weight_max | 0.999 |
| average_clustering_unweighted | 0.4356 |
| edges_in_physical_subnetwork | 14 |
| edges_with_experiments_channel_score | 43 |

## Limitations
- STRING edges are functional associations (combined evidence channels), not physical interactions. `in_physical_subnetwork` marks pairs in STRING's physical-complex subnetwork, which is itself a scored prediction; `has_experiments_channel_score` reflects STRING's experiments channel, which can include evidence transferred from other organisms. Neither is equivalent to an experimentally confirmed physical interaction (e.g. from BioGRID/IntAct).
- Seeds are defined by GO annotation in UniProtKB, so the network covers only proteins already annotated to the configured terms; poorly annotated proteomes yield smaller seed sets.
- With string.add_nodes = 0 only interactions among seeds are retrieved; associations to non-seed proteins are outside the graph by design.
- STRING represents Deinococcus radiodurans by taxon 243230 (a single strain), not the configured species taxon 1299; results describe that strain's proteome.
- 11 of 52 seed proteins have no association at combined score >= 0.7 with another seed.
- 4 valid seed records were not carried into the graph at the identifier-mapping stage (see excluded_records.tsv).
