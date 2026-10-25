# Data-quality report: Deinococcus radiodurans

Run `20261025T174738Z_3dda43d0e49e` · STRING 12.5 · edge threshold 0.7 (retrieved at STRING required_score 400)

## Species assignment
| item | value |
|---|---|
| NCBI taxon (configured) | 1299 |
| STRING taxon (expected) | 243230 |
| STRING taxa returned by identifier mapping | [243230] |
| UniProt organism ids among seeds | ['243230'] |
| all retained records match expected STRING taxon | True |

## Seeds (UniProtKB, GO-defined)
Query: `((go:0006979) OR (go:0098869) OR (go:0016209)) AND (taxonomy_id:1299)`

| item | value |
|---|---|
| UniProtKB entries for organism (context) | 3104 |
| seed records retrieved | 16 |
| seed fraction of organism entries | 0.0052 |
| seeds passing validation | 16 |
| seeds mapped to STRING | 16 |
| mapping rate | 1.0 |

Excluded seed records:
_none_

Excluded mapping records:
_none_

## Edges
| item | value |
|---|---|
| functional records retrieved | 38 |
| functional edges retained | 11 |
| physical-subnetwork records retrieved | 4 |
| physical-subnetwork edges retained | 4 |
| physical edges without functional counterpart | 0 |

Excluded functional records:
- `below_threshold`: 27

Excluded physical records:
_none_

## Network coverage
| item | value |
|---|---|
| nodes in graph | 16 |
| nodes with >=1 retained edge | 7 |
| isolated nodes (kept, flagged) | 9 |
| fraction of mapped seeds connected | 0.4375 |

## Annotation coverage
| item | value |
|---|---|
| nodes | 16 |
| with_uniprot_accession | 16 (1.0) |
| uniprot_reviewed | 4 (0.25) |
| with_gene_name | 4 (0.25) |
| with_protein_name | 16 (1.0) |
| with_string_annotation | 16 (1.0) |
| with_seed_go_term_listed_verbatim | 6 (0.375) |
| uniprot_xref_disagrees_with_string_mapping | 0 (0.0) |

## Graph summary statistics
| item | value |
|---|---|
| nodes | 16 |
| edges | 11 |
| isolated_nodes | 9 |
| density | 0.091667 |
| connected_components | 11 |
| largest_component_nodes | 5 |
| largest_component_fraction | 0.3125 |
| degree_mean | 1.375 |
| degree_max | 4 |
| weighted_degree_mean | 1.1591 |
| edge_weight_min | 0.733 |
| edge_weight_median | 0.806 |
| edge_weight_max | 0.986 |
| average_clustering_unweighted | 0.3125 |
| edges_in_physical_subnetwork | 4 |
| edges_with_experiments_channel_score | 9 |

## Limitations
- STRING edges are functional associations (combined evidence channels), not physical interactions. `in_physical_subnetwork` marks pairs in STRING's physical-complex subnetwork, which is itself a scored prediction; `has_experiments_channel_score` reflects STRING's experiments channel, which can include evidence transferred from other organisms. Neither is equivalent to an experimentally confirmed physical interaction (e.g. from BioGRID/IntAct).
- Seeds are defined by GO annotation in UniProtKB, so the network covers only proteins already annotated to the configured terms; poorly annotated proteomes yield smaller seed sets.
- With string.add_nodes = 0 only interactions among seeds are retrieved; associations to non-seed proteins are outside the graph by design.
- STRING represents Deinococcus radiodurans by taxon 243230 (a single strain), not the configured species taxon 1299; results describe that strain's proteome.
- 9 of 16 seed proteins have no association at combined score >= 0.7 with another seed.
