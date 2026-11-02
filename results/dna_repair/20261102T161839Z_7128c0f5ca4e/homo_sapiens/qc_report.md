# Data-quality report: Homo sapiens

Run `20261102T161839Z_7128c0f5ca4e` · STRING 12.5 · edge threshold 0.7 (retrieved at STRING required_score 400)

## Species assignment
| item | value |
|---|---|
| NCBI taxon (configured) | 9606 |
| STRING taxon (expected) | 9606 |
| STRING taxa returned by identifier mapping | [9606] |
| UniProt organism ids among seeds | ['9606'] |
| all retained records match expected STRING taxon | True |

## Seeds (UniProtKB, GO-defined)
Query: `((go:0006281)) AND (organism_id:9606) AND (reviewed:true)`

| item | value |
|---|---|
| UniProtKB entries for organism (context) | 20431 |
| seed records retrieved | 432 |
| seed fraction of organism entries | 0.0211 |
| seeds passing validation | 432 |
| seeds mapped to STRING | 426 |
| mapping rate | 0.9861 |

Excluded seed records:
_none_

Excluded mapping records:
- `unmapped_in_string`: 6

## Edges
| item | value |
|---|---|
| functional records retrieved | 13052 |
| functional edges retained | 5285 |
| physical-subnetwork records retrieved | 5450 |
| physical-subnetwork edges retained | 2393 |
| physical edges without functional counterpart | 0 |

Excluded functional records:
- `below_threshold`: 7767

Excluded physical records:
- `below_threshold`: 3057

## Network coverage
| item | value |
|---|---|
| nodes in graph | 426 |
| nodes with >=1 retained edge | 403 |
| isolated nodes (kept, flagged) | 23 |
| fraction of mapped seeds connected | 0.946 |

## Annotation coverage
| item | value |
|---|---|
| nodes | 426 |
| with_uniprot_accession | 426 (1.0) |
| uniprot_reviewed | 426 (1.0) |
| with_gene_name | 426 (1.0) |
| with_protein_name | 426 (1.0) |
| with_string_annotation | 426 (1.0) |
| with_seed_go_term_listed_verbatim | 167 (0.392) |
| uniprot_xref_disagrees_with_string_mapping | 0 (0.0) |

## Graph summary statistics
| item | value |
|---|---|
| nodes | 426 |
| edges | 5285 |
| isolated_nodes | 23 |
| density | 0.058382 |
| connected_components | 25 |
| largest_component_nodes | 400 |
| largest_component_fraction | 0.939 |
| degree_mean | 24.8122 |
| degree_max | 161 |
| weighted_degree_mean | 21.6861 |
| edge_weight_min | 0.7 |
| edge_weight_median | 0.879 |
| edge_weight_max | 0.999 |
| average_clustering_unweighted | 0.4902 |
| edges_in_physical_subnetwork | 2393 |
| edges_with_experiments_channel_score | 4106 |

## Limitations
- STRING edges are functional associations (combined evidence channels), not physical interactions. `in_physical_subnetwork` marks pairs in STRING's physical-complex subnetwork, which is itself a scored prediction; `has_experiments_channel_score` reflects STRING's experiments channel, which can include evidence transferred from other organisms. Neither is equivalent to an experimentally confirmed physical interaction (e.g. from BioGRID/IntAct).
- Seeds are defined by GO annotation in UniProtKB, so the network covers only proteins already annotated to the configured terms; poorly annotated proteomes yield smaller seed sets.
- With string.add_nodes = 0 only interactions among seeds are retrieved; associations to non-seed proteins are outside the graph by design.
- 23 of 426 seed proteins have no association at combined score >= 0.7 with another seed.
- 6 valid seed records were not carried into the graph at the identifier-mapping stage (see excluded_records.tsv).
