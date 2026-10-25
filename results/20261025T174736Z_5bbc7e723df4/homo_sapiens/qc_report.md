# Data-quality report: Homo sapiens

Run `20261025T174736Z_5bbc7e723df4` · STRING 12.5 · edge threshold 0.7 (retrieved at STRING required_score 400)

## Species assignment
| item | value |
|---|---|
| NCBI taxon (configured) | 9606 |
| STRING taxon (expected) | 9606 |
| STRING taxa returned by identifier mapping | [9606] |
| UniProt organism ids among seeds | ['9606'] |
| all retained records match expected STRING taxon | True |

## Seeds (UniProtKB, GO-defined)
Query: `((go:0006979) OR (go:0098869) OR (go:0016209)) AND (organism_id:9606) AND (reviewed:true)`

| item | value |
|---|---|
| UniProtKB entries for organism (context) | 20431 |
| seed records retrieved | 385 |
| seed fraction of organism entries | 0.0188 |
| seeds passing validation | 385 |
| seeds mapped to STRING | 382 |
| mapping rate | 0.9922 |

Excluded seed records:
_none_

Excluded mapping records:
- `collapsed_to_same_string_id`: 1
- `unmapped_in_string`: 2

## Edges
| item | value |
|---|---|
| functional records retrieved | 5550 |
| functional edges retained | 2210 |
| physical-subnetwork records retrieved | 1979 |
| physical-subnetwork edges retained | 1152 |
| physical edges without functional counterpart | 0 |

Excluded functional records:
- `below_threshold`: 3340

Excluded physical records:
- `below_threshold`: 827

## Network coverage
| item | value |
|---|---|
| nodes in graph | 382 |
| nodes with >=1 retained edge | 312 |
| isolated nodes (kept, flagged) | 70 |
| fraction of mapped seeds connected | 0.8168 |

## Annotation coverage
| item | value |
|---|---|
| nodes | 382 |
| with_uniprot_accession | 382 (1.0) |
| uniprot_reviewed | 382 (1.0) |
| with_gene_name | 382 (1.0) |
| with_protein_name | 382 (1.0) |
| with_string_annotation | 382 (1.0) |
| with_seed_go_term_listed_verbatim | 155 (0.4058) |
| uniprot_xref_disagrees_with_string_mapping | 0 (0.0) |

## Graph summary statistics
| item | value |
|---|---|
| nodes | 382 |
| edges | 2210 |
| isolated_nodes | 70 |
| density | 0.030369 |
| connected_components | 79 |
| largest_component_nodes | 295 |
| largest_component_fraction | 0.7723 |
| degree_mean | 11.5707 |
| degree_max | 90 |
| weighted_degree_mean | 10.3283 |
| edge_weight_min | 0.7 |
| edge_weight_median | 0.92 |
| edge_weight_max | 0.999 |
| average_clustering_unweighted | 0.3606 |
| edges_in_physical_subnetwork | 1152 |
| edges_with_experiments_channel_score | 1659 |

## Limitations
- STRING edges are functional associations (combined evidence channels), not physical interactions. `in_physical_subnetwork` marks pairs in STRING's physical-complex subnetwork, which is itself a scored prediction; `has_experiments_channel_score` reflects STRING's experiments channel, which can include evidence transferred from other organisms. Neither is equivalent to an experimentally confirmed physical interaction (e.g. from BioGRID/IntAct).
- Seeds are defined by GO annotation in UniProtKB, so the network covers only proteins already annotated to the configured terms; poorly annotated proteomes yield smaller seed sets.
- With string.add_nodes = 0 only interactions among seeds are retrieved; associations to non-seed proteins are outside the graph by design.
- 70 of 382 seed proteins have no association at combined score >= 0.7 with another seed.
- 3 valid seed records were not carried into the graph at the identifier-mapping stage (see excluded_records.tsv).
