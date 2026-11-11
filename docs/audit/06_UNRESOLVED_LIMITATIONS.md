# 6. Unresolved scientific limitations

## Must be resolved before submission

1. **Human verification.** The first author has not yet re-derived the key results personally.
   The manuscript does not claim otherwise.
Resolved in v0.4.0: IntAct counts (C1), coverage boundary (C2), accurate AI-use acknowledgment,
and citations for the background claims (report 3); the v0.4.0 manuscript compiles in Overleaf. Resolved in v0.4.1: the KEGG responses are
no longer redistributed (report 1); as a consequence, outside reproduction is byte-exact only
while KEGG serves the snapshot data.

## Limitations of the study design (to be stated, not fixed)

6. **Functional associations, not physical interactions.** STRING edges combine co-expression, text
   mining, genomic context and transferred evidence. There are no experimentally determined physical
   interactions for the bacterial seed proteins in IntAct, BioGRID or MINT.
7. **Orthology transfer inside STRING.** Edge conservation between orthologs is partly built into
   the input scores, so homology-anchored alignment results are not independent evidence of
   conserved wiring.
8. **Unequal coverage.** 16 vs 382 proteins for oxidative stress. GO annotation depth, not only
   biology, shapes the networks.
9. **No ground truth for analogy.** Candidate analogues cannot be evaluated for accuracy; the
   ortholog proxy is evaluable only for DNA repair (16 positives; 14 for role similarity).
10. **A priori weights and thresholds.** The score weights and category thresholds are chosen,
    not learned. Top-20 lists depend moderately on the weights.
11. **Limited statistical resolution.** 100 rewired samples bound p at about 0.01; nine alignment
    tests are not adjusted for multiplicity.
12. **Curated evidence is selective.** 12 bacterial proteins, abstract level, targeted search; human
    perturbation evidence only from GO codes.
13. **Non-protein defences are outside the model.** The Mn²⁺–metabolite antioxidants suggested by
    the sodA studies cannot be represented in a protein network.
14. **No experimental validation.** All proposed analogues are hypotheses. The SOD1 complementation
    proposal must also deal with SOD1's dependence on the copper chaperone CCS.

## Parts of the research plan not done

15. Graphlet degree vectors; Leiden communities; module-preservation statistics; MAGNA++ or other
    aligners; cross-species embeddings; transcriptomic overlays; the human-genetics cross-check.
