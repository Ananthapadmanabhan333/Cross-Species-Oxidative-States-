# Phase 2 findings: oxidative-stress defence networks, D. radiodurans vs human

Reference run: `results/phase2/20261025T203119Z_863b77a88678/` (generated report:
`technical_report.md`). Seed-sensitivity runs: `20261025T203256Z_bdecaafcd1e1` (seed 1) and
`20261025T203431Z_d5645afdbed0` (seed 987654). Input: Phase 1 run `20261025T174736Z_5bbc7e723df4`
(STRING 12.5, UniProt 2026_03), re-verified by checksum before analysis.

This document interprets those outputs. Every number below comes from them. **These are
computational candidate relationships ranked by available evidence; none demonstrates biological
equivalence, and none has been experimentally tested here.**

## 1. What the data can and cannot support

The comparison is strongly limited by coverage:

- The *D. radiodurans* network has **16 proteins and 11 edges** at combined score ≥ 0.7, and 9 of
  the 16 proteins are isolated. The human network has 382 proteins and 2,210 edges.
- Only 2 of the 16 *D. radiodurans* proteins (katA, sodA) have experimental GO evidence. The
  other 14 have electronic-only (IEA; 4 proteins) or phylogenetic-plus-electronic (IBA/IEA;
  10 proteins) GO evidence.
- **The bacterial sparsity reflects seed-set size, not unusually poor connectivity.** Random
  16-protein subsets of the human seeds have a mean of 3.5 edges (95th percentile 9); the bacterial
  seed set has 11, more than 97.7 % of the subsets. STRING's PPI-enrichment test is significant for
  both seed sets at every threshold (*D. radiodurans* at 0.7: 11 observed vs ~0 expected edges,
  p = 1.2 × 10⁻¹⁷). Because the seeds were selected by shared GO annotation, this connectivity is
  expected and is not evidence of a new pathway.

The consequence is that **annotation-level evidence is informative, and topology-level evidence for
*D. radiodurans* is weak**. The analyses below are read with that in mind.

## 2. Topology and modules

- **Human.** The network has one large component (295 of 382 nodes) and moderate clustering (0.36).
  Louvain at resolution 1.0 gives about 85 modules, 15 of them non-singleton (modularity 0.54), with
  moderate stability across 50 seeds (mean pairwise ARI 0.74). Seven modules have FDR-significant
  enrichment against the seed background, among them:
  - an antioxidant/detoxification module of 69 proteins (SOD1–3, CAT, PRDX1–6, GPX2–8, TXNRD1–3,
    glutathione synthesis and transferases, NFE2L2/KEAP1); its top terms are oxidoreductase activity
    (FDR 2 × 10⁻²⁷) and cellular oxidant detoxification;
  - a proteasome module (45 proteins);
  - a mitochondrial inner-membrane module (35 proteins);
  - a peroxisomal-import module (7 proteins);
  - a 133-protein module enriched for broad regulatory terms. Its size and generic terms suggest a
    heterogeneous grouping rather than one mechanism.
- ***D. radiodurans*.** Louvain simply recovers the two non-trivial components: {katA, DR_A0259,
  sodA, sodC, DR_A0202}, the catalase/SOD module, and {DR_1208, DR_1209}, two peroxiredoxin-family
  proteins. That partition is unique, hence ARI = 1. The only significant term (FDR 0.017) is the
  KEGG map "Peroxisome" (dra04146). This is an artefact of ortholog-based KEGG map membership —
  bacteria have no peroxisomes — and is not organelle evidence.
- **Module correspondence.** Both bacterial modules match the human antioxidant module best
  (core-term similarity 0.49 and 0.47). This is descriptive: it shows that the same enzyme classes
  dominate both, not that the modules are organised alike.

## 3. Network alignment: no evidence of a conserved topology

IsoRank with topology only (α = 1) maps 9 of the 11 *D. radiodurans* edges onto human edges.
That beats random one-to-one mappings (p = 0.001) but **not** degree-preserving rewired bacterial
graphs aligned the same way (p = 0.10). The run also did not converge. When the homology prior is
added (α = 0.8 and 0.5), all 5 reciprocal-best-hit pairs are recovered, but edge conservation is no
better than for rewired graphs (p = 0.63 and 0.88).

**Conclusion:** with 11 edges, there is no evidence that the bacterial association topology
matches the human one beyond what its degree sequence and the density of the human network allow.
A dense target network lets almost any small graph be embedded well; the rewired null is what
exposes this.

## 4. Candidate relationships

All 6,112 protein pairs were scored. Category counts (eligible / all): 5 / 5 RBH orthologs,
134 / 168 homologs, 56 / 56 candidate analogues, 150 / 232 shared-function pairs, 41 / 41
similar-pathway pairs; the rest are speculative or unsupported.

### 4.1 Orthologs: the expected backbone

The five reciprocal best hits are:

| D. radiodurans | human | rank | shared EC |
|---|---|---|---|
| katA | CAT | 5 | 1.11.1.6 |
| sodA | SOD2 | 11 | 1.15.1.1 |
| msrB | MSRB3 | 3 | 1.8.4.12 |
| DR_A0202 | SOD1 | 20 | 1.15.1.1 |
| DR_2242 | PRDX2 | 30 | only sub-subclass 1.11.1.- (EC 1.11.1.29 vs 1.11.1.24) |

These are well-known conserved defence enzymes; recovering them is a sanity check on the
pipeline, not a discovery.

### 4.2 Candidate analogues: same reaction, no detected homology

These pairs share catalytic activity (EC) and have independent pathway, GO or experimental
support, but show no homology signal in best-hit search, PANTHER, Pfam or InterPro.

**Superoxide dismutation by two unrelated enzyme families.** *D. radiodurans* sodA is a Mn/Fe-type
SOD (Pfam PF00081/PF02777); human SOD1 and SOD3 are Cu/Zn-type (PF00080). The reverse holds for
the bacterial Cu/Zn proteins sodC and DR_A0202 versus human SOD2, a Mn-SOD. The relevant pairs are:

| pair | rank | 5–95 % rank interval |
|---|---|---|
| sodA–SOD3 | 4 | 1–17 |
| sodA–SOD1 | 33 | – |
| sodC–SOD2 | 44 | – |
| DR_A0202–SOD2 | 46 | – |

This matches the textbook view that Cu/Zn and Mn/Fe SODs are structurally unrelated families
catalysing the same reaction. The pipeline arrives at that classification from the data alone,
which supports the analogue logic. sodA's oxidative-stress role in *D. radiodurans* is backed by
mutant evidence: according to PubMed, katA and sodA insertion mutants are more sensitive to
ionising radiation than the wild type (Markillie et al. 1999, J Bacteriol,
[doi:10.1128/JB.181.2.666-669.1999](https://doi.org/10.1128/JB.181.2.666-669.1999)).

**Peroxide removal by catalase versus peroxiredoxins.** katA (catalase, EC 1.11.1.6) is paired with
human PRDX2–5 (EC 1.11.1.24) at ranks 42–65. They share only the EC sub-subclass 1.11.1 (peroxidases),
have different folds, and show no detected homology. They are functional analogues for H₂O₂
detoxification at the level of the reaction class only: the chemistry differs (dismutation versus
thiol-based reduction), so this is a weaker analogy than the SOD case. According to PubMed, katA
activity in *D. radiodurans* is linked to H₂O₂ tolerance (Li et al. 2013, PLoS One,
[doi:10.1371/journal.pone.0054420](https://doi.org/10.1371/journal.pone.0054420)).

**Cytoglobin as a non-canonical SOD.** sodA–CYGB ranks 26; the bacterial Cu/Zn SODs pair with
CYGB at ranks 36–39. The shared EC 1.15.1.1 for CYGB rests on one study, which according to PubMed
reports superoxide dismutation by cytoglobin via iron redox cycling (Zweier et al. 2021, PNAS,
[doi:10.1073/pnas.2105053118](https://doi.org/10.1073/pnas.2105053118)). This is the most
hypothesis-generating candidate and also the least established.

### 4.3 Homologs with shared function

The best-scoring pairs overall are homologs, for example:

- DR_A0202/sodC–SOD3 (shared PANTHER family and EC 1.15.1.1; ranks 1–2, intervals 1–3);
- DR_1982–TXNRD1–3 (EC 1.8.1.9, shared Pfam);
- DR_0846/DR_1208–PRDX1–5 (EC 1.11.1.24, shared Pfam).

They identify conserved enzyme families rather than analogues.

### 4.4 Network-role analogues

Fifty-eight eligible non-homologous pairs with functional or pathway evidence also have role
similarity ≥ 0.75 (e.g. katA–PRDX2/PRDX5 ≈ 0.83). Across all eligible pairs, 663 non-homologous
pairs reach that role similarity. Role similarity rests on 7 connected bacterial proteins and showed
no ability to separate orthologs from other pairs (descriptive AUC 0.39, 3 positives). It is
therefore reported but **not** used as evidence of analogous mechanism.

### 4.5 Proteins without a supported partner

DR_1538, DR_A0145, DR_A0301 and DR_1857 have no EC number. Their best *eligible* pairs score
0.32–0.44 and are GO-only (`shared_function_annotation`, DR_1538–PRDX5) or `unsupported`. Their
higher-scoring pairs (0.52–0.58, e.g. DR_1538–GPX8) fall below the evidence-coverage minimum (0.40
vs 0.45) and are therefore not ranked; this is the eligibility gate working as intended.

The only homology signal for DR_1538 and DR_1857 is a shared InterPro entry with ribosomal protein
S3 (IPR015946, a K-homology-like α/β fold). It is a fold-level signature with no functional bearing
here. This illustrates that the InterPro tier is deliberately broad: it errs towards calling pairs
homologs, which is the conservative direction for analogue claims. These four proteins are the
least characterised bacterial seeds and the natural targets for Phase 3 characterisation.

## 5. Robustness

| perturbation | result |
|---|---|
| Ranking weights (500 Dirichlet samples) | median Spearman 0.96 (5 %: 0.91); median top-20 overlap 0.90 (5 %: 0.80) |
| Leave one component out | the ranking depends most on `ec` (top-20 overlap drops to 0.45 without it) and `literature` (Spearman drops to 0.82); every other component: top-20 overlap ≥ 0.85, Spearman ≥ 0.91 |
| Similarity measure | top-20 overlap with Ochiai: 0.90 (simGIC), 0.75 (min); candidate-analogue counts 32 / 56 / 120 |
| STRING threshold 0.4–0.9 | top-20 overlap with the 0.7 ranking 0.90–1.00 |
| Random seed (two alternative master seeds) | top-20 overlap 1.0, score Spearman ≥ 0.999, category agreement ≥ 99.7 % |
| Offline re-run | all 33 result tables byte-identical |

The candidate-analogue count is sensitive to the similarity measure (32–120). That category
should be read as "pairs worth examining", not as a fixed set.

## 6. Embeddings

NetMF link prediction on the human network (5 held-out splits): AUC 0.933 ± 0.012, slightly below
Adamic–Adar (0.949), Jaccard (0.946) and common neighbours (0.945). Neighbourhood heuristics are at
least as good here, so embeddings add nothing and are not used in the ranking. The bacterial
network (11 edges) is too small to evaluate, and no cross-species embedding comparison was made.

## 7. Evaluation and its limits

There is no ground truth for cross-species analogy, so no accuracy is claimed. The proxy check
(how well each component separates the 5 RBH ortholog pairs from all other pairs) has too few
positives to be evaluable. Its annotation-component AUCs (0.95–0.98) are also inflated, because GO
and EC annotations are often transferred between homologs.

## 8. Limitations specific to these findings

- **Perturbation evidence is undercounted.** The katA/sodA mutant study above is coded IDA in GO,
  not IMP, so the perturbation component records no *D. radiodurans* protein. Missing evidence does
  not lower scores, but the component should not be read as "no phenotype data exist".
- **GO evidence is mostly inferred.** Most bacterial GO evidence is electronic (IEA) or
  phylogenetic (IBA), and IBA annotations are themselves inferred from homologs.
- **"No detected homology" does not exclude remote homology**, although the InterPro check reduces
  this risk (it reclassified, for example, thioredoxin-fold peroxiredoxin pairs as homologs).
- **The seed definition omits known mechanisms.** Seeds are GO-defined, so *D. radiodurans* defences
  that are not protein-centred or not GO-annotated are outside the network — notably Mn²⁺–metabolite
  antioxidant complexes and DNA-repair-linked recovery.
- **STRING edges are functional associations**, not physical interactions.

## 9. Proposed next steps (Phase 3, pending your review)

1. **Increase bacterial coverage before any further topology work.** Options: STRING `add_nodes`
   neighbours, a lower threshold, or orthology-projected seeds (bacterial proteins sharing families
   with human seeds), each reported as a separate, labelled network variant.
2. **Add curated perturbation and interaction evidence**, e.g. mutant phenotypes from the
   *D. radiodurans* literature, and BioGRID/IntAct for physical interactions.
3. **Prioritise experimental tests by evidence and uncertainty.** The strongest analogue pair is
   sodA vs SOD1/SOD3 (cross-family complementation of a *D. radiodurans* sodA deletion under
   H₂O₂/paraquat). The most hypothesis-generating is CYGB as an SOD analogue.

DNA repair remains out of scope until you have reviewed these results.
