# Phase 2 findings: oxidative-stress defence networks, D. radiodurans vs human

Reference run: `results/phase2/20261025T222114Z_a6bb58dd4682/` (generated report:
`technical_report.md`). Seed-sensitivity runs: `20261025T222736Z_58b73da6abf5` (seed 1) and
`20261025T222911Z_073fbaf0033f` (seed 987654). Input: Phase 1 run `20261025T174736Z_5bbc7e723df4`
(STRING 12.5, UniProt 2026_03), re-verified by checksum before analysis.

This revision adds three things:

- **KEGG Orthology (KO)** as a curated, cross-kingdom orthology resource;
- **literature-curated perturbation evidence** (`config/curated_evidence.tsv`; every row
  verified against its PubMed abstract), including one finding that contradicts a simple
  protective role;
- **a KO-projected coverage variant** of the bacterial network.

It also fixes two problems:

- **Perturbation scoring:** the earlier rule halved the score when only one protein had evidence,
  which quietly treated missing evidence as negative.
- **Row order:** ties are now broken deterministically, so outputs are byte-identical across
  Python hash seeds.

Every number below comes from the run outputs. **These are computational candidate relationships
ranked by available evidence. None demonstrates biological equivalence, and none has been tested
experimentally here.**

## 1. What the data can and cannot support

- The *D. radiodurans* network has **16 proteins and 11 edges** at combined score ≥ 0.7, and 9 of
  the 16 proteins are isolated. The human network has 382 proteins and 2,210 edges.
- Five of the 16 bacterial proteins now have experimental-grade support (GO experimental codes or
  curated mutant/biochemical studies): katA, sodA, DRA0259, DR1538 and DR0846. The other 11 rest
  on electronic or phylogenetic GO evidence.
- **The sparsity reflects seed-set size, not poor connectivity.** Random 16-protein subsets of the
  human seeds have a mean of 3.5 edges (95th percentile 9); the bacterial set has 11, more than
  97.7 % of the subsets. STRING's PPI-enrichment test is significant for both seed sets at every
  threshold (*D. radiodurans* at 0.7: 11 vs ~0 expected edges, p = 1.2 × 10⁻¹⁷). The seeds share GO
  annotation, so this connectivity is expected.
- **Enlarging coverage via orthology adds little.** The KO-projected variant takes every bacterial
  gene in a KO group carried by a human seed protein. Of 345 human seed KOs, only 16 *D. radiodurans*
  genes match (original seeds included). The variant has 27 proteins (11 new), 16 edges, and
  14 isolated proteins, so the bacterial network stays sparse.

Annotation-level evidence is therefore informative, while topology-level evidence for
*D. radiodurans* remains weak.

## 2. Topology and modules

- **Human.** The network has one large component (295 of 382 nodes). Louvain at resolution 1.0
  gives about 85 modules, 15 of them non-singleton (modularity 0.54), with moderate stability across
  50 seeds (mean ARI 0.74). Seven modules have FDR-significant enrichment against the seed
  background:
  - an antioxidant/detoxification module of 69 proteins (SOD1–3, CAT, PRDX1–6, GPX2–8, TXNRD1–3,
    glutathione synthesis and transferases, NFE2L2/KEAP1);
  - a proteasome module (45 proteins);
  - a mitochondrial inner-membrane module (35 proteins);
  - a peroxisomal-import module (7 proteins);
  - a heterogeneous 133-protein module with broad regulatory terms.
- ***D. radiodurans*.** Louvain recovers the two non-trivial components: {katA, DRA0259, sodA,
  sodC, DRA0202} and {DR1208, DR1209}, hence ARI = 1. Its only significant term, the KEGG map
  "Peroxisome" (dra04146, FDR 0.017), is an artefact of ortholog-based map membership: bacteria
  have no peroxisomes.
- **Module correspondence.** Both bacterial modules match the human antioxidant module best
  (core-term similarity 0.49 and 0.47). This is descriptive only.

## 3. Network alignment: no evidence of a conserved topology

| network | prior | edges conserved | p vs random mapping | p vs degree-preserving rewired null |
|---|---|---|---|---|
| GO seeds (11 edges) | none (α = 1) | 9 | 0.001 | 0.10 (100 samples) |
| GO seeds | homology (α = 0.8 / 0.5) | 8 / 7 | 0.001 | 0.63 / 0.88 |
| KO-projected (16 edges) | none (α = 1) | 9 | 0.001 | 1.0 (200 samples) |
| KO-projected | shared KO (α = 0.8) | 13 | 0.001 | 0.02 (200 samples) |

- **Topology alone shows no conservation** beyond what degree sequences and the density of the
  human network allow, whether coverage is low or modestly enlarged.
- **The one low p-value (0.02) comes from the KO-prior alignment of the enlarged network.** It is
  not independent evidence: STRING transfers association evidence between orthologs of different
  organisms, so conservation between KO-orthologous pairs is partly built into the input scores.
- **Beating random mappings (p = 0.001) only shows that the optimiser works.**

## 4. Candidate relationships

All 6,112 pairs were scored. Eligible / all pairs per category:

| category | eligible / all |
|---|---|
| reciprocal-best-hit orthologs | 5 / 5 |
| KO-only orthologs | 3 / 3 |
| homologs | 135 / 165 |
| candidate analogues | 81 / 81 |
| shared function annotation | 152 / 207 |
| similar pathway role | 41 / 41 |
| speculative or unsupported | the remainder |

### 4.1 Orthologs: the expected backbone

The five reciprocal best hits, with rank and 5–95 % rank interval under weight perturbation:

| D. radiodurans | human | rank | interval |
|---|---|---|---|
| katA | CAT | 7 | 1–16 |
| msrB | MSRB3 | 10 | 1–370 |
| sodA | SOD2 | 16 | 10–84 |
| DR_2242 | PRDX2 | 26 | – |
| DR_A0202 | SOD1 | 42 | – |

KEGG KO adds three orthologs that are not reciprocal best hits: DRA0259–CAT (K03781, rank 9),
msrB–MSRB2 (K07305, rank 33) and sodC–SOD1 (K04565, rank 46). No reciprocal-best-hit pair
contradicts its KO assignment. Recovering these conserved enzymes is a sanity check, not a discovery.

### 4.2 The best-supported bacterial protein: peroxiredoxin Q (DR0846)

The top of the ranking is now DR0846 paired with human PRDX5, PRDX1, PRDX4 and PRDX3 (ranks 1–5;
intervals within 1–10). These are homolog pairs. DR0846 rises because it has direct curated
evidence. According to PubMed, it is a peroxiredoxin Q with peroxidase and chaperone activity, and
a DR0846/DR1998 double deletion is less tolerant of oxidative stress than either single mutant
(Cho et al. 2018, FEBS Lett, [doi:10.1002/1873-3468.13302](https://doi.org/10.1002/1873-3468.13302)).

### 4.3 Candidate analogues: same reaction, no detected homology

These pairs share catalytic activity (EC) and have independent pathway, GO or experimental support,
but show no homology signal in best-hit search, KO, PANTHER, Pfam or InterPro.

**Superoxide dismutation by two unrelated enzyme families.** Mn/Fe SOD sodA pairs with Cu/Zn SOD3
(rank 19, interval 9–74) and SOD1 (rank 123). The bacterial Cu/Zn proteins sodC and DRA0202 pair
with human Mn-SOD2 (ranks 81 and 86). This matches the established view that Cu/Zn and Mn/Fe SODs
are structurally unrelated families catalysing the same reaction.

**The sodA evidence is mixed, and the pipeline now says so.** According to PubMed, sodA insertion
mutants are more sensitive to ionising radiation (Markillie et al. 1999, J Bacteriol,
[doi:10.1128/JB.181.2.666-669.1999](https://doi.org/10.1128/JB.181.2.666-669.1999)). However,
ΔMnSOD mutants showed *increased* gamma-radiation survival while becoming sensitive to atmospheric
oxygen during desiccation. The authors attribute this to small-molecule Mn antioxidants that can
supplant MnSOD (Gaidamakova et al. 2022, mBio,
[doi:10.1128/mbio.03394-21](https://doi.org/10.1128/mbio.03394-21)). sodA therefore scores 0.5 on
perturbation, and every sodA pair lists this as contradictory evidence. This is why sodA–SOD3 fell
from rank 4 to rank 19. **The finding matters for the research question:** in *D. radiodurans*,
the functional analogue of human enzymatic SOD may be partly non-protein (Mn²⁺–metabolite
complexes), which a protein network cannot represent.

**Peroxide removal by catalase versus peroxiredoxins.** katA and DRA0259 (catalases, EC 1.11.1.6)
pair with human PRDX2–5 (EC 1.11.1.24) at ranks 73–99. They share only the EC sub-subclass 1.11.1,
have different folds, and their chemistry differs (dismutation vs thiol-based reduction). This is an
analogy of reaction class only. According to PubMed, Δdr1998 (katA) shows the greatest loss of H₂O₂
resistance of the three catalase mutants (Jeong et al. 2015, BBRC,
[doi:10.1016/j.bbrc.2015.12.017](https://doi.org/10.1016/j.bbrc.2015.12.017)), and disruption of
DRA0259 reduces antioxidant activity (Wang et al. 2016, Int J Radiat Biol,
[doi:10.3109/09553002.2016.1150618](https://doi.org/10.3109/09553002.2016.1150618)).

**Cytoglobin as a non-canonical SOD.** sodA–CYGB ranks 70 (interval 35–449). CYGB's EC 1.15.1.1
rests on one study: according to PubMed, cytoglobin dismutates superoxide via iron redox cycling
(Zweier et al. 2021, PNAS, [doi:10.1073/pnas.2105053118](https://doi.org/10.1073/pnas.2105053118)).
This is hypothesis-generating and weakly established.

### 4.4 Shared function without shared family: OsmC (DR1538) vs glutathione peroxidases

DR1538 (OsmC) pairs with GPX8, GPX7 and GPX2 (ranks 11–13) on GO function similarity plus its
curated disruption phenotype (Wang et al. 2016, above). No EC number is annotated for DR1538, so
these remain `shared_function_annotation`, not candidate analogues. OsmC-family proteins and GPXs
are unrelated families that both reduce hydroperoxides. **This pair deserves an activity assay**,
because it would become a candidate analogue if a shared catalytic activity were shown.

### 4.5 Network roles and proteins without support

- **Network roles.** 58 eligible non-homologous pairs with functional or pathway evidence also have
  role similarity ≥ 0.75. Role similarity rests on 7 connected bacterial proteins and does not
  separate orthologs from other pairs (descriptive AUC 0.35, 5 positives), so it is reported but
  not used as evidence of mechanism.
- **Proteins without a supported partner.** DRA0145, DRA0301 and DR1857 have no EC number and no
  informative homology; their best eligible pairs are `unsupported` (0.42–0.51). DR1538 and DR1857
  share only a fold-level InterPro signature with ribosomal protein S3 (IPR015946), which is not
  functionally informative. These are the natural targets for Phase 3 characterisation.

## 5. Robustness

| perturbation | result |
|---|---|
| Ranking weights (500 Dirichlet samples) | median Spearman 0.98 (5 %: 0.93); median top-20 overlap 0.80 (5 %: 0.65) |
| Leave one component out | top-20 overlap falls most without `ec` (0.60), `homology` (0.60), `go` (0.65) or `literature` (0.70); Spearman falls most without `perturbation` (0.89) |
| Similarity measure | top-20 overlap with Ochiai: 0.90 (simGIC), 0.60 (min); candidate analogues 65 / 81 / 120 |
| STRING threshold 0.4–0.9 | top-20 overlap with the 0.7 ranking 0.80–1.00 |
| Random seed (two alternatives) | top-50 overlap 1.0, score Spearman ≥ 0.9995, category agreement ≥ 99.7 % |
| Python hash seed (3 alternatives, offline) | all 37 result tables byte-identical |

**The top of the ranking is less stable than before** (top-20 overlap 0.80 vs 0.90). The new
curated and KO evidence is concentrated on a few bacterial proteins, so the weight given to that
evidence matters more. The candidate-analogue set also depends on the similarity measure (65–120).
Read it as "pairs worth examining", not as a fixed list.

## 6. Embeddings

NetMF link prediction on the human network: AUC 0.933 ± 0.012, slightly below Adamic–Adar (0.949),
Jaccard (0.946) and common neighbours (0.945). Embeddings add nothing here and are not used in the
ranking. The bacterial network is too small to evaluate, and no cross-species embedding comparison
is made.

## 7. Evaluation and its limits

There is no ground truth for cross-species analogy, and no accuracy is claimed. The orthology
proxy, which now includes KO orthologs, has 8 positive pairs, still below the 10 needed to be
evaluable. The annotation-component AUCs (0.96–0.99) are inflated, because annotations are often
transferred between homologs.

## 8. Limitations specific to these findings

- **Curated evidence is incomplete.** It covers 5 bacterial proteins, was curated from PubMed
  abstracts only, and comes from a targeted search. A katE2 deletion study (PMID:32366051) was
  excluded because its abstract names no locus tag. Human perturbation evidence still comes only
  from GO evidence codes.
- **Most bacterial GO evidence is inferred** (IEA/IBA). "No detected homology" does not exclude
  remote homology, although KO and InterPro checks reduce this risk.
- **Non-protein defences are outside the network.** These include the Mn²⁺–metabolite antioxidants
  that the sodA evidence points to.
- **STRING edges are functional associations,** partly transferred between organisms through
  orthology.

## 9. Proposed next steps (Phase 3, pending your review)

1. **Test the highest-value pairs experimentally.**
   - sodA vs SOD1/SOD3, by cross-family complementation of a sodA deletion under H₂O₂/paraquat
     and desiccation. The mixed sodA evidence makes the condition choice critical.
   - An activity assay for DR1538 (OsmC) to resolve the GPX comparison.
2. **Represent non-protein antioxidants explicitly** (Mn²⁺–metabolite complexes), e.g. as
   annotated network context. Without them, the bacterial–human comparison is incomplete.
3. **Broaden curated evidence systematically.** Full-text curation, BioGRID/IntAct for physical
   interactions, and human mutant phenotypes beyond GO.

DNA repair remains out of scope until you have reviewed these results.
