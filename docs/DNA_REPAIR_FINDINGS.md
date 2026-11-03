# DNA-repair findings: D. radiodurans vs human

Reference run: `results/dna_repair/phase2/20261103T112433Z_73788e2787df/` (generated report:
`technical_report.md`). Seed-sensitivity runs: `20261103T114216Z_0e41d8f5b4be` (seed 1) and
`20261103T114851Z_41c672396042` (seed 987654). Input: Phase 1 run `20261102T161839Z_7128c0f5ca4e`
(`config/phase1_dna_repair.yaml`; STRING 12.5), re-verified by checksum before analysis.

The code and analysis parameters are identical to the oxidative-stress analysis
([PHASE2_FINDINGS.md](PHASE2_FINDINGS.md)); see [PHASE2_METHODS.md](PHASE2_METHODS.md). Every
number below comes from the run outputs. **These are computational candidate relationships. None
demonstrates biological equivalence, and none has been tested experimentally here.**

## 1. Data and coverage

Seeds are UniProtKB entries annotated with GO:0006281 (DNA repair), including descendant terms.
This gives a far better-balanced comparison than oxidative stress:

| | *D. radiodurans* (STRING 243230) | *H. sapiens* |
|---|---|---|
| seed records / mapped to STRING | 56 / 52 | 432 / 426 |
| edges at combined score ≥ 0.7 | 119 | 5,285 |
| isolated proteins | 11 | 23 |
| largest component | 41 | 400 |
| proteins with experimental-grade evidence (GO or curated) | 10 | 258 |

- **Unmapped seeds are reported, not substituted.** Four *D. radiodurans* and six human seeds have
  no STRING entry.
- **The bacterial seed set is more connected than size alone predicts.** Random 52-protein human
  subsets have a mean of 77 edges (95th percentile 116); the bacterial set has 119, more than
  95.3 % of the subsets. STRING's PPI-enrichment test is significant for both (the seeds share GO
  annotation, so this is expected).
- **Enlarging the network via orthology adds almost nothing.** The KO-projected coverage variant
  adds a single protein, so the GO-defined bacterial DNA-repair network is already near the
  coverage that orthology can provide.

## 2. Bacterial modules recover known repair pathways

Louvain at resolution 1.0 gives 5 non-singleton modules (mean ARI across 50 seeds 0.76; modularity
0.39). All 5 have FDR-significant enrichment against the seed background:

| module | members | top enrichment |
|---|---|---|
| base-excision repair | ung, nth, DR_0354, DR_2074, DR_2285, DR_2438, DR_2584, … | base-excision repair (FDR 3 × 10⁻⁴) |
| RecFOR recombination | recF, recO, recR, recG, recN, radA, ruvA, ruvC | DNA recombination; homologous recombination (dra03440) |
| nucleotide-excision repair | uvrA, uvrB, uvrC, mfd, ruvB, DR_A0188, DR_0856 | excinuclease repair complex; NER (dra03420) |
| replication / mismatch repair | polA, ligA, mutL, mutS2, mutM, ssb, sbcC, sbcD, … | DNA replication (FDR 0.023) |
| radiation response | recA, lexA, ddrA, ddrB, pprA | response to ionising radiation; response to desiccation |

The human network (9 non-singleton modules, mean ARI 0.64) separates into modules for the
Fanconi anaemia pathway and D-loop resolution, DNA replication, base-excision repair,
nucleotide-excision repair, translesion synthesis, and SMC5/6 cohesion, plus regulatory and
proteostasis groupings.

**Module correspondence is weakly discriminating.** The bacterial base-excision module matches the
human base-excision module best (core-term similarity 0.28). But the bacterial excision-repair and
recombination modules *also* match the human base-excision module best, not the human
nucleotide-excision or recombination modules. Shared generic repair terms dominate the module
profiles, so this layer says "DNA repair resembles DNA repair" more than it maps pathways one to one.

## 3. Network alignment

| prior | edges conserved / 119 | rewired-null mean | p vs rewired null |
|---|---|---|---|
| none (α = 1) | 114 | 110.0 | 0.08 (100 samples) |
| homology (α = 0.8) | 87 | 77.7 | 0.02 |
| homology (α = 0.5) | 89 | 75.3 | 0.01 |

- **Topology alone shows no significant conservation** (p = 0.08; KO-projected variant
  p = 0.10). A dense human network lets most small graphs align well, so 114 of 119 edges are
  conserved even for degree-matched random graphs (null mean 110.0).
- **Anchored on homologs, the bacterial association structure is conserved more than chance.**
  This holds against degree-preserving rewired graphs (p = 0.01–0.02; KO variant p = 0.005). With
  100 rewired samples the smallest attainable p is 1/101 ≈ 0.01 (1/201 with 200), so these values
  sit at the resolution limit of the test.
- **Correction (2026-11-03).** An earlier version of this analysis reported p = 0.05 for topology
  alone (KO variant 0.08). Near-ties in IsoRank's greedy matching were then decided by
  floating-point rounding, which varies between CPUs; ties are now broken deterministically by
  node order, and the result is identical on different hardware.
- **This is not independent evidence of conserved wiring.** STRING transfers association evidence
  between orthologs of different organisms, so part of this conservation is built into the input
  scores. It is consistent with conserved pathway structure but cannot demonstrate it.

## 4. Candidate relationships

Eligible / all pairs per category:

| category | eligible / all |
|---|---|
| reciprocal-best-hit orthologs | 13 / 14 |
| KO-only orthologs | 2 / 2 |
| homologs | 267 / 302 |
| distant homologs (common fold only) | 1,083 / 1,157 |
| candidate analogues | 90 / 90 |
| shared function annotation | 305 / 328 |
| similar pathway role | 1,006 / 1,006 |
| speculative or unsupported | the remainder |

### 4.1 Orthologs and homologs: the conserved core

- **Reciprocal best hits include well-known ortholog pairs:**

  | D. radiodurans | human | rank |
  |---|---|---|
  | nth | NTHL1 | 2 |
  | DR_0354 | APEX1 | 4 |
  | DR_0715 | TDG | 6 |
  | DR_2285 | MUTYH | 7 |
  | DR_0467 | POLL | 11 |
  | recA | DMC1 | 30 |
  | DR_2074 | MPG | 33 |
  | ung | UNG | 82 |
  | polA | POLN | 91 |
  | mutL | MLH1 | 127 |

- **Homologs** include mutM with NEIL1/NEIL3 (ranks 3–5; the Fpg/Nei family), recA with RAD51
  (rank 169), DR_1289 (RecQ) with WRN (rank 40), and ddrA with RAD52 (shared Pfam domain). The last
  is consistent with the PubMed-recorded description of DR0423/DdrA as a homologue of eukaryotic
  Rad52 (Harris et al. 2004, [doi:10.1371/journal.pbio.0020304](https://doi.org/10.1371/journal.pbio.0020304)).
- **Two reciprocal best hits should not be read as orthologs.** ligA–RFC1 (bitscore 65.5; rank 478)
  most likely shares only a BRCT domain: ligA is a DNA ligase, RFC1 a clamp loader. The evidence
  integration ranks it low, but its category label is misleading. In addition, 9 of the 14
  reciprocal best hits fall in different KEGG KO groups. These include true orthologs such as
  mutL–MLH1, because KEGG often splits bacterial and eukaryotic members of one family. This is
  recorded as an uncertainty, not a contradiction.

### 4.2 Candidate analogues: convergent DNA glycosylases

These pairs have no homology signal at any tier (best hit, KO, PANTHER, Pfam, InterPro, common
fold) and share catalytic activity with independent support.

- **Oxidised-base glycosylases with different folds:**

  | D. radiodurans | human | rank |
  |---|---|---|
  | mutM (Fpg family) | NTHL1 (HhH-GPD family) | 26 |
  | mutM | OGG1 (HhH-GPD) | 47 |
  | nth (HhH-GPD) | NEIL3 (Fpg/Nei family) | 31 |
  | nth | NEIL1 | 52 |

  All share EC 4.2.99.18 (DNA-(apurinic or apyrimidinic site) lyase). The bacterial Fpg/MutM and
  the eukaryotic OGG1 both excise 8-oxoguanine, which makes them the clearest example of
  convergent function the pipeline recovers.
- **3-methyladenine glycosylases.** DR_2584 (EC 3.2.2.21) vs human MPG ranks 62.
- **Weak analogues, driven by broad EC numbers.** DR_0467 pairs with POLK, POLN and POLG (ranks
  79–118) on EC 2.7.7.7 (any DNA-directed DNA polymerase). mutM and nth pair with POLB on AP-lyase
  activity. These share an activity class, not a mechanism.
- **Correction: the ligases are not analogues.** I had expected ligA (NAD⁺-dependent) vs LIG1
  (ATP-dependent) to be the textbook case of unrelated isofunctional enzymes. The pipeline
  classifies the pair as `distant_homolog` (shared common-fold signature; rank 514), which is
  consistent with both ligase families sharing adenylation and OB-fold domains. The data were right
  and my prior framing was not.

### 4.3 Pathway-level analogy: bacterial vs eukaryotic excision repair

- **NER.** uvrA–XPC is `similar_pathway_role` (shared KEGG NER map, no shared function annotation;
  rank 2,771), and uvrB–ERCC2 is a distant homolog (both helicases; rank 621). The bacterial UvrABC
  and eukaryotic XPC/TFIIH systems perform the same pathway with largely different proteins. The
  pipeline places that analogy at the pathway level, where it belongs, not at the level of single
  enzyme pairs.
- **Recombination mediators.** recO–RAD52 is `similar_pathway_role` (rank 176), consistent with
  both acting in recombination without detectable homology.

### 4.4 Deinococcus-specific proteins and the curated evidence

The table `config/curated_evidence_dna_repair.tsv` holds 13 PubMed-verified findings for 7
proteins, curated at abstract level. All of the following is according to PubMed:

- **recA.** The rec30 mutant repairs few or no double-strand breaks (Daly et al. 1994,
  [doi:10.1128/jb.176.12.3508-3517.1994](https://doi.org/10.1128/jb.176.12.3508-3517.1994)).
  RecA's co-protease activity, not its recombination activity, may underlie radiation resistance
  (Satoh et al. 2002, [doi:10.1093/oxfordjournals.jbchem.a003066](https://doi.org/10.1093/oxfordjournals.jbchem.a003066)).
- **ddrA and ddrB.** Mutants are radiation- and/or mitomycin-sensitive (Harris et al. 2004; Selvam
  et al. 2013, [doi:10.1371/journal.pone.0069007](https://doi.org/10.1371/journal.pone.0069007)).
  A ddrB/recA double mutant is defective in fragment rejoining (Xu et al. 2010,
  [doi:10.1016/j.dnarep.2010.04.006](https://doi.org/10.1016/j.dnarep.2010.04.006)).
- **recN and uvsE.** Mutants are sensitive to DNA-damaging agents (Funayama et al. 1999,
  [doi:10.1016/s0921-8777(99)00044-0](https://doi.org/10.1016/s0921-8777(99)00044-0)) and to UV in
  a uvrA background (Earl et al. 2002, [doi:10.1128/jb.184.4.1003-1009.2002](https://doi.org/10.1128/jb.184.4.1003-1009.2002)).
- **pprA: mixed evidence.** PprA stimulates DNA end-joining by ligases, and pprA loss reduces
  radiation resistance (Narumi et al. 2004, [doi:10.1111/j.1365-2958.2004.04272.x](https://doi.org/10.1111/j.1365-2958.2004.04272.x)).
  However, ΔpprA repairs double-strand breaks efficiently, and its radiosensitivity is attributed to
  defective DNA segregation and cell division (Devigne et al. 2013,
  [doi:10.1016/j.dnarep.2013.01.004](https://doi.org/10.1016/j.dnarep.2013.01.004)). The pipeline
  records both and flags the contradiction. pprA has no supported human counterpart (best pairs are
  `speculative`), which fits its genus-specific status.
- **lexA: contradicting evidence.** Neither lexA deletion nor non-cleavable LexA changes survival
  or double-strand-break repair after gamma irradiation (Jolivet et al. 2006,
  [doi:10.1111/j.1365-2958.2005.04946.x](https://doi.org/10.1111/j.1365-2958.2005.04946.x)). Its GO
  "DNA repair" annotation is not borne out by mutant data in this organism, and its pairs carry
  that contradiction.

## 5. Experimentally confirmed physical interactions (IntAct)

| | *D. radiodurans* | *H. sapiens* |
|---|---|---|
| physical interaction pairs between seeds | 0 | 599 |
| STRING edges (≥ 0.7) confirmed | 0 of 119 | 431 of 5,285 (8.2 %) |
| STRING physical-subnetwork edges confirmed | 0 of 14 | 395 of 2,393 (16.5 %) |
| IntAct pairs absent from the STRING graph | – | 168 (111 absent even at score ≥ 0.4) |

Most human STRING edges are functional associations without a deposited physical interaction.
**No experimentally confirmed physical interactions exist for the bacterial repair proteins in
IntAct, BioGRID or MINT**, so the bacterial network is entirely functional-association based.

## 6. Robustness and evaluation

| check | result |
|---|---|
| Ranking weights (500 Dirichlet samples) | median Spearman 0.98 (5 %: 0.93); median top-20 overlap 0.70 (5 %: 0.50) |
| Leave one component out | top-20 overlap falls most without `homology` (0.35) and `ec` (0.55) |
| Similarity measure | top-20 overlap with Ochiai: 0.85 (simGIC), 0.50 (min); candidate analogues 67 / 90 / 108 |
| STRING threshold 0.4–0.9 | top-20 overlap with the 0.7 ranking 0.60–1.00 |
| Random seed (two alternatives) | top-20 overlap 1.0; score Spearman ≥ 0.9988 |
| Python hash seed (offline) | all 42 result tables byte-identical |
| Embeddings (NetMF) | link-prediction AUC 0.879 (human) and 0.728 (*D. radiodurans*), below Adamic–Adar (0.952 / 0.931). Not used. |
| Proxy (16 ortholog pairs; evaluable) | `role` (topology only, no leakage) AUC 0.46: topology does not recover orthologs. Annotation components 0.83–0.97 (inflated by annotation transfer) |

**The top of this ranking is less stable than for oxidative stress.** Without the homology
component, only 35 % of the top 20 remain, because the top of the list is mostly orthologs. Read
individual ranks within their reported 5–95 % intervals.

## 7. Limitations

- **The curated evidence is selective.** It covers 7 bacterial proteins, at abstract level, from a
  targeted search. Findings whose abstracts give no locus tag (e.g. uvrA1/uvrA2) were not used.
- **Broad EC numbers inflate some "analogues"** (e.g. EC 2.7.7.7 for all DNA polymerases).
- **"No detected homology" cannot exclude very remote homology**, although common-fold signatures
  are now checked.
- **STRING edges are functional associations,** partly transferred between organisms through
  orthology. This affects both alignment p-values.

## 8. Proposed next steps

1. **Test the convergent glycosylase pairs.** Cross-complementation of a *D. radiodurans* mutM
   deletion with human OGG1, and of nth with NEIL1, under oxidative or ionising damage, measuring
   8-oxoG removal.
2. **Characterise the Deinococcus-specific factors against human end-processing proteins**
   (pprA, ddrA, ddrB). The current data give no supported human counterpart for pprA, and its
   mutant phenotype points to a segregation role.
3. **Extend curation to full texts and to human mutant phenotypes beyond GO.**
