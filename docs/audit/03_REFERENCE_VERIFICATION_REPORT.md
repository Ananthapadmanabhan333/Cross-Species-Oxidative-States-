# 3. Reference verification report

All 43 references in `paper.tex` were checked for existence and bibliographic accuracy (authors,
title, journal, volume, pages, year, DOI): 37 against their PubMed record (PMID below) and 6, which
PubMed does not index, against Crossref by DOI. No reference was found to be invented or
mis-attributed.

**Claim support.** Bibliographic accuracy is not the same as support for the sentence citing it.
The curated mutant studies were checked against their PubMed abstracts when the evidence tables
were built. Methods and database references are cited for definitions or data sources. The four
background references on *D. radiodurans* (daly2004, daly2007, daly2009, slade2011) were checked at
abstract level only: daly2004 supports "above 10 kGy" and growth under chronic irradiation; daly2007
supports Mn-dependent protection of proteins. Full-text checks were not made for any reference.

| Key | Verified against | Key | Verified against |
|---|---|---|---|
| daly2004 | PMID 15459345 | maslov2002 | PMID 11988575 |
| daly2007 | PMID 17373858 | pesquita2008 | PMID 18460186 |
| daly2009 | PMID 19172147 | guimera2005 | PMID 15729348 |
| slade2011 | PMID 21372322 | netmf | Crossref 10.1145/3159652.3159706 |
| orthofinder | PMID 31727128 | adamic2003 | Crossref 10.1016/s0378-8733(03)00009-1 |
| isorank | PMID 18725631 | networkx | Crossref 10.25080/tcwv9851 |
| magna | PMID 25792552 | harris2004 | PMID 15361932 |
| guzzi2018 | PMID 28062413 | cho2018 | PMID 30488429 |
| uniprot2025 | PMID 39552041 | markillie1999 | PMID 9882685 |
| go2023 | PMID 36866529 | jeong2015 | PMID 26692481 |
| string2025 | PMID 39558183 | wang2016 | PMID 26948123 |
| intact2022 | PMID 34761267 | daly1994 | PMID 8206827 |
| psicquic | PMID 21716279 | selvam2013 | PMID 23840905 |
| quickgo | PMID 19744993 | funayama1999 | PMID 10556595 |
| kegg2025 | PMID 39417505 | earl2002 | PMID 11807060 |
| panther2022 | PMID 34717010 | xu2010 | PMID 20451472 |
| interpro2025 | PMID 39565202 | gaidamakova2022 | PMID 35012337 |
| yu2007 | PMID 17447836 | jolivet2006 | PMID 16359339 |
| louvain | Crossref 10.1088/1742-5468/2008/10/P10008 | narumi2004 | PMID 15458422 |
| ari | Crossref 10.1007/BF01908075 | devigne2013 | PMID 23403184 |
| bh | Crossref 10.1111/j.2517-6161.1995.tb02031.x | leiden | PMID 30914743 |
| | | langfelder2011 | PMID 21283776 |

## Issues

- **Missing DOIs.** `ari`, `bh`, `adamic2003` and `networkx` have no DOI in the manuscript. The
  verified DOIs are listed above.
- **Year conventions.** `cho2018` and `jeong2015` are cited with their print-issue years (2019 and
  2016); PubMed lists the online years (2018 and 2015). The keys are internal labels and do not
  appear in the PDF.

## Claims that need a citation

Candidates were found in PubMed and checked at abstract level. Each must be read before use; none
has been added to the manuscript.

| Claim in the manuscript | Candidate source | What the abstract supports |
|---|---|---|
| Fpg/MutM (bacteria) and OGG1 (eukaryotes) both initiate 8-oxoguanine repair | Boiteux, Coste & Castaing 2017, *Free Radic Biol Med* 107:179–201, PMID 27903453, doi:10.1016/j.freeradbiomed.2016.11.042 | function in both kingdoms; structural comparison stated as covered but not in the abstract |
| Oxidative-base glycosylases form Endo III (HhH-GPD) and Fpg families | Ide 2001, *Prog Nucleic Acid Res Mol Biol* 68:207–221, PMID 11554298, doi:10.1016/s0079-6603(01)68101-7 | the two subgroups |
| NAD⁺- and ATP-dependent ligases share nucleotidyltransferase and OB-fold domains | Shuman & Lima 2004, *Curr Opin Struct Biol* 14:757–764, PMID 15582400, doi:10.1016/j.sbi.2004.10.006 | shared core of a nucleotidyltransferase domain fused to an OB-fold domain |
| Eukaryotic NER uses proteins unrelated to UvrABC | Hoeijmakers 1991, *J Cell Sci* 100:687–691, PMID 1814927, doi:10.1242/jcs.100.4.687 | "no eukaryotic homologs of the UvrABC proteins have been identified" (as of 1991; a recent source is preferable) |
| Cu/Zn and Mn/Fe SODs are unrelated families | none found yet | author to supply |
| Bacterial interaction data are sparse and biased | none searched | author to supply, or support with the IntAct counts in this paper |
| DdrA, DdrB and PprA have no human ortholog (Introduction) | none searched | author to supply; note that harris2004 describes DdrA as a Rad52 homolog (report 2, B1) |
| SOD1 copper loading needs CCS (feasibility of the proposed experiment) | Culotta et al. 1997, *J Biol Chem* 272:23469–23472, PMID 9295278, doi:10.1074/jbc.272.38.23469 | copper delivery to SOD1 by CCS; yeast without it make inactive SOD1 |
