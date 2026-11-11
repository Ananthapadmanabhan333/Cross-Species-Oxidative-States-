# KEGG responses (not redistributed)

The pipeline uses two KEGG REST responses: the KEGG Orthology (KO) assignments of every gene of
*D. radiodurans* (`link/ko/dra`) and of human (`link/ko/hsa`). KEGG's terms
(https://www.kegg.jp/kegg/legal.html) do not allow its data to be redistributed, so this
repository keeps only the provenance record of each response (`*.meta.json`: URL, retrieval time,
size and SHA-256), not the data.

Restore the data from KEGG with:

    python scripts/fetch_kegg.py

The script stores a response only if its SHA-256 matches the record, that is, if KEGG still serves
exactly the data the reference runs used (retrieved 2026-10-25). KEGG is updated regularly, and by
2026-11-11 both responses had already changed slightly (dra: 2 lines removed, 1 added; hsa:
8 lines added), so the reference runs can no longer be reproduced byte for byte from KEGG. To run
the analysis on the current KEGG release instead:

    python scripts/fetch_kegg.py --accept-current

KEGG is free for academic use; commercial use needs a licence. Please cite KEGG:

- M. Kanehisa and S. Goto, "KEGG: Kyoto Encyclopedia of Genes and Genomes," *Nucleic Acids Res.*,
  28(1):27-30, 2000. doi:10.1093/nar/28.1.27
- M. Kanehisa *et al.*, "KEGG: biological systems database as a model of the real world,"
  *Nucleic Acids Res.*, 53(D1):D672-D677, 2025. doi:10.1093/nar/gkae909
