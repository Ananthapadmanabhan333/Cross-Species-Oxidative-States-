# 5. AI-use disclosure

This declaration describes the role of Claude (Anthropic), used through Claude Code, in this
project. It is based on the repository history and the working sessions. Fourteen of the fifteen
commits carry the trailer `Co-Authored-By: Claude Opus 5.5`; the commits were made under the
researcher's git identity, but their content was produced by Claude.

## Who did what

| Activity | Claude | Researcher (Ananthapadmanabhan DS) | Prof. Pietro Hiram Guzzi |
|---|---|---|---|
| Research question and plan | — | brought the project brief | wrote the research plan (main.tex in Overleaf) |
| Code generation | wrote all pipeline code, tests, scripts, configuration and CI | specified goals and constraints; approved changes | — |
| Data acquisition | issued all API requests (STRING, UniProt, QuickGO, KEGG, IntAct) and cached the responses | ran the sessions on their computer | — |
| Computational analysis | executed every run, including reference and sensitivity runs | directed which analyses to run | — |
| Literature curation | searched PubMed, extracted the curated findings from abstracts, verified reference metadata | — | — |
| Data interpretation | wrote the findings documents and all interpretations in the manuscript | reviewed outputs during the sessions; independent verification in progress | — |
| Manuscript drafting | drafted the full manuscript and the Overleaf progress section, and edited them (including with the avoid-ai-writing skill as an editorial checklist) | chose format, placement and authorship; edited the author name | — |
| Error detection | found and fixed the CPU-dependence bugs; found bugs C1 and C2 in this audit | — | — |
| Verification | wrote the independent audit script and these reports | to re-derive key results personally (pending) | — |

No wet-lab or other experimental work was performed by anyone in this project; all results are
computational.

## Limits that readers should know

- The independent audit was also written by Claude. It re-implements the methods separately from
  the pipeline but is not a human check.
- Claude can make errors that look plausible. Two such errors (C1, C2) survived until this audit.
- The literature findings were curated from abstracts, not full texts.

## Proposed acknowledgment for the manuscript

The current acknowledgment ("prepared with the assistance of Claude ... under the direction of the
authors") understates Claude's role. A more accurate wording, for the authors to approve or adapt:

> Claude (Anthropic), an AI model used through Claude Code, wrote the analysis software, retrieved
> the data, ran all analyses, curated the literature evidence from PubMed abstracts and drafted this
> manuscript, under the direction of the first author. The authors reviewed the results, are
> responsible for the content, and have released the code, cached data and an independent
> re-computation of the reported numbers for verification.

Use this wording only once its claims are true: in particular, "the authors reviewed the results"
should follow the section-by-section review.
