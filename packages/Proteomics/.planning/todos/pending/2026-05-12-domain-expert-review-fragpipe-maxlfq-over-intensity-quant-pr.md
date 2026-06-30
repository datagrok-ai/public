---
created: 2026-05-12T19:23:45.471Z
title: Domain-expert review FragPipe MaxLFQ-over-Intensity quant preference
area: import
files:
  - src/parsers/fragpipe-parser.ts
---

## Problem

The FragPipe parser added during Phase 2 of the `prepare-package-for-claude`
harness loop picks one quant column per sample from FragPipe's
`combined_protein.tsv`, preferring `<sample> MaxLFQ Intensity` over
`<sample> Intensity`, then falling back to `<sample> Razor Intensity`.

That choice was made by a coding agent reasoning from FragPipe documentation
("MaxLFQ is the recommended LFQ quant"), not by a proteomics domain expert
familiar with how the rest of *this* package's pipeline behaves
(normalize -> impute -> differential expression). The agent self-flagged this
during the Phase 2 self-report as a judgment call worth confirming.

If MaxLFQ-first is wrong for the downstream pipeline (e.g. it inflates
missingness, distorts variance for DEqMS / limma, or interacts badly with
quantile / VSN normalization), users will silently get wrong analyses on
FragPipe imports.

## Solution

1. Show a proteomics domain expert the FragPipe parser column selection logic
   and the downstream pipeline (normalization.ts, imputation.ts,
   differential-expression.ts).
2. Confirm one of:
   - **MaxLFQ-first stays.** Document the rationale in
     `src/parsers/fragpipe-parser.ts` as a comment.
   - **Different default** (raw Intensity / Razor first). Update the parser
     and tests in `src/tests/fragpipe-parser.ts`.
   - **Make it user-selectable** at import time, similar to how the generic
     parser surfaces intensity column selection. Adds a small dialog before
     the table opens.
3. Whatever the outcome, capture it in the package conventions section of
   `CLAUDE.md` so the next vendor parser added (TMT, DIA-NN, ...) doesn't
   re-litigate the same decision.
