---
created: 2026-05-16T00:00:00.000Z
title: DEqMS result rows are misaligned to proteins — silently scrambles the volcano
area: analysis
files:
  - scripts/deqms_de.R (outputResult(fit, coef=2) returns rows sorted by significance, not input order)
  - src/analysis/differential-expression.ts:139 (copyDEResultsToFrame: positional getRawData copy, no protein-ID join)
  - src/analysis/differential-expression.ts:209 (runDeqmsDE)
  - src/tests/ (no test pins DEqMS protein↔stat alignment)
---

## Problem

Selecting **DEqMS** in the Differential Expression dialog produces a
plausible-looking but **scientifically wrong** volcano: every protein is assigned
another protein's log2FC / p-value / significance. The plot keeps the correct
marginal distributions (rabbit-ears shape, a top-left outlier) so the corruption
is **silent** — a user would trust it.

Confirmed empirically (2026-05-16) on the BP DMD-vs-WT client dataset:
log2FC correlation limma~DEqMS = **−0.001** (must be ~1.0 — fold change is the
same deterministic quantity on the same expression matrix regardless of test).
Concordance with the client's own Spectronaut-Candidates hit list at the client
thresholds (Q<0.05, |log2FC|>0.58): limma κ=0.90 / 429-of-465 shared, **DEqMS
κ=0.017 / 35-of-465 shared, only 2 of 126 DMD-enriched shared**. DMD itself —
the obvious #1 hit, correct in limma — is absent from the DEqMS top hits.

Root cause (confirmed in source):
- `deqms_de.R` limma fallback uses `topTable(fit, coef=2, sort.by="none")` →
  preserves input row order → limma path is correct.
- `deqms_de.R` DEqMS path uses `outputResult(fit, coef=2)`, whose default is to
  **sort rows by significance**; nothing reorders them back.
- `copyDEResultsToFrame` (differential-expression.ts:139) copies the R result
  columns **positionally by row index** with no join on protein ID, so row *i*
  of the protein frame receives the *i*-th most-significant protein's stats.

## Solution

Minimal R-side fix in `deqms_de.R` — reorder back to input order before building
`result` (rownames were set to `seq_len(nrow(exprMat))`):

```r
tt <- outputResult(fit, coef = 2)
tt <- tt[order(as.numeric(rownames(tt))), ]
```

Harden the TS side as defence-in-depth: `copyDEResultsToFrame` should not assume
positional parity — carry an explicit key/order column from the R result and
align on it (also audit `runLimmaDE` / VSN paths for the same latent assumption,
even though limma currently happens to preserve order via `sort.by="none"`).

Add a regression test that pins protein↔stat alignment for DEqMS (e.g. a tiny
fixture where the most-significant protein is NOT row 0, asserting its log2FC
lands on the right protein). Until fixed, **DEqMS must not be used** for any
deliverable; limma is the validated path. Consider disabling/guarding the DEqMS
option in the DE dialog until the fix and test land.

## Resolution (2026-05-16)

Fixed via the key-alignment approach (chosen over a deqms-only reorder so the
same silent-corruption class can't recur for limma/future R DE methods):

- `scripts/deqms_de.R` + `scripts/limma_de.R` now emit a `row` column =
  1-based input index (from `rownames`), in every result branch. limma also
  sets `rownames(exprMat) <- seq_len(nRows)` so the key is well-defined.
- `copyDEResultsToFrame` (now exported) aligns result→df by the `row` key, not
  by position; falls back to positional only if `row` is absent. Output buffers
  pre-filled with FLOAT_NULL so uncovered rows are null, not spurious 0.
- Regression test added: `src/tests/analysis.ts` category "DE result
  alignment" — a shuffled mock result (most-significant NOT row 0) must realign
  by key; plus a positional-fallback test. Pure TS, no R env needed.

**VERIFIED END-TO-END (2026-05-16).** `tsc` clean; regression tests "DE result
alignment" 2/2 PASS on localhost (full build — `--skip-build` reuses a stale
test bundle, must rebuild). Fixed package published to localhost; a real DEqMS
re-run on it, reconciled against candidates, confirms the scramble is gone:
log2FC corr limma~DEqMS 0.992 (was −0.001), cand~DEqMS 0.951 (was −0.005);
DEqMS vs candidates at client thresholds κ=0.899 / 428-of-465 (was κ=0.017 /
35); DEqMS top hits now DMD/ADH1B/SQOR/… (biologically correct). DEqMS is now a
valid, usable DE method. Resolved.
