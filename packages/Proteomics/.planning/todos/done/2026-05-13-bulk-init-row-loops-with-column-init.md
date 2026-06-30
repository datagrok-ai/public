---
created: 2026-05-13T22:51:00.000Z
title: Replace row-by-row col.set() loops with Column.init across the package
area: performance
files:
  - src/analysis/imputation.ts
  - src/analysis/normalization.ts
  - src/analysis/differential-expression.ts
  - src/analysis/enrichment.ts
  - src/parsers/shared-utils.ts
  - src/parsers/maxquant-parser.ts
  - src/parsers/fragpipe-parser.ts
  - src/parsers/spectronaut-candidates-parser.ts
  - src/viewers/heatmap.ts
  - src/viewers/qc-computations.ts
  - src/viewers/enrichment-viewers.ts
---

## Problem

The pattern `df.columns.addNew*(name); for (let i=0; i<df.rowCount; i++) col.set(i, ...)`
appears at ~57 call sites across 10 files in the package. Each `col.set(i, ...)`
hops the JS → native bridge once. On small dev fixtures (8–100 rows) the overhead
is imperceptible; on real proteomics data (a typical Spectronaut Candidates
report is 5k–10k rows; a MaxQuant proteinGroups file can be similar) it becomes
the dominant cost.

Demonstrated on `createVolcanoPlot` (commit 2c54c8eaf4): two row-loops there
accounted for 99% of volcano open time on a real 8329-row Candidates report
(11.7s total). Swapping to `Column.init(callback)` reclaimed 255× — total dropped
from 11750ms to 46ms.

The remaining call sites have the same shape and will have the same problem
the moment they hit a similarly-sized DataFrame. Several of them (the heatmap
top-N filter, the enrichment-volcano wiring, the per-step imputation /
normalization passes) run during the protein-table pipeline on the same
8k-row inputs that triggered the volcano regression.

## Solution

Sweep all 10 files and convert eligible row loops to `Column.init(callback)`.
The shape to look for:

```typescript
const col = df.columns.addNewFloat('foo'); // or addNewString / addNewBool
for (let i = 0; i < df.rowCount; i++)
  col.set(i, /* expr in terms of i */);
```

becomes:

```typescript
const col = df.columns.addNewFloat('foo');
col.init((i) => /* expr in terms of i */);
```

Some loops can't be converted naively — they mutate state across iterations
(e.g. `runDifferentialExpression` collects p-values and FDR-corrects across
all rows before writing back; `imputeKnn` uses scanned neighbor sets). Those
are fine to leave, but the in-place writeback at the end of those functions
can still be a `Column.init` over a precomputed array.

`Column.init` works on already-existing columns too — it doesn't require
addNew. So in-place updates like the `direction` column in volcano.ts (which
needs to preserve the column reference for downstream viewers) can also use it.

Suggested order, by likely impact:

1. **`src/parsers/shared-utils.ts`** — `log2TransformColumns` and
   `copyAsLog2Columns` run on every parser path (MaxQuant, Spectronaut,
   FragPipe, generic) over the full protein × sample matrix. Highest leverage.
2. **`src/analysis/imputation.ts`** — `imputeZero`, `imputeMean`,
   `imputeMedian`, `imputeMinProb` are all pure per-cell rewrites; one-line
   swaps. `imputeKnn` only the final writeback.
3. **`src/analysis/normalization.ts`** — `medianNormalize`, `quantileNormalize`
   final writeback.
4. **`src/viewers/heatmap.ts`** and **`src/viewers/qc-computations.ts`** —
   builds z-score / CV / MA columns over the protein DataFrame.
5. **`src/parsers/maxquant-parser.ts`**, **`fragpipe-parser.ts`**,
   **`spectronaut-candidates-parser.ts`** — `filterContaminantRows` style
   helpers and `significant` column builds.
6. **`src/analysis/enrichment.ts`** and **`src/viewers/enrichment-viewers.ts`** —
   smaller DataFrames (one row per GO term), so impact is bounded but the
   pattern is the same.
7. **`src/analysis/differential-expression.ts`** — only the final p-value /
   significance writeback; the inter-row stats collection has to stay JS-side.

## Measurement

Before and after each batch, time the affected entry point (parse, normalize,
impute, DE, volcano, heatmap) on the real 8329-row Spectronaut Candidates
fixture or the 1569-row CPTAC MaxQuant fixture in `files/demo/`. A
package-wide `Column.init`-only sweep should make every step of the pipeline
sub-100ms on those fixtures.

## Verification

- Full test suite remains 120/120 green.
- No behavior change is intended; this is purely a runtime swap. Visual
  smoke check on the HYE candidates fixture (volcano, heatmap, QC dashboard)
  should be visually identical before and after.

## Origin

Surfaced by user-reported slow volcano on a real 8329-row Spectronaut
Candidates file (2026-05-13). Root-caused, fixed for volcano in commit
2c54c8eaf4, then audited the rest of the package and found the same
pattern everywhere.
