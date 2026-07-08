# Proteomics — Numerical & Mathematical Correctness Audit (Merged)

## 1. Summary

This is the reconciled result of two independent audits of the numerical methods in the
Proteomics plugin. Both audits agree the **core math is correct** — PCA (Jacobi
eigendecomposition + small-covariance trick), differential expression (Welch t-test,
BH-FDR), CV (Welford), MA-plot, MinProb parameters, and the R limma/DEqMS/VSN mechanics
are all sound.

The defects live at the **edges and in the SPC / enrichment glue code**, not in the
central pipeline:

- **1 confirmed code bug both audits found:** Nelson Rule 4 uses the wrong pattern
  (alternation about the centerline instead of alternation of direction).
- **3 real issues only the external audit caught:** the enrichment member-gene alignment
  (highest-impact of the three), the Box–Muller `u1 = 0` edge, and Nelson Rule 8 missing
  the "both sides" condition.
- **1 issue only the internal audit caught:** quantile normalization emits `NaN` for a
  column with a single non-null value.
- **Deviations-from-reference (not hard bugs):** unweighted kNN vs. Troyanskaya's weighted
  average; the "LOESS" trend is a moving average; DEqMS multi-TMT minimum-count rule.

Nothing in the central statistics is wrong. The priority order for fixes is Rule 4 →
enrichment alignment → Rule 8 → Box–Muller → quantile guard.

## 2. Source documents

- [`proteomics-audit-in.md`](proteomics-audit-in.md) — internal audit (line-level read of
  every kernel; caught the quantile `NaN` edge).
- [`proteomics-audit-out.md`](proteomics-audit-out.md) — external audit (verified against
  primary literature and reference implementations; caught the enrichment, Rule 8, and
  Box–Muller issues, with richer deviation-from-reference notes).

## 3. Per-file check

Legend: ✅ correct · ⚠️ deviation / caveat · ❌ bug · ℹ️ informational

| File | Status | Notes |
|------|:------:|-------|
| `src/analysis/spc.ts` | ❌ | Rule 4 wrong pattern (bug); Rule 8 missing "both sides" (deviation). Rules 1–3, 5–7, iterative baseline, meanSd, median, Pearson all ✅ |
| `src/analysis/enrichment.ts` | ❌ | Member-gene recovery indexes `intersections` against the raw query instead of the mapped ENSG list. g:GOSt params, FDR-from-p_value, direction split all ✅ |
| `src/analysis/imputation.ts` | ⚠️ | Box–Muller `u1 = 0` → `Infinity`/`NaN` (bug, rare); kNN unweighted vs. Troyanskaya (deviation). MinProb params ✅ |
| `src/analysis/normalization.ts` | ⚠️ | Quantile emits `NaN` for a single-non-null column (edge bug); interpolation is a valid limma-style variant ℹ️. Median centering ✅, VSN ✅ |
| `src/analysis/differential-expression.ts` | ✅ | log2FC, Welch t-test, BH-FDR (alpha-independent), `row`-key realignment all correct. DEqMS multi-TMT min-count is a doc note ℹ️ |
| `src/analysis/pca.ts` | ✅ | Small-covariance trick, √λ score scaling, Jacobi rotation, centering, variance-explained all verified correct |
| `src/viewers/qc-computations.ts` | ✅ | MA-plot and Welford CV (with `2^log2` back-transform) correct; "LOESS" trend is a moving average ℹ️ |

## 4. Issues (severity-descending)

### 4.1 Nelson Rule 4 uses the wrong pattern — **High**

- **Where:** [`src/analysis/spc.ts:291-301`](src/analysis/spc.ts#L291-L301)
- **Summary:** Canonical Nelson Test 4 (1984) is "14 points in a row alternating in
  **direction**" — a zigzag detected from the sign of successive differences
  `series[i] − series[i-1]`; the centerline is irrelevant. The code instead classifies each
  point as above/below the mean and requires that label to alternate. A true sawtooth that
  stays entirely above the mean is missed; a slow oscillation crossing the mean every point
  false-fires. Disabled by default (`nelson_4: false`), so no live impact until enabled, but
  it is a defect in a named canonical algorithm. Both audits agree.
- **How to fix:** compute `d[i] = sign(series[i] − series[i-1])` and require `d[i] = −d[i-1]`
  for 13 consecutive steps across the last 14 points; handle zero-differences explicitly.

### 4.2 Enrichment member-gene names can be misaligned — **Medium-High**

- **Where:** [`src/analysis/enrichment.ts:213-222`](src/analysis/enrichment.ts#L213-L222)
- **Summary:** `buildEnrichmentDf` recovers the member genes per term by indexing
  `r.intersections[j]` against `queryGenes[j]` (raw input order). Per the g:Profiler API,
  the `intersections` array is positionally aligned to the **mapped** Ensembl gene list
  (`genes_metadata.query.<name>.ensgs`), not the raw query. This is only correct when every
  input maps 1:1 to a single ENSG with no unmapped, multi-mapped, or background-filtered
  genes. Otherwise the "Intersection" column shows the wrong gene names. Missed by the
  internal audit (enrichment was wrongly dismissed as "fully delegated"); confirmed real.
- **How to fix:** align `intersections[j]` to `genes_metadata.query.<name>.ensgs[j]` and
  reverse-map ENSG → original symbol via the response's `mapping` table (as the gprofiler2
  R client does), instead of indexing against `queryGenes`.

### 4.3 Nelson Rule 8 missing the "both sides" condition — **Medium**

- **Where:** [`src/analysis/spc.ts:334-338`](src/analysis/spc.ts#L334-L338)
- **Summary:** Nelson Test 8 is "eight points in a row on **both sides** of the centerline
  with none in Zone C" (a mixture pattern). The code only checks `|z| > 1` for all 8 points
  and does not require at least one above and one below the mean, so it also fires on eight
  consecutive points beyond 1σ on a single side — which is a different (shift) signal. This
  matches the Minitab "either side" convention but not Nelson's literal definition. The
  internal audit marked this correct; the external audit is right.
- **How to fix:** additionally require the 8-point window to contain ≥1 point above and ≥1
  below the mean. If the Minitab convention is intended, document it explicitly instead.

### 4.4 Box–Muller can produce Infinity/NaN — **Medium**

- **Where:** [`src/analysis/imputation.ts:9-14`](src/analysis/imputation.ts#L9-L14)
- **Summary:** `z = √(−2·ln u1)·cos(2π·u2)` is the standard transform, but `Math.random()`
  returns `[0, 1)`, so `u1` can be exactly 0 → `ln(0) = −∞` → an `Infinity`/`NaN` draw. Over
  many imputed cells this eventually occurs, corrupting MinProb output. The internal audit
  noted only MinProb non-determinism; the external audit caught the concrete edge.
- **How to fix:** draw `u1` in `(0, 1]` via `u1 = 1 − Math.random()`, or resample while
  `u1 === 0`.

### 4.5 kNN uses an unweighted neighbor average — **Medium (deviation) / Low**

- **Where:** [`src/analysis/imputation.ts:163-180`](src/analysis/imputation.ts#L163-L180)
- **Summary:** Troyanskaya (2001) KNNimpute imputes each missing value with an
  inverse-distance-**weighted** average of the k neighbors; the code uses an **unweighted**
  mean, biasing estimates slightly toward more-distant neighbors. Orientation (impute per
  protein/row using nearest rows) and the overlap-normalized distance are both fine. Note:
  unweighted kNN imputation is itself a legitimate method and the code does not claim
  Troyanskaya specifically, so the real severity is closer to informational.
- **How to fix:** weight each neighbor by `1/distance` (or a Gaussian kernel) to match the
  cited reference; or document that an unweighted variant is intended.

### 4.6 Quantile normalization emits NaN for a single-non-null column — **Low**

- **Where:** [`src/analysis/normalization.ts:90-97`](src/analysis/normalization.ts#L90-L97)
  (and the `maxValid === 1` case at [`normalization.ts:74`](src/analysis/normalization.ts#L74))
- **Summary:** When a column has exactly one non-null value (`n === 1`) while `maxValid > 1`,
  `globalPos = (r / (n - 1)) * (maxValid - 1)` evaluates `0 / 0 = NaN`, and `NaN` is written
  back as the normalized value. Degenerate for real proteomics (thousands of proteins), but
  it produces a wrong value instead of a no-op. Caught only by the internal audit.
- **How to fix:** guard with `if (n < 2) continue;` in the assignment loop and skip the
  `maxValid < 2` case up front.

### 4.7 DEqMS multi-TMT minimum-count rule not applied — **Low**

- **Where:** [`src/analysis/differential-expression.ts:236-267`](src/analysis/differential-expression.ts#L236-L267)
- **Summary:** For inputs spanning multiple TMT experiments, Zhu et al. (2020) recommend
  using the **minimum** PSM/peptide count across experiments as the variance-regression
  covariate (data set D5). The code passes the user-selected count column straight through,
  which is fine for a single set but suboptimal for multi-set inputs. Not a bug for the
  common single-set case.
- **How to fix:** for multi-set inputs, apply (or at least document) the
  minimum-count-across-experiments rule before calling `DeqmsDE`.

### 4.8 "LOESS" trend is a moving average — **Informational**

- **Where:** [`src/viewers/qc-computations.ts:57-106`](src/viewers/qc-computations.ts#L57-L106)
- **Summary:** `computeLoessTrend` is a window-mean over points sorted by A — no local
  weighted regression, no tricube weighting. The docstring's "approximates loess" understates
  the difference. Both audits agree; internally sound as a visual trend line.
- **How to fix:** rename to `computeMovingAverageTrend` to avoid implying local regression;
  implement local linear regression with tricube weights if a true loess line is wanted.

### 4.9 Quantile NA-handling reference and MNAR caveat — **Informational**

- **Where:** [`src/analysis/normalization.ts:37-103`](src/analysis/normalization.ts#L37-L103)
- **Summary:** The fractional-rank interpolation for unequal non-missing counts is a valid
  variant matching limma `normalizeQuantiles` (not identical to preprocessCore). Results will
  differ slightly between the two references, and with heavy/non-random (MNAR) missingness —
  common in proteomics — neither variant is guaranteed unbiased.
- **How to fix:** state in the docstring that the implementation follows the limma
  interpolation convention, and note the MNAR caveat.

### 4.10 Inconsistent null tests — **Low**

- **Where:** `src/analysis/normalization.ts`, `src/viewers/qc-computations.ts`
- **Summary:** Some paths test `raw[i] === DG.FLOAT_NULL` (sentinel on a materialized array)
  while others use `col.isNone(i)`. If the two ever disagree (or `col.length` differs from
  `raw.length`), a cell could be treated as present/absent inconsistently. Not a math error,
  but it can corrupt which cells are normalized/counted.
- **How to fix:** standardize on a single null test throughout.
