# Statistics

Statistical methods for two-group comparison, rank correlation, multi-group
contrasts, dose-response analysis, exact categorical tests, and covariate
adjustment. Reference values are validated against scipy, R, SAS PROC GLM,
and the original publications. The mathematical backend is [jstat](https://github.com/jstat/jstat) —
distributions, special functions, and matrix inverse — wrapped in a typed
TypeScript API.

```typescript
import {stats} from '@datagrok-libraries/sci-comp';

const {
  welchTTest, mannWhitneyU, hedgesG,
  spearman, severityTrend,
  welchPairwise, dunnettPairwise, bonferroniCorrect,
  jonckheere, cochranArmitage, thresholdTest,
  fisherExact2x2,
  pavaIncreasing, pavaDecreasing, williamsTest,
  runAncova,
} = stats;
```

| Method | Use case | Reference |
|--------|----------|-----------|
| `welchTTest` | Two-sample mean comparison, unequal variances | [Welch (1947)](https://doi.org/10.1093/biomet/34.1-2.28) |
| `mannWhitneyU` | Distribution-free two-sample test | [Mann & Whitney (1947)](https://doi.org/10.1214/aoms/1177730491) |
| `hedgesG` | Bias-corrected effect size (Cohen's d × J) | [Hedges (1981)](https://doi.org/10.3102/10769986006002107) |
| `spearman` | Rank correlation | [Spearman (1904)](https://doi.org/10.2307/1412159) |
| `severityTrend` | Dose-severity rank correlation (with constant-y guard) | — |
| `welchPairwise` | Per-treated-group raw p-values vs control | — |
| `dunnettPairwise` | Many-to-one comparison vs control (FWER-controlled) | [Dunnett (1955)](https://doi.org/10.1080/01621459.1955.10501294) |
| `bonferroniCorrect` | Multiplicity correction | [Bonferroni (1936)](https://en.wikipedia.org/wiki/Bonferroni_correction) |
| `jonckheere` | Trend test for ordered groups (continuous) | [Jonckheere (1954)](https://doi.org/10.1093/biomet/41.1-2.133) |
| `cochranArmitage` | Trend test for proportions (incidence) | [Cochran (1954)](https://doi.org/10.2307/3001616), [Armitage (1955)](https://doi.org/10.2307/3001775) |
| `thresholdTest` | Sequential threshold test for proportions | [Young (1985)](https://www.sra.org/) |
| `fisherExact2x2` | Exact test for 2×2 contingency tables | [Fisher (1935)](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) |
| `pavaIncreasing/Decreasing` | Isotonic regression (PAVA) | [Barlow et al. (1972)](https://en.wikipedia.org/wiki/Isotonic_regression) |
| `williamsTest` | Step-down dose-response with PAVA | [Williams (1971)](https://doi.org/10.2307/2528930), [Williams (1972)](https://doi.org/10.2307/2556895) |
| `runAncova` | Covariate-adjusted group comparison (OLS) | Standard ANCOVA |

All methods accept a flexible numeric input type:

```typescript
type NumericInput =
  | readonly number[]
  | Int8Array | Uint8Array | Int16Array | Uint16Array
  | Int32Array | Uint32Array | Float32Array | Float64Array;
```

Inputs are normalised to `Float64Array` internally. NaN is the standard
"missing" sentinel — methods that take samples strip NaN before computing
(or do pairwise NaN removal for paired inputs like Spearman and ANCOVA).

## Comparing two groups

Three views on the same data: parametric significance (Welch), distribution-free
significance (MWU), and the magnitude of the difference (Hedges' g).

```typescript
const a = [69, 70, 66, 63, 68, 70, 69, 67, 62];
const b = [68, 62, 67, 68, 69, 67, 61, 59];

const t = welchTTest(a, b);
// {statistic: ~1.31, pValue: ~0.20}

const u = mannWhitneyU(a, b, 'two-sided');
// {statistic: ~52, pValue: ~0.18}

const g = hedgesG(a, b);
// ~0.43  (small effect)
```

Mann-Whitney uses the **exact** distribution when `min(n1, n2) ≤ 8` without
ties (matching scipy's `method='auto'`), and the asymptotic normal
approximation with continuity correction and tie adjustment otherwise.

## Correlation

```typescript
const x = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110];
const y = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17];

const r = spearman(x, y);
// {rho: -29/165 ≈ -0.176, pValue: ~0.63}
```

`severityTrend` is a thin wrapper that explicitly returns `null` when the
severities are constant (correlation is undefined for a flat response):

```typescript
severityTrend([0, 10, 50, 100], [0.1, 0.5, 1.2, 2.0]);
// {rho: 1.0, pValue: 0}

severityTrend([0, 10, 50, 100], [1.0, 1.0, 1.0, 1.0]);
// {rho: null, pValue: null}
```

## Pairwise vs control (FWER control)

Compare each treated group against the control while keeping the
family-wise error rate at α. Two strategies on the same data:

```typescript
const control = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32];
const treated = [
  {doseLevel: 1, values: [9.76, 8.80, 7.68, 9.36]},
  {doseLevel: 2, values: [12.80, 9.68, 12.16, 9.20, 10.55]},
];

// Strategy 1: Welch + Bonferroni
const raw = welchPairwise(control, treated).map((r) => r.pValueWelch);
const adj = bonferroniCorrect(raw);
// adj ≈ [0.621, 0.030]

// Strategy 2: Dunnett's many-to-one
const dun = dunnettPairwise(control, treated);
// [{doseLevel: 1, statistic: ~0.857, pValueAdj: ~0.620, ...},
//  {doseLevel: 2, statistic: ~3.694, pValueAdj: ~0.006, ...}]
```

Dunnett uses the joint distribution of the contrasts via 2D numerical
integration of the Dunnett (1955) factorisation — uniformly more powerful
than Welch + Bonferroni when within-group variance is homogeneous.

## Trend tests

Both reduce to a Z-statistic and a two-sided p-value. Continuous data uses
**Jonckheere-Terpstra** (sum of pairwise Mann-Whitney counts); incidence
data uses **Cochran-Armitage** (linear contrast over scores):

```typescript
// Continuous trend across ordered groups
jonckheere([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
// {statistic: ~3.32, pValue: ~0.0009}

// Incidence trend (Tang 2006: 0/12, 0/12, 1/12, 3/12)
cochranArmitage([0, 0, 1, 3], [12, 12, 12, 12]);
// {zStatistic: ~2.34, chi2Statistic: ~5.45, pValue: ~0.020, ...}
```

`cochranArmitage` accepts custom scores, alternative direction, variance
method (binomial / hypergeometric), and an optional Buonaccorsi-style
modified statistic:

```typescript
cochranArmitage([0, 0, 1, 3], [12, 12, 12, 12], {
  scores: [0, 1, 4, 8],          // dose-proportional scores
  alternative: 'increasing',
  variance: 'hypergeometric',    // matches R prop.trend.test
  modified: true,                // adds zModified, pValueModified
});
```

## Threshold test (sequential dose-response)

Walks from the lowest dose, pooling each non-significant group into a
cumulative "control". Stops at the first significant comparison —
that group is the Effect Level (EL); prior groups are NOELs.

```typescript
// Young (1985) Table 5(c): NCI mice DDT data
const counts = [4, 1, 1, 7];
const totals = [56, 49, 48, 41];

const steps = thresholdTest(counts, totals, {
  alpha: 0.05,
  adjustAlpha: true,             // Šidák correction across the k-1 comparisons
});
// 3 steps; last step has effectGroup = 3, noelGroups = [0, 1, 2]
```

## Fisher exact 2×2

Two-sided p-value uses the "minlike" definition (matches scipy's default).
Returns one-sided p-values and the sample odds ratio as well.

```typescript
// Lady Tasting Tea (Fisher 1935)
const r = fisherExact2x2([[3, 1], [1, 3]]);
// {oddsRatio: 9, pValue: 34/70, pGreater: 17/70, pLess: 69/70}

// Perfect separation → Infinity odds ratio, valid p-value
fisherExact2x2([[5, 0], [0, 5]]);
// {oddsRatio: Infinity, pValue: 1/126, ...}
```

## Williams step-down (dose-response trend)

Two ingredients:

```typescript
// 1. Pool-Adjacent-Violators isotonic regression
pavaIncreasing([10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7], [1, 1, 1, 1, 1, 1, 1]);
// Float64Array [10.1, 10.1, 10.1, 10.6, 11.4, 11.8, 11.8]
```

PAVA is also exposed standalone (`pavaIncreasing`, `pavaDecreasing`) for
isotonic-regression use cases outside Williams' test.

```typescript
// 2. Step-down test using critical-value tables from Williams 1971/1972
const w = williamsTest(
  [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7],   // means (index 0 = control)
  Array(7).fill(Math.sqrt(1.16)),               // standard deviations
  [8, 8, 8, 8, 8, 8, 8],                        // sample sizes
  ['ctrl', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'], // labels
  {direction: 'increase', alpha: 0.05},
);
// w.minimumEffectiveDose === 'd4'
// w.stepDownResults: [d6 sig, d5 sig, d4 sig, d3 NS] — stops at first non-sig
```

`direction: 'auto'` (default) infers the sense of the trend from the
highest-dose mean vs the control mean. Critical values come from
[`williams-tables.ts`](./tests/williams-tables.ts) — the digitised
Williams (1971, 1972) tables with the 1972 extrapolation formula for
unequal control-vs-dose replication (`w = c/r`).

## ANCOVA

Fits `y ~ C(group) + covariate` via OLS and reports adjusted (LS) means,
pairwise contrasts vs control, slope homogeneity F-test, and an effect
decomposition (total = direct + indirect):

```typescript
// Montgomery (2012) Table 15.10 — fiber strength by machine, diameter as covariate
const y = [36, 41, 39, 42, 49, 40, 48, 39, 45, 44, 35, 37, 42, 34, 32];
const x = [20, 25, 24, 25, 32, 22, 28, 22, 30, 28, 21, 23, 26, 21, 15];
const g = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3];

const r = runAncova(y, x, g, {controlGroup: 1, alpha: 0.05});
// r.modelRSquared        ≈ 0.9192   (matches SAS: 0.9192)
// r.adjustedMeans[0]     ≈ {group: 1, adjustedMean: 40.382, ...}
// r.slope                ≈ {tStatistic: 8.365, pValue: 4e-6}
// r.slopeHomogeneity     ≈ {fStatistic: 0.488, pValue: 0.629, homogeneous: true}
```

The slope-homogeneity F-test is a precondition: when slopes differ
significantly across groups (`homogeneous: false`), the basic ANCOVA model
is misspecified and the adjusted means should not be interpreted at face
value.

For Lazic-style organ-free covariates (Lazic et al. 2020), set
`useOrganFreeBw: true` — the covariate becomes `BW − organ_value`.

Returns `null` when there are too few observations (`n < k + 2` or `k < 2`).

## Running examples

Runnable examples are located in `examples/`:

```bash
# Welch t-test + Mann-Whitney U + Hedges' g on the same data
npx tsx src/stats/examples/compare-two-groups.ts

# Spearman rank correlation + severity-trend (with constant-y edge case)
npx tsx src/stats/examples/correlation.ts

# Welch + Bonferroni vs Dunnett's many-to-one (FWER comparison)
npx tsx src/stats/examples/pairwise-vs-control.ts

# Jonckheere-Terpstra + Cochran-Armitage trend tests
npx tsx src/stats/examples/trend-tests.ts

# Sequential threshold test for proportions (Young 1985)
npx tsx src/stats/examples/threshold-test.ts

# Fisher 2×2 exact test (Lady Tasting Tea + edge cases)
npx tsx src/stats/examples/fisher-exact.ts

# PAVA + Williams' step-down test (Williams 1971/1972 numerical examples)
npx tsx src/stats/examples/williams-step-down.ts

# ANCOVA on Montgomery Table 15.10 (reproduces SAS PROC GLM)
npx tsx src/stats/examples/ancova.ts
```

## Validation

Each method has a Jest test suite in `__tests__/` that loads JSON fixtures
from `__tests__/fixtures/`. The fixtures are generated by
[`scripts/generate-fixtures.py`](../../scripts/generate-fixtures.py), which
runs the SENDEX Python validation suite with scipy and serialises the
input data + reference values. A total of 179 cases across 14 fixtures —
all pass with the same tolerances as the Python suite (1e-3 for p-values,
1e-2 for statistics, tighter for closed-form quantities).

To regenerate fixtures (only when scipy versions change or a new test case
is added), use the SEND-TEST conda environment:

```bash
"$USERPROFILE/anaconda3/envs/SEND-TEST/python.exe" scripts/generate-fixtures.py
```
