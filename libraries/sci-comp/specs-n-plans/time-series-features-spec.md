# Time Series Feature Extraction — Implementation Spec

## Overview

TypeScript implementation of 34 Tier 1 tsfresh-compatible feature calculators for short time series (N=10–15).
Part of `@datagrok-libraries/sci-comp`, located in `src/time-series/feature-extraction/`.

## Dependencies

- **jstat** `^1.9.6` — for `jStat.studentt.cdf()` in `linear_trend` pvalue calculation.

No other external dependencies.

## File Structure

```
libraries/sci-comp/src/
  time-series/
    feature-extraction/
      types.ts
      dataframe.ts
      stats-utils.ts
      linear-trend.ts
      calculators.ts
      extract.ts
      index.ts
```

---

## Module: `types.ts`

```typescript
export type NumericArray =
  | Int32Array<ArrayBufferLike>
  | Uint32Array<ArrayBufferLike>
  | Float32Array<ArrayBufferLike>
  | Float64Array<ArrayBufferLike>;

export interface TimeSeriesColumn {
  name: string;
  data: NumericArray;
}

export interface TimeSeriesDataFrame {
  ids: Uint32Array<ArrayBufferLike>;
  time: NumericArray;
  rowCount: number;
  columns: TimeSeriesColumn[];
}

export interface FeatureColumn {
  name: string;
  data: Float64Array;
}

export interface FeatureMatrix {
  sampleIds: Uint32Array<ArrayBufferLike>;
  columns: FeatureColumn[];
  nSamples: number;
}

export interface ExtractOptions {
  validate?: boolean;  // default: false
}
```

---

## Module: `dataframe.ts`

### `SampleRange`

```typescript
interface SampleRange {
  start: number;  // index in flat arrays
  length: number; // number of points in this sample
}
```

### `buildIndex(df: TimeSeriesDataFrame): Map<number, SampleRange>`

Single pass over `df.ids` (O(n)). Assumes data is sorted by id (all rows for one id are contiguous).
Returns map from id → {start, length}.

Algorithm:
```
prevId = ids[0], start = 0
for i = 1..rowCount-1:
  if ids[i] !== prevId:
    map.set(prevId, {start, length: i - start})
    prevId = ids[i]
    start = i
map.set(prevId, {start, length: rowCount - start})
```

### `validate(df: TimeSeriesDataFrame): void`

Checks:
1. `ids.length === time.length === rowCount`
2. Every `column.data.length === rowCount`
3. For Float32Array/Float64Array columns: no NaN or ±Inf values
4. ids are contiguous (sorted by id)

Throws `Error` with descriptive message on failure.

---

## Module: `stats-utils.ts`

### `tdistCdf(t: number, df: number): number`

Wrapper around jStat:

```typescript
import { jStat } from 'jstat';

export function tdistCdf(t: number, df: number): number {
  return jStat.studentt.cdf(t, df);
}
```

### `twoTailPvalue(tStat: number, df: number): number`

```typescript
export function twoTailPvalue(tStat: number, df: number): number {
  return 2 * (1 - tdistCdf(Math.abs(tStat), df));
}
```

---

## Module: `linear-trend.ts`

### `linearTrend(x: NumericArray, n: number): LinearTrendResult`

```typescript
export interface LinearTrendResult {
  slope: number;
  intercept: number;
  rvalue: number;
  pvalue: number;
  stderr: number;
}
```

Single-pass OLS regression against t = 0, 1, 2, ..., n-1.

**Formulas** (one pass, accumulate `sumX`, `sumX2`, `sumT`, `sumT2`, `sumTX`):

Since t = 0..n-1:
- `sumT = n*(n-1)/2`
- `sumT2 = n*(n-1)*(2*n-1)/6`

These are computed analytically, no loop needed. Only `sumX`, `sumX2`, `sumTX` require a loop.

```
meanT = sumT / n
meanX = sumX / n

SS_tt = sumT2 - n * meanT^2        // Σ(t - meanT)^2
SS_xx = sumX2 - n * meanX^2        // Σ(x - meanX)^2
SS_tx = sumTX - n * meanT * meanX  // Σ(t - meanT)(x - meanX)

slope = SS_tx / SS_tt              // SS_tt > 0 when n >= 2
intercept = meanX - slope * meanT
rvalue = SS_tx / sqrt(SS_tt * SS_xx)  // 0 if SS_xx == 0

RSS = SS_xx - slope * SS_tx        // residual sum of squares
if n > 2 and SS_tt > 0:
  stderr = sqrt(RSS / (n - 2)) / sqrt(SS_tt)
else:
  stderr = 0

if stderr > 0:
  tStat = slope / stderr
  pvalue = twoTailPvalue(tStat, n - 2)  // uses jStat
else:
  pvalue = 1.0
```

Edge cases:
- `n < 2`: return all zeros (should not happen per spec, but guard)
- `SS_xx == 0` (constant series): rvalue = 0, slope = 0, pvalue = 1.0
- `n == 2`: stderr = 0 (no residual df), pvalue = 1.0 if slope is finite

Match: `scipy.stats.linregress` behavior.

---

## Module: `calculators.ts`

### Design: Three Passes + Median + Linear Trend

For each sample slice, computations are grouped into passes to minimize iterations.

### Pass 1: Basic statistics

**Input:** `data: NumericArray`, `n: number`
**Output:** `Pass1Result`
**Iterations:** one loop over `data[0..n-1]`

```typescript
export interface Pass1Result {
  n: number;
  sum: number;
  sumSq: number;
  min: number;
  max: number;
  firstMinIdx: number;
  firstMaxIdx: number;
  lastMinIdx: number;
  lastMaxIdx: number;
  minCount: number;
  maxCount: number;
  mean: number;     // derived: sum / n
  variance: number; // derived: sumSq/n - mean^2  (population, ddof=0)
  std: number;      // derived: sqrt(variance)
  range: number;    // derived: max - min
}
```

Algorithm:
```
sum = 0, sumSq = 0
min = +Inf, max = -Inf
firstMinIdx = 0, firstMaxIdx = 0
lastMinIdx = 0, lastMaxIdx = 0
minCount = 0, maxCount = 0

for i = 0..n-1:
  v = data[i]
  sum += v
  sumSq += v * v
  if v < min:
    min = v; firstMinIdx = i; lastMinIdx = i; minCount = 1
  else if v == min:
    lastMinIdx = i; minCount++
  if v > max:
    max = v; firstMaxIdx = i; lastMaxIdx = i; maxCount = 1
  else if v == max:
    lastMaxIdx = i; maxCount++

mean = sum / n
variance = sumSq / n - mean * mean
// guard against floating-point negative variance:
if variance < 0: variance = 0
std = sqrt(variance)
range = max - min
```

**Features extracted from Pass 1:**

| Feature | Expression |
|---------|------------|
| `mean` | `p1.mean` |
| `minimum` | `p1.min` |
| `maximum` | `p1.max` |
| `length` | `p1.n` |
| `sum_values` | `p1.sum` |
| `abs_energy` | `p1.sumSq` |
| `root_mean_square` | `sqrt(p1.sumSq / n)` |
| `standard_deviation` | `p1.std` |
| `variance` | `p1.variance` |
| `variance_larger_than_standard_deviation` | `p1.variance > 1 ? 1.0 : 0.0` |
| `first_location_of_minimum` | `p1.firstMinIdx / n` |
| `first_location_of_maximum` | `p1.firstMaxIdx / n` |
| `last_location_of_minimum` | `p1.lastMinIdx / n` |
| `last_location_of_maximum` | `p1.lastMaxIdx / n` |
| `has_duplicate_max` | `p1.maxCount >= 2 ? 1.0 : 0.0` |
| `has_duplicate_min` | `p1.minCount >= 2 ? 1.0 : 0.0` |

### Pass 2: Differences, uniqueness, crossings

**Input:** `data: NumericArray`, `n: number`, `mean: number`, `std: number`
**Output:** `Pass2Result`
**Iterations:** one loop over `data[0..n-1]`, tracking previous value

```typescript
export interface Pass2Result {
  sumAbsDiff: number;
  sumSqDiff: number;
  uniqueCount: number;
  hasDuplicate: boolean;
  reoccurringDatapointCount: number;
  crossingCount: number;       // crossings of M_VALUE
  longestStrikeAbove: number;
  longestStrikeBelow: number;
}
```

Algorithm:
```
sumAbsDiff = 0, sumSqDiff = 0
crossingCount = 0

// Uniqueness via Map<number, number> (value → count)
valueCounts = new Map()

// Strikes (depend on mean, computed in pass 1)
curAbove = 0, maxAbove = 0
curBelow = 0, maxBelow = 0

prev = data[0]
valueCounts.set(prev, 1)
// init strikes for first element
if prev > mean: curAbove = 1
if prev < mean: curBelow = 1

for i = 1..n-1:
  v = data[i]
  d = v - prev

  // Differences
  sumAbsDiff += abs(d)
  sumSqDiff += d * d

  // Crossings of m=0 (configurable)
  if (prev - M) * (v - M) < 0: crossingCount++

  // Uniqueness
  valueCounts.set(v, (valueCounts.get(v) ?? 0) + 1)

  // Strikes
  if v > mean:
    curAbove++; maxAbove = max(maxAbove, curAbove)
    curBelow = 0
  else if v < mean:
    curBelow++; maxBelow = max(maxBelow, curBelow)
    curAbove = 0
  else:
    curAbove = 0; curBelow = 0

  prev = v

maxAbove = max(maxAbove, curAbove)
maxBelow = max(maxBelow, curBelow)

// Uniqueness summary
uniqueCount = valueCounts.size
hasDuplicate = uniqueCount != n
reoccurringDatapointCount = 0
for count of valueCounts.values():
  if count > 1: reoccurringDatapointCount += count
```

**Features extracted from Pass 2:**

| Feature | Expression |
|---------|------------|
| `mean_change` | `(data[n-1] - data[0]) / (n - 1)` (direct access, not from pass) |
| `mean_abs_change` | `p2.sumAbsDiff / (n - 1)` |
| `absolute_sum_of_changes` | `p2.sumAbsDiff` |
| `mean_second_derivative_central` | `(data[n-1] - data[n-2] - data[1] + data[0]) / (2 * (n - 2))` (direct) |
| `cid_ce (normalize=false)` | `sqrt(p2.sumSqDiff)` |
| `has_duplicate` | `p2.hasDuplicate ? 1.0 : 0.0` |
| `percentage_of_reoccurring_datapoints_to_all_datapoints` | `p2.reoccurringDatapointCount / n` |
| `ratio_value_number_to_time_series_length` | `p2.uniqueCount / n` |
| `number_crossing_m` | `p2.crossingCount` |
| `longest_strike_above_mean` | `p2.longestStrikeAbove` |
| `longest_strike_below_mean` | `p2.longestStrikeBelow` |

### `cid_ce` with normalize=true

Requires separate computation when std > 0:
```
z = (data - mean) / std
cidNorm = sqrt( Σ(z[i+1] - z[i])^2 )
       = sqrt( Σ((data[i+1] - data[i]) / std)^2 )
       = sqrt( sumSqDiff ) / std
       = cid_ce_false / std
```

So: `cidNormalized = std > 0 ? sqrt(p2.sumSqDiff) / p1.std : 0`

No extra loop needed.

### Pass 3: Threshold-based features

**Input:** `data: NumericArray`, `n: number`, `mean: number`, `std: number`, `median: number`, `range: number`
**Output:** `Pass3Result`
**Iterations:** one loop over `data[0..n-1]`

```typescript
export interface Pass3Result {
  countAboveMean: number;
  countBelowMean: number;
  // ratio_beyond_r_sigma for each r:
  beyondRSigmaCounts: number[];  // one per R_SIGMA_VALUES entry
}
```

Algorithm:
```
countAbove = 0, countBelow = 0
beyondCounts = [0, 0, 0]  // for r = 1, 1.5, 2

for i = 0..n-1:
  v = data[i]
  if v > mean: countAbove++
  if v < mean: countBelow++
  absDeviation = abs(v - mean)
  for j = 0..len(R_SIGMA_VALUES)-1:
    if absDeviation > R_SIGMA_VALUES[j] * std:
      beyondCounts[j]++
```

**Features extracted from Pass 3:**

| Feature | Expression |
|---------|------------|
| `count_above_mean` | `p3.countAboveMean` |
| `count_below_mean` | `p3.countBelowMean` |
| `ratio_beyond_r_sigma` (per r) | `p3.beyondRSigmaCounts[j] / n` |

**Features computed from pass 1 + median (no extra loop):**

| Feature | Expression |
|---------|------------|
| `symmetry_looking` (per r) | `abs(mean - median) < r * range ? 1.0 : 0.0` |
| `large_standard_deviation` (per r) | `std > r * range ? 1.0 : 0.0` |

### `computeMedian(data: NumericArray, n: number, buf: Float64Array): number`

Copy to buffer, sort, return middle element.

```
for i = 0..n-1: buf[i] = data[i]
sort buf[0..n-1]
if n % 2 == 0: return (buf[n/2 - 1] + buf[n/2]) / 2
else: return buf[(n-1) / 2]
```

`buf` is pre-allocated once in the orchestrator (size = max sample length).

---

## Module: `extract.ts`

### `extractFeatures(df: TimeSeriesDataFrame, options?: ExtractOptions): FeatureMatrix`

Orchestrator. Steps:

1. **Validate** (if `options?.validate`).
2. **Build index** — `buildIndex(df)` → `Map<number, SampleRange>`.
3. **Compute feature names** — for each column × each feature → name string.
4. **Allocate result** — `FeatureColumn[]`, each with `Float64Array(nSamples)`.
5. **Allocate reusable buffer** — `Float64Array(maxSampleLength)` for median.
6. **Main loop** — for each sample (si = 0..nSamples-1):
   - For each column (col):
     - `slice = col.data.subarray(range.start, range.start + range.length)`
     - `n = range.length`
     - Run pass1, computeMedian, pass2, pass3, linearTrend
     - Write results to `resultColumns[offset + featureIdx].data[si]`

### Feature naming convention

Pattern: `{columnName}__{featureName}` or `{columnName}__{featureName}__{paramName}_{paramValue}`

Examples:
```
pH__mean
pH__cid_ce__normalize_true
pH__cid_ce__normalize_false
pH__symmetry_looking__r_0.05
pH__linear_trend__attr_slope
pH__ratio_beyond_r_sigma__r_2
pH__number_crossing_m__m_0
concentration__mean
```

### Features per column

| Group | Count |
|-------|-------|
| Simple parameterless (mean, median, ..., longest_strike_below_mean) | 32 |
| cid_ce × 2 (normalize true/false) | 2 |
| symmetry_looking × len(R_VALUES) | 3 |
| large_standard_deviation × len(R_VALUES) | 3 |
| linear_trend × 5 (slope, intercept, rvalue, pvalue, stderr) | 5 |
| ratio_beyond_r_sigma × len(R_SIGMA_VALUES) | 3 |
| number_crossing_m × 1 | 1 |
| **Total per column** | **~49** |

For 2 columns (pH + concentration): **~98 features**.

### Configuration

Default parameters are constants:

```typescript
const R_VALUES = [0.05, 0.25, 0.45];
const R_SIGMA_VALUES = [1, 1.5, 2];
const M_VALUES = [0];
```

These match the Python test config. Future: accept via `ExtractOptions`.

---

## Module: `index.ts`

Public API re-exports:

```typescript
export { extractFeatures } from './extract';
export { validate } from './dataframe';
export type {
  NumericArray,
  TimeSeriesColumn,
  TimeSeriesDataFrame,
  FeatureColumn,
  FeatureMatrix,
  ExtractOptions,
} from './types';
```

---

## Edge Cases

| Condition | Behavior |
|-----------|----------|
| `n == 0` | Skip sample, all features = NaN |
| `n == 1` | mean/min/max/sum valid; diff-based features = 0 or NaN; linear_trend = all 0 |
| `n == 2` | All features valid except: `mean_second_derivative_central` = 0 (division by 0 guard); linear_trend stderr = 0, pvalue = 1.0 |
| `std == 0` (constant series) | `cid_ce(normalize=true)` = 0; `ratio_beyond_r_sigma` = 0; `linear_trend` rvalue = 0, slope = 0 |
| All same value | `has_duplicate` = true; `percentage_of_reoccurring_datapoints` = 1.0; `variance` = 0 |
| Negative values | All formulas work; `number_crossing_m` with m=0 counts zero-crossings |

---

## Boolean Feature Encoding

tsfresh returns Python `bool` for some features. In `FeatureColumn.data` (Float64Array), encode as:
- `true` → `1.0`
- `false` → `0.0`

Affected features: `variance_larger_than_standard_deviation`, `has_duplicate`, `has_duplicate_max`, `has_duplicate_min`, `symmetry_looking`, `large_standard_deviation`.

---

## Performance Notes

- All passes are O(n) per sample. Total: O(S × C × n) where S = samples, C = columns.
- Median requires O(n log n) sort but n ≤ 15, so negligible.
- One `Float64Array` allocation for sort buffer, reused across all samples.
- `subarray()` returns a view (no copy) for all TypedArray types.
- Map allocation in Pass 2 for uniqueness: at most n=15 entries, trivial.
- jStat `studentt.cdf` called once per column per sample for linear_trend pvalue.
