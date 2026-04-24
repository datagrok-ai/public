# Time Series Feature Extraction — Test Spec

## Overview

Two test strategies: **unit tests** (per-calculator) and **integration tests** (full pipeline).

## Test Framework

Standard Jest/Vitest (whichever sci-comp already uses). No Python dependency at runtime —
golden values are pre-computed and stored as constants.

## Tolerances

| Category | Tolerance | Rationale |
|----------|-----------|-----------|
| Exact integer features | `=== expected` | length, counts, boolean-as-float |
| Basic statistics | `atol = 1e-10` | Simple sums, direct formulas |
| Derived statistics (std, variance, rms) | `atol = 1e-10` | One sqrt, minimal error accumulation |
| linear_trend slope, intercept, rvalue, stderr | `atol = 1e-7` | Multiple divisions, accumulated rounding |
| linear_trend pvalue | `atol = 1e-6` | jStat t-CDF vs scipy — implementations may differ slightly |
| cid_ce, mean_abs_change | `atol = 1e-10` | Simple sums + one sqrt |
| percentage_of_reoccurring_datapoints | `atol = 1e-10` | Integer counting / n |

Helper:
```typescript
function assertClose(actual: number, expected: number, atol: number, label: string): void {
  if (Math.abs(actual - expected) > atol)
    throw new Error(`${label}: expected ${expected}, got ${actual}, diff=${Math.abs(actual - expected)}`);
}

function assertExact(actual: number, expected: number, label: string): void {
  if (actual !== expected)
    throw new Error(`${label}: expected ${expected}, got ${actual}`);
}
```

---

## Golden Test Data

### Dataset 1: Normal random (seed=42)

Pre-generated in Python, stored as typed array literals in test file.

```python
# Generation script (Python):
import numpy as np
rng = np.random.default_rng(42)
for n in [10, 12, 15]:
    x = rng.normal(loc=7, scale=1, size=n)
    print(f"n={n}: [{', '.join(f'{v:.16e}' for v in x)}]")
```

Three series: `n=10`, `n=12`, `n=15`.

### Dataset 2: Constant series

```
x = [5.0, 5.0, 5.0, ..., 5.0]  (n=10)
```

Tests: std=0, variance=0, cid_ce=0, has_duplicate=true, all diffs=0,
percentage_of_reoccurring_datapoints=1.0, linear_trend slope=0.

### Dataset 3: Monotonic ramp

```
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  (n=10)
```

Tests: linear_trend slope=1.0, rvalue=1.0, pvalue≈0, mean_change=1.0,
first_location_of_minimum=0, last_location_of_maximum=0.9.

### Dataset 4: Alternating

```
x = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1]  (n=10)
```

Tests: mean=0, number_crossing_m(m=0)=9, high cid_ce, symmetry.

### Dataset 5: Step function

```
x = [0, 0, 0, 0, 0, 5, 5, 5, 5, 5]  (n=10)
```

Tests: has_duplicate=true, has_duplicate_max=true, has_duplicate_min=true,
percentage_of_reoccurring_datapoints=1.0, longest_strike_above_mean=5.

### Dataset 6: Single spike

```
x = [0, 0, 0, 0, 10, 0, 0, 0, 0, 0]  (n=10)
```

Tests: max=10, first_location_of_maximum=4/10=0.4, large variance,
ratio_beyond_r_sigma high.

### Dataset 7: Negative values

```
x = [-3.5, -1.2, 0.0, 2.1, -4.0, 1.5, -0.3, 3.3, -2.2, 0.8]  (n=10)
```

Tests: number_crossing_m(m=0) counts zero-crossings, negative min, mixed signs.

### Dataset 8: Near-constant with one outlier

```
x = [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 14.0]  (n=10)
```

Tests: ratio_beyond_r_sigma, skewed distribution, last_location_of_maximum = 0.9.

---

## Golden Value Generation

Python script to pre-compute all expected values for each dataset:

```python
import numpy as np
import json
from scipy.stats import linregress

R_VALUES = [0.05, 0.25, 0.45]
R_SIGMA_VALUES = [1, 1.5, 2]
M_VALUE = 0

def longest_strike(cond):
    max_len = cur = 0
    for v in cond:
        if v:
            cur += 1
            max_len = max(max_len, cur)
        else:
            cur = 0
    return max_len

def compute_golden(x):
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    mean = np.mean(x)
    std = np.std(x)       # ddof=0
    var = np.var(x)        # ddof=0
    med = np.median(x)
    diff = np.diff(x)
    rng = np.max(x) - np.min(x)

    # linear trend via scipy (exact reference)
    if n >= 2:
        lr = linregress(np.arange(n), x)
        slope, intercept, rvalue, pvalue, stderr = lr
    else:
        slope = intercept = rvalue = stderr = 0.0
        pvalue = 1.0

    g = {}

    # Basic
    g["mean"] = mean
    g["median"] = med
    g["minimum"] = float(np.min(x))
    g["maximum"] = float(np.max(x))
    g["length"] = n
    g["sum_values"] = float(np.sum(x))
    g["abs_energy"] = float(np.sum(x**2))
    g["root_mean_square"] = float(np.sqrt(np.mean(x**2)))
    g["standard_deviation"] = float(std)
    g["variance"] = float(var)

    # Changes
    g["mean_change"] = float((x[-1] - x[0]) / (n - 1)) if n > 1 else 0.0
    g["mean_abs_change"] = float(np.mean(np.abs(diff))) if n > 1 else 0.0
    g["absolute_sum_of_changes"] = float(np.sum(np.abs(diff))) if n > 1 else 0.0
    g["mean_second_derivative_central"] = (
        float((x[-1] - x[-2] - x[1] + x[0]) / (2 * (n - 2))) if n > 2 else 0.0
    )

    # Boolean (as float 1.0 / 0.0)
    g["variance_larger_than_standard_deviation"] = 1.0 if var > 1 else 0.0
    g["has_duplicate"] = 1.0 if len(np.unique(x)) != n else 0.0
    g["has_duplicate_max"] = 1.0 if np.sum(x == np.max(x)) >= 2 else 0.0
    g["has_duplicate_min"] = 1.0 if np.sum(x == np.min(x)) >= 2 else 0.0

    # Counts
    g["count_above_mean"] = int(np.sum(x > mean))
    g["count_below_mean"] = int(np.sum(x < mean))

    # Locations
    g["first_location_of_minimum"] = float(np.argmin(x) / n)
    g["first_location_of_maximum"] = float(np.argmax(x) / n)
    g["last_location_of_minimum"] = float((n - 1 - np.argmin(x[::-1])) / n)
    g["last_location_of_maximum"] = float((n - 1 - np.argmax(x[::-1])) / n)

    # Uniqueness
    unique, counts = np.unique(x, return_counts=True)
    g["percentage_of_reoccurring_datapoints_to_all_datapoints"] = float(
        np.sum(counts[counts > 1]) / n
    )
    g["ratio_value_number_to_time_series_length"] = float(len(unique) / n)

    # CID
    g["cid_ce__normalize_false"] = float(np.sqrt(np.sum(diff**2))) if n > 1 else 0.0
    if std > 0 and n > 1:
        g["cid_ce__normalize_true"] = float(np.sqrt(np.sum(diff**2)) / std)
    else:
        g["cid_ce__normalize_true"] = 0.0

    # Strikes
    g["longest_strike_above_mean"] = longest_strike(x > mean)
    g["longest_strike_below_mean"] = longest_strike(x < mean)

    # Linear trend (scipy reference)
    g["linear_trend__attr_slope"] = float(slope)
    g["linear_trend__attr_intercept"] = float(intercept)
    g["linear_trend__attr_rvalue"] = float(rvalue)
    g["linear_trend__attr_pvalue"] = float(pvalue)
    g["linear_trend__attr_stderr"] = float(stderr)

    # Parametric: symmetry_looking
    for r in R_VALUES:
        g[f"symmetry_looking__r_{r}"] = 1.0 if abs(mean - med) < r * rng else 0.0

    # Parametric: large_standard_deviation
    for r in R_VALUES:
        g[f"large_standard_deviation__r_{r}"] = 1.0 if std > r * rng else 0.0

    # Parametric: ratio_beyond_r_sigma
    for r in R_SIGMA_VALUES:
        if std > 0:
            g[f"ratio_beyond_r_sigma__r_{r}"] = float(np.sum(np.abs(x - mean) > r * std) / n)
        else:
            g[f"ratio_beyond_r_sigma__r_{r}"] = 0.0

    # number_crossing_m
    g[f"number_crossing_m__m_{M_VALUE}"] = int(np.sum(
        (x[:-1] - M_VALUE) * (x[1:] - M_VALUE) < 0
    )) if n > 1 else 0

    return g

# Generate and print for all datasets
datasets = {
    "normal_n10": np.random.default_rng(42).normal(7, 1, 10),
    "normal_n12": np.random.default_rng(42).normal(7, 1, 12),
    "normal_n15": np.random.default_rng(42).normal(7, 1, 15),
    "constant": np.full(10, 5.0),
    "ramp": np.arange(10, dtype=float),
    "alternating": np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=float),
    "step": np.array([0,0,0,0,0,5,5,5,5,5], dtype=float),
    "spike": np.array([0,0,0,0,10,0,0,0,0,0], dtype=float),
    "negative": np.array([-3.5,-1.2,0.0,2.1,-4.0,1.5,-0.3,3.3,-2.2,0.8]),
    "near_constant_outlier": np.array([7,7,7,7,7,7,7,7,7,14], dtype=float),
}

for name, data in datasets.items():
    golden = compute_golden(data)
    print(f"\n// === {name} (n={len(data)}) ===")
    print(f"const {name}_data = new Float64Array([{', '.join(str(v) for v in data)}]);")
    print(f"const {name}_expected = {json.dumps(golden, indent=2)};")
```

---

## Unit Tests

### Test group: `calculators.test.ts`

#### Pass 1 tests

```
test("pass1 — normal data n=10")
  input: normal_n10_data
  assert: mean, min, max, sum, sumSq, firstMinIdx, lastMaxIdx, minCount, maxCount
  tolerance: 1e-10

test("pass1 — constant series")
  input: constant_data
  assert: mean=5, min=5, max=5, std=0, variance=0, range=0
  assert: minCount=10, maxCount=10

test("pass1 — ramp")
  input: ramp_data
  assert: min=0, max=9, firstMinIdx=0, firstMaxIdx=9
  assert: lastMinIdx=0, lastMaxIdx=9
```

#### Pass 2 tests

```
test("pass2 — differences and uniqueness")
  input: normal_n10_data, mean from pass1
  assert: sumAbsDiff, crossingCount, uniqueCount
  assert: hasDuplicate=false (float data, extremely unlikely)

test("pass2 — constant series")
  assert: sumAbsDiff=0, sumSqDiff=0
  assert: hasDuplicate=true, uniqueCount=1
  assert: reoccurringDatapointCount=10

test("pass2 — alternating")
  assert: crossingCount=9 (every consecutive pair crosses 0)
  assert: sumAbsDiff = 9 * 2 = 18
```

#### Pass 3 tests

```
test("pass3 — count above/below mean")
  input: ramp_data, mean=4.5
  assert: countAboveMean=5 (5,6,7,8,9), countBelowMean=5 (0,1,2,3,4)
  note: 4.5 is not equal to any element, so no "equal to mean" cases

test("pass3 — ratio_beyond_r_sigma")
  input: spike_data, mean=1.0, std=...
  assert: ratio at r=1 > 0 (spike is far from mean)

test("pass3 — constant series")
  assert: countAboveMean=0, countBelowMean=0 (all equal to mean)
  assert: all ratio_beyond = 0 (std=0 guard)
```

#### Median tests

```
test("computeMedian — even n")
  input: [1, 2, 3, 4], n=4
  assert: 2.5

test("computeMedian — odd n")
  input: [1, 2, 3, 4, 5], n=5
  assert: 3.0

test("computeMedian — does not mutate input")
  input: Int32Array([3, 1, 2])
  assert: result=2, input unchanged
```

#### Linear trend tests

```
test("linearTrend — perfect line")
  input: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  assert: slope=1.0, intercept=0.0, rvalue=1.0, stderr≈0, pvalue≈0
  tolerance: slope/intercept/rvalue 1e-10, pvalue 1e-6

test("linearTrend — constant")
  input: [5, 5, 5, 5, 5]
  assert: slope=0, intercept=5.0, rvalue=0, pvalue=1.0, stderr=0

test("linearTrend — negative slope")
  input: [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
  assert: slope=-1.0, rvalue=-1.0

test("linearTrend — noisy data")
  input: normal_n10_data
  assert: matches scipy.stats.linregress golden values
  tolerance: slope/intercept/rvalue/stderr 1e-7, pvalue 1e-6

test("linearTrend — n=2")
  input: [1, 3]
  assert: slope=2.0, intercept=1.0, rvalue=1.0, stderr=0.0, pvalue=1.0
  note: n-2=0 degrees of freedom, so stderr=0 and pvalue=1.0

test("linearTrend — n=3")
  input: [1, 2, 4]
  assert: matches scipy golden values
  note: first non-trivial case for pvalue (df=1)
```

---

## Integration Tests

### Test group: `extract.test.ts`

#### Single sample, single column

```
test("extractFeatures — single sample, all features")
  input: DataFrame with 1 sample (id=1), 1 column ("value"), n=10
  data: normal_n10_data
  assert: FeatureMatrix has nSamples=1
  assert: each column name matches "{colName}__{featureName}" pattern
  assert: each feature value matches golden expected value
  tolerance: per-feature as defined above
```

#### Multiple samples, variable lengths

```
test("extractFeatures — 3 samples with different lengths")
  input: DataFrame with ids=[1,1,...,2,2,...,3,3,...], lengths 10, 12, 15
  data: normal_n10, normal_n12, normal_n15
  assert: FeatureMatrix has nSamples=3
  assert: sampleIds = Uint32Array([1, 2, 3])
  assert: each sample's features match its golden values
```

#### Multiple columns

```
test("extractFeatures — 2 columns (pH + concentration)")
  input: 1 sample, 2 columns
  assert: feature count = 2 × features_per_column
  assert: column names prefixed with "pH__" and "concentration__"
  assert: values match independently computed golden values
```

#### Edge case datasets

```
test("extractFeatures — constant series")
  data: [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
  assert: std=0, variance=0, cid_ce_norm=0, has_duplicate=1.0
  assert: percentage_of_reoccurring_datapoints=1.0
  assert: linear_trend slope=0, pvalue=1.0

test("extractFeatures — step function")
  data: [0,0,0,0,0,5,5,5,5,5]
  assert: has_duplicate=1.0, has_duplicate_max=1.0, has_duplicate_min=1.0
  assert: longest_strike_above_mean=5

test("extractFeatures — alternating")
  data: [1,-1,1,-1,1,-1,1,-1,1,-1]
  assert: mean=0.0
  assert: number_crossing_m__m_0 = 9
```

---

## Input Type Tests

### Test group: `types.test.ts`

Verify that all NumericArray types produce identical results:

```
test("Int32Array input produces same features as Float64Array")
  input: same integer data as Int32Array and Float64Array
  assert: all features match within 1e-10

test("Uint32Array input")
  input: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] as Uint32Array
  assert: same features as Float64Array version

test("Float32Array input")
  input: normal data as Float32Array
  assert: features match Float64 within Float32 precision (atol=1e-5)
  note: Float32 has ~7 decimal digits precision
```

---

## Validation Tests

### Test group: `validate.test.ts`

```
test("validate — passes on clean data")
  assert: no throw

test("validate — throws on NaN in Float64Array column")
  input: column with data[3] = NaN
  assert: throws Error mentioning "NaN" and index 3

test("validate — throws on Inf")
  input: column with data[0] = Infinity
  assert: throws Error

test("validate — skips NaN check for Int32Array")
  input: Int32Array column (cannot contain NaN by definition)
  assert: no throw

test("validate — throws on mismatched lengths")
  input: ids.length=10, time.length=10, column.data.length=9
  assert: throws Error

test("validate — throws on non-contiguous ids")
  input: ids = [1, 2, 1, 2] (interleaved)
  assert: throws Error about contiguous grouping

test("validate — not called when options.validate=false")
  input: data with NaN
  assert: no throw from extractFeatures (NaN propagates to results)
```

---

## Naming Convention Tests

### Test group: `naming.test.ts`

```
test("feature names follow tsfresh convention")
  assert: every FeatureColumn.name matches regex:
    /^[a-zA-Z0-9_]+__[a-zA-Z0-9_]+(__[a-zA-Z0-9_.]+)?$/

test("parameterized feature names include params")
  assert: "value__cid_ce__normalize_true" in names
  assert: "value__cid_ce__normalize_false" in names
  assert: "value__symmetry_looking__r_0.05" in names
  assert: "value__linear_trend__attr_slope" in names
  assert: "value__ratio_beyond_r_sigma__r_1.5" in names

test("multi-column names are prefixed correctly")
  input: columns named "pH" and "conc"
  assert: names starting with "pH__" and "conc__"
  assert: no names without prefix
```

---

## Performance Tests (optional)

```
test("extractFeatures — 1000 samples × 2 columns × n=15 under 100ms")
  input: synthetic DataFrame with 1000 samples
  assert: execution time < 100ms
  note: 1000 × 2 × ~49 features = ~98,000 feature values

test("extractFeatures — no memory leaks on repeated calls")
  run 100 times
  assert: no significant heap growth
```

---

## Cross-Validation Against tsfresh

For release validation (not CI), a Python script generates golden values using
actual tsfresh `extract_features()` and dumps to JSON. TypeScript tests load
this JSON and compare.

```python
# generate_golden.py
from tsfresh.feature_extraction import extract_features
import json, numpy as np

# ... same tier1_fc_params as in the Python test ...

datasets = { ... }
golden = {}
for name, x in datasets.items():
    df = pd.DataFrame({"id": 1, "time": np.arange(len(x)), "value": x})
    feats = extract_features(df, column_id="id", column_sort="time",
                             default_fc_parameters=tier1_fc_params,
                             disable_progressbar=True)
    golden[name] = feats.iloc[0].to_dict()

with open("golden_tsfresh.json", "w") as f:
    json.dump(golden, f, indent=2)
```

TypeScript test:
```typescript
import golden from './golden_tsfresh.json';

test.each(Object.keys(golden))("matches tsfresh for dataset %s", (name) => {
  const data = datasets[name];
  const result = extractFeatures(makeDataFrame(data));
  for (const col of result.columns) {
    const tsfreshKey = toTsfreshKey(col.name);
    assertClose(col.data[0], golden[name][tsfreshKey], tolerance(col.name), col.name);
  }
});
```
