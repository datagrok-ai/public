# Time Series Feature Extraction

Pure TypeScript implementation of feature calculators
for [time series](https://en.wikipedia.org/wiki/Time_series). Part of `@datagrok-libraries/sci-comp`. Feature names follow naming conventions compatible with the [tsfresh](https://tsfresh.readthedocs.io/en/latest/index.html) library, enabling seamless interoperability.

Extracts **45 features** per value column:

### Basic statistics

| # | Feature | Description |
|---|--------|------------|
| 1 | mean | Average value of the time series |
| 2 | median | Median value |
| 3 | minimum | Minimum value |
| 4 | maximum | Maximum value |
| 5 | length | Length of the time series |
| 6 | sum_values | Sum of all values |
| 7 | abs_energy | Sum of squared values |
| 8 | root_mean_square | Square root of mean squared values |
| 9 | standard_deviation | Standard deviation |
| 10 | variance | Variance |

### Boolean flags

| # | Feature | Description |
|---|--------|------------|
| 11 | variance_larger_than_standard_deviation | Whether variance > standard deviation |
| 12 | has_duplicate | Whether there are duplicate values |
| 13 | has_duplicate_max | Whether maximum value appears more than once |
| 14 | has_duplicate_min | Whether minimum value appears more than once |

### Location

| # | Feature | Description |
|---|--------|------------|
| 15 | first_location_of_minimum | Index of first occurrence of minimum |
| 16 | first_location_of_maximum | Index of first occurrence of maximum |
| 17 | last_location_of_minimum | Index of last occurrence of minimum |
| 18 | last_location_of_maximum | Index of last occurrence of maximum |

### Change metrics

| # | Feature | Description |
|---|--------|------------|
| 19 | mean_change | Mean of consecutive differences |
| 20 | mean_abs_change | Mean of absolute consecutive differences |
| 21 | absolute_sum_of_changes | Sum of absolute differences |
| 22 | mean_second_derivative_central | Mean second derivative (central difference) |

### Uniqueness

| # | Feature | Description |
|---|--------|------------|
| 23 | percentage_of_reoccurring_datapoints_to_all_datapoints | Fraction of repeating values |
| 24 | ratio_value_number_to_time_series_length | Number of unique values divided by length |

### Counting

| # | Feature | Description |
|---|--------|------------|
| 25 | count_above_mean | Number of values above mean |
| 26 | count_below_mean | Number of values below mean |
| 27 | longest_strike_above_mean | Longest consecutive sequence above mean |
| 28 | longest_strike_below_mean | Longest consecutive sequence below mean |

### Complexity

| # | Feature | Description |
|---|--------|------------|
| 29 | cid_ce (normalize=false) | Complexity estimate (raw) |
| 30 | cid_ce (normalize=true) | Complexity estimate (normalized) |

### Symmetry & spread

| # | Feature | Description |
|---|--------|------------|
| 31 | symmetry_looking (r=0.05) | Symmetry check with tolerance r=0.05 |
| 32 | symmetry_looking (r=0.25) | Symmetry check with tolerance r=0.25 |
| 33 | symmetry_looking (r=0.45) | Symmetry check with tolerance r=0.45 |
| 34 | large_standard_deviation (r=0.05) | Std deviation larger than r * range |
| 35 | large_standard_deviation (r=0.25) | Same with r=0.25 |
| 36 | large_standard_deviation (r=0.45) | Same with r=0.45 |

### Linear trend

| # | Feature | Description |
|---|--------|------------|
| 37 | slope | Slope of linear regression |
| 38 | intercept | Intercept of regression line |
| 39 | rvalue | Correlation coefficient |
| 40 | pvalue | p-value of slope |
| 41 | stderr | Standard error of slope |

### Threshold

| # | Feature | Description |
|---|--------|------------|
| 42 | ratio_beyond_r_sigma (r=1) | Fraction beyond 1σ |
| 43 | ratio_beyond_r_sigma (r=1.5) | Fraction beyond 1.5σ |
| 44 | ratio_beyond_r_sigma (r=2) | Fraction beyond 2σ |

### Crossings

| # | Feature | Description |
|---|--------|------------|
| 45 | number_crossing_m (m=0) | Number of crossings of level m=0 |

**Total number of features: 45**

## Installation

```bash
npm install @datagrok-libraries/sci-comp
```

## Usage

```typescript
import {timeSeries} from '@datagrok-libraries/sci-comp';

const {extractFeatures, validate} = timeSeries.featureExtraction;
```

### Input format

Data is packed into a `TimeSeriesDataFrame` — a columnar structure where
`ids` groups rows into samples and `time` provides the time axis:

```typescript
import type {TimeSeriesDataFrame} from '@datagrok-libraries/sci-comp';

const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1]),
  time: new Float64Array([0, 1, 2, 3, 4]),
  rowCount: 5,
  columns: [{
    name: 'temperature',
    data: new Float64Array([20.1, 20.5, 21.0, 20.8, 21.3]),
  }],
};
```

All rows for the same `id` must be contiguous. Multiple value columns can be
provided — each produces its own set of features.

### Single sample, single column

```typescript
const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
  time: new Float64Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
  rowCount: 10,
  columns: [{
    name: 'temperature',
    data: new Float64Array([20.1, 20.5, 21.0, 20.8, 21.3, 21.7, 22.0, 21.5, 21.8, 22.2]),
  }],
};

const result = extractFeatures(df);
// result.nSamples === 1
// result.columns.length === 45
// result.columns[0].name === 'temperature__mean'
// result.columns[0].data[0] === 21.29
```

### Multiple samples

Pack several time series into one DataFrame using different ids.
Features are extracted for each sample independently:

```typescript
const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]),
  time: new Float64Array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4]),
  rowCount: 15,
  columns: [{
    name: 'pH',
    data: new Float64Array([
      7.0, 7.1, 7.2, 7.3, 7.4,   // sample 1: rising
      7.0, 7.0, 7.0, 7.0, 7.0,   // sample 2: constant
      6.8, 7.3, 6.9, 7.5, 7.1,   // sample 3: noisy
    ]),
  }],
};

const result = extractFeatures(df);
// result.nSamples === 3
// result.sampleIds === Uint32Array([1, 2, 3])
// result.columns[0].data[0] — mean for sample 1
// result.columns[0].data[1] — mean for sample 2
// result.columns[0].data[2] — mean for sample 3
```

### Multiple columns

When several measurement channels are provided, each column produces its own
set of features prefixed with the column name:

```typescript
const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
  time: new Float64Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
  rowCount: 10,
  columns: [
    {name: 'pH', data: new Float64Array([7.0, 7.1, 7.2, 7.3, 7.4, 7.3, 7.2, 7.1, 7.0, 6.9])},
    {name: 'concentration', data: new Float64Array([0.5, 0.8, 1.2, 1.5, 2.0, 2.3, 2.8, 3.1, 3.5, 4.0])},
  ],
};

const result = extractFeatures(df);
// result.columns.length === 90  (45 per column)
// 'pH__mean', 'pH__linear_trend__attr_slope', ...
// 'concentration__mean', 'concentration__linear_trend__attr_slope', ...
```

### Validation

Pass `validate: true` to check the input for NaN, Infinity, length mismatches,
and non-contiguous ids before extraction:

```typescript
extractFeatures(df, {validate: true});
// throws Error if data contains NaN, Inf, or structural issues

validate(df);
// standalone validation without extraction
```

### Feature naming convention

Feature names follow the tsfresh convention:

```
{column}__{feature}
{column}__{feature}__{param}_{value}
```

Examples:

```
temperature__mean
temperature__cid_ce__normalize_true
temperature__symmetry_looking__r_0.05
temperature__linear_trend__attr_slope
temperature__ratio_beyond_r_sigma__r_2
temperature__number_crossing_m__m_0
```

### Supported input types

Value columns accept any typed numeric array:

| Type | Precision |
|------|-----------|
| `Float64Array` | Full precision |
| `Float32Array` | ~7 decimal digits |
| `Int32Array` | Integer data |
| `Uint32Array` | Unsigned integer data |

## Running examples

```bash
# Single sample, single column
npx tsx src/time-series/feature-extraction/examples/single-sample.ts

# Multiple samples — side-by-side comparison
npx tsx src/time-series/feature-extraction/examples/multiple-samples.ts

# Multiple columns — cross-channel comparison
npx tsx src/time-series/feature-extraction/examples/multiple-columns.ts
```
## References

1. Christ, M., Braun, N., Neuffer, J. and Kempa-Liehr A.W. (2018). Time Series FeatuRe Extraction on basis of Scalable Hypothesis tests (tsfresh – A Python package). Neurocomputing 307 (2018) 72-77, https://doi.org/10.1016/j.neucom.2018.03.067.
  
2. Herff, C., Krusienski, D.J. (2019). Extracting Features from Time Series. In: Kubben, P., Dumontier, M., Dekker, A. (eds) Fundamentals of Clinical Data Science. Springer, Cham. https://doi.org/10.1007/978-3-319-99713-1_7
