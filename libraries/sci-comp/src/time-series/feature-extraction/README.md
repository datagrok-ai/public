# Time Series Feature Extraction

Pure TypeScript implementation of feature calculators
for [time series](https://en.wikipedia.org/wiki/Time_series). Part of `@datagrok-libraries/sci-comp`. Feature names follow naming conventions compatible with the [tsfresh](https://tsfresh.readthedocs.io/en/latest/index.html) library, enabling seamless interoperability.

Extracts **45 features** per value column:

| # | Category | Features | Count |
|---|----------|----------|-------|
| 1 | Basic statistics | mean, median, minimum, maximum, length, sum_values, abs_energy, root_mean_square, standard_deviation, variance | 10 |
| 2 | Boolean flags | variance_larger_than_standard_deviation, has_duplicate, has_duplicate_max, has_duplicate_min | 4 |
| 3 | Location | first_location_of_minimum, first_location_of_maximum, last_location_of_minimum, last_location_of_maximum | 4 |
| 4 | Change metrics | mean_change, mean_abs_change, absolute_sum_of_changes, mean_second_derivative_central | 4 |
| 5 | Uniqueness | percentage_of_reoccurring_datapoints_to_all_datapoints, ratio_value_number_to_time_series_length | 2 |
| 6 | Counting | count_above_mean, count_below_mean, longest_strike_above_mean, longest_strike_below_mean | 4 |
| 7 | Complexity | cid_ce (normalize=false), cid_ce (normalize=true) | 2 |
| 8 | Symmetry & spread | symmetry_looking (r=0.05, 0.25, 0.45), large_standard_deviation (r=0.05, 0.25, 0.45) | 6 |
| 9 | Linear trend | slope, intercept, rvalue, pvalue, stderr | 5 |
| 10 | Threshold | ratio_beyond_r_sigma (r=1, 1.5, 2) | 3 |
| 11 | Crossings | number_crossing_m (m=0) | 1 |
| | | **Total** | **45** |

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
