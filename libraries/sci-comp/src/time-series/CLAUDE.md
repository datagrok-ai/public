# CLAUDE.md — time-series

Currently a single submodule: `feature-extraction`. No dependency on `datagrok-api`.

## Architecture

```
src/time-series/
  index.ts                          # Re-exports feature-extraction
  feature-extraction/
    types.ts                        # NumericArray, TimeSeriesDataFrame, TimeSeriesColumn, FeatureColumn, FeatureMatrix, ExtractOptions
    extract.ts                      # extractFeatures() — main entry, orchestrates passes + median + linear trend
    calculators.ts                  # pass1 (basic stats), pass2 (diffs, uniqueness, strikes), pass3 (threshold), computeMedian
    linear-trend.ts                 # linearTrend() — slope, intercept, rvalue, pvalue, stderr via least-squares
    dataframe.ts                    # buildIndex(), validate() — sample grouping and input validation
    stats-utils.ts                  # Statistical helpers
    __tests__/                      # Tests by module: calculators, extract, linear-trend, naming, types, validate
    examples/                       # single-sample.ts, multiple-samples.ts, multiple-columns.ts
```

## Key design patterns

- **tsfresh-compatible naming**: produces ~45 features per value column following the tsfresh convention — `{col}__{feature}` or `{col}__{feature}__{param}_{value}`. When adding a feature, match the upstream tsfresh name exactly so downstream tooling can interop. Don't invent new naming styles.

- **Multi-pass calculator architecture** (`calculators.ts`): computations are split into three passes to minimize iterations over each sample slice — pass1 (basic stats), pass2 (diffs / uniqueness / crossings / strikes), pass3 (threshold-based counts), plus separate median and linear trend. New features go in the pass that already iterates over the data they need; only add a fourth pass if the access pattern is genuinely different.

- **Columnar I/O**: input is a `TimeSeriesDataFrame` with `ids` (sample grouping), `time`, and value `columns`. Output is a `FeatureMatrix` — one row per sample, one `FeatureColumn` per feature. Don't introduce row-oriented helpers.

- **`NumericArray`**: input data may be `Int32Array | Uint32Array | Float32Array | Float64Array`. Internal accumulators stay in `Float64Array` regardless.

- **Contiguous id grouping invariant**: all rows for the same sample id MUST be contiguous in input (sorted by id). `validate()` enforces this — do not relax. The pass loops depend on it for O(n) slicing.

## Tests

Organized by module rather than by feature: `calculators.test.ts`, `extract.test.ts`, `linear-trend.test.ts`, `naming.test.ts`, `types.test.ts`, `validate.test.ts`.
