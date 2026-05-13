/**
 * Feature extraction orchestrator.
 *
 * Wires together index building, multi-pass calculators, median, and linear
 * trend into a single `extractFeatures()` entry point.
 */

import {TimeSeriesDataFrame, FeatureColumn, FeatureMatrix, ExtractOptions} from './types';
import {buildIndex} from './dataframe';
import {validate} from './dataframe';
import {pass1, pass2, pass3, computeMedian} from './calculators';
import {linearTrend} from './linear-trend';

/* ================================================================== */
/*  Configuration                                                      */
/* ================================================================== */

/** Thresholds for symmetry_looking and large_standard_deviation. */
const R_VALUES = [0.05, 0.25, 0.45];

/** Thresholds for ratio_beyond_r_sigma. */
const R_SIGMA_VALUES = [1, 1.5, 2];

/** Crossing thresholds for number_crossing_m. */
const M_VALUES = [0];

/* ================================================================== */
/*  Feature naming                                                     */
/* ================================================================== */

/**
 * Build the ordered list of feature names for a single value column.
 *
 * Names follow the tsfresh convention: `{col}__{feature}` or
 * `{col}__{feature}__{param}_{value}`.
 */
function buildFeatureNames(colName: string): string[] {
  const names: string[] = [];

  // Simple parameterless features
  names.push(`${colName}__mean`);
  names.push(`${colName}__median`);
  names.push(`${colName}__minimum`);
  names.push(`${colName}__maximum`);
  names.push(`${colName}__length`);
  names.push(`${colName}__sum_values`);
  names.push(`${colName}__abs_energy`);
  names.push(`${colName}__root_mean_square`);
  names.push(`${colName}__standard_deviation`);
  names.push(`${colName}__variance`);
  names.push(`${colName}__variance_larger_than_standard_deviation`);
  names.push(`${colName}__first_location_of_minimum`);
  names.push(`${colName}__first_location_of_maximum`);
  names.push(`${colName}__last_location_of_minimum`);
  names.push(`${colName}__last_location_of_maximum`);
  names.push(`${colName}__has_duplicate_max`);
  names.push(`${colName}__has_duplicate_min`);
  names.push(`${colName}__mean_change`);
  names.push(`${colName}__mean_abs_change`);
  names.push(`${colName}__absolute_sum_of_changes`);
  names.push(`${colName}__mean_second_derivative_central`);
  names.push(`${colName}__has_duplicate`);
  names.push(`${colName}__percentage_of_reoccurring_datapoints_to_all_datapoints`);
  names.push(`${colName}__ratio_value_number_to_time_series_length`);
  names.push(`${colName}__longest_strike_above_mean`);
  names.push(`${colName}__longest_strike_below_mean`);
  names.push(`${colName}__count_above_mean`);
  names.push(`${colName}__count_below_mean`);

  // cid_ce parametric
  names.push(`${colName}__cid_ce__normalize_false`);
  names.push(`${colName}__cid_ce__normalize_true`);

  // symmetry_looking parametric
  for (const r of R_VALUES)
    names.push(`${colName}__symmetry_looking__r_${r}`);

  // large_standard_deviation parametric
  for (const r of R_VALUES)
    names.push(`${colName}__large_standard_deviation__r_${r}`);

  // linear_trend parametric
  names.push(`${colName}__linear_trend__attr_slope`);
  names.push(`${colName}__linear_trend__attr_intercept`);
  names.push(`${colName}__linear_trend__attr_rvalue`);
  names.push(`${colName}__linear_trend__attr_pvalue`);
  names.push(`${colName}__linear_trend__attr_stderr`);

  // ratio_beyond_r_sigma parametric
  for (const r of R_SIGMA_VALUES)
    names.push(`${colName}__ratio_beyond_r_sigma__r_${r}`);

  // number_crossing_m parametric
  for (const m of M_VALUES)
    names.push(`${colName}__number_crossing_m__m_${m}`);

  return names;
}

/* ================================================================== */
/*  Per-sample computation                                             */
/* ================================================================== */

/**
 * Compute all feature values for a single sample slice.
 *
 * Runs pass 1–3, median, and linear trend, then assembles the results
 * into a flat array matching the order of {@link buildFeatureNames}.
 */
function computeSampleFeatures(
  data: Float64Array | Float32Array | Int32Array | Uint32Array,
  n: number,
  medianBuf: Float64Array,
): number[] {
  if (n === 0)
    return new Array(buildFeatureNames('_').length).fill(NaN);

  const p1 = pass1(data, n);
  const median = computeMedian(data, n, medianBuf);

  const p2Result = n > 1 ? pass2(data, n, p1.mean, M_VALUES[0]) : null;
  const p3Result = pass3(data, n, p1.mean, p1.std, R_SIGMA_VALUES);
  const lt = linearTrend(data, n);

  const values: number[] = [];

  // Simple parameterless
  values.push(p1.mean);
  values.push(median);
  values.push(p1.min);
  values.push(p1.max);
  values.push(p1.n);
  values.push(p1.sum);
  values.push(p1.sumSq);
  values.push(Math.sqrt(p1.sumSq / n));
  values.push(p1.std);
  values.push(p1.variance);
  values.push(p1.variance > 1 ? 1.0 : 0.0);
  values.push(p1.firstMinIdx / n);
  values.push(p1.firstMaxIdx / n);
  values.push((p1.lastMinIdx + 1) / n);
  values.push((p1.lastMaxIdx + 1) / n);
  values.push(p1.maxCount >= 2 ? 1.0 : 0.0);
  values.push(p1.minCount >= 2 ? 1.0 : 0.0);

  // mean_change
  if (n > 1)
    values.push((data[n - 1] - data[0]) / (n - 1));
  else
    values.push(0);

  // mean_abs_change
  values.push(p2Result ? p2Result.sumAbsDiff / (n - 1) : 0);

  // absolute_sum_of_changes
  values.push(p2Result ? p2Result.sumAbsDiff : 0);

  // mean_second_derivative_central
  if (n > 2)
    values.push((data[n - 1] - data[n - 2] - data[1] + data[0]) / (2 * (n - 2)));
  else
    values.push(0);

  // has_duplicate
  values.push(p2Result ? (p2Result.hasDuplicate ? 1.0 : 0.0) : 0.0);

  // percentage_of_reoccurring_datapoints_to_all_datapoints
  values.push(p2Result ? p2Result.reoccurringDatapointCount / n : 0);

  // ratio_value_number_to_time_series_length
  values.push(p2Result ? p2Result.uniqueCount / n : 1);

  // longest_strike_above_mean, longest_strike_below_mean
  values.push(p2Result ? p2Result.longestStrikeAbove : 0);
  values.push(p2Result ? p2Result.longestStrikeBelow : 0);

  // count_above_mean, count_below_mean
  values.push(p3Result.countAboveMean);
  values.push(p3Result.countBelowMean);

  // cid_ce normalize=false, normalize=true
  const cidFalse = p2Result ? Math.sqrt(p2Result.sumSqDiff) : 0;
  values.push(cidFalse);
  values.push(p1.std > 0 ? cidFalse / p1.std : 0);

  // symmetry_looking
  for (const r of R_VALUES)
    values.push(Math.abs(p1.mean - median) < r * p1.range ? 1.0 : 0.0);

  // large_standard_deviation
  for (const r of R_VALUES)
    values.push(p1.std > r * p1.range ? 1.0 : 0.0);

  // linear_trend
  values.push(lt.slope);
  values.push(lt.intercept);
  values.push(lt.rvalue);
  values.push(lt.pvalue);
  values.push(lt.stderr);

  // ratio_beyond_r_sigma
  for (let j = 0; j < R_SIGMA_VALUES.length; j++)
    values.push(p3Result.beyondRSigmaCounts[j] / n);

  // number_crossing_m
  for (let _j = 0; _j < M_VALUES.length; _j++)
    values.push(p2Result ? p2Result.crossingCount : 0);

  return values;
}

/* ================================================================== */
/*  Public API                                                         */
/* ================================================================== */

/**
 * Extract tsfresh-compatible features from a time-series DataFrame.
 *
 * For each sample (grouped by `df.ids`) and each value column, runs
 * three passes plus median and linear trend, producing ~45 features
 * per column.
 *
 * @param df      Input DataFrame with ids, time, and value columns.
 * @param options Optional settings (e.g. `{validate: true}`).
 * @returns Feature matrix with one row per sample and one column per feature.
 */
export function extractFeatures(
  df: TimeSeriesDataFrame, options?: ExtractOptions,
): FeatureMatrix {
  if (options?.validate) validate(df);

  const index = buildIndex(df);
  const nSamples = index.size;
  const sampleIds = new Uint32Array(nSamples);

  // Build feature names for all columns
  const allFeatureNames: string[] = [];
  for (const col of df.columns) {
    const names = buildFeatureNames(col.name);
    allFeatureNames.push(...names);
  }

  const featuresPerColumn = buildFeatureNames('_').length;

  // Allocate result columns
  const resultColumns: FeatureColumn[] = allFeatureNames.map((name) => ({
    name,
    data: new Float64Array(nSamples),
  }));

  // Find max sample length for median buffer
  let maxLen = 0;
  for (const range of index.values())
    if (range.length > maxLen) maxLen = range.length;
  const medianBuf = new Float64Array(maxLen);

  let si = 0;
  for (const [id, range] of index) {
    sampleIds[si] = id;

    for (let ci = 0; ci < df.columns.length; ci++) {
      const col = df.columns[ci];
      const slice = col.data.subarray(range.start, range.start + range.length);
      const featureValues = computeSampleFeatures(slice, range.length, medianBuf);
      const offset = ci * featuresPerColumn;

      for (let fi = 0; fi < featureValues.length; fi++)
        resultColumns[offset + fi].data[si] = featureValues[fi];
    }

    si++;
  }

  return {sampleIds, columns: resultColumns, nSamples};
}
