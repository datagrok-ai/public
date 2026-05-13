import {TimeSeriesDataFrame} from '../types';

// === Test datasets ===

export const constantData = new Float64Array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);
export const rampData = new Float64Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
export const alternatingData = new Float64Array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1]);
export const stepData = new Float64Array([0, 0, 0, 0, 0, 5, 5, 5, 5, 5]);
export const spikeData = new Float64Array([0, 0, 0, 0, 10, 0, 0, 0, 0, 0]);
export const negativeData = new Float64Array([-3.5, -1.2, 0.0, 2.1, -4.0, 1.5, -0.3, 3.3, -2.2, 0.8]);
export const nearConstantOutlierData = new Float64Array([7, 7, 7, 7, 7, 7, 7, 7, 7, 14]);

// === Golden expected values ===
// Computed using the same formulas as tsfresh 0.21.x / scipy.stats.linregress

export const constantExpected: Record<string, number> = {
  'mean': 5, 'median': 5, 'minimum': 5, 'maximum': 5, 'length': 10,
  'sum_values': 50, 'abs_energy': 250, 'root_mean_square': 5,
  'standard_deviation': 0, 'variance': 0,
  'variance_larger_than_standard_deviation': 0,
  'first_location_of_minimum': 0, 'first_location_of_maximum': 0,
  'last_location_of_minimum': 1.0, 'last_location_of_maximum': 1.0,
  'has_duplicate_max': 1, 'has_duplicate_min': 1,
  'mean_change': 0, 'mean_abs_change': 0, 'absolute_sum_of_changes': 0,
  'mean_second_derivative_central': 0,
  'has_duplicate': 1,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 1,
  'ratio_value_number_to_time_series_length': 0.1,
  'longest_strike_above_mean': 0, 'longest_strike_below_mean': 0,
  'count_above_mean': 0, 'count_below_mean': 0,
  'cid_ce__normalize_false': 0, 'cid_ce__normalize_true': 0,
  'symmetry_looking__r_0.05': 0, 'symmetry_looking__r_0.25': 0, 'symmetry_looking__r_0.45': 0,
  'large_standard_deviation__r_0.05': 0, 'large_standard_deviation__r_0.25': 0,
  'large_standard_deviation__r_0.45': 0,
  'linear_trend__attr_slope': 0, 'linear_trend__attr_intercept': 5,
  'linear_trend__attr_rvalue': 0, 'linear_trend__attr_pvalue': 1,
  'linear_trend__attr_stderr': 0,
  'ratio_beyond_r_sigma__r_1': 0, 'ratio_beyond_r_sigma__r_1.5': 0, 'ratio_beyond_r_sigma__r_2': 0,
  'number_crossing_m__m_0': 0,
};

export const rampExpected: Record<string, number> = {
  'mean': 4.5, 'median': 4.5, 'minimum': 0, 'maximum': 9, 'length': 10,
  'sum_values': 45, 'abs_energy': 285, 'root_mean_square': 5.338539126015656,
  'standard_deviation': 2.8722813232690143, 'variance': 8.25,
  'variance_larger_than_standard_deviation': 1,
  'first_location_of_minimum': 0, 'first_location_of_maximum': 0.9,
  'last_location_of_minimum': 0.1, 'last_location_of_maximum': 1.0,
  'has_duplicate_max': 0, 'has_duplicate_min': 0,
  'mean_change': 1, 'mean_abs_change': 1, 'absolute_sum_of_changes': 9,
  'mean_second_derivative_central': 0,
  'has_duplicate': 0,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 0,
  'ratio_value_number_to_time_series_length': 1,
  'longest_strike_above_mean': 5, 'longest_strike_below_mean': 5,
  'count_above_mean': 5, 'count_below_mean': 5,
  'cid_ce__normalize_false': 3, 'cid_ce__normalize_true': 1.044465935734187,
  'symmetry_looking__r_0.05': 1, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 0,
  'linear_trend__attr_slope': 1, 'linear_trend__attr_intercept': 0,
  'linear_trend__attr_rvalue': 1, 'linear_trend__attr_pvalue': 0,
  'linear_trend__attr_stderr': 0,
  'ratio_beyond_r_sigma__r_1': 0.4, 'ratio_beyond_r_sigma__r_1.5': 0.2,
  'ratio_beyond_r_sigma__r_2': 0,
  'number_crossing_m__m_0': 1,
};

export const alternatingExpected: Record<string, number> = {
  'mean': 0, 'median': 0, 'minimum': -1, 'maximum': 1, 'length': 10,
  'sum_values': 0, 'abs_energy': 10, 'root_mean_square': 1,
  'standard_deviation': 1, 'variance': 1,
  'variance_larger_than_standard_deviation': 0,
  'first_location_of_minimum': 0.1, 'first_location_of_maximum': 0,
  'last_location_of_minimum': 1.0, 'last_location_of_maximum': 0.9,
  'has_duplicate_max': 1, 'has_duplicate_min': 1,
  'mean_change': -0.2222222222222222, 'mean_abs_change': 2,
  'absolute_sum_of_changes': 18, 'mean_second_derivative_central': 0,
  'has_duplicate': 1,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 1,
  'ratio_value_number_to_time_series_length': 0.2,
  'longest_strike_above_mean': 1, 'longest_strike_below_mean': 1,
  'count_above_mean': 5, 'count_below_mean': 5,
  'cid_ce__normalize_false': 6, 'cid_ce__normalize_true': 6,
  'symmetry_looking__r_0.05': 1, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 1,
  'linear_trend__attr_slope': -0.06060606060606061,
  'linear_trend__attr_intercept': 0.2727272727272727,
  'linear_trend__attr_rvalue': -0.17407765595569785,
  'linear_trend__attr_pvalue': 0.6305360755570297,
  'linear_trend__attr_stderr': 0.12121212121212122,
  'ratio_beyond_r_sigma__r_1': 0, 'ratio_beyond_r_sigma__r_1.5': 0,
  'ratio_beyond_r_sigma__r_2': 0,
  'number_crossing_m__m_0': 9,
};

export const stepExpected: Record<string, number> = {
  'mean': 2.5, 'median': 2.5, 'minimum': 0, 'maximum': 5, 'length': 10,
  'sum_values': 25, 'abs_energy': 125, 'root_mean_square': 3.5355339059327378,
  'standard_deviation': 2.5, 'variance': 6.25,
  'variance_larger_than_standard_deviation': 1,
  'first_location_of_minimum': 0, 'first_location_of_maximum': 0.5,
  'last_location_of_minimum': 0.5, 'last_location_of_maximum': 1.0,
  'has_duplicate_max': 1, 'has_duplicate_min': 1,
  'mean_change': 0.5555555555555556, 'mean_abs_change': 0.5555555555555556,
  'absolute_sum_of_changes': 5, 'mean_second_derivative_central': 0,
  'has_duplicate': 1,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 1,
  'ratio_value_number_to_time_series_length': 0.2,
  'longest_strike_above_mean': 5, 'longest_strike_below_mean': 5,
  'count_above_mean': 5, 'count_below_mean': 5,
  'cid_ce__normalize_false': 5, 'cid_ce__normalize_true': 2,
  'symmetry_looking__r_0.05': 1, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 1,
  'linear_trend__attr_slope': 0.7575757575757576,
  'linear_trend__attr_intercept': -0.9090909090909092,
  'linear_trend__attr_rvalue': 0.8703882797784891,
  'linear_trend__attr_pvalue': 0.0010528257934054874,
  'linear_trend__attr_stderr': 0.1515151515151515,
  'ratio_beyond_r_sigma__r_1': 0, 'ratio_beyond_r_sigma__r_1.5': 0,
  'ratio_beyond_r_sigma__r_2': 0,
  'number_crossing_m__m_0': 1,
};

export const spikeExpected: Record<string, number> = {
  'mean': 1, 'median': 0, 'minimum': 0, 'maximum': 10, 'length': 10,
  'sum_values': 10, 'abs_energy': 100, 'root_mean_square': 3.1622776601683795,
  'standard_deviation': 3, 'variance': 9,
  'variance_larger_than_standard_deviation': 1,
  'first_location_of_minimum': 0, 'first_location_of_maximum': 0.4,
  'last_location_of_minimum': 1.0, 'last_location_of_maximum': 0.5,
  'has_duplicate_max': 0, 'has_duplicate_min': 1,
  'mean_change': 0, 'mean_abs_change': 2.2222222222222223,
  'absolute_sum_of_changes': 20, 'mean_second_derivative_central': 0,
  'has_duplicate': 1,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 0.9,
  'ratio_value_number_to_time_series_length': 0.2,
  'longest_strike_above_mean': 1, 'longest_strike_below_mean': 5,
  'count_above_mean': 1, 'count_below_mean': 9,
  'cid_ce__normalize_false': 14.142135623730951,
  'cid_ce__normalize_true': 4.714045207910317,
  'symmetry_looking__r_0.05': 0, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 0,
  'linear_trend__attr_slope': -0.06060606060606061,
  'linear_trend__attr_intercept': 1.2727272727272727,
  'linear_trend__attr_rvalue': -0.058025885318565944,
  'linear_trend__attr_pvalue': 0.873494892371369,
  'linear_trend__attr_stderr': 0.3686522745635285,
  'ratio_beyond_r_sigma__r_1': 0.1, 'ratio_beyond_r_sigma__r_1.5': 0.1,
  'ratio_beyond_r_sigma__r_2': 0.1,
  'number_crossing_m__m_0': 2,
};

export const negativeExpected: Record<string, number> = {
  'mean': -0.35, 'median': -0.15, 'minimum': -4, 'maximum': 3.3, 'length': 10,
  'sum_values': -3.5, 'abs_energy': 52.81000000000001,
  'root_mean_square': 2.2980426453832403,
  'standard_deviation': 2.271233145231902, 'variance': 5.158499999999999,
  'variance_larger_than_standard_deviation': 1,
  'first_location_of_minimum': 0.4, 'first_location_of_maximum': 0.7,
  'last_location_of_minimum': 0.5, 'last_location_of_maximum': 0.8,
  'has_duplicate_max': 0, 'has_duplicate_min': 0,
  'mean_change': 0.47777777777777775, 'mean_abs_change': 3.4555555555555557,
  'absolute_sum_of_changes': 31.1, 'mean_second_derivative_central': 0.04375000000000001,
  'has_duplicate': 0,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 0,
  'ratio_value_number_to_time_series_length': 1,
  'longest_strike_above_mean': 3, 'longest_strike_below_mean': 2,
  'count_above_mean': 6, 'count_below_mean': 4,
  'cid_ce__normalize_false': 11.577996372429903,
  'cid_ce__normalize_true': 5.097669693988084,
  'symmetry_looking__r_0.05': 1, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 0,
  'linear_trend__attr_slope': 0.2818181818181818,
  'linear_trend__attr_intercept': -1.6181818181818182,
  'linear_trend__attr_rvalue': 0.3563971853322638,
  'linear_trend__attr_pvalue': 0.3120888302970459,
  'linear_trend__attr_stderr': 0.26121141812462506,
  'ratio_beyond_r_sigma__r_1': 0.4, 'ratio_beyond_r_sigma__r_1.5': 0.2,
  'ratio_beyond_r_sigma__r_2': 0,
  'number_crossing_m__m_0': 7,
};

export const nearConstantOutlierExpected: Record<string, number> = {
  'mean': 7.7, 'median': 7, 'minimum': 7, 'maximum': 14, 'length': 10,
  'sum_values': 77, 'abs_energy': 637, 'root_mean_square': 7.981227975693966,
  'standard_deviation': 2.1, 'variance': 4.41,
  'variance_larger_than_standard_deviation': 1,
  'first_location_of_minimum': 0, 'first_location_of_maximum': 0.9,
  'last_location_of_minimum': 0.9, 'last_location_of_maximum': 1.0,
  'has_duplicate_max': 0, 'has_duplicate_min': 1,
  'mean_change': 0.7777777777777778, 'mean_abs_change': 0.7777777777777778,
  'absolute_sum_of_changes': 7, 'mean_second_derivative_central': 0.4375,
  'has_duplicate': 1,
  'percentage_of_reoccurring_datapoints_to_all_datapoints': 0.9,
  'ratio_value_number_to_time_series_length': 0.2,
  'longest_strike_above_mean': 1, 'longest_strike_below_mean': 9,
  'count_above_mean': 1, 'count_below_mean': 9,
  'cid_ce__normalize_false': 7, 'cid_ce__normalize_true': 3.333333333333333,
  'symmetry_looking__r_0.05': 0, 'symmetry_looking__r_0.25': 1, 'symmetry_looking__r_0.45': 1,
  'large_standard_deviation__r_0.05': 1, 'large_standard_deviation__r_0.25': 1,
  'large_standard_deviation__r_0.45': 0,
  'linear_trend__attr_slope': 0.38181818181818183,
  'linear_trend__attr_intercept': 5.981818181818182,
  'linear_trend__attr_rvalue': 0.5222329678670934,
  'linear_trend__attr_pvalue': 0.1215029188171235,
  'linear_trend__attr_stderr': 0.22044283005422083,
  'ratio_beyond_r_sigma__r_1': 0.1, 'ratio_beyond_r_sigma__r_1.5': 0.1,
  'ratio_beyond_r_sigma__r_2': 0.1,
  'number_crossing_m__m_0': 0,
};

// === Tolerance helpers ===

const TREND_FEATURES = new Set([
  'linear_trend__attr_slope', 'linear_trend__attr_intercept',
  'linear_trend__attr_rvalue', 'linear_trend__attr_stderr',
]);
const TREND_PVALUE = 'linear_trend__attr_pvalue';
const EXACT_FEATURES = new Set(['length']);
const BOOL_FEATURES = new Set([
  'variance_larger_than_standard_deviation', 'has_duplicate',
  'has_duplicate_max', 'has_duplicate_min',
]);

export function toleranceFor(featureName: string): number | 'exact' {
  // Strip column prefix to get bare feature name
  const bare = featureName.includes('__') ?
    featureName.substring(featureName.indexOf('__') + 2) : featureName;

  if (EXACT_FEATURES.has(bare)) return 'exact';
  if (BOOL_FEATURES.has(bare)) return 'exact';
  if (bare.startsWith('symmetry_looking') || bare.startsWith('large_standard_deviation'))
    return 'exact';
  if (bare === TREND_PVALUE) return 1e-6;
  if (TREND_FEATURES.has(bare)) return 1e-7;
  return 1e-10;
}

export function assertClose(actual: number, expected: number, atol: number, label: string): void {
  if (Math.abs(actual - expected) > atol) {
    throw new Error(
      `${label}: expected ${expected}, got ${actual}, diff=${Math.abs(actual - expected)}`,
    );
  }
}

export function assertExact(actual: number, expected: number, label: string): void {
  if (actual !== expected)
    throw new Error(`${label}: expected ${expected}, got ${actual}`);
}

// === DataFrame builder ===

export function makeDataFrame(
  data: Float64Array, colName = 'value',
): TimeSeriesDataFrame {
  const n = data.length;
  return {
    ids: new Uint32Array(n).fill(1),
    time: new Float64Array(Array.from({length: n}, (_, i) => i)),
    rowCount: n,
    columns: [{name: colName, data}],
  };
}

export function makeMultiSampleDataFrame(
  samples: {id: number; data: Float64Array}[],
  colName = 'value',
): TimeSeriesDataFrame {
  const totalRows = samples.reduce((s, d) => s + d.data.length, 0);
  const ids = new Uint32Array(totalRows);
  const time = new Float64Array(totalRows);
  const values = new Float64Array(totalRows);

  let offset = 0;
  for (const sample of samples) {
    for (let i = 0; i < sample.data.length; i++) {
      ids[offset + i] = sample.id;
      time[offset + i] = i;
      values[offset + i] = sample.data[i];
    }
    offset += sample.data.length;
  }

  return {
    ids, time, rowCount: totalRows,
    columns: [{name: colName, data: values}],
  };
}
