import * as DG from 'datagrok-api/dg';

/** The summary statistics compared between the actual and the predicted values, in the built-in tool's order. */
export const STATISTIC_IDS =
  ['min', 'max', 'med', 'avg', 'stdev', 'variance', 'skew', 'kurt', 'q1', 'q2', 'q3'] as const;
type StatisticId = typeof STATISTIC_IDS[number];

export const STATISTIC_DESCRIPTIONS: Record<StatisticId, string> = {
  min: 'Minimum value in the column. It represents the smallest value and should align with domain knowledge.',
  max: 'Maximum value in the column. It represents the largest value and should align with domain knowledge.',
  med: 'Median value in the column. It represents the middle value and is robust to outliers.',
  avg: 'Average value in the column. Large deviations from the median may indicate skewed data.',
  stdev: 'Standard deviation of values in the column. High values indicate large variability.',
  variance: 'Variance of values in the column. It measures the dispersion of the dataset.',
  skew: 'Skewness of values in the column. Positive skew indicates a long right tail; negative skew indicates a long ' +
    'left tail.',
  kurt: 'Kurtosis of values in the column. High kurtosis indicates heavy tails or outliers.',
  q1: 'First quartile (25th percentile) value in the column. It should be less than the median.',
  q2: 'Second quartile (median or 50th percentile) value in the column. It represents the middle of the dataset.',
  q3: 'Third quartile (75th percentile) value in the column. It should be greater than the median.',
};

interface StatisticRow { id: StatisticId; actual: number; predicted: number }

/** Each statistic of [actual] and of [predicted]; `q2` is the median, listed again as the built-in tool did. */
export function statisticsOf(actual: DG.Column, predicted: DG.Column): StatisticRow[] {
  const actualStats = DG.Stats.fromColumn(actual);
  const predictedStats = DG.Stats.fromColumn(predicted);
  return STATISTIC_IDS.map((id) => ({id, actual: actualStats[id], predicted: predictedStats[id]}));
}

/** A float column `residuals`: [prediction] minus [target] per row, empty where either value is. */
export function residualsOf(target: DG.Column, prediction: DG.Column): DG.Column {
  const values = Float32Array.from({length: target.length}, (_, i) =>
    target.isNone(i) || prediction.isNone(i) ? DG.FLOAT_NULL : prediction.getNumber(i) - target.getNumber(i));
  return DG.Column.fromFloat32Array('residuals', values);
}
