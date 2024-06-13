// Regression tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_fitLinearRegressionParams, _fitLinearRegressionParamsWithDataNormalizing} from '../wasm/EDAAPI';

/** Computes coefficients of linear regression */
export function computeLinRegressionCoefs(features: DG.ColumnList, targets: DG.Column): Float32Array {
  const featuresCount = features.length;

  const yAvg = targets.stats.avg;
  const yStdev = targets.stats.stdev;

  const params = new Float32Array(featuresCount + 1).fill(0);
  params[featuresCount] = yAvg;

  if ((targets.length < featuresCount) || (yStdev === 0))
    return params;

  const nonConstFeatureColsIndeces: number[] = [];
  const nonConstFeatureCols: DG.Column[] = [];
  const nonConstFeatureAvgs = new Float32Array(featuresCount);
  const nonConstFeatureStdevs = new Float32Array(featuresCount);

  let idx = 0;
  let nonConstFeaturesCount = 0;
  for (const col of features) {
    const stats = col.stats;

    if (stats.stdev > 0) {
      nonConstFeatureColsIndeces.push(idx);
      nonConstFeatureCols.push(col);
      nonConstFeatureAvgs[nonConstFeaturesCount] = stats.avg;
      nonConstFeatureStdevs[nonConstFeaturesCount] = stats.stdev;
      ++nonConstFeaturesCount;
    }

    ++idx;
  }

  if (nonConstFeaturesCount === 0)
    return params;

  const tempParams = _fitLinearRegressionParamsWithDataNormalizing(
    DG.DataFrame.fromColumns(nonConstFeatureCols).columns,
    DG.Column.fromFloat32Array('xAvgs', nonConstFeatureAvgs, nonConstFeaturesCount),
    DG.Column.fromFloat32Array('xStdevs', nonConstFeatureStdevs, nonConstFeaturesCount),
    targets,
    yAvg,
    yStdev,
    nonConstFeaturesCount + 1,
  ).getRawData();

  for (let i = 0; i < nonConstFeaturesCount; ++i)
    params[nonConstFeatureColsIndeces[i]] = tempParams[i];

  params[featuresCount] = tempParams[nonConstFeaturesCount];

  return params;
} // computeLinRegressionCoefs

/** */
export function getPredictionByLinearRegression(features: DG.ColumnList, params: Float32Array): DG.Column {
  const featuresCount = features.length;
  if (featuresCount !== params.length - 1)
    throw new Error('Incorrect parameters count');

  const col = features.byIndex(0);
  const samplesCount = col.length;
  const prediction = new Float32Array(samplesCount);

  let rawData = col.getRawData();
  const bias = params[featuresCount];
  let weight = params[0];

  for (let i = 0; i < samplesCount; ++i)
    prediction[i] = bias + weight * rawData[i];

  for (let j = 1; j < featuresCount; ++j) {
    rawData = features.byIndex(j).getRawData();
    weight = params[j];

    for (let i = 0; i < samplesCount; ++i)
      prediction[i] += weight * rawData[i];
  }

  return DG.Column.fromFloat32Array('prediction', prediction, samplesCount);
}
