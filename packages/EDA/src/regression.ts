// Regression tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_fitLinearRegressionParams, _fitLinearRegressionParamsWithDataNormalizing} from '../wasm/EDAAPI';

/** Compute coefficients of linear regression */
export function getLinearRegressionParams(features: DG.ColumnList, targets: DG.Column): Float32Array {
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

/** Return prediction of linear regression model */
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
} // getPredictionByLinearRegression

/** Generate test dataset */
export function getTestDatasetForLinearRegression(rowCount: number, colCount: number,
  featuresScale: number, featuresBias: number, paramsScale: number, paramsBias: number): DG.DataFrame {
  const df = grok.data.demo.randomWalk(rowCount, colCount + 1);
  const cols = df.columns;
  const noiseCol = cols.byIndex(colCount);
  noiseCol.name = 'y (noisy)';
  const yNoisy = noiseCol.getRawData();
  const y = new Float32Array(rowCount).fill(paramsBias);

  let idx = 0;
  let scale = 0;
  let bias = 0;
  let weight = 0;

  for (const col of cols) {
    col.name = `x${idx}`;
    scale = Math.random() * featuresScale;
    bias = Math.random() * featuresBias;
    const arr = col.getRawData();
    weight = Math.random() * paramsScale;

    for (let j = 0; j < rowCount; ++j) {
      arr[j] = scale * arr[j] + bias;
      y[j] += arr[j] * weight;
    }

    ++idx;

    if (idx === colCount)
      break;
  }

  scale = Math.random() * featuresScale;
  bias = Math.random() * featuresBias;

  for (let j = 0; j < rowCount; ++j)
    yNoisy[j] = scale * yNoisy[j] + y[j];

  cols.add(DG.Column.fromFloat32Array('y', y, rowCount));

  return df;
} // getTestDatasetForLinearRegression
