// Regression tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Linear regression migrated to Rust + WASM (sci-comp-ml); now async.
import {_fitLinearRegressionParamsWithDataNormalizing} from '../wasm/eda-api';
import {getPlsAnalysis} from './pls/pls-tools';

export const TOLERANCE = 1e-7;

// Default PLS components count
const PLS_COMPONENTS_COUNT = 10;

/** Default OLS epochs (OLS_EPOCHS in eda-api); worst case, loss tol may early-stop sooner. */
const GD_EPOCHS = 1000;

/** Applicability and interactivity thresholds for linear regression. */
enum LIN_REG_LIMITS {
  // isApplicable: hard bounds — fit must complete without exhausting browser memory
  MAX_FEATURES = 1000,
  MAX_FEATURES_X_SAMPLES = 1e7, // ≈80 MB Float64 flatX (~160 MB peak); ~10 s synchronous at 1000 GD epochs
  // isInteractive: responsive (~0.5 s) full-batch GD budgets
  INTERACTIVE_MEMORY_BUDGET = 2.5e6, // Float64 flatX elements (≈20 MB TS-side, ~40 MB peak)
  INTERACTIVE_TIME_BUDGET = 2.5e8, // N·M·epochs work units (≈0.5 s at conservative 5e8 units/s)
}

/** Check whether linear regression can be applied to the given data */
export function isLinearRegressionApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  for (const col of features) {
    if (!col.matches('numerical'))
      return false;
  }
  if (!target.matches('numerical'))
    return false;

  const featuresCount = features.length;
  const samplesCount = target.length;

  if (featuresCount > LIN_REG_LIMITS.MAX_FEATURES)
    return false;

  if (featuresCount * samplesCount > LIN_REG_LIMITS.MAX_FEATURES_X_SAMPLES)
    return false;

  return true;
}

/** Check whether linear regression training is fast enough for interactive use */
export function isLinearRegressionInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  const M = features.length;
  const N = target.length;

  // Standardised design matrix N×M (Float64, TS flatX); GD forms no M×M normal matrix. Real peak ~2×.
  const memory = N * M;
  if (memory > LIN_REG_LIMITS.INTERACTIVE_MEMORY_BUDGET)
    return false;

  // Full-batch GD: ~N·M work per epoch × epochs; budget = main-thread freeze duration
  const time = GD_EPOCHS * N * M;
  if (time > LIN_REG_LIMITS.INTERACTIVE_TIME_BUDGET)
    return false;

  return true;
}

/** Compute coefficients of linear regression.*/
export async function getLinearRegressionParams(features: DG.ColumnList, targets: DG.Column,
  l1: number = 0, l2: number = 0,
  lr?: number, epochs?: number, tol?: number): Promise<Float32Array> {
  const featuresCount = features.length;
  const samplesCount = targets.length;

  const yAvg = targets.stats.avg;
  const yStdev = targets.stats.stdev;

  const params = new Float32Array(featuresCount + 1).fill(0);
  params[featuresCount] = yAvg;

  // The trivial case
  if ((yStdev === 0) || (samplesCount === 1))
    return params;

  try {
    // Analyze inputs sizes

    // Non-constant columns data
    const nonConstFeatureColsIndeces: number[] = [];
    const nonConstFeatureCols: DG.Column[] = [];
    const nonConstFeatureAvgs = new Float32Array(featuresCount);
    const nonConstFeatureStdevs = new Float32Array(featuresCount);

    let idx = 0;
    let nonConstFeaturesCount = 0;

    // Extract non-constant columns data
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

    // The trivial case
    if (nonConstFeaturesCount === 0)
      return params;

    // Compute parameters of linear regression
    const tempParams = (await _fitLinearRegressionParamsWithDataNormalizing(
      DG.DataFrame.fromColumns(nonConstFeatureCols).columns,
      DG.Column.fromFloat32Array('xAvgs', nonConstFeatureAvgs, nonConstFeaturesCount),
      DG.Column.fromFloat32Array('xStdevs', nonConstFeatureStdevs, nonConstFeaturesCount),
      targets,
      yAvg,
      yStdev,
      nonConstFeaturesCount + 1,
      l1,
      l2,
      lr,
      epochs,
      tol,
    )).getRawData();

    // Extract params taking into account non-constant columns
    for (let i = 0; i < nonConstFeaturesCount; ++i)
      params[nonConstFeatureColsIndeces[i]] = tempParams[i];

    params[featuresCount] = tempParams[nonConstFeaturesCount];
  } catch (e) {
    // Apply PLS regression if regular linear regression failed
    const paramsByPLS = await getLinearRegressionParamsUsingPLS(
      features,
      targets,
      componentsCount(features.length, targets.length),
    );

    let tmpSum = 0;

    // Compute bias (due to the centering feature of PLS)
    for (let i = 0; i < featuresCount; ++i) {
      params[i] = paramsByPLS[i];
      tmpSum += paramsByPLS[i] * features.byIndex(i).stats.avg;
    }

    params[featuresCount] -= tmpSum;
  }

  return params;
} // computeLinRegressionCoefs

/** Structure of the serialized linear regression model */
enum LIN_REG_MODEL {
  FEATURES = 'features',
  VALUES = 'values',
  BIAS = 'bias',
}

/** Pack the model: a two-column dataframe whose last row carries the bias */
export function packLinearRegressionModel(params: Float32Array, featureNames: string[]): Uint8Array {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(LIN_REG_MODEL.FEATURES, featureNames.concat(LIN_REG_MODEL.BIAS)),
    DG.Column.fromFloat32Array(LIN_REG_MODEL.VALUES, params),
  ]).toByteArray();
}

/** Unpack the model. Feature names come back without the trailing bias item,
 * and are absent for models saved before the names were stored. */
export function unpackLinearRegressionModel(model: Uint8Array): {params: Float32Array, names?: string[]} {
  try {
    const df = DG.DataFrame.fromByteArray(model);
    const featuresCol = df.col(LIN_REG_MODEL.FEATURES);
    const valuesCol = df.col(LIN_REG_MODEL.VALUES);

    if ((featuresCol === null) || (valuesCol === null))
      throw new Error('unexpected columns');

    return {
      params: new Float32Array(valuesCol.getRawData().subarray(0, df.rowCount)),
      names: (featuresCol.toList() as string[]).slice(0, -1),
    };
  } catch (_) {
    // Legacy model: a bare Float32Array buffer of coefficients plus the bias.
    // slice() copies into a fresh 4-aligned buffer — the incoming view may sit
    // at an offset Float32Array cannot address.
    if (model.byteLength % 4 !== 0)
      throw new Error('Failed to load model: unexpected format');

    return {params: new Float32Array(model.slice().buffer)};
  }
}

/** Return the model's features in the training order. */
function orderedFeatures(features: DG.ColumnList, params: Float32Array, names?: string[]): DG.Column[] {
  const expected = params.length - 1;

  if (names === undefined) {
    if (features.length !== expected)
      throw new Error('Incorrect parameters count');

    return features.toList();
  }

  if (names.length !== expected)
    throw new Error('Incorrect parameters count');

  const byName = new Map<string, DG.Column>();
  features.toList().forEach((col) => byName.set(col.name, col));

  return names.map((name) => {
    const col = byName.get(name);

    if (col === undefined)
      throw new Error(`Model feature is missing in the table: "${name}"`);

    return col;
  });
}

/** Return prediction of linear regression model */
export function getPredictionByLinearRegression(features: DG.ColumnList, params: Float32Array,
  featureNames?: string[]): DG.Column {
  const cols = orderedFeatures(features, params, featureNames);
  const featuresCount = cols.length;

  const samplesCount = cols[0].length;
  const prediction = new Float32Array(samplesCount);

  let rawData = cols[0].getRawData();
  const bias = params[featuresCount];
  let weight = params[0];

  for (let i = 0; i < samplesCount; ++i)
    prediction[i] = bias + weight * rawData[i];

  for (let j = 1; j < featuresCount; ++j) {
    rawData = cols[j].getRawData();
    weight = params[j];

    for (let i = 0; i < samplesCount; ++i)
      prediction[i] += weight * rawData[i];
  }

  return DG.Column.fromFloat32Array(
    features.getUnusedName('prediction'),
    prediction,
    samplesCount,
  );
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

/** Reteurn linear regression params using the PLS method */
async function getLinearRegressionParamsUsingPLS(features: DG.ColumnList,
  targets: DG.Column, components: number): Promise<Float32Array> {
  const plsAnalysis = await getPlsAnalysis({
    table: DG.DataFrame.fromColumns([targets]),
    features: features,
    predict: targets,
    components: components,
    names: undefined,
    isQuadratic: false,
  });

  return plsAnalysis.regressionCoefficients.getRawData() as Float32Array;
}

/** Return number of PLS components to be used */
const componentsCount = (featuresCount: number, samplesCount: number) => {
  if (samplesCount <= featuresCount)
    return Math.min(PLS_COMPONENTS_COUNT, samplesCount);

  return Math.min(PLS_COMPONENTS_COUNT, featuresCount);
};
