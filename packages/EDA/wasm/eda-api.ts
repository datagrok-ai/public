// New JS API for the EDA ML methods, backed by Rust + WebAssembly
// (the `sci-comp-ml` crate) instead of the retired C++/Emscripten module.
//
// This file replaces, export by export, the matching functions of the
// generated `EDAAPI.js`. Each export keeps the exact name, signature and
// return shape of its `EDAAPI` counterpart, so the `src/*.ts` consumers
// only repoint their import (`../wasm/EDAAPI` -> `../wasm/eda-api`); their
// bodies do not change. The old `EDAAPI` stays in place until every method
// is migrated and green, then it is deleted.
//
// "Variant B": the kernels never standardise. Standardisation happens here
// (the adapter), the caller keeps passing raw columns. Columns from
// Datagrok are read with `getRawData()`, whose typed array may be longer
// than the row count — we trim to `col.length` and convert to f64.

import * as DG from 'datagrok-api/dg';

import {WasmElasticNet, WasmPca, WasmPls, WasmSoftmax} from './sci_comp_ml.js';
import {ensureEdaMlInit} from '../src/wasm-loader';

// Convergence controls of the retired C++ NIPALS (PCA.h: TOL, MAX_ITER).
// Reproduced here so the Rust kernel iterates identically.
const PCA_TOL = 1e-6;
const PCA_MAX_ITER = 100;

// OLS-via-gradient-descent controls for the linear-regression path. The
// retired C++ solved the normal equations in closed form; here full-batch
// GD on standardised data (l1 = l2 = 0) approaches the same least-squares
// solution. Data is unit-variance, so a modest rate is stable; epochs are
// capped to keep the large-sample Performance run within budget, and the
// loss `tol` early-stops once converged.
const OLS_LR = 0.1;
const OLS_EPOCHS = 1000;
const OLS_TOL = 1e-7;

/**
 * Principal component analysis via NIPALS.
 *
 * Drop-in for the `EDAAPI` export of the same name: takes the feature
 * columns and the requested component count, returns a `DG.DataFrame` of
 * `componentsCount` score columns named `"1".."k"` (Float32, length
 * `rowCount`) — the shape `computePCA` expects.
 *
 * The C++ NIPALS standardised each column to `(x - mean) / sigma_pop`
 * internally (`X - means`, then `colwise().normalized() * sqrt(n)`); the
 * Rust kernel does not standardise, so we do it here with the population
 * standard deviation, trimming each raw array to the real row count.
 *
 * (The name keeps the `InWebWorker` suffix for contract compatibility;
 * this first cut runs on the main thread. A worker can be reintroduced
 * later without changing the signature.)
 */
export async function _principalComponentAnalysisNipalsInWebWorker(
  table: DG.DataFrame, columns: DG.ColumnList, componentsCount: number): Promise<DG.DataFrame> {
  await ensureEdaMlInit();

  const cols = columns.toList();
  const m = cols.length;
  const nRows = table.rowCount;

  // Standardise (population) into the flat [col0, col1, ...] boundary form.
  const flat = new Float64Array(m * nRows);
  for (let j = 0; j < m; ++j) {
    const raw = cols[j].getRawData(); // capacity may exceed nRows
    let sum = 0;
    for (let i = 0; i < nRows; ++i)
      sum += raw[i];
    const mean = sum / nRows;

    let ss = 0;
    for (let i = 0; i < nRows; ++i) {
      const d = raw[i] - mean;
      ss += d * d;
    }
    let sigma = Math.sqrt(ss / nRows); // population sigma, matching C++
    if (sigma === 0)
      sigma = 1; // constant column guard

    const base = j * nRows;
    for (let i = 0; i < nRows; ++i)
      flat[base + i] = (raw[i] - mean) / sigma;
  }

  const pca = new WasmPca(componentsCount, PCA_TOL, PCA_MAX_ITER);
  try {
    pca.fit(flat, nRows);

    const k = pca.nComponents();
    const scores = pca.scores(); // flat k x nRows, row-major (component by component)

    const outCols: DG.Column[] = [];
    for (let a = 0; a < componentsCount; ++a) {
      const data = new Float32Array(nRows);
      // Components beyond the extracted rank stay zero (C++ returned all
      // componentsCount columns; the exhausted ones are ~0 anyway).
      if (a < k) {
        const base = a * nRows;
        for (let i = 0; i < nRows; ++i)
          data[i] = scores[base + i];
      }
      outCols.push(DG.Column.fromFloat32Array((a + 1).toString(), data));
    }

    return DG.DataFrame.fromColumns(outCols);
  } finally {
    pca.free();
  }
}

/**
 * Partial least squares regression (PLS1).
 *
 * Drop-in for the `EDAAPI` export of the same name. Returns an array
 * indexed by `WASM_OUTPUT_IDX`:
 *   [0] prediction       — `D · b` (raw features, **no intercept**), length rowCount
 *   [1] regr_coeffs      — raw-space coefficients `b`, length = #features
 *   [2] tScores          — X-scores, `componentsCount` columns of length rowCount
 *   [3] uScores          — Y-scores, `componentsCount` columns of length rowCount
 *   [4] xLoadings        — X-loadings, `componentsCount` columns of length #features
 *   [5] yLoadings        — Y-loadings `q`, one column of length componentsCount
 *
 * Matches the retired C++ `partialLeastSquareExtended`: predictors and
 * response are standardised here with the population std (`||x - mean|| /
 * sqrt(n)`), the kernel trains on the standardised data, coefficients are
 * unwound to raw space via the stats, and `prediction = D · b` is computed
 * on the raw predictors without the intercept (the consumer adds the bias
 * via `debiasedPrediction`).
 */
export async function _partialLeastSquareRegressionInWebWorker(
  table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, componentsCount: number,
): Promise<[DG.Column, DG.Column, DG.Column[], DG.Column[], DG.Column[], DG.Column]> {
  await ensureEdaMlInit();

  const cols = features.toList();
  const m = cols.length;
  const nRows = table.rowCount;

  // Raw predictor arrays (kept for the D·b prediction) + population stats,
  // building the standardised flat [col0, col1, ...] boundary form.
  const raw = cols.map((c) => c.getRawData()); // each capacity may exceed nRows
  const xMeans = new Float64Array(m);
  const xStds = new Float64Array(m);
  const flatX = new Float64Array(m * nRows);
  for (let j = 0; j < m; ++j) {
    const r = raw[j];
    let sum = 0;
    for (let i = 0; i < nRows; ++i)
      sum += r[i];
    const mean = sum / nRows;

    let ss = 0;
    for (let i = 0; i < nRows; ++i) {
      const d = r[i] - mean;
      ss += d * d;
    }
    let sigma = Math.sqrt(ss / nRows);
    if (sigma === 0)
      sigma = 1;

    xMeans[j] = mean;
    xStds[j] = sigma;
    const base = j * nRows;
    for (let i = 0; i < nRows; ++i)
      flatX[base + i] = (r[i] - mean) / sigma;
  }

  // Response stats + standardised response.
  const ry = predict.getRawData();
  let ysum = 0;
  for (let i = 0; i < nRows; ++i)
    ysum += ry[i];
  const yMean = ysum / nRows;
  let yss = 0;
  for (let i = 0; i < nRows; ++i) {
    const d = ry[i] - yMean;
    yss += d * d;
  }
  let yStd = Math.sqrt(yss / nRows);
  if (yStd === 0)
    yStd = 1;
  const stdY = new Float64Array(nRows);
  for (let i = 0; i < nRows; ++i)
    stdY[i] = (ry[i] - yMean) / yStd;

  const pls = new WasmPls(componentsCount);
  try {
    pls.setFeatureStats(xMeans, xStds, yMean, yStd);
    pls.fit(flatX, nRows, stdY);

    const b = pls.regressionCoefficients(); // raw space, length m

    // prediction = D · b on raw predictors, no intercept (C++ parity).
    const predData = new Float32Array(nRows);
    for (let i = 0; i < nRows; ++i) {
      let s = 0;
      for (let j = 0; j < m; ++j)
        s += b[j] * raw[j][i];
      predData[i] = s;
    }

    const tFlat = pls.tScores(); // A x nRows
    const uFlat = pls.uScores(); // A x nRows
    const pFlat = pls.xLoadings(); // A x m
    const qArr = pls.yLoadings(); // length A
    const a = qArr.length;

    const tScores: DG.Column[] = [];
    const uScores: DG.Column[] = [];
    for (let k = 0; k < a; ++k) {
      const tcol = new Float32Array(nRows);
      const ucol = new Float32Array(nRows);
      const base = k * nRows;
      for (let i = 0; i < nRows; ++i) {
        tcol[i] = tFlat[base + i];
        ucol[i] = uFlat[base + i];
      }
      tScores.push(DG.Column.fromFloat32Array((k + 1).toString(), tcol));
      uScores.push(DG.Column.fromFloat32Array((k + 1).toString(), ucol));
    }

    const xLoadings: DG.Column[] = [];
    for (let k = 0; k < a; ++k) {
      const pcol = new Float32Array(m);
      const base = k * m;
      for (let j = 0; j < m; ++j)
        pcol[j] = pFlat[base + j];
      xLoadings.push(DG.Column.fromFloat32Array((k + 1).toString(), pcol));
    }

    const prediction = DG.Column.fromFloat32Array('0', predData);
    const regrCoeffs = DG.Column.fromFloat32Array('0', Float32Array.from(b));
    const yLoadings = DG.Column.fromFloat32Array('0', Float32Array.from(qArr));

    return [prediction, regrCoeffs, tScores, uScores, xLoadings, yLoadings];
  } finally {
    pls.free();
  }
}

/**
 * Softmax (multinomial logistic) classifier training.
 *
 * Drop-in for the `EDAAPI` export of the same name. Returns a
 * `DG.DataFrame` whose `.columns` are `classesCount` parameter columns of
 * length `featuresCount + 1` — the `[W | B]` layout (per class: weights
 * then bias) the `SoftmaxClassifier` consumer slices by index.
 *
 * Matches the retired C++ `fitSoftmax`: the predictors are standardised
 * here with the caller-supplied per-feature `featureAvgs`/`featureStdDevs`
 * (`(x - avg) / stdev`, leaving a column merely centred when `stdev <= 0`,
 * as the C++ did); the kernel trains on the standardised columns with the
 * reference class-frequency weighting (the `WasmSoftmax` default). The
 * parameters come back in standardised space — the consumer's `predict`
 * standardises its input with the same stats before applying them.
 *
 * Note: unlike the C++ `_fitSoftmax`, this is async (the wasm needs a
 * one-time init), so the call site awaits it.
 */
export async function _fitSoftmax(
  features: DG.ColumnList, featureAvgs: DG.Column, featureStdDevs: DG.Column, targets: DG.Column,
  classesCount: number, iterCount: number, learningRate: number, penalty: number, tolerance: number,
  paramsRows: number, _paramsCols: number,
): Promise<DG.DataFrame> {
  await ensureEdaMlInit();

  const cols = features.toList();
  const m = cols.length;
  const nRows = targets.length;
  const avgs = featureAvgs.getRawData();
  const stdevs = featureStdDevs.getRawData();

  // Standardise predictors with the supplied stats into the flat boundary form.
  const flatX = new Float64Array(m * nRows);
  for (let j = 0; j < m; ++j) {
    const r = cols[j].getRawData(); // capacity may exceed nRows
    const avg = avgs[j];
    const sd = stdevs[j] > 0 ? stdevs[j] : 1; // C++ leaves column centred when stdev <= 0
    const base = j * nRows;
    for (let i = 0; i < nRows; ++i)
      flatX[base + i] = (r[i] - avg) / sd;
  }

  const labels = (targets.getRawData() as Int32Array).subarray(0, nRows);

  const sm = new WasmSoftmax(classesCount, learningRate, iterCount, penalty, tolerance);
  try {
    sm.fit(flatX, nRows, labels);

    const params = sm.params(); // flat c x (n+1) row-major [W|B]
    const outCols: DG.Column[] = [];
    for (let c = 0; c < classesCount; ++c) {
      const colData = new Float32Array(paramsRows);
      const base = c * paramsRows;
      for (let r = 0; r < paramsRows; ++r)
        colData[r] = params[base + r];
      outCols.push(DG.Column.fromFloat32Array((c + 1).toString(), colData));
    }

    return DG.DataFrame.fromColumns(outCols);
  } finally {
    sm.free();
  }
}

/**
 * Linear-regression coefficients with internal data normalisation.
 *
 * Drop-in for the `EDAAPI` export of the same name. Returns a single
 * `DG.Column` of length `paramsCount` (= #features + 1): the raw-space
 * weights `[w0 … w_{m-1}]` followed by the bias.
 *
 * Matches the retired C++ `fitLinearRegressionParamsWithDataNormalizing`:
 * predictors and response are standardised with the caller-supplied stats
 * (`(x - avg) / stdev`, `(y - yAvg) / yStdev`), OLS is solved on the
 * standardised data, then coefficients are rescaled back to raw space:
 *   w_raw_j = w_std_j · yStdev / stdev_j
 *   bias    = yAvg − Σ w_raw_j · avg_j
 * The standardised intercept is ~0 (centred data) and, as in the C++, is
 * not carried into the raw bias.
 *
 * Note: unlike the C++ `_fit...`, this is async (the wasm needs a one-time
 * init), so the call site awaits it.
 */
export async function _fitLinearRegressionParamsWithDataNormalizing(
  features: DG.ColumnList, featureAvgs: DG.Column, featureStdDevs: DG.Column, targets: DG.Column,
  targetsAvg: number, targetsStdDev: number, paramsCount: number,
): Promise<DG.Column> {
  await ensureEdaMlInit();

  const cols = features.toList();
  const m = cols.length;
  const nRows = targets.length;
  const xAvgs = featureAvgs.getRawData();
  const xStdevs = featureStdDevs.getRawData();

  // Standardise predictors and response with the supplied stats.
  const flatX = new Float64Array(m * nRows);
  for (let j = 0; j < m; ++j) {
    const r = cols[j].getRawData(); // capacity may exceed nRows
    const avg = xAvgs[j];
    const sd = xStdevs[j] > 0 ? xStdevs[j] : 1;
    const base = j * nRows;
    for (let i = 0; i < nRows; ++i)
      flatX[base + i] = (r[i] - avg) / sd;
  }

  const ry = targets.getRawData();
  const ySd = targetsStdDev > 0 ? targetsStdDev : 1;
  const stdY = new Float64Array(nRows);
  for (let i = 0; i < nRows; ++i)
    stdY[i] = (ry[i] - targetsAvg) / ySd;

  // OLS via full-batch GD; no feature stats -> standardised-space coeffs.
  const model = new WasmElasticNet(OLS_LR, OLS_EPOCHS, 0, 0, OLS_TOL);
  try {
    model.fit(flatX, nRows, stdY);
    const cwb = model.coefficientsWithBias(); // [w_std_0 … w_std_{m-1}, b_std]

    const out = new Float32Array(paramsCount);
    let sum = 0;
    for (let j = 0; j < m; ++j) {
      const sd = xStdevs[j] > 0 ? xStdevs[j] : 1;
      const wRaw = cwb[j] * targetsStdDev / sd;
      out[j] = wRaw;
      sum += wRaw * xAvgs[j];
    }
    out[m] = targetsAvg - sum;

    return DG.Column.fromFloat32Array('params', out, paramsCount);
  } finally {
    model.free();
  }
}
