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

import {WasmPca} from './sci_comp_ml.js';
import {ensureEdaMlInit} from '../src/wasm-loader';

// Convergence controls of the retired C++ NIPALS (PCA.h: TOL, MAX_ITER).
// Reproduced here so the Rust kernel iterates identically.
const PCA_TOL = 1e-6;
const PCA_MAX_ITER = 100;

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
