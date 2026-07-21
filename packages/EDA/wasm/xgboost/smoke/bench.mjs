// Fit/predict benchmark suite for the XGBoost wasm module (node, V8 = Chrome engine).
//
// Usage:
//   node bench.mjs <XGBoostAPI-node.js> <label>            # new wrapper API
//   node bench.mjs <old-XGBoostAPI.js>  <label> --legacy   # pre-1.7.0 module
//
// Scenarios (regression ones also run in --legacy mode for comparison):
//   shapes     - fit+predict across data shapes (interactive to stress)
//   objectives - binary:logistic, multi:softmax K=4 / K=10
//   depth/iter - hyperparameter scaling at 50k x 20
//   batch      - predict latency for 1 / 100 / 1000 rows (tooltip scoring)
//   load       - model deserialization time + serialized size
//   missing    - 20% missing values in features
//
// Output: RESULT <label> <scenario> <metric> <ms|bytes>
import {createRequire} from 'node:module';
import {readFileSync} from 'node:fs';
import {resolve, dirname, join} from 'node:path';

const [modulePath, label, legacyFlag] = process.argv.slice(2);
const legacy = legacyFlag === '--legacy';
if (!modulePath || !label) {
  console.error('usage: node bench.mjs <module.js> <label> [--legacy]');
  process.exit(2);
}

const require = createRequire(import.meta.url);
const factory = require(resolve(modulePath));
const wasmPath = join(dirname(resolve(modulePath)),
  resolve(modulePath).replace(/\.js$/, '.wasm').split(/[\\/]/).pop());
const m = await factory({wasmBinary: readFileSync(wasmPath)});

const ITER = 20, ETA = 0.3, DEPTH = 6, LAMBDA = 1, ALPHA = 0;
const MISSING = -2147483648;
const F = 4;
const LEGACY_RESERVE = 10000000; // ints, as in the old package glue

let s = 42 >>> 0;
function rnd() {
  s = (Math.imul(s, 1103515245) + 12345) >>> 0;
  return ((s >>> 16) & 0x7fff) / 32768;
}

function makeData(nRows, nCols, nClass = 0, missingFrac = 0) {
  const x = new Float32Array(nRows * nCols); // col-major
  const y = new Float32Array(nRows);
  for (let j = 0; j < nCols; ++j) {
    for (let i = 0; i < nRows; ++i) {
      const v = rnd();
      x[j * nRows + i] = (missingFrac > 0 && rnd() < missingFrac) ? MISSING : v;
      if (j < 8) y[i] += v * (j + 1);
    }
  }
  if (nClass > 0) {
    const lim = Math.min(nCols, 8) * (Math.min(nCols, 8) + 1) / 2;
    for (let i = 0; i < nRows; ++i)
      y[i] = Math.min(nClass - 1, Math.floor(y[i] * nClass / lim));
  }
  return {x, y};
}

function timeMs(fn, runs = 1, warmup = 0) {
  for (let k = 0; k < warmup; ++k) fn();
  const t0 = performance.now();
  for (let k = 0; k < runs; ++k) fn();
  return (performance.now() - t0) / runs;
}

function report(scenario, metric, value, digits = 1) {
  console.log(`RESULT ${label} ${scenario} ${metric} ${value.toFixed(digits)}`);
}

// ------------------------------------------------------------- new API ----

function newFitPredict(scenario, nRows, nCols, opts = {}) {
  const {objective = 0, numClass = 0, depth = DEPTH, iters = ITER,
    missingFrac = 0, predictRuns = 10} = opts;
  const {x, y} = makeData(nRows, nCols, objective === 0 ? 0 : numClass, missingFrac);
  const xPtr = m._malloc(x.length * F), yPtr = m._malloc(y.length * F);
  const outPtr = m._malloc(nRows * F);
  m.HEAPF32.set(x, xPtr / F);
  m.HEAPF32.set(y, yPtr / F);

  const t0 = performance.now();
  const h = m._xgbTrain(xPtr, nRows, nCols, MISSING, yPtr, objective, numClass,
    iters, ETA, depth, LAMBDA, ALPHA);
  report(scenario, 'fit', performance.now() - t0, 0);
  if (h <= 0) throw new Error(`train failed: ${scenario}`);

  report(scenario, 'predict', timeMs(() => {
    if (m._xgbPredict(h, xPtr, nRows, nCols, MISSING, outPtr, nRows) !== 0)
      throw new Error(`predict failed: ${scenario}`);
  }, predictRuns, 3));

  [xPtr, yPtr, outPtr].forEach((p) => m._free(p));
  return h;
}

function newBatchAndLoad() {
  // Model trained on 50k x 20; predict tiny batches (tooltip scoring).
  const nRows = 50000, nCols = 20;
  const h = newFitPredict('batch-src-50kx20', nRows, nCols);

  const {x} = makeData(1000, nCols);
  const xPtr = m._malloc(x.length * F), outPtr = m._malloc(1000 * F);
  m.HEAPF32.set(x, xPtr / F);
  for (const [rows, runs] of [[1, 1000], [100, 200], [1000, 50]]) {
    report(`batch-${rows}row`, 'predict', timeMs(() => {
      m._xgbPredict(h, xPtr, rows, nCols, MISSING, outPtr, rows);
    }, runs, 10), 3);
  }

  // Serialization size + load (deserialization) time.
  const size = m._xgbModelSize(h);
  report('model-50kx20', 'bytes', size, 0);
  const bytesPtr = m._malloc(size);
  m._xgbModelCopy(h, bytesPtr, size);
  const handles = [];
  report('model-50kx20', 'load', timeMs(() => {
    handles.push(m._xgbLoadModel(bytesPtr, size));
  }, 20, 2));
  handles.forEach((hh) => m._xgbFreeModel(hh));
  m._xgbFreeModel(h);
  [xPtr, outPtr, bytesPtr].forEach((p) => m._free(p));
}

// ---------------------------------------------------------- legacy API ----

function legacyFitPredict(scenario, nRows, nCols, opts = {}) {
  const {missingFrac = 0, predictRuns = 10} = opts;
  const {x, y} = makeData(nRows, nCols, 0, missingFrac);
  const xPtr = m._malloc(x.length * F), yPtr = m._malloc(y.length * F);
  const outPtr = m._malloc(nRows * F);
  const modelPtr = m._malloc(LEGACY_RESERVE * 4), utilsPtr = m._malloc(4);
  m.HEAPF32.set(x, xPtr / F);
  m.HEAPF32.set(y, yPtr / F);

  const t0 = performance.now();
  m._train(xPtr, nRows, nCols, MISSING, yPtr, nRows,
    ITER, ETA, DEPTH, LAMBDA, ALPHA, utilsPtr, 1, modelPtr, LEGACY_RESERVE);
  report(scenario, 'fit', performance.now() - t0, 0);
  const modelBytes = m.HEAP32[utilsPtr / 4];
  if (modelBytes <= 0) throw new Error(`legacy train failed: ${scenario}`);

  // Shipped legacy behavior: the model is re-parsed on EVERY predict call.
  report(scenario, 'predict', timeMs(() => {
    m._predict(xPtr, nRows, nCols, MISSING, modelPtr, modelBytes, outPtr, nRows);
  }, predictRuns, 3));

  const result = {modelPtr, modelBytes, nCols};
  [xPtr, yPtr, outPtr, utilsPtr].forEach((p) => m._free(p));
  return result;
}

function legacyBatch() {
  const {modelPtr, modelBytes, nCols} = legacyFitPredict('batch-src-50kx20', 50000, 20);
  const {x} = makeData(1000, nCols);
  const xPtr = m._malloc(x.length * F), outPtr = m._malloc(1000 * F);
  m.HEAPF32.set(x, xPtr / F);
  for (const [rows, runs] of [[1, 200], [100, 100], [1000, 50]]) {
    report(`batch-${rows}row`, 'predict', timeMs(() => {
      m._predict(xPtr, rows, nCols, MISSING, modelPtr, modelBytes, outPtr, rows);
    }, runs, 5), 3);
  }
  report('model-50kx20', 'bytes', modelBytes, 0);
  [xPtr, outPtr, modelPtr].forEach((p) => m._free(p));
}

// ------------------------------------------------------------- suite ------

const SHAPES = [[1000, 10], [10000, 10], [10000, 100], [50000, 20],
  [100000, 10], [100000, 50]];

if (legacy) {
  legacyFitPredict('warmup-1kx10', 1000, 10);
  for (const [r, c] of SHAPES)
    legacyFitPredict(`reg-${r / 1000}kx${c}`, r, c, {predictRuns: r >= 100000 ? 5 : 10});
  legacyBatch();
  legacyFitPredict('missing20-50kx20', 50000, 20, {missingFrac: 0.2});
} else {
  newFitPredict('warmup-1kx10', 1000, 10);
  for (const [r, c] of SHAPES)
    newFitPredict(`reg-${r / 1000}kx${c}`, r, c, {predictRuns: r >= 100000 ? 5 : 10});

  newFitPredict('binary-10kx100', 10000, 100, {objective: 1, numClass: 2});
  newFitPredict('binary-50kx20', 50000, 20, {objective: 1, numClass: 2});
  newFitPredict('multi4-10kx100', 10000, 100, {objective: 2, numClass: 4});
  newFitPredict('multi4-50kx20', 50000, 20, {objective: 2, numClass: 4});
  newFitPredict('multi10-50kx20', 50000, 20, {objective: 2, numClass: 10});

  newFitPredict('depth3-50kx20', 50000, 20, {depth: 3});
  newFitPredict('depth9-50kx20', 50000, 20, {depth: 9});
  newFitPredict('iter100-50kx20', 50000, 20, {iters: 100});

  newBatchAndLoad();
  newFitPredict('missing20-50kx20', 50000, 20, {missingFrac: 0.2});
}

report('heap', 'bytes', m.HEAPU8.length, 0);
