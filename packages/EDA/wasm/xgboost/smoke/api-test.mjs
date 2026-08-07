// Node test for the wasm wrapper API (xgboost-api.cpp exports).
// Usage: node api-test.mjs <path-to-XGBoostAPI-node.js>
//
// Drives the module exactly like the TypeScript glue will: column-major
// float32 data written into the wasm heap, train -> save -> load -> predict.
import {createRequire} from 'node:module';
import {resolve} from 'node:path';

const modulePath = process.argv[2];
if (!modulePath) {
  console.error('usage: node api-test.mjs <path to XGBoostAPI.js>');
  process.exit(2);
}
// The module is a classic MODULARIZE build with a UMD tail: under node it
// is a CommonJS export of the factory.
const require = createRequire(import.meta.url);
const XGBoost = require(resolve(modulePath));
const m = await XGBoost();

const ROWS = 120, COLS = 4, MISSING = -2147483648;
const OBJ = {regression: 0, binary: 1, multiclass: 2};

// Same deterministic generator as smoke.cpp.
function makeData(nClass) {
  const x = new Float32Array(ROWS * COLS);  // column-major
  const y = new Float32Array(ROWS);
  let s = 42 >>> 0;
  for (let i = 0; i < ROWS; ++i) {
    let sum = 0;
    for (let j = 0; j < COLS; ++j) {
      s = (Math.imul(s, 1103515245) + 12345) >>> 0;
      const v = ((s >>> 16) & 0x7fff) / 32768;
      x[j * ROWS + i] = v;  // column j, row i
      sum += v;
    }
    y[i] = nClass <= 0 ? sum : Math.floor(sum * nClass / COLS) % nClass;
  }
  return {x, y};
}

const FLOAT_BYTES = 4;
function alloc(len) { return m._malloc(len * FLOAT_BYTES); }
function writeF32(ptr, arr) { m.HEAPF32.set(arr, ptr / FLOAT_BYTES); }
function readF32(ptr, len) { return m.HEAPF32.slice(ptr / FLOAT_BYTES, ptr / FLOAT_BYTES + len); }

let failures = 0;
function check(cond, label) {
  console.log(`${cond ? 'OK  ' : 'FAIL'} ${label}`);
  if (!cond) ++failures;
}

function runObjective(name, objective, numClass) {
  const {x, y} = makeData(objective === OBJ.regression ? 0 : numClass);
  const xBuf = alloc(x.length), yBuf = alloc(y.length), pBuf = alloc(ROWS);
  writeF32(xBuf, x);
  writeF32(yBuf, y);

  const h = m._xgbTrain(xBuf, ROWS, COLS, MISSING, yBuf, objective, numClass, 10, 0.3, 3, 1.0, 0.0);
  check(h > 0, `${name}: train -> handle ${h}`);

  // Predict with the live handle.
  let rc = m._xgbPredict(h, xBuf, ROWS, COLS, MISSING, pBuf, ROWS);
  const pred = readF32(pBuf, ROWS);
  const inRange = [...pred].every((p) =>
    Number.isFinite(p) &&
    (objective !== OBJ.binary || (p >= 0 && p <= 1)) &&
    (objective !== OBJ.multiclass || (Number.isInteger(p) && p >= 0 && p < numClass)));
  check(rc === 0 && inRange, `${name}: predict values sane`);

  // Save -> load roundtrip must reproduce predictions exactly.
  const size = m._xgbModelSize(h);
  check(size > 0, `${name}: model size ${size} bytes`);
  const bytesBuf = m._malloc(size);
  check(m._xgbModelCopy(h, bytesBuf, size) === size, `${name}: model copy`);
  const h2 = m._xgbLoadModel(bytesBuf, size);
  check(h2 > 0 && h2 !== h, `${name}: load -> new handle ${h2}`);

  const pBuf2 = alloc(ROWS);
  rc = m._xgbPredict(h2, xBuf, ROWS, COLS, MISSING, pBuf2, ROWS);
  const pred2 = readF32(pBuf2, ROWS);
  check(rc === 0 && pred.every((p, i) => p === pred2[i]), `${name}: roundtrip predictions identical`);

  m._xgbFreeModel(h);
  m._xgbFreeModel(h2);
  check(m._xgbPredict(h, xBuf, ROWS, COLS, MISSING, pBuf, ROWS) === -1, `${name}: freed handle rejected`);

  [xBuf, yBuf, pBuf, pBuf2, bytesBuf].forEach((b) => m._free(b));
}

runObjective('reg:squarederror', OBJ.regression, 0);
runObjective('binary:logistic', OBJ.binary, 2);
runObjective('multi:softmax', OBJ.multiclass, 3);

// Bad handle paths must fail cleanly, not abort.
check(m._xgbModelSize(0) === -1, 'zero handle rejected');
check(m._xgbModelSize(999) === -1, 'unknown handle rejected');

console.log(failures ? `\nAPI TEST FAILED (${failures})` : '\nAPI TEST PASSED');
process.exit(failures ? 1 : 0);
