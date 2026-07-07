// JS/TS glue for the XGBoost wasm module (wasm/XGBoostAPI.js + .wasm).
//
// Design (see xgboost-update-plan.md, phases 2-3):
//  - DG-free layer: callers pass typed-array views + dimensions; all
//    DG.Column access happens in src/xgbooster.ts. Column views MUST already
//    be sliced to the column length (col.getRawData() may be longer than
//    col.length - capacity mechanism).
//  - Training runs in a long-lived web worker (transferable buffers, promise
//    queue); prediction is synchronous on the main thread over a cached
//    model handle.
//  - The wasm build has C++ exceptions disabled: a library abort kills the
//    instance. Inputs are validated upstream; on a crash the worker (or the
//    main-thread module) is re-created and model handles become invalid.
//
// The module is imported by webpack in both bundles (classic MODULARIZE +
// UMD tail); the .wasm binary is emitted to dist/ via webpack config and
// located relative to the running bundle, like sci_comp_ml_bg.wasm.

// The .js extension is REQUIRED: webpack's resolve.extensions puts .wasm
// first, so an extensionless './XGBoostAPI' would resolve to XGBoostAPI.wasm
// (a file-loader URL string), not the loader script.
import XGBoostFactory, {XgbWasmModule} from './XGBoostAPI.js';

/** Objective codes shared with the C wrapper (xgboost-api.cpp). */
export enum XgbObjective {
  Regression = 0, // reg:squarederror -> predicted value
  Binary = 1, // binary:logistic  -> probability of class 1
  Multiclass = 2, // multi:softmax    -> class index
}

export interface XgbHyperParams {
  iterations: number;
  eta: number;
  maxDepth: number;
  lambda: number;
  alpha: number;
}

/** Numeric column view: already sliced to the column length. */
export type ColumnView = Float32Array | Float64Array | Int32Array;

const FLOAT_BYTES = 4;

// ---------------------------------------------------------------- module --

const bundleUrl: string =
  (typeof document !== 'undefined' && document.currentScript) ?
    (document.currentScript as HTMLScriptElement).src :
    self.location.href;

/** Absolute URL of the wasm binary in `dist/`, resolved from the bundle URL. */
export function xgbWasmUrl(): string {
  return new URL('XGBoostAPI.wasm', bundleUrl).href;
}

let module_: XgbWasmModule | null = null;
let initPromise: Promise<void> | undefined;

/** Initialize the main-thread module (idempotent). Called from package init. */
export async function initXgboost(): Promise<void> {
  initPromise ??= XGBoostFactory({locateFile: () => xgbWasmUrl()})
    .then((m) => {module_ = m;});
  return initPromise;
}

/** Drop the crashed module instance; the next initXgboost() re-instantiates. */
function resetModule(): void {
  module_ = null;
  initPromise = undefined;
  void initXgboost(); // fire-and-forget re-init so the next predict succeeds
}

function getModule(): XgbWasmModule {
  if (module_ === null)
    throw new Error('XGBoost wasm module is not initialized (call initXgboost first)');
  return module_;
}

// ----------------------------------------------------------------- arena --

/** Growing wasm-heap buffer reused across calls (geometric growth, no shrink). */
class WasmArena {
  private ptr = 0;
  private capacity = 0; // in bytes

  ensureBytes(m: XgbWasmModule, bytes: number): number {
    if (bytes > this.capacity) {
      if (this.ptr !== 0) m._free(this.ptr);
      this.capacity = Math.ceil(bytes * 1.5 / 16) * 16;
      this.ptr = m._malloc(this.capacity);
    }
    return this.ptr;
  }

  invalidate(): void {
    this.ptr = 0;
    this.capacity = 0;
  }
}

// Main-thread arenas: features, prediction output, model bytes.
const featuresArena = new WasmArena();
const outputArena = new WasmArena();
const bytesArena = new WasmArena();

/** Copy column views into the wasm heap col-major; returns the base pointer.
 *  Views of any numeric type are converted by TypedArray.set(). */
function writeColumns(m: XgbWasmModule, arena: WasmArena, cols: ColumnView[], nRows: number): number {
  const ptr = arena.ensureBytes(m, cols.length * nRows * FLOAT_BYTES);
  // Re-read the heap view after ensureBytes: malloc may have grown memory.
  const heap = m.HEAPF32;
  for (let j = 0; j < cols.length; ++j)
    heap.set(cols[j], ptr / FLOAT_BYTES + j * nRows);
  return ptr;
}

// ---------------------------------------------------------------- models --

/** Load a model from serialized bytes; returns a live handle (> 0). */
export function loadXgbModel(modelBytes: Uint8Array): number {
  const m = getModule();
  const ptr = bytesArena.ensureBytes(m, modelBytes.length);
  m.HEAPU8.set(modelBytes, ptr);
  const handle = m._xgbLoadModel(ptr, modelBytes.length);
  if (handle <= 0)
    throw new Error('XGBoost: failed to load model from bytes');
  return handle;
}

export function freeXgbModel(handle: number): void {
  if (module_ !== null && handle > 0)
    module_._xgbFreeModel(handle);
}

/** Synchronous prediction over a live handle. Returns nRows floats
 *  (value / probability of class 1 / class index, by objective). */
export function predictXgb(handle: number, cols: ColumnView[], nRows: number,
  missing: number): Float32Array {
  const m = getModule();
  try {
    const dataPtr = writeColumns(m, featuresArena, cols, nRows);
    const outPtr = outputArena.ensureBytes(m, nRows * FLOAT_BYTES);
    const rc = m._xgbPredict(handle, dataPtr, nRows, cols.length, missing, outPtr, nRows);
    if (rc !== 0)
      throw new Error('XGBoost: prediction failed (invalid model handle or data shape)');
    // Re-read the view: _xgbPredict may have grown memory.
    return m.HEAPF32.slice(outPtr / FLOAT_BYTES, outPtr / FLOAT_BYTES + nRows);
  } catch (err) {
    if (err instanceof WebAssembly.RuntimeError) {
      // The instance is dead (exceptions are disabled in this build).
      [featuresArena, outputArena, bytesArena].forEach((a) => a.invalidate());
      resetModule();
      throw new Error('XGBoost module crashed and was re-initialized; ' +
        'cached model handles are invalid - retry the prediction');
    }
    throw err;
  }
}

// ---------------------------------------------------------------- worker --

interface FitRequest {
  data: Float32Array; // col-major features (transferable)
  labels: Float32Array; // (transferable)
  nRows: number;
  nCols: number;
  missing: number;
  objective: XgbObjective;
  numClass: number;
  hyper: XgbHyperParams;
}

let worker: Worker | null = null;
let workerQueue: Promise<unknown> = Promise.resolve();

function getWorker(): Worker {
  worker ??= new Worker(new URL('./workers/xgboostWorker', import.meta.url));
  return worker;
}

function killWorker(): void {
  worker?.terminate();
  worker = null;
}

function runFitJob(req: FitRequest): Promise<Uint8Array> {
  return new Promise<Uint8Array>((resolve, reject) => {
    const w = getWorker();
    const fail = (message: string) => {
      // The worker's wasm instance may be dead (no exceptions) - recreate.
      killWorker();
      reject(new Error(message));
    };
    w.onmessage = (e) => {
      if (e.data.error !== undefined)
        fail(`XGBoost training failed: ${e.data.error}`);
      else
        resolve(new Uint8Array(e.data.model));
    };
    w.onerror = (e) => fail(`XGBoost training worker crashed: ${e.message ?? 'unknown error'}`);
    w.postMessage({wasmUrl: xgbWasmUrl(), ...req}, [req.data.buffer, req.labels.buffer]);
  });
}

/** Train in the long-lived worker; jobs are queued (one wasm instance there).
 *  Returns the serialized model (UBJSON). */
export function fitXgb(cols: ColumnView[], labels: ColumnView, nRows: number,
  missing: number, objective: XgbObjective, numClass: number,
  hyper: XgbHyperParams): Promise<Uint8Array> {
  // One contiguous col-major copy: it is transferred (zero-copy) to the
  // worker, so DG's own buffers are never detached. TypedArray.set converts
  // int/double views to float32 on the fly.
  const data = new Float32Array(cols.length * nRows);
  for (let j = 0; j < cols.length; ++j)
    data.set(cols[j], j * nRows);
  const labelsF32 = new Float32Array(labels.length);
  labelsF32.set(labels);

  const job = workerQueue.then(() => runFitJob({
    data, labels: labelsF32, nRows, nCols: cols.length, missing, objective, numClass, hyper,
  }));
  workerQueue = job.catch(() => undefined); // keep the queue alive after failures
  return job;
}
