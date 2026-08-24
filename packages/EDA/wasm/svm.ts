// JS/TS glue for the SVM wasm module (wasm/SVMAPI.js + .wasm).
//
// Design (mirrors wasm/xgbooster.ts):
//  - DG-free layer: callers pass typed-array views + dimensions; all
//    DG.Column access happens in src/svm.ts. Column views MUST already be
//    sliced to the column length (col.getRawData() may be longer than
//    col.length - capacity mechanism).
//  - Training runs in a long-lived web worker (transferable buffers, promise
//    queue); prediction is synchronous on the main thread over a cached
//    model handle.
//  - The wasm build has C++ exceptions disabled: a library abort kills the
//    instance. Inputs are validated upstream; on a crash the worker (or the
//    main-thread module) is re-created and model handles become invalid.
//
// LIBSVM has no `missing` sentinel (unlike XGBoost); missing values are
// rejected in TypeScript before any wasm call.

// The .js extension is REQUIRED: webpack's resolve.extensions puts .wasm
// first, so an extensionless './SVMAPI' would resolve to SVMAPI.wasm (a
// file-loader URL string), not the loader script. Types come from the sidecar
// SVMAPI.d.ts.
import SVMFactory, {SvmWasmModule} from './SVMAPI.js';

/** svm_type codes, shared with the C wrapper (svm-api.cpp) and svm.h. */
export enum SvmType {
  CSvc = 0,
  NuSvc = 1,
  OneClass = 2,
  EpsilonSvr = 3,
  NuSvr = 4,
}

/** kernel_type codes, shared with svm.h. */
export enum KernelType {
  Linear = 0,
  Poly = 1,
  Rbf = 2,
  Sigmoid = 3,
  Precomputed = 4,
}

export interface SvmHyperParams {
  svmType: SvmType;
  kernelType: KernelType;
  C: number;
  gamma: number; // 0 = 1/features
  degree: number;
  coef0: number;
  epsilon: number; // SVR loss epsilon (param.p)
  nu: number;
  eps: number; // stopping tolerance
  cacheSize: number; // MB
  shrinking: number; // 0/1
  probability: number; // 0/1
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
export function svmWasmUrl(): string {
  return new URL('SVMAPI.wasm', bundleUrl).href;
}

let module_: SvmWasmModule | null = null;
let initPromise: Promise<void> | undefined;

/** Initialize the main-thread module (idempotent). Called from package init. */
export async function initSvm(): Promise<void> {
  initPromise ??= SVMFactory({locateFile: () => svmWasmUrl()})
    .then((m) => {module_ = m;});
  return initPromise;
}

/** Drop the crashed module instance; the next initSvm() re-instantiates. */
function resetModule(): void {
  module_ = null;
  initPromise = undefined;
  void initSvm(); // fire-and-forget re-init so the next predict succeeds
}

function getModule(): SvmWasmModule {
  if (module_ === null)
    throw new Error('SVM wasm module is not initialized (call initSvm first)');
  return module_;
}

// ----------------------------------------------------------------- arena --

/** Growing wasm-heap buffer reused across calls (geometric growth, no shrink). */
class WasmArena {
  private ptr = 0;
  private capacity = 0; // in bytes

  ensureBytes(m: SvmWasmModule, bytes: number): number {
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
function writeColumns(m: SvmWasmModule, arena: WasmArena, cols: ColumnView[], nRows: number): number {
  const ptr = arena.ensureBytes(m, cols.length * nRows * FLOAT_BYTES);
  // Re-read the heap view after ensureBytes: malloc may have grown memory.
  const heap = m.HEAPF32;
  for (let j = 0; j < cols.length; ++j)
    heap.set(cols[j], ptr / FLOAT_BYTES + j * nRows);
  return ptr;
}

// ---------------------------------------------------------------- models --

/** Load a model from serialized bytes; returns a live handle (> 0). */
export function loadSvmModel(modelBytes: Uint8Array): number {
  const m = getModule();
  const ptr = bytesArena.ensureBytes(m, modelBytes.length);
  m.HEAPU8.set(modelBytes, ptr);
  const handle = m._svmLoadModel(ptr, modelBytes.length);
  if (handle <= 0)
    throw new Error('SVM: failed to load model from bytes');
  return handle;
}

export function freeSvmModel(handle: number): void {
  if (module_ !== null && handle > 0)
    module_._svmFreeModel(handle);
}

/** Synchronous prediction over a live handle. Returns nRows floats
 *  (classification: the class label as a category code; regression: value). */
export function predictSvm(handle: number, cols: ColumnView[], nRows: number): Float32Array {
  const m = getModule();
  try {
    const dataPtr = writeColumns(m, featuresArena, cols, nRows);
    const outPtr = outputArena.ensureBytes(m, nRows * FLOAT_BYTES);
    const rc = m._svmPredict(handle, dataPtr, nRows, cols.length, outPtr, nRows);
    if (rc !== 0)
      throw new Error('SVM: prediction failed (invalid model handle or data shape)');
    // Re-read the view: _svmPredict may have grown memory.
    return m.HEAPF32.slice(outPtr / FLOAT_BYTES, outPtr / FLOAT_BYTES + nRows);
  } catch (err) {
    if (err instanceof WebAssembly.RuntimeError) {
      // The instance is dead (exceptions are disabled in this build).
      [featuresArena, outputArena, bytesArena].forEach((a) => a.invalidate());
      resetModule();
      throw new Error('SVM module crashed and was re-initialized; ' +
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
  hyper: SvmHyperParams;
}

let worker: Worker | null = null;
let workerQueue: Promise<unknown> = Promise.resolve();

function getWorker(): Worker {
  worker ??= new Worker(new URL('./workers/svmWorker', import.meta.url));
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
        fail(`SVM training failed: ${e.data.error}`);
      else
        resolve(new Uint8Array(e.data.model));
    };
    w.onerror = (e) => fail(`SVM training worker crashed: ${e.message ?? 'unknown error'}`);
    w.postMessage({wasmUrl: svmWasmUrl(), ...req}, [req.data.buffer, req.labels.buffer]);
  });
}

/** Train in the long-lived worker; jobs are queued (one wasm instance there).
 *  Returns the serialized model (LIBSVM text format). */
export function fitSvm(cols: ColumnView[], labels: ColumnView, nRows: number,
  hyper: SvmHyperParams): Promise<Uint8Array> {
  // One contiguous col-major copy: it is transferred (zero-copy) to the
  // worker, so DG's own buffers are never detached. TypedArray.set converts
  // int/double views to float32 on the fly.
  const data = new Float32Array(cols.length * nRows);
  for (let j = 0; j < cols.length; ++j)
    data.set(cols[j], j * nRows);
  const labelsF32 = new Float32Array(labels.length);
  labelsF32.set(labels);

  const job = workerQueue.then(() => runFitJob({
    data, labels: labelsF32, nRows, nCols: cols.length, hyper,
  }));
  workerQueue = job.catch(() => undefined); // keep the queue alive after failures
  return job;
}
