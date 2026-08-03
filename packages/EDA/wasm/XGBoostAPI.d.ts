// Hand-written types for the generated Emscripten module wasm/XGBoostAPI.js
// (classic MODULARIZE + UMD tail; the factory is the CommonJS export).
// Exported functions are defined in wasm/xgboost/xgboost-api.cpp.

export interface XgbWasmModule {
  _malloc(bytes: number): number;
  _free(ptr: number): void;

  /** Returns a model handle (> 0) or 0 on failure. */
  _xgbTrain(dataPtr: number, nRows: number, nCols: number, missing: number,
    labelsPtr: number, objective: number, numClass: number,
    iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number): number;
  /** Serializes the model into an internal cache; returns byte count or -1. */
  _xgbModelSize(handle: number): number;
  /** Copies the serialized model into dst; returns byte count or -1. */
  _xgbModelCopy(handle: number, dstPtr: number, dstSize: number): number;
  /** Returns a model handle (> 0) or 0 on failure. */
  _xgbLoadModel(bytesPtr: number, size: number): number;
  /** Returns 0 on success, -1 on failure; out receives nRows floats. */
  _xgbPredict(handle: number, dataPtr: number, nRows: number, nCols: number,
    missing: number, outPtr: number, outLen: number): number;
  _xgbFreeModel(handle: number): void;

  /** Heap views: re-read after any call that may grow memory. */
  HEAPF32: Float32Array;
  HEAPU8: Uint8Array;
}

export interface XgbModuleOverrides {
  locateFile?: (file: string, scriptDirectory: string) => string;
  onAbort?: (reason: unknown) => void;
}

declare function XGBoost(overrides?: XgbModuleOverrides): Promise<XgbWasmModule>;
export default XGBoost;
