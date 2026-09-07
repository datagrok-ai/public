// Hand-written types for the generated Emscripten module wasm/SVMAPI.js
// (classic MODULARIZE + UMD tail; the factory is the CommonJS export).
// Exported functions are defined in wasm/libsvm/svm-api.cpp.

export interface SvmWasmModule {
  _malloc(bytes: number): number;
  _free(ptr: number): void;

  /** Returns a model handle (> 0) or 0 on failure.
   *  colMajor: nRows*nCols float32 (column j at [j*nRows, (j+1)*nRows)).
   *  labels: nRows float32 (classification: category codes; regression: y).
   *  svmType/kernelType: enum codes from svm.h. */
  _svmTrain(dataPtr: number, nRows: number, nCols: number, labelsPtr: number,
    svmType: number, kernelType: number, C: number, gamma: number, degree: number,
    coef0: number, epsilon: number, nu: number, eps: number,
    cacheSize: number, shrinking: number, probability: number): number;
  /** Serializes the model into an internal cache; returns byte count or -1. */
  _svmModelSize(handle: number): number;
  /** Copies the serialized model into dst; returns byte count or -1. */
  _svmModelCopy(handle: number, dstPtr: number, dstSize: number): number;
  /** Returns a model handle (> 0) or 0 on failure. */
  _svmLoadModel(bytesPtr: number, size: number): number;
  /** Returns 0 on success, -1 on failure; out receives nRows floats. */
  _svmPredict(handle: number, dataPtr: number, nRows: number, nCols: number,
    outPtr: number, outLen: number): number;
  _svmFreeModel(handle: number): void;

  /** Heap views: re-read after any call that may grow memory. */
  HEAPF32: Float32Array;
  HEAPU8: Uint8Array;
}

export interface SvmModuleOverrides {
  locateFile?: (file: string, scriptDirectory: string) => string;
  onAbort?: (reason: unknown) => void;
}

declare function SVM(overrides?: SvmModuleOverrides): Promise<SvmWasmModule>;
export default SVM;
