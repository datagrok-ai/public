// JavaScript API for call wasm-functions from the XGBoostAPI module

const INT_BYTES = 4;
const FLOAT_BYTES = 4;

export async function initXgboost() {
  await initXGBoostModule();
}

/** Returns buffer & offset: the Int32Array case */
export function memAllocInt32Arr(size) {
  return {
    buf: XGBoostModule.HEAP32.buffer,
    off: XGBoostModule._malloc(size * INT_BYTES),
  };
}

/** Returns buffer & offset: the Float32Array case */
export function memAllocFloat32Arr(size) {
  return {
    buf: XGBoostModule.HEAPF32.buffer,
    off: XGBoostModule._malloc(size * FLOAT_BYTES),
  };
}

/** Free memory */
export function memFree(ptr) {
  XGBoostModule._free(ptr);
}
