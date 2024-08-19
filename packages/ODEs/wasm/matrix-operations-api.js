// Matrix operations by wasm-computations

/** Initialize wasm-module */
export async function initMatrOperApi() {
  await initMatrOper();
}

/** Returns buffer & 2 offsets */
export function memAlloc(size) {
  return {
    buf: matrOperLib.HEAPF64.buffer,
    off1: matrOperLib._malloc(size * 8),
    off2: matrOperLib._malloc(size * 8),
  };
}

/** Free memory */
export function memFree(ptr) {
  matrOperLib._free(ptr);
}

/** Computes inverse matrix: previously allocated wasm-memory is used. */
export function inverseMatrix(source, rowCount, output) {
  if (rowCount === 1)
    output[0] = 1.0 / source[0];
  else if (rowCount === 2) {
    const det = source[0] * source[3] - source[1] * source[2];

    output[0] = source[3] / det;
    output[1] = -source[1] / det;
    output[2] = -source[2] / det;
    output[3] = source[0] / det;
  }
  else     
    matrOperLib._invMatrixD(source.byteOffset, rowCount, output.byteOffset);
}
