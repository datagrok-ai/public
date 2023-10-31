// Matrix operations by wasm-computations

/** Initialize wasm-module */
export async function initMatrOperApi() {
  await initMatrOper();
}

/** Returns inverse matrix */
export function getInverseMatrix(source, rowCount) {
  if (rowCount === 1)
    return new Float64Array([1.0 / source[0]]);
  
  if (rowCount === 2) {
    const det = source[0] * source[3] - source[1] * source[2];

    return new Float64Array([
      source[3] / det,
      -source[1] / det,
      -source[2] / det,
      source[0] / det
    ]);
  }
  
  const size = source.length;  
  const sourceBuf = matrOperLib._malloc(size * 8);
  const invBuf = matrOperLib._malloc(size * 8);
  matrOperLib.HEAPF64.set(source, sourceBuf >> 3);  
  matrOperLib._invMatrixD(sourceBuf, rowCount, invBuf);  
  const inv = new Float64Array(size);
  
  for (let i = 0; i < size; ++i)
    inv[i] = matrOperLib.HEAPF64[invBuf / 8 + i];

  matrOperLib._free(sourceBuf);
  matrOperLib._free(invBuf);

  return inv;
}

/** Returns inverse matrix */
export function invMatrix(source, rowCount, output) {
  if (rowCount === 1)
    output[0] = 1.0 / source[0];
  else if (rowCount === 2) {
    const det = source[0] * source[3] - source[1] * source[2];

    output[0] = source[3] / det;
    output[1] = -source[1] / det;
    output[2] = -source[2] / det;
    output[3] = source[0] / det;
  }
  else {  
    const size = source.length;  
    const sourceBuf = matrOperLib._malloc(size * 8);
    const invBuf = matrOperLib._malloc(size * 8);
    matrOperLib.HEAPF64.set(source, sourceBuf >> 3);  
    matrOperLib._invMatrixD(sourceBuf, rowCount, invBuf);   
  
    for (let i = 0; i < size; ++i)
      output[i] = matrOperLib.HEAPF64[invBuf / 8 + i];

    matrOperLib._free(sourceBuf);
    matrOperLib._free(invBuf);
  }
}
