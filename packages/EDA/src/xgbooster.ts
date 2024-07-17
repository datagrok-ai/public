import {memAllocInt32Arr, memAllocFloat32Arr, memFree} from '../wasm/xgbooster';

export function testXGBoost() {
  const iterations = 20;
  const eta = 0.3;
  const maxDepth = 6;
  const lambda = 1;
  const alpha = 0;

  const ROWS = 5;
  const COLS = 3;

  const featuresMemory = memAllocFloat32Arr(ROWS * COLS);
  const features = new Float32Array(featuresMemory.buf, featuresMemory.off, ROWS * COLS);

  const labelsMemory = memAllocFloat32Arr(ROWS);
  const labels = new Float32Array(labelsMemory.buf, labelsMemory.off, ROWS);

  const predictionsMemory = memAllocFloat32Arr(ROWS);
  const predictions = new Float32Array(predictionsMemory.buf, predictionsMemory.off, ROWS);

  const MISSING_VALUE = -1;

  const modelSizeMemory = memAllocInt32Arr(1);
  const modelSizePtr = new Int32Array(modelSizeMemory.buf, modelSizeMemory.off, 1);

  const reserved = 1000000;

  const modelBytesMemory = memAllocInt32Arr(1);
  const modelBytes = new Int32Array(modelBytesMemory.buf, modelBytesMemory.off, 1);

  memFree(featuresMemory.off);
  memFree(labelsMemory.off);
  memFree(predictionsMemory.off);
  memFree(modelSizeMemory.off);
  memFree(modelBytesMemory.off);

  console.log('Done!');
}
