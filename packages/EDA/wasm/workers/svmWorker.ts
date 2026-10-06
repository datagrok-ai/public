// Worker that runs LIBSVM training off the UI thread: receives a col-major
// Float32Array + labels (transferables), returns serialized model bytes.

// .js extension required.
import SVMFactory from '../SVMAPI.js';
import type {SvmWasmModule} from '../SVMAPI.js';
import type {SvmHyperParams} from '../svm';

const FLOAT_BYTES = 4;

let module_: SvmWasmModule | undefined;

async function ensureInit(wasmUrl: string): Promise<SvmWasmModule> {
  module_ ??= await SVMFactory({locateFile: () => wasmUrl});
  return module_;
}

onmessage = async (e: MessageEvent) => {
  const msg = e.data;
  try {
    const m = await ensureInit(msg.wasmUrl);

    const data: Float32Array = msg.data;
    const labels: Float32Array = msg.labels;
    const h: SvmHyperParams = msg.hyper;

    const dataPtr = m._malloc(data.length * FLOAT_BYTES);
    const labelsPtr = m._malloc(labels.length * FLOAT_BYTES);
    try {
      // Re-read heap views after each malloc: memory may have grown.
      m.HEAPF32.set(data, dataPtr / FLOAT_BYTES);
      m.HEAPF32.set(labels, labelsPtr / FLOAT_BYTES);

      const handle = m._svmTrain(
        dataPtr, msg.nRows, msg.nCols, labelsPtr,
        h.svmType, h.kernelType, h.C, h.gamma, h.degree, h.coef0,
        h.epsilon, h.nu, h.eps, h.cacheSize, h.shrinking, h.probability);
      if (handle <= 0)
        throw new Error('svmTrain failed');

      try {
        const size = m._svmModelSize(handle);
        if (size <= 0)
          throw new Error('svmModelSize failed');
        const bytesPtr = m._malloc(size);
        try {
          if (m._svmModelCopy(handle, bytesPtr, size) !== size)
            throw new Error('svmModelCopy failed');
          // slice() copies out of the heap so the buffer is transferable.
          const model = m.HEAPU8.slice(bytesPtr, bytesPtr + size);
          postMessage({model: model.buffer}, {transfer: [model.buffer]});
        } finally {
          m._free(bytesPtr);
        }
      } finally {
        m._svmFreeModel(handle);
      }
    } finally {
      m._free(dataPtr);
      m._free(labelsPtr);
    }
  } catch (err) {
    postMessage({error: err instanceof Error ? err.message : String(err)});
  }
};
