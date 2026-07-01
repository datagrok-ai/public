// Web worker that runs the heavy sci-comp-ml fits (PCA / PLS) off the UI
// thread. The C++/Emscripten PCA and PLS used to run in workers; this
// restores that for the Rust + WASM kernels.
//
// The main thread does all DG.Column access (standardisation, building the
// result columns) and ships plain typed arrays here; the worker only owns
// the wasm `fit`. The wasm URL is passed in (a worker has no
// `document.currentScript`); the module is initialised once per worker.

// @ts-ignore - wasm-pack glue resolved via its sidecar .d.ts
import initEdaMl, {WasmPca, WasmPls} from '../../wasm/sci_comp_ml.js';

let initPromise: Promise<unknown> | undefined;

function ensureInit(wasmUrl: string): Promise<unknown> {
  if (initPromise === undefined)
    initPromise = initEdaMl(wasmUrl);
  return initPromise;
}

onmessage = async (e: MessageEvent) => {
  const msg = e.data;
  try {
    await ensureInit(msg.wasmUrl);

    if (msg.method === 'pca') {
      const pca = new WasmPca(msg.componentsCount, msg.tol, msg.maxIter);
      try {
        pca.fit(msg.flatX, msg.nRows);
        const nComponents = pca.nComponents();
        const scores = pca.scores(); // Float64Array, nComponents x nRows
        postMessage({nComponents, scores});
      } finally {
        pca.free();
      }
    } else if (msg.method === 'pls') {
      const pls = new WasmPls(msg.componentsCount);
      try {
        pls.setFeatureStats(msg.xMeans, msg.xStds, msg.yMean, msg.yStd);
        pls.fit(msg.flatX, msg.nRows, msg.stdY);
        const b = pls.regressionCoefficients(); // raw space, length m
        const t = pls.tScores(); // A x nRows
        const u = pls.uScores(); // A x nRows
        const p = pls.xLoadings(); // A x m
        const q = pls.yLoadings(); // length A
        postMessage({b, t, u, p, q});
      } finally {
        pls.free();
      }
    } else {
      postMessage({error: `Unknown method: ${msg.method}`});
    }
  } catch (err) {
    postMessage({error: err instanceof Error ? err.message : String(err)});
  }
};
