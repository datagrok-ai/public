# SVMAPI wasm module (LIBSVM)

Builds `wasm/SVMAPI.js` + `wasm/SVMAPI.wasm` — the Emscripten wasm module that
wraps [LIBSVM](https://github.com/cjlin1/libsvm) for the EDA package. Glue:
`svm-api.cpp` (this folder). TypeScript side: `wasm/svm.ts` (DG-free) and
`src/svm.ts` (DG-facing).

## Provenance

- Upstream: `cjlin1/libsvm`, tag **v3.37** (`LIBSVM_VERSION 337`), BSD-3-Clause.
- The source is **not** vendored here. `build.ps1` compiles `svm.cpp` straight
  from an external clone (`-LibsvmDir`); only the built artifacts and this
  wrapper live in the EDA package. The clone is never modified.

## Prerequisites

- A clone of `cjlin1/libsvm` at tag `v3.37` (any location).
- An `emsdk` install (any recent Emscripten; `emsdk_env.ps1` on PATH root).
- `node` on PATH (for the mandatory smoke gate).

## Rebuild

```powershell
powershell -File build.ps1 -LibsvmDir <libsvm-clone> -EmsdkDir <emsdk>
```

The script:
1. Sources the Emscripten environment.
2. Compiles two link targets from `<clone>/svm.cpp` + `svm-api.cpp`:
   - `SVMAPI` (`-sENVIRONMENT=web,worker`) — the shipped artifact.
   - `SVMAPI-node` (`...,node`) — used only by the smoke test, never copied.
3. Runs `smoke/api-test.mjs` against `SVMAPI-node.js`. **Artifacts are copied
   into `wasm/` only if it passes.**

Then re-bundle the loader and re-emit the `.wasm` to `dist/`, and run the suite:

```
npx webpack
grok test --host local --no-retry --category SVM
```

## emcc flags — why

| Flag | Reason |
|---|---|
| `-sMODULARIZE=1 -sEXPORT_NAME=SVM` | Classic MODULARIZE + UMD tail → the same file is importable by webpack in both the main bundle and the worker, and `require()`-able under node for the smoke test. |
| `-sENVIRONMENT=web,worker` (shipped) | Excludes node: a node-enabled loader emits `require("node:fs")`, which breaks webpack bundling. The `web,worker,node` variant is built separately for the smoke test only. |
| `-sGROWABLE_ARRAYBUFFERS=0` | Mandatory. Emscripten 6.x otherwise backs the heap with a resizable ArrayBuffer, and `TextDecoder` throws per spec over such a buffer on the first string LIBSVM would decode. Node lacks the resizable path, so tests pass while real browsers crash. |
| `-sFORCE_FILESYSTEM=1` | LIBSVM (de)serializes only through files; the wrapper round-trips the model through Emscripten MEMFS (`/svm.model`). |
| `-sALLOW_MEMORY_GROWTH=1` | Training builds an `nRows*(nCols+1)` `svm_node` array. |
| (no `-fexceptions`) | LIBSVM's core never throws and never calls `exit()`. Validation is done in TypeScript and by `svm_check_parameter`. |

## Wrapper design notes

Two LIBSVM-specific points handled in `svm-api.cpp`:

1. **No in-memory serialization.** `svm_save_model` / `svm_load_model` are
   file-based → round-trip through MEMFS.
2. **`svm_train` does not copy the support vectors** (`model->SV` point into the
   input buffers). Right after training the wrapper saves+loads through MEMFS to
   get a self-contained model (`free_sv=1`) and frees the input buffers. The
   handle table therefore always holds self-contained models.
