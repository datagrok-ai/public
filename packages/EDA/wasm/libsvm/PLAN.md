# SVM (libsvm) integration plan for the EDA package

**Ticket:** GROK-20754 · **Branch:** `vmakarichev/grok-20754-eda-svm`
**Base:** LIBSVM v3.37 (`cjlin1/libsvm`, BSD-3-Clause)
**Integration model:** mirrors the **XGBoost** stack (C++/Emscripten + worker + packed container); feature standardization is borrowed from **Softmax**.

---

## Fixed decisions

| # | Decision | Choice |
|---|---|---|
| A | libsvm source location | **External clone** `ml-cpp/libsvm` (like XGBoost), `build.ps1 -LibsvmDir`; svm.cpp is NOT vendored into the package |
| B | v1 scope | **Two tasks**: classification (`C_SVC`) + regression (`EPSILON_SVR`), selected by target type |
| C | Hyperparameter tooltip style | **Style B** (positive "Used by… / … only") |
| D | Probability estimates, nu-variants, one-class | Out of v1 (possible extensions) |
| E | Kernels | Linear / Polynomial / RBF / Sigmoid; default RBF |

---

## Point 0. Architecture

```
[UI: ML | Models | Train Model / Apply Model]
        │
        ▼
src/svm.ts            ← DG-facing SVM class: DG.Column ⇄ typed arrays,
   (class SVM)           feature standardization, model packing (container v2)
        │
        ▼
wasm/svm.ts           ← DG-free glue: wasm-heap arenas, worker queue,
   (glue)                synchronous predict over a cached handle
        │
   ┌────┴─────┐
   ▼          ▼
svmWorker.ts   SVMAPI.js/.wasm   ← training in a worker (SMO is slow),
 (training)     (C++ core)          predict synchronous on the main thread
```

**Separation of concerns (as in XGBoost):**
- `src/svm.ts` — the only place that touches `DG.Column`; also does standardization (variant B — the kernel does NOT normalize).
- `wasm/svm.ts` — typed arrays + dimensions only, no DG.
- Training runs in a long-lived worker (transferable buffers, job queue). Prediction is synchronous on the main thread over a cached model handle.

**Two differences from XGBoost that drive the wrapper design:**

1. **No in-memory serialization.** libsvm only has file-based `svm_save_model(filename)` / `svm_load_model(filename)`. → Save/load via **Emscripten MEMFS** (`/svm.model`), reusing libsvm's native text format with no upstream edits.
2. **`svm_train` does NOT copy the support vectors.** `model->SV[j]` point into the input `svm_node` buffers. → Right after training do a **round-trip save+load via MEMFS** → a self-contained model (`free_sv=1`), after which the input buffers are freed. Invariant: the handle table always holds self-contained models.

**Two tasks (selected by target column type, like XGBoost):**

| Task | Target type | Internal libsvm type |
|---|---|---|
| Classification (binary + multiclass one-vs-one) | `string` / `bool` | `C_SVC` |
| Regression | `numerical` (int/float) | `EPSILON_SVR` |

`bool` → 2-class string (`'true'`/`'false'`) with BOOL reconstruction at predict time (like softmax/xgboost).

---

## Point 1. File layout (Variant A — external clone)

**New — C++/build:**
```
wasm/libsvm/
  svm-api.cpp          # thin extern "C" wrapper over libsvm
  build.ps1            # emcc build from ml-cpp/libsvm/svm.cpp + smoke gate
  README.md            # self-contained rebuild instructions
  smoke/api-test.mjs   # mandatory gate: node train→save→load→predict
```
**New — committed build artifacts:**
```
wasm/SVMAPI.js
wasm/SVMAPI.wasm
wasm/SVMAPI.d.ts       # hand-written module types (like XGBoostAPI.d.ts)
```
**New — TS:**
```
wasm/svm.ts                 # DG-free glue          (like wasm/xgbooster.ts)
wasm/workers/svmWorker.ts   # training worker       (like xgboostWorker.ts)
src/svm.ts                  # DG-facing SVM class   (like src/xgbooster.ts)
```
**Edits:**
```
src/package.ts                          # +4 functions + initSvm() in init
webpack.config.js                       # +'./wasm/SVMAPI.wasm' in entry.package[]
src/tests/classifiers-tests.ts          # +category('SVM')
src/tests/model-serialization-tests.ts  # +SVM round-trip
CHANGELOG.md, CREDITS.md, README.md, CLAUDE.md
```
**Not touched:** `src/wasm-loader.ts` (Rust sci-comp-ml only; the SVM glue resolves its own URL via `document.currentScript`); `ml-cpp/libsvm/*` (the upstream clone is read but never modified).

---

## Point 2. C++ wrapper `wasm/libsvm/svm-api.cpp`

Surface 1:1 with `wasm/xgboost/xgboost-api.cpp`.

**2.1 Slot table** (copy of the XGBoost scheme): `struct ModelSlot { svm_model* model; std::string bytes; }`; slot 0 is unused (handle 0 = error); `GetSlot`, `StoreModel`.

**2.2 Input → svm_node**: columnar float32 (column j at `[j*nRows,(j+1)*nRows)`). Internally build dense `svm_node` rows: `nCols` nodes `{index=1..nCols, value=(double)}` + terminator `{index=-1}`; `x_space` = `nRows*(nCols+1)` nodes. Indices are 1-based, dense (zeros are not skipped).

**2.3 `svmTrain(...) -> handle`:**
```cpp
int svmTrain(const float* colMajor, int nRows, int nCols, const float* labels,
             int svmType, int kernelType, double C, double gamma, int degree,
             double coef0, double epsilon, double nu, double eps,
             double cacheSize, int shrinking, int probability);
```
Steps: fill `svm_problem`/`svm_parameter` (`gamma==0` → `1.0/nCols`) → `svm_check_parameter()` (string → free + return 0) → `svm_train()` → **round-trip via MEMFS** → free `x_space`/`prob` → `StoreModel` → handle.

**2.4 Serialization via MEMFS (libsvm-specific):**
```cpp
saveToBytes(model): svm_save_model("/svm.model", m); read the file into std::string
loadFromBytes(b,n): write b to "/svm.model"; return svm_load_model("/svm.model")
```
- `svmModelSize(h)` → `slot.bytes = saveToBytes(...)`, return length.
- `svmModelCopy(h,dst,n)` → memcpy from `slot.bytes`.
- `svmLoadModel(bytes,n)` → `loadFromBytes` → `StoreModel`.
- Fixed path, strictly sequential access (guaranteed by the worker).

**2.5 `svmPredict(h, colMajor, nRows, nCols, out, outLen)`**: per row — build `svm_node`, `out[i] = (float)svm_predict(model, row)` (class label as category code / regression value). `svmFreeModel(h)` → `svm_free_and_destroy_model`.

**2.6 Printing**: in `svmTrain`, call `svm_set_print_string_function(quiet)` — silence the SMO progress output.

---

## Point 3. Build pipeline `wasm/libsvm/build.ps1`

Shorter than the XGBoost one (no CMake/worktree/patch — a single self-contained `svm.cpp`).

**3.1 Params/env:** `-LibsvmDir` (default `ml-cpp/libsvm`), `-EmsdkDir`, `-EdaWasmDir`, `-Tag v3.37`; `& emsdk_env.ps1`.

**3.2 Compile (one command):**
```
emcc -O3 <clone>/svm.cpp svm-api.cpp -I <clone> -o SVMAPI.js \
  -sMODULARIZE=1 -sEXPORT_NAME=SVM -sENVIRONMENT=web,worker \
  -sALLOW_MEMORY_GROWTH=1 -sGROWABLE_ARRAYBUFFERS=0 -sFORCE_FILESYSTEM=1 \
  -sEXPORTED_FUNCTIONS=_svmTrain,_svmModelSize,_svmModelCopy,_svmLoadModel,_svmPredict,_svmFreeModel,_malloc,_free \
  -sEXPORTED_RUNTIME_METHODS=HEAPF32,HEAPU8
```

| Flag | Why |
|---|---|
| `-sMODULARIZE=1 -sEXPORT_NAME=SVM` | UMD tail → importable by webpack in the bundle and worker, require()-able under node |
| `-sENVIRONMENT=web,worker` | WITHOUT node — otherwise `require("node:fs")` breaks bundling |
| `-sGROWABLE_ARRAYBUFFERS=0` | **mandatory** — TextDecoder throws on libsvm strings in new Chrome |
| `-sFORCE_FILESYSTEM=1` | MEMFS for save/load |
| `-sALLOW_MEMORY_GROWTH=1` | large `x_space` |
| no `-fexceptions` | libsvm core is exception-free; validation in TS + `svm_check_parameter` |

SIMD/LTO — not in v1 (no obvious gain; benchmark later if desired).

**3.3 Gate + copy:** node `smoke/api-test.mjs` (train→save→load→predict + negative `C<=0`→0); artifacts are copied into `wasm/` **only on success**. For the smoke test — a separate link `SVMAPI-node.js` (`-sENVIRONMENT=web,worker,node`), never copied into `wasm/`.

**3.4 `smoke/api-test.mjs`:** linearly separable 2D set; round-trip matches; reasonable accuracy; negative check.

**3.5 Rebuild (README):** `build.ps1` → `npx webpack` → `grok test --host local --no-retry --category SVM`.

---

## Point 4. DG-free glue `wasm/svm.ts`

Copy of `wasm/xgbooster.ts`, with SVM signatures/enums.

**4.1 Enums + hyperparameters (C++ contract):**
```ts
enum SvmType    { CSvc=0, NuSvc=1, OneClass=2, EpsilonSvr=3, NuSvr=4 }
enum KernelType { Linear=0, Poly=1, Rbf=2, Sigmoid=3, Precomputed=4 }
interface SvmHyperParams { svmType; kernelType; C; gamma; degree; coef0;
                           epsilon; nu; eps; cacheSize; shrinking; probability }
type ColumnView = Float32Array | Float64Array | Int32Array;
```
**4.2 URL/init:** `bundleUrl` from `document.currentScript`; `svmWasmUrl()`; `initSvm()` idempotent; `resetModule()` after a crash.
**4.3 Arenas:** `WasmArena` + `writeColumns` (re-read the heap view after malloc); features/output/modelBytes arenas.
**4.4 load/predict (synchronous):** `loadSvmModel`, `freeSvmModel`, `predictSvm` (slice `HEAPF32.slice`; on `WebAssembly.RuntimeError` — invalidate arenas + `resetModule`). Difference: **no `missing` parameter** (missing values are rejected in TS).
**4.5 Training in the worker:** long-lived Worker + `workerQueue` (one instance → serialized jobs), transferable col-major `Float32Array` + labels; `fitSvm(cols, labels, nRows, hyper)`.
**4.6 Not carried over:** `numClass`/`objective` — libsvm infers the class count from the labels.

---

## Point 5. DG-facing class `src/svm.ts`

Skeleton — `src/xgbooster.ts`; standardization — from `src/softmax-classifier.ts`.

**5.1 Fields:** `modelBytes?`, `taskType?`, `targetCategories?`, `targetWasBool`, `avgs`/`stdevs` (standardization is part of the model), `handle`.
**5.2 `isApplicable`** (like xgboost): all features numerical; target numerical/string/bool. **`isInteractive`** — thresholds **lower** than xgboost (SMO is superlinear): ~5,000 rows / ~50 features (numbers to be refined).
**5.3 Standardization:** `extractStats` (avg/stdev from `col.stats`) + `normalized` (`(raw−avg)/stdev`, `stdev==0`→0). One path for fit and predict.
**5.4 `fit`:** validate → `initSvm()` → bool→string → pick `svmType` by target type → category codes / numeric y → `extractStats` + `normalized` → `gamma=0`→`1/features` → `fitSvm(...)` → `invalidateHandle`.
**5.5 `predict`:** `normalized` (SAME avgs/stdevs) → `loadSvmModel` (handle cache) → `predictSvm` (crash → invalidate + rethrow) → decode: classification "code→category/bool" (simpler than xgboost — libsvm returns a ready label, no 0.5 threshold), regression — restore the target type (int/bigint/float).
**5.6 `toBytes`/constructor (container v2):** `[headerSize(4)|headerDf|categoriesDf|avgs/stdevs Df|libsvm-bytes]`; header: `Type, ModelSize, Categories size, StatsSize, WasBool, SvmType, KernelType, Version`; missing `Version` → "retrain" error. Standardization stored as float32 columns `Avg-s`/`Stddev-s` (softmax pattern).
**5.7 `dispose`/handle cache:** `invalidateHandle` → `freeSvmModel`; repeated predicts reuse the handle, after dispose — reload from bytes.

> Hyperparameters are NOT duplicated in the header — `svm_type`/`kernel_type`/`gamma`/`degree`/`coef0` already live inside the libsvm bytes; the container stores only what the TS layer needs (task type, categories, standardization, bool flag).

---

## Point 6. Function registration `src/package.ts`

A block of 4 functions with `meta.mlname: 'SVM'` (modeled on XGBoost). The model is auto-listed in the **ML | Models | Train Model / Apply Model** UI; no separate top-menu items are added.

**6.1 `trainSVM` (mlrole: train)** — `outputs: [{name:'model', type:'dynamic'}]`, returns `Uint8Array`. Parameters and **final tooltips (style B):**

| Parameter | type / options | description |
|---|---|---|
| kernel | string, `choices:['Linear','Polynomial','RBF','Sigmoid']`, `initialValue:'RBF'` | `Kernel function: linear, polynomial, RBF, or sigmoid.` |
| cost | double, `min:0.001 max:1000 initialValue:1` | `Regularization strength C. Higher C fits the training data more tightly.` |
| gamma | double, `min:0 max:100 initialValue:0` | `Kernel coefficient. Used by the RBF, polynomial, and sigmoid kernels; 0 = 1/features.` |
| degree | int, `min:1 max:10 initialValue:3` | `Degree of the polynomial kernel. Used by the polynomial kernel only.` |
| coef0 | double, `min:-100 max:100 initialValue:0` | `Independent kernel term. Used by the polynomial and sigmoid kernels only.` |
| epsilon | double, `min:0 max:10 initialValue:0.1` | `Epsilon in the SVR loss. Used for regression targets only.` |

**6.2 `applySVM` (mlrole: apply)** — copy of `applyXGBooster`: `new SVM(model).predict(df.columns)` → DataFrame.
**6.3 `isApplicableSVM` / `isInteractiveSVM`** — delegate to the class's static methods.
**6.4 init:** `await initSvm()` in `PackageFunctions.init()` next to `initXgboost()`.
**6.5** After edits — `grok api` (regenerates `package.g.ts`/`package-api.ts`) + `grok check --soft`.

> UI note: all parameters are always visible; libsvm ignores `degree`/`coef0`/`epsilon` for an unsuitable kernel/task (safe). Style-B tooltips make this clear. Conditional visibility — optional, later.

---

## Point 7. webpack / package.json / init edits

**7.1 `webpack.config.js`:** add `'./wasm/SVMAPI.wasm'` to `entry.package` (next to `XGBoostAPI.wasm`) — file-loader emits the binary into `dist/`. Nothing else.
**7.2 `package.json`:** EDA has no `sources`/`browserFeatures` fields — the mechanism is the same as XGBoost; nothing to add. `scripts` untouched (wasm build is a separate `build.ps1`).
**7.3 init:** `await initSvm()` (see 6.4); the training worker is created lazily on the first `fitSvm`.
**7.4 `src/wasm-loader.ts`** — not touched (Rust sci-comp-ml only).

---

## Point 8. Tests

**8.1 `category('SVM')` in classifiers-tests.ts** (mirrors XGBoost): Correctness (acc>0.9, seeded), Bool target (round-trip, BOOL), Multiclass (one-vs-one, categories ⊆ target), Regression float (madError, round-trip), Int target, Capacity tail (slice to `col.length`), Validation rejects bad inputs (throw before wasm), Performance (benchmark).
**8.2 SVM-specific:** standardization matters (feature ×1000 — accuracy does not degrade); kernel round-trip (RBF and Linear); `gamma=auto` (`0`→`1/nfeatures`).
**8.3 `model-serialization-tests.ts`:** SVM on iris/cars/winequality — pack→unpack identical.
**8.4 Sizes/thresholds:** SMO is superlinear — perf benchmark ~1–2K rows × 10–20 features (NOT 50K like xgboost); document why.
**8.5 Native gate:** `smoke/api-test.mjs` under node — mandatory in `build.ps1`.
**8.6 Run:** `grok test --host local --no-retry --category SVM` + full `classifiers`.

> Determinism: C_SVC/EPSILON_SVR are deterministic (SMO is not randomized) → seeded datasets are reproducible. If probability is enabled later (Platt over 5-fold CV with `rand()`), write probability tests against a range.

---

## Point 9. Documentation

**9.1 `CHANGELOG.md` → `## v.next`:**
```
* GROK-20754: Added SVM (libsvm) support for classification and regression
```
**9.2 `CREDITS.md`** (required by BSD-3):
```
* [LIBSVM](https://github.com/cjlin1/libsvm) v3.37 — C.-C. Chang, C.-J. Lin. BSD-3-Clause.
```
**9.3 `wasm/libsvm/README.md`:** prerequisites (clone `cjlin1/libsvm` tag v3.37, emsdk, node), `build.ps1` command, smoke gate, flags and why, version pin.
**9.4 `packages/EDA/CLAUDE.md`** — an "SVM (libsvm)" subsection, do-not-revert notes: MEMFS serialization; round-trip after `svm_train` (SVs are not copied); standardization in TS; `isInteractive` thresholds lower than XGBoost; link to the README; `category 'SVM'` under Testing.
**9.5 Fix the wording** "XGBoost is the only remaining C++/Emscripten kernel" in `packages/EDA/CLAUDE.md` and `.claude/rules/wasm.md` → **two** modules (XGBoost + libsvm/SVM).
**9.6 `packages/EDA/README.md`** — in the "Supervised machine learning" section:
```
* [Support vector machine](https://en.wikipedia.org/wiki/Support-vector_machine) (SVM) classification and regression
```

---

## Work order (stages)

1. **C++ wrapper + build**: `svm-api.cpp`, `build.ps1`, smoke gate → `SVMAPI.js/.wasm`; verify train/save/load/predict under node.
2. **DG-free layer + worker**: `wasm/svm.ts`, `svmWorker.ts`.
3. **DG class**: `src/svm.ts` (standardization, container, decoding).
4. **Registration + init**: `package.ts`, `webpack`, `grok api && grok check`.
5. **Tests + docs**: `npm run build`, then `grok test --category SVM`.

## Possible extensions (after v1)

- Probability estimates (`svm_predict_probability`, `-b 1`; Platt + pairwise coupling) — behind a flag; account for CV nondeterminism.
- nu-variants (`NU_SVC` / `NU_SVR`) — alternative parameterization.
- One-class SVM (`ONE_CLASS`) — a separate task (novelty detection, target-free training), its own input.
- Conditional hyperparameter visibility in the Train Model form.
- SIMD/LTO in the build based on benchmark results.
