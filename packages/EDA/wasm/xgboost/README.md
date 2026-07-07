# XGBoost wasm module: build assets

Everything needed to (re)build `wasm/XGBoostAPI.js` + `wasm/XGBoostAPI.wasm`
lives in this folder: the C wrapper, the source patch, the build script and
the tests. The XGBoost sources themselves are NOT stored in the repo - they
are fetched from the public upstream (see Prerequisites); build intermediates
go to the system TEMP directory and never enter the package.

Upstream: [dmlc/xgboost](https://github.com/dmlc/xgboost), tag **v3.3.0**,
commit `d5cd2b40725d55747447f66e4a24f9a2c341b0bf`, plus the minimal-build patch.

Key benchmark-driven decisions (2026-07-07, see CHANGELOG 1.7.0): flags
`-O3 -msimd128 -flto` - fastest fit AND predict; SIMD without LTO regresses
fit on wide data by up to 45%; `-Os` makes predict ~5x slower; disabling C++
exceptions saves 531 KB (validation moved to TypeScript); vs the pre-1.7.0
module: training 1.3-2.2x faster, prediction 1.4-3.4x faster.

## Prerequisites (one-time)

Clone the two public dependencies anywhere you like:

```powershell
# 1. XGBoost sources (the build checks out tag v3.3.0 into a temporary
#    git worktree - your clone is never modified)
git clone --recursive https://github.com/dmlc/xgboost <xgboost-dir>

# 2. Emscripten SDK 6.x
git clone https://github.com/emscripten-core/emsdk <emsdk-dir>
cd <emsdk-dir>
./emsdk install latest
./emsdk activate latest
```

`node` must be in PATH (it runs the api-test gate). No other tools are needed:
cmake and ninja come with emsdk.

## Build

From this folder:

```powershell
powershell -File build.ps1 -XgboostDir <xgboost-dir> -EmsdkDir <emsdk-dir>
```

(The script also has defaults matching the standard Datagrok developer layout,
so the parameters may be omitted if they resolve for you.)

The script: creates a temporary worktree at the pinned tag, applies
`patches/*.patch`, builds with emcmake in TEMP, runs `smoke/api-test.mjs`
against the fresh module and - ONLY if it passes - copies `XGBoostAPI.js` +
`XGBoostAPI.wasm` into the package's `wasm/`. Then:

```powershell
npx webpack        # re-bundles the loader and re-emits the .wasm to dist/
grok test --host local --no-retry --category XGBoost
```

## Contents

| What | Why |
|---|---|
| `xgboost-api.cpp` | C wrapper (6 functions: `xgbTrain`, `xgbModelSize`, `xgbModelCopy`, `xgbLoadModel`, `xgbPredict`, `xgbFreeModel`): columnar input via the array interface (no transposition), in-place predict without DMatrix, a table of long-lived model handles. C++ exceptions are disabled in the wasm build - all input validation MUST happen on the TypeScript side (`src/xgbooster.ts`). |
| `CMakeLists.txt` | Wrapper project for emcmake (`-DXGBOOST_SOURCE_DIR=<clone>`); production flags `-O3 -msimd128 -flto` (chosen by benchmark; use SIMD and LTO only together), `-sGROWABLE_ARRAYBUFFERS=0` is mandatory (Emscripten 6.x otherwise breaks TextDecoder in new Chrome). Two link targets: production `web,worker` and test-only `XGBoostAPI-node` (`web,worker,node`). |
| `build.ps1` | Full build cycle with the api-test gate, as described above. |
| `patches/xgboost-wasm-minimal-v3.3.0.patch` | Library trimming: the CMake flag `XGBOOST_WASM_MINIMAL` excludes ranking/survival/quantile/hinge objectives, AUC/rank/survival metrics, `exact`/`approx` tree methods and `gblinear`; plus a `bst_idx_t` to `size_t` narrowing fix for wasm32 (`quantile.cc`). What remains: `gbtree` + `hist`, reg/binary/multi objectives, rmse/logloss/mlogloss metrics. |
| `smoke/smoke.cpp` | Smoke test via the C API: 3 objectives + 6 negative trimming checks. The negative checks need C++ exceptions - compile the wasm variant with `-DSMOKE_SKIP_NEGATIVE`; the full 9/9 run is done on a native build (see `smoke/CMakeLists.txt`). |
| `smoke/api-test.mjs` | Node test of the wasm wrapper exports (23 checks, including a byte-exact model roundtrip). Built into `build.ps1`: if it fails, artifacts are NOT copied. |
| `smoke/bench.mjs` | Wasm-level benchmark matrix: shapes, objectives, depth/iterations, small batches, model load, missing values; `--legacy` mode for the pre-1.7.0 module. |

## Updating the XGBoost version

1. Fetch the new tag in your clone; `build.ps1 -Tag vX.Y.Z` will fail on the
   patch - regenerate it first: branch off the new tag, `git apply` the current
   patch, fix conflicts by the rule: undefined symbol
   `__dmlc_registry_file_tag_X` means drop tag X from the registry; an
   `Unknown objective/metric/updater` error in smoke means return the source
   file to the build. Then `git diff <tag>` into a new file under `patches/`.
2. Run the native smoke (9/9) and `build.ps1` (api-test is the gate).
3. Update the version/commit in this README, in `CREDITS.md` and the package
   CHANGELOG.

## Testing rule

Every module rebuild must pass `api-test.mjs` (built into build.ps1); when the
patch or the XGBoost version changes, additionally run the full native
`smoke.cpp`. Platform tests: `grok test --host local --no-retry --category
XGBoost`.
