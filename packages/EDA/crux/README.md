# Crux modules

Crux WebAssembly modules used by `src/crux/` (the Crux bridge). `rspack.config.js` copies this folder (without this
README and `engine.toml`) into `dist/crux/`, which `grok publish` uploads; this source folder is not published itself
(`.npmignore`).

| Path | What | In git |
|---|---|---|
| `crux-module-worker.js` | Byte-for-byte copy of crux-core `crates/crux-wasm/js/crux-module-worker.mjs` | yes |
| `eda/crux_wasm.js`, `eda/crux_wasm_bg.wasm` | `module-eda`, wasm-bindgen `--target web` | yes |
| `eda/ENGINE_REV` | The first 12 hex digits of the crux-core commit the module was built from | yes |
| `engine.toml` | The pin: what `module-eda` was built from and how | yes |

The module files are committed (plain git, no LFS), so every checkout publishes a working package. They change only
when the pin changes.

The worker has the `.js` extension, not `.mjs`: Datlas serves package files by the `mime` package's extension
map, which has no `mjs`, and a module worker whose script arrives without a JavaScript MIME type does not start.

## The pin (`engine.toml`)

| Field | Meaning |
|---|---|
| `branch`, `engine_rev` | The crux-core branch and full commit hash the module was built from |
| `crate`, `features` | The crate and its feature set (`--no-default-features --features ...`) |
| `profile`, `rustflags`, `target_dir` | The cargo profile, the wasm32 rustflags and the cargo target directory |
| `wasm_bindgen`, `toolchain` | The wasm-bindgen CLI and the Rust toolchain versions |

## Building `module-eda`

Source: crux-core, branch `eda/module-eda`, commit `d23ccae34`. Toolchain: Rust 1.95.0 (rustup installs it and the
`wasm32-unknown-unknown` target on the first call), wasm-bindgen CLI 0.2.129
(`cargo install wasm-bindgen-cli --version 0.2.129 --locked`), `node`. crux-core's
`tools/parity/wasm-modules/modules.py` does not run on Windows (it imports `fcntl`), so the commands are run by hand
from Git Bash, from the crux-core root:

```bash
CARGO_TARGET_DIR=target/wasm-eda CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS="-C target-feature=+simd128" cargo build --profile wasm-ship -p crux-wasm --target wasm32-unknown-unknown --no-default-features --features "crux-fn/ml,crux-fn/ml-model,crux-fn/xgboost,crux-fn/catboost,crux-fn/conformal,crux-fn/hpo,crux-fn/shap,crux-fn/imblearn"
wasm-bindgen --target nodejs --out-name crux_wasm --out-dir target/wasm-eda/pkg-node target/wasm-eda/wasm32-unknown-unknown/wasm-ship/crux_wasm.wasm
wasm-bindgen --target web --out-name crux_wasm --out-dir target/wasm-eda/pkg-web target/wasm-eda/wasm32-unknown-unknown/wasm-ship/crux_wasm.wasm
```

The feature set has no `crux-fn/lightgbm`, so the module registers 49 functions, without `ml.lightgbm_train` /
`ml.lightgbm_predict`, and `ml.regress` / `ml.classify` reject `model: "lightgbm"`. The `lightgbm` crate is still
compiled in as a dependency of the `shap` crate (its tree explainer).

The nodejs target is only for node probes. Sizes at `d23ccae34`: raw `crux_wasm.wasm` 16,406,303 B,
`crux_wasm_bg.wasm` 15,038,039 B (5,047,201 B gzipped).

Then copy `target/wasm-eda/pkg-web/{crux_wasm.js,crux_wasm_bg.wasm}` (not the `.d.ts` files) into `eda/`, write
`eda/ENGINE_REV` = the 12-hex commit, update `engine.toml`, and run `grok build --force` in EDA (Turborepo does not
see a swapped module and would replay the old `dist/`). If the worker changed, copy
`crates/crux-wasm/js/crux-module-worker.mjs` to `crux-module-worker.js`.
