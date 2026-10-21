# Crux modules

Crux WebAssembly modules used by `src/crux/` (the Crux bridge). `rspack.config.js` copies this folder (without this
README) into `dist/crux/`, which `grok publish` uploads; this source folder is not published itself (`.npmignore`).

| Path | What | In git |
|---|---|---|
| `crux-module-worker.js` | Byte-for-byte copy of crux-core `crates/crux-wasm/js/crux-module-worker.mjs` | yes |
| `ml/crux_wasm.js`, `ml/crux_wasm_bg.wasm` | `module-ml`, wasm-bindgen `--target web` | no (gitignored) |

`ml/` is gitignored: build it locally with the recipe below. A fresh checkout without it builds and publishes,
but the `Crux bridge` tests fail with the loader's "failed to load" error. `module-eda` replaces `ml/` after
Crux stage 0.

The worker has the `.js` extension, not `.mjs`: Datlas serves package files by the `mime` package's extension
map, which has no `mjs`, and a module worker whose script arrives without a JavaScript MIME type does not start.

## Building `module-ml`

Below, `<crux-core>` is the crux-core checkout and `<EDA>` is `public/packages/EDA`.

Source: `<crux-core>`, commit `44203d6a6`. `tools/parity/wasm-modules/modules.py`
does not run on Windows (it imports `fcntl`), so its two commands are run by hand from Git Bash. Toolchain:
`cargo` (rustup installs the pinned `1.95.0` and the `wasm32-unknown-unknown` target on the first call), `node`.

1. wasm-bindgen CLI at the version `Cargo.lock` pins (the CLI must match the crate):

   ```bash
   cargo install wasm-bindgen-cli --version 0.2.129 --locked
   ```

2. Build the module (profile `wasm-ship`: fat LTO, one codegen unit; `+simd128` only for wasm32):

   ```bash
   cd <crux-core> && CARGO_TARGET_DIR=target/wasm-modules CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS="-C target-feature=+simd128" cargo build --profile wasm-ship -p crux-wasm --target wasm32-unknown-unknown --no-default-features --features module-ml
   ```

   Measured on 2026-10-21, clean `target/wasm-modules`: 3 min 14 s, peak
   ~3.7 GB for the rustc/linker processes. Output `target/wasm-modules/wasm32-unknown-unknown/wasm-ship/crux_wasm.wasm`,
   19.7 MB.

3. wasm-bindgen, web target (the nodejs target is only for node probes):

   ```bash
   cd <crux-core> && wasm-bindgen --target web --out-name crux_wasm --out-dir target/crux-modules/web/ml target/wasm-modules/wasm32-unknown-unknown/wasm-ship/crux_wasm.wasm
   cd <crux-core> && wasm-bindgen --target nodejs --out-name crux_wasm --out-dir target/crux-modules/nodejs/ml target/wasm-modules/wasm32-unknown-unknown/wasm-ship/crux_wasm.wasm
   ```

4. Copy into this folder (not the `.d.ts` files):

   ```bash
   cp <crux-core>/target/crux-modules/web/ml/crux_wasm.js <crux-core>/target/crux-modules/web/ml/crux_wasm_bg.wasm <EDA>/crux/ml/
   cp <crux-core>/crates/crux-wasm/js/crux-module-worker.mjs <EDA>/crux/crux-module-worker.js
   ```

   `ml/crux_wasm_bg.wasm` is 18.2 MB after wasm-bindgen.

Then `grok build` in EDA emits the folder into `dist/crux/`.
