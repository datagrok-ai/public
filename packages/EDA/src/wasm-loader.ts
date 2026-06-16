// Lazy initialisation + cache for the sci-comp-ml (Rust + WASM) module.
//
// The module is a wasm-pack `--target web` build copied into `wasm/`
// (`sci_comp_ml.js` glue + `sci_comp_ml_bg.wasm` binary). It needs one
// async `init` before any class is constructed; we cache that promise so
// repeated calls (and concurrent ones) share a single instantiation.
//
// The binary is emitted next to the bundle in `dist/` (via the `.wasm`
// file-loader rule in webpack.config.js). We resolve it relative to the
// served bundle URL rather than `_package.webRoot` - the latter reads the
// package handle, which is not initialised in the separate
// `package-test.js` bundle, so it would throw under `grok test`.
//
// The bundle URL is captured from `document.currentScript` at module-load
// time (it is only valid during the bundle's initial synchronous
// execution) — the same mechanism the Emscripten EDA loader uses for
// EDA.wasm, so it works in both `package.js` and `package-test.js`.

import initEdaMl from '../wasm/sci_comp_ml.js';

const bundleUrl: string =
  (typeof document !== 'undefined' && document.currentScript)
    ? (document.currentScript as HTMLScriptElement).src
    : self.location.href;

let initPromise: Promise<unknown> | undefined;

/** Absolute URL of the wasm binary in `dist/`, resolved from the bundle URL. */
export function edaMlWasmUrl(): string {
  return new URL('sci_comp_ml_bg.wasm', bundleUrl).href;
}

/** Ensure the sci-comp-ml wasm is instantiated; safe to call repeatedly. */
export function ensureEdaMlInit(): Promise<unknown> {
  if (initPromise === undefined)
    initPromise = initEdaMl(edaMlWasmUrl());

  return initPromise;
}
