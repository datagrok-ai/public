import {PatinaeModule} from './patinae-types';

// Resolved against the served bundle (dist/package.js) rather than _package.webRoot so that the
// same code works from the package-test.js bundle under `grok test`.
const bundleUrl: string = (typeof document !== 'undefined' && document.currentScript) ?
  (document.currentScript as HTMLScriptElement).src : self.location.href;

let viewerModule: Promise<PatinaeModule> | null = null;

// Imports the vendored Patinae ES module (wasm/) once; the wasm binary is fetched by its glue.
export function getPatinae(): Promise<PatinaeModule> {
  const url = new URL('../wasm/patinae-viewer.js', bundleUrl).href;
  return viewerModule ??= import(/* webpackIgnore: true */ url);
}
