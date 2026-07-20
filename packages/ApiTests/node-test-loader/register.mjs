// Stubs browser asset imports so the Node test runner can load the same browser-oriented
// js-api / test TypeScript sources webpack bundles, unchanged. Loaded via
// `node --import tsx --import ./node-test-loader/register.mjs` (after tsx, so the CJS
// extension overrides below win over tsx's transformer).
import {register} from 'node:module';
import Module from 'node:module';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc.js';
import advancedFormat from 'dayjs/plugin/advancedFormat.js';

// ESM side: intercept `import './x.css'` from ESM-graph modules.
register('./hooks.mjs', import.meta.url);

// Browser env shim: the test/js-api sources assume dayjs plugins are already loaded.
dayjs.extend(utc);
dayjs.extend(advancedFormat);

// CJS side: intercept `require('./x.css')` from modules tsx treats as CommonJS.
const ASSET_EXTENSIONS = [
  '.css', '.scss', '.sass', '.less',
  '.svg', '.png', '.jpg', '.jpeg', '.gif', '.webp',
  '.woff', '.woff2', '.ttf', '.eot',
];
const stub = (module) => { module.exports = {}; };
for (const ext of ASSET_EXTENSIONS)
  Module._extensions[ext] = stub;

// Single-instance js-api: `startDatagrok` (webpack CommonJS bundle) defines the runtime
// classes and sets globalThis.DG/grok/ui. Without this, test files would load a second
// tsc-compiled copy of the same classes via 'datagrok-api/dg' — structurally identical
// but failing every `instanceof` against runtime-created objects. Alias the subpath
// imports to the runtime globals instead (tests run after startDatagrok, so the lazy
// property reads resolve).
const RUNTIME_GLOBALS = {
  'datagrok-api/dg': () => globalThis.DG,
  'datagrok-api/grok': () => globalThis.grok,
  'datagrok-api/ui': () => globalThis.ui,
};
const origResolve = Module._resolveFilename;
Module._resolveFilename = function(request, ...args) {
  if (request in RUNTIME_GLOBALS) return `\0dg-runtime:${request}`;
  return origResolve.call(this, request, ...args);
};
const origLoad = Module._load;
Module._load = function(request, ...args) {
  if (request in RUNTIME_GLOBALS) {
    const value = RUNTIME_GLOBALS[request]();
    if (value === undefined)
      throw new Error(`${request} requested before startDatagrok() initialized the runtime`);
    return value;
  }
  return origLoad.call(this, request, ...args);
};
