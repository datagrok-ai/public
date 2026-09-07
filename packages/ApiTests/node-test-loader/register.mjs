// Stubs browser asset imports so the Node test runner can load the same browser-oriented
// js-api / test TypeScript sources webpack bundles, unchanged. Loaded via
// `node --import tsx --import ./node-test-loader/register.mjs` (after tsx, so the CJS
// extension overrides below win over tsx's transformer).
import {register} from 'node:module';
import Module from 'node:module';
import {mkdirSync, writeFileSync} from 'node:fs';
import {dirname, join} from 'node:path';
import {fileURLToPath} from 'node:url';
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
// imports to the runtime globals instead.
const RUNTIME_GLOBALS = {
  'datagrok-api/dg': 'DG',
  'datagrok-api/grok': 'grok',
  'datagrok-api/ui': 'ui',
};

let stubModules = null;

const bindable = (name) => {
  try {
    new Function(`var ${name};`);
    return true;
  }
  catch (_) {
    return false;
  }
};

/// Generates one real module per runtime global and points the resolver at it. The runner
/// calls this once startDatagrok() has produced the globals, before it imports any test
/// file. It has to be a file on disk rather than a virtual module: a test file may be
/// resolved by our own hook, by tsx's, or by Node's, and only a real path survives all
/// three — returning a synthetic specifier made Node read it as a directory import.
export function bindRuntimeGlobals() {
  const dir = join(dirname(fileURLToPath(import.meta.url)), '.runtime');
  mkdirSync(dir, {recursive: true});
  stubModules = {};
  for (const request of Object.keys(RUNTIME_GLOBALS)) {
    const global = RUNTIME_GLOBALS[request];
    const value = globalThis[global];
    if (value === undefined)
      throw new Error(`${request} bound before startDatagrok() initialized the runtime`);
    // ESM named exports are static, so mirror the keys the runtime object carries now.
    const names = Object.keys(value).filter(bindable);
    const file = join(dir, `${global}.mjs`);
    writeFileSync(file, `const g = globalThis.${global};\n` +
      (names.length ? `export const {${names.join(', ')}} = g;\n` : '') +
      'export default g;\n');
    stubModules[request] = file;
  }
}

const origResolve = Module._resolveFilename;
Module._resolveFilename = function(request, ...args) {
  if (stubModules && request in stubModules) return stubModules[request];
  return origResolve.call(this, request, ...args);
};
const origLoad = Module._load;
Module._load = function(request, ...args) {
  if (request in RUNTIME_GLOBALS) {
    const value = globalThis[RUNTIME_GLOBALS[request]];
    if (value === undefined)
      throw new Error(`${request} requested before startDatagrok() initialized the runtime`);
    return value;
  }
  return origLoad.call(this, request, ...args);
};
