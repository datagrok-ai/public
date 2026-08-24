// ESM module-customization hook used by the Node test runner (package-test-node.ts).
// The js-api / test sources are written for the browser and statically import CSS
// (and occasionally other assets), which webpack handles via loaders but Node cannot
// parse. This hook returns an empty module for any asset import so the runner can load
// the same TypeScript sources the browser bundles, unchanged.
const ASSET = /\.(css|scss|sass|less|svg|png|jpe?g|gif|webp|woff2?|ttf|eot)(\?.*)?$/;

// Same single-instance aliasing register.mjs installs for the CommonJS graph, for modules
// the ESM loader owns instead. Those never reach Module._load, so without this they get the
// CJS resolver's `dg-runtime:` sentinel back and Node imports it as a directory.
//
// Hooks run on their own thread and cannot read globalThis.DG, so the runner re-registers
// this module with the runtime's export names once startDatagrok() has produced them; the
// registration register.mjs makes earlier carries no names and only stubs assets.
const RUNTIME_GLOBALS = {
  'datagrok-api/dg': 'DG',
  'datagrok-api/grok': 'grok',
  'datagrok-api/ui': 'ui',
};
const RUNTIME_SCHEME = 'dg-runtime:';

let runtimeExports = null;

const bindable = (name) => {
  try {
    new Function(`var ${name};`);
    return true;
  }
  catch (_) {
    return false;
  }
};

export async function initialize(data) {
  runtimeExports = data?.runtimeExports ?? null;
}

export async function resolve(specifier, context, nextResolve) {
  if (runtimeExports && specifier in RUNTIME_GLOBALS)
    return {url: RUNTIME_SCHEME + specifier, shortCircuit: true};
  return nextResolve(specifier, context);
}

export async function load(url, context, nextLoad) {
  if (ASSET.test(url))
    return {format: 'module', source: 'export default {};', shortCircuit: true};
  if (url.startsWith(RUNTIME_SCHEME)) {
    const global = RUNTIME_GLOBALS[url.slice(RUNTIME_SCHEME.length)];
    const names = (runtimeExports[global] ?? []).filter(bindable);
    const source = `const g = globalThis.${global};\n` +
      (names.length ? `export const {${names.join(', ')}} = g;\n` : '') +
      'export default g;\n';
    return {format: 'module', source, shortCircuit: true};
  }
  return nextLoad(url, context);
}
