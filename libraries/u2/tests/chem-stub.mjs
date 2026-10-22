/* Module hooks serving the two `datagrok-api` modules `src/dg/inputs/molecule.ts` imports: a
   `grok.chem` whose `drawMolecule` hands back a host naming the molecule it was asked for, and an
   empty `ui`. Not a test file: `node --test tests/*.test.js` skips it. */

const GROK = `export const chem = {
  drawMolecule(molecule, w, h) {
    const host = document.createElement('div');
    host.className = 'chem-stub-depiction';
    host.dataset.molecule = molecule;
    return host;
  },
};`;

const URLS = {'datagrok-api/grok': 'u2-test:chem-stub-grok', 'datagrok-api/ui': 'u2-test:chem-stub-ui'};
const SOURCES = {'u2-test:chem-stub-grok': GROK, 'u2-test:chem-stub-ui': 'export {};'};

export async function resolve(specifier, context, next) {
  const url = URLS[specifier];
  return url === undefined ? next(specifier, context) : {url, format: 'module', shortCircuit: true};
}

export async function load(url, context, next) {
  const source = SOURCES[url];
  return source === undefined ? next(url, context) : {source, format: 'module', shortCircuit: true};
}
