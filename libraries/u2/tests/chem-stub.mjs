/* Module hooks serving the two `datagrok-api` modules `src/dg/inputs/molecule.ts` imports: a
   `grok.chem` whose `drawMolecule` hands back a host naming the molecule it was asked for and whose
   `Sketcher` draws a text only once `setMolecule` was given it (a SMILES without balanced brackets
   reads as none), and a `ui` whose `input.molecule` is a platform molecule input with its
   click-to-edit host (`platformOpened` counts the clicks its own handler got) and whose `dialog`
   records each one shown. Not a test file: `node --test tests/*.test.js` skips it. */

const GROK = `export const chem = {
  drawMolecule(molecule, w, h) {
    const host = document.createElement('div');
    host.className = 'chem-stub-depiction';
    host.dataset.molecule = molecule;
    return host;
  },
  isMolBlock: (s) => s.includes('M  END'),
  checkSmiles: (s) => /^[A-Za-z0-9()=#]+$/.test(s) && s.split('(').length === s.split(')').length,
  Sketcher: class {
    static made = [];
    constructor() {
      chem.Sketcher.made.push(this);
      this.root = document.createElement('div');
      this.molInput = document.createElement('input');
      this.molInput.value = '';
      this.drawn = '';
    }
    setMolecule(molecule) {
      this.molInput.value = molecule;
      this.drawn = chem.checkSmiles(molecule) ? molecule : '';
    }
    isEmpty() {
      return this.drawn === '';
    }
    getSmiles() {
      return this.drawn;
    }
  },
};`;

const UI = `export const shown = [];
export function dialog() {
  const d = {content: [], ok: null,
    add(x) { this.content.push(x); return this; },
    onOK(fn) { this.ok = fn; return this; },
    show() { shown.push(this); return this; }};
  return d;
}
export const input = {
  molecule(caption, options) {
    const root = document.createElement('div');
    root.className = 'd4-input-molecule ui-input-root';
    const label = document.createElement('label');
    label.className = 'ui-label ui-input-label';
    label.textContent = caption;
    const host = document.createElement('div');
    host.className = 'd4-input-molecule-canvas-host';
    const canvas = document.createElement('canvas');
    canvas.className = 'ui-input-editor d4-input-molecule-canvas';
    host.append(canvas);
    root.append(label, host);
    const listeners = [];
    let value = options?.value ?? '';
    const dg = {root, caption, platformOpened: 0,
      get value() { return value; },
      set value(v) { value = v; for (const fn of listeners.slice()) fn(v); },
      onChanged: {subscribe(fn) {
        listeners.push(fn);
        return {unsubscribe() { listeners.splice(listeners.indexOf(fn), 1); }};
      }}};
    host.addEventListener('click', () => dg.platformOpened++);
    root.platform = dg;
    return dg;
  },
};`;

const URLS = {'datagrok-api/grok': 'u2-test:chem-stub-grok', 'datagrok-api/ui': 'u2-test:chem-stub-ui'};
const SOURCES = {'u2-test:chem-stub-grok': GROK, 'u2-test:chem-stub-ui': UI};

export async function resolve(specifier, context, next) {
  const url = URLS[specifier];
  return url === undefined ? next(specifier, context) : {url, format: 'module', shortCircuit: true};
}

export async function load(url, context, next) {
  const source = SOURCES[url];
  return source === undefined ? next(url, context) : {source, format: 'module', shortCircuit: true};
}
