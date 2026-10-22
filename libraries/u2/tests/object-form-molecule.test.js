/* A readonly Molecule field is drawn by the platform's molecule renderer, loaded on first use and
   kept across a refresh while the value stands; every other readonly field is its text. The
   `grok.chem` the renderer draws with is tests/chem-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Access} from '../src/core/access.js';
import {EMPTY} from '../src/core/text.js';
import {propertyForm} from '../src/dg/forms/object-form.js';

register('./chem-stub.mjs', import.meta.url);

function form(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function prop(name, propertyType, options = {}) {
  return {name, propertyType, get: (x) => x[name], set: (x, v) => x[name] = v, ...options};
}

const READONLY = {access: Access.from({can: {}, fields: {structure: 'readonly', name: 'readonly'}})};
const value = (f, name) => f.root.querySelector(`[data-u2-name="${name}"] [data-u2-part="readonly-value"]`);

/** The depiction, once the renderer's first import has landed — a module load, not a microtask. */
async function depiction(f) {
  for (let i = 0; i < 100 && value(f, 'structure').querySelector('.chem-stub-depiction') === null; i++)
    await flush();
  return value(f, 'structure').querySelector('.chem-stub-depiction');
}

form('a readonly molecule is drawn, and the drawing stays through a refresh while the value does', async () => {
  const row = {structure: 'CCO', name: 'Ethanol'};
  const f = propertyForm([prop('structure', 'string', {semType: 'Molecule'}), prop('name', 'string')], row, READONLY);
  document.body.append(f.root);
  assert.equal(value(f, 'structure').textContent, 'CCO', 'the string until the drawing lands');
  assert.equal(value(f, 'name').textContent, 'Ethanol', 'any other field is its text');
  row.structure = 'CCN';
  f.refresh();
  assert.equal(value(f, 'structure').textContent, 'CCN');
  const drawn = await depiction(f);
  assert.notEqual(drawn, null, 'the platform\'s depiction');
  assert.equal(drawn.dataset.molecule, 'CCN', 'the value it stands for now, never one it no longer does');

  f.refresh();
  await flush();
  assert.equal(value(f, 'structure').querySelector('.chem-stub-depiction'), drawn, 'not redrawn for the same value');

  row.structure = 'CCC';
  f.refresh();
  await flush();
  assert.equal(value(f, 'structure').querySelector('.chem-stub-depiction').dataset.molecule, 'CCC');

  row.structure = null;
  f.refresh();
  assert.equal(value(f, 'structure').textContent, '', 'an empty molecule is an empty field');
  assert.equal(value(f, 'structure').querySelector('.chem-stub-depiction'), null);
  f.dispose();
});

form('a required structure is a required editor; text the sketcher drew nothing for stays, flagged unless it is SMILES', async () => {
  const grok = await import('datagrok-api/grok');
  const ui = await import('datagrok-api/ui');
  await import('../src/dg/inputs/molecule.js');
  const row = {structure: ''};
  const f = propertyForm([prop('structure', 'string', {semType: 'Molecule', nullable: false})], row);
  document.body.append(f.root);
  const input = f.input('structure');
  assert.equal(input.nullable, false, 'the property says it must hold a value');
  assert.equal(input.required, true);
  const canvas = f.root.querySelector('.d4-input-molecule-canvas-host canvas');
  const sketch = (enter) => {
    fire(canvas, 'click');
    const sketcher = grok.chem.Sketcher.made.at(-1);
    enter(sketcher);
    ui.shown.at(-1).ok();
  };
  sketch((sketcher) => sketcher.molInput.value = 'C1CC(');
  assert.equal(row.structure, 'C1CC(', 'what was typed is not dropped');
  assert.equal(input.validity.value, 'Cannot read the structure');
  assert.equal(input.root.querySelector('[data-u2-part="error"]').textContent, 'Cannot read the structure',
    'and the field says so');
  sketch((sketcher) => {
    assert.equal(sketcher.molInput.value, 'C1CC(', 'reopened, the text is there to fix');
    sketcher.molInput.value = 'CCN';
  });
  assert.equal(row.structure, 'CCN', 'a SMILES typed without Enter is kept');
  assert.equal(input.validity.value, null);
  sketch((sketcher) => sketcher.setMolecule('CCO'));
  assert.equal(row.structure, 'CCO', 'a drawing is what the sketcher answers');
  assert.equal(input.validity.value, null);
  sketch((sketcher) => sketcher.setMolecule(''));
  assert.equal(input.validity.value, EMPTY, 'cleared, it is required again');
  assert.equal(f.root.querySelector('.d4-input-molecule').platform.platformOpened, 0,
    'the platform\'s own dialog never opened');
  f.dispose();
});
