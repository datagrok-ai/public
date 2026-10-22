/* A readonly Molecule field is drawn by the platform's molecule renderer, loaded on first use and
   kept across a refresh while the value stands; every other readonly field is its text. The
   `grok.chem` the renderer draws with is tests/chem-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Access} from '../src/core/access.js';
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
