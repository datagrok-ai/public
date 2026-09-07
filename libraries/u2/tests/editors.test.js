/* The editor registry: a semType rule supplies the input, override beats it, unmatched props
   fall back to the generated editor, and unregister restores the default. Structural — the
   molecule seed lives in the platform-bound module and is not loaded here. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, TextInput} from '../src/index.js';
import {Editors} from '../src/dg/forms/editors.js';
import {propertyForm} from '../src/dg/forms/object-form.js';

function smoke(name, body) {
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

const structureProp = {name: 'smiles', caption: 'Structure', type: 'string', semType: 'Test',
  get: (o) => o.smiles, set: (o, v) => o.smiles = v};
const nameProp = {name: 'name', type: 'string', get: (o) => o.name, set: (o, v) => o.name = v};

function structureRule() {
  return {
    match: (p) => p.semType === 'Test',
    create: (p, options) => {
      const input = new TextInput(options);
      input.root.dataset.test = 'registered';
      return input;
    },
  };
}

smoke('a matching rule supplies the editor, wired to the property both ways', async () => {
  const unregister = Editors.register(structureRule());
  try {
    const target = {name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O'};
    const form = propertyForm([nameProp, structureProp], target);
    const input = form.input('smiles');
    assert.equal(input.root.dataset.test, 'registered');
    assert.equal(input.value.value, target.smiles, 'seeded from the property');

    input.value.value = 'c1ccccc1';
    await flush();
    assert.equal(target.smiles, 'c1ccccc1', 'edit writes through');
    assert.equal(form.input('name').root.dataset.test, undefined, 'unmatched prop stays generated');
    form.dispose();
  } finally {
    unregister();
  }
});

smoke('an override input beats the registry', async () => {
  const unregister = Editors.register(structureRule());
  const custom = new TextInput({label: 'Custom'});
  try {
    const form = propertyForm([structureProp], {smiles: 'C'}, {overrides: {smiles: {input: custom}}});
    assert.equal(form.input('smiles'), custom);
    form.dispose();
  } finally {
    unregister();
    custom.dispose();
  }
});

smoke('unregister restores the generated editor', async () => {
  Editors.register(structureRule())();
  const form = propertyForm([structureProp], {smiles: 'C'});
  assert.equal(form.input('smiles').root.dataset.test, undefined);
  form.dispose();
});
