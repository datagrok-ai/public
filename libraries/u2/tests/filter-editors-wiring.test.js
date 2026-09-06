/* Loading the dg filter entry alone wires the hint router: `FilterBuilder.defaultEditors` honours a
   property's `inputType`/`editor` hints without anything importing the object form first. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextArea} from '../src/components/inputs/text-input.js';
import {FilterBuilder} from '../src/components/filter/filter-builder.js';

register('./dg-stub.mjs', import.meta.url);
await import('../src/dg/filter/index.js');

test('the dg filter entry routes editor hints on its own', () => {
  const scope = new Scope();
  try {
    const prop = {name: 'notes', type: 'string', inputType: 'TextArea'};
    const input = Scope.runWith(scope, () => FilterBuilder.defaultEditors(prop, {}));
    assert.ok(input instanceof TextArea);
  } finally {
    scope.dispose();
    resetDom();
  }
});
