/* ListInput and MapInput: the two collection editors ported from Dart — the quote-aware CSV
   contract plus the expand toggle, and the key/value rows with their duplicate-key validity. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {ListInput} from '../src/components/list-input.js';
import {MapInput} from '../src/components/map-input.js';

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

function mount(component) {
  document.body.append(component.root);
  return component;
}

function type(el, text) {
  el.value = text;
  fire(el, 'input');
}

function field(input) {
  return input.root.querySelector('input.u2-list-text');
}

function area(input) {
  return input.root.querySelector('textarea.u2-list-text');
}

function rows(input) {
  return input.root.querySelectorAll('.u2-map-row');
}

function keys(input) {
  return input.root.querySelectorAll('.u2-map-key');
}

function values(input) {
  return input.root.querySelectorAll('.u2-map-value');
}

function rowButton(input, index, action) {
  return rows(input)[index].querySelectorAll('.u2-map-btn')[action === 'add' ? 0 : 1];
}

smoke('ListInput: the quote-aware split round-trips', async () => {
  assert.deepEqual(ListInput.parse('a,"b,c",d'), ['a', 'b,c', 'd']);
  assert.deepEqual(ListInput.parse(''), []);
  assert.deepEqual(ListInput.parse(' a , b '), ['a', 'b'], 'items are trimmed');
  assert.deepEqual(ListInput.parse('a,\nb'), ['a', 'b'], 'line breaks never split an item');
  assert.equal(ListInput.join(['a', 'b,c', 'd']), 'a,"b,c",d');
  assert.equal(ListInput.join([]), '');

  const input = mount(new ListInput({label: 'Synonyms'}));
  assert.equal(input.root.dataset.u2, 'list-input');
  assert.deepEqual(input.value.value, []);
  assert.equal(field(input).placeholder, 'Comma-separated values');

  type(field(input), 'a,"b,c",d');
  assert.deepEqual(input.value.value, ['a', 'b,c', 'd']);

  input.value.value = ['x', 'y,z'];
  assert.equal(field(input).value, 'x,"y,z"', 'a programmatic write re-quotes the text');
  assert.equal(area(input).value, 'x,"y,z"', 'and reaches the collapsed textarea too');

  type(field(input), '');
  assert.deepEqual(input.value.value, []);
  input.dispose();
});

smoke('ListInput: the toggle swaps field and textarea, keeping the text', async () => {
  const input = mount(new ListInput({label: 'Synonyms', value: ['ASA', 'aspirin']}));
  const toggle = input.root.querySelector('.u2-list-toggle');
  assert.equal(input.expanded, false);
  assert.equal(field(input).hidden, false);
  assert.equal(area(input).hidden, true);
  assert.equal(toggle.getAttribute('aria-label'), 'Expand');

  type(field(input), 'ASA,aspirin,2-acetoxybenzoic acid');
  fire(toggle, 'click');
  assert.equal(input.expanded, true);
  assert.equal(field(input).hidden, true);
  assert.equal(area(input).hidden, false);
  assert.equal(area(input).value, 'ASA,aspirin,2-acetoxybenzoic acid', 'expanding keeps the text');
  assert.equal(toggle.getAttribute('aria-label'), 'Collapse');

  type(area(input), 'ASA,aspirin');
  assert.deepEqual(input.value.value, ['ASA', 'aspirin'], 'the textarea edits the same value');
  fire(toggle, 'click');
  assert.equal(input.expanded, false);
  assert.equal(field(input).value, 'ASA,aspirin', 'collapsing keeps it too');
  input.dispose();
});

smoke('ListInput: two inputs on one signal; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const first = mount(new ListInput({label: 'Tags', value: ['a']}));
  const second = mount(new ListInput({label: 'Also tags', bind: first.value}));
  type(field(first), 'a,b');
  assert.equal(field(second).value, 'a,b');

  second.dispose();
  first.dispose();
  assert.equal(Scope.liveCount, before);

  type(field(first), 'a,b,c');
  assert.deepEqual(first.value.value, ['a', 'b'], 'listeners died with the scope');
});

smoke('MapInput: rows write through, and an empty key stays out of the value', async () => {
  const changes = [];
  const input = mount(new MapInput({label: 'Params', onChanged: (v) => changes.push(v)}));
  assert.equal(input.root.dataset.u2, 'map-input');
  assert.equal(rows(input).length, 1, 'an empty map still offers one row to type into');
  assert.deepEqual(input.value.value, {});

  type(values(input)[0], 'DMSO');
  assert.deepEqual(input.value.value, {}, 'a value without a key is not an entry');
  assert.deepEqual(changes, []);

  type(keys(input)[0], 'solvent');
  assert.deepEqual(input.value.value, {solvent: 'DMSO'});
  assert.equal(changes.length, 1);

  type(keys(input)[0], '');
  assert.deepEqual(input.value.value, {}, 'clearing the key drops the entry again');
  input.dispose();
});

smoke('MapInput: add and remove rows', async () => {
  const input = mount(new MapInput({label: 'Params', value: {solvent: 'DMSO'}}));
  assert.equal(rows(input).length, 1);
  assert.equal(keys(input)[0].value, 'solvent');
  assert.equal(values(input)[0].value, 'DMSO');

  fire(rowButton(input, 0, 'add'), 'click');
  assert.equal(rows(input).length, 2, 'the new row lands below');
  assert.equal(keys(input)[1].value, '');
  assert.deepEqual(input.value.value, {solvent: 'DMSO'}, 'an empty row is not a change');

  type(keys(input)[1], 'temperature');
  type(values(input)[1], '25');
  assert.deepEqual(input.value.value, {solvent: 'DMSO', temperature: '25'});

  fire(rowButton(input, 0, 'remove'), 'click');
  assert.equal(rows(input).length, 1);
  assert.deepEqual(input.value.value, {temperature: '25'});

  fire(rowButton(input, 0, 'remove'), 'click');
  assert.equal(rows(input).length, 1, 'the last row is never removed, only emptied');
  assert.deepEqual(input.value.value, {});
  input.dispose();
});

smoke('MapInput: duplicate keys mark the cells and the input invalid', async () => {
  const input = mount(new MapInput({label: 'Params', value: {a: '1', b: '2'}}));
  assert.equal(input.validity.value, null);

  type(keys(input)[1], 'a');
  assert.equal(input.validity.value, 'Duplicate keys are not allowed');
  assert.equal(keys(input)[0].classList.contains('u2-invalid'), true);
  assert.equal(keys(input)[1].classList.contains('u2-invalid'), true);
  assert.equal(input.root.querySelector('.u2-input-error').textContent,
    'Duplicate keys are not allowed');

  type(keys(input)[1], 'c');
  assert.equal(input.validity.value, null);
  assert.equal(keys(input)[0].classList.contains('u2-invalid'), false);
  assert.deepEqual(input.value.value, {a: '1', c: '2'});
  input.dispose();
});

smoke('MapInput: a programmatic write rebuilds the rows; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const input = mount(new MapInput({label: 'Params'}));

  input.value.value = {solvent: 'DMSO', temperature: '25'};
  assert.equal(rows(input).length, 2);
  assert.deepEqual(keys(input).map((k) => k.value), ['solvent', 'temperature']);
  assert.deepEqual(values(input).map((v) => v.value), ['DMSO', '25']);

  input.value.value = {};
  assert.equal(rows(input).length, 1);
  assert.equal(keys(input)[0].value, '');

  input.dispose();
  assert.equal(Scope.liveCount, before);
  type(keys(input)[0], 'solvent');
  assert.deepEqual(input.value.value, {}, 'listeners died with the scope');
});

smoke('a disabled map input stays disabled through a programmatic write', () => {
  const input = mount(new MapInput({value: {a: '1'}}));
  input.enabled = false;
  input.value.value = {a: '1', b: '2'};
  const cells = input.root.querySelectorAll('.u2-map-key, .u2-map-value, .u2-map-btn');
  assert.ok(cells.length > 0);
  for (let i = 0; i < cells.length; i++)
    assert.equal(cells[i].disabled, true, 'rebuilt rows keep the disabled state');
  input.dispose();
});
