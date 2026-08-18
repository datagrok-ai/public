/* MultiSelect: the popup machine, tag rendering in the closed field, and the value contract —
   every write is a fresh array in items order, and programmatic writes re-render both surfaces. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Form, TextInput} from '../src/index.js';
import {MultiSelect} from '../src/components/multi-select.js';

const ITEMS = ['MW', 'LogP', 'TPSA', 'HBA', {value: 'hbd', label: 'H-bond donors'}];
const MANY = ['MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'Rotatable bonds', 'Rings', 'Formal charge',
  'Heavy atoms', 'Fraction Csp3'];

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

function texts(root, selector) {
  return root.querySelectorAll(selector).map((el) => el.textContent);
}

function select(options = {}) {
  return mount(new MultiSelect({items: ITEMS, placeholder: 'Pick properties…', ...options}));
}

function field(ms) {
  return ms.root.querySelector('.u2-multi-select-field');
}

function popup() {
  return document.body.querySelector('.u2-multi-select-popup');
}

function rows() {
  return popup().querySelectorAll('.u2-multi-select-option');
}

function tags(ms) {
  return texts(ms.root, '.u2-tag-text');
}

smoke('machine: click opens, Esc closes, state and ARIA follow', async () => {
  const ms = select();
  assert.equal(ms.machineState.value, 'idle');

  field(ms).focus();
  assert.equal(ms.machineState.value, 'focused');

  fire(field(ms), 'click');
  await flush();
  assert.equal(ms.isOpen.value, true);
  assert.equal(ms.machineState.value, 'open');
  assert.equal(field(ms).getAttribute('aria-expanded'), 'true');
  assert.equal(popup().querySelector('.u2-multi-select-options').getAttribute('role'), 'listbox');
  assert.deepEqual(texts(popup(), '.u2-multi-select-option .u2-multi-select-text'),
    ['MW', 'LogP', 'TPSA', 'HBA', 'H-bond donors']);

  fire(popup(), 'keydown', {key: 'Escape'});
  assert.equal(ms.isOpen.value, false);
  assert.equal(ms.machineState.value, 'focused');
  assert.equal(popup(), null, 'the overlay is gone');
  assert.equal(field(ms).getAttribute('aria-expanded'), 'false');

  fire(field(ms), 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(ms.isOpen.value, true, 'ArrowDown reopens');
  ms.dispose();
});

smoke('toggling keeps the popup open and writes a fresh array in items order', async () => {
  const ms = select();
  fire(field(ms), 'click');
  await flush();

  const empty = ms.value.value;
  fire(rows()[4], 'pointerdown');
  fire(rows()[1], 'pointerdown');
  assert.equal(ms.isOpen.value, true, 'a multi-select stays open while picking');
  assert.equal(popup().isConnected, true);
  assert.deepEqual(ms.value.value, ['LogP', 'hbd'], 'items order, not click order');
  assert.notEqual(ms.value.value, empty, 'every write is a fresh array');
  assert.deepEqual(rows().map((r) => r.getAttribute('aria-selected')),
    ['false', 'true', 'false', 'false', 'true']);
  assert.deepEqual(tags(ms), ['LogP', 'H-bond donors']);

  const before = ms.value.value;
  fire(rows()[1], 'pointerdown');
  assert.deepEqual(ms.value.value, ['hbd'], 'a second click deselects');
  assert.notEqual(ms.value.value, before);
  assert.deepEqual(before, ['LogP', 'hbd'], 'the previous array is never mutated');
  ms.dispose();
});

smoke('keyboard: Down moves the active row, Space toggles it, Enter closes', async () => {
  const ms = select();
  fire(field(ms), 'keydown', {key: ' '});
  await flush();

  fire(popup(), 'keydown', {key: 'ArrowDown'});
  fire(popup(), 'keydown', {key: 'ArrowDown'});
  fire(popup(), 'keydown', {key: 'ArrowUp'});
  const active = rows()[0];
  assert.equal(active.classList.contains('u2-multi-select-option-active'), true);
  assert.equal(popup().getAttribute('aria-activedescendant'), active.id);

  fire(popup(), 'keydown', {key: ' '});
  assert.deepEqual(ms.value.value, ['MW']);
  assert.equal(ms.isOpen.value, true);

  fire(popup(), 'keydown', {key: 'Enter'});
  assert.equal(ms.isOpen.value, false);
  assert.deepEqual(ms.value.value, ['MW']);
  ms.dispose();
});

smoke('a tag ✕ removes its value without opening the popup', async () => {
  const ms = select({value: ['MW', 'TPSA']});
  const remove = ms.root.querySelectorAll('.u2-tag-remove')[0];
  assert.equal(remove.getAttribute('aria-label'), 'Remove MW');

  fire(remove, 'click');
  await flush();
  assert.deepEqual(ms.value.value, ['TPSA']);
  assert.deepEqual(tags(ms), ['TPSA']);
  assert.equal(ms.isOpen.value, false, 'removing never opens the popup');
  assert.equal(popup(), null);
  ms.dispose();
});

smoke('select all / none toggles the whole list from one row', async () => {
  const ms = select();
  fire(field(ms), 'click');
  await flush();

  const all = popup().querySelector('.u2-multi-select-all');
  assert.equal(all.textContent, 'Select all');
  fire(all, 'pointerdown');
  assert.deepEqual(ms.value.value, ['MW', 'LogP', 'TPSA', 'HBA', 'hbd']);
  assert.equal(all.textContent, 'Select none');
  assert.equal(all.classList.contains('u2-multi-select-selected'), true);

  fire(all, 'pointerdown');
  assert.deepEqual(ms.value.value, []);
  assert.equal(all.textContent, 'Select all');
  ms.dispose();
});

smoke('filter: on by default above eight items, narrows the rows, select-all follows it', async () => {
  const few = select();
  fire(field(few), 'click');
  await flush();
  assert.equal(popup().querySelector('.u2-multi-select-filter-input'), null, 'five items need no filter');
  few.dispose();

  const ms = mount(new MultiSelect({items: MANY}));
  fire(field(ms), 'click');
  await flush();
  const filter = popup().querySelector('.u2-multi-select-filter-input');
  assert.equal(rows().length, 10);

  filter.value = 'on';
  fire(filter, 'input');
  await flush();
  assert.deepEqual(texts(popup(), '.u2-multi-select-option .u2-multi-select-text'),
    ['Rotatable bonds', 'Fraction Csp3']);

  fire(popup().querySelector('.u2-multi-select-all'), 'pointerdown');
  assert.deepEqual(ms.value.value, ['Rotatable bonds', 'Fraction Csp3'],
    'select all acts on the filtered rows');

  filter.value = 'zzz';
  fire(filter, 'input');
  await flush();
  assert.equal(rows().length, 0);
  assert.equal(popup().querySelector('.u2-multi-select-empty').textContent, 'No matches');

  filter.value = '';
  fire(filter, 'input');
  await flush();
  assert.equal(rows().length, 10);
  ms.dispose();
});

smoke('closed field: placeholder, maxTags overflow, programmatic writes re-render both surfaces', async () => {
  const ms = select();
  assert.equal(ms.root.querySelector('.u2-multi-select-placeholder').textContent, 'Pick properties…');

  ms.value.value = ['MW', 'LogP', 'TPSA', 'HBA', 'hbd'];
  assert.deepEqual(tags(ms), ['MW', 'LogP', 'TPSA']);
  assert.equal(ms.root.querySelector('.u2-multi-select-more').textContent, '+2 more');
  assert.equal(ms.root.querySelector('.u2-multi-select-placeholder'), null);

  fire(field(ms), 'click');
  await flush();
  assert.deepEqual(rows().map((r) => r.getAttribute('aria-selected')),
    ['true', 'true', 'true', 'true', 'true'], 'a programmatic write reaches the checkboxes too');

  ms.value.value = [];
  assert.deepEqual(tags(ms), []);
  assert.deepEqual(rows().map((r) => r.getAttribute('aria-selected')),
    ['false', 'false', 'false', 'false', 'false']);
  assert.equal(ms.root.querySelector('.u2-multi-select-placeholder').textContent, 'Pick properties…');

  const wide = select({maxTags: 1, value: ['MW', 'LogP']});
  assert.deepEqual(tags(wide), ['MW']);
  assert.equal(wide.root.querySelector('.u2-multi-select-more').textContent, '+1 more');
  wide.dispose();
  ms.dispose();
});

smoke('setItems prunes vanished values and stays silent when nothing vanished', async () => {
  const changes = [];
  const ms = select({value: ['LogP', 'hbd'], onChanged: (v) => changes.push(v)});
  assert.deepEqual(tags(ms), ['LogP', 'H-bond donors']);

  ms.setItems(['MW', 'TPSA', {value: 'hbd', label: 'H-bond donors'}]);
  assert.deepEqual(ms.value.value, ['hbd'], 'LogP is gone from the items, so it leaves the value');
  assert.deepEqual(tags(ms), ['H-bond donors']);
  assert.equal(changes.length, 1);

  ms.setItems(['TPSA', {value: 'hbd', label: 'H-bond donors'}]);
  assert.deepEqual(ms.value.value, ['hbd']);
  assert.equal(changes.length, 1, 'nothing vanished — no value churn');

  fire(field(ms), 'click');
  await flush();
  assert.deepEqual(texts(popup(), '.u2-multi-select-option .u2-multi-select-text'),
    ['TPSA', 'H-bond donors'], 'the popup shows the new items');
  ms.dispose();
});

smoke('validators and Form validity come from the input base', async () => {
  const ms = select({label: 'Properties'});
  ms.addValidator((v) => v.length < 2 ? 'Pick at least 2' : null);
  assert.equal(ms.validity.value, 'Pick at least 2');
  assert.equal(field(ms).classList.contains('u2-invalid'), true);
  assert.equal(ms.root.querySelector('.u2-input-error').textContent, 'Pick at least 2');

  const form = mount(new Form());
  const name = new TextInput({label: 'Name', value: 'Series A'});
  form.addAll([name, ms]);
  await flush();
  assert.equal(form.validity.value, 'Pick at least 2');
  assert.deepEqual(form.getValues(), {Name: 'Series A', Properties: []});

  ms.value.value = ['MW', 'TPSA'];
  assert.equal(ms.validity.value, null);
  assert.equal(form.validity.value, null);
  assert.equal(field(ms).classList.contains('u2-invalid'), false);
  form.dispose();
  name.dispose();
  ms.dispose();
});

smoke('dispose: closes the popup and kills the listeners', async () => {
  const before = Scope.liveCount;
  const ms = select();
  fire(field(ms), 'click');
  await flush();

  const open = popup();
  assert.ok(Scope.liveCount > before);
  ms.dispose();
  assert.equal(open.isConnected, false);
  assert.equal(Scope.liveCount, before);

  fire(field(ms), 'click');
  await flush();
  assert.equal(popup(), null, 'listeners died with the scope');
});
