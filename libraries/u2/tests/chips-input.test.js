/* ChipsInput: one pressable chip per item, the value in item order, setItems pruning as a
   system write, a disabled input ignoring clicks, chips kept by value across changes and addressable
   as items; registered as u2-chips-input. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Input} from '../src/index.js';
import {ChipsInput} from '../src/components/inputs/chips-input.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';

function ui(name, body) {
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

const chips = (input) => input.root.querySelectorAll('.u2-chips-chip');
const pressed = (input) => chips(input).filter((c) => c.getAttribute('aria-pressed') === 'true').map((c) => c.textContent);

ui('a chip per item, pressed where the value says; a click toggles, the value keeps item order', () => {
  const input = new ChipsInput({label: 'Visible to', items: ['Sales', 'Developers', 'Chemists'], value: ['Chemists']});
  assert.equal(input.root.dataset.u2, 'chips-input');
  assert.deepEqual(chips(input).map((c) => c.textContent), ['Sales', 'Developers', 'Chemists']);
  assert.deepEqual(pressed(input), ['Chemists']);
  fire(chips(input)[0], 'click');
  assert.deepEqual(input.value.value, ['Sales', 'Chemists']);
  assert.deepEqual(pressed(input), ['Sales', 'Chemists']);
  fire(chips(input)[2], 'click');
  assert.deepEqual(input.value.value, ['Sales']);
  input.dispose();
});

ui('setItems prunes a pick whose item vanished as a system write; an empty list says so', () => {
  const input = new ChipsInput({label: 'Visible to', items: ['Sales', 'Chemists'], value: ['Sales', 'Chemists']});
  let system = null;
  input.effect(() => {
    input.value.value;
    system = Input.isSystemWrite;
  });
  input.setItems(['Chemists']);
  assert.deepEqual(input.value.value, ['Chemists']);
  assert.equal(system, true);
  input.setItems([]);
  assert.equal(input.root.querySelector('.u2-chips-empty').textContent, 'No items');
  input.dispose();
});

ui('the pressed chip is the same element, still focused, so a second Space toggles it back', () => {
  const input = new ChipsInput({label: 'Visible to', items: ['Sales', 'Developers'], value: []});
  document.body.append(input.root);
  const chip = chips(input)[1];
  chip.focus();
  fire(chip, 'click');
  assert.deepEqual(input.value.value, ['Developers']);
  assert.equal(chips(input)[1] === chip, true, 'kept, not rebuilt');
  assert.equal(chip.classList.contains('u2-chips-chip-on'), true);
  assert.equal(document.activeElement === chip, true);
  fire(chip, 'click');
  assert.deepEqual(input.value.value, []);
  assert.equal(chip.getAttribute('aria-pressed'), 'false');
  input.setItems(['Chemists', 'Developers', {value: 'Sales', label: 'Sales team'}]);
  assert.equal(chips(input)[1] === chip, true, 'a chip whose item stays keeps its element across setItems');
  assert.equal(document.activeElement === chip, true);
  assert.deepEqual(chips(input).map((c) => c.textContent), ['Chemists', 'Developers', 'Sales team'],
    'a relabelled item is patched in place');
  input.dispose();
});

ui('chips are addressable items: data-u2-item, role, item-label; the pressed ones are the selection', async () => {
  const input = new ChipsInput({label: 'Visible to', items: ['Sales', {value: 'dev', label: 'Developers'}],
    value: ['dev']});
  document.body.append(input.root);
  assert.deepEqual(chips(input).map((c) => [c.dataset.u2Item, c.getAttribute('role'),
    c.querySelector('[data-u2-part="item-label"]').textContent]), [['Sales', 'button', 'Sales'], ['dev', 'button', 'Developers']]);
  assert.equal(input.root.querySelector('.u2-chip'), null, 'the entity chip class is not used');
  const state = input.getState();
  assert.deepEqual(state.items.list.map((i) => [i.key, i.label, i.selected]), [['Sales', 'Sales', false], ['dev', 'Developers', true]]);
  assert.deepEqual(state.selection, ['dev']);
  assert.equal(state.text, 'Developers');
  assert.equal(state.readOnly, false);
  assert.equal(state.required, false);
  assert.equal(await input.revealItem('developers') === chips(input)[1], true);
  input.dispose();
});

ui('a disabled input ignores clicks', () => {
  const input = new ChipsInput({label: 'Visible to', items: ['Sales'], enabled: false});
  fire(chips(input)[0], 'click');
  assert.deepEqual(input.value.value, []);
  assert.equal(chips(input)[0].disabled, true);
  input.dispose();
});

ui('registered as u2-chips-input', () => {
  const reg = new Registry();
  registerAll(reg);
  const built = reg.get('u2-chips-input').create({label: 'Visible to', items: ['A', 'B'], value: ['B']});
  assert.deepEqual(pressed(built), ['B']);
  built.dispose();
});
