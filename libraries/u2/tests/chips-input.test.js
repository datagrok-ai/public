/* ChipsInput: one pressable chip per item, the value in item order, setItems pruning as a
   system write, a disabled input ignoring clicks; registered as u2-chips-input. */

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

const chips = (input) => input.root.querySelectorAll('.u2-chip');
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
