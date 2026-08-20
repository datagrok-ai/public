/* RadioInput: the native group (one name, checked ⇄ signal), setItems, and the button variant. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {RadioInput} from '../src/components/radio-input.js';

const STAGES = ['Discovery', 'Preclinical', 'Clinical'];

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

function radios(input) {
  return input.root.querySelectorAll('.u2-radio-item input');
}

function checked(input) {
  return radios(input).filter((r) => r.checked).map((r) => r.value);
}

/** What a click on a native radio leaves behind: one checked button, then a change event. */
function pick(input, value) {
  const all = radios(input);
  for (const radio of all)
    radio.checked = radio.value === value;
  fire(all.find((r) => r.value === value), 'change');
}

smoke('one group name, labels, and no selection until something is picked', async () => {
  const input = mount(new RadioInput({label: 'Stage', items: STAGES}));
  assert.equal(input.root.dataset.u2, 'radio-input');
  assert.equal(input.root.querySelector('.u2-radio').getAttribute('role'), 'radiogroup');
  assert.deepEqual(radios(input).map((r) => r.value), STAGES);
  assert.deepEqual(input.root.querySelectorAll('.u2-radio-item span').map((s) => s.textContent), STAGES);
  assert.equal(new Set(radios(input).map((r) => r.name)).size, 1, 'one shared name makes it one group');
  assert.equal(input.value.value, null);
  assert.deepEqual(checked(input), []);

  const other = mount(new RadioInput({label: 'Also stage', items: STAGES}));
  assert.notEqual(radios(other)[0].name, radios(input)[0].name, 'two inputs never share a name');
  input.dispose();
  other.dispose();
});

smoke('picking writes the value, and a programmatic write checks the button', async () => {
  const changes = [];
  const input = mount(new RadioInput({label: 'Stage', items: STAGES, value: 'Discovery',
    onChanged: (v) => changes.push(v)}));
  assert.deepEqual(checked(input), ['Discovery'], 'the initial value reaches the buttons');

  pick(input, 'Clinical');
  assert.equal(input.value.value, 'Clinical');
  assert.deepEqual(changes, ['Clinical']);

  input.value.value = 'Preclinical';
  assert.deepEqual(checked(input), ['Preclinical'], 'a programmatic write checks the button');

  input.value.value = null;
  assert.deepEqual(checked(input), [], 'and null clears the group');
  input.dispose();
});

smoke('setItems keeps the value while the item survives, drops it otherwise', async () => {
  const input = mount(new RadioInput({label: 'Stage', items: STAGES, value: 'Preclinical'}));

  input.setItems(['Preclinical', 'Clinical']);
  assert.equal(input.value.value, 'Preclinical');
  assert.deepEqual(checked(input), ['Preclinical'], 'the rebuilt buttons carry the selection');

  input.setItems(['Discovery', 'Clinical']);
  assert.equal(input.value.value, null, 'a vanished item cannot stay selected');
  assert.deepEqual(checked(input), []);

  pick(input, 'Clinical');
  assert.equal(input.value.value, 'Clinical', 'the group listener survives a rebuild');
  input.dispose();
});

smoke('buttons variant marks the picked item; tooltips land on the items', async () => {
  const input = mount(new RadioInput({label: 'Stage', items: STAGES, buttons: true, value: 'Discovery',
    itemTooltips: {Discovery: 'Hit finding'}}));
  const items = input.root.querySelectorAll('.u2-radio-item');
  assert.equal(input.root.querySelector('.u2-radio').classList.contains('u2-radio-buttons'), true);
  assert.equal(items[0].classList.contains('u2-radio-pressed'), true);
  assert.equal(items[2].classList.contains('u2-radio-pressed'), false);
  assert.equal(items[0].title, 'Hit finding');
  assert.equal(items[1].title, '');

  pick(input, 'Clinical');
  assert.equal(items[0].classList.contains('u2-radio-pressed'), false);
  assert.equal(items[2].classList.contains('u2-radio-pressed'), true);
  input.dispose();
});

smoke('validators come from the input base; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const input = mount(new RadioInput({label: 'Stage', items: STAGES}));
  input.addValidator((v) => v === null ? 'Pick a stage' : null);
  assert.equal(input.validity.value, 'Pick a stage');
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Pick a stage');

  pick(input, 'Discovery');
  assert.equal(input.validity.value, null);

  input.dispose();
  assert.equal(Scope.liveCount, before);
  pick(input, 'Clinical');
  assert.equal(input.value.value, 'Discovery', 'listeners died with the scope');
});

smoke('a disabled group stays disabled through setItems', () => {
  const input = mount(new RadioInput({items: STAGES}));
  input.enabled = false;
  input.setItems(['Phase I', 'Phase II']);
  const radios = input.root.querySelectorAll('input');
  assert.equal(radios.length, 2);
  for (let i = 0; i < radios.length; i++)
    assert.equal(radios[i].disabled, true, 'rebuilt radios keep the disabled state');
  input.dispose();
});
