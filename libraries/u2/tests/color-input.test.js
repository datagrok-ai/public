/* ColorInput: the hex text contract (parse, validity, last-good swatch), the native picker
   behind the swatch, and the swatch-only variant. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {ColorInput} from '../src/components/color-input.js';

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

function text(input) {
  return input.root.querySelector('.u2-color-text');
}

function swatch(input) {
  return input.root.querySelector('.u2-color-swatch');
}

function picker(input) {
  return input.root.querySelector('.u2-color-picker');
}

function type(input, value) {
  const el = text(input);
  el.value = value;
  fire(el, 'input');
}

smoke('typing a hex color paints the swatch and writes the value', async () => {
  const changes = [];
  const input = mount(new ColorInput({label: 'Color', onChanged: (v) => changes.push(v)}));
  assert.equal(input.root.dataset.u2, 'color-input');
  assert.equal(input.value.value, '#000000', 'black by default, as Dart ColorInput is');
  assert.equal(text(input).value, '#000000');
  assert.equal(swatch(input).style.background, '#000000');
  assert.equal(picker(input).type, 'color');

  type(input, '#ff0000');
  assert.equal(input.value.value, '#ff0000');
  assert.equal(swatch(input).style.background, '#ff0000');
  assert.equal(picker(input).value, '#ff0000', 'the native picker follows the text');
  assert.equal(input.validity.value, null);
  assert.deepEqual(changes, ['#ff0000']);

  type(input, '#0F0');
  assert.equal(input.value.value, '#0F0', 'the short form is a color too');
  assert.equal(picker(input).value, '#00ff00', '…expanded for the picker, which takes no short form');
  input.dispose();
});

smoke('text that is not a color marks the input invalid and keeps the last good color', async () => {
  const input = mount(new ColorInput({label: 'Color', value: '#123456'}));

  type(input, 'crimson');
  assert.equal(input.validity.value, 'Not a color');
  assert.equal(input.value.value, '#123456', 'the value signal keeps the last good color');
  assert.equal(swatch(input).style.background, '#123456');
  assert.equal(text(input).value, 'crimson', 'the typed text stays on screen');
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Not a color');

  type(input, '#12345');
  assert.equal(input.validity.value, 'Not a color', 'five digits is not a hex color');

  type(input, '#abcdef');
  assert.equal(input.validity.value, null);
  assert.equal(input.value.value, '#abcdef');
  input.dispose();
});

smoke('the swatch opens the native picker, and the picker writes text, swatch and value', async () => {
  const input = mount(new ColorInput({label: 'Color'}));
  const opened = [];
  picker(input).addEventListener('click', () => opened.push(1));

  fire(swatch(input), 'click');
  assert.equal(opened.length, 1, 'the swatch clicks the hidden picker');

  picker(input).value = '#336699';
  fire(picker(input), 'input');
  assert.equal(input.value.value, '#336699');
  assert.equal(text(input).value, '#336699', 'the text follows the picker');
  assert.equal(swatch(input).style.background, '#336699');
  input.dispose();
});

smoke('the hex field is the editor and the swatch rides the options rail, right of it', async () => {
  const input = mount(new ColorInput({label: 'Color'}));
  assert.equal(text(input).classList.contains('u2-input-editor'), true);
  const rail = input.root.querySelector('.u2-input-options');
  assert.equal(rail.contains(swatch(input)), true, 'as Dart puts colorDiv in the options div');
  assert.equal(text(input).nextSibling, rail);
  input.dispose();
});

smoke('swatchOnly hides the text field; a programmatic write still repaints', async () => {
  const input = mount(new ColorInput({label: 'Color', value: '#ff0000', swatchOnly: true}));
  assert.equal(text(input).hidden, true);
  assert.equal(input.root.classList.contains('u2-color-box-only'), true);

  input.value.value = '#00ff00';
  assert.equal(swatch(input).style.background, '#00ff00');
  assert.equal(text(input).value, '#00ff00');
  input.dispose();
});

smoke('two inputs on one signal stay in step; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const first = mount(new ColorInput({label: 'Fill', value: '#111111'}));
  const second = mount(new ColorInput({label: 'Also fill', bind: first.value}));
  assert.equal(text(second).value, '#111111');

  type(first, '#222222');
  assert.equal(text(second).value, '#222222');
  assert.equal(swatch(second).style.background, '#222222');

  second.dispose();
  first.dispose();
  assert.equal(Scope.liveCount, before);

  type(first, '#333333');
  assert.equal(first.value.value, '#222222', 'listeners died with the scope');
});
