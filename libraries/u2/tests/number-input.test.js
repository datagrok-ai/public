/* NumberInput parity with the platform's `NumberInput<T>` (core d4 number_input.dart /
   float_input.dart): min/max validate instead of clamping typed text, the slider materializes only
   within real bounds, the clicker steps inside them, and floats pick their own precision. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {NumberInput} from '../src/components/number-input.js';

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

function editor(input) {
  return input.root.querySelector('.u2-number-box input');
}

function slider(input) {
  return input.root.querySelector('.u2-number-slider');
}

function clickers(input) {
  return input.root.querySelectorAll('.u2-number-click');
}

function type(input, text) {
  const el = editor(input);
  el.value = text;
  fire(el, 'input');
  return el;
}

smoke('out-of-range text stays on screen, reaches the value, and only marks the input invalid', () => {
  const input = mount(new NumberInput({label: 'Dose', min: 0, max: 100}));
  type(input, '1000');
  assert.equal(input.value.value, 1000, 'the typed number reaches the signal, as it does in Dart');
  assert.equal(input.validity.value, 'Value must be at most 100');
  assert.equal(editor(input).value, '1000');

  fire(editor(input), 'blur');
  assert.equal(input.value.value, 1000, 'blur does not clamp');
  assert.equal(editor(input).value, '1000');

  type(input, '-1');
  assert.equal(input.validity.value, 'Value must be at least 0');

  type(input, '50');
  assert.equal(input.validity.value, null);
  input.dispose();
});

smoke('the stepper is bounded: it walks back into range and clears the validity', () => {
  const input = mount(new NumberInput({mode: 'int', min: 0, max: 10, step: 2, spinner: true}));
  type(input, '99');
  assert.equal(input.validity.value, 'Value must be at most 10');
  fire(input.root.querySelector('.u2-number-spin'), 'click');
  assert.equal(input.value.value, 10);
  assert.equal(editor(input).value, '10');
  assert.equal(input.validity.value, null);
  input.dispose();
});

smoke('the slider materializes only within real bounds', () => {
  const bounded = mount(new NumberInput({min: 0, max: 100, slider: true}));
  assert.notEqual(slider(bounded), null);
  bounded.dispose();

  const unbounded = mount(new NumberInput({min: 0, slider: true}));
  assert.equal(slider(unbounded), null, 'one bound is not enough');
  unbounded.dispose();

  const empty = mount(new NumberInput({min: 5, max: 5, slider: true}));
  assert.equal(slider(empty), null, 'an empty range gets no slider');
  empty.dispose();

  const off = mount(new NumberInput({min: 0, max: 10}));
  assert.equal(slider(off), null, 'and it is opt-in');
  off.dispose();
});

smoke('the slider step defaults per mode and drags as a user edit', () => {
  const float = mount(new NumberInput({min: 0, max: 100, slider: true, value: 10}));
  assert.equal(slider(float).min, '0');
  assert.equal(slider(float).max, '100');
  assert.equal(slider(float).step, '1', 'floats default to a hundredth of the range');
  assert.equal(slider(float).value, '10', 'and the thumb starts on the value');

  slider(float).value = '75';
  fire(slider(float), 'input');
  assert.equal(float.value.value, 75);
  assert.equal(editor(float).value, '75', 'the text box follows the thumb');

  float.value.value = 20;
  assert.equal(slider(float).value, '20', 'and the thumb follows a programmatic write');
  type(float, '30');
  assert.equal(slider(float).value, '30', 'and typing');
  float.dispose();

  const int = mount(new NumberInput({mode: 'int', min: 0, max: 10, slider: true}));
  assert.equal(slider(int).step, '1');
  int.dispose();

  const stepped = mount(new NumberInput({min: 0, max: 1, step: 0.5, slider: true}));
  assert.equal(slider(stepped).step, '0.5', 'an explicit step wins');
  stepped.dispose();
});

smoke('the clicker sits on the options rail and steps inside the bounds', () => {
  const input = mount(new NumberInput({mode: 'int', min: 1, max: 10, clicker: true, value: 9}));
  const [minus, plus] = clickers(input);
  assert.equal(clickers(input).length, 2);
  assert.equal(minus.parentNode, input.root.querySelector('.u2-input-options'),
    'on the rail, as Dart adds it through addOptions');
  assert.equal(input.root.querySelector('.u2-number-spinner'), null, 'and the spinner stays opt-in');

  fire(plus, 'click');
  assert.equal(input.value.value, 10);
  fire(plus, 'click');
  assert.equal(input.value.value, 10, 'the clicker clamps at max');
  fire(minus, 'click');
  assert.equal(input.value.value, 9);
  input.dispose();

  const stepped = mount(new NumberInput({min: 0, max: 10, step: 0.5, clicker: true, value: 1}));
  fire(clickers(stepped)[1], 'click');
  assert.equal(stepped.value.value, 1.5, 'stepping by `step`');
  stepped.dispose();
});

smoke('format drives the display but never rewrites text under the cursor', () => {
  const input = mount(new NumberInput({format: (v) => `${v.toFixed(1)}`, value: 2}));
  assert.equal(editor(input).value, '2.0');

  type(input, '3.14159');
  assert.equal(editor(input).value, '3.14159', 'typing is left alone');
  assert.equal(input.value.value, 3.14159);

  fire(editor(input), 'blur');
  assert.equal(editor(input).value, '3.1', 'and formatted on commit');

  input.value.value = 5;
  assert.equal(editor(input).value, '5.0');
  input.dispose();
});

smoke('floats pick their precision off value, bounds and step', () => {
  const stepped = mount(new NumberInput({step: 0.25, value: 0.5}));
  assert.equal(editor(stepped).value, '0.50', 'as many decimals as the most precise of them');
  stepped.dispose();

  const whole = mount(new NumberInput({min: 0, max: 10, value: 2}));
  assert.equal(editor(whole).value, '2', 'nothing fractional, nothing added');
  whole.dispose();

  const tenths = mount(new NumberInput({step: 0.1, value: 2}));
  assert.equal(editor(tenths).value, '2.0');
  tenths.dispose();

  const small = mount(new NumberInput({value: 1e-7}));
  assert.equal(editor(small).value, '1e-7', 'and tiny values keep their exponent');
  small.dispose();

  const int = mount(new NumberInput({mode: 'int', step: 2, value: 5}));
  assert.equal(editor(int).value, '5', 'ints are never given decimals');
  int.dispose();
});

smoke('postfix renders on the options rail, right of the editor', () => {
  const input = mount(new NumberInput({label: 'Dose', postfix: 'mg', value: 250}));
  const postfix = input.root.querySelector('.u2-input-postfix');
  const rail = input.root.querySelector('.u2-input-options');
  assert.equal(postfix.textContent, 'mg');
  assert.equal(postfix.parentNode, rail);
  assert.equal(input.root.querySelector('.u2-input-editor').nextSibling, rail);
  assert.equal(rail.parentNode, input.root.querySelector('.u2-input-box'));
  input.dispose();
});

smoke('a bare number input is a bare field: no spinner, no clicker, no rail', () => {
  const input = mount(new NumberInput({label: 'Dose', min: 0, max: 10}));
  assert.equal(input.root.querySelector('.u2-number-spinner'), null);
  assert.equal(clickers(input).length, 0);
  assert.equal(input.root.querySelector('.u2-input-options'), null, 'nothing to put there');
  input.dispose();

  const opted = mount(new NumberInput({label: 'Dose', spinner: true}));
  assert.notEqual(opted.root.querySelector('.u2-number-spinner'), null);
  fire(opted.root.querySelector('.u2-number-spin'), 'click');
  assert.equal(opted.value.value, null, 'and it is the wave-1 spinner, empty field included');
  opted.dispose();
});

smoke('the slider lives inside the editor box, under the field', () => {
  const input = mount(new NumberInput({label: 'Dose', min: 0, max: 100, slider: true, value: 5}));
  assert.equal(slider(input).parentNode, input.root.querySelector('.u2-number'),
    'positioned against the field box, as sliderDiv holds it in Dart');
  input.dispose();
});

smoke('the stepper and the slider follow the value, not the text a format rewrote', () => {
  const input = mount(new NumberInput({min: 0, max: 5000, step: 1, slider: true, clicker: true,
    format: (v) => v.toLocaleString('en-US'), value: 1234}));
  assert.equal(editor(input).value, '1,234', 'a format the parser would reject');
  assert.equal(slider(input).value, '1234');

  fire(clickers(input)[1], 'click');
  assert.equal(input.value.value, 1235, 'the stepper reads the signal');
  assert.equal(slider(input).value, '1235');
  input.dispose();
});

smoke('the stepper does nothing on an empty field, as the platform clicker does', () => {
  const input = mount(new NumberInput({min: 0, max: 10, clicker: true}));
  fire(clickers(input)[1], 'click');
  assert.equal(input.value.value, null);
  assert.equal(editor(input).value, '');

  fire(clickers(input)[0], 'click');
  assert.equal(input.value.value, null);
  input.dispose();
});

smoke('an exponent carries no decimals to count', () => {
  const input = mount(new NumberInput({value: 1.5e-7}));
  assert.equal(editor(input).value, '1.5e-7');
  input.dispose();
});
