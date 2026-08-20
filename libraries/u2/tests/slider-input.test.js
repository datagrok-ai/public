/* SliderInput: the track contract (bounds, default step, drag ⇄ signal) and the drag tooltip —
   a bare track, as Dart's `SliderInput` is, with the value only in the tooltip. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Tooltip} from '../src/index.js';
import {SliderInput} from '../src/components/slider-input.js';

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

function track(input) {
  return input.root.querySelector('.u2-slider-track');
}

function readout(input) {
  return input.text;
}

function drag(input, value) {
  const el = track(input);
  el.value = String(value);
  fire(el, 'input');
}

smoke('bounds and the default step come from min/max', async () => {
  const input = mount(new SliderInput({label: 'Opacity', min: 0, max: 50, value: 10}));
  assert.equal(track(input).type, 'range');
  assert.equal(track(input).min, '0');
  assert.equal(track(input).max, '50');
  assert.equal(track(input).step, '0.5', '(max - min) / 100');
  assert.equal(track(input).value, '10');

  const stepped = mount(new SliderInput({label: 'Count', min: 0, max: 10, step: 2}));
  assert.equal(track(stepped).step, '2', 'an explicit step wins');

  const bare = mount(new SliderInput({label: 'Bare'}));
  assert.equal(track(bare).min, '0', '0…100 by default, as a bare range element is');
  assert.equal(track(bare).max, '100');
  assert.equal(track(bare).step, '1');
  assert.equal(bare.value.value, null);
  assert.equal(track(bare).value, '0', 'a null value parks the thumb at the minimum');
  assert.equal(readout(bare), '', '…and has no number to show');
  assert.equal(bare.root.querySelector('.u2-slider-value'), null,
    'no readout beside the track — the platform has none either');
  input.dispose();
  stepped.dispose();
  bare.dispose();
});

smoke('dragging writes the signal, a programmatic write moves the thumb', async () => {
  const changes = [];
  const input = mount(new SliderInput({label: 'Dose', min: 0, max: 100, value: 25,
    onChanged: (v) => changes.push(v)}));
  assert.equal(readout(input), '25');

  drag(input, 60);
  assert.equal(input.value.value, 60);
  assert.equal(readout(input), '60');
  assert.equal(Tooltip.root.textContent, '60', 'the value surfaces in the tooltip while dragging');
  assert.deepEqual(changes, [60]);
  Tooltip.hide();

  input.value.value = 5;
  assert.equal(track(input).value, '5', 'a programmatic write echoes into the track');
  assert.equal(readout(input), '5');
  assert.deepEqual(changes, [60, 5]);
  input.dispose();
});

smoke('format drives the readout only; the value stays a number', async () => {
  const input = mount(new SliderInput({label: 'Threshold', min: 0, max: 1, step: 0.01, value: 0.5,
    format: (v) => `${Math.round(v * 100)}%`}));
  assert.equal(readout(input), '50%');

  drag(input, 0.25);
  assert.equal(input.value.value, 0.25);
  assert.equal(readout(input), '25%');
  assert.equal(Tooltip.root.textContent, '25%', 'the format drives the tooltip');
  Tooltip.hide();
  input.dispose();
});

smoke('two sliders on one signal stay in step; dispose kills the listeners', async () => {
  const before = Scope.liveCount;
  const first = mount(new SliderInput({label: 'From', min: 0, max: 10, value: 3}));
  const second = mount(new SliderInput({label: 'Also from', min: 0, max: 10, bind: first.value}));
  assert.equal(track(second).value, '3');

  drag(first, 8);
  assert.equal(second.value.value, 8);
  assert.equal(track(second).value, '8');
  assert.equal(readout(second), '8');
  Tooltip.hide();

  second.dispose();
  first.dispose();
  assert.equal(Scope.liveCount, before);

  drag(first, 1);
  assert.equal(first.value.value, 8, 'listeners died with the scope');
});
