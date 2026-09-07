/* Tour: the walkthrough service contract — one masked overlay plus an anchored popup in the shared
   layer, spec-name targeting, skip-forward over targets that are not there, the NEXT/SKIP/DONE
   keyboard surface, auto-advance, and a teardown that leaves neither DOM nor live scopes behind. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Overlay, Scope, signal} from '../src/index.js';
import {Tour} from '../src/components/display/tour.js';

// the shim has no SVG namespace factory; the mask is plain elements here, attributes and all
document.createElementNS ??= (_ns, tag) => document.createElement(tag);

let running = null;

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      running?.finish();
      running = null;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function start(options) {
  running = Tour.run(options);
  return running;
}

/** A visible element the tour can spotlight: the shim answers rects from `.rect`. */
function target(name, x = 10, y = 20, width = 100, height = 30, parent = document.body) {
  const el = document.createElement('div');
  el.dataset.u2Name = name;
  el.rect = new DOMRect(x, y, width, height);
  parent.append(el);
  return el;
}

const box = (rect) => ['x', 'y', 'width', 'height'].map((a) => Number(rect.getAttribute(a)));
const holes = (tour) => tour.overlay.querySelectorAll('.u2-tour-hole');
const counter = (tour) => tour.popup.querySelector('.u2-tour-counter').textContent;
const nextButton = (tour) => tour.popup.querySelector('.u2-tour-next');
const skipButton = (tour) => tour.popup.querySelector('.u2-tour-skip');

smoke('run mounts the overlay and the popup in the layer host and punches the target hole', async () => {
  target('nameInput', 10, 20, 100, 30);
  const tour = start({steps: [{target: 'nameInput', content: 'Name your model here.'}]});

  assert.equal(tour.overlay.parentNode, Overlay.host);
  assert.equal(tour.popup.parentNode, Overlay.host);
  assert.equal(tour.active.value, true);
  assert.ok(Number(tour.popup.style.zIndex) > Number(tour.overlay.style.zIndex),
    'the popup paints above the dim');
  assert.equal(tour.popup.querySelector('.u2-tour-content').textContent, 'Name your model here.');
  // 4px of padding on every side
  assert.deepEqual(box(holes(tour)[0]), [6, 16, 108, 38]);
  assert.deepEqual(box(holes(tour)[1]), [0, 0, 0, 0], 'no extra: the second hole is empty');
});

smoke('a string target is a spec name looked up inside root; an element target is used as is', async () => {
  const root = document.createElement('div');
  document.body.append(root);
  target('field', 200, 300, 40, 40, root);
  target('field', 1, 1, 5, 5);
  const outside = target('other', 500, 600, 60, 20);

  const tour = start({root, steps: [
    {target: 'field', content: 'Inside the root.'},
    {target: outside, content: 'By element, wherever it lives.', extra: 'field'},
  ]});
  assert.deepEqual(box(holes(tour)[0]), [196, 296, 48, 48], 'the root-scoped element won');

  tour.next();
  assert.deepEqual(box(holes(tour)[0]), [496, 596, 68, 28]);
  assert.deepEqual(box(holes(tour)[1]), [196, 296, 48, 48], 'extra punches a second hole');
});

smoke('a target that is not there is skipped forward, warned about once, and never blocks', async () => {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(message);
  try {
    target('last', 40, 50, 20, 20);
    // a present but unrendered element is as good as a missing one
    target('zeroSized', 0, 0, 0, 0);
    const tour = start({steps: [
      {target: 'gone', content: 'Never shown.'},
      {target: 'zeroSized', content: 'Never shown either.'},
      {target: 'last', content: 'Shown.'},
    ]});

    assert.equal(tour.step.value, 2);
    assert.deepEqual(warnings, ['u2: tour target "gone" not found — step skipped',
      'u2: tour target "zeroSized" not found — step skipped']);

    tour.next();
    assert.equal(tour.active.value, false, 'nothing renderable left: the tour is done');
    assert.equal(warnings.length, 2, 'one warning per target, not per attempt');
  } finally {
    console.warn = original;
  }
});

smoke('the counter walks the steps and the last one says DONE', async () => {
  target('a', 10, 10, 30, 30);
  target('b', 60, 10, 30, 30);
  const results = [];
  const tour = start({steps: [
    {target: 'a', content: 'First.'},
    {target: 'b', content: 'Second.'},
  ], onFinish: (r) => results.push(r)});

  assert.equal(counter(tour), '1 / 2');
  assert.equal(nextButton(tour).textContent, 'NEXT');
  assert.equal(skipButton(tour).textContent, 'SKIP');

  fire(nextButton(tour), 'click');
  assert.equal(tour.step.value, 1);
  assert.equal(counter(tour), '2 / 2');
  assert.equal(nextButton(tour).textContent, 'DONE');
  assert.deepEqual(box(holes(tour)[0]), [56, 6, 38, 38], 're-anchored on the new target');

  fire(nextButton(tour), 'click');
  assert.deepEqual(results, ['done']);
});

smoke('custom button labels override the defaults', async () => {
  target('a', 10, 10, 30, 30);
  const tour = start({steps: [{target: 'a', content: 'Only step.'}],
    buttons: {next: 'Continue', skip: 'Not now', done: 'Got it'}});
  assert.equal(skipButton(tour).textContent, 'Not now');
  assert.equal(nextButton(tour).textContent, 'Got it', 'a single step is the last one');
});

smoke('Esc skips, Enter advances, and Enter inside a text field does not', async () => {
  target('a', 10, 10, 30, 30);
  target('b', 60, 10, 30, 30);
  const input = document.createElement('input');
  document.body.append(input);
  const results = [];
  const options = {steps: [{target: 'a', content: 'First.'}, {target: 'b', content: 'Second.'}],
    onFinish: (r) => results.push(r)};

  const tour = start(options);
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(tour.step.value, 0, 'Enter in a field types, it does not advance');

  fire(document.body, 'keydown', {key: 'Enter'});
  assert.equal(tour.step.value, 1);

  fire(document.body, 'keydown', {key: 'Escape'});
  assert.deepEqual(results, ['skipped']);
  assert.equal(tour.overlay.isConnected, false);

  const second = start(options);
  fire(document.body, 'keydown', {key: 'Escape'});
  assert.equal(second.active.value, false);
  assert.deepEqual(results, ['skipped', 'skipped']);
});

smoke('Tab cycles inside the popup and the step opens with NEXT focused', async () => {
  target('a', 10, 10, 30, 30);
  const tour = start({steps: [{target: 'a', content: 'Only step.'}]});
  assert.equal(document.activeElement, nextButton(tour));

  fire(nextButton(tour), 'keydown', {key: 'Tab'});
  assert.equal(document.activeElement, skipButton(tour), 'wrapped forward to the first control');

  fire(skipButton(tour), 'keydown', {key: 'Tab', shiftKey: true});
  assert.equal(document.activeElement, nextButton(tour), 'and back again');
});

smoke('advanceOn advances once, and its step scope dies when the step is left', async () => {
  target('a', 10, 10, 30, 30);
  target('b', 60, 10, 30, 30);
  target('c', 110, 10, 30, 30);
  const ran = signal(false);
  const tour = start({steps: [
    {target: 'a', content: 'Click it.', advanceOn: ran},
    {target: 'b', content: 'Second.'},
    {target: 'c', content: 'Third.'},
  ]});

  assert.equal(tour.step.value, 0, 'a falsy signal keeps the step');
  ran.value = true;
  assert.equal(tour.step.value, 1);

  ran.value = false;
  ran.value = true;
  assert.equal(tour.step.value, 1, 'the left step no longer listens');
});

smoke('advanceOn also takes a promise, and a stale resolve is ignored', async () => {
  target('a', 10, 10, 30, 30);
  target('b', 60, 10, 30, 30);
  let arrive = () => {};
  const tour = start({steps: [
    {target: 'a', content: 'Waiting.', advanceOn: new Promise((resolve) => arrive = resolve)},
    {target: 'b', content: 'Second.'},
  ]});

  tour.next();
  arrive();
  await flush();
  assert.equal(tour.step.value, 1, 'the promise of a step already left does not advance again');
});

smoke('starting a second tour finishes the running one', async () => {
  target('a', 10, 10, 30, 30);
  const results = [];
  const first = start({steps: [{target: 'a', content: 'First tour.'}],
    onFinish: (r) => results.push(`first:${r}`)});
  const second = start({steps: [{target: 'a', content: 'Second tour.'}],
    onFinish: (r) => results.push(`second:${r}`)});

  assert.equal(first.active.value, false);
  assert.equal(first.overlay.isConnected, false);
  assert.equal(second.active.value, true);
  assert.deepEqual(results, ['first:skipped']);
});

smoke('finish is idempotent, cleans the layer host and gives focus back', async () => {
  const before = document.createElement('button');
  document.body.append(before);
  before.focus();
  target('a', 10, 10, 30, 30);

  const results = [];
  const tour = start({steps: [{target: 'a', content: 'Only step.'}], onFinish: (r) => results.push(r)});
  assert.notEqual(document.activeElement, before, 'the popup took focus');

  tour.finish();
  tour.finish();
  tour.next();
  assert.equal(tour.overlay.isConnected, false);
  assert.equal(tour.popup.isConnected, false);
  assert.equal(Overlay.host.querySelector('.u2-tour-overlay'), null);
  assert.equal(document.activeElement, before, 'focus is back where the tour found it');
  assert.deepEqual(results, ['done'], 'onFinish fires once');
  await flush();
});
