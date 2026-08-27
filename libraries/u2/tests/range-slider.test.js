/* RangeSlider: handle and band dragging over a faked track rect, snapping and clamping, the
   vertical orientation, keyboard steps, two-way signals, and listener death on dispose. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, RangeSlider, signal} from '../src/index.js';

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

/** Mounts a slider whose 200px track starts at 0 — clientX 2×value for a 0…100 range. */
function mount(options) {
  const slider = new RangeSlider(options);
  document.body.append(slider.root);
  const track = slider.root.querySelector('.u2-range-track');
  track.rect = options?.vertical ? new DOMRect(0, 0, 16, 200) : new DOMRect(0, 0, 200, 16);
  return slider;
}

function handle(slider, kind) {
  return slider.root.querySelector(`.u2-range-handle-${kind}`);
}

function move(clientX, clientY = 8) {
  window.dispatchEvent(new PointerEvent('pointermove', {clientX, clientY}));
}

function up() {
  window.dispatchEvent(new PointerEvent('pointerup', {}));
}

smoke('defaults and rendering: ends park at the bounds, styles follow the signals', async () => {
  const slider = mount({});
  assert.equal(slider.lo.value, 0);
  assert.equal(slider.hi.value, 100);

  slider.setValues(25, 75);
  await flush();
  assert.equal(handle(slider, 'lo').style.left, '25%');
  assert.equal(handle(slider, 'hi').style.left, '75%');
  const band = slider.root.querySelector('.u2-range-band');
  assert.equal(band.style.left, '25%');
  assert.equal(band.style.width, '50%');
  assert.equal(handle(slider, 'lo').getAttribute('aria-valuenow'), '25');
  assert.equal(handle(slider, 'lo').getAttribute('role'), 'slider');
  slider.dispose();
});

smoke('a grabbed handle drags relative to the grab point; the other end clamps it', async () => {
  const changes = [];
  const slider = mount({lo: 20, hi: 60, onChanged: (lo, hi) => changes.push([lo, hi])});
  const events = [];
  slider.root.addEventListener('input', () => events.push('input'));
  slider.root.addEventListener('change', () => events.push('change'));

  fire(handle(slider, 'lo'), 'pointerdown', {clientX: 50, clientY: 8});
  assert.equal(slider.lo.value, 20, 'grabbing a handle never teleports it');
  move(70);
  assert.equal(slider.lo.value, 30, 'the grab offset rides along');
  move(180);
  assert.equal(slider.lo.value, 60, 'lo never passes hi');
  up();
  assert.deepEqual(events.filter((e) => e === 'change'), ['change']);
  assert.ok(events.includes('input'));
  assert.deepEqual(changes.at(-1), [60, 60]);

  move(20);
  assert.equal(slider.lo.value, 60, 'listeners are gone after pointerup');
  slider.dispose();
});

smoke('step snaps dragging; minRange keeps the ends apart', async () => {
  const slider = mount({step: 10, minRange: 20, lo: 20, hi: 60});
  fire(handle(slider, 'lo'), 'pointerdown', {clientX: 40, clientY: 8});
  move(94);
  assert.equal(slider.lo.value, 40, '20 + 27 = 47 snaps to 50, minRange clamps it down to hi − 20');
  up();
  slider.dispose();
});

smoke('a press between the handles pans, even off the thin band line', async () => {
  const slider = mount({lo: 20, hi: 60});
  fire(slider.root, 'pointerdown', {clientX: 80, clientY: 1});
  assert.equal(slider.lo.value, 20, 'no jump on grab');
  move(100);
  assert.equal(slider.lo.value, 30, 'the whole selection panned');
  assert.equal(slider.hi.value, 70);
  up();
  slider.dispose();
});

smoke('band pans snap on a step slider', async () => {
  const slider = mount({step: 10, lo: 20, hi: 60});
  fire(slider.root, 'pointerdown', {clientX: 80, clientY: 8});
  move(94);
  assert.equal(slider.lo.value, 30, '20 + 7 snaps to 30');
  assert.equal(slider.hi.value, 70);
  up();
  slider.dispose();
});

smoke('default formatting keeps at most two decimals; an explicit format wins', async () => {
  const slider = mount({lo: 20, hi: 60});
  slider.lo.value = 33.333333;
  await flush();
  assert.equal(handle(slider, 'lo').getAttribute('aria-valuetext'), '33.33');
  slider.lo.value = 25;
  await flush();
  assert.equal(handle(slider, 'lo').getAttribute('aria-valuetext'), '25', 'no forced zeros');
  slider.dispose();

  const formatted = mount({lo: 0.5, hi: 0.9, min: 0, max: 1, format: (v) => `${v * 100}%`});
  await flush();
  assert.equal(handle(formatted, 'lo').getAttribute('aria-valuetext'), '50%');
  formatted.dispose();
});

smoke('band drag moves the window and preserves the span', async () => {
  const slider = mount({lo: 20, hi: 60});
  const band = slider.root.querySelector('.u2-range-band');
  fire(band, 'pointerdown', {clientX: 80, clientY: 8});
  assert.equal(slider.lo.value, 20, 'grabbing the band moves nothing yet');
  move(120);
  assert.equal(slider.lo.value, 40);
  assert.equal(slider.hi.value, 80);
  move(200);
  assert.equal(slider.lo.value, 60, 'the window stops at the edge');
  assert.equal(slider.hi.value, 100);
  up();
  slider.dispose();
});

smoke('track click pulls the nearest handle', async () => {
  const slider = mount({lo: 20, hi: 60});
  const track = slider.root.querySelector('.u2-range-track');
  fire(track, 'pointerdown', {clientX: 180, clientY: 8});
  assert.equal(slider.hi.value, 90, 'the click was nearer hi');
  assert.equal(slider.lo.value, 20);
  up();
  fire(track, 'pointerdown', {clientX: 10, clientY: 8});
  assert.equal(slider.lo.value, 5);
  up();
  slider.dispose();
});

smoke('vertical: bottom is min, top is max; ArrowUp raises the value', async () => {
  const slider = mount({vertical: true, lo: 40, hi: 80});
  assert.ok(slider.root.classList.contains('u2-range-slider-vertical'));
  assert.equal(handle(slider, 'lo').getAttribute('aria-orientation'), 'vertical');

  const track = slider.root.querySelector('.u2-range-track');
  fire(track, 'pointerdown', {clientX: 8, clientY: 190});
  assert.equal(slider.lo.value, 5, '190px down a 200px track is 5% up — below lo, so lo jumps');
  up();
  await flush();
  assert.equal(handle(slider, 'lo').style.bottom, '5%');
  assert.ok(!handle(slider, 'lo').style.left, 'vertical positioning never touches left');

  fire(handle(slider, 'lo'), 'keydown', {key: 'ArrowUp'});
  assert.equal(slider.lo.value, 6);
  slider.dispose();
});

smoke('keyboard: arrows step, Home/End go to the reachable bounds', async () => {
  const slider = mount({lo: 20, hi: 60, minRange: 10});
  const lo = handle(slider, 'lo');
  fire(lo, 'keydown', {key: 'ArrowRight'});
  assert.equal(slider.lo.value, 21);
  fire(lo, 'keydown', {key: 'ArrowLeft'});
  assert.equal(slider.lo.value, 20);
  fire(lo, 'keydown', {key: 'PageUp'});
  assert.equal(slider.lo.value, 30);
  fire(lo, 'keydown', {key: 'Home'});
  assert.equal(slider.lo.value, 0);
  fire(lo, 'keydown', {key: 'End'});
  assert.equal(slider.lo.value, 50, 'End stops minRange short of hi');
  fire(handle(slider, 'hi'), 'keydown', {key: 'End'});
  assert.equal(slider.hi.value, 100);
  slider.dispose();
});

smoke('external signals are adopted, not copied; dispose severs a drag in flight', async () => {
  const lo = signal(10);
  const hi = signal(90);
  const slider = mount({lo, hi});
  assert.equal(slider.lo, lo, 'the slider works on the caller\'s signal');

  lo.value = 30;
  await flush();
  assert.equal(handle(slider, 'lo').style.left, '30%');

  fire(handle(slider, 'hi'), 'pointerdown', {clientX: 170, clientY: 8});
  move(150);
  assert.equal(hi.value, 80, 'the grab kept its 5-unit offset: 75 + 5');
  slider.dispose();
  move(100);
  assert.equal(hi.value, 80, 'window listeners died with the scope');
});

smoke('setValues clamps to the bounds and keeps the order', async () => {
  const slider = mount({minRange: 10});
  slider.setValues(-5, 200);
  assert.equal(slider.lo.value, 0);
  assert.equal(slider.hi.value, 100);
  slider.setValues(50, 20);
  assert.equal(slider.lo.value, 50);
  assert.equal(slider.hi.value, 60, 'hi is pushed up to lo + minRange');
  slider.dispose();
});
