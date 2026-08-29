/* Form layouts: the layout enum's class emission, the captionAlign live toggle, the
   Form.add(input, host) seam, and the auto machinery — _pickLayout's platform truth table
   (threshold + 10px hysteresis, js-api ui.ts:1443-1448) and the class toggle over stubbed
   widths. Layout metrics are dom-shim fakes; nothing here asserts pixels. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {Form} from '../src/components/forms/form.js';
import {TextInput} from '../src/components/inputs/text-input.js';

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

function mount(f) {
  document.body.append(f.root);
  return f;
}

/** A labeled input whose label measures `width` in the fake layout. */
function input(label, width) {
  const it = new TextInput({label, name: label});
  it.root.querySelector('.u2-input-label').offsetWidth = width;
  return it;
}

const tall = (f) => f.root.classList.contains('u2-form-tall');
const wide = (f) => f.root.classList.contains('u2-form-wide');

test('_pickLayout: the fit truth table — tall when squeezed, WIDE when roomy, band and hidden host keep current', () => {
  const pick = Form._pickLayout;
  // hidden host (clientWidth 0) keeps the current state — never a wrong flip
  assert.equal(pick(0, 120, 150, 'wide'), 'wide');
  assert.equal(pick(0, 120, 150, 'tall'), 'tall');
  assert.equal(pick(0, 120, 150, 'normal'), 'normal', 'the never-measured state holds too');
  // too narrow: width - label < minEditor → tall, from any state
  assert.equal(pick(200, 120, 150, 'wide'), 'tall');
  assert.equal(pick(200, 120, 150, 'tall'), 'tall');
  assert.equal(pick(200, 120, 150, 'normal'), 'tall');
  assert.equal(pick(269, 120, 150, 'wide'), 'tall');
  // roomy: width - label > minEditor + 10 → wide (the user-ruled roomy state), from any state
  assert.equal(pick(400, 120, 150, 'wide'), 'wide');
  assert.equal(pick(400, 120, 150, 'tall'), 'wide');
  assert.equal(pick(400, 120, 150, 'normal'), 'wide');
  assert.equal(pick(281, 120, 150, 'tall'), 'wide');
  // the ±10 hysteresis band [minEditor, minEditor + 10] keeps the current state — no flicker
  assert.equal(pick(270, 120, 150, 'wide'), 'wide');
  assert.equal(pick(270, 120, 150, 'tall'), 'tall');
  assert.equal(pick(275, 120, 150, 'wide'), 'wide');
  assert.equal(pick(275, 120, 150, 'tall'), 'tall');
  assert.equal(pick(280, 120, 150, 'wide'), 'wide');
  assert.equal(pick(280, 120, 150, 'tall'), 'tall');
  assert.equal(pick(275, 120, 150, 'normal'), 'normal');
});

ui('layout classes: wide and tall are emitted, auto and normal are not, old options are ignored', () => {
  const auto = new Form();
  assert.equal(auto.root.classList.contains('u2-form-wide'), false);
  assert.equal(tall(auto), false);
  auto.dispose();

  const wide = new Form({layout: 'wide'});
  assert.equal(wide.root.classList.contains('u2-form-wide'), true);
  wide.dispose();

  const forced = new Form({layout: 'tall'});
  assert.equal(tall(forced), true);
  forced.dispose();

  const normal = new Form({layout: 'normal'});
  assert.equal(normal.root.classList.contains('u2-form-wide'), false);
  assert.equal(tall(normal), false);
  normal.dispose();
});

ui('captionAlign: left emits the class, right (the default) does not, a signal toggles it live', () => {
  const right = new Form({captionAlign: 'right'});
  assert.equal(right.root.classList.contains('u2-form-captions-left'), false);
  right.dispose();

  const left = new Form({captionAlign: 'left'});
  assert.equal(left.root.classList.contains('u2-form-captions-left'), true);
  left.dispose();

  const align = signal('right');
  const live = new Form({captionAlign: align});
  assert.equal(live.root.classList.contains('u2-form-captions-left'), false);
  align.value = 'left';
  assert.equal(live.root.classList.contains('u2-form-captions-left'), true);
  align.value = 'right';
  assert.equal(live.root.classList.contains('u2-form-captions-left'), false);
  live.dispose();
  align.value = 'left';
  assert.equal(live.root.classList.contains('u2-form-captions-left'), false, 'the effect died with the form');
});

ui('add(input, host) lays the row out in the host while the input still joins the form', async () => {
  const f = mount(new Form());
  const host = document.createElement('div');
  f.root.append(host);
  const name = input('name', 60);
  f.add(name, host);
  const plain = input('plain', 40);
  f.add(plain);
  assert.equal(name.root.parentNode, host);
  assert.equal(plain.root.parentNode, f.root.querySelector('.u2-form-rows'));
  assert.deepEqual(f.inputs, [name, plain], 'both joined the form');
  name.value.value = 'Aspirin';
  plain.value.value = 'x';
  assert.deepEqual(f.getValues(), {name: 'Aspirin', plain: 'x'});
  await flush();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '60px',
    'hosted labels join the measured column');
  f.dispose();
  name.dispose();
  plain.dispose();
});

ui('auto: the measure pass flips to tall when squeezed and to WIDE when roomy', async () => {
  const f = mount(new Form());
  const name = input('Compound name', 120);
  f.root.clientWidth = 200;
  f.add(name);
  await flush();
  assert.equal(tall(f), true, '200 - 120 < 150 → tall');
  assert.equal(wide(f), false);
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '120px',
    'label measuring keeps running while auto-tall');

  f.root.clientWidth = 400;
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), false, '400 - 120 > 160 → the roomy state');
  assert.equal(wide(f), true, 'which is wide, by user ruling — not the plain caption-left look');
  f.dispose();
  name.dispose();
});

ui('auto: no flip inside the 10px hysteresis band, and a hidden host keeps the current state', async () => {
  const f = mount(new Form());
  const name = input('Compound name', 120);
  f.root.clientWidth = 275; // 155 sits inside [150, 160]
  f.add(name);
  await flush();
  assert.equal(tall(f) || wide(f), false, 'the band keeps the never-measured plain state');

  f.root.clientWidth = 200;
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), true);
  f.root.clientWidth = 275;
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), true, 'the band keeps tall too — no flicker');

  f.root.clientWidth = 0;
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), true, 'clientWidth 0 (hidden host) never rewrites the state');

  f.root.clientWidth = 400;
  f._scheduleAlign();
  await flush();
  assert.equal(wide(f), true);
  f.root.clientWidth = 275;
  f._scheduleAlign();
  await flush();
  assert.equal(wide(f), true, 'the band keeps wide from the wide side too');
  f.dispose();
  name.dispose();
});

ui('explicit normal never switches; explicit tall skips label measuring entirely', async () => {
  const never = mount(new Form({layout: 'normal'}));
  const a = input('Compound name', 120);
  never.root.clientWidth = 200;
  never.add(a);
  await flush();
  assert.equal(tall(never), false);
  assert.equal(wide(never), false, 'explicit normal never goes wide either');
  never.dispose();
  a.dispose();

  const forced = mount(new Form({layout: 'tall'}));
  const b = input('Compound name', 120);
  forced.root.clientWidth = 200;
  forced.add(b);
  await flush();
  assert.equal(tall(forced), true);
  assert.equal(forced.root.style.getPropertyValue('--u2-form-label-width'), '',
    'no measure pass under explicit tall');
  forced.dispose();
  b.dispose();
});

ui('dialogs are exempt: an auto form inside .u2-dialog never flips tall', async () => {
  const dialog = document.createElement('div');
  dialog.className = 'u2-dialog';
  document.body.append(dialog);
  const f = new Form();
  dialog.append(f.root);
  const name = input('Compound name', 120);
  f.root.clientWidth = 200;
  f.add(name);
  await flush();
  assert.equal(tall(f), false, 'the platform exempts dialog forms (ui.ts:1437-1441)');
  assert.equal(wide(f), false, 'a dialog form keeps the plain aligned look — the platform\'s dialog skin');

  // a form that flipped wide BEFORE moving into a dialog loses the class on the next pass
  f.root.classList.add('u2-form-wide');
  f._scheduleAlign();
  await flush();
  assert.equal(wide(f), false, 'the exemption clears a stale wide state');
  f.dispose();
  name.dispose();
});

ui('the --dg-form-min-editor-width token drives the threshold, and the 150 fallback is never cached', async () => {
  const f = mount(new Form());
  const name = input('Compound name', 120);
  // token not served yet (detached-stylesheet boot): 300 - 120 = 180 > 160 → wide on the fallback
  f.root.clientWidth = 300;
  f.add(name);
  await flush();
  assert.equal(tall(f), false, 'the 150 fallback keeps caption-left');
  assert.equal(wide(f), true);

  // the token appears (stylesheet applied): 180 < 200 → tall — a cached fallback would ignore it
  f.root.style.setProperty('--dg-form-min-editor-width', '200px');
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), true, 'the late token is honored — the fallback was not cached');
  f.dispose();
  name.dispose();
});

ui('clientWidth is read AFTER the tall class comes off — the flip-back decision uses untall geometry', async () => {
  const f = mount(new Form());
  const name = input('Compound name', 120);
  const widths = {normal: 200, tall: 180};
  Object.defineProperty(f.root, 'clientWidth', {
    configurable: true,
    get: () => f.root.classList.contains('u2-form-tall') ? widths.tall : widths.normal,
  });
  f.add(name);
  await flush();
  assert.equal(tall(f), true, '200 - 120 < 150 → tall');

  widths.normal = 400;
  f._scheduleAlign();
  await flush();
  assert.equal(tall(f), false,
    'the pass measures with the class removed (400) — reading the tall-state 180 would stick tall');
  assert.equal(wide(f), true, 'and the roomy read lands on wide');
  f.dispose();
  name.dispose();
});

ui('folding a section (labels measuring 0) never reflows the shared label column', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const long = input('A very long caption', 88);
  const short = input('Short', 48);
  f.addAll([long, short]);
  await flush();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '88px');

  // fold: the widest label goes display:none and measures 0 — its remembered width stands in
  const label = long.root.querySelector('.u2-input-label');
  label.offsetWidth = 0;
  f._scheduleAlign();
  await flush();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '88px',
    'the column held through the fold');

  label.offsetWidth = 88;
  f._scheduleAlign();
  await flush();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '88px', 'and after expand');
  f.dispose();
  long.dispose();
  short.dispose();
});

ui('a nested form leaves the auto switching to the outer one', async () => {
  const outer = mount(new Form());
  const inner = new Form();
  outer.root.append(inner.root);
  const name = input('Compound name', 120);
  inner.root.clientWidth = 200;
  inner.add(name);
  await flush();
  assert.equal(tall(inner), false, 'the nested form never flips itself (ui.ts:1419-1424 parity)');
  assert.equal(wide(inner), false, 'exempt forms keep the plain caption-left look, never wide');
  outer.dispose();
  inner.dispose();
  name.dispose();
});

ui('dispose disconnects the ResizeObserver', () => {
  const f = new Form();
  const observer = ResizeObserver.instances.find((o) => o.targets.includes(f.root));
  assert.notEqual(observer, undefined, 'an auto form observes its root');
  f.dispose();
  assert.equal(observer.disconnected, true);

  const fixed = new Form({layout: 'normal'});
  assert.equal(ResizeObserver.instances.some((o) => o.targets.includes(fixed.root)), false,
    'a fixed layout attaches no observer');
  fixed.dispose();
});

// the shim loads no stylesheets — the level-heading cascade is pinned at the source: the l3
// sub-head style must carry the doubled class, or .u2-section-header's bold wins by load order
test('form.css: the l3 sub-head selector outranks the section-header bold', () => {
  const css = readFileSync(new URL('../css/form.css', import.meta.url), 'utf8');
  assert.match(css, /\.u2-form-category\.u2-form-category-l3 \{/,
    'l3 needs the doubled class to win the weight cascade');
  assert.match(css, /\.u2-form-category\.u2-form-category-l1 \{/,
    'l1 (h1 scale) is doubled the same way');
  assert.doesNotMatch(css, /^\.u2-form-category-l2/m,
    'l2 IS the base category style — no bare l2 block to lose the cascade with');
});

ui('auto survives a missing ResizeObserver (the headless guard)', async () => {
  const saved = globalThis.ResizeObserver;
  delete globalThis.ResizeObserver;
  try {
    const f = mount(new Form());
    const name = input('Compound name', 120);
    f.root.clientWidth = 200;
    f.add(name);
    await flush();
    assert.equal(tall(f), true, 'the rAF pass still runs without an observer');
    f.dispose();
    name.dispose();
  } finally {
    globalThis.ResizeObserver = saved;
  }
});
