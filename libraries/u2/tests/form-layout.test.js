/* Form layouts: the layout enum's class emission, the captionAlign live toggle, the
   Form.add(input, host) seam, and the auto machinery — _pickLayout's platform truth table
   (threshold + 10px hysteresis, js-api ui.ts:1443-1448) and the class toggle over stubbed
   widths. Layout metrics are dom-shim fakes; nothing here asserts pixels. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {Form, LabelColumn} from '../src/components/forms/form.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {ChoiceInput} from '../src/components/inputs/choice-input.js';

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

ui('alignment lands before the next task: a form built and attached in one go is aligned after a microtask', async () => {
  const f = mount(new Form({layout: 'normal'}));
  f.addAll([input('Name', 60), input('Natural analog', 97)]);
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '', 'not while the building code runs');
  await Promise.resolve();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '97px');
  for (const i of f.inputs.slice())
    i.dispose();
  f.dispose();
});

ui('labelWidth: a fixed column — nothing measured, the auto layout decides with it', async () => {
  const f = mount(new Form({labelWidth: 28}));
  f.add(input('R1', 20));
  f.root.clientWidth = 400;
  f.layoutNow();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '28px');
  assert.equal(wide(f), true);
  f.inputs[0].dispose();
  f.dispose();
});

ui('labelColumn: the forms sharing a column take its widest caption, a hidden one\'s included', async () => {
  const column = new LabelColumn();
  const a = mount(new Form({layout: 'normal', labelColumn: column}));
  const b = mount(new Form({layout: 'normal', labelColumn: column}));
  a.add(input('Name', 60));
  b.add(input('Natural analog', 97));
  await flush();
  assert.equal(column.width.value, 97);
  assert.equal(a.root.style.getPropertyValue('--u2-form-label-width'), '97px', 'the narrower form follows');
  assert.equal(b.root.style.getPropertyValue('--u2-form-label-width'), '97px');
  b.inputs[0].dispose();
  b.dispose();
  assert.equal(column.width.value, 60, 'a disposed form leaves the column');
  assert.equal(a.root.style.getPropertyValue('--u2-form-label-width'), '60px');
  a.inputs[0].dispose();
  a.dispose();
});

ui('error: the form\'s first shown message — none on a fresh form, validity has it from the start', async () => {
  const f = mount(new Form());
  const name = new TextInput({label: 'Name', required: true});
  f.add(name);
  assert.notEqual(f.validity.value, null);
  assert.equal(f.error.value, null);
  f.validate();
  assert.equal(f.error.value, f.validity.value);
  name.dispose();
  f.dispose();
});

ui('width: an input\'s own editor width beats the widths a form fits', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const cap = new ChoiceInput({label: 'R1', items: ['H', 'OH'], width: 110});
  const name = new TextInput({label: 'Name', width: '12em'});
  f.addAll([cap, name]);
  await flush();
  assert.equal(cap.root.querySelector('select').style.width, '110px');
  assert.equal(name.root.querySelector('input').style.width, '12em');
  cap.dispose();
  name.dispose();
  f.dispose();
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

  // a fixed caption-left layout observes too: built hidden, it aligns its labels once it shows
  const fixed = new Form({layout: 'normal'});
  const fixedObserver = ResizeObserver.instances.find((o) => o.targets.includes(fixed.root));
  assert.notEqual(fixedObserver, undefined, 'a fixed caption-left layout observes its root');
  fixed.dispose();
  assert.equal(fixedObserver.disconnected, true);

  const tall = new Form({layout: 'tall'});
  assert.equal(ResizeObserver.instances.some((o) => o.targets.includes(tall.root)), false,
    'a tall form has no label column to align');
  tall.dispose();
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

ui('validity skips hidden and disabled inputs (P79); validate() shows every counted verdict', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const hidden = new TextInput({label: 'Code', name: 'code', required: true});
  const off = new TextInput({label: 'Note', name: 'note', required: true, enabled: false});
  f.addAll([hidden, off]);
  hidden.visible = false;
  assert.equal(f.validity.peek(), null);
  assert.equal(f.validate(), true);
  hidden.visible = true;
  assert.notEqual(f.validity.peek(), null, 'shown again, it fails the form');
  assert.equal(f.validate(), false);
  f.dispose();
  hidden.dispose();
  off.dispose();
});

ui('choices take the platform width: the longest option plus 30px, measured into --u2-choice-width', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const role = new ChoiceInput({label: 'Role', name: 'role', items: ['Viewer', 'Administrator'], value: 'Viewer'});
  f.add(role);
  const probes = [];
  const append = f.root.append.bind(f.root);
  f.root.append = (el) => {
    probes.push(el);
    el.getBoundingClientRect = () => ({width: el.textContent.length * 7});
    append(el);
  };
  await flush();
  assert.equal(role.root.style.getPropertyValue('--u2-choice-width'), `${'Administrator'.length * 7 + 30}px`);
  assert.equal(probes.length, 1, 'one probe per pass, removed after');
  assert.equal(probes[0].parentNode, null);
  f.dispose();
  role.dispose();
});

test('form.css: the auto form keeps ui.form widths in its wide state; the error line stays inside the row', () => {
  const form = readFileSync(new URL('../css/form.css', import.meta.url), 'utf8');
  assert.match(form, /\.u2-form-auto\.u2-form-wide \.u2-input-root:not\(\.u2-input-fill\)\[data-u2='choice-input'\] select\.u2-input-editor/);
  assert.match(form, /width: var\(--u2-choice-width, 200px\);/);
  assert.match(form, /\.u2-form-auto\.u2-form-wide \.u2-input-root:not\(\.u2-input-fill\) \.u2-input-editor\.u2-number/);
  const inputs = readFileSync(new URL('../css/inputs.css', import.meta.url), 'utf8');
  const at = inputs.indexOf('\n.u2-input-error {');
  const error = inputs.slice(at, inputs.indexOf('}', at));
  assert.match(error, /box-sizing: border-box;/, 'the form indents it by padding, inside its 100% basis');
});

ui('Form.fitChoices: a dialog body that lays inputs out as a form measures its choices the same way', async () => {
  const host = document.createElement('div');
  document.body.append(host);
  const stage = new ChoiceInput({label: 'Stage', name: 'stage', items: ['I', 'Preclinical'], value: 'I'});
  const name = new TextInput({label: 'Name', name: 'name'});
  host.append(stage.root, name.root);
  const append = host.append.bind(host);
  host.append = (el) => {
    el.getBoundingClientRect = () => ({width: el.textContent.length * 7});
    append(el);
  };
  Form.fitChoices(host, [stage.root, name.root]);
  assert.equal(stage.root.style.getPropertyValue('--u2-choice-width'), `${'Preclinical'.length * 7 + 30}px`);
  assert.equal(name.root.style.getPropertyValue('--u2-choice-width'), '', 'only choices are measured');
  stage.dispose();
  name.dispose();
});

// the shim cannot lay out: the platform's form widths are pinned at the source
test('form.css / choice.css / inputs.css: the platform form widths hold in every host', () => {
  const form = readFileSync(new URL('../css/form.css', import.meta.url), 'utf8');
  const choice = readFileSync(new URL('../css/choice.css', import.meta.url), 'utf8');
  const inputs = readFileSync(new URL('../css/inputs.css', import.meta.url), 'utf8');
  const block = (sheet, selector) => {
    const at = sheet.indexOf(selector);
    assert.notEqual(at, -1, `${selector} is a rule`);
    return sheet.slice(sheet.indexOf('{', at) + 1, sheet.indexOf('}', at));
  };
  const grow = ':not(.u2-multi-choice):not(.u2-radio):not(.u2-image)';
  assert.match(block(form, `.u2-form .u2-input-box > div.u2-input-editor${grow}`), /flex-grow: 1;/,
    'a normal-state form (dialog-hosted, nested) grows its div editors, as ui.form does');
  assert.ok(form.includes(`.u2-dialog-content > .u2-input-root > .u2-input-box > div.u2-input-editor${grow}`),
    'and so does a u2 Dialog body');
  assert.match(block(form, '.ui-form .ui-input-root > .u2-input-box:has(> input.u2-input-editor)'), /flex: 1;/,
    'a box bridged into a platform form or dialog grows like the platform text editor');
  assert.match(block(form, '.u2-form-wide .u2-input-root[data-u2=\'multi-choice-input\'] .u2-input-editor'), /flex: none;/,
    'the wide form\'s list exception keeps no 0 basis');
  assert.match(block(form, '.u2-form:not(.u2-form-tall) .u2-input-error'),
    /padding-left: calc\(var\(--u2-form-label-width, 100px\) \+ var\(--dg-space-m\)\);/, 'the verdict starts under the editor');
  assert.match(block(choice, ':where([data-u2=\'choice-input\'] > .u2-input-box) > select.u2-input-editor'),
    /width: var\(--u2-choice-width, 200px\);/, 'the measured choice width travels with the input');
  assert.match(block(choice, 'div.u2-input-editor.u2-multi-choice'), /min-width: 150px;/);
  assert.match(block(inputs, '\n.u2-input-editor.u2-columns {'), /min-width: 100px;/);
});

ui('the auto layout marks itself, so the platform widths hold in its wide state', async () => {
  const auto = new Form();
  assert.equal(auto.root.classList.contains('u2-form-auto'), true);
  const wide = new Form({layout: 'wide'});
  assert.equal(wide.root.classList.contains('u2-form-auto'), false);
  auto.dispose();
  wide.dispose();
});

ui('addGroup: one caption over several inputs — the group caption sizes the column, its inputs do not', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const name = input('Name', 40);
  const line = input('Line', 30);
  const background = input('A long background caption', 160);
  f.add(name).addGroup('Colors', [line, background], {captionsAfter: true});
  f.root.querySelector('.u2-form-group-label').offsetWidth = 52;
  f.layoutNow();
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '52px', 'the group caption, not its inputs\'');
  assert.deepEqual(f.inputs.map((i) => i.name), ['Name', 'Line', 'A long background caption']);
  const group = f.root.querySelector('[role="group"]');
  assert.equal(document.getElementById(group.getAttribute('aria-labelledby')).textContent, 'Colors');
  assert.equal(line.root.classList.contains('u2-input-caption-after'), true);
  for (const i of [name, line, background])
    i.dispose();
  f.dispose();
});

ui('addColumns: the inputs in one grid, at the minimum and column count asked, still the form\'s own', async () => {
  const f = mount(new Form({layout: 'normal'}));
  const symbol = input('Symbol', 44);
  const name = input('Name', 38);
  const polymer = input('Polymer', 52);
  f.addColumns([symbol, name, polymer], {minWidth: 176});
  f.layoutNow();
  const grid = f.root.querySelector('.u2-form-columns');
  assert.deepEqual([...grid.children], [symbol.root, name.root, polymer.root]);
  assert.equal(grid.style.getPropertyValue('--u2-form-column-min'), '176px');
  assert.equal(grid.style.getPropertyValue('--u2-form-columns'), '', 'two columns by default');
  assert.equal(f.root.style.getPropertyValue('--u2-form-label-width'), '52px', 'its captions size the one label column');
  assert.deepEqual(f.inputs.map((i) => i.name), ['Symbol', 'Name', 'Polymer']);
  const host = document.createElement('div');
  f.root.append(host);
  const three = [input('A', 10), input('B', 10), input('C', 10)];
  f.addColumns(three, {columns: 3, host});
  assert.equal(host.querySelector('.u2-form-columns').style.getPropertyValue('--u2-form-columns'), '3');
  for (const i of [symbol, name, polymer, ...three])
    i.dispose();
  f.dispose();
});

test('addColumns: a track never goes under the minimum, so a narrow row holds fewer columns rather than squeezing them', () => {
  const sheet = readFileSync(new URL('../css/form.css', import.meta.url), 'utf8');
  const at = sheet.indexOf('.u2-form-columns {');
  const rule = sheet.slice(at, sheet.indexOf('}', at));
  assert.match(rule, /repeat\(auto-fill, minmax\(min\(100%, max\(var\(--u2-form-column-min, 200px\)/);
  assert.match(rule, /\/ var\(--u2-form-columns, 2\)\)\), 1fr\)\)/);
  assert.match(sheet, /\.u2-form\.u2-form \.u2-form-columns > \.u2-input-root:not\(\.u2-input-fill\) \.u2-input-editor \{\s+flex-shrink: 1;/,
    'a fixed-width editor (a number\'s 150px) gives way to its column');
});

ui('quiet: the form says so, and an input says while it is empty', async () => {
  const f = mount(new Form({quiet: true}));
  const name = new TextInput({label: 'Name', name: 'name'});
  f.add(name);
  assert.equal(f.root.classList.contains('u2-quiet'), true);
  assert.equal(name.root.classList.contains('u2-input-empty'), true, 'nothing typed: it keeps its underline');
  name.value.value = 'Cysteine';
  assert.equal(name.root.classList.contains('u2-input-empty'), false);
  name.value.value = '';
  assert.equal(name.root.classList.contains('u2-input-empty'), true);
  const sheet = readFileSync(new URL('../css/inputs.css', import.meta.url), 'utf8');
  assert.match(sheet, /\.u2-quiet \{\s+--dg-input-border: var\(--dg-border-width\) solid transparent;/);
  assert.match(sheet, /\.u2-quiet \.u2-input-root\.u2-input-empty,/);
  name.dispose();
  f.dispose();
});

ui('addButtons: the row goes where `host` says, at the start when asked, and keeps the form-button look there', async () => {
  const f = mount(new Form());
  const footer = document.createElement('div');
  document.body.append(footer);
  const save = document.createElement('button');
  f.addButtons((row) => row.append(save), {align: 'start', host: footer});
  const row = footer.firstElementChild ?? footer.children[0];
  assert.equal(row.classList.contains('u2-form-buttons'), true);
  assert.equal(row.classList.contains('u2-form-buttons-start'), true);
  assert.equal(row.contains(save), true);
  assert.equal(f.root.contains(row), false, 'placed in the host, not under the rows');
  const sheet = readFileSync(new URL('../css/form.css', import.meta.url), 'utf8');
  assert.match(sheet, /\.u2-form-buttons \.u2-btn-primary:not\(:disabled\)[^{]*\{\s+background: var\(--dg-blue-2\);/,
    'a primary button in the row is the platform form button outside the form too');
  f.dispose();
});
