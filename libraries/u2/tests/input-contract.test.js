/* One input contract on the Input base (WO-C7): visible, readOnly, required (touch-gated), focus(),
   onEnter/onEscape that keep the key to themselves, fill, undefined as empty (P74), commitOn on the
   numeric editors, the `change` event (R18), the core inputs' state and items, the form skipping
   hidden and disabled inputs (P79) — and no definite-assignment field in a class whose
   createEditor runs inside the base constructor (P18). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {state} from 'datagrok-api/u2core';
import {signal, Scope, Input, TextInput, TextArea, NumberInput, BoolInput, ChoiceInput, MultiChoiceInput,
  RadioInput, SliderInput, QNumInput, Form, EMPTY, registry, registerAll, renderSpec, SpecContext, each} from '../src/index.js';
import {BigIntInput} from '../src/components/inputs/bigint-input.js';

registerAll();

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

const editorOf = (input) => input.root.querySelector('input, textarea, select');

smoke('visible: the option hides the root, live, and state.visible follows', async () => {
  const shown = signal(false);
  const input = mount(new TextInput({label: 'Name', visible: shown}));
  assert.equal(input.root.hidden, true);
  assert.equal(input.visible, false);
  assert.equal(state(input.root).visible, false);
  shown.value = true;
  assert.equal(input.root.hidden, false);
  assert.equal(state(input.root).visible, true);
  input.visible = false;
  assert.equal(input.root.hidden, true);
  input.dispose();
});

smoke('readOnly: native readOnly and aria-readonly on a text editor; the class; live; state', async () => {
  const ro = signal(true);
  const input = mount(new TextInput({label: 'Name', value: 'Ada', readOnly: ro}));
  const editor = editorOf(input);
  assert.equal(editor.readOnly, true);
  assert.equal(editor.getAttribute('aria-readonly'), 'true');
  assert.equal(input.root.classList.contains('u2-input-readonly'), true);
  assert.equal(state(input.root).readOnly, true);
  ro.value = false;
  assert.equal(editor.readOnly, false);
  assert.equal(editor.hasAttribute('aria-readonly'), false);
  assert.equal(input.readOnly, false);
  assert.equal(state(input.root).readOnly, false);
  input.dispose();
});

smoke('readOnly: a checkbox, a switch and a number stepper ignore the user', async () => {
  const box = mount(new BoolInput({label: 'Active', value: false, readOnly: true}));
  fire(box.root.querySelector('input'), 'click');
  box.root.querySelector('input').click();
  assert.equal(box.value.peek(), false, 'the box does not toggle');
  const sw = mount(new BoolInput({label: 'On', value: false, switch: true, readOnly: true}));
  fire(sw.root.querySelector('[role="switch"]'), 'click');
  fire(sw.root.querySelector('[role="switch"]'), 'keydown', {key: ' '});
  fire(sw.root.querySelector('[role="switch"]'), 'keydown', {key: 'Enter'});
  assert.equal(sw.value.peek(), false, 'the switch does not toggle');
  const n = mount(new NumberInput({label: 'Count', value: 5, readOnly: true}));
  const field = n.root.querySelector('.u2-number-box input');
  assert.equal(field.readOnly, true);
  const step = fire(field, 'keydown', {key: 'ArrowUp'});
  assert.equal(step, false, 'the stepping key is cancelled');
  assert.equal(n.value.peek(), 5);
  box.dispose();
  sw.dispose();
  n.dispose();
});

smoke('required: aria-required; an empty value is invalid, shown only once touched', async () => {
  const input = mount(new TextInput({label: 'Name', required: true}));
  const editor = editorOf(input);
  assert.equal(editor.getAttribute('aria-required'), 'true');
  assert.equal(input.required, true);
  assert.equal(input.validity.peek(), EMPTY, 'the verdict is there');
  assert.equal(input.root.querySelector('[data-u2-part="error"]').textContent, '', 'but not shown yet');
  assert.equal(editor.hasAttribute('aria-invalid'), false);
  let s = state(input.root);
  assert.equal(s.valid, false);
  assert.equal(s.error, null, 'state.error is what the user sees');
  assert.equal(s.required, true);
  fire(editor, 'focusout');
  assert.equal(input.root.querySelector('[data-u2-part="error"]').textContent, EMPTY, 'leaving it shows it');
  assert.equal(editor.getAttribute('aria-invalid'), 'true');
  s = state(input.root);
  assert.equal(s.error, EMPTY);
  input.value.value = 'Ada';
  assert.equal(input.validity.peek(), null);
  assert.equal(Input.isEmpty([]), true);
  assert.equal(Input.isEmpty(undefined), true);
  assert.equal(Input.isEmpty(0), false);
  input.dispose();
});

smoke('required: an added validator speaks first; the required verdict only when it has nothing to say', async () => {
  const input = mount(new TextInput({label: 'Code', required: true, value: 'x'}));
  input.addValidator((v) => v === 'x' ? 'Not x' : null);
  input.touch();
  assert.equal(input.validity.peek(), 'Not x');
  input.value.value = '';
  assert.equal(input.validity.peek(), EMPTY);
  input.dispose();
});

smoke('nullable: false is required: a cleared NumberInput fails after touch', async () => {
  const input = mount(new NumberInput({label: 'Age', value: 3, nullable: false}));
  const field = input.root.querySelector('.u2-number-box input');
  field.value = '';
  fire(field, 'input');
  assert.equal(input.value.peek(), null);
  assert.equal(input.validity.peek(), EMPTY);
  assert.equal(input.root.querySelector('[data-u2-part="error"]').textContent, EMPTY, 'an edit touches it');
  input.dispose();
});

smoke('focus() lands on the focus target; a multi choice focuses its first box', async () => {
  const input = mount(new TextInput({label: 'Name'}));
  input.focus();
  assert.equal(document.activeElement, editorOf(input));
  const multi = mount(new MultiChoiceInput({label: 'Tags', items: ['a', 'b']}));
  multi.focus();
  assert.equal(document.activeElement, multi.root.querySelector('input[type="checkbox"]'));
  input.dispose();
  multi.dispose();
});

smoke('onEnter/onEscape: called with the event and the input, the key kept from the document', async () => {
  const seen = [];
  const onDocument = (e) => seen.push(`document ${e.key}`);
  document.addEventListener('keydown', onDocument);
  try {
    const input = mount(new TextInput({label: 'Name', value: 'Ada',
      onEnter: (e, it) => seen.push(`enter ${it.value.peek()} ${e.defaultPrevented}`),
      onEscape: (e, it) => seen.push(`escape ${it === input}`)}));
    const n = mount(new NumberInput({label: 'Count', value: 1, onEnter: (_, it) => seen.push(`number ${it.value.peek()}`)}));
    fire(editorOf(input), 'keydown', {key: 'Enter'});
    fire(editorOf(input), 'keydown', {key: 'Escape'});
    fire(editorOf(input), 'keydown', {key: 'a'});
    const field = n.root.querySelector('.u2-number-box input');
    field.value = '7';
    fire(field, 'keydown', {key: 'Enter'});
    assert.deepEqual(seen, ['enter Ada true', 'escape true', 'document a', 'number 7']);
    input.dispose();
    n.dispose();
  } finally {
    document.removeEventListener('keydown', onDocument);
  }
});

smoke('fill: the class the skin stretches; TextArea rows', async () => {
  const input = mount(new TextInput({label: 'Name', fill: true}));
  assert.equal(input.root.classList.contains('u2-input-fill'), true);
  const area = mount(new TextArea({label: 'Notes', rows: 5}));
  assert.equal(area.root.querySelector('textarea').rows, 5);
  const css = readFileSync(new URL('../css/inputs.css', import.meta.url), 'utf8');
  assert.match(css, /\.u2-input-root\.u2-input-fill > \.u2-input-box > \.u2-input-editor/);
  input.dispose();
  area.dispose();
});

smoke('P74: undefined is empty in every core editor', async () => {
  const none = () => signal(undefined);
  const number = mount(new NumberInput({label: 'N', bind: none(), min: 0, max: 10, slider: true}));
  assert.equal(number.root.querySelector('.u2-number-box input').value, '');
  assert.equal(state(number.root).value, null);
  assert.equal(number.validity.peek(), null);
  const big = mount(new BigIntInput({label: 'B', bind: none()}));
  assert.equal(editorOf(big).value, '');
  const q = mount(new QNumInput({label: 'Q', bind: none()}));
  assert.equal(editorOf(q).value, '');
  const slider = mount(new SliderInput({label: 'S', bind: none()}));
  assert.equal(slider.text, '');
  const sw = mount(new BoolInput({label: 'On', switch: true, bind: none()}));
  assert.equal(sw.root.querySelector('[role="switch"]').getAttribute('aria-checked'), 'false');
  const multi = mount(new MultiChoiceInput({label: 'M', items: ['a'], bind: none()}));
  assert.equal(state(multi.root).text, '');
  for (const x of [number, big, q, slider, sw, multi])
    x.dispose();
});

smoke('commitOn change: a NumberInput writes on blur or Enter only, and shows the text meanwhile', async () => {
  const input = mount(new NumberInput({label: 'N', value: 1, commitOn: 'change'}));
  const field = input.root.querySelector('.u2-number-box input');
  field.value = '42';
  fire(field, 'input');
  assert.equal(input.value.peek(), 1, 'typing does not write');
  assert.equal(field.value, '42');
  fire(field, 'blur');
  assert.equal(input.value.peek(), 42, 'leaving writes');
  field.value = '43';
  fire(field, 'input');
  fire(field, 'keydown', {key: 'Enter'});
  assert.equal(input.value.peek(), 43, 'Enter writes');
  input.dispose();
});

smoke('commitOn input: a NumberInput keeps the text being typed ("1." stays)', async () => {
  const input = mount(new NumberInput({label: 'N'}));
  const field = input.root.querySelector('.u2-number-box input');
  field.value = '1.';
  fire(field, 'input');
  await flush();
  assert.equal(input.value.peek(), 1);
  assert.equal(field.value, '1.', 'the echo of the edit never rewrites the field');
  input.dispose();
});

smoke('change event: on(\'change\') gets the new value and the input; the initial value is no change', async () => {
  const input = mount(new TextInput({label: 'Name', value: 'a'}));
  const got = [];
  input.on('change', (x, it) => {
    if (typeof x === 'string')
      got.push([x, it === input]);
  });
  input.value.value = 'b';
  assert.deepEqual(got, [['b', true]]);
  input.dispose();
});

smoke('state: bool checked; choice, radio and multi choice items; multi choice text is the checked captions', async () => {
  const bool = mount(new BoolInput({label: 'Active', value: true}));
  assert.equal(state(bool.root).checked, true);
  const choice = mount(new ChoiceInput({label: 'Stage', items: [{value: 'p1', label: 'Phase I'}, 'Discovery'],
    value: 'p1'}));
  let s = state(choice.root);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.label, i.selected]),
    [['p1', 'Phase I', true], ['Discovery', 'Discovery', false]]);
  assert.deepEqual(s.selection, ['p1']);
  assert.equal(s.open, undefined, 'a native select reports no open state');
  const radio = mount(new RadioInput({label: 'Mode', items: ['fast', 'slow'], value: 'slow'}));
  s = state(radio.root);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.checked]), [['fast', false], ['slow', true]]);
  assert.equal(s.text, 'slow');
  const row = radio.root.querySelector('[data-u2-item="slow"]');
  assert.equal(row.querySelector('[data-u2-part="item-label"]').textContent, 'slow');
  assert.equal(row.querySelector('input').dataset.value, 'slow');
  assert.equal(await radio.revealItem('fast'), radio.root.querySelector('[data-u2-item="fast"]'));
  const multi = mount(new MultiChoiceInput({label: 'Targets', items: [{value: 'gpcr', label: 'GPCR'},
    {value: 'kin', label: 'Kinase'}, {value: 'ion', label: 'Ion channel'}], value: ['kin', 'ion']}));
  s = state(multi.root);
  assert.equal(s.text, 'Kinase, Ion channel');
  assert.deepEqual(s.selection, ['kin', 'ion']);
  assert.deepEqual(s.items.list.map((i) => i.checked), [false, true, true]);
  const item = multi.root.querySelector('[data-u2-item="kin"]');
  assert.equal(item.dataset.value, 'kin');
  assert.equal(item.querySelector('[data-u2-part="item-label"]').textContent, 'Kinase');
  assert.equal(await multi.revealItem('Ion channel'), multi.root.querySelector('[data-u2-item="ion"]'));
  for (const x of [bool, choice, radio, multi])
    x.dispose();
});

smoke('form: hidden and disabled invalid inputs do not fail it; validate() touches and focuses the first', async () => {
  const form = mount(new Form({layout: 'normal'}));
  const hidden = new TextInput({label: 'Hidden', name: 'hidden', required: true});
  const disabled = new TextInput({label: 'Disabled', name: 'disabled', required: true});
  const name = new TextInput({label: 'Name', name: 'name', required: true});
  form.addAll([hidden, disabled, name]);
  hidden.visible = false;
  disabled.enabled = false;
  assert.equal(form.validity.peek(), EMPTY, 'the visible empty field fails the form');
  assert.equal(name.root.querySelector('[data-u2-part="error"]').textContent, '', 'untouched, nothing shown');
  name.value.value = 'Ada';
  assert.equal(form.validity.peek(), null, 'the hidden and the disabled ones do not');
  assert.equal(form.validate(), true);
  name.value.value = '';
  assert.equal(form.validate(), false);
  assert.equal(name.root.querySelector('[data-u2-part="error"]').textContent, EMPTY, 'validate() touched it');
  assert.equal(document.activeElement, editorOf(name));
  assert.equal(hidden.root.querySelector('[data-u2-part="error"]').textContent, '', 'a hidden input stays untouched');
  const s = state(form.root);
  assert.equal(s.valid, false);
  assert.equal(s.error, EMPTY);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.label, i.invalid ?? false, i.disabled ?? false]),
    [['hidden', 'Hidden', true, false], ['disabled', 'Disabled', true, true], ['name', 'Name', true, false]]);
  hidden.visible = true;
  name.value.value = 'Ada';
  assert.equal(form.validity.peek(), EMPTY, 'shown again, it counts');
  // an ancestor hidden by other means is seen on the next read
  const section = document.createElement('div');
  form.root.append(section);
  section.append(hidden.root);
  section.hidden = true;
  assert.equal(form.validate(), true);
  form.dispose();
  for (const x of [hidden, disabled, name])
    x.dispose();
});

smoke('form: an input built in a keyed row leaves the form with its row; remove() takes one out', async () => {
  const form = mount(new Form({layout: 'normal'}));
  const rows = document.createElement('div');
  form.addElement(rows);
  const keys = signal(['R1', 'R2']);
  const owner = new Scope();
  Scope.runWith(owner, () => each(rows, keys, (k) => k, (k) => {
    const row = document.createElement('div');
    form.add(new TextInput({label: k, name: k, required: true}), row);
    return row;
  }));
  assert.deepEqual(form.inputs.map((i) => i.name), ['R1', 'R2']);
  assert.equal(form.validity.peek(), EMPTY);
  keys.value = ['R1'];
  assert.deepEqual(form.inputs.map((i) => i.name), ['R1'], 'the dropped row took its input out');
  form.inputs[0].value.value = 'H';
  assert.equal(form.validity.peek(), null, 'the gone input no longer fails the form');
  const extra = new TextInput({label: 'Extra', name: 'extra', required: true});
  form.add(extra);
  assert.equal(form.validity.peek(), EMPTY);
  form.remove(extra);
  assert.equal(extra.root.isConnected, false);
  assert.deepEqual(form.inputs.map((i) => i.name), ['R1']);
  assert.equal(form.validity.peek(), null);
  owner.dispose();
  assert.deepEqual(form.inputs, []);
  extra.dispose();
  form.dispose();
});

smoke('reset: a fresh start for another record — no waiting verdict, no unreadable text, no onChanged', async () => {
  const changes = [];
  const name = mount(new TextInput({label: 'Name', required: true, onChanged: (v) => changes.push(v)}));
  name.touch();
  assert.notEqual(name.root.querySelector('[data-u2-part="error"]').textContent, '', 'touched and empty: shown');
  name.reset('');
  assert.equal(name.root.querySelector('[data-u2-part="error"]').textContent, '', 'the new record has not been touched');
  name.reset('Ada');
  assert.equal(name.value.value, 'Ada');
  assert.deepEqual(changes, [], 'a reset is no edit');
  const id = mount(new NumberInput({label: 'ID', value: 0}));
  const editor = id.root.querySelector('input');
  editor.value = 'abc';
  fire(editor, 'input');
  assert.equal(id.validity.value, 'Not a number');
  id.reset(0);
  assert.equal(id.validity.value, null, 'the same value, its text read again');
  assert.equal(editor.value, '0');
  name.dispose();
  id.dispose();
});

smoke('addOptions: a caller puts its own icon on an input\'s options rail', async () => {
  const input = mount(new ChoiceInput({label: 'R1', items: ['H', 'OH']}));
  const pencil = document.createElement('button');
  input.addOptions(pencil);
  assert.equal(pencil.parentElement, input.root.querySelector('[data-u2-part="options"]'));
  input.dispose();
});

smoke('spec: readOnly and required are registry props every input tag takes', async () => {
  const tags = registry.metas().filter((m) => m.category === 'Inputs' && m.visual !== false);
  for (const tag of ['u2-text-input', 'u2-number-input', 'u2-bool-input', 'u2-choice-input', 'u2-radio-input'])
    assert.ok(tags.find((m) => m.tag === tag).props.some((p) => p.name === 'readOnly'), `${tag} has readOnly`);
  const ctx = new SpecContext();
  const built = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-text-input',
    props: {label: 'Name', readOnly: true, required: true}}}, ctx, registry);
  const root = built.root.querySelector('[data-u2="text-input"]');
  const s = state(root);
  assert.equal(s.readOnly, true);
  assert.equal(s.required, true);
  built.dispose();
});

/** Classes whose createEditor still declares `!:` fields — their owners convert them and empty this
 * list (C10: the popup and picker inputs; C6: chips, access grid; C13: from-dart-input). */
const P18_PENDING = [];

test('P18: a class with createEditor declares no `!:` field (the base constructor runs it first)', () => {
  const src = new URL('../src/', import.meta.url);
  const offenders = [];
  for (const dir of ['components/inputs', 'components/forms', 'dg/inputs']) {
    for (const name of readdirSync(new URL(dir + '/', src))) {
      if (!name.endsWith('.ts') || name.endsWith('.d.ts'))
        continue;
      const path = `${dir}/${name}`;
      const text = readFileSync(new URL(path, src), 'utf8');
      if (!text.includes('createEditor(') || P18_PENDING.includes(path))
        continue;
      const fields = text.split('\n').filter((line) =>
        /^\s+(private |protected |public )?(readonly )?\w+!:/.test(line) && !line.includes('declare '));
      if (fields.length > 0)
        offenders.push(`${path}: ${fields.map((f) => f.trim()).join(' | ')}`);
    }
  }
  assert.deepEqual(offenders, [], 'use `private declare _x: T;` for fields createEditor assigns');
});
