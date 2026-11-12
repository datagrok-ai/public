/* The u2 ↔ bdd contract on the base classes: label/editor/error wiring and ARIA state on every
   input, the state block `DG.U2.state(el)` reads, the element → control door across portaled
   popups, and the registry walk (every visual tag stamps `data-u2`, every input example is named).
   Contract: core/docs/features/ui2/CONTRACT.md §1–§2. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {state, status, control, items, revealItem, props, setProp} from 'datagrok-api/u2core';
import {Scope, Control, Overlay, TextInput, NumberInput, BoolInput, ChoiceInput, registry, registerAll,
  renderSpec, SpecContext} from '../src/index.js';

registerAll();

/** Input tags whose registry example builds an unnamed input (no label, aria-label or placeholder). */
const UNNAMED_INPUT_EXAMPLES = [];

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

function mounted(input) {
  document.body.append(input.root);
  return input;
}

const INPUTS = [
  ['text-input', () => new TextInput({label: 'Name', value: 'Ada'}), 'input'],
  ['number-input', () => new NumberInput({label: 'Age', value: 42}), 'input'],
  ['bool-input', () => new BoolInput({label: 'Active', value: true}), 'input'],
  ['choice-input', () => new ChoiceInput({label: 'Stage', items: ['Discovery', 'Phase I'], value: 'Phase I'}), 'select'],
];

for (const [kind, make, editorTag] of INPUTS) {
  smoke(`${kind}: the label names the editor, ids derive from the root id`, async () => {
    const input = mounted(make());
    const root = input.root;
    const labelEl = root.querySelector('[data-u2-part="label"]');
    const target = root.querySelector(editorTag);
    assert.match(root.id, /^u2-\d+$/);
    assert.equal(labelEl.id, `${root.id}-label`);
    assert.equal(target.id, `${root.id}-editor`);
    assert.equal(labelEl.getAttribute('for'), target.id);
    assert.equal(root.querySelector('[data-u2-part="error"]').id, `${root.id}-error`);
    input.dispose();
  });

  smoke(`${kind}: state carries kind, label, value, text, validity and enabled`, async () => {
    const input = mounted(make());
    const s = state(input.root);
    assert.equal(s.kind, kind);
    assert.equal(s.label, input.label);
    assert.equal(s.enabled, true);
    assert.equal(s.visible, true);
    assert.equal(s.focused, false);
    assert.equal(s.valid, true);
    assert.equal(s.error, null);
    assert.deepEqual(s.value, Control.jsonValue(input.peek()));
    assert.equal(typeof s.text, 'string');
    assert.deepEqual(JSON.parse(JSON.stringify(s)), s, 'JSON-safe');
    assert.deepEqual(status(input.root).state, s, 'the status carries the same block');
    input.dispose();
  });
}

smoke('text: value as displayed; choice: the caption; bool: true/false', async () => {
  const text = mounted(new TextInput({label: 'Name', value: 'Ada'}));
  assert.equal(state(text.root).text, 'Ada');
  text.value = 'Grace';
  assert.equal(state(text.root).text, 'Grace');
  const choice = mounted(new ChoiceInput({label: 'Stage', items: [{value: 'p1', label: 'Phase I'}], value: 'p1'}));
  assert.equal(state(choice.root).value, 'p1');
  assert.equal(state(choice.root).text, 'Phase I');
  const bool = mounted(new BoolInput({label: 'Active', value: true}));
  assert.equal(state(bool.root).value, true);
  assert.equal(state(bool.root).text, 'true');
  text.dispose();
  choice.dispose();
  bool.dispose();
});

smoke('invalid: aria-invalid and aria-describedby point at the error the user reads', async () => {
  const input = mounted(new TextInput({label: 'Code', value: ''}));
  input.addValidator((v) => v ? null : 'Code is required');
  const editor = input.root.querySelector('input');
  const error = input.root.querySelector('[data-u2-part="error"]');
  assert.equal(state(input.root).valid, false, 'the verdict is there from the start');
  assert.equal(editor.hasAttribute('aria-invalid'), false, 'and shows once the input is touched');
  assert.equal(state(input.root).error, null);
  input.touch();
  assert.equal(editor.getAttribute('aria-invalid'), 'true');
  assert.equal(editor.getAttribute('aria-describedby'), error.id);
  assert.equal(error.textContent, 'Code is required');
  assert.equal(state(input.root).valid, false);
  assert.equal(state(input.root).error, 'Code is required');

  input.value = 'X1';
  assert.equal(editor.hasAttribute('aria-invalid'), false);
  assert.equal(editor.hasAttribute('aria-describedby'), false);
  assert.equal(state(input.root).valid, true);
  assert.equal(state(input.root).error, null);
  input.dispose();
});

smoke('disabled: aria-disabled on the root, state.enabled false', async () => {
  const input = mounted(new TextInput({label: 'Name'}));
  input.disabled = 'Not editable in this test';
  assert.equal(input.root.getAttribute('aria-disabled'), 'true');
  assert.equal(input.root.querySelector('input').disabled, true);
  assert.equal(state(input.root).enabled, false);
  input.disabled = null;
  assert.equal(input.root.hasAttribute('aria-disabled'), false);
  assert.equal(state(input.root).enabled, true);
  input.dispose();
});

smoke('visible = false sets hidden; state.visible follows hidden ancestors and detachment', async () => {
  const input = mounted(new TextInput({label: 'Name'}));
  input.visible = false;
  assert.equal(input.root.hasAttribute('hidden'), true);
  assert.equal(input.visible, false);
  assert.equal(state(input.root).visible, false);
  input.visible = true;
  assert.equal(state(input.root).visible, true);

  const host = document.createElement('div');
  document.body.append(host);
  host.append(input.root);
  host.hidden = true;
  assert.equal(state(input.root).visible, false, 'a hidden ancestor hides it');
  host.hidden = false;
  input.root.remove();
  assert.equal(state(input.root).visible, false, 'detached');
  input.dispose();
});

smoke('focused: focus on the editor or inside a popup the control owns', async () => {
  const input = mounted(new TextInput({label: 'Name'}));
  input.root.querySelector('input').focus();
  assert.equal(state(input.root).focused, true);

  const content = document.createElement('div');
  const option = document.createElement('button');
  content.append(option);
  Overlay.show(input.box, content, {owner: input.scope});
  assert.equal(content.dataset.u2OwnerId, input.rootId);
  option.focus();
  assert.equal(state(input.root).focused, true, 'focus in the owned popup');
  input.dispose();
});

smoke('inline: no label part; the label text or the explicit name becomes the aria-label', async () => {
  const labelled = mounted(new TextInput({inline: true, label: 'Filter'}));
  assert.equal(labelled.root.querySelector('[data-u2-part="label"]'), null);
  assert.equal(labelled.root.querySelector('input').getAttribute('aria-label'), 'Filter');
  assert.equal(state(labelled.root).label, 'Filter');

  const named = mounted(new TextInput({inline: true, name: 'query'}));
  assert.equal(named.root.querySelector('input').getAttribute('aria-label'), 'query');

  const placeholder = mounted(new TextInput({inline: true, placeholder: 'Search…'}));
  assert.equal(state(placeholder.root).label, 'Search…', 'the placeholder names it last');
  assert.equal(state(placeholder.root).placeholder, 'Search…');
  labelled.dispose();
  named.dispose();
  placeholder.dispose();
});

smoke('forElement resolves a node in a portaled popup to the popup owner', async () => {
  const input = mounted(new TextInput({label: 'Name'}));
  const content = document.createElement('div');
  const row = document.createElement('span');
  content.append(row);
  Overlay.show(input.box, content, {owner: input.scope});
  assert.equal(content.parentElement.classList.contains('u2-overlay'), true, 'portaled to the body host');
  assert.equal(Control.forElement(row), input);
  assert.equal(control(row), input);
  assert.equal(state(row).kind, 'text-input');

  const unowned = document.createElement('div');
  document.body.append(unowned);
  assert.equal(Control.forElement(unowned), undefined);
  assert.equal(state(unowned), null, 'null outside any u2 control');
  input.dispose();
});

smoke('onChanged runs owned by the builder, not by the effect run, and gets the input', async () => {
  const owner = new Scope();
  const calls = [];
  let opened;
  const input = owner.build(() => new TextInput({label: 'Name', onChanged: (v, source) => {
    calls.push([v, source]);
    opened = Scope.ambient;
  }}));
  input.value = 'a';
  assert.deepEqual(calls, [['a', input]]);
  assert.equal(opened, owner, 'what onChanged builds belongs to the builder');
  owner.dispose();
});

smoke('items: the default reveal searches rendered items by key, then label; ambiguity names the keys', async () => {
  const list = new Control();
  list.root.dataset.u2 = 'sample-list';
  for (const [key, text] of [['a', 'Alpha'], ['b', 'Beta'], ['c', 'beta']]) {
    const item = document.createElement('div');
    item.dataset.u2Item = key;
    item.textContent = ` ${text} `;
    list.root.append(item);
  }
  document.body.append(list.root);
  assert.equal(items(list.root), null, 'a plain control lists no items');
  assert.equal((await revealItem(list.root, 'a')).dataset.u2Item, 'a');
  assert.equal((await revealItem(list.root, 'alpha')).dataset.u2Item, 'a');
  assert.equal(await revealItem(list.root, 'Gamma'), null);
  await assert.rejects(revealItem(list.root, 'BETA'), {message: 'u2: 2 items match "BETA": b, c'});

  const popup = document.createElement('div');
  popup.dataset.u2OwnerId = list.rootId;
  const option = document.createElement('div');
  option.dataset.u2Item = 'd';
  option.textContent = 'Delta';
  popup.append(option);
  document.body.append(popup);
  assert.equal(await revealItem(list.root, 'Delta'), option, 'items in an owned popup');
  list.dispose();
});

smoke('props / setProp: declared properties, read-only and unknown ones throw', async () => {
  class Sample extends Control {
    constructor() {
      super();
      this.fixed = 'construction-only';
      this._title = 'untitled';
      this.root.dataset.u2 = 'sample';
    }

    get title() {
      return this._title;
    }

    set title(x) {
      this._title = x;
    }
  }
  const sample = new Sample();
  sample.componentMeta = {props: [{name: 'title', type: 'string'}, {name: 'fixed', type: 'string'}]};
  document.body.append(sample.root);
  assert.deepEqual(props(sample.root), [
    {name: 'title', type: 'string', value: 'untitled', writable: true},
    {name: 'fixed', type: 'string', value: 'construction-only', writable: false},
  ]);
  setProp(sample.root, 'title', 'renamed');
  assert.equal(sample.title, 'renamed');
  assert.throws(() => setProp(sample.root, 'fixed', 'x'), {message: 'u2: sample.fixed is read-only'});
  assert.throws(() => setProp(sample.root, 'nope', 1), {message: 'u2: sample has no property nope'});
  sample.dispose();
});

test('jsonValue: the projection behind state.value', () => {
  assert.equal(Control.jsonValue(new Date(0)), '1970-01-01T00:00:00.000Z');
  assert.equal(Control.jsonValue({name: 'demog', rowCount: 5}), 'demog');
  assert.deepEqual(Control.jsonValue([1, 'a', null, undefined, () => 1]), [1, 'a', null, null, null]);
  assert.deepEqual(Control.jsonValue({a: 1, b: {c: [true]}, f: () => 1}), {a: 1, b: {c: [true]}});
  assert.deepEqual(Control.jsonValue(new Set(['x', 'y'])), ['x', 'y']);
  assert.equal(Control.jsonValue(NaN), 'NaN');
  assert.equal(Control.jsonValue(10n), '10');
  assert.equal(Control.jsonValue(() => 1), undefined);
  assert.equal(Control.jsonValue(document.createElement('div')), undefined);
  assert.equal(Control.jsonValue(null), null);
});

smoke('registry walk: every visual tag stamps data-u2; every input example has an accessible name', async () => {
  const unnamed = [];
  for (const meta of registry.metas()) {
    if (meta.visual === false)
      continue;
    const instance = renderSpec({$schema: 'dg-ui/1', root: meta.example}, new SpecContext(), registry);
    document.body.append(instance.root);
    const root = instance.root.children[0];
    assert.notEqual(root?.dataset.u2, undefined, `${meta.tag}: its example root carries no data-u2`);
    if (meta.category === 'Inputs') {
      const s = state(root);
      assert.equal(s.kind, root.dataset.u2, `${meta.tag}: state answers its root`);
      if (!s.label)
        unnamed.push(meta.tag);
    }
    instance.dispose();
  }
  assert.deepEqual(unnamed.sort(), UNNAMED_INPUT_EXAMPLES, 'input examples without an accessible name');
});
