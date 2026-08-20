/* PropertyEditor — the metadata editor's visibility matrix, its record value model and its
   re-targeting. The control mirrors `DG.Property.propertyOptions` rather than reading it, so the
   whole suite runs off-platform, like the object-form one. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {propertyEditor} from '../src/dg/property-editor.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function editorTest(name, body) {
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

function labels(editor) {
  const found = editor.root.querySelectorAll('.u2-input-label');
  return Array.from(found).map((l) => l.textContent);
}

function fieldOf(editor, label, selector = 'input') {
  const rows = editor.root.querySelectorAll('.u2-input-root');
  const row = Array.from(rows).find((r) => r.querySelector('.u2-input-label')?.textContent === label);
  return row.querySelector(selector);
}

function typeSelect(editor) {
  return fieldOf(editor, 'Type', 'select');
}

editorTest('shows the identity fields for every property, the numerical ones for a number', () => {
  const editor = mount(propertyEditor({name: 'dose', type: 'double', min: 0, max: 1000}));
  assert.deepEqual(labels(editor), ['Name', 'Type', 'Friendly name', 'Description', 'Nullable',
    'Semantic type', 'Input type', 'Editor', 'Units', 'Format', 'Category',
    'Min', 'Max', 'Step', 'Show slider', 'Show plus/minus']);
  assert.equal(fieldOf(editor, 'Name').value, 'dose');
  assert.equal(fieldOf(editor, 'Min').value, '0');
  assert.equal(fieldOf(editor, 'Description', 'textarea') != null, true,
    'the description carries the platform\'s own TextArea editor hint');
  assert.deepEqual(Array.from(typeSelect(editor).children).map((o) => o.value),
    ['string', 'int', 'double', 'num', 'bigint', 'qnum', 'bool', 'datetime', 'list', 'map', 'file']);
  editor.dispose();
});

editorTest('switching the type swaps the type-dependent fields', () => {
  const records = [];
  const editor = mount(propertyEditor({name: 'dose', type: 'double', min: 0, max: 1000},
    {onChanged: (options) => records.push(options)}));

  const type = typeSelect(editor);
  type.value = 'string';
  fire(type, 'change');
  const shown = labels(editor);
  assert.equal(shown.includes('Min'), false, 'min is gone');
  assert.equal(shown.includes('Max'), false);
  assert.equal(shown.includes('Choices'), true, 'and choices took its place');

  type.value = 'bool';
  fire(type, 'change');
  assert.deepEqual(labels(editor).slice(11), [], 'a bool asks for nothing extra');

  type.value = 'int';
  fire(type, 'change');
  assert.deepEqual(labels(editor).slice(11), ['Min', 'Max', 'Step', 'Show slider', 'Show plus/minus']);
  assert.equal(fieldOf(editor, 'Min').value, '0', 'the value the record still carries is shown again');
  editor.dispose();
});

editorTest('every edit reports the whole record, and the options signal carries it', () => {
  const records = [];
  const editor = mount(propertyEditor({name: 'dose', type: 'double', min: 0, units: 'mg'},
    {onChanged: (options) => records.push(options)}));
  assert.deepEqual(records, [], 'building the editor is not an edit');
  assert.deepEqual(editor.options.value, {name: 'dose', type: 'double', min: 0, units: 'mg'});

  const max = fieldOf(editor, 'Max');
  max.value = '5';
  fire(max, 'input');
  assert.deepEqual(records, [{name: 'dose', type: 'double', min: 0, units: 'mg', max: 5}],
    'the untouched fields come along');
  assert.deepEqual(editor.options.value, records[0]);

  const name = fieldOf(editor, 'Name');
  name.value = 'weight';
  fire(name, 'change');
  assert.equal(records.length, 2);
  assert.deepEqual(records[1], {name: 'weight', type: 'double', min: 0, units: 'mg', max: 5});
  assert.notEqual(records[1], records[0], 'a fresh record every time, so a reader sees the change');
  editor.dispose();
});

editorTest('the name is reported on commit, not per keystroke', () => {
  const records = [];
  const editor = mount(propertyEditor({name: 'dose', type: 'double'},
    {onChanged: (options) => records.push(options)}));
  const name = fieldOf(editor, 'Name');
  // a caller renames the property it keys by this name, so 'dose' → 'name' must not rename
  // through 'n', 'na' and 'nam' on the way
  for (const typed of ['n', 'na', 'nam', 'name']) {
    name.value = typed;
    fire(name, 'input');
  }
  assert.deepEqual(records, [], 'typing is not an edit');
  assert.equal(editor.options.value.name, 'dose', 'and the record still carries the old name');

  fire(name, 'change');
  assert.deepEqual(records, [{name: 'name', type: 'double'}], 'blur or Enter reports it once');

  // every other field keeps reporting as it is typed — the sync counters of the A/B demo count
  // those keystrokes
  const units = fieldOf(editor, 'Units');
  units.value = 'mg';
  fire(units, 'input');
  assert.equal(records.length, 2);
  assert.equal(records[1].units, 'mg');
  editor.dispose();
});

editorTest('choices edit as a comma-separated list', () => {
  const records = [];
  const editor = mount(propertyEditor({name: 'series', type: 'string'},
    {onChanged: (options) => records.push(options)}));
  const choices = fieldOf(editor, 'Choices');
  choices.value = 'a,"b,c"';
  fire(choices, 'input');
  assert.deepEqual(records[records.length - 1].choices, ['a', 'b,c']);
  editor.dispose();
});

editorTest('the caller\'s record is never touched — the editor works on a copy', () => {
  const source = {name: 'dose', type: 'double'};
  const editor = mount(propertyEditor(source));
  const name = fieldOf(editor, 'Name');
  name.value = 'weight';
  fire(name, 'change');
  assert.deepEqual(source, {name: 'dose', type: 'double'});
  assert.equal(editor.options.value.name, 'weight');
  editor.dispose();
});

editorTest('a caller\'s validator reports on the field, and survives a re-target', async () => {
  const taken = ['weight'];
  const editor = mount(propertyEditor({name: 'dose', type: 'double'},
    {validators: {name: (x) => taken.includes(x) ? 'Name is already taken' : null}}));
  const error = () => editor.root.querySelector('.u2-input-error').textContent;
  assert.equal(error(), '');

  const name = fieldOf(editor, 'Name');
  name.value = 'weight';
  fire(name, 'change');
  assert.equal(error(), 'Name is already taken');

  editor.setTarget({name: 'weight', type: 'int'});
  await flush();
  assert.equal(error(), 'Name is already taken', 'the rebuilt field validates the same way');
  editor.dispose();
});

editorTest('re-targeting shows the new record and leaks nothing', async () => {
  const changes = [];
  const editor = mount(propertyEditor({name: 'dose', type: 'double'},
    {onChanged: (options) => changes.push(options)}));
  const live = [];
  for (let i = 0; i < 10; i++) {
    editor.setTarget({name: `p${i}`, type: i % 2 ? 'string' : 'int', min: i});
    assert.equal(editor.options.value.name, `p${i}`, 'the record changes at once');
    await flush();                                    // the fields are rebuilt on a microtask
    assert.equal(fieldOf(editor, 'Name').value, `p${i}`);
    assert.equal(labels(editor).includes('Choices'), i % 2 === 1);
    // the two shapes carry a different number of fields; each must cost the same every time
    live.push(Scope.liveCount);
    assert.equal(live[i], live[i % 2], 'one generation at a time');
  }
  assert.deepEqual(changes, [], 're-targeting is not an edit');
  editor.dispose();
});

editorTest('re-targeting from an onChanged handler rebuilds once, on the last target', async () => {
  const seen = [];
  let editor = null;
  // the handler runs inside the effect the current generation owns, so rebuilding there would
  // dispose the effect that is running
  editor = mount(propertyEditor({name: 'dose', type: 'double'}, {onChanged: (options) => {
    seen.push(options);
    editor.setTarget({name: 'series', type: 'string'});
    editor.setTarget({name: 'weight', type: 'int', min: 1});
  }}));
  const name = fieldOf(editor, 'Name');
  name.value = 'x';
  fire(name, 'change');
  await flush();

  assert.deepEqual(seen, [{name: 'x', type: 'double'}], 'the edit was reported, and nothing threw');
  assert.equal(fieldOf(editor, 'Name').value, 'weight', 'the last target won');
  assert.deepEqual(labels(editor).slice(11), ['Min', 'Max', 'Step', 'Show slider', 'Show plus/minus']);
  assert.equal(fieldOf(editor, 'Min').value, '1');
  const settled = Scope.liveCount;

  const max = fieldOf(editor, 'Max');
  max.value = '5';
  fire(max, 'input');
  await flush();
  assert.equal(seen.length, 2, 'the rebuilt fields edit as the first ones did');
  assert.equal(seen[1].max, 5);
  assert.equal(Scope.liveCount, settled, 'and the round left no generation behind');
  editor.dispose();
});
