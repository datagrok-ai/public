/* Platform pickers (src/dg/inputs/pickers.ts) and the column renderer (src/dg/entities/column-renderer.ts) over
   the shell stub in tests/platform-stub.mjs: the live item lists, the rich column picker, and the
   renderer's glyph/semType/unknown-name rows. The tables are the getter-backed frames of
   tests/platform-doubles.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {DataFrame} from './platform-doubles.mjs';

register('./platform-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {Input} = await import('../src/core/input-base.js');
const {columnInput, tableInput, tablesInput, ColumnPicker} = await import('../src/dg/inputs/pickers.js');
const {columnRenderer} = await import('../src/dg/entities/column-renderer.js');

function picker(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      grok.resetShell();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const DEMO = () => new DataFrame([
  {name: 'age', type: 'int'}, {name: 'height', type: 'double'},
  {name: 'sex', type: 'string', semType: 'Sex'}, {name: 'started', type: 'datetime'},
]);

function options(input) {
  return input.root.querySelectorAll('option').map((o) => o.value).filter((v) => v !== '');
}

/** A local file as the import action reads it: a name plus the two readers it may reach for. */
function localFile(name, text = 'a,b\n1,2') {
  return {name, text: async () => text, arrayBuffer: async () => new Uint8Array([1, 2, 3]).buffer};
}

/** Picks `file` through the import action's own hidden picker, as the OS dialog would. */
async function importFile(input, file) {
  const picker = input.root.querySelector('.u2-input-options input');
  picker.files = [file];
  fire(picker, 'change');
  await flush();
}

picker('columnInput: names follow the table, a dropped column clears the value', () => {
  const table = DEMO();
  const input = columnInput('Column', table);
  assert.deepEqual(options(input), ['age', 'height', 'sex', 'started']);
  input.value.value = 'sex';

  table.columns.remove('sex');
  assert.deepEqual(options(input), ['age', 'height', 'started']);
  assert.equal(input.value.value, null);

  input.dispose();
  assert.equal(table.liveSubscriptions(), 0, 'disposal takes the subscription with it');
});

picker('columnInput: a renamed column keeps the value and re-lists under the new name', () => {
  const table = DEMO();
  const input = columnInput('Column', table);
  input.value.value = 'sex';

  table.columns.byName('sex').name = 'gender';
  assert.deepEqual(options(input), ['age', 'height', 'gender', 'started']);
  assert.equal(input.value.value, 'gender', 'the value follows the column, it is not dropped');

  table.columns.byName('age').name = 'ageYears';
  assert.equal(input.value.value, 'gender', 'renaming another column leaves the value alone');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0, 'both subscriptions go with the input');
});

picker('columnInput: the filter decides what is offered', () => {
  const input = columnInput('Number', DEMO(), {filter: (c) => c.isNumerical});
  assert.deepEqual(options(input), ['age', 'height', 'started'], 'datetime is numerical to the platform');
  input.dispose();
});

picker('tableInput: the item list follows the open tables', () => {
  grok.openTable('demog');
  const input = tableInput('Table');
  assert.deepEqual(options(input), ['demog']);

  grok.openTable('cars');
  assert.deepEqual(options(input), ['demog', 'cars']);
  input.value.value = 'cars';

  grok.closeTable('cars');
  assert.deepEqual(options(input), ['demog']);
  assert.equal(input.value.value, null, 'a closed table takes the value with it');

  input.dispose();
  assert.equal(grok.liveSubscriptions(), 0);
});

picker('tablesInput: a checkbox per open table, the value is their names', () => {
  grok.openTable('demog');
  grok.openTable('cars');
  const input = tablesInput('Tables');
  const boxes = () => input.root.querySelectorAll('input[type=checkbox]');
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars']);

  const box = boxes()[1];
  box.checked = true;
  fire(box, 'change');
  assert.deepEqual(input.value.value, ['cars']);

  grok.openTable('geo');
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars', 'geo']);
  assert.deepEqual(input.value.value, ['cars'], 'opening a table keeps what was checked');
  assert.equal(boxes()[1].checked, true);

  grok.closeTable('cars');
  assert.deepEqual(input.value.value, [], 'a closed table leaves the value with it');

  input.dispose();
  assert.equal(grok.liveSubscriptions(), 0);
});

picker('tableInput: the import action opens a file, adds the table and picks it', async () => {
  const input = tableInput('Table');
  let systemic = null;
  input.effect(() => {
    input.value.value;
    systemic = Input.isSystemWrite;
  });

  await importFile(input, localFile('demog.csv'));
  assert.deepEqual(grok.shell.tableNames, ['demog'], 'the file joins the workspace, named after it');
  assert.deepEqual(options(input), ['demog']);
  assert.equal(input.value.value, 'demog');
  assert.equal(systemic, false, 'the imported table is picked the way a user picks one');

  grok.closeTable('demog');
  assert.equal(input.value.value, null);
  assert.equal(systemic, true, 'while the prune a close leaves behind is the input speaking');
  input.dispose();
});

picker('tableInput: a d42 blob comes in through the binary reader', async () => {
  const input = tableInput('Table');
  await importFile(input, localFile('keys.d42'));
  assert.deepEqual(options(input), ['keys']);
  assert.equal(input.value.value, 'keys');
  input.dispose();
});

picker('tableInput: an extension neither reader takes is refused, and the input carries on',
  async () => {
    const errors = [];
    grok.shell.error = (m) => errors.push(m);
    grok.openTable('demog');
    const input = tableInput('Table');
    input.value.value = 'demog';

    await importFile(input, localFile('report.xlsx'));
    assert.deepEqual(errors, ['File extension .xlsx is not supported.']);
    assert.deepEqual(options(input), ['demog'], 'nothing was added');
    assert.equal(input.value.value, 'demog', 'and nothing was taken away');

    await importFile(input, localFile('cars.csv'));
    assert.equal(input.value.value, 'cars');
    input.dispose();
  });

picker('tableInput: a file the reader chokes on is reported, not thrown', async () => {
  const errors = [];
  grok.shell.error = (m) => errors.push(m);
  grok.data.parseCsv = () => {
    throw new Error('unexpected end of input');
  };

  const input = tableInput('Table');
  await importFile(input, localFile('broken.csv'));
  assert.deepEqual(errors, ['broken.csv could not be opened.']);
  assert.equal(input.value.value, null);
  input.dispose();
});

picker('tablesInput: an imported table joins the checks', async () => {
  grok.openTable('demog');
  const input = tablesInput('Tables');
  const boxes = () => input.root.querySelectorAll('input[type=checkbox]');

  await importFile(input, localFile('cars.csv'));
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars']);
  assert.deepEqual(input.value.value, ['cars'], 'the imported table comes in checked');
  assert.equal(boxes()[1].checked, true);

  await importFile(input, localFile('geo.csv'));
  assert.deepEqual(input.value.value, ['cars', 'geo'], 'and joins what was checked already');
  input.dispose();
});

picker('rich columnInput: suggestions filter as typed and a pick becomes the value', async () => {
  const table = DEMO();
  const input = columnInput('Column', table, {rich: true});
  assert.ok(input instanceof ColumnPicker);
  assert.equal(input.root.dataset.u2, 'column-picker');
  document.body.append(input.root);

  const box = input.root.querySelector('input');
  box.focus();
  fire(box, 'keydown', {key: 'ArrowDown'});
  await flush();
  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.deepEqual(popup.querySelectorAll('.u2-column-name').map((el) => el.textContent),
    ['age', 'height', 'sex', 'started']);

  box.value = 'he';
  fire(box, 'input');
  await flush();
  assert.deepEqual(popup.querySelectorAll('.u2-column-name').map((el) => el.textContent), ['height']);

  fire(popup.querySelector('.u2-typeahead-option'), 'pointerdown');
  assert.equal(input.value.value, 'height');
  assert.equal(box.value, 'height');

  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

picker('rich columnInput: two-way with the value signal; a dropped column clears it', async () => {
  const table = DEMO();
  const input = new ColumnPicker({label: 'Column', table, value: 'sex'});
  const box = input.root.querySelector('input');
  assert.equal(box.value, 'sex', 'the constructed value shows in the box');

  input.value.value = 'age';
  assert.equal(box.value, 'age');

  table.columns.remove('age');
  assert.equal(input.value.value, null);
  assert.equal(box.value, '');
  input.dispose();
});

picker('rich columnInput: a renamed column carries the value and the box text', () => {
  const table = DEMO();
  const input = new ColumnPicker({label: 'Column', table, value: 'sex'});
  const box = input.root.querySelector('input');

  table.columns.byName('sex').name = 'gender';
  assert.equal(input.value.value, 'gender');
  assert.equal(box.value, 'gender');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

picker('tablesInput: an empty shell reads as empty, not as a broken box', () => {
  const input = tablesInput('Tables');
  assert.equal(input.root.querySelector('.u2-multi-choice-empty').textContent, 'No open tables');

  grok.openTable('demog');
  assert.equal(input.root.querySelector('.u2-multi-choice-empty'), null);
  grok.closeTable('demog');
  assert.equal(input.root.querySelector('.u2-multi-choice-empty').textContent, 'No open tables');
  input.dispose();
});

picker('columnRenderer: glyph per type, semType hint, and a name the table lost', () => {
  const renderer = columnRenderer(DEMO());
  const glyph = (name) => renderer.icon(name).className;
  assert.ok(glyph('age').includes('fa-hashtag'));
  assert.ok(glyph('sex').includes('fa-font'));
  assert.ok(glyph('started').includes('fa-calendar-alt'));
  assert.ok(glyph('gone').includes('fa-question-circle'), 'an unknown name still renders');

  const row = renderer.listItem('sex');
  assert.equal(row.querySelector('.u2-column-name').textContent, 'sex');
  assert.equal(row.querySelector('.u2-column-semtype').textContent, 'Sex');
  assert.equal(renderer.listItem('age').querySelector('.u2-column-semtype'), null,
    'no semType, no hint');
  assert.equal(renderer.tooltip('sex'), 'sex · string · Sex');
  assert.equal(renderer.tooltip('gone'), 'gone');
  assert.equal(renderer.caption('gone'), 'gone');
});
