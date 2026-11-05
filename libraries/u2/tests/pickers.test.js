/* Platform pickers (src/dg/inputs/pickers.ts) and the column renderer (src/dg/entities/column-renderer.ts) over
   the shell stub in tests/platform-stub.mjs: the live item lists, the rich column picker, and the
   renderer's glyph/semType/unknown-name rows. The tables are the getter-backed frames of
   tests/platform-doubles.mjs. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Usage} from '../src/core/scope.js';
import {DataFrame, strayTables} from './platform-doubles.mjs';

register('./platform-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {Input} = await import('../src/core/input-base.js');
const {columnInput, tableInput, tablesInput, ColumnPicker, TableInput} = await import('../src/dg/inputs/pickers.js');
const {ColumnInput} = await import('../src/dg/inputs/column-combo.js');
const {columnsInput, columnsMapInput, aggregatedColumnsInput} = await import('../src/dg/inputs/columns.js');
const {columnRenderer} = await import('../src/dg/entities/column-renderer.js');

/** The names of the tables a test opened: the shell's workspace starts with its stray ones (a nameless among them). */
const opened = () => grok.shell.tableNames.slice(strayTables().length);

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

// `~height.cache` is hidden, as Chem's `~<column>.canonicalSmiles` is: the platform lists neither
const DEMO = () => new DataFrame([
  {name: 'age', type: 'int'}, {name: 'height', type: 'double'}, {name: '~height.cache', type: 'double'},
  {name: 'sex', type: 'string', semType: 'Sex'}, {name: 'started', type: 'datetime'},
]);

function options(input) {
  return [...input.root.querySelectorAll('option')].map((o) => o.value).filter((v) => v !== '');
}

/** A local file as the import action reads it: a name plus the two readers it may reach for. */
function localFile(name, text = 'a,b\n1,2') {
  return {name, text: async () => text, arrayBuffer: async () => new Uint8Array([1, 2, 3]).buffer};
}

/** Picks `file` through the import action's own hidden picker, as the OS dialog would. Under the
 * action API the picker is transient (FP-P-2) — the open-file icon mints it — so when the rail
 * carries no input yet, the drive goes through the icon first. */
async function importFile(input, file) {
  let picker = input.root.querySelector('input[type=file]');
  if (picker === null) {
    fire(input.root.querySelector('.u2-input-options [aria-label="Open file"]'), 'click');
    await flush();
    picker = input.root.querySelector('input[type=file]') ?? document.querySelector('input[type=file]');
  }
  picker.files = [file];
  fire(picker, 'change');
  await flush();
}

picker('columnInput: names follow the table, a dropped column clears the value', () => {
  const table = DEMO();
  const input = columnInput('Column', table);
  assert.deepEqual(options(input), ['age', 'height', 'sex', 'started']);
  input.value = 'sex';

  table.columns.remove('sex');
  assert.deepEqual(options(input), ['age', 'height', 'started']);
  assert.equal(input.value, null);

  input.dispose();
  assert.equal(table.liveSubscriptions(), 0, 'disposal takes the subscription with it');
});

picker('columnInput: a renamed column keeps the value and re-lists under the new name', () => {
  const table = DEMO();
  const input = columnInput('Column', table);
  input.value = 'sex';

  table.columns.byName('sex').name = 'gender';
  assert.deepEqual(options(input), ['age', 'height', 'gender', 'started']);
  assert.equal(input.value, 'gender', 'the value follows the column, it is not dropped');

  table.columns.byName('age').name = 'ageYears';
  assert.equal(input.value, 'gender', 'renaming another column leaves the value alone');
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
  input.value = 'cars';

  grok.closeTable('cars');
  assert.deepEqual(options(input), ['demog']);
  assert.equal(input.value, null, 'a closed table takes the value with it');

  input.dispose();
  assert.equal(grok.liveSubscriptions(), 0);
});

picker('tableInput, tablesInput: a table the shell holds without a name is not offered', () => {
  grok.openTable('demog');
  grok.openTable('nameless').name = undefined;
  const input = tableInput('Table');
  const many = tablesInput('Tables');
  assert.deepEqual(options(input), ['demog']);
  assert.deepEqual([...many.root.querySelectorAll('input[type=checkbox]')].map((b) => b.value), ['demog']);

  grok.openTable('cars');
  assert.deepEqual(options(input), ['demog', 'cars']);

  input.dispose();
  many.dispose();
});

picker('tablesInput: a checkbox per open table, the value is their names', () => {
  grok.openTable('demog');
  grok.openTable('cars');
  const input = tablesInput('Tables');
  const boxes = () => [...input.root.querySelectorAll('input[type=checkbox]')];
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars']);

  const box = boxes()[1];
  box.checked = true;
  fire(box, 'change');
  assert.deepEqual(input.value, ['cars']);

  grok.openTable('geo');
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars', 'geo']);
  assert.deepEqual(input.value, ['cars'], 'opening a table keeps what was checked');
  assert.equal(boxes()[1].checked, true);

  grok.closeTable('cars');
  assert.deepEqual(input.value, [], 'a closed table leaves the value with it');

  input.dispose();
  assert.equal(grok.liveSubscriptions(), 0);
});

picker('tableInput: the import action opens a file, adds the table and picks it', async () => {
  const input = tableInput('Table');
  let systemic = null;
  input.effect(() => {
    input.value;
    systemic = Input.isSystemWrite;
  });

  await importFile(input, localFile('demog.csv'));
  assert.deepEqual(opened(), ['demog'], 'the file joins the workspace, named after it');
  assert.deepEqual(options(input), ['demog']);
  assert.equal(input.value, 'demog');
  assert.equal(systemic, false, 'the imported table is picked the way a user picks one');

  grok.closeTable('demog');
  assert.equal(input.value, null);
  assert.equal(systemic, true, 'while the prune a close leaves behind is the input speaking');
  input.dispose();
});

picker('tableInput: a d42 blob comes in through the binary reader', async () => {
  const input = tableInput('Table');
  await importFile(input, localFile('keys.d42'));
  assert.deepEqual(options(input), ['keys']);
  assert.equal(input.value, 'keys');
  input.dispose();
});

picker('tableInput: an extension neither reader takes is refused, and the input carries on',
  async () => {
    const errors = [];
    grok.shell.error = (m) => errors.push(m);
    grok.openTable('demog');
    const input = tableInput('Table');
    input.value = 'demog';

    await importFile(input, localFile('report.xlsx'));
    assert.deepEqual(errors, ['File extension .xlsx is not supported.']);
    assert.deepEqual(options(input), ['demog'], 'nothing was added');
    assert.equal(input.value, 'demog', 'and nothing was taken away');

    await importFile(input, localFile('cars.csv'));
    assert.equal(input.value, 'cars');
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
  assert.equal(input.value, null);
  input.dispose();
});

picker('tablesInput: an imported table joins the checks', async () => {
  grok.openTable('demog');
  const input = tablesInput('Tables');
  const boxes = () => [...input.root.querySelectorAll('input[type=checkbox]')];

  await importFile(input, localFile('cars.csv'));
  assert.deepEqual(boxes().map((b) => b.value), ['demog', 'cars']);
  assert.deepEqual(input.value, ['cars'], 'the imported table comes in checked');
  assert.equal(boxes()[1].checked, true);

  await importFile(input, localFile('geo.csv'));
  assert.deepEqual(input.value, ['cars', 'geo'], 'and joins what was checked already');
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
  assert.deepEqual([...popup.querySelectorAll('.u2-column-name')].map((el) => el.textContent),
    ['age', 'height', 'sex', 'started']);

  box.value = 'he';
  fire(box, 'input');
  await flush();
  assert.deepEqual([...popup.querySelectorAll('.u2-column-name')].map((el) => el.textContent), ['height']);

  fire(popup.querySelector('.u2-typeahead-option'), 'pointerdown');
  assert.equal(input.value, 'height');
  assert.equal(box.value, 'height');

  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

picker('rich columnInput: two-way with the value signal; a dropped column clears it', async () => {
  const table = DEMO();
  const input = new ColumnPicker({label: 'Column', table, value: 'sex'});
  const box = input.root.querySelector('input');
  assert.equal(box.value, 'sex', 'the constructed value shows in the box');

  input.value = 'age';
  assert.equal(box.value, 'age');

  table.columns.remove('age');
  assert.equal(input.value, null);
  assert.equal(box.value, '');
  input.dispose();
});

picker('rich columnInput: a renamed column carries the value and the box text', () => {
  const table = DEMO();
  const input = new ColumnPicker({label: 'Column', table, value: 'sex'});
  const box = input.root.querySelector('input');

  table.columns.byName('sex').name = 'gender';
  assert.equal(input.value, 'gender');
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

// --- the TableInput action rail (plan-pickers FP-P-2), written blind to the implementation ---

function railIcons(input) {
  return [...input.root.querySelectorAll('.u2-input-options [data-u2="icon-button"]')];
}

/** Installs the fake platform doors for the body, sync or async, and removes them after. */
function withGlobals(globals, body) {
  Object.assign(globalThis, globals);
  const clean = () => {
    for (const name of Object.keys(globals))
      delete globalThis[name];
  };
  let result;
  try {
    result = body();
  } catch (e) {
    clean();
    throw e;
  }
  if (result instanceof Promise)
    return result.finally(clean);
  clean();
  return result;
}

picker('tableInput: the headless rail is the open-file action alone', () => {
  const input = tableInput('Table');
  assert.ok(input instanceof TableInput, 'the factory narrows to the class');
  assert.deepEqual(input.actions.map((a) => a.id), ['open-file'],
    'no platform doors, no platform icons — every existing DOM stays');
  const icons = railIcons(input);
  assert.equal(icons.length, 1);
  assert.ok(icons[0].querySelector('.fa-folder-open'));
  assert.equal(icons[0].getAttribute('aria-label'), 'Open file');
  input.dispose();
});

picker('tableInput: the platform doors put their icons on, in the platform order', () => {
  withGlobals({grok_UI_PickTableFromFiles: async () => null, grok_UI_PickTableFromQuery: async () => null},
    () => {
      const input = tableInput('Table');
      assert.deepEqual(input.actions.map((a) => a.id), ['open-file', 'add-from-files', 'query-db']);
      const icons = railIcons(input);
      assert.deepEqual(icons.map((el) => el.getAttribute('aria-label')),
        ['Open file', 'Add file from Files', 'Query database'], 'the Dart tooltips, verbatim (P-P-1)');
      assert.ok(icons[1].querySelector('.fa-folder-tree'));
      assert.ok(icons[2].querySelector('.fa-database'));
      input.dispose();
    });
});

picker('tableInput: each door gates its own icon', () => {
  withGlobals({grok_UI_PickTableFromFiles: async () => null}, () => {
    const input = tableInput('Table');
    assert.deepEqual(input.actions.map((a) => a.id), ['open-file', 'add-from-files']);
    input.dispose();
  });
  withGlobals({grok_UI_PickTableFromQuery: async () => null}, () => {
    const input = tableInput('Table');
    assert.deepEqual(input.actions.map((a) => a.id), ['open-file', 'query-db']);
    input.dispose();
  });
});

picker('tableInput: options.actions replaces the defaults wholesale at construction', () => {
  const input = tableInput('Table',
    {actions: [{id: 'only', icon: 'bolt', tooltip: 'Only', run: () => {}}]});
  assert.deepEqual(input.actions.map((a) => a.id), ['only']);
  const icons = railIcons(input);
  assert.equal(icons.length, 1);
  assert.equal(icons[0].getAttribute('aria-label'), 'Only');
  assert.equal(input.root.querySelector('.fa-folder-open'), null, 'no open-file default rides along');
  input.dispose();
});

picker('TableInput.actions: the setter re-renders the rail; the getter answers a copy', () => {
  const input = tableInput('Table');
  const runs = [];
  input.actions = [{id: 'mine', icon: 'star', tooltip: 'Mine', run: (i) => runs.push(i)}];
  assert.deepEqual(railIcons(input).map((el) => el.getAttribute('aria-label')), ['Mine'],
    'a set replaces the rail whole');
  assert.equal(input.root.querySelector('.fa-folder-open'), null);
  assert.ok(railIcons(input)[0].querySelector('.fa-star'));

  fire(railIcons(input)[0], 'click');
  assert.deepEqual(runs, [input], 'run receives the input itself');

  input.actions.push({id: 'sneaky', icon: 'bug', tooltip: 'Nope', run: () => {}});
  assert.deepEqual(input.actions.map((a) => a.id), ['mine'], 'the getter hands out a copy');
  assert.equal(railIcons(input).length, 1);

  input.actions = [...input.actions, {id: 'more', icon: 'plus', tooltip: 'More', run: () => {}}];
  assert.deepEqual(railIcons(input).map((el) => el.getAttribute('aria-label')), ['Mine', 'More'],
    'read-modify-write is the composition point');

  input.actions = [];
  assert.equal(railIcons(input).length, 0, 'and an empty set clears the rail');
  input.dispose();
});

picker('TableInput.actions: two inputs keep their own sets', () => {
  const a = tableInput('A');
  const b = tableInput('B');
  b.actions = [{id: 'b-only', icon: 'star', tooltip: 'B only', run: () => {}}];
  assert.deepEqual(a.actions.map((x) => x.id), ['open-file'], 'the sibling is untouched');
  assert.equal(a.root.querySelector('.fa-star'), null);
  assert.ok(b.root.querySelector('.fa-star'));
  a.dispose();
  b.dispose();
});

picker('tableInput: the Files action adds the resolved frame to the workspace and picks it as a user edit',
  async () => {
    let resolve;
    await withGlobals({grok_UI_PickTableFromFiles: () => new Promise((r) => resolve = r)},
      async () => {
        const input = tableInput('Table');
        let systemic = null;
        input.effect(() => {
          input.value;
          systemic = Input.isSystemWrite;
        });
        fire(input.root.querySelector('.u2-input-options [aria-label="Add file from Files"]'), 'click');
        resolve(new DataFrame([{name: 'a', type: 'int'}], [], 'plates'));
        await flush();
        assert.deepEqual(opened(), ['plates'],
          'the frame joins the workspace (divergence #16 — u2, unlike Dart, adds it)');
        assert.deepEqual(options(input), ['plates']);
        assert.equal(input.value, 'plates');
        assert.equal(systemic, false, 'picked the way a user picks one');
        input.dispose();
      });
  });

picker('tableInput: a cancelled dialog leaves the input alone', async () => {
  await withGlobals({grok_UI_PickTableFromQuery: async () => null}, async () => {
    grok.openTable('demog');
    const input = tableInput('Table');
    input.value = 'demog';
    fire(input.root.querySelector('.u2-input-options [aria-label="Query database"]'), 'click');
    await flush();
    assert.deepEqual(opened(), ['demog'], 'a null resolution adds nothing');
    assert.equal(input.value, 'demog', 'and takes nothing away');
    input.dispose();
  });
});

picker('tableInput: a pick landing after dispose is dropped', async () => {
  let resolve;
  await withGlobals({grok_UI_PickTableFromFiles: () => new Promise((r) => resolve = r)},
    async () => {
      const input = tableInput('Table');
      fire(input.root.querySelector('.u2-input-options [aria-label="Add file from Files"]'), 'click');
      input.dispose();
      resolve(new DataFrame([], [], 'late'));
      await flush();
      assert.deepEqual(opened(), [], 'a disposed input adds nothing');
    });
});

picker('tableInput: detection happens at construction — the setter renders what it was given', () => {
  const input = withGlobals(
    {grok_UI_PickTableFromFiles: async () => null, grok_UI_PickTableFromQuery: async () => null},
    () => tableInput('Table'));
  input.actions = input.actions;
  assert.deepEqual(input.actions.map((a) => a.id), ['open-file', 'add-from-files', 'query-db'],
    'the doors leaving does not thin what was already given');
  assert.equal(railIcons(input).length, 3);

  const late = tableInput('Late');
  assert.deepEqual(late.actions.map((a) => a.id), ['open-file'], 'while a new input detects afresh');
  input.dispose();
  late.dispose();
});

picker('tableInput: a rejecting picker dialog is reported at the boundary, not unhandled', async () => {
  const errors = [];
  grok.shell.error = (m) => errors.push(m);
  await withGlobals({grok_UI_PickTableFromFiles: async () => {
    throw new Error('boom');
  }}, async () => {
    const input = tableInput('Table');
    fire(input.root.querySelector('.u2-input-options [aria-label="Add file from Files"]'), 'click');
    await flush();
    assert.deepEqual(errors, ['Error: boom'], 'the rejection surfaces as a shell error');
    assert.equal(input.peek(), null, 'a failed pick writes nothing');
    assert.equal(input.scope.isDisposed, false, 'the input survives');
    input.dispose();
  });
});

picker('TableInput.actions: ten re-sets leave the input scope flat', () => {
  const input = tableInput('Table');
  const disposers = input.scope._disposers?.size ?? 0;
  for (let i = 0; i < 10; i++)
    input.actions = [{id: 'mine', icon: 'star', tooltip: 'Mine', run: () => {}}];
  assert.equal(input.scope._disposers?.size ?? 0, disposers,
    'no per-render cleanup piles up on the input');
  assert.equal(railIcons(input).length, 1);
  input.dispose();
});

// --- what the platform lists (src/dg/shell/platform.ts) over the hostile doubles ---

/** A campaign table read back from Chem: `~Core.canonicalSmiles` first, `~Scaffold.canonicalSmiles` added by the
 * frame double, a hidden numeric cache and a column with no name before `MW`. */
const CAMPAIGN = () => new DataFrame([
  {name: '~Core.canonicalSmiles', type: 'string', semType: 'Molecule'}, {name: 'ID', type: 'string'},
  {name: 'Core', type: 'string', semType: 'Molecule'}, {name: 'Scaffold', type: 'string', semType: 'Molecule'},
  {name: '~MW.cache', type: 'double'}, {name: '', type: 'int'}, {name: 'MW', type: 'double'},
], [], 'campaign');
const LISTED = ['ID', 'Core', 'Scaffold', 'MW'];
const isMolecule = (c) => c.semType === 'Molecule';
const keys = (input) => input.getItems({all: true}).list.map((i) => i.key);

picker('columnInput: no variant lists a hidden ~ column or a nameless one; the filter narrows what is listed', () => {
  for (const variant of [{}, {rich: true}, {grid: true}]) {
    const table = CAMPAIGN();
    const all = columnInput('Column', table, variant);
    const molecules = columnInput('Molecules', table, {...variant, filter: isMolecule});
    assert.deepEqual(keys(all), LISTED, JSON.stringify(variant));
    assert.deepEqual(keys(molecules), ['Core', 'Scaffold'], JSON.stringify(variant));
    all.dispose();
    molecules.dispose();
  }
  const plain = columnInput('Column', CAMPAIGN());
  assert.deepEqual(options(plain), LISTED, 'the select itself');
  plain.dispose();
});

picker('columns, columns map and aggregations list no hidden column, and default to none', () => {
  const table = CAMPAIGN();
  const many = columnsInput('Columns', table, {value: LISTED});
  assert.deepEqual(keys(many), LISTED);
  assert.equal(many.root.querySelector('.u2-columns-summary').textContent, '(4) All',
    'every listed column picked is all of them');

  const map = columnsMapInput('Map', ['structure'], table);
  assert.deepEqual(keys(map.input('structure')), LISTED);

  const aggregations = aggregatedColumnsInput('Aggregations', table);
  fire(aggregations.root.querySelector('[aria-label="Add aggregation"]'), 'click');
  assert.deepEqual(aggregations.getState().value, [{column: 'MW', aggregation: 'avg'}],
    'a new row starts on the first listed numerical column');
  assert.deepEqual([...aggregations.root.querySelector('select').querySelectorAll('option')].map((o) => o.value), LISTED,
    'the row\'s column choice');

  const combo = new ColumnInput({label: 'Column', table});
  fire(combo.root.querySelector('.u2-columns'), 'keydown', {key: 'ArrowDown'});
  assert.equal(combo.getState().value, 'ID', 'an arrow from nothing steps onto the first listed column');
  for (const input of [many, map, aggregations, combo])
    input.dispose();
});

picker('columnInput: takes value, nullable and required, as every input does — each variant', () => {
  for (const variant of [{}, {rich: true}, {grid: true}]) {
    const input = columnInput('Molecules', CAMPAIGN(),
      {...variant, filter: isMolecule, value: 'Core', nullable: false, required: true, tooltipText: 'Structures'});
    const state = input.getState();
    assert.deepEqual([state.value, state.required, input.nullable, input.tooltipText],
      ['Core', true, false, 'Structures'], JSON.stringify(variant));
    input.dispose();
  }
  const strict = columnInput('Column', CAMPAIGN(), {nullable: false});
  assert.deepEqual([...strict.root.querySelectorAll('option')].map((o) => o.value), LISTED, 'not nullable: no empty choice');
  strict.dispose();
});

picker('tableInput: table is the picked frame itself — of two open under one name, the one given', () => {
  const first = grok.openTable('demog');
  const cars = grok.openTable('cars');
  const second = grok.openTable('copy');
  second.name = 'demog';
  const input = tableInput('Table');
  const many = tablesInput('Tables');
  assert.deepEqual(options(input), ['demog', 'cars'], 'a name two tables share is offered once');
  assert.deepEqual([...many.root.querySelectorAll('input[type=checkbox]')].map((b) => b.value), ['demog', 'cars']);
  const seen = [];
  input.effect(() => seen.push(input.table));
  assert.equal(input.table, null);

  input.table = second;
  assert.equal(input.getState().value, 'demog');
  assert.ok(input.table === second, 'the frame given, not the first of its name');
  assert.ok(seen.at(-1) === second, 'tracked');

  input.table = cars;
  assert.ok(input.table === cars);
  const select = input.root.querySelector('select');
  select.value = 'demog';
  fire(select, 'change');
  assert.ok(input.table === first, 'a pick by name is the first table of that name');

  input.table = second;
  grok.shell.dart.tables = grok.shell.dart.tables.filter((t) => t !== second);
  grok.shell.dart.tableRemoved.fire();
  assert.ok(input.table === first, 'once the frame given is closed, the name answers');
  input.dispose();
  many.dispose();
});

picker('tableInput: a table it does not offer — nameless, or not open — throws in strict mode; in production it ' +
  'is reported, and clears the pick', () => {
  const demog = grok.openTable('demog');
  const input = tableInput('Table');
  input.table = demog;
  assert.throws(() => input.table = new DataFrame([], [], 'loose'), {name: 'U2UsageError', code: 'table-not-offered'});
  const mode = Usage.mode;
  Usage.mode = 'report';
  try {
    const reports = () => Usage.reports().filter((r) => r.code === 'table-not-offered').length;
    const before = reports();
    input.table = grok.shell.tables[0];
    assert.equal(grok.shell.tables[0].name, undefined, 'the stray table renamed to null');
    assert.deepEqual([input.table, input.getState().value, reports()], [null, null, before + 1]);
    input.table = demog;
    input.table = new DataFrame([], [], 'loose');
    assert.deepEqual([input.table, input.getState().value, reports()], [null, null, before + 2]);
    input.dispose();
  } finally {
    Usage.mode = mode;
  }
});
