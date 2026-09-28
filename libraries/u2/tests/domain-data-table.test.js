/* `domains.dataTable` (WO 3-10) over the memory backend: the columns come from the table's
   schema the way the grid's decoration derives them (the name column first, the system and `~`
   service ones out, a hidden column absent), the cells carry the source writer's verdicts by ROW
   KEY, and an edit repaints the window that is already rendered without swapping the items. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {backends} from '../src/sources/backends.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {timestamp} from '../src/core/elements.js';
import {Rows} from '../src/sources/rows-like.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {DomainDataTable} = await import('../src/dg/domain/grid.js');
const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const headers = (dt) => dt.root.querySelectorAll('.u2-data-table-head').map((el) => el.textContent);
const row = (dt, i) => dt.root.querySelector(`.u2-data-table-row[data-index="${i}"]`);
const cell = (dt, i, column) => row(dt, i).children[dt.columns().indexOf(column)];

async function dataTable(options = {}) {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const dt = domains.dataTable(src, options);
  document.body.append(dt.root);
  await flush();
  dt.table.value.root.clientHeight = 400;
  await flush();
  return {table, src, dt};
}

scoped('the columns are the schema\'s: the name column first, system and service columns out', async () => {
  const {src, dt} = await dataTable();
  assert.equal(dt.root.dataset.u2, 'domain-data-table');
  assert.equal(dt.columns()[0], 'title', 'the name column leads');
  for (const system of ['id', 'version', 'created_on', 'updated_on', 'author_id'])
    assert.equal(dt.columns().includes(system), false, `${system} is not the table's to show`);
  assert.equal(dt.columns().some((c) => c.startsWith('~')), false, 'no service column');
  assert.deepEqual(dt.columns(), ['title', 'project_id', 'number', 'description', 'done', 'reporter',
    'tags', 'due', 'weight', 'priority']);
  assert.equal(headers(dt)[0], 'title');
  assert.equal(headers(dt).includes('Priority'), true, 'the schema\'s caption is the header');
  assert.equal(cell(dt, 0, 'title').textContent, 'Aspirin');
  assert.equal(cell(dt, 0, 'weight').textContent, '1.5');
  dt.dispose();
  src.dispose();
});

scoped('columns and hiddenColumns narrow the selection', async () => {
  const {src, dt} = await dataTable({columns: ['title', 'id', 'priority'], hiddenColumns: ['priority']});
  assert.deepEqual(dt.columns(), ['title'], 'a system column and a hidden one are dropped from the request too');
  assert.deepEqual(headers(dt), ['title']);
  dt.dispose();
  src.dispose();

  const second = await dataTable({hiddenColumns: ['project_id']});
  assert.equal(second.dt.columns().includes('project_id'), false);
  second.dt.dispose();
  second.src.dispose();
});

scoped('an edit through the source paints its cell, and a refusal paints and explains it', async () => {
  const {src, dt} = await dataTable();
  const changed = () => dt.root.querySelectorAll('.u2-cell-changed').map((c) => c.textContent);
  assert.deepEqual(changed(), []);

  src.rows.items.peek()[1].priority = 'high';
  await flush();
  assert.deepEqual(changed(), ['high'], 'the pending cell, keyed by the row it belongs to');
  assert.equal(cell(dt, 1, 'priority').className.includes('u2-cell-error'), false);

  src.rows.items.peek()[1].title = '';
  await flush();
  const invalid = cell(dt, 1, 'title');
  assert.equal(invalid.className.includes('u2-cell-error'), true);
  assert.equal(invalid.title, 'Value can\'t be empty', 'the refusal is the cell\'s tooltip');
  assert.equal(cell(dt, 0, 'title').className.includes('u2-cell-error'), false, 'and only that cell');

  src.discard();
  await flush();
  assert.deepEqual(changed(), [], 'a discard takes the paint back');
  dt.dispose();
  src.dispose();
});

scoped('a cell the caller may not write reads as read-only', async () => {
  backends.domain = backend({access: {can: {view: true, edit: true},
    fields: {id: 'readonly', title: 'editable', priority: 'readonly'}}});
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const dt = domains.dataTable(src, {columns: ['title', 'priority']});
  document.body.append(dt.root);
  await flush();
  dt.table.value.root.clientHeight = 400;
  await flush();
  assert.equal(cell(dt, 0, 'priority').className.includes('u2-cell-readonly'), true);
  assert.equal(cell(dt, 0, 'title').className.includes('u2-cell-readonly'), false);
  dt.dispose();
  src.dispose();
});

scoped('the selection and the source\'s current row are one thing', async () => {
  const {src, dt} = await dataTable();
  dt.table.value.selectedIndex.value = 2;
  await flush();
  assert.equal(src.currentRow.value.title, 'Naproxen');
  src.currentRow.value = src.rows.byKey('i1');
  await flush();
  assert.equal(dt.table.value.selectedIndex.value, 0);
  dt.dispose();
  src.dispose();
});

scoped('a reference cell shows the row it points at, not its uuid', async () => {
  // the caption the query projected rides with the row: the first paint is the name, with no
  // lookup at all — the target's table is held open here, as a round trip holds it, and the cell
  // is right anyway
  const inner = backend();
  backends.domain = {saveAll: (edits) => inner.saveAll(edits),
    table: (address) => address === 'grit.project' ? new Promise(() => {}) : inner.table(address)};
  const early = (await domains.table('grit.issue')).source({pageSize: 10});
  const first = domains.dataTable(early, {columns: ['project_id']});
  document.body.append(first.root);
  await flush();
  first.table.value.root.clientHeight = 400;
  // a scroller with no height defers its first window to the next frame (list.ts `_later`), and
  // the shim only settles those on a flush
  await flush();
  assert.equal(cell(first, 0, 'project_id').textContent, 'Grit');
  first.dispose();

  // with the captions off the one lookup stands in for them, and until it answers a placeholder,
  // never the id — a uuid flashing in the cell reads as the value of the column
  const off = (await domains.table('grit.issue')).source({pageSize: 10, captions: []});
  const none = domains.dataTable(off, {columns: ['project_id']});
  document.body.append(none.root);
  await flush();
  none.table.value.root.clientHeight = 400;
  await flush();
  assert.equal(cell(none, 0, 'project_id').textContent, '…');
  none.dispose();
  off.dispose();
  early.dispose();

  const {src, dt} = await dataTable({columns: ['title', 'project_id']});
  await flush();
  assert.equal(cell(dt, 0, 'project_id').textContent, 'Grit', 'the name column of the row it points at');
  assert.equal(cell(dt, 2, 'project_id').textContent, 'Datagrok');
  assert.equal(cell(dt, 0, 'project_id').textContent.includes('p1'), false, 'never the raw id');

  // an id the target does not hold answers as itself, as the picker's resolve does — the lookup
  // path is the only one that ever sees a raw id
  const loose = (await domains.table('grit.issue')).source({pageSize: 10, captions: []});
  const lt = domains.dataTable(loose, {columns: ['project_id']});
  document.body.append(lt.root);
  await flush();
  lt.table.value.root.clientHeight = 400;
  await flush();
  await flush();
  assert.equal(cell(lt, 0, 'project_id').textContent, 'Grit', 'the lookup answers the target\'s name');
  loose.rows.items.peek()[1].project_id = 'nope';
  await flush();
  await flush();
  assert.equal(cell(lt, 1, 'project_id').textContent, 'nope');
  lt.dispose();
  loose.dispose();
  dt.dispose();
  src.dispose();
});

scoped('a datetime cell is formatted, and a float loses its packing noise', async () => {
  const {src, dt} = await dataTable({columns: ['due', 'weight']});
  const due = cell(dt, 0, 'due');
  assert.notEqual(due.querySelector('.u2-timestamp'), null, 'the shared timestamp, as every other u2 surface');
  assert.equal(due.textContent.includes('2026-01-01T00:00:00Z'), false, 'never the raw ISO string');
  // a value stamped at midnight UTC is a DATE: read in the local zone it would be the day before
  // everywhere west of Greenwich (the runner sits at UTC-5), so it is printed from the UTC parts
  assert.equal(due.textContent, 'Jan 1, 2026');
  assert.equal(due.textContent, timestamp('2026-01-01T00:00:00Z', undefined, {utcDates: true}).textContent);
  assert.notEqual(due.querySelector('.u2-timestamp').title, '', 'the day is the tooltip too');
  assert.equal(cell(dt, 1, 'due').textContent, '', 'an empty cell stays empty');

  src.rows.items.peek()[0].weight = 1.7999999999999998;
  await flush();
  assert.equal(cell(dt, 0, 'weight').textContent, '1.8');
  assert.equal(cell(dt, 1, 'weight').textContent, '0.5');
  dt.dispose();
  src.dispose();
});

scoped('the column tracks follow the column types, so a wide table is not a wall of truncation', async () => {
  const {src, dt} = await dataTable();
  const tracks = dt.root.querySelector('.u2-data-table-header').style.gridTemplateColumns.split(') ')
    .map((t) => t.endsWith(')') ? t : `${t})`);
  const of = (column) => tracks[dt.columns().indexOf(column)];
  assert.equal(of('title'), 'minmax(120px, 1.5fr)', 'the name column takes the most room');
  assert.equal(of('done'), 'minmax(48px, 0.5fr)', 'a flag the least');
  assert.equal(of('number'), 'minmax(48px, 0.5fr)');
  assert.equal(of('weight'), 'minmax(48px, 0.5fr)');
  assert.equal(of('due'), 'minmax(96px, 0.8fr)');
  assert.equal(of('description'), 'minmax(72px, 1fr)');
  assert.equal(of('project_id'), 'minmax(96px, 1.5fr)', 'a reference reads as a name, not as an id');
  // ten columns of minimums must still fit a pane, or the last one is off-screen behind a
  // scrollbar the rows (absolutely positioned) never paint
  const floor = tracks.reduce((sum, t) => sum + Number(/minmax\((\d+)px/.exec(t)[1]), 0);
  assert.equal(floor <= 800, true, `the minimums sum to ${floor}px`);
  const numeric = dt.root.querySelector('.u2-data-table-row').children[dt.columns().indexOf('weight')];
  assert.equal(numeric.className.includes('u2-data-table-align-right'), true, 'numbers read right-aligned');
  assert.equal(dt.root.querySelectorAll('.u2-data-table-head')[0].title, 'title',
    'a truncated header says what it is');
  dt.dispose();
  src.dispose();
});

scoped('spec: u2-domain-data-table over a bound source', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-data-table');
  assert.equal(meta.category, 'Display');
  assert.deepEqual(meta.props.map((p) => p.name).slice(0, 4),
    ['source', 'columns', 'hiddenColumns', 'rowHeight']);
  assert.match(meta.usage, /u2-domain-grid/);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-data-table', name: 'rows',
    bind: {source: '$.issues'}, props: {columns: ['title', 'priority']}}},
  new SpecContext({data: {issues: signal(src)}}), reg);
  await flush();
  const dt = instance.node('rows');
  assert.equal(dt instanceof DomainDataTable, true);
  assert.equal(dt.source, src);
  assert.deepEqual(dt.columns(), ['title', 'priority']);
  instance.dispose();
  src.dispose();
});

scoped('a ref cell takes the caption the query projected; only a frame without one costs a lookup', async () => {
  const {table, src, dt} = await dataTable({columns: ['title', 'project_id']});
  const prop = table.properties.find((p) => p.name === 'project_id');

  const carried = DomainDataTable.caption(prop, {project_id: 'p1', [Rows.caption('project_id')]: 'Grit'},
    'project_id', 'p1');
  assert.equal(carried, 'Grit', 'the wire\'s answer on the first paint — a string, so nothing to repaint');

  assert.equal(DomainDataTable.caption(prop, {project_id: 'p1', [Rows.caption('project_id')]: null},
    'project_id', 'p1'), '', 'a null caption is an empty cell, not a miss — the caller may not see the target');

  // a draft, a `User` and a `Group` are never projected: the one lookup stays, and never the uuid
  const cache = new Map();
  const looked = DomainDataTable.caption(prop, {project_id: 'p1'}, 'project_id', 'p1', cache);
  assert.equal(looked.textContent, '…');
  await flush();
  await flush();
  assert.equal(looked.textContent, 'Grit');
  assert.equal(DomainDataTable.caption(prop, {project_id: 'p1'}, 'project_id', 'p1', cache), 'Grit',
    'a second render reads what the lookup answered — one lookup per id, not per paint');
  dt.dispose();
  src.dispose();
});
