/* `domains.grid` (2-4b) over the platform doubles: the source's frame is the grid's, the js-api
   editor is attached with the source's edit state — decorated first — re-attached over the next
   frame and released when the source lets go, and the `u2-domain-grid` tag. The real grid runs
   only in the platform — U2Demo `U2: domain session` covers that. `DG` comes from tests/dg-stub.mjs. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, signal, Registry, SpecContext, renderSpec, registerAll} from '@datagrok-libraries/u2';
import {DomainSource} from '../src/core/domain-source.js';
import {DataFrame, Property, Stream, WidgetDescriptor, platform} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainGrid} = await import('../src/dg/grid.js');
const {EditorEditState} = await import('../src/dg/editor-state.js');
const {registerDomainComponents} = await import('../src/dg/registrations.js');
const DG = await import('datagrok-api/dg');

/** The look the grid host reads to fit its header, columns and height. */
const LOOK = [new Property('colHeaderHeight', 'int', {defaultValue: 22}),
  new Property('colHeaderFont', 'string', {defaultValue: 'bold 13px Roboto'}), new Property('rowHeight', 'int', {defaultValue: 28})];

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', [new Property('allowEdit', 'bool', {defaultValue: false}),
      ...LOOK])];
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      WidgetDescriptor.registry = [];
      DG.DomainObjectHandler.decorated.length = 0;
      DG.DomainObjectHandler.names = {};
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const env = {designTime: false, subBinds: {}, resolve: () => null};

/** The columns the frame double carries — `project_id` is the ref column the platform's ref-cell
 * renderer keys on (its semType is the target's row type). */
const COLUMNS = [{name: 'id', type: 'string'}, {name: 'created_on', type: 'string'},
  {name: 'title', type: 'string'}, {name: 'priority', type: 'string'}, {name: 'due', type: 'datetime'},
  {name: 'project_id', type: 'string', semType: 'grit.project'}];

/** The js-api editor's surface the grid host and the edit state touch, over a frame double. */
function fakeEditor(rows, columns = COLUMNS) {
  const df = new DataFrame(columns, rows);
  return {
    df,
    editor: {
      table: 'grit.issue', isDirty: false, changeCount: 0, isSaving: false, detached: 0,
      onChanged: new Stream(), onDirtyChanged: new Stream(), onSavingChanged: new Stream(), onSaved: new Stream(),
      onRefused: new Stream(),
      get dataFrame() { return df; },
      stateOf: () => '', errorsOf: () => ({}), isChanged: () => false, errorOf: () => null,
      setValue() {}, addRow: () => -1, markDeleted() {}, unmarkDeleted() {}, discard() {},
      save: async () => true,
      detach() { this.detached++; },
    },
  };
}

/** The memory backend's tables, their frames replaced by a frame double with a fake js-api editor
 * attached — what the dg backend hands a source. */
function editorBackend(memory, columns = COLUMNS) {
  const editors = [];
  const cell = (row, c) => {
    const value = c.from === undefined ? row[c.name] : c.from(row);
    return c.type === 'datetime' && value != null ? new Date(value) : value ?? null;
  };
  return {editors, table: async (address) => {
    const t = await memory.table(address);
    return {
      ...t, access: () => t.access(), count: (f, s) => t.count(f, s), query: (s) => t.query(s),
      transaction: (ops) => t.transaction(ops),
      frame: async (spec) => {
        const rows = await t.query(spec);
        const {df, editor} = fakeEditor(rows.map((r) =>
          Object.fromEntries(columns.map((c) => [c.name, cell(r, c)]))), columns);
        editors.push(editor);
        const edit = new EditorEditState(editor);
        return {df, edit, append: async () => 0, dispose: () => edit.dispose()};
      },
    };
  }};
}

scoped('over a memory source: the frame is the grid\'s; the memory writer is not a platform editor', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  assert.equal(grid.root.dataset.u2, 'domain-grid');
  assert.equal(grid.root.classList.contains('u2-domain-grid'), true);
  assert.equal(grid.grid.type, 'Grid');
  assert.equal(grid.grid.root.parentElement, grid.root);
  assert.equal(grid.grid.dataFrame, src.df.value, 'the source\'s frame reaches the viewer');
  assert.equal(grid.grid.editor, null);
  assert.equal(DG.DomainObjectHandler.decorated.length, 0, 'nothing to decorate for');
  grid.dispose();
  src.dispose();
});

scoped('the column header fits two lines of caption, unless the caller sets its own height', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  const own = domains.grid(src, {look: {colHeaderHeight: 60}});
  await flush();
  assert.equal(grid.grid.props.colHeaderHeight, 44, 'a unit on the second line is not cut');
  assert.equal(own.grid.props.colHeaderHeight, 60);
  own.grid.onColumnResized.fire();
  assert.equal(own.grid.props.colHeaderHeight, 60, 'a resize does not refit a height the caller set');
  grid.dispose();
  own.dispose();
  src.dispose();
});

scoped('the header grows to the longest caption at its column\'s width, between two and four lines', () => {
  const measure = (text) => text.length * 7;
  assert.equal(DomainGrid.captionLines('Title', 92, measure), 1);
  assert.equal(DomainGrid.captionLines('Dosed Concentration (uM)', 92, measure), 3);
  assert.equal(DomainGrid.captionLines('Reported CLint (uL/min/10^6 cells)', 92, measure), 4);
  assert.equal(DomainGrid.captionLines('Reported CLint (uL/min/10^6 cells)', 300, measure), 1);
  assert.equal(DomainGrid.captionLines('Supercalifragilistic', 92, measure), 1, 'a long word takes its own line');
  assert.equal(DomainGrid.headerHeight(0), 44, 'never under two lines');
  assert.equal(DomainGrid.headerHeight(2), 44);
  assert.equal(DomainGrid.headerHeight(3), 61, 'a third line of the 13px header font');
  assert.equal(DomainGrid.headerHeight(9), 78, 'never over four');
  assert.equal(DomainGrid.headerHeight(3, 'bold 11px Roboto'), 55, 'the grid\'s own header font');
});

scoped('a cell keeps to one line, cut with an ellipsis', async () => {
  emsBackends.domain = editorBackend(backend());
  const src = new DomainSource({table: 'grit.issue', pageSize: 10}, env);
  const grid = domains.grid(src);
  src.start();
  await flush();
  const columns = grid.grid.columns.toList();
  assert.equal(columns.length, COLUMNS.length);
  for (const column of columns)
    assert.equal(column.contentCellStyle.textWrap, 'ellipsis', column.name);
  src.dispose();
  grid.dispose();
});

scoped('the js-api editor is attached with the source\'s edit state, decorated first, re-attached per frame, released with it', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  const src = new DomainSource({table: 'grit.issue', pageSize: 10}, env);
  const grid = domains.grid(src);
  assert.equal(grid.grid.editor, null, 'nothing before the first frame');
  src.start();
  await flush();
  assert.equal(src.state.value, 'ready');
  assert.equal(grid.grid.editor, inner.editors[0], 'the source\'s own editor');
  assert.equal(grid.grid.dataFrame, inner.editors[0].dataFrame, 'over the editor\'s frame');
  assert.equal(DG.DomainObjectHandler.decorated.length, 1);
  const [d] = DG.DomainObjectHandler.decorated;
  assert.deepEqual([d.grid, d.table, d.dataFrame], [grid.grid, 'grit.issue', inner.editors[0].dataFrame]);
  await src.refresh();
  await flush();
  assert.equal(inner.editors.length, 2);
  assert.equal(inner.editors[0].detached, 1, 'the old writer left with its frame');
  assert.equal(grid.grid.editor, inner.editors[1], 're-attached over the next frame');
  assert.equal(grid.grid.dart.attached.length, 2);
  src.dispose();
  assert.equal(grid.grid.editor, null, 'released when the source lets the frame go');
  assert.equal(inner.editors[1].detached, 1);
  grid.dispose();
});

scoped('spec: u2-domain-grid over a bound source', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-grid');
  assert.equal(meta.category, 'Display');
  assert.deepEqual(meta.props.map((p) => p.name), ['source']);
  assert.match(meta.usage, /Save and Discard are the session's/);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-grid', name: 'grid',
    bind: {source: '$.issues'}}}, new SpecContext({data: {issues: signal(src)}}), reg);
  await flush();
  const grid = instance.node('grid');
  assert.equal(grid instanceof DomainGrid, true);
  assert.equal(grid.source, src);
  assert.equal(grid.grid.dataFrame, src.df.value);
  instance.dispose();
  src.dispose();
});

scoped('on attach: the schema\'s captions as friendly names, the system columns hidden', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  const columns = grid.grid.columns;
  assert.deepEqual([columns.byName('id').visible, columns.byName('created_on').visible], [false, false],
    'the system columns are not the grid\'s to show');
  assert.equal(columns.byName('title').visible, true);
  assert.equal(columns.byName('priority').caption, 'Priority', 'the schema\'s caption');
  assert.equal(columns.byName('title').caption, 'title', 'a column the schema does not caption stays raw');
  assert.equal(columns.byName('priority').name, 'priority', 'the grid column keeps the column\'s name');
  assert.equal(src.df.value.columns.byName('due').getTag('format'), 'MMM d, yyyy',
    'a datetime column whose values all fall on midnight is a date');
  grid.dispose();
  src.dispose();
});

scoped('a column is as wide as its data and its caption in two lines, over a floor by type and under a cap', () => {
  const measure = (text) => text.length * 7;
  assert.equal(DomainGrid.captionWidth('Text', measure), 28);
  assert.equal(DomainGrid.captionWidth('Dosed Concentration (uM)', measure), 126, '"Dosed" over "Concentration (uM)"');
  const text = {name: 'text', type: 'string'};
  assert.equal(DomainGrid.columnWidth(text, false, 20, 36), 96, 'an empty text column still reads');
  assert.equal(DomainGrid.columnWidth(text, false, 210, 36), 210, 'the data\'s own width');
  assert.equal(DomainGrid.columnWidth(text, false, 900, 36), 320, 'one long text does not take the pane');
  assert.equal(DomainGrid.columnWidth({name: 'n', type: 'double'}, false, 40, 141), 141, 'the caption\'s width');
  assert.equal(DomainGrid.columnWidth({name: 'b', type: 'bool'}, false, 30, 20), 48, 'a flag the least');
  assert.equal(DomainGrid.columnWidth({name: 'p', type: 'string', semType: 'grit.project'}, false, 280, 50), 160,
    'a reference draws a name, not the id the data width measures');
  assert.equal(DomainGrid.columnWidth(text, true, 20, 36), 160, 'the name column the most');
});

scoped('the shown columns are fitted once per frame; the height follows the rows', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  const columns = grid.grid.columns;
  assert.equal(columns.byName('title').width, 160, 'the name column');
  assert.equal(columns.byName('priority').width, 96);
  assert.equal(columns.byName('id').width, undefined, 'a hidden column is left alone');
  columns.byName('priority').width = 300;
  grid.grid.onDataFrameChanged.fire();
  assert.equal(columns.byName('priority').width, 300, 'a width the user set stays over the same frame');
  const rows = src.df.value.rowCount;
  const height = () => grid.root.style.getPropertyValue('--u2-domain-grid-height');
  assert.equal(height(), `${44 + rows * 28 + 12}px`, 'the header, the rows and a scrollbar');
  src.df.value.dart.rows.push({});
  src.df.value.onRowsAdded.fire();
  assert.equal(height(), `${44 + (rows + 1) * 28 + 12}px`);
  grid.dispose();
  src.dispose();
});

scoped('hiddenColumns: what the caller hides goes with the system columns; a value with a time keeps it', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src, {hiddenColumns: ['priority']});
  await flush();
  assert.equal(grid.grid.columns.byName('priority').visible, false, 'a child grid hides its FK this way');
  assert.equal(grid.grid.columns.byName('title').visible, true);
  const due = src.df.value.columns.byName('due');
  assert.equal(DomainGrid.isDateOnly(due), true);
  src.df.value.dart.rows[1].due = new Date('2026-01-02T09:30:00Z');
  assert.equal(DomainGrid.isDateOnly(due), false, 'one value with a time makes the column a datetime');
  grid.dispose();
  src.dispose();
});

scoped('ref cells draw the target row\'s caption; an id the platform cannot name keeps the id', async () => {
  emsBackends.domain = editorBackend(backend());
  DG.DomainObjectHandler.names['grit.project|p1'] = 'Grit';
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  const refs = grid.grid.columns.byName('project_id');
  assert.deepEqual([refs.cellText(0), refs.cellText(1)], ['Grit', 'Grit'], 'the target row\'s name, not its id');
  assert.equal(refs.cellText(2), 'p2', 'an id with no name falls back to the id');
  assert.equal(grid.grid.columns.byName('title').cellText(0), 'Aspirin', 'a plain column draws its value');
  grid.dispose();
  src.dispose();
});

scoped('the frame is decorated once the GRID holds it, and again over every frame it repoints to', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  DG.DomainObjectHandler.names['grit.project|p1'] = 'Grit';
  const src = new DomainSource({table: 'grit.issue', pageSize: 10}, env);
  const grid = domains.grid(src);
  src.start();
  await flush();
  const held = () => DG.DomainObjectHandler.decorated.every((d) => d.held);
  assert.equal(DG.DomainObjectHandler.decorated.length, 1);
  assert.equal(held(), true, 'never over a frame the grid does not hold');
  await src.refresh();
  await flush();
  assert.equal(DG.DomainObjectHandler.decorated.length, 2, 'the next frame is decorated too');
  assert.equal(DG.DomainObjectHandler.decorated[1].dataFrame, inner.editors[1].dataFrame);
  assert.equal(held(), true);
  assert.equal(grid.grid.columns.byName('project_id').cellText(0), 'Grit', 'the new frame\'s ref cells too');
  grid.dispose();
  src.dispose();
});

scoped('a join table\'s two ref columns both draw captions', async () => {
  emsBackends.domain = editorBackend(backend(), [{name: 'id', type: 'string'},
    {name: 'project_id', type: 'string', semType: 'grit.project'},
    {name: 'issue_id', type: 'string', semType: 'grit.issue', from: (row) => row.id}]);
  Object.assign(DG.DomainObjectHandler.names,
    {'grit.project|p1': 'Grit', 'grit.issue|i1': 'Aspirin', 'grit.issue|i2': 'Ibuprofen'});
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  const columns = grid.grid.columns;
  assert.deepEqual([columns.byName('project_id').cellText(0), columns.byName('issue_id').cellText(0)],
    ['Grit', 'Aspirin']);
  assert.deepEqual([columns.byName('project_id').cellText(1), columns.byName('issue_id').cellText(1)],
    ['Grit', 'Ibuprofen']);
  grid.dispose();
  src.dispose();
});

scoped('a refusal the grid raised over a read-only cell reaches the status bar, and no second balloon', async () => {
  const inner = editorBackend(backend());
  emsBackends.domain = inner;
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const grid = domains.grid(src);
  await flush();
  const [editor] = inner.editors;
  assert.equal(editor.onRefused.count, 1, 'the grid host follows the editor it attached');
  editor.onRefused.fire({row: 0, column: 'title', message: 'Label is read-only'});
  await flush();
  assert.equal(src.summary.value, 'Cannot save: Label is read-only');
  assert.equal(document.body.querySelectorAll('.u2-notify-error').length, 0,
    'the grid balloons it already — the source does not say it twice');
  editor.onChanged.fire(editor);
  await flush();
  assert.equal(src.error.value, undefined, 'the next edit takes the refusal back');
  await src.refresh();
  await flush();
  assert.equal(editor.onRefused.count, 0, 'the writer that left the frame is not followed');
  assert.equal(inner.editors[1].onRefused.count, 1, 're-subscribed over the next frame');
  grid.dispose();
  assert.equal(inner.editors[1].onRefused.count, 0);
  src.dispose();
});
