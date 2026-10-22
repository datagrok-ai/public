/* The final UX round over the memory backend: a record page's pending delete, the list's order, the
   captions a frame carries in the target's format, the grid's number formats, the folded filters
   pane, the platform's current object, a child tab's New, the poll's count and the host answers a
   load reuses. */

import {emsBackends} from '../src/core/domain-backend.js';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {DataFrame, Property, Stream, WidgetDescriptor, platform} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';
import {backends, Scope, Rows, notify, Filters, allowedActions} from '@datagrok-libraries/u2';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains, DomainTable} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');
const {EditorEditState} = await import('../src/dg/editor-state.js');
const grok = await import('datagrok-api/grok');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    const dg = globalThis.DG;
    try {
      globalThis.DG = {...dg, format: (x, format) => `${format.replace(/0+$/, '')}${String(x).padStart(6, '0')}`};
      await body();
    } finally {
      globalThis.DG = dg;
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
    assert.equal(DomainApp.live.size, 0, 'no app left behind');
  });
}

const BASE = '/apps/T/Issues';
const titles = (a) => a.listSource.rows.items.value.map((r) => r.title);

async function app(options = {}) {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const a = domains.app({table, base: BASE, pageSize: 10, ...options});
  document.body.append(a.root);
  await flush();
  return {table, app: a};
}

/** Registration-shaped tables: a key autonumbered in a display format, a reference to one, a row keyed by one. */
const REG = {
  name: 'reg',
  tables: {
    compound: {businessKey: ['regno'], columns: {regno: {type: 'int', format: 'DG-000000'}, label: {type: 'string'}}},
    batch: {businessKey: ['batch_regno'],
      columns: {compound_id: {type: 'ref', ref: 'compound'}, batch_regno: {type: 'int', format: 'DGB-000000'},
        made: {type: 'datetime', format: 'yyyy-MM-dd'}}},
    note: {businessKey: ['batch_id'], columns: {batch_id: {type: 'ref', ref: 'batch'}, text: {type: 'string'}}},
  },
};
const REG_ROWS = {
  compound: [{id: 'c1', regno: 7, label: 'Aspirin'}, {id: 'c2', regno: 123, label: 'Caffeine'}],
  batch: [{id: 'b1', compound_id: 'c1', batch_regno: 45}, {id: 'b2', compound_id: 'c2', batch_regno: 46}],
  note: [{id: 'n1', batch_id: 'b1', text: 'dried'}],
};

scoped('a delete from the record page shows until the save, Restore takes it back, and the save lands on the list',
  async () => {
    const {app: a} = await app();
    await a.goTo('entity', 'i2');
    await flush();
    const offered = allowedActions(a.menuActions(), {access: a.entitySource.value.access.value}).map((x) => x.name);
    assert.equal(offered.includes('Open'), false, 'the row open is not offered again');
    assert.equal(offered.includes('Delete'), true);
    const banner = a.root.querySelector('[data-u2-part="deleting"]');
    assert.equal(banner.hidden, true);
    a.menuActions().find((x) => x.name === 'Delete').run();
    await flush();
    assert.equal(banner.hidden, false, 'the pending delete is on the page, not only in the status bar');
    assert.match(banner.textContent, /Ibuprofen will be deleted when you save/);
    assert.equal(a.root.classList.contains('u2-domain-app-deleting-row'), true);
    fire(banner.querySelector('button'), 'click');
    await flush();
    assert.equal(banner.hidden, true, 'Restore takes it back');
    assert.equal(a.session.isDirty.value, false);
    a.menuActions().find((x) => x.name === 'Delete').run();
    assert.equal(await a.session.save(), true);
    await flush();
    await flush();
    assert.equal(a.page.value, 'list', 'a saved deletion leaves its page for the list');
    assert.equal(a.listSource.rows.byKey('i2'), undefined);
    assert.equal(await a.open(`${BASE}?entity=i2`), true);
    await flush();
    const empty = a.form.value.root.querySelector('[data-u2-part="empty"]').textContent;
    assert.equal(empty, 'This issue does not exist or was deleted.', 'and its address no longer echoes the id');
    a.dispose();
  });

scoped('the list opens in creation order; Sort offers the table\'s orders, the one in force checked', async () => {
  const {table, app: a} = await app();
  assert.deepEqual(DomainApp.sortChoices(table), [['Oldest first', 'created_on,project_id,number'],
    ['Newest first', '!created_on,!project_id,!number'], ['By title', 'title']]);
  assert.equal(a.listSource.sort.value, 'created_on,project_id,number', 'not the table\'s arbitrary order');
  const other = table.source();
  await flush();
  other.newRow({project_id: 'p2', number: 2, title: 'Aardvark'});
  assert.equal(await other.save(), true);
  await flush();
  assert.deepEqual(titles(a), ['Aspirin', 'Ibuprofen', 'Naproxen', 'Aardvark']);
  assert.equal(await a.setSort('!created_on,!project_id,!number'), true);
  await flush();
  assert.equal(titles(a)[0], 'Aardvark', 'newest first: what was just registered is on top');
  const [, , , more] = a.ribbon().main;
  document.body.append(more.root);
  fire(more.root, 'click');
  await flush();
  const labels = [...document.querySelectorAll('[role="menuitem"] .u2-menu-label')].map((el) => el.textContent);
  assert.equal(labels.includes('Sort'), true, 'the orders are a ⋯ submenu');
  fire(document.body, 'keydown', {key: 'Escape'});
  await a.setSort('title');
  await flush();
  assert.deepEqual(titles(a), ['Aardvark', 'Aspirin', 'Ibuprofen', 'Naproxen']);
  other.dispose();
  a.dispose();
  const b = domains.app({table, base: BASE, sort: '!title'});
  document.body.append(b.root);
  await flush();
  assert.deepEqual(titles(b), ['Naproxen', 'Ibuprofen', 'Aspirin', 'Aardvark'], 'an app sets its own');
  b.dispose();
});

scoped('a frame carries its captions in the target\'s format: the schema renderer, a page loaded next', async () => {
  emsBackends.domain = new MemoryDomainBackend(REG, {rows: structuredClone(REG_ROWS)});
  const notes = (await domains.table('reg.note')).source({pageSize: 10});
  await flush();
  const note = notes.rows.byKey('n1');
  assert.equal(note[Rows.caption('batch_id')], 'DGB-000045');
  assert.equal(DomainTable.schemaRenderer(() => notes.schema).caption(note), 'DGB-000045',
    'a row keyed by a numbered reference is named by it');
  notes.dispose();
  const batches = (await domains.table('reg.batch')).source({pageSize: 1});
  await flush();
  await batches.loadMore();
  await flush();
  assert.deepEqual(batches.rows.items.value.map((r) => r[Rows.caption('compound_id')]), ['DG-000007', 'DG-000123'],
    'the page appended later too');
  batches.dispose();
});

scoped('the platform grid reads a numbered key, a date and a reference to a numbered row in their display formats',
  async () => {
    const memory = new MemoryDomainBackend(REG, {rows: structuredClone(REG_ROWS)});
    const columns = [{name: 'id', type: 'string'}, {name: 'batch_regno', type: 'int'}, {name: 'made', type: 'datetime'},
      {name: 'compound_id', type: 'string', semType: 'reg.compound'}, {name: Rows.caption('compound_id'), type: 'string'}];
    emsBackends.domain = {saveAll: (edits) => memory.saveAll(edits), table: async (address) => {
      const t = await memory.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {frame: async (spec) => {
        const df = new DataFrame(columns, (await t.query(spec)).map((r) =>
          Object.fromEntries(columns.map((c) => [c.name, r[c.name] ?? null]))));
        const editor = {table: address, isDirty: false, changeCount: 0, isSaving: false,
          onChanged: new Stream(), onDirtyChanged: new Stream(), onSavingChanged: new Stream(), onSaved: new Stream(),
          onRefused: new Stream(), get dataFrame() { return df; }, stateOf: () => '', errorsOf: () => ({}),
          isChanged: () => false, errorOf: () => null, detach() {}};
        const edit = new EditorEditState(editor);
        return {df, edit, append: async () => 0, dispose: () => edit.dispose()};
      }});
    }};
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', [new Property('allowEdit', 'bool', {defaultValue: false})])];
    try {
      const source = (await domains.table('reg.batch')).source({pageSize: 10});
      const grid = domains.grid(source);
      await flush();
      const df = source.df.value;
      assert.equal(df.columns.byName('batch_regno').getTag('format'), 'DGB-000000', 'DGB-000045, as on the form');
      assert.equal(df.columns.byName('made').getTag('format'), 'yyyy-MM-dd', 'a date in the schema\'s format');
      assert.equal(df.get(Rows.caption('compound_id'), 0), 'DG-000007',
        'the caption the platform\'s reference cell draws, in the target\'s format');
      grid.dispose();
      source.dispose();
    } finally {
      WidgetDescriptor.registry = [];
      platform.reset();
    }
  });

scoped('the left filters pane is folded while it holds no condition; a filter applied from outside opens it',
  async () => {
    const {app: a} = await app({filters: {mode: 'builder', placement: 'left'}});
    const pane = a.root.querySelector('.u2-domain-app-filters-pane');
    const strip = a.root.querySelector('[data-u2-part="filters-strip"]');
    const fold = () => fire(pane.querySelector('.u2-domain-app-filters-head button'), 'click');
    assert.equal(pane.hidden, true, 'an empty pane leaves its width to the list');
    assert.equal(strip.hidden, false);
    assert.equal(a.root.querySelector('.u2-splitter').classList.contains('u2-domain-app-filters-folded'), true);
    fire(strip.querySelector('button'), 'click');
    await flush();
    assert.equal(pane.hidden, false, 'the strip opens it');
    fold();
    await flush();
    assert.equal(pane.hidden, true, 'and its head folds it');
    assert.equal(await a.open(`${BASE}?q=${encodeURIComponent('done = true')}`), true);
    await flush();
    assert.equal(pane.hidden, false, 'a filter from the address shows where it is');
    fold();
    await flush();
    assert.equal(strip.querySelector('.u2-domain-app-filters-count').textContent, '1',
      'folded, the strip says how many conditions are in force');
    a.dispose();
  });

scoped('the row the user moves to is the platform\'s current object; the row the list opens on is not', async () => {
  const writes = grok.shell.dart.writes.length;
  const {table, app: a} = await app();
  assert.equal(a.listSource.currentRow.value.title, 'Aspirin');
  assert.equal(grok.shell.dart.writes.length, writes, 'opening the list moves nothing');
  a.list.list.root.focus();
  a.list.list.selectedIndex.value = 1;
  await flush();
  assert.equal(grok.shell.o.typeName, 'grit.issue', 'the row as the platform\'s DomainRow');
  assert.equal(grok.shell.o.values.title, 'Ibuprofen');
  a.dispose();
  const b = domains.app({table, base: BASE, pageSize: 10, context: (row) => row.number === 1 ? `issue ${row.title}` : null});
  document.body.append(b.root);
  await flush();
  b.list.list.root.focus();
  b.list.list.selectedIndex.value = 2;
  await flush();
  assert.equal(grok.shell.o, 'issue Naproxen', 'what the app maps the row to');
  b.list.list.selectedIndex.value = 1;
  await flush();
  assert.equal(grok.shell.o, 'issue Naproxen', 'null shows nothing new');
  b.dispose();
});

scoped('a child tab\'s New opens the child\'s create page where the app routes it; a draft\'s children stay in the tab',
  async () => {
    emsBackends.domain = backend();
    const [projects, issues] = await Promise.all([domains.table('grit.project'), domains.table('grit.issue')]);
    const issueApp = domains.app({table: issues, base: BASE, pageSize: 10});
    document.body.append(issueApp.root);
    const asked = [];
    const a = domains.app({table: projects, base: '/apps/T/Projects', pageSize: 10, history: false,
      children: {tables: [issues], mode: 'list', onNew: (entry, values) => {
        asked.push([entry.fk, values]);
        DomainApp.createAt(BASE, values);
      }}});
    document.body.append(a.root);
    await flush();
    await a.goTo('entity', 'p1');
    await flush();
    await flush();
    fire(a.panes.querySelector('button[data-u2="new-button"]'), 'click');
    await flush();
    assert.deepEqual(asked, [['project_id', {project_id: 'p1'}]]);
    assert.equal(issueApp.entity.value, DomainApp.NEW, 'the issue app\'s create page');
    assert.equal(issueApp.entitySource.value.currentRow.value.project_id, 'p1', 'with its project filled in');
    assert.equal(a.session.isDirty.value, false, 'and no blank row in the tab');
    await a.create();
    await flush();
    await flush();
    fire(a.panes.querySelector('button[data-u2="new-button"]'), 'click');
    await flush();
    assert.equal(asked.length, 1, 'a draft\'s child is added in the tab, to be saved with it');
    assert.equal(a.session.isDirty.value, true);
    a.session.discard();
    await flush();
    a.dispose();
    issueApp.dispose();
  });

scoped('a count the poll\'s reload asked for never lands over a newer load\'s', async () => {
  const inner = backend();
  const gate = {hold: false, release: undefined};
  emsBackends.domain = {saveAll: (edits) => inner.saveAll(edits), table: async (address) => {
    const t = await inner.table(address);
    return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {count: (scope) => {
      const answer = t.count(scope);
      if (!gate.hold)
        return answer;
      gate.hold = false;
      return new Promise((resolve) => gate.release = () => resolve(answer));
    }});
  }};
  const source = (await domains.table('grit.issue')).source({pageSize: 10});
  await flush();
  assert.equal(source.total.value, 3);
  gate.hold = true;
  const reloading = source._reload();
  await flush();
  assert.equal(typeof gate.release, 'function', 'the reload waits on its count');
  source.query.value = 'done = true';
  await flush();
  assert.equal(source.total.value, 1);
  gate.release();
  await reloading;
  assert.equal(source.total.value, 1, 'the reload\'s count is about a query no longer in force');
  source.dispose();
});

/** `has` answers by substring; the calls are counted. */
async function compounds() {
  const schema = {name: 'chem', tables: {compound: {businessKey: ['code'], columns: {
    code: {type: 'int'}, structure: {type: 'string', semType: 'Molecule'}, label: {type: 'string', searchable: true}}}}};
  emsBackends.domain = new MemoryDomainBackend(schema, {rows: {compound: [
    {id: 'c1', code: 1, structure: 'CCO', label: 'ethanol'}, {id: 'c2', code: 2, structure: 'CC(=O)O', label: 'acetic acid'},
    {id: 'c3', code: 3, structure: 'c1ccccc1O', label: 'phenol'}]}});
  const table = await domains.table('chem.compound');
  const calls = [];
  const off = table.operators.add({id: 'has', label: 'has', semType: 'Molecule', rows: async (c) => {
    calls.push(c);
    return (await emsBackends.domain.tableSync('chem.compound').query({})).filter((r) =>
      r.structure.includes(String(c.value))).map((r) => r.id);
  }});
  return {table, calls, off};
}

scoped('a host condition is answered again only once the table moves: a new order or search reuses its ids',
  async () => {
    const {table, calls, off} = await compounds();
    try {
      const source = table.source({query: 'structure has "CC"'});
      await flush();
      assert.equal(calls.length, 1);
      source.sort.value = '!code';
      await flush();
      source.search.value = 'acid';
      await flush();
      assert.equal(calls.length, 1, 'the order and the search changed, the table did not');
      assert.deepEqual(source.rows.items.value.map((r) => r.id), ['c2']);
      const other = table.source();
      await flush();
      other.newRow({code: 4, structure: 'CCC(=O)O', label: 'propionic acid'});
      assert.equal(await other.save(), true);
      await flush();
      assert.equal(calls.length, 2, 'a write moves the table: the condition is answered again');
      assert.deepEqual(source.rows.items.value.map((r) => r.label).sort(), ['acetic acid', 'propionic acid']);
      other.dispose();
      source.dispose();
    } finally {
      off();
    }
  });

scoped('a refusal over a derived column names the input it is derived from — the field the user can change',
  async () => {
    emsBackends.domain = backend();
    const table = await domains.table('grit.issue');
    table.derived.add({from: ['title'], to: ['description'], derive: (row) => ({description: `${row.title}!`})});
    table.validators.add('description', (value) => String(value).includes('?') ? 'This title is taken' : null);
    const source = table.source({pageSize: 10});
    await flush();
    source.rows.byKey('i1').title = 'Aspirin?';
    await source.settled();
    assert.equal(await source.save(), false);
    assert.equal(source.error.value.message, 'Cannot save: Aspirin?: Title: This title is taken');
    source.dispose();
  });
