/* The four platform-backed sources (WO-14) over FAKE backends — never the platform. What is
   pinned here: a bound param makes a func-source reactive, its outputs read undefined until the
   run lands, the design-time policies decide whether it runs at all and with what cap, paging and
   filtering reset a collection, a table source follows the workspace, and an entity ref answers
   its properties as binding leaves. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {backends} from '../src/sources/backends.js';
import {FuncSource, DESIGN_ROWS} from '../src/sources/func-source.js';
import {EntitySource} from '../src/sources/entity-source.js';
import {TableSource} from '../src/sources/table-source.js';
import {EntityRef} from '../src/sources/entity-ref.js';
import {DataFrame, Entity, Property, Stream, User} from './platform-doubles.mjs';

/** Every test runs over its own fake backends and must leave the live-scope count where it was. */
function source(name, body) {
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

/** Past the AsyncValue debounce a source keeps by default — entity-ref declares no `debounce`
 * prop, so its lookups always ride the 300 ms one. */
const settle = () => new Promise((resolve) => setTimeout(resolve, 400));

function env(options = {}) {
  return {designTime: options.designTime ?? false, subBinds: options.subBinds ?? {},
    resolve: options.resolve ?? (() => null)};
}

function fakeRunner(descriptor, result) {
  return {
    calls: [],
    find(name) {
      return name === descriptor.name ? descriptor : null;
    },
    run(name, params, abort, options) {
      this.calls.push({name, params: {...params}, options});
      return Promise.resolve(result(params, this.calls.length));
    },
  };
}

/** A VISUAL query — the one kind that is known to read, and the only one design-time auto-run and
 * the row cap are allowed to reach. */
const QUERY = {name: 'Orders', kind: 'table-query', inputs: [new Property('days', 'int')],
  outputs: [new Property('rows', 'int')]};
/** Hand-written SQL: `DG.DataQuery` covers DML and DDL, so opening a spec must not run it. */
const SQL_QUERY = {...QUERY, kind: 'query'};
const SCRIPT = {name: 'Send', kind: 'script', inputs: [], outputs: [new Property('sent', 'bool')]};
const FRAME_QUERY = {name: 'Frames', kind: 'query', inputs: [],
  outputs: [new Property('orders', 'dataframe')]};

source('func-source: a bound param makes the run reactive, and outputs read undefined until ready',
  async () => {
    const runner = fakeRunner(QUERY, (params) => ({rows: params.days * 2}));
    backends.funcRunner = runner;
    const days = signal(7);
    const src = new FuncSource({func: 'Orders', debounce: 0},
      env({subBinds: {'params.days': '$.daysInput'}, resolve: () => days}));

    const rows = src.bindStep('rows');
    assert.equal(src.bindStep(''), rows, 'the single output is the default binding');
    assert.equal(rows.value, undefined, 'undefined until the run lands');
    assert.equal(src.bindStep('state').value, 'loading');

    src.start();
    await flush();
    assert.equal(runner.calls.length, 1, 'construction and start coalesce into one run');
    assert.deepEqual(runner.calls[0].params, {days: 7});
    assert.equal(rows.value, 14);
    assert.equal(src.bindStep('state').value, 'ready');

    days.value = 30;
    await flush();
    assert.equal(runner.calls.length, 2, 'a bound param changing re-runs the source');
    assert.deepEqual(runner.calls[1].params, {days: 30});
    assert.equal(rows.value, 60);

    assert.equal(src.bindStep('nope'), null, 'a known function offers only what it declares');
    assert.deepEqual(src.bindProps().map((p) => p.name), ['state', 'error', 'rows']);
    src.dispose();
  });

source('func-source: a DataFrame output walks through dfBindings', async () => {
  const df = new DataFrame([{name: 'name', type: 'string'}], [{name: 'aspirin'}, {name: 'ibuprofen'}]);
  backends.funcRunner = fakeRunner(FRAME_QUERY, () => ({orders: df}));
  const src = new FuncSource({func: 'Frames', debounce: 0}, env());
  const orders = src.bindStep('orders');
  assert.equal(orders.bindStep('currentRow').bindStep('name').value, undefined, 'no frame yet');

  await flush();
  assert.equal(orders.bindStep('rowCount').value, 2);
  const name = orders.bindStep('currentRow').bindStep('name');
  assert.equal(name.value, 'aspirin');
  orders.bindStep('currentRowIdx').value = 1;
  assert.equal(name.value, 'ibuprofen', 'stepping the row index moves the bound field');
  name.value = 'naproxen';
  assert.equal(df.dart.rows[1].name, 'naproxen', 'and editing the field writes the frame');
  assert.deepEqual(src.bindProps().map((p) => p.name), ['state', 'error', 'df', 'currentRowIdx', 'currentRow',
    'selection', 'filter', 'rowCount', 'columns'],
  'a single-table source enumerates the frame\'s own steps — no output name the author never chose');
  assert.equal(src.bindProps().find((p) => p.name === 'df').default, true, 'and a path stopping at it binds the frame');
  assert.equal(src.bindStep('currentRow').bindStep('name'), name,
    'a single-table source answers the frame\'s own steps too');
  assert.equal(src.bindStep('rowCount').value, 2);
  src.dispose();
});

/* A platform `DG.Property` carries name and type on PROTOTYPE getters, so a spread of one is an
   empty object — an output enumerated that way is nameless, and the binding picker shows a blank
   row nothing can be picked from. Plain-object fakes hid this until the designer was driven
   against a real package function; every descriptor here declares getter-backed doubles now. */
source('func-source: outputs enumerate by name even when the platform keeps them on getters', async () => {
  const descriptor = {name: 'Frames', kind: 'query', inputs: [],
    outputs: [new Property('orders', 'dataframe'), new Property('total', 'int')]};
  backends.funcRunner = fakeRunner(descriptor, () => ({}));
  const src = new FuncSource({func: 'Frames', designData: 'schema', debounce: 0},
    env({designTime: true}));
  assert.deepEqual(src.bindProps().map((p) => p.name), ['state', 'error', 'orders', 'total']);
  assert.deepEqual(src.bindProps().map((p) => p.type), ['string', 'string', 'dataframe', 'int']);
  assert.equal(src.bindProps()[2].walkable, true, 'the frame output is still walkable');
  src.dispose();
});

source('func-source: the design-time cap is passed exactly when the designer is up', async () => {
  backends.funcRunner = fakeRunner(QUERY, () => ({rows: 1}));
  const running = new FuncSource({func: 'Orders', debounce: 0}, env());
  await flush();
  assert.equal(backends.funcRunner.calls[0].options, undefined, 'a running form is never capped');
  running.dispose();

  backends.funcRunner = fakeRunner(QUERY, () => ({rows: 1}));
  const designed = new FuncSource({func: 'Orders', debounce: 0}, env({designTime: true}));
  await flush();
  assert.equal(designed.designData, 'live', 'a query defaults to live at design time');
  assert.equal(backends.funcRunner.calls.length, 1,
    'a visual query is the one kind that auto-runs at design time');
  assert.deepEqual(backends.funcRunner.calls[0].options, {maxRows: DESIGN_ROWS});
  designed.dispose();
});

source('func-source: the schema policy runs nothing by itself — an explicit Refresh does',
  async () => {
    backends.funcRunner = fakeRunner(QUERY, () => ({rows: 1}));
    const src = new FuncSource({func: 'Orders', designData: 'schema', debounce: 0},
      env({designTime: true}));
    await flush();
    assert.equal(src.bindStep('state').value, 'idle');
    assert.equal(src.bindStep('rows').value, undefined);
    assert.deepEqual(src.bindProps().map((p) => p.name), ['state', 'error', 'rows'],
      'the outputs still enumerate — that is what schema is for');
    assert.equal(backends.funcRunner.calls.length, 0, 'nothing runs on its own');

    await src.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.equal(backends.funcRunner.calls.length, 1,
      'DD9 forbids the AUTO-run, not the one the user asked for — a Refresh that did nothing is ' +
      'what made the source look broken');
    assert.equal(src.bindStep('state').value, 'ready');
    assert.equal(src.bindStep('rows').value, 1, 'and what it brought back outranks the policy');
    src.dispose();
  });

source('func-source: a refresh runs immediately, whatever the debounce', async () => {
  backends.funcRunner = fakeRunner(SCRIPT, () => ({sent: true}));
  const src = new FuncSource({func: 'Send', designData: 'live', debounce: 5000},
    env({designTime: true}));
  await src.getFunctions().find((f) => f.name === 'refresh').apply();
  assert.equal(backends.funcRunner.calls.length, 1,
    'the debounce coalesces dependency changes — it must not sit on the user\'s own ask');
  assert.equal(src.bindStep('state').value, 'ready');
  src.dispose();
});

source('func-source: a frame\'s columns enumerate once a run — or a sample — has landed',
  async () => {
    backends.funcRunner = fakeRunner(FRAME_QUERY,
      () => ({orders: new DataFrame([{name: 'customer', type: 'string'}], [{customer: 'Bayer'}])}));
    const src = new FuncSource({func: 'Frames', designData: 'schema', debounce: 0},
      env({designTime: true}));
    await flush();
    assert.deepEqual(src.bindStep('orders').bindStep('currentRow').bindProps(), [],
      'no frame, no columns — the picker says so instead of showing an empty branch');

    await src.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.deepEqual(src.bindStep('orders').bindStep('currentRow').bindProps().map((p) => p.name),
      ['customer'], 'a source that HAS run knows its columns, exactly as a table source does');
    src.dispose();
  });

source('func-source: the sample keeps one frame — a re-read must not repoint every binding',
  async () => {
    backends.funcRunner = fakeRunner(FRAME_QUERY, () => ({orders: null}));
    backends.tableFromRows = (rows) => new DataFrame([{name: 'name', type: 'string'}], rows);
    const src = new FuncSource({func: 'Frames', designData: 'sample', sample: [{name: 'row'}]},
      env({designTime: true}));
    const first = src.bindStep('orders').bindStep('df').value;
    // any write anywhere invalidates a dependency-free computed: without the cache the sample is
    // rebuilt here, and every field bound under it is looking at a frame nothing points to
    signal(0).value = 1;
    assert.equal(src.bindStep('orders').bindStep('df').value, first, 'the same frame, read after read');
    src.dispose();
  });

source('func-source: the sample policy builds its outputs through tableFromRows', async () => {
  backends.funcRunner = fakeRunner(FRAME_QUERY, () => ({orders: null}));
  const built = [];
  backends.tableFromRows = (rows) => {
    built.push(rows);
    return new DataFrame([{name: 'name', type: 'string'}], rows);
  };
  const src = new FuncSource({func: 'Frames', designData: 'sample', sample: [{name: 'sample row'}]},
    env({designTime: true}));

  assert.equal(src.bindStep('state').value, 'ready');
  assert.deepEqual(built, [[{name: 'sample row'}]], 'the bare array is the single output\'s rows');
  assert.equal(src.bindStep('orders').bindStep('currentRow').bindStep('name').value, 'sample row');
  assert.equal(backends.funcRunner.calls.length, 0);
  src.dispose();
});

source('func-source: a hand-written query never auto-runs at design time, but refresh does',
  async () => {
    backends.funcRunner = fakeRunner(SQL_QUERY, () => ({rows: 1}));
    const src = new FuncSource({func: 'Orders', debounce: 0}, env({designTime: true}));
    src.start();
    await flush();
    assert.equal(src.designData, 'live', 'the policy is still live — it is the auto-kick that is off');
    assert.equal(backends.funcRunner.calls.length, 0,
      'DG.DataQuery is hand-written SQL: GrokConnect executes DML and DDL as written');

    src.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.equal(backends.funcRunner.calls.length, 1, 'the explicit Refresh is allowed');
    src.dispose();
  });

source('func-source: a side-effectful function never auto-runs at design time, but refresh does',
  async () => {
    backends.funcRunner = fakeRunner(SCRIPT, () => ({sent: true}));
    const src = new FuncSource({func: 'Send', designData: 'live', debounce: 0},
      env({designTime: true}));
    src.start();
    await flush();
    assert.equal(backends.funcRunner.calls.length, 0, 'DD9: never AUTO-run');

    src.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.equal(backends.funcRunner.calls.length, 1, 'the explicit Refresh is allowed');
    assert.equal(src.bindStep('sent').value, true);
    src.dispose();
  });

source('func-source: a failed run says why — the message is a binding step', async () => {
  backends.funcRunner = {find: (name) => name === 'Orders' ? QUERY : null,
    run: () => Promise.reject(new Error('connection refused'))};
  const src = new FuncSource({func: 'Orders', debounce: 0}, env());
  await flush();
  assert.equal(src.bindStep('state').value, 'error');
  assert.equal(src.bindStep('error').value, 'connection refused',
    'without this a failed query and an empty one are the same empty form');
  assert.equal(src.bindStep('rows').value, undefined);
  src.dispose();
});

source('func-source: a field bound before the first run shows nothing, not "undefined"', async () => {
  backends.funcRunner = fakeRunner(FRAME_QUERY, () => new Promise(() => {}));
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-func-source', name: 'orders', props: {func: 'Frames'}}],
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-text-input', name: 'customer', props: {label: 'Customer'},
        bind: {value: '$.orders.currentRow.customer'}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);
  await flush();
  assert.equal(instance.root.querySelector('input').value, '');
  instance.dispose();
});

source('func-source: without a backend it refuses to build', async () => {
  assert.throws(() => new FuncSource({func: 'Orders'}, env()), /no platform backend/);
});

function fakeCollection(items) {
  const collection = {
    pages: [],
    filtered: undefined,
    list(options) {
      collection.pages.push(options);
      const rows = options.filter === undefined ? items :
        items.filter((i) => i.name.includes(options.filter));
      const from = (options.pageNumber ?? 0) * options.pageSize;
      return Promise.resolve(rows.slice(from, from + options.pageSize));
    },
    count: () => Promise.resolve(items.length),
    filter(query) {
      collection.filtered = query;
      return collection;
    },
    order: () => collection,
  };
  return collection;
}

const PEOPLE = [new Entity('alice', {friendlyName: 'Alice M'}), new Entity('bob'), new Entity('carol'),
  new Entity('alicia'), new Entity('dave')];

source('entity-source: pages accumulate, names label them, and a filter change reloads from page 0',
  async () => {
    const collection = fakeCollection(PEOPLE);
    backends.dapi = (name) => name === 'users' ? collection : null;
    const src = new EntitySource({entity: 'users', pageSize: 2}, env());
    src.start();
    await flush();

    assert.equal(src.bindStep('items').value.length, 2);
    assert.deepEqual(src.bindStep('names').value, ['Alice M', 'bob'],
      'friendlyName wins, name is the fallback');
    assert.equal(src.bindStep('total').value, 5);
    assert.equal(src.bindStep(''), src.bindStep('items'), 'items is the default binding');

    src.getFunctions().find((f) => f.name === 'loadMore').apply();
    await flush();
    assert.equal(src.bindStep('items').value.length, 4);
    src.getFunctions().find((f) => f.name === 'loadMore').apply();
    await flush();
    assert.equal(src.bindStep('items').value.length, 5);
    assert.equal(src.bindStep('state').value, 'done');

    src.filter.value = 'ali';
    await flush();
    assert.deepEqual(src.bindStep('names').value, ['Alice M', 'alicia'],
      'a new filter is a new collection, from the first page');
    assert.equal(collection.pages[collection.pages.length - 1].pageNumber, 0);
    src.dispose();
  });

source('entity-source: a bound filter resolves at start and drives the reload', async () => {
  const collection = fakeCollection(PEOPLE);
  backends.dapi = () => collection;
  const query = signal('bob');
  const src = new EntitySource({entity: 'users', pageSize: 10},
    env({subBinds: {filter: '$.queryInput'}, resolve: () => query}));
  src.start();
  await flush();
  assert.deepEqual(src.bindStep('names').value, ['bob']);

  query.value = 'carol';
  await flush();
  assert.deepEqual(src.bindStep('names').value, ['carol']);

  // whatever the user types reaches the server's smart-filter grammar: quotes and backslashes go
  query.value = 'car"ol\\';
  await flush();
  assert.equal(src.filter.value, 'carol');
  assert.equal(collection.pages[collection.pages.length - 1].filter, 'carol');
  src.dispose();
});

source('entity-source: a spec feeds a choice input straight from the tray', async () => {
  backends.dapi = () => fakeCollection(PEOPLE);
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-entity-source', name: 'people', props: {entity: 'users', pageSize: 10}}],
    root: {tag: 'u2-form', name: 'form', children: [
      {tag: 'u2-choice-input', name: 'pick', props: {label: 'User'}, bind: {items: '$.people.names'}},
    ]},
  }, new SpecContext(), reg);
  document.body.append(instance.root);
  await flush();

  const options = [...instance.root.querySelectorAll('option')].map((o) => o.textContent);
  assert.equal(options.includes('Alice M'), true, options.join('|'));
  assert.equal(options.includes('dave'), true, options.join('|'));
  instance.dispose();
});

source('entity-source: an unknown collection and a missing backend both contain', async () => {
  assert.throws(() => new EntitySource({entity: 'users'}, env()), /no platform backend/);

  backends.dapi = () => null;
  const src = new EntitySource({entity: 'nope'}, env());
  src.start();
  await flush();
  assert.equal(src.bindStep('state').value, 'error');
  assert.deepEqual(src.bindStep('items').value, []);
  src.dispose();
});

source('table-source: the frame repoints when the workspace changes', async () => {
  const demog = new DataFrame([{name: 'age', type: 'int'}], [{age: 40}, {age: 55}]);
  const tables = {};
  const changed = new Stream();
  backends.workspace = {
    table: (name) => tables[name] ?? null,
    tableNames: () => Object.keys(tables),
    onTablesChanged: changed,
  };
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(message);
  let src;
  try {
    src = new TableSource({table: 'demog'});
  } finally {
    console.warn = original;
  }
  assert.equal(warnings.length, 1, warnings.join('|'));
  assert.match(warnings[0], /no open table named "demog"/);
  assert.equal(src.bindStep('rowCount').value, 0, 'an absent table contains: no rows, no break');
  const age = src.bindStep('currentRow').bindStep('age');
  assert.equal(age.value, undefined);

  tables.demog = demog;
  changed.fire({});
  assert.equal(src.bindStep('rowCount').value, 2);
  assert.equal(age.value, 40, 'the pre-bound field fills when the table opens');
  assert.equal(src.bindStep('df').value, demog);

  delete tables.demog;
  changed.fire({});
  assert.equal(age.value, undefined, 'and empties when it closes');
  src.dispose();
  assert.equal(changed.count, 0, 'disposal drops the workspace subscription');
});

source('table-source: without a backend it refuses to build', async () => {
  assert.throws(() => new TableSource({table: 'demog'}), /no platform backend/);
});

source('entity-ref: a found entity answers its properties as binding leaves', async () => {
  const alice = new User('alice', {email: 'alice@datagrok.ai'});
  backends.dapiFind = (collection, id) =>
    Promise.resolve(collection === 'users' && id === 'u1' ? alice : null);
  const src = new EntityRef({entityType: 'users', id: 'u1'}, env());
  const name = src.bindStep('name');
  assert.equal(name.value, undefined, 'undefined until the lookup lands');

  await settle();
  assert.equal(src.bindStep('state').value, 'ready');
  assert.equal(name.value, 'alice');
  assert.equal(src.bindStep('email').value, 'alice@datagrok.ai');
  assert.equal(src.bindStep(''), src.bindStep('entity'));
  const props = src.bindProps().map((p) => p.name);
  assert.deepEqual(props.slice(0, 3), ['state', 'error', 'entity']);
  assert.equal(props.includes('email'), true, 'prototype getters enumerate; `_`-prefixed names never do');
  assert.equal(props.includes('dart'), false, 'the raw Dart handle an entity keeps as its one own field is no leaf');
  src.dispose();
});

source('entity-ref: a missing entity ends in the error state and its leaves read undefined',
  async () => {
    backends.dapiFind = (collection, entityId) =>
      Promise.resolve(entityId === 'u2' ? new Entity('u2') : null);
    const id = signal('u1');
    const src = new EntityRef({entityType: 'users'}, env({subBinds: {id: '$.idInput'},
      resolve: () => id}));
    src.start();
    await settle();
    assert.equal(src.bindStep('state').value, 'error');
    assert.equal(src.bindStep('name').value, undefined);
    assert.deepEqual(src.bindProps().map((p) => p.name), ['state', 'error', 'entity']);
    assert.match(src.bindStep('error').value, /was not found/,
      'the failure is reachable — an empty form says why it is empty');

    id.value = 'u2';
    await settle();
    assert.equal(src.bindStep('state').value, 'ready');
    assert.equal(src.bindStep('name').value, 'u2', 'a bound id repoints the lookup');
    src.dispose();
  });

source('entity-ref: the schema policy never looks anything up', async () => {
  let calls = 0;
  backends.dapiFind = () => {
    calls++;
    return Promise.resolve(new Entity('alice'));
  };
  const src = new EntityRef({entityType: 'users', id: 'u1', designData: 'schema'},
    env({designTime: true}));
  await settle();
  assert.equal(calls, 0);
  assert.equal(src.bindStep('state').value, 'idle');
  src.dispose();
});
