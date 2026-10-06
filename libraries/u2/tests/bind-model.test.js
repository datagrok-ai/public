/* WO-15 — the binding picker's model: what the tree offers, in what order, and what a leaf's path
   assembles to. The laziness contract is the load-bearing assertion here — enumerating a frame's
   columns must not allocate a single column signal — so the fake frame counts its reads. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Component} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {backends} from '../src/sources/backends.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';
import {brokenCount} from '../src/dg/designer/node-ref.js';
import {bindTree, NOTHING_YET} from '../src/dg/designer/bind-model.js';

register('./dg-stub.mjs', import.meta.url);
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');

/** A frame whose handle counts its reads: allocating one column signal costs a `get`. */
const frame = () => new DataFrame([{name: 'name', type: 'string', semType: 'Text'},
  {name: 'Mol Weight', type: 'double', semType: 'Molecule Weight'}], [{'name': 'aspirin', 'Mol Weight': 180.16}]);

/** A tray source over a frame — what a WO-14 table source will be, reduced to the protocol. */
class TableSource extends Component {
  constructor(df, scope) {
    super();
    this._df = dfBindings(signal(df), scope);
  }

  bindStep(name) {
    return this._df.bindStep(name);
  }

  bindProps() {
    return this._df.bindProps();
  }
}

const BROKEN = {
  tag: 'u2-e-broken-source',
  visual: false,
  createComponent: () => {
    throw new Error('no platform backend');
  },
  description: 'A source that cannot be built',
  props: [],
  example: {tag: 'u2-e-broken-source', name: 'nope'},
};

const SPEC = {
  $schema: 'dg-ui/1',
  components: [
    {tag: 'u2-state', name: 'draft', props: {type: 'string', initial: 'seed'}},
    {tag: 'u2-e-broken-source', name: 'ghost'},
  ],
  root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}},
    {tag: 'div', name: 'note', props: {text: 'plain html'}},
  ]},
};

function model(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    const scope = new Scope();
    const reg = new Registry();
    registerAll(reg);
    reg.register(BROKEN);
    const df = frame();
    const orders = Scope.runWith(scope, () => new TableSource(df, scope));
    const ctx = new SpecContext({data: {reagent: 'Ethanol', orders}});
    const instance = renderSpec(JSON.parse(JSON.stringify(SPEC)), ctx, reg);
    try {
      await body({instance, df, roots: bindTree(instance)});
    } finally {
      instance.dispose();
      scope.dispose();
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** A branch, by the name its label starts with. */
const at = (nodes, name) => nodes.find((n) => n.label === name || n.label.startsWith(`${name} `));

model('roots: context data first, then the tray, then the named controls', ({instance, roots}) => {
  assert.deepEqual(roots.map((n) => n.label),
    ['reagent : string ⇄', 'orders', 'draft (u2-state)', 'nameInput (u2-text-input)']);
  assert.deepEqual(roots.map((n) => n.path), ['$.reagent', '$.orders', '$.draft', '$.nameInput'],
    'every root is bindable on its own — a source answers its default step');
  assert.equal(at(roots, 'reagent').children, undefined, 'a plain data signal is a leaf');
  assert.equal(brokenCount(instance), 1, 'the tray carries a component that could not be built');
  assert.equal(roots.some((n) => n.label.startsWith('ghost')), false,
    'and it is skipped, not walked — its map entry is the placeholder element');
  assert.equal(roots.some((n) => n.label.startsWith('note')), false, 'a plain HTML node is no bind source');
  assert.equal(roots.some((n) => n.label.startsWith('form')), false, 'nor is a container without signals');
});

model('a component expands through bindProps, and a control offers its signal-backed props',
  ({roots}) => {
    const draft = at(roots, 'draft').children();
    assert.deepEqual(draft.map((n) => [n.label, n.path]), [['value : string ⇄', '$.draft.value']]);

    const input = at(roots, 'nameInput').children();
    const value = at(input, 'value');
    assert.equal(value.path, '$.nameInput.value');
    assert.equal(value.label.endsWith('⇄'), true, 'a writable leaf carries the two-way affordance');
  });

model('a walkable step is a group with no path, and its leaves assemble bracketed paths',
  ({roots}) => {
    const steps = at(roots, 'orders').children();
    assert.deepEqual(steps.map((n) => n.label), ['df : dataframe', 'currentRowIdx : int ⇄',
      'currentRow : object', 'selection : object', 'filter : object', 'rowCount : int',
      'columns : string_list']);
    const row = at(steps, 'currentRow');
    assert.equal(row.path, null, 'a step that is nothing on its own commits no path');
    assert.equal(at(steps, 'rowCount').path, '$.orders.rowCount');

    const columns = row.children();
    assert.deepEqual(columns.map((n) => [n.label, n.path]), [
      ['name : string · Text ⇄', '$.orders.currentRow.name'],
      ['Mol Weight : double · Molecule Weight ⇄', '$.orders.currentRow[\'Mol Weight\']'],
    ]);
  });

model('enumeration allocates nothing: walking the whole tree never reads the frame', ({df, roots}) => {
  const walk = (nodes, depth) => {
    for (const node of nodes) {
      if (node.children && depth > 0)
        walk(node.children(), depth - 1);
    }
  };
  walk(roots, 4);
  assert.equal(df.dart.reads, 0, 'a column signal is created by bindStep, never by bindProps');

  const columns = at(at(roots, 'orders').children(), 'currentRow');
  assert.equal(columns.path, null);
  assert.equal(df.dart.reads, 0);
});

/* The acceptance blocker: a source that has not run knows no columns, and an expander that opens
   onto nothing is indistinguishable from a broken one — "I genuinely thought my data source was
   broken". Whatever cannot be enumerated yet says why. */
test('a walkable step with nothing under it explains itself, and the note is not pickable',
  async () => {
    const live = Scope.liveCount;
    const scope = new Scope();
    const reg = new Registry();
    registerAll(reg);
    const empty = Scope.runWith(scope, () => new TableSource(undefined, scope));
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'form'}},
      new SpecContext({data: {orders: empty}}), reg);
    try {
      const row = at(at(bindTree(instance), 'orders').children(), 'currentRow');
      const children = row.children();
      assert.deepEqual(children.map((n) => n.label), [NOTHING_YET]);
      assert.equal(children[0].path, null, 'the picker only commits a path — a note cannot be picked');
    } finally {
      instance.dispose();
      scope.dispose();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });

/* WO-V5 — a viewer node is a bind source with no picker code of its own: `bindProps()` answers the
   walkable `table` step ahead of the look properties, and a property with a setter is a ⇄ leaf. */
test('a viewer node walks through table to the seven frame steps; its writable look props carry ⇄',
  async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', [
      new Property('allowEdit', 'bool'), new Property('rowHeight', 'int', {set: null})])];
    const scope = new Scope();
    const reg = new Registry();
    registerPlatformComponents(reg);
    const df = frame();
    const orders = Scope.runWith(scope, () => new TableSource(df, scope));
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
      {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.orders'}}]}},
    new SpecContext({data: {orders}}), reg);
    try {
      assert.equal(brokenCount(instance), 0);
      const grid = at(bindTree(instance), 'grid');
      assert.equal(grid.label, 'grid (u2-viewer-grid)');
      assert.equal(grid.path, null, 'a viewer has no default step: a group to walk, not a pick');
      const props = grid.children();
      assert.deepEqual(props.map((n) => [n.label, n.path]), [
        ['table : dataframe', null],
        ['allowEdit : bool ⇄', '$.grid.allowEdit'],
        ['rowHeight : int', '$.grid.rowHeight'],
      ], 'the table step first, then every look property — ⇄ exactly where a setter exists');

      const steps = at(props, 'table').children();
      assert.deepEqual(steps.map((n) => n.label), ['df : dataframe', 'currentRowIdx : int ⇄',
        'currentRow : object', 'selection : object', 'filter : object', 'rowCount : int',
        'columns : string_list'], 'the same seven steps a table source answers');
      assert.equal(at(steps, 'rowCount').path, '$.grid.table.rowCount');
      assert.deepEqual(at(steps, 'currentRow').children().map((n) => n.path),
        ['$.grid.table.currentRow.name', '$.grid.table.currentRow[\'Mol Weight\']']);
      assert.equal(df.dart.reads, 0, 'enumeration through the viewer allocates nothing either');
    } finally {
      instance.dispose();
      scope.dispose();
      WidgetDescriptor.registry = [];
      platform.reset();
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });

/* U1 of the viewers acceptance pass: the picker offered a function source as `orders › orders › df`,
   and OK on the source's own row did nothing — the source IS its one table (func-source's
   bindStep), so its row is the pick (`$.orders`, the sample's own path) and the frame's steps sit
   right under it. A function with several outputs keeps their names; neither probes a source to
   find out — `default` is metadata. */
test('a single-table function source is pickable as itself and lists the frame\'s steps; several outputs keep their names',
  async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    const saved = {...backends};
    const outputs = {
      Frames: [new Property('orders', 'dataframe')],
      Both: [new Property('orders', 'dataframe'), new Property('total', 'int')],
    };
    backends.funcRunner = {
      find: (name) => name in outputs ? {name, kind: 'query', inputs: [], outputs: outputs[name]} : null,
      run: () => Promise.resolve({}),
    };
    const reg = new Registry();
    registerAll(reg);
    const instance = renderSpec({$schema: 'dg-ui/1', components: [
      {tag: 'u2-func-source', name: 'orders', props: {func: 'Frames', auto: false}},
      {tag: 'u2-func-source', name: 'multi', props: {func: 'Both', auto: false}},
    ], root: {tag: 'u2-form', name: 'form'}}, new SpecContext(), reg);
    try {
      const roots = bindTree(instance);
      const orders = at(roots, 'orders');
      assert.equal(orders.path, '$.orders', 'the source itself is the pick — its default step is the frame');
      const steps = orders.children();
      assert.deepEqual(steps.map((n) => n.label), ['state : string', 'error : string', 'df : dataframe',
        'currentRowIdx : int ⇄', 'currentRow : object', 'selection : object', 'filter : object',
        'rowCount : int', 'columns : string_list'], 'the frame\'s steps right under the source');
      assert.equal(at(steps, 'df').path, '$.orders.df');
      assert.equal(at(steps, 'currentRow').path, null);

      const multi = at(roots, 'multi');
      assert.equal(multi.path, null, 'no default output — a group that needs one more step');
      const named = multi.children();
      assert.deepEqual(named.map((n) => [n.label, n.path]), [['state : string', '$.multi.state'],
        ['error : string', '$.multi.error'], ['orders : dataframe', null], ['total : int', '$.multi.total']]);
      assert.deepEqual(at(named, 'orders').children().map((n) => n.path).slice(0, 2),
        ['$.multi.orders.df', '$.multi.orders.currentRowIdx'], 'a named output walks as a frame');
    } finally {
      instance.dispose();
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
