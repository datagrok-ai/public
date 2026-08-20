/* WO-15 — the binding picker's model: what the tree offers, in what order, and what a leaf's path
   assembles to. The laziness contract is the load-bearing assertion here — enumerating a frame's
   columns must not allocate a single column signal — so the fake frame counts its reads. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Component} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame} from './platform-doubles.mjs';
import {brokenCount} from '../src/dg/designer/node-ref.js';
import {bindTree, NOTHING_YET} from '../src/dg/designer/bind-model.js';

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
