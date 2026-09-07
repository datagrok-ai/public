/* DG.Widget extends Control (V-4): the platform doubles ARE u2 controls by inheritance, over the
   kill-walk globals of tests/dg-stub.mjs. Pinned: identity and `Control.is`, the property tier
   over a Dart viewer — the look as `propertyTarget`, the echo suppression and the walkable `table`
   step — the two-site property-change seam, the lifecycle by kind both ways, and the
   ambient-adoption contract: a bare construction is never adopted, `viewerControl`'s build is —
   and the WO-U7 laziness: construction mints no scope; the first use does, wiring the lifecycle
   at mint. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Signal, isWritableSignal} from '../src/core/signals.js';
import {Component, Control} from '../src/core/component.js';
import {isBindSource} from '../src/spec/bind-source.js';
import {Registry} from '../src/spec/registry.js';
import {DartWidget, DataFrame, JsViewer, Property, Viewer, ViewerMetaHelper,
  WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {viewerControl} = await import('../src/dg/viewers/viewer-control.js');

const demog = (rows = 2) => new DataFrame([{name: 'age', type: 'int'}, {name: 'sex'}],
  Array.from({length: rows}, (_, i) => ({age: 30 + i, sex: i % 2 ? 'M' : 'F'})), 'demog');

const descriptors = () => [
  new WidgetDescriptor('Grid', [
    new Property('allowEdit', 'bool', {userEditable: true, category: 'Data', defaultValue: true}),
    new Property('columnNames', 'list', {subType: 'string'}),
    new Property('rowHeight', 'int', {set: null, defaultValue: 20}),
  ]),
  new WidgetDescriptor('Filters', [new Property('filters', 'list', {subType: 'map'})]),
  new WidgetDescriptor('Scatter plot', [new Property('xColumnName', 'string')]),
];

/** Every test runs over the fixture descriptors, a clean platform record and a clean document, and
 * must leave the live-scope count where it was. */
function inherited(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    WidgetDescriptor.registry = descriptors();
    try {
      await body();
    } finally {
      WidgetDescriptor.registry = [];
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

class Gauge extends JsViewer {
  constructor() {
    super();
    this.detaches = 0;
    this.changes = [];
    this.level = this.addProperty('level', 'int', 1);
  }

  onPropertyChanged(p) {
    this.changes.push(p?.name ?? null);
  }

  detach() {
    this.detaches++;
    super.detach();
  }
}

class PoliteGauge extends JsViewer {
  constructor() {
    super();
    this.level = this.addProperty('level', 'int', 1);
  }

  onPropertyChanged(p) {
    super.onPropertyChanged(p);
  }
}

inherited('a platform viewer is a u2 control by inheritance', () => {
  const df = demog();
  const v = Viewer.grid(df);
  assert.ok(v instanceof Viewer && v instanceof Control);
  assert.equal(Control.is(v), true);
  assert.ok(v.root.classList.contains('u2-component'), 'stamped by the Control constructor');
  assert.equal(v.root.classList.contains('u2-viewer'), false, 'u2-viewer arrives from viewerControl');

  assert.equal(v.getProperties, Viewer.prototype.getProperties, 'the platform override wins');
  assert.equal(v.root, v.dart.root);
  assert.equal(v.propertyTier, true);
  for (const member of ['run', 'effect', 'own', 'dispose', 'link', 'bindStep', 'bindProps'])
    assert.equal(typeof v[member], 'function', member);
  assert.equal(v.propertyTarget, globalThis.grok_Viewer_Get_Look(v.dart), 'a viewer\'s properties are over its look');
  assert.ok(v.scope instanceof Scope && v._u2 !== undefined);
  v.dispose();
});

inherited('construction mints no scope; the first use mints one and wires the lifecycle', () => {
  const live = Scope.liveCount;
  const v = Viewer.grid(demog());
  assert.equal(v._u2.scope, undefined, 'a fresh wrapper owns nothing');
  assert.equal(Scope.liveCount, live, 'no scope minted at construction');
  assert.equal(platform.cleanups.length, 0, 'no kill cleanup registered at construction');

  v.bindStep('allowEdit');
  assert.ok(v._u2.scope instanceof Scope);
  assert.equal(Scope.liveCount, live + 1, 'first use mints exactly one scope');
  assert.equal(platform.cleanups.length, 1, 'the Dart lifecycle is wired at mint');
  const seen = [];
  v.own(() => seen.push(platform.kills.length));
  v.dispose();
  assert.deepEqual(platform.kills, [v.root], 'the mint-wired kill fires');
  assert.deepEqual(seen, [1], 'the lifecycle disposer registered first: the kill precedes later owns');
});

inherited('Component.is / Control.is never mint', () => {
  const v = Viewer.grid(demog());
  const live = Scope.liveCount;
  assert.equal(Control.is(v), true);
  assert.equal(Component.is(v), true);
  assert.equal(v._u2.scope, undefined, 'probing a wrapper mints nothing');
  assert.equal(Scope.liveCount, live);
  v.dispose();
});

inherited('the registry stamp lands on componentMeta, apart from the platform meta helper', () => {
  const v = Viewer.grid(demog());
  assert.ok(v.meta instanceof ViewerMetaHelper);
  const reg = new Registry();
  reg.register({tag: 'u2-fake-grid', create: () => v, description: 'Fake', props: [], example: {tag: 'u2-fake-grid'}});
  const built = reg.get('u2-fake-grid').create({});
  assert.equal(built, v);
  assert.equal(v.componentMeta.tag, 'u2-fake-grid');
  assert.ok(v.meta instanceof ViewerMetaHelper, 'the stamp never touches the platform helper');
  assert.deepEqual(v.specProps, {});
  v.dispose();
});

inherited('dispose() kills a Dart-owned viewer through the kill-walk, once', () => {
  const v = Viewer.grid(demog());
  let detached = 0;
  v.onDetached.subscribe(() => detached++);
  assert.equal(v._u2.scope, undefined, 'never engaged — dispose mints, wires and kills');
  v.dispose();
  assert.deepEqual(platform.kills, [v.root]);
  assert.equal(v.scope.isDisposed, true);
  assert.equal(v.isDetached, true);
  assert.equal(detached, 1);
  v.dispose();
  assert.equal(platform.kills.length, 1, 'a second dispose kills nothing');
});

inherited('the platform killing a parent disposes the scope, and u2 adds no kill of its own', () => {
  const v = Viewer.grid(demog());
  let detached = 0;
  v.onDetached.subscribe(() => detached++);
  v.bindStep('allowEdit');
  const parent = document.createElement('div');
  parent.append(v.root);
  globalThis.grok_Widget_Kill(parent);
  assert.equal(v.scope.isDisposed, true);
  assert.deepEqual(platform.kills, [parent], 'the cleanup path kills nothing more');
  assert.equal(detached, 1);
  assert.equal(platform.cleanups.length, 0, 'the cleanup ran once and is gone');
});

inherited('a JS-owned viewer detaches both ways, once', () => {
  const g = new Gauge();
  assert.ok(g instanceof JsViewer && Control.is(g));
  g.dispose();
  assert.equal(g.detaches, 1, 'dispose → the class\'s own detach');
  assert.equal(g.isDetached, true);
  assert.equal(platform.kills.length, 0);

  const h = new Gauge();
  h.effect(() => {});
  h.detach();
  assert.equal(h.scope.isDisposed, true, 'detach → the scope');
  assert.equal(h.detaches, 1, 'no second detach');
  h.dispose();
  assert.equal(h.detaches, 1);

  const k = new Gauge();
  k.effect(() => {});
  const parent = document.createElement('div');
  parent.append(k.root);
  globalThis.grok_Widget_Kill(parent);
  assert.equal(k.scope.isDisposed, true, 'the kill-walk reaches a JS-owned viewer through its detach');
  assert.equal(k.detaches, 1);

  const quiet = new Gauge();
  quiet.detach();
  assert.equal(quiet._u2.scope, undefined, 'detach never mints');
});

inherited('ambient adoption is the factory\'s, not the constructor\'s', () => {
  const s = new Scope();
  const bare = Scope.runWith(s, () => Viewer.grid(demog()));
  const built = Scope.runWith(s, () => viewerControl('Grid', {table: demog()}));
  assert.ok(built.root.classList.contains('u2-viewer'));
  s.dispose();
  assert.equal(bare.scope.isDisposed, false, 'a bare construction is never ambient-adopted');
  assert.equal(built.scope.isDisposed, true, 'viewerControl\'s build re-adds the adoption');
  assert.deepEqual(platform.kills, [built.root]);
  bare.dispose();
});

inherited('propertyTarget is the look for a Dart viewer, the handle for a DartWidget, the widget itself for a JsViewer', () => {
  const v = Viewer.grid(demog());
  assert.equal(v.propertyTarget, v.dart.look);
  assert.equal(v.bindStep('allowEdit').value, true, 'read off the look');

  const w = new DartWidget({n: 1, properties: [Property.create('n', 'int', (d) => d.n, (d, x) => d.n = x, 1)]});
  assert.equal(w.propertyTarget, w.dart);
  const n = w.bindStep('n');
  assert.equal(n.value, 1, 'read off the handle');
  n.value = 2;
  assert.equal(w.dart.n, 2, 'written to the handle');

  const j = new Gauge();
  assert.equal(j.propertyTarget, j);
  assert.equal(j.bindStep('level').value, 1, 'read off the widget');
  v.dispose();
  w.dispose();
  j.dispose();
});

inherited('the property tier reads and writes the look through the props bag', async () => {
  const v = Viewer.grid(demog());
  const events = [];
  v.onPropertyValueChanged.subscribe((e) => events.push(e.args.property.name));
  assert.equal(v.props.allowEdit, true);

  const allowEdit = v.bindStep('allowEdit');
  assert.ok(allowEdit instanceof Signal && isWritableSignal(allowEdit));
  assert.equal(allowEdit.value, true);
  allowEdit.value = false;
  assert.equal(v.props.allowEdit, false, 'written through to the look');
  assert.deepEqual(events, ['allowEdit'], 'one platform event per write');

  const runs = [];
  const s = new Scope();
  s.effect(() => runs.push(allowEdit.value));
  v.props.allowEdit = true;
  await flush();
  assert.deepEqual(runs, [false, true], 'an external write updates the signal once');
  assert.deepEqual(events, ['allowEdit', 'allowEdit'], 'and echoes no write back');
  v.props.allowEdit = true;
  await flush();
  assert.deepEqual(runs, [false, true], 'a same-value event re-fires nobody');

  const rowHeight = v.bindStep('rowHeight');
  assert.ok(rowHeight instanceof Signal && !isWritableSignal(rowHeight), 'no setter → read-only');
  assert.equal(rowHeight.value, 20);
  assert.equal(v.bindStep('nope'), null);
  s.dispose();
  v.dispose();
});

inherited('filters on a FilterGroup are written to the look like any property', async () => {
  const fg = Viewer.filters(demog());
  const events = [];
  fg.onPropertyValueChanged.subscribe((e) => events.push(e.args.property.name));
  const filters = fg.bindStep('filters');
  assert.equal(filters.value, null);
  filters.value = [{type: 'categorical', column: 'sex'}];
  assert.deepEqual(fg.props.filters, [{type: 'categorical', column: 'sex'}]);
  filters.value = [{type: 'histogram', column: 'age'}];
  assert.deepEqual(fg.props.filters, [{type: 'histogram', column: 'age'}], 'the list replaces, never merges');
  assert.deepEqual(events, ['filters', 'filters'], 'one platform event per write');
  await flush();
  assert.deepEqual(filters.value, fg.props.filters, 'the step equals the look after the refresh');
  assert.deepEqual(events, ['filters', 'filters'], 'and the refresh echoes no write back');
  fg.dispose();
});

inherited('the table step is the frame\'s binding surface and follows a repoint', () => {
  const df = demog();
  const v = Viewer.grid(df);
  const table = v.bindStep('table');
  assert.ok(isBindSource(table));
  assert.equal(v.bindStep('table'), table, 'cached');
  const rowCount = table.bindStep('rowCount');
  assert.equal(rowCount.value, 2);
  assert.equal(table.bindStep('currentRow').bindStep('age').value, 30);
  v.dataFrame = demog(5);
  assert.equal(rowCount.value, 5);

  const reads = v.dart.propertyReads;
  const props = v.bindProps();
  assert.deepEqual(props[0], {name: 'table', type: 'dataframe', walkable: true});
  assert.deepEqual(props.slice(1).map((p) => [p.name, p.writable]),
    [['allowEdit', true], ['columnNames', true], ['rowHeight', false]]);
  assert.equal(v._u2.propertySignals.size, 0, 'enumeration allocates no signal');
  v.bindProps();
  assert.equal(v.dart.propertyReads, reads + 1, 'the property list is indexed once');
  v.dispose();
});

inherited('two grids keep their own table step', () => {
  const a = Viewer.grid(demog(2));
  const b = Viewer.grid(demog(5));
  const tableA = a.bindStep('table');
  const tableB = b.bindStep('table');
  assert.notEqual(tableA, tableB);
  assert.equal(tableA.bindStep('rowCount').value, 2);
  assert.equal(tableB.bindStep('rowCount').value, 5);
  a.dataFrame = demog(3);
  assert.equal(tableA.bindStep('rowCount').value, 3);
  assert.equal(tableB.bindStep('rowCount').value, 5, 'a repoint of one leaves the other');
  a.dispose();
  b.dispose();
});

inherited('the property-change seam notifies from both sites', async () => {
  const g = new Gauge();
  const level = g.bindStep('level');
  assert.equal(level.value, 1);
  g.props.level = 7;
  await flush();
  assert.equal(level.value, 7, 'the setter site covers an override that never calls super');
  assert.deepEqual(g.changes, ['level'], 'the class\'s own hook still runs');
  g.dispose();

  const p = new PoliteGauge();
  const step = p.bindStep('level');
  let runs = 0;
  p.effect(() => {
    step.value;
    runs++;
  });
  p.props.level = 5;
  await flush();
  assert.deepEqual([step.value, runs], [5, 2], 'the double notify collapses to one refresh');
  p.level = 9;
  p.onPropertyChanged(null);
  await flush();
  assert.equal(step.value, 9, 'a direct onPropertyChanged call refreshes through the virtual site');
  p.dispose();
});
