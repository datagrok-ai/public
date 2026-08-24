/* adopt(): a js-api widget made a u2 Control in place (VP-6, VP-7) over the getter-backed doubles
   of tests/platform-doubles.mjs and the kill-walk globals of tests/dg-stub.mjs. Pinned: identity
   and the per-class mixed prototype, the platform members winning, the `meta` accessor, the
   lifecycle by kind both ways, the ambient-scope adoption, and the property tier over a Dart viewer
   — the look as `propertyTarget`, `filters` written to the look like any property, the echo
   suppression and the walkable `table` step. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Signal, isWritableSignal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {isBindSource} from '../src/spec/bind-source.js';
import {Registry} from '../src/spec/registry.js';
import {DartWidget, DataFrame, FilterGroup, Grid, JsViewer, Property, Viewer, ViewerMetaHelper,
  WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {adopt, dartOwned} = await import('../src/dg/viewers/adopt.js');

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
function adopted(name, body) {
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

adopted('the adopted widget is the same object, under one mixed prototype per js-api class', () => {
  const df = demog();
  const grid = Viewer.grid(df);
  const proto = Object.getPrototypeOf(grid);
  const v = adopt(grid);
  assert.equal(v, grid);
  assert.ok(v instanceof Grid && v instanceof Viewer);
  assert.equal(Control.is(v), true);
  assert.equal(adopt(v), v, 'idempotent');
  assert.equal(Object.getPrototypeOf(Object.getPrototypeOf(v)), proto, 'mixed in between, the class untouched');
  const twin = adopt(Viewer.grid(df));
  assert.equal(Object.getPrototypeOf(twin), Object.getPrototypeOf(v), 'two grids share one mixed prototype');
  const filters = adopt(Viewer.filters(df));
  assert.notEqual(Object.getPrototypeOf(filters), Object.getPrototypeOf(v), 'a FilterGroup gets its own');
  assert.ok(v.root.classList.contains('u2-component') && v.root.classList.contains('u2-viewer'));
  assert.equal(dartOwned(v), true);

  assert.equal(v.getProperties, Viewer.prototype.getProperties, 'the platform member wins');
  assert.equal(v.onEvent, Viewer.prototype.onEvent);
  assert.equal(v.root, grid.dart.root);
  assert.equal(v.propertyTier, true);
  for (const member of ['run', 'effect', 'own', 'dispose', 'link', 'bindStep', 'bindProps'])
    assert.equal(typeof v[member], 'function', member);
  assert.equal(v.propertyTarget, globalThis.grok_Viewer_Get_Look(v.dart), 'a viewer\'s properties are over its look');
  assert.ok(v.scope instanceof Scope && v._u2 !== undefined);

  v.dispose();
  twin.dispose();
  filters.dispose();
});

adopted('meta answers the platform helper before the stamp and the registry stamp after', () => {
  const v = adopt(Viewer.grid(demog()));
  assert.ok(v.meta instanceof ViewerMetaHelper);
  const reg = new Registry();
  reg.register({tag: 'u2-fake-grid', create: () => v, description: 'Fake', props: [], example: {tag: 'u2-fake-grid'}});
  const built = reg.get('u2-fake-grid').create({});
  assert.equal(built, v);
  assert.equal(v.meta.tag, 'u2-fake-grid');
  assert.deepEqual(v.specProps, {});
  v.dispose();
});

adopted('dispose() kills a Dart-owned viewer through the kill-walk, once', () => {
  const v = adopt(Viewer.grid(demog()));
  let detached = 0;
  v.onDetached.subscribe(() => detached++);
  v.dispose();
  assert.deepEqual(platform.kills, [v.root]);
  assert.equal(v.scope.isDisposed, true);
  assert.equal(v.isDetached, true);
  assert.equal(detached, 1);
  v.dispose();
  assert.equal(platform.kills.length, 1, 'a second dispose kills nothing');
});

adopted('the platform killing a parent disposes the scope, and u2 adds no kill of its own', () => {
  const v = adopt(Viewer.grid(demog()));
  let detached = 0;
  v.onDetached.subscribe(() => detached++);
  const parent = document.createElement('div');
  parent.append(v.root);
  globalThis.grok_Widget_Kill(parent);
  assert.equal(v.scope.isDisposed, true);
  assert.deepEqual(platform.kills, [parent], 'the cleanup path kills nothing more');
  assert.equal(detached, 1);
  assert.equal(platform.cleanups.length, 0, 'the cleanup ran once and is gone');
});

adopted('a JS-owned widget detaches both ways, once', () => {
  const g = adopt(new Gauge());
  assert.equal(dartOwned(g), false);
  assert.ok(g.root.classList.contains('u2-viewer'));
  g.dispose();
  assert.equal(g.detaches, 1, 'dispose → the class\'s own detach');
  assert.equal(g.isDetached, true);
  assert.equal(platform.kills.length, 0);

  const h = adopt(new Gauge());
  h.detach();
  assert.equal(h.scope.isDisposed, true, 'detach → the scope');
  assert.equal(h.detaches, 1, 'no second detach');
  h.dispose();
  assert.equal(h.detaches, 1);

  const k = adopt(new Gauge());
  const parent = document.createElement('div');
  parent.append(k.root);
  globalThis.grok_Widget_Kill(parent);
  assert.equal(k.scope.isDisposed, true, 'the kill-walk reaches a JS-owned widget through its detach');
  assert.equal(k.detaches, 1);
});

adopted('adopting inside a scope hands the widget to it', () => {
  const s = new Scope();
  const v = Scope.runWith(s, () => adopt(Viewer.grid(demog())));
  assert.equal(v.scope.isDisposed, false);
  s.dispose();
  assert.equal(v.scope.isDisposed, true);
  assert.deepEqual(platform.kills, [v.root]);
});

adopted('without the kill globals a Dart-owned widget is adopted, warned about once, and never killed', () => {
  const kill = globalThis.grok_Widget_Kill;
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    delete globalThis.grok_Widget_Kill;
    const v = adopt(Viewer.grid(demog()));
    const w = adopt(Viewer.grid(demog()));
    assert.equal(warnings.filter((m) => m.includes('grok_Widget_Kill')).length, 1);
    v.dispose();
    w.dispose();
    assert.equal(platform.kills.length, 0);
  } finally {
    globalThis.grok_Widget_Kill = kill;
    console.warn = warn;
  }
});

adopted('propertyTarget is the look for a Dart viewer, the handle for a DartWidget, the widget itself for a JsViewer', () => {
  const v = adopt(Viewer.grid(demog()));
  assert.equal(v.propertyTarget, v.dart.look);
  assert.equal(v.bindStep('allowEdit').value, true, 'read off the look');

  const w = adopt(new DartWidget({n: 1, properties: [Property.create('n', 'int', (d) => d.n, (d, x) => d.n = x, 1)]}));
  assert.equal(w.propertyTarget, w.dart);
  const n = w.bindStep('n');
  assert.equal(n.value, 1, 'read off the handle');
  n.value = 2;
  assert.equal(w.dart.n, 2, 'written to the handle');

  const j = adopt(new Gauge());
  assert.equal(j.propertyTarget, j);
  assert.equal(j.bindStep('level').value, 1, 'read off the widget');
  v.dispose();
  w.dispose();
  j.dispose();
});

adopted('without grok_Viewer_Get_Look a Dart viewer is adopted, warned about once, and keeps the default target', () => {
  const look = globalThis.grok_Viewer_Get_Look;
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  class Lone extends Grid {}
  const lone = () => new Lone({type: 'Grid', descriptor: WidgetDescriptor.getByName('Grid'), dataFrame: demog()});
  try {
    delete globalThis.grok_Viewer_Get_Look;
    const v = adopt(lone());
    const w = adopt(lone());
    assert.equal(warnings.filter((m) => m.includes('grok_Viewer_Get_Look')).length, 1);
    assert.equal(v.propertyTarget, v);
    assert.equal(w.propertyTarget, w);
    v.dispose();
    w.dispose();
  } finally {
    globalThis.grok_Viewer_Get_Look = look;
    console.warn = warn;
  }
});

adopted('the property tier reads and writes the look through the props bag', async () => {
  const v = adopt(Viewer.grid(demog()));
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

adopted('filters on a FilterGroup are written to the look like any property', async () => {
  const fg = adopt(Viewer.filters(demog()));
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

adopted('the table step is the frame\'s binding surface and follows a repoint', () => {
  const df = demog();
  const v = adopt(Viewer.grid(df));
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

adopted('two grids under the one mixed prototype keep their own table step', () => {
  const a = adopt(Viewer.grid(demog(2)));
  const b = adopt(Viewer.grid(demog(5)));
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

adopted('a JS-owned widget\'s properties announce through onPropertyChanged', async () => {
  const g = adopt(new Gauge());
  const level = g.bindStep('level');
  assert.equal(level.value, 1);
  level.value = 3;
  assert.equal(g.level, 3);
  assert.deepEqual(g.changes, ['level'], 'the class\'s own hook still runs');
  g.props.level = 7;
  await flush();
  assert.equal(level.value, 7, 'the tier follows an outside write');
  g.dispose();
});
