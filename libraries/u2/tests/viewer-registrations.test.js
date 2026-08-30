/* registerPlatformViewers / registerPlatformComponents (VP-9) over a private registry and the
   WO-V3 descriptor doubles, a spec through the viewer tags (bind, two-way look props, events
   through the viewer's own bus, dump, the Design/Run re-creation), and the viewer sample rendering
   through the `sample` policy (VP-15/16). The core registry stays platform-free: the global one
   after `registerAll()` alone carries no viewer tag, which is what keeps the manifest gate green. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Registry, registry as globalRegistry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {backends} from '../src/sources/backends.js';
import {brokenCount} from '../src/dg/designer/node-ref.js';
import {DataFrame, Grid, Property, Viewer, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {registerPlatformViewers, registerPlatformComponents, kebab, toSpecProp, VIEWER_USAGE} =
  await import('../src/dg/viewers/registrations.js');
const {VIEWER_SAMPLES, platformSamples} = await import('../src/dg/viewers/samples.js');

const descriptors = () => [
  new WidgetDescriptor('Grid', [
    new Property('table', 'dataframe'),
    new Property('viewer', 'object'),
    new Property('allowEdit', 'bool', {userEditable: true, category: 'Data', defaultValue: false, description: 'Editable cells',
      nullable: false}),
    new Property('columnNames', 'list', {subType: 'string'}),
    new Property('filters', 'list', {subType: 'map', userEditable: null}),
    new Property('controlsFont', 'string', {userEditable: false}),
    new Property('rowHeight', 'int', {set: null}),
    new Property('onCellClick', 'string', {category: 'Events'}),
    new Property('initializationFunction', 'string', {category: 'Data'}),
    new Property('onInitializedScript', 'string', {category: 'Data'}),
    new Property('onNodeExpandFunction', 'string'),
    new Property('look', 'lookandfeel'),
    new Property('defaultCellStyle', 'gridcellstyle'),
  ], {events: [{name: 'd4-grid-cell-click', eventName: 'OnCellClicked'}]}),
  new WidgetDescriptor('_makeInspectorPanel', [new Property('xColumnName', 'string')]),
  new WidgetDescriptor('Filters', [new Property('filters', 'list', {subType: 'map'}),
    new Property('columnNames', 'list', {subType: 'string'})]),
  new WidgetDescriptor('Form', []),
  new WidgetDescriptor('Scatter plot', [
    new Property('xColumnName', 'string'), new Property('yColumnName', 'string'),
    new Property('markerMinSize', 'num', {min: 1, max: 50}),
    new Property('markerType', 'string', {choices: ['circle', 'square']}),
  ], {description: 'Dots'}),
  new WidgetDescriptor('3D scatter plot', [new Property('xColumnName', 'string')]),
];

const rows = (n) => Array.from({length: n}, (_, i) =>
  ({orderId: 1001 + i, customer: `C${i}`, city: 'Kyiv', total: 100 * i}));
const frame = (records) => new DataFrame(Object.keys(records[0] ?? {}).map((name) => ({name})), records, 'orders');

/** The platform seam as the sample needs it: `demoOrders` answers a 4-row frame when run, and
 * sample rows become a frame without a run. */
function fakeBackends() {
  const runner = {
    calls: [],
    find: (name) => name === 'demoOrders' ?
      {name, kind: 'function', inputs: [{name: 'days', type: 'int'}], outputs: [{name: 'orders', type: 'dataframe'}]} :
      null,
    run(name, params) {
      this.calls.push({name, params});
      return Promise.resolve({orders: frame(rows(4))});
    },
  };
  backends.funcRunner = runner;
  backends.tableFromRows = (records) => frame(records);
  return runner;
}

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    WidgetDescriptor.registry = descriptors();
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      WidgetDescriptor.registry = [];
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const settle = () => new Promise((resolve) => setTimeout(resolve, 200));

scoped('kebab', () => {
  assert.deepEqual(['Grid', 'Scatter plot', '3D scatter plot', 'PC Plot', 'Tile Viewer', 'Heat map', 'Heatmap'].map(kebab),
    ['grid', 'scatter-plot', '3d-scatter-plot', 'pc-plot', 'tile-viewer', 'heat-map', 'heatmap']);
});

scoped('every descriptor becomes a u2-viewer-* tag with table, look and its user-editable properties', () => {
  const reg = new Registry();
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    registerPlatformViewers(reg);
  } finally {
    console.warn = warn;
  }
  assert.deepEqual(reg.metas().map((m) => m.tag),
    ['u2-viewer-grid', 'u2-viewer-filters', 'u2-viewer-form', 'u2-viewer-scatter-plot', 'u2-viewer-3d-scatter-plot'],
    'a package-internal descriptor (_makeInspectorPanel) is no viewer');
  assert.deepEqual(warnings, [], 'and passing it over is silent');
  const grid = reg.get('u2-viewer-grid');
  assert.equal(grid.category, 'Viewers');
  assert.equal(grid.label, 'Grid');
  assert.equal(grid.usage, VIEWER_USAGE.Grid);
  assert.deepEqual(grid.props.map((p) => p.name), ['table', 'look', 'allowEdit', 'columnNames', 'filters', 'rowHeight'],
    'controlsFont (userEditable false), the descriptor\'s own table/viewer, the Events category, the script hooks ' +
    '(initializationFunction, onInitializedScript, onNodeExpandFunction) and the Dart-internal look structures ' +
    'absent; filters (null) present');
  assert.deepEqual(grid.props[0], {name: 'table', type: 'dataframe', bindable: true,
    description: 'The DataFrame this viewer shows — bind it to a source (`$.orders`).'}, 'bind-only, never required');
  assert.equal(grid.props[1].type, 'object');
  assert.deepEqual(grid.props[2], {name: 'allowEdit', type: 'bool', bindable: true, twoWay: true,
    category: 'Data', description: 'Editable cells'}, 'nullable is not copied');
  assert.equal(grid.props[3].type, 'string_list');
  assert.equal(grid.props[4].type, 'object');
  assert.equal(grid.props[5].twoWay, false, 'twoWay iff the property has a setter');
  assert.ok(grid.props.every((p) => !('propertyType' in p)), 'the platform type never reaches checkProp');
  assert.deepEqual(grid.events, ['d4-data-event', 'd4-grid-cell-click']);
  assert.deepEqual(grid.example, {tag: 'u2-viewer-grid', bind: {table: '$.orders'}});
  assert.equal(grid.designerActions, undefined);

  const plot = reg.get('u2-viewer-scatter-plot');
  assert.equal(plot.description, 'Dots');
  assert.equal(plot.label, 'Scatter plot');
  const glyph = document.createElement('i');
  WidgetDescriptor.prototype.createIcon = function() { return this.name === 'Scatter plot' ? glyph : null; };
  try {
    assert.equal(plot.icon(), glyph, 'the platform\'s own icon, made on demand');
  } finally {
    delete WidgetDescriptor.prototype.createIcon;
  }
  assert.deepEqual(plot.props.find((p) => p.name === 'markerMinSize'),
    {name: 'markerMinSize', type: 'double', bindable: true, twoWay: true, min: 1, max: 50});
  assert.deepEqual(plot.props.find((p) => p.name === 'markerType').choices, ['circle', 'square']);
  assert.deepEqual(plot.events, ['d4-data-event']);
  assert.ok(plot.props.slice(2).every((p) => p.bindable), 'every property bindable; look is the literal escape hatch');
  assert.equal(reg.get('u2-viewer-filters').designerActions[0].name, 'Add filter for column…');
  assert.equal(reg.get('u2-viewer-filters').designerActions[0].icon, 'filter');
  assert.equal(reg.get('u2-viewer-form').props.length, 2);

  registerPlatformViewers(reg);
  assert.equal(reg.metas().length, 5, 'a second call is a no-op');
  const manifest = JSON.stringify(reg.manifest());
  assert.ok(manifest.includes('u2-viewer-grid') && manifest.includes('"label":"Scatter plot"') && !manifest.includes('"icon"'),
    'the label is plain JSON and stays; the icon closure is stripped');
});

scoped('toSpecProp reads field by field and never spreads', () => {
  const p = new Property('x', 'num', {min: 0, friendlyName: 'X', nullable: true, choices: []});
  assert.deepEqual(toSpecProp(p), {name: 'x', type: 'double', bindable: true, twoWay: false,
    friendlyName: 'X', min: 0});
  assert.deepEqual(Object.keys({...p}), ['dart']);
});

scoped('a taken tag is skipped with one warning', () => {
  WidgetDescriptor.registry.push(new WidgetDescriptor('Scatter Plot', []), new WidgetDescriptor('Scatter plot', []));
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    const reg = new Registry();
    registerPlatformViewers(reg);
    assert.equal(reg.metas().length, 5);
    assert.deepEqual(warnings, ['u2 spec: u2-viewer-scatter-plot: "Scatter Plot" skipped — the tag is taken by "Scatter plot"',
      'u2 spec: u2-viewer-scatter-plot: "Scatter plot" is already registered']);
  } finally {
    console.warn = warn;
  }
});

scoped('registerPlatformComponents is the core registry plus the viewers; the core registry alone has none', () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  const core = reg.metas().filter((m) => !m.tag.startsWith('u2-viewer-'));
  assert.equal(core.length, 40);
  assert.equal(reg.metas().length, 45);
  registerPlatformComponents(reg);
  assert.equal(reg.metas().length, 45);

  const fch = reg.get('u2-func-call-history-browser');
  assert.equal(fch.props.find((p) => p.name === 'functionName').bindable, true);
  const fci = reg.get('u2-func-call-input');
  assert.equal(fci.props.find((p) => p.name === 'value').twoWay, true);
  assert.equal(fci.props.find((p) => p.name === 'functionName').bindable, true);
  const ff = reg.get('u2-func-form');
  assert.notEqual(ff.props.find((p) => p.name === 'functionName').bindable, true, 're-render tier');
  assert.equal(ff.props.find((p) => p.name === 'callId').bindable, true);
  assert.equal(ff.props.find((p) => p.name === 'showHistory').bindable, true);
  assert.equal(ff.props.find((p) => p.name === 'showRun').bindable, true);

  registerAll();
  assert.ok(globalRegistry.metas().every((m) => !m.tag.startsWith('u2-viewer-')), 'the manifest stays platform-free');
  assert.deepEqual(platformSamples(globalRegistry), []);
  assert.deepEqual(platformSamples(undefined), []);
  assert.equal(platformSamples(reg), VIEWER_SAMPLES);
});

scoped('a spec renders a viewer over a bound table, links a look prop two-way and routes its events', async () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  const scope = new Scope();
  const df = frame(rows(3));
  const x = signal('total');
  let hits = 0;
  const ctx = new SpecContext({data: {orders: dfBindings(signal(df), scope), x}, commands: {hit: () => hits++}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [
    {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.orders'}, props: {allowEdit: true}},
    {tag: 'u2-viewer-scatter-plot', name: 'plot', bind: {table: '$.orders', xColumnName: '$.x'},
      props: {yColumnName: 'orderId', look: {markerMinSize: 3}}, on: {'d4-data-event': 'cmd:hit'}},
  ]}}, ctx, reg);
  assert.equal(brokenCount(instance), 0);
  const grid = instance.node('grid');
  assert.ok(grid instanceof Grid);
  assert.equal(grid.dataFrame, df);
  assert.equal(grid.props.allowEdit, true);
  assert.equal(grid.root.dataset.u2Name, 'grid');
  assert.equal(grid.componentMeta.tag, 'u2-viewer-grid');

  const plot = instance.node('plot');
  assert.deepEqual([plot.props.xColumnName, plot.props.yColumnName, plot.props.markerMinSize], ['total', 'orderId', 3]);
  x.value = 'customer';
  assert.equal(plot.props.xColumnName, 'customer');
  plot.props.xColumnName = 'orderId';
  await flush();
  assert.equal(x.value, 'orderId', 'two-way through link');

  plot.onEvent('d4-data-event').fire({});
  assert.equal(hits, 1, 'the viewer\'s own event bus');

  assert.deepEqual(instance.dump().root.children[1], {tag: 'u2-viewer-scatter-plot', name: 'plot',
    props: {yColumnName: 'orderId', look: {markerMinSize: 3}}, bind: {table: '$.orders', xColumnName: '$.x'},
    on: {'d4-data-event': 'cmd:hit'}}, 'bind and literal props only');

  instance.dispose();
  assert.ok(platform.kills.includes(grid.root) && platform.kills.includes(plot.root));
  assert.equal(grid.isDetached, true);
  scope.dispose();
});

scoped('the Design/Run toggle re-creates a viewer bound to a source; the old one is killed', async () => {
  fakeBackends();
  const reg = new Registry();
  registerPlatformComponents(reg);
  const instance = renderSpec({$schema: 'dg-ui/1',
    components: [{tag: 'u2-func-source', name: 'orders',
      props: {func: 'demoOrders', designData: 'sample', debounce: 0, sample: rows(3)}}],
    root: {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.orders'}}},
  new SpecContext(), reg, {designTime: true});
  const before = instance.node('grid');
  assert.equal(before.dataFrame.rowCount, 3, 'the sample rows, without a run');
  assert.equal(backends.funcRunner.calls.length, 0);

  instance.setDesignTime(false);
  await settle();
  const after = instance.node('grid');
  assert.notEqual(after, before);
  assert.equal(before.isDetached, true);
  assert.ok(platform.kills.includes(before.root));
  assert.equal(backends.funcRunner.calls.length, 1);
  assert.equal(after.dataFrame.rowCount, 4, 'the live run');
  assert.equal(brokenCount(instance), 0);
  instance.dispose();
});

scoped('the viewer sample renders zero-broken at design time over real sample rows, nothing run', () => {
  const runner = fakeBackends();
  const reg = new Registry();
  registerPlatformComponents(reg);
  assert.deepEqual(VIEWER_SAMPLES.map((s) => s.name), ['Master–detail (grid + filters)']);
  const instance = renderSpec(JSON.parse(JSON.stringify(VIEWER_SAMPLES[0].spec)), new SpecContext(), reg, {designTime: true});
  assert.equal(brokenCount(instance), 0);
  const grid = instance.node('grid');
  assert.ok(grid instanceof Grid);
  assert.equal(grid.dataFrame.rowCount, 3);
  assert.equal(grid.props.allowEdit, false);
  assert.deepEqual(instance.node('filters').props.columnNames, ['city', 'total']);
  assert.equal(instance.node('form').dataFrame, grid.dataFrame);
  assert.equal(instance.node('customerField').value.peek(), 'Aspirin Labs');
  assert.equal(runner.calls.length, 0, 'DD9: nothing ran');
  instance.dispose();
  assert.equal(platform.kills.length, 3);
});
