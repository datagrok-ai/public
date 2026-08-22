/* The contract of tests/platform-doubles.mjs: every double keeps its fields on prototype getters
   over the dart handle, as the real entity does — so a spread copies nothing but `dart` — and the
   loader hooks serve these very classes, never a second copy. The viewer doubles also pin the
   probed platform facts the viewer integration designs around (viewers/plan.md §1: P6, P7, P9,
   P11, P14) and the kill-walk globals dg-stub installs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import './dom-shim.js';
import {BitSet, Column, DartWidget, DataFrame, DataQuery, Entity, EventType, FileInfo, FilterGroup, Func,
  Grid, JsViewer, Package, Property, Script, Shell, TableQuery, UnreadableFileInfo, User, Viewer,
  ViewerMetaHelper, Widget, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const DG = await import('datagrok-api/dg');
const grok = await import('datagrok-api/grok');

const demog = () => new DataFrame([{name: 'age', type: 'int', semType: 'Age'}, {name: 'sex'}],
  [{age: 40, sex: 'F'}, {age: 55, sex: 'M'}], 'demog');

/** The WO-V3 fixture: a descriptor per type the `Viewer` statics build, declared as the live list is. */
const descriptors = () => [
  new WidgetDescriptor('Grid', [
    new Property('allowEdit', 'bool', {userEditable: true, category: 'Data'}),
    new Property('columnNames', 'list', {subType: 'string'}),
    new Property('filters', 'list', {subType: 'map', userEditable: null}),
  ], {events: [{name: 'd4-grid-cell-click', eventName: 'OnCellClicked'}]}),
  new WidgetDescriptor('Filters', [new Property('filters', 'list', {subType: 'map'}),
    new Property('columnNames', 'list', {subType: 'string'})]),
  new WidgetDescriptor('Form', []),
  new WidgetDescriptor('Scatter plot', [
    new Property('xColumnName', 'string'), new Property('yColumnName', 'string'),
    new Property('markerMinSize', 'num', {min: 1, max: 50}),
    new Property('showRegressionLine', 'bool', {defaultValue: false}),
  ], {description: 'Dots', synonyms: ['xy plot']}),
];

const doubles = () => {
  const frame = demog();
  const [descriptor] = descriptors();
  const grid = new Grid({type: 'Grid', descriptor, dataFrame: frame});
  return [
    new Property('days', 'int', {description: 'How many'}),
    new Property('size', 'num', {min: 1, max: 50, step: 0.5, inputType: 'Float', editor: 'slider', format: '0.0',
      units: 'px', showSlider: true, showPlusMinus: false}),
    new Func('orders'), new DataQuery('sql'), new Script('py'), new TableQuery('visual'),
    new FileInfo('Demo:Files/geo/cities.csv'), new UnreadableFileInfo('keys/server.pem'),
    frame, frame.columns, frame.columns.byName('age'), frame.selection, new BitSet(3),
    new Entity('alice'), new User('alice'), new Package('Chem'), new Shell(),
    descriptor, descriptor.events[0], grid, grid.meta, new FilterGroup({type: 'Filters', descriptor}),
    new Viewer({type: 'Form', descriptor}), new JsViewer(), new Widget(document.createElement('div')),
    new DartWidget({type: 'Legend'}),
  ];
};

// stricter than the real DataFrame, which also keeps columns/rows/filter/temp/tags as own fields — none read by u2
test('spreading any double copies the dart handle and nothing else', () => {
  for (const double of doubles()) {
    assert.deepEqual(Object.keys({...double}), ['dart'], double.constructor.name);
    assert.deepEqual(Object.keys(double), ['dart'], double.constructor.name);
  }
});

test('a property answers what it was created with, and what fromOptions refined it with', () => {
  const get = () => 1;
  const p = Property.create('days', 'int', get, null, 7).fromOptions({description: 'How many',
    choices: ['1', '7'], caption: 'Days', nullable: false, category: 'Misc', semType: 'Span'});
  assert.deepEqual([p.name, p.propertyType, p.type, p.get, p.set, p.defaultValue],
    ['days', 'int', 'int', get, null, 7]);
  assert.deepEqual([p.description, p.choices, p.caption, p.nullable, p.category, p.semType],
    ['How many', ['1', '7'], 'Days', false, 'Misc', 'Span']);
  assert.equal(new Property('x', 'string').choices, null, 'unset choices are null, as the Dart field is');
  assert.equal(new Property('x', 'string').caption, undefined, 'what was not given is not there');
  const size = new Property('size', 'num', {min: 1, max: 50, step: 0.5, inputType: 'Float', editor: 'slider',
    format: '0.0', units: 'px', showSlider: true, showPlusMinus: false});
  assert.deepEqual([size.min, size.max, size.step, size.inputType, size.editor, size.format, size.units,
    size.showSlider, size.showPlusMinus], [1, 50, 0.5, 'Float', 'slider', '0.0', 'px', true, false],
  'the numeric-editor fields read through getters too (property.ts:179-237)');
  assert.deepEqual(Object.keys({...size}), ['dart']);
});

test('a function answers its declaration, its kind and the call it prepares', async () => {
  const chem = new Package('Chem');
  const inputs = [new Property('smiles', 'string')];
  const f = new Func('save', {package: chem, friendlyName: 'Save it', description: 'Saves',
    inputs, options: {role: 'file-handler'}, run: async (params) => ({ok: params.smiles})});
  assert.deepEqual([f.name, f.friendlyName, f.description, f.package, f.inputs, f.outputs, f.options],
    ['save', 'Save it', 'Saves', chem, inputs, [], {role: 'file-handler'}]);
  assert.deepEqual(await f.prepare({smiles: 'CCO'}).call(), {outputs: {ok: 'CCO'}});
  assert.equal(new Func('orders', {namespace: 'admin'}).nqName, 'admin:orders');
  assert.equal(new Func('orders').nqName, 'orders');

  const query = new TableQuery('visual');
  assert.ok(query instanceof DataQuery && query instanceof Func && query instanceof Entity);
  query.limit = 100;
  assert.equal(query.limit, 100);
  assert.deepEqual(Object.keys(query), ['dart'], 'the limit went behind the handle');

  Func.registry = [f, new Func('save'), new Func('other', {options: {role: 'file-handler'}})];
  assert.equal(Func.find({name: 'save'}).length, 2);
  assert.deepEqual(Func.find({name: 'save', package: 'Chem'}), [f]);
  assert.equal(Func.find({meta: {role: 'file-handler'}}).length, 2);
  Func.registry = [];
});

test('a file reads its names off the path, and its share off the connection', () => {
  const home = new Entity('Files', {nqName: 'Demo:Files'});
  const file = new FileInfo('Demo:Files/geo/cities.csv', {connection: home});
  assert.deepEqual([file.fullPath, file.path, file.fileName, file.extension, file.name],
    ['Demo:Files/geo/cities.csv', 'geo/cities.csv', 'cities.csv', 'csv', 'cities.csv']);
  assert.deepEqual([file.connection, file.isFile, file.isDirectory, file.data], [home, true, false, null]);

  const folder = new FileInfo('a/tables', {directory: true});
  assert.deepEqual([folder.path, folder.extension, folder.isFile, folder.isDirectory], ['a/tables', '', false, true]);
  assert.equal(new FileInfo('Demo:Files').path, '', 'the share itself has no path under it');
  const local = FileInfo.fromBytes('server.pem', new Uint8Array([1]));
  assert.deepEqual([local.connection, local.data.length], [null, 1]);
  assert.throws(() => new UnreadableFileInfo('keys/server.pem').connection, RangeError);
});

test('a frame answers its cells, counts its reads and writes, and fires the platform events', () => {
  const df = demog();
  assert.deepEqual([df.name, df.rowCount, df.currentRowIdx, df.columns.names()], ['demog', 2, 0, ['age', 'sex']]);
  assert.ok(df.selection instanceof BitSet && df.filter instanceof BitSet);
  assert.equal(df.selection.length, 2);
  const age = df.columns.byName('age');
  assert.ok(age instanceof Column);
  assert.deepEqual([age.name, age.type, age.semType, age.isNumerical], ['age', 'int', 'Age', true]);
  assert.deepEqual([df.columns.byName('sex').type, df.columns.byName('sex').semType], ['string', null]);
  assert.equal(df.columns.byName('gone'), null);
  assert.equal(df.columns.contains('sex'), true);
  assert.equal(new Column('when', 'datetime').isNumerical, true, 'as in ddt');

  const seen = [];
  const watched = ['onCurrentRowChanged', 'onValuesChanged', 'onColumnsChanged', 'onColumnNameChanged'];
  const subs = watched.map((event) => df[event].subscribe((x) => seen.push([event, x])));
  assert.equal(df.liveSubscriptions(), 4);

  assert.equal(df.get('age', 1), 55);
  assert.equal(df.dart.reads, 1);
  df.set('age', 1, 56);
  assert.deepEqual(df.dart.writes, [['age', 1, 56]]);
  assert.equal(df.dart.rows[1].age, 56);
  df.currentRowIdx = 0;
  df.currentRowIdx = 1;
  assert.equal(df.dart.idxWrites, 2, 'every write counts');
  age.name = 'years';
  df.columns.remove('sex');
  df.columns.add(new Column('flag', 'bool'));
  assert.deepEqual(seen, [['onValuesChanged', {column: 'age', row: 1}], ['onCurrentRowChanged', 1],
    ['onColumnNameChanged', {args: {oldName: 'age', newName: 'years'}}],
    ['onColumnsChanged', undefined], ['onColumnsChanged', undefined]],
  'a rename fires onColumnNameChanged alone; an unchanged row index fires nothing');
  assert.deepEqual(df.columns.names(), ['years', 'flag']);
  df.columns.remove('gone');
  assert.deepEqual(df.columns.names(), ['years', 'flag'], 'an unknown name removes nothing');
  subs[0].unsubscribe();
  assert.equal(df.liveSubscriptions(), 3, 'an unsubscribed stream is let go');
  const bytes = new Uint8Array([1, 2]);
  const read = DataFrame.fromByteArray(bytes);
  assert.deepEqual([read.rowCount, read.dart.bytes], [0, bytes], 'the bytes stay behind the handle');
});

test('an entity row answers its names; a user adds an email', () => {
  const alice = new User('alice', {id: 'u1', friendlyName: 'Alice M', email: 'alice@datagrok.ai'});
  assert.deepEqual([alice.id, alice.name, alice.friendlyName, alice.nqName, alice.email],
    ['u1', 'alice', 'Alice M', 'alice', 'alice@datagrok.ai']);
  assert.equal(new Entity('bob').friendlyName, undefined,
    'no fallback here: the platform derives camelCaseToWords(name), and every u2 read goes friendlyName ?? name');
});

test('the shell records every assignment of the current object, and names tables apart', () => {
  const shell = new Shell();
  const ref = {};
  shell.o = ref;
  shell.o = null;
  assert.deepEqual(shell.dart.writes, [ref, null]);
  assert.equal(shell.o, null);
  assert.equal(shell.windows.showContextPanel, false);

  const added = [];
  shell.dart.tableAdded.subscribe(() => added.push(shell.tableNames));
  assert.equal(shell.addTable(demog()).name, 'demog');
  assert.equal(shell.addTable(demog()).name, 'demog (2)');
  assert.deepEqual(added, [['demog'], ['demog', 'demog (2)']]);
  shell.closeTable('demog');
  assert.deepEqual(shell.tableNames, ['demog (2)']);
  assert.doesNotThrow(() => shell.info('x') ?? shell.warning('x') ?? shell.error('x'));
});

test('a descriptor answers its declaration; its properties are defined over the look (P2, P4)', () => {
  WidgetDescriptor.registry = descriptors();
  const grid = WidgetDescriptor.getByName('Grid');
  assert.ok(grid instanceof WidgetDescriptor);
  assert.equal(WidgetDescriptor.getDescriptors().length, 4);
  assert.equal(WidgetDescriptor.getByName('Nope'), null);
  assert.deepEqual([grid.name, grid.description, grid.synonyms], ['Grid', '', []]);
  assert.deepEqual(WidgetDescriptor.getByName('Scatter plot').synonyms, ['xy plot']);
  assert.deepEqual(grid.properties.map((p) => [p.name, p.propertyType, p.userEditable, p.propertySubType, p.category]),
    [['allowEdit', 'bool', true, null, 'Data'], ['columnNames', 'list', undefined, 'string', undefined],
      ['filters', 'list', null, 'map', undefined]]);
  assert.deepEqual(grid.properties.map((p) => globalThis.grok_Property_Get_PropertySubType(p.dart)),
    [null, 'string', 'map'], 'the subtype is read off the handle, as the global does');
  assert.ok(grid.properties.every((p) => typeof p.get === 'function' && typeof p.set === 'function'),
    'get and set present on every descriptor property, as live (P1)');
  const [event] = grid.events;
  assert.ok(event instanceof EventType);
  assert.deepEqual([event.name, event.eventName], ['d4-grid-cell-click', 'OnCellClicked']);
  assert.deepEqual(WidgetDescriptor.getByName('Form').events, []);

  const given = new Property('x', 'int', {get: () => 1, set: null});
  new WidgetDescriptor('Own', [given]);
  assert.deepEqual([given.get(), given.set], [1, null], 'accessors given are kept');
  WidgetDescriptor.registry = [];
});

test('a viewer builds from its descriptor; its properties read the look, never the viewer (P5, P6)', () => {
  WidgetDescriptor.registry = descriptors();
  const df = demog();
  const grid = Viewer.fromType('Grid', df, {allowEdit: false});
  assert.ok(grid instanceof Grid && grid instanceof Viewer && grid instanceof Widget);
  assert.deepEqual([grid.type, grid.descriptor.name, grid.dataFrame, grid.table], ['Grid', 'Grid', df, df]);
  assert.equal(grid.root.getAttribute('data-widget'), 'true');
  assert.ok(platform.widgets.has(grid), 'a built viewer is registered with the platform');
  const look = globalThis.grok_Viewer_Get_Look(grid.dart);
  assert.deepEqual(look, {allowEdit: false, columnNames: null, filters: null},
    'every property on the look, defaults stamped, options applied');

  const reads = grid.dart.propertyReads;
  const [allowEdit] = grid.getProperties();
  assert.equal(grid.dart.propertyReads, reads + 1, 'every getProperties call counts on the handle');
  assert.ok(allowEdit instanceof Property);
  assert.throws(() => allowEdit.get(grid.dart), /NoSuchMethodError: method not found: 'allowEdit'/);
  assert.throws(() => allowEdit.get(grid), /Receiver: Instance of 'Grid'/);
  assert.throws(() => allowEdit.set(grid.dart, true), /NoSuchMethodError/);
  assert.equal(allowEdit.get(look), false);
  assert.equal(grid.props.allowEdit, false, 'the bag is over the look');
  assert.throws(() => grid.props.get('allowEdit'), /NoSuchMethodError/,
    'the bag\'s own get/set take the viewer, as the platform\'s do');
  assert.throws(() => grid.props.nope, /Property not found: nope/);
  assert.ok('allowEdit' in grid.props && !('nope' in grid.props));
  assert.deepEqual(Object.keys(grid.props), ['allowEdit', 'columnNames', 'filters']);

  assert.ok(Viewer.grid(df) instanceof Grid && Grid.create(df) instanceof Grid);
  assert.ok(Viewer.filters(df) instanceof FilterGroup);
  assert.ok(Viewer.scatterPlot(df) instanceof Viewer && !(Viewer.scatterPlot(df) instanceof Grid));
  assert.deepEqual([Viewer.form(df).type, Viewer.form(df).getProperties()], ['Form', []]);
  assert.throws(() => Viewer.fromType('Nope', df), /Unknown viewer type: Nope/);
  assert.throws(() => { grid.root = null; }, TypeError, 'a Dart viewer root is read-only, as in js-api');
  WidgetDescriptor.registry = [];
});

test('a property write fires onPropertyValueChanged with the property, a same value too; setOptions fires per key (P7)', () => {
  WidgetDescriptor.registry = descriptors();
  const plot = Viewer.scatterPlot(demog(), {xColumnName: 'age'});
  const [xColumnName] = plot.descriptor.properties;
  const seen = [];
  plot.onPropertyValueChanged.subscribe((e) => seen.push([e.args.property.name, e.args.property === xColumnName, e.data]));
  plot.props.xColumnName = 'sex';
  plot.props.xColumnName = 'sex';
  plot.setOptions({yColumnName: 'age', markerMinSize: 3, type: 'Scatter plot', nope: 1});
  assert.deepEqual(seen, [['xColumnName', true, undefined], ['xColumnName', true, undefined],
    ['yColumnName', false, undefined], ['markerMinSize', false, undefined]]);
  assert.deepEqual(plot.getOptions(), {type: 'Scatter plot', look: {xColumnName: 'sex', yColumnName: 'age', markerMinSize: 3}});
  assert.deepEqual(plot.getOptions(true).look,
    {xColumnName: 'sex', yColumnName: 'age', markerMinSize: 3, showRegressionLine: false});
  assert.ok(!('nope' in plot.props) && !('nope' in globalThis.grok_Viewer_Get_Look(plot.dart)),
    'an unknown option is ignored, not stored');
  assert.equal(plot.onEvent('d4-data-event'), plot.onDataEvent, 'one stream per id');
  assert.deepEqual([plot.onPropertyValueChanged.count, plot.liveSubscriptions()], [1, 1]);
  WidgetDescriptor.registry = [];
});

test('repointing the frame fires onDataFrameChanged (P11); meta is getter-only (P8); an empty frame is creatable', () => {
  WidgetDescriptor.registry = descriptors();
  const grid = Viewer.grid(DataFrame.create(0));
  assert.deepEqual([grid.dataFrame.rowCount, grid.dataFrame.columns.names()], [0, []]);
  let fired = 0;
  grid.onDataFrameChanged.subscribe(() => fired++);
  const df = demog();
  grid.dataFrame = df;
  assert.deepEqual([fired, grid.dataFrame, grid.dataFrame.dart === df.dart], [1, df, true]);
  assert.ok(grid.meta instanceof ViewerMetaHelper);
  assert.equal(grid.meta, grid.meta);
  assert.throws(() => { grid.meta = {}; }, TypeError);
  assert.equal(DataFrame.create(3).rowCount, 3);
  assert.equal(DataFrame.create().rowCount, 0);
  WidgetDescriptor.registry = [];
});

test('a filter group rebuilds through updateOrAdd, which writes the filters property (P14)', () => {
  WidgetDescriptor.registry = descriptors();
  const fg = Viewer.filters(demog(), {filters: [{type: 'categorical', column: 'sex'}]});
  assert.ok(fg instanceof FilterGroup);
  const seen = [];
  fg.onPropertyValueChanged.subscribe((e) => seen.push(e.args.property.name));
  fg.updateOrAdd({type: 'histogram', column: 'age'});
  assert.deepEqual(fg.props.filters, [{type: 'categorical', column: 'sex'}, {type: 'histogram', column: 'age'}]);
  fg.updateOrAdd({type: 'histogram', column: 'age', min: 1});
  assert.deepEqual(fg.props.filters, [{type: 'categorical', column: 'sex'}, {type: 'histogram', column: 'age', min: 1}],
    'the same column and type is updated in place');
  fg.props.filters = [];
  assert.deepEqual([seen, fg.props.filters], [['filters', 'filters', 'filters'], []],
    'a plain write fires and stores; only updateOrAdd would rebuild the live panel');
  WidgetDescriptor.registry = [];
});

test('a JS viewer owns its properties: addProperty declares over a field and its setter tells onPropertyChanged', () => {
  class Gauge extends JsViewer {
    constructor() {
      super();
      this.changed = [];
      this.level = this.addProperty('level', 'int', 5, {category: 'Data'});
      this.unit = this.addProperty('units', 'string', 'mg', {fieldName: 'unit'});
    }

    onPropertyChanged(p) { this.changed.push(p?.name ?? null); }
  }
  const g = new Gauge();
  assert.ok(g instanceof Viewer && g instanceof Widget && !(g instanceof Grid));
  assert.deepEqual([g.type, g.level, g.unit, g.descriptor, g.dataFrame], ['JsViewer', 5, 'mg', null, null]);
  assert.equal(g.getProperties(), g._properties);
  const [level, units] = g.getProperties();
  assert.deepEqual([level.name, level.propertyType, level.defaultValue, level.category, units.dart.fieldName],
    ['level', 'int', 5, 'Data', undefined]);
  assert.deepEqual([level.get(g), g.props.level, g.props.units], [5, 5, 'mg']);
  g.props.level = 7;
  g.props.units = 'g';
  assert.deepEqual([g.level, g.unit, g.changed], [7, 'g', ['level', 'units']]);
  level.set(g, 7);
  assert.deepEqual(g.changed, ['level', 'units', 'level'], 'a same value tells too');
  assert.equal(g.root.getAttribute('data-widget'), null);
  g.root = document.createElement('div');
  assert.ok(platform.widgets.has(g));

  let detached = 0;
  let cancelled = 0;
  g.onDetached.subscribe(() => detached++);
  g.sub({unsubscribe: () => cancelled++});
  g.detach();
  assert.deepEqual([g.isDetached, cancelled, detached], [true, 1, 0],
    'JS detach cancels the subs and does nothing Dart-side (P9)');
});

test('a widget registers through toDart; a Dart widget\'s properties take the handle', () => {
  const w = Widget.fromRoot(document.createElement('div'));
  assert.ok(!platform.widgets.has(w));
  assert.equal(w.toDart(), w.dart);
  assert.ok(platform.widgets.has(w), 'toDart registers');
  w.temp.x = 1;
  assert.deepEqual([w.type, w.dart.temp.x, w.getProperties(), w.getFunctions(), w.aiDescription], ['Unknown', 1, [], [], null]);
  w.aiDescription = 'briefed';
  assert.equal(w.dart.aiDescription, 'briefed');

  const x = Property.create('x', 'int', (d) => d.x, (d, v) => { d.x = v; }, 0);
  const legend = new DartWidget({type: 'Legend', properties: [x], x: 1});
  assert.deepEqual([legend.type, legend.getProperties(), legend.props.x], ['Legend', [x], 1]);
  legend.props.x = 2;
  assert.equal(legend.dart.x, 2);
  assert.equal(legend.root.getAttribute('data-widget'), null);
});

test('the kill-walk runs the cleanups under the element, then detaches every widget under it (P9)', () => {
  WidgetDescriptor.registry = descriptors();
  platform.reset();
  const parent = document.createElement('div');
  const pane = document.createElement('div');
  parent.appendChild(pane);
  const grid = Viewer.grid(demog());
  pane.appendChild(grid.root);
  const js = new JsViewer();
  parent.appendChild(js.root);
  const hosted = Widget.fromRoot(document.createElement('div'));
  parent.appendChild(hosted.root);
  hosted.toDart();
  const outside = Viewer.grid(demog());
  const log = [];
  globalThis.grok_Widget_RegisterCleanup(grid.root, () => log.push(['cleanup grid', grid.isDetached]));
  globalThis.grok_Widget_RegisterCleanup(pane, () => log.push(['cleanup pane', grid.isDetached]));
  globalThis.grok_Widget_RegisterCleanup(outside.root, () => log.push(['cleanup outside']));
  grid.onDetached.subscribe(() => log.push(['grid detached']));
  js.onDetached.subscribe(() => log.push(['js detached']));
  hosted.subs.push({unsubscribe: () => log.push(['hosted unsubscribed'])});
  assert.deepEqual([platform.cleanups.length, platform.widgets.size], [3, 4]);

  globalThis.grok_Widget_Kill(parent);
  assert.deepEqual(platform.kills, [parent]);
  assert.deepEqual(log, [['cleanup grid', false], ['cleanup pane', false], ['grid detached'], ['js detached'],
    ['hosted unsubscribed']], 'cleanups first, with the viewer still attached; then the detaches');
  assert.deepEqual([grid.isDetached, js.isDetached, hosted.isDetached, outside.isDetached], [true, true, true, false]);
  assert.equal(platform.cleanups.length, 1, 'the cleanups that ran are gone; the one outside stays');

  globalThis.grok_Widget_Kill(parent);
  assert.deepEqual([log.length, platform.kills], [5, [parent, parent]], 'a second kill is recorded and does nothing');
  platform.reset();
  assert.deepEqual([platform.kills, platform.cleanups, platform.widgets.size], [[], [], 0]);
  WidgetDescriptor.registry = [];
});

test('the loader hook serves the same classes, not a second copy', () => {
  assert.equal(DG.Property, Property);
  assert.equal(DG.Func, Func);
  assert.equal(DG.DataQuery, DataQuery);
  assert.equal(DG.Script, Script);
  assert.equal(DG.TableQuery, TableQuery);
  assert.equal(DG.FileInfo, FileInfo);
  assert.equal(DG.DataFrame, DataFrame);
  assert.equal(DG.Widget, Widget);
  assert.equal(DG.DartWidget, DartWidget);
  assert.equal(DG.Viewer, Viewer);
  assert.equal(DG.Grid, Grid);
  assert.equal(DG.FilterGroup, FilterGroup);
  assert.equal(DG.JsViewer, JsViewer);
  assert.equal(DG.WidgetDescriptor, WidgetDescriptor);
  assert.equal(DG.EventType, EventType);
  assert.equal(DG.platform, platform);
  assert.ok(['grok_Widget_Kill', 'grok_Widget_RegisterCleanup', 'grok_Property_Get_PropertySubType',
    'grok_Viewer_Get_Look'].every((name) => typeof globalThis[name] === 'function'), 'the globals are installed with DG');
  assert.ok(grok.shell instanceof Shell);
  assert.ok(new DG.FileInfo('a.csv') instanceof FileInfo, 'so an instanceof in the module under test holds');
});
