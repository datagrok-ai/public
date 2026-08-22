/* viewerOf / viewerControl / viewers.* / viewerSettings (VP-8, VP-18, VP-19) over the viewer
   doubles: the empty-frame start and the repoint by frame identity, named keys over `look`, the
   renderer's signals peeked and never linked, the fluent factories' two-way links, and the settings
   form over the look. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal, computed} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {DataFrame, FilterGroup, Grid, Property, Viewer, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {viewerOf, viewerControl} = await import('../src/dg/viewers/viewer-control.js');
const {viewers, viewerSettings} = await import('../src/dg/viewers/viewers.js');

const demog = () => new DataFrame([{name: 'age', type: 'int'}, {name: 'weight', type: 'double'}],
  [{age: 40, weight: 70}, {age: 55, weight: 80}], 'demog');

const descriptors = () => [
  new WidgetDescriptor('Grid', [new Property('allowEdit', 'bool', {defaultValue: false})]),
  new WidgetDescriptor('Filters', [new Property('columnNames', 'list', {subType: 'string'})]),
  new WidgetDescriptor('Form', []),
  new WidgetDescriptor('Scatter plot', [
    new Property('xColumnName', 'string'), new Property('yColumnName', 'string'),
    new Property('markerMinSize', 'num', {defaultValue: 1}),
  ]),
  new WidgetDescriptor('Summary', []),
];

function scoped(name, body) {
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

scoped('an unbound table starts the viewer over an empty frame; the first frame repoints it, a same frame does not', () => {
  const table = signal(undefined);
  const v = viewerControl('Grid', {table});
  assert.ok(v instanceof Grid && Control.is(v));
  assert.equal(v.dataFrame.rowCount, 0);
  let repoints = 0;
  v.onDataFrameChanged.subscribe(() => repoints++);

  const df = demog();
  table.value = df;
  assert.equal(v.dataFrame, df);
  assert.equal(repoints, 1);
  const twin = new DataFrame();
  twin.dart = df.dart;
  table.value = twin;
  assert.equal(repoints, 1, 'another wrapper over the same frame repoints nothing');
  table.value = undefined;
  assert.equal(v.dataFrame, df, 'an unbound frame leaves the viewer where it was');
  table.value = demog();
  assert.equal(repoints, 2);

  const plain = viewerControl('Grid', {});
  assert.equal(plain.dataFrame.rowCount, 0);
  plain.dispose();
  v.dispose();
  assert.equal(platform.kills.length, 2);
});

scoped('a type that cannot be built over the empty placeholder fails with what to do, the platform\'s error as the cause', () => {
  const fromType = Viewer.fromType;
  // Summary indexes the frame's columns at construction and throws on an empty one (M2)
  Viewer.fromType = (type, table, options) => {
    if (type === 'Summary' && table.rowCount === 0)
      throw new RangeError('Invalid value: Not in inclusive range 0..0: 0');
    return fromType(type, table, options);
  };
  try {
    assert.throws(() => viewerControl('Summary', {table: signal(undefined)}),
      (e) => e instanceof Error && e.message === 'Summary needs a table — bind `table` to a data source' &&
        e.cause instanceof RangeError);
    assert.throws(() => viewerControl('Summary', {}), {message: 'Summary needs a table — bind `table` to a data source'});
    const empty = new DataFrame([{name: 'age', type: 'int'}], [], 'empty');
    assert.throws(() => viewerControl('Summary', {table: signal(empty)}), RangeError,
      'a bound frame that fails is the platform\'s own error');
    const v = viewerControl('Summary', {table: signal(demog())});
    assert.equal(v.type, 'Summary');
    v.dispose();
  } finally {
    Viewer.fromType = fromType;
  }
});

scoped('a repoint re-applies the construction-time literal look, never a bound option', async () => {
  const table = signal(undefined);
  const fg = viewerControl('Filters', {table, columnNames: ['city', 'total']});
  let repoints = 0;
  // the platform resolves column-naming props against the frame and drops what the placeholder lacks
  fg.onDataFrameChanged.subscribe(() => {
    repoints++;
    fg.props.columnNames = [];
  });
  const df = demog();
  table.value = df;
  assert.deepEqual(fg.props.columnNames, ['city', 'total']);
  assert.equal(repoints, 1);
  const twin = new DataFrame();
  twin.dart = df.dart;
  table.value = twin;
  assert.equal(repoints, 1, 'the same frame repoints nothing');

  const x = signal('age');
  const v = viewerControl('Scatter plot', {table, xColumnName: x, look: {xColumnName: 'ignored', yColumnName: 'weight'}});
  v.props.xColumnName = 'weight';
  v.props.yColumnName = 'age';
  table.value = demog();
  assert.equal(v.props.xColumnName, 'weight', 'a renderer-bound prop is the link\'s, not re-applied');
  assert.equal(v.props.yColumnName, 'weight', 'a literal one is');

  const y = signal('weight');
  const w = viewerOf(Viewer.scatterPlot, table, {xColumnName: 'age', yColumnName: y});
  y.value = 'age';
  await flush();
  w.props.xColumnName = 'weight';
  table.value = demog();
  assert.deepEqual([w.props.xColumnName, w.props.yColumnName], ['age', 'age']);
  fg.dispose();
  v.dispose();
  w.dispose();
});

scoped('viewerControl: a named key wins over look, and a bound signal is peeked, never linked', () => {
  const x = signal('age');
  const v = viewerControl('Scatter plot', {table: signal(demog()), xColumnName: x,
    look: {yColumnName: 'weight', xColumnName: 'ignored'}});
  assert.equal(v.props.xColumnName, 'age');
  assert.equal(v.props.yColumnName, 'weight');
  x.value = 'weight';
  assert.equal(v.props.xColumnName, 'age', 'the renderer bridges bound props; the factory does not');
  v.props.xColumnName = 'age';
  assert.equal(x.value, 'weight');
  v.dispose();
});

scoped('viewerOf links a signal option two-way where writable, one-way where not', async () => {
  const x = signal('age');
  const y = signal('weight');
  const v = viewerOf(Viewer.scatterPlot, demog(), {xColumnName: x, yColumnName: computed(() => y.value), markerMinSize: 4});
  assert.deepEqual([v.props.xColumnName, v.props.yColumnName, v.props.markerMinSize], ['age', 'weight', 4]);
  x.value = 'weight';
  assert.equal(v.props.xColumnName, 'weight', 'follows');
  v.props.xColumnName = 'age';
  await flush();
  assert.equal(x.value, 'age', 'writes back');
  y.value = 'age';
  assert.equal(v.props.yColumnName, 'age', 'a read-only signal is followed');
  v.props.yColumnName = 'weight';
  await flush();
  assert.equal(y.value, 'age', 'and never written');
  v.dispose();
  x.value = 'weight';
  assert.equal(v.props.xColumnName, 'age', 'the links die with the viewer');
});

scoped('viewerOf takes a frame signal and the viewer is owned by the ambient scope', () => {
  const table = signal(undefined);
  const s = new Scope();
  const v = Scope.runWith(s, () => viewerOf(Viewer.grid, table, {allowEdit: true}));
  assert.equal(v.dataFrame.rowCount, 0);
  assert.equal(v.props.allowEdit, true);
  table.value = demog();
  assert.equal(v.dataFrame.rowCount, 2);
  s.dispose();
  assert.equal(v.scope.isDisposed, true);
  assert.deepEqual(platform.kills, [v.root]);
});

scoped('viewers.*: one-liners over the js-api statics', () => {
  const df = demog();
  const grid = viewers.grid(df, {allowEdit: true});
  assert.ok(grid instanceof Grid && Control.is(grid));
  assert.equal(grid.props.allowEdit, true);
  const filters = viewers.filters(df, {columnNames: ['age']});
  assert.ok(filters instanceof FilterGroup);
  assert.deepEqual(filters.props.columnNames, ['age']);
  const form = viewers.form(df);
  assert.equal(form.type, 'Form');
  const byType = viewers.byType('Scatter plot', df, {xColumnName: 'age'});
  assert.equal(byType.props.xColumnName, 'age');
  for (const v of [grid, filters, form, byType])
    v.dispose();
  assert.equal(platform.kills.length, 4);
});

scoped('viewerSettings builds a form whose fields read and write the look', () => {
  const v = viewers.scatterPlot(demog(), {xColumnName: 'age'});
  const form = viewerSettings(v);
  assert.equal(form.target, globalThis.grok_Viewer_Get_Look(v.dart));
  assert.equal(form.input('xColumnName').value.peek(), 'age');
  form.input('xColumnName').value.value = 'weight';
  assert.equal(v.props.xColumnName, 'weight');
  form.dispose();
  v.dispose();
});
