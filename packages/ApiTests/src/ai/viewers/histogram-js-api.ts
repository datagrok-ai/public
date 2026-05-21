import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectPropAndLook, expectRoundTripPropAndLook, findProp, look, subscribeAll, withAttachedViewer} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:219 (DG.Viewer.histogram),
// public/js-api/src/viewer.ts:671 (HistogramViewer),
// public/js-api/src/dataframe/data-frame.ts:565 (df.plot.histogram).
// Histogram JS-only surface: typed factory + df.plot.histogram, the
// includeDefaults flag on getOptions, setOptions round-trip, dataFrame
// getter/setter swap, getInfo/getProperties shape, the four event
// Observables (onBinsSelected / onLineSelected / onMouseOverBins /
// onMouseOverLine), and close() on attached. Note: close() on a
// never-attached viewer (factory but no addViewer) throws inside
// grok_Viewer_Close — that's why detached viewers are not close()-d.
category('AI: Viewers: Histogram JS API', () => {
  test('factory typed', async () => {
    const v = DG.Viewer.histogram(demog(100), {value: 'age', bins: 15});
    expect(v instanceof DG.HistogramViewer, true);
    expectPropAndLook(v, {bins: 15, valueColumnName: 'age'});
  });

  test('factory via DataFrame.plot.histogram', async () => {
    const df = demog();
    const v = df.plot.histogram({value: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.dataFrame === df, true);
    expectPropAndLook(v, {valueColumnName: 'height'});
  });

  test('getOptions excludes defaults by default', async () => {
    const v = DG.Viewer.histogram(demog());
    const lean = v.getOptions().look;
    const full = look(v);
    expect(Object.keys(full).length > Object.keys(lean).length, true);
    expect(Object.keys(full).length >= 10, true);
    expect('bins' in full, true);
  });

  test('setOptions round-trip', async () => {
    const v = DG.Viewer.histogram(demog(), {value: 'age'});
    expectRoundTripPropAndLook(v, {bins: 7, showXAxis: true, splitColumnName: 'race'});
    const opts = v.getOptions(true);
    expect(opts.type, DG.VIEWER.HISTOGRAM);
    expect(typeof opts.id, 'string');
  });

  test('dataFrame getter setter swap', async () => {
    const df1 = demog(20);
    const v = DG.Viewer.histogram(df1, {value: 'age'});
    expect(v.dataFrame === df1, true);
    const df2 = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.INT, 'x', [1, 2, 3, 4, 5, 6, 7, 8])]);
    v.dataFrame = df2;
    expect(v.dataFrame.rowCount, 8);
    expect(v.dataFrame.columns.contains('x'), true);
  });

  test('getInfo getProperties descriptor shape', async () => {
    const v = DG.Viewer.histogram(demog(20), {value: 'age'});
    expect(typeof v.getInfo(), 'object');
    expect(findProp(v, 'bins') != null, true);
    expect(v.descriptor.name, DG.VIEWER.HISTOGRAM);
  });

  test('event streams are rxjs Observables', async () => {
    const v = DG.Viewer.histogram(demog(20), {value: 'age'});
    subscribeAll([v.onBinsSelected, v.onLineSelected, v.onMouseOverBins, v.onMouseOverLine])();
  });

  test('close detaches without throwing on attached viewer', async () => {
    await withAttachedViewer<DG.HistogramViewer>(demog(20), DG.VIEWER.HISTOGRAM, {value: 'age'},
      (v) => v.close());
  });
});
