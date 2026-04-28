import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: Histogram JS API', () => {
  test('factory typed', async () => {
    const df = grok.data.demo.demog(100);
    const v = DG.Viewer.histogram(df, {value: 'age', bins: 15});
    expect(v instanceof DG.HistogramViewer, true);
    expect(v.type, DG.VIEWER.HISTOGRAM);
    expect(v.props['bins'], 15);
    expect(v.props['valueColumnName'], 'age');
    const look = v.getOptions(true).look;
    expect(look['bins'], 15);
    expect(look['valueColumnName'], 'age');
  });

  test('factory via DataFrame.plot.histogram', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.HISTOGRAM);
    expect(v.dataFrame === df, true);
    expect(v.props['valueColumnName'], 'height');
  });

  test('getOptions excludes defaults by default', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df);
    const lean = v.getOptions().look;
    const full = v.getOptions(true).look;
    expect(Object.keys(full).length > Object.keys(lean).length, true);
    expect(Object.keys(full).length >= 10, true);
    expect('bins' in full, true);
  });

  test('setOptions round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    v.setOptions({bins: 7, showXAxis: true, splitColumnName: 'race'});
    expect(v.props['bins'], 7);
    expect(v.props['showXAxis'], true);
    expect(v.props['splitColumnName'], 'race');
    const opts = v.getOptions(true);
    expect(opts.type, DG.VIEWER.HISTOGRAM);
    expect(typeof opts.id, 'string');
    expect(opts.look['bins'], 7);
    expect(opts.look['showXAxis'], true);
    expect(opts.look['splitColumnName'], 'race');
  });

  test('dataFrame getter setter swap', async () => {
    const df1 = grok.data.demo.demog(20);
    const v = DG.Viewer.histogram(df1, {value: 'age'});
    expect(v.dataFrame === df1, true);
    const df2 = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'x', [1, 2, 3, 4, 5, 6, 7, 8]),
    ]);
    v.dataFrame = df2;
    expect(v.dataFrame.rowCount, 8);
    expect(v.dataFrame.columns.contains('x'), true);
  });

  test('getInfo getProperties descriptor shape', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const info = v.getInfo();
    expect(info != null, true);
    expect(typeof info, 'object');
    const props = v.getProperties();
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    var hasBins = false;
    for (var p of props)
      if (p.name === 'bins') { hasBins = true; break; }
    expect(hasBins, true);
    expect(v.descriptor.name, DG.VIEWER.HISTOGRAM);
  });

  test('event streams are rxjs Observables', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const subs: Subscription[] = [];
    try {
      const streams = [v.onBinsSelected, v.onLineSelected, v.onMouseOverBins, v.onMouseOverLine];
      for (var s of streams) {
        expect(s != null, true);
        expect(typeof s.subscribe, 'function');
        const sub = s.subscribe(() => {});
        expect(typeof sub.unsubscribe, 'function');
        subs.push(sub);
      }
    }
    finally {
      for (var sub of subs)
        sub.unsubscribe();
    }
  });

  test('close detaches without throwing on attached viewer', async () => {
    const df = grok.data.demo.demog(20);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'}) as DG.HistogramViewer;
      expect(v.type, DG.VIEWER.HISTOGRAM);
      var threw = false;
      try {
        v.close();
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
    }
    finally {
      tv.close();
    }
  });
});
