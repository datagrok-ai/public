import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: LineChart JS API', () => {
  test('factory friendly-key aliasing via DG.Viewer.lineChart', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height'], split: 'race'});
    expect(v instanceof DG.LineChartViewer, true);
    expect(v.type, DG.VIEWER.LINE_CHART);
    expect(v.props['xColumnName'], 'age');
    const ys = v.props['yColumnNames'];
    expect(Array.isArray(ys), true);
    expect(ys.indexOf('height') >= 0, true);
    expect(v.props['splitColumnName'], 'race');
    const look = v.getOptions(true).look;
    expect(look['xColumnName'], 'age');
    expect(Array.isArray(look['yColumnNames']), true);
    expect((look['yColumnNames'] as string[]).indexOf('height') >= 0, true);
  });

  test('factory via df.plot.line shorthand', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.line({xColumnName: 'age', yColumnNames: ['weight']});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.LINE_CHART);
    expect(v.dataFrame === df, true);
    expect(v.props['xColumnName'], 'age');
    const ys = v.props['yColumnNames'];
    expect(Array.isArray(ys), true);
    expect(ys.indexOf('weight') >= 0, true);
  });

  test('yColumnNames and splitColumnNames array round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {
      xColumnName: 'age',
      yColumnNames: ['height', 'weight'],
      splitColumnNames: ['race', 'sex'],
    });
    const ys = v.props['yColumnNames'];
    expect(Array.isArray(ys), true);
    expectArray(ys, ['height', 'weight']);
    const splits = v.props['splitColumnNames'];
    expect(Array.isArray(splits), true);
    expectArray(splits, ['race', 'sex']);
    const look = v.getOptions(true).look;
    expect(Array.isArray(look['yColumnNames']), true);
    expectArray(look['yColumnNames'] as string[], ['height', 'weight']);
    expect(Array.isArray(look['splitColumnNames']), true);
    expectArray(look['splitColumnNames'] as string[], ['race', 'sex']);
  });

  test('multiAxis and splineTension round-trip via setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {xColumnName: 'age', yColumnNames: ['height', 'weight']});
    expect(v.props['multiAxis'], false);
    v.setOptions({multiAxis: true, splineTension: 0.5});
    expect(v.props['multiAxis'], true);
    expect(v.props['splineTension'], 0.5);
    const look = v.getOptions(true).look;
    expect(look['multiAxis'], true);
    expect(look['splineTension'], 0.5);
    const props = v.getProperties();
    var multiAxisProp: DG.Property | undefined;
    var splineProp: DG.Property | undefined;
    for (var p of props) {
      if (p.name === 'multiAxis') multiAxisProp = p;
      if (p.name === 'splineTension') splineProp = p;
    }
    expect(multiAxisProp != null, true);
    expect(splineProp != null, true);
  });

  test('onLineSelected and onZoomed are rxjs Observables', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.lineChart(df, {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
    expect(v instanceof DG.LineChartViewer, true);
    const subs: Subscription[] = [];
    try {
      const streams: Observable<any>[] = [v.onLineSelected, v.onZoomed, v.onResetView];
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

  test('view.addViewer attaches a typed LineChartViewer with activeFrame', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.LINE_CHART, {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      expect(v.type, DG.VIEWER.LINE_CHART);
      expect(v instanceof DG.LineChartViewer, true);
      var found: DG.Viewer | undefined;
      for (var x of tv.viewers)
        if (x.type === DG.VIEWER.LINE_CHART) { found = x; break; }
      expect(found != null, true);
      expect(found instanceof DG.LineChartViewer, true);
      const af = v.activeFrame;
      expect(af != null, true);
      const afJs = DG.toJs(af);
      expect(afJs instanceof DG.DataFrame, true);
    }
    finally {
      tv.close();
    }
  });
});
