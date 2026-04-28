import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: BoxPlot JS API', () => {
  test('factory DG.Viewer.boxPlot returns typed DG.BoxPlot; df.plot.box returns base Viewer', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'race'});
    expect(v instanceof DG.BoxPlot, true);
    expect(v.type, DG.VIEWER.BOX_PLOT);
    expect(v.dataFrame === df, true);
    expect(v.props['valueColumnName'], 'age');
    expect(v.props['category1ColumnName'], 'race');
    const v2 = df.plot.box({valueColumnName: 'age', category1ColumnName: 'race'});
    expect(v2 instanceof DG.Viewer, true);
    expect(v2.type, DG.VIEWER.BOX_PLOT);
  });

  test('friendly-key aliases value/category1/category2 collapse to canonical names', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race', category2: 'sex'});
    expect(v instanceof DG.BoxPlot, true);
    expect(v.props['valueColumnName'], 'age');
    expect(v.props['category1ColumnName'], 'race');
    expect(v.props['category2ColumnName'], 'sex');
    const look = v.getOptions(true).look;
    expect(look['valueColumnName'], 'age');
    expect(look['category1ColumnName'], 'race');
    expect(look['category2ColumnName'], 'sex');
  });

  test('combined category1ColumnName + category2ColumnName JSON envelope round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'race'});
    v.setOptions({category1ColumnName: 'race', category2ColumnName: 'sex'});
    expect(v.props['category1ColumnName'], 'race');
    expect(v.props['category2ColumnName'], 'sex');
    const look = v.getOptions(true).look;
    expect(look['category1ColumnName'], 'race');
    expect(look['category2ColumnName'], 'sex');
  });

  test('showStatistics + statisticsFormat combined round-trip with getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'race'});
    v.setOptions({showStatistics: true, statisticsFormat: 'two digits after comma'});
    expect(v.props['showStatistics'], true);
    expect(v.props['statisticsFormat'], 'two digits after comma');
    const look = v.getOptions(true).look;
    expect(look['showStatistics'], true);
    expect(look['statisticsFormat'], 'two digits after comma');
    const props = v.getProperties();
    var fmtProp: DG.Property | undefined;
    for (var p of props)
      if (p.name === 'statisticsFormat') fmtProp = p;
    expect(fmtProp != null, true);
    expect(Array.isArray(fmtProp!.choices), true);
    expect(fmtProp!.choices.length > 0, true);
    expect(fmtProp!.choices.indexOf('two digits after comma') >= 0, true);
  });

  test('onResetView, onAfterDrawScene, onBeforeDrawScene, onPointClicked are rxjs Observables', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'sex'});
    expect(v instanceof DG.BoxPlot, true);
    const subs: Subscription[] = [];
    try {
      const streams: Observable<any>[] = [v.onResetView, v.onAfterDrawScene, v.onBeforeDrawScene, v.onPointClicked];
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

  test('view.addViewer attaches a typed DG.BoxPlot', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BOX_PLOT,
        {valueColumnName: 'age', category1ColumnName: 'race'}) as DG.BoxPlot;
      expect(v.type, DG.VIEWER.BOX_PLOT);
      expect(v instanceof DG.BoxPlot, true);
      var found: DG.Viewer | undefined;
      for (var x of tv.viewers)
        if (x.type === DG.VIEWER.BOX_PLOT) { found = x; break; }
      expect(found != null, true);
      expect(found instanceof DG.BoxPlot, true);
    }
    finally {
      tv.close();
    }
  });
});
