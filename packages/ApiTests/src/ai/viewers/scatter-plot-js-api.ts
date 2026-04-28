import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: ScatterPlot JS API', () => {
  test('df.plot.scatter shorthand returns typed ScatterPlotViewer', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({xColumnName: 'age', yColumnName: 'height', colorColumnName: 'race'});
    expect(v instanceof DG.ScatterPlotViewer, true);
    expect(v.type, DG.VIEWER.SCATTER_PLOT);
    expect(v.dataFrame === df, true);
    expect(v.props['xColumnName'], 'age');
    expect(v.props['yColumnName'], 'height');
    expect(v.props['colorColumnName'], 'race');
    const look = v.getOptions(true).look;
    expect(look['xColumnName'], 'age');
    expect(look['yColumnName'], 'height');
    expect(look['colorColumnName'], 'race');
  });

  test('meta.formulaLines helper round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    const fl = v.meta.formulaLines;
    expect(fl != null, true);
    expect(Array.isArray(fl.items), true);
    expect(fl.items.length, 0);
    fl.add({title: 'AI test line', formula: '${height} = ${age} + 100', type: 'line'});
    const items = fl.items;
    expect(items.length, 1);
    expect(items[0].title, 'AI test line');
    expect(items[0].formula, '${height} = ${age} + 100');
    const raw = v.props['formulaLines'];
    expect(typeof raw, 'string');
    const parsed = JSON.parse(raw);
    expect(Array.isArray(parsed), true);
    expect(parsed.length, 1);
    expect(parsed[0].title, 'AI test line');
    fl.clear();
    expect(fl.items.length, 0);
  });

  test('meta.annotationRegions helper round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    const ar = v.meta.annotationRegions;
    expect(ar != null, true);
    expect(Array.isArray(ar.items), true);
    expect(ar.items.length, 0);
    ar.add({header: 'AI test region', description: 'helper round-trip', opacity: 50});
    const items = ar.items;
    expect(items.length, 1);
    expect(items[0].header, 'AI test region');
    expect(items[0].description, 'helper round-trip');
    const raw = v.props['annotationRegions'];
    expect(typeof raw, 'string');
    const parsed = JSON.parse(raw);
    expect(Array.isArray(parsed), true);
    expect(parsed.length, 1);
    expect(parsed[0].header, 'AI test region');
    ar.clear();
    expect(ar.items.length, 0);
  });

  test('zoom() reflects xMin/xMax/yMin/yMax in attached look', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.SCATTER_PLOT,
        {xColumnName: 'age', yColumnName: 'height'}) as DG.ScatterPlotViewer;
      expect(v instanceof DG.ScatterPlotViewer, true);
      v.zoom(20, 150, 60, 200);
      const look = v.getOptions(true).look;
      expect('xMin' in look, true);
      expect('xMax' in look, true);
      expect('yMin' in look, true);
      expect('yMax' in look, true);
    }
    finally {
      tv.close();
    }
  });

  test('Scatter event Observables are subscribable', async () => {
    const df = grok.data.demo.demog(20);
    const v = df.plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    const subs: Subscription[] = [];
    try {
      const streams: Observable<any>[] = [
        v.onZoomed,
        v.onResetView,
        v.onViewportChanged,
        v.onAfterDrawScene,
        v.onBeforeDrawScene,
        v.onPointClicked,
        v.onPointDoubleClicked,
      ];
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

  test('view.addViewer typed instance + viewBox/getInfo geometry', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.SCATTER_PLOT,
        {xColumnName: 'age', yColumnName: 'height'}) as DG.ScatterPlotViewer;
      expect(v.type, DG.VIEWER.SCATTER_PLOT);
      expect(v instanceof DG.ScatterPlotViewer, true);
      var found: DG.Viewer | undefined;
      for (var x of tv.viewers)
        if (x.type === DG.VIEWER.SCATTER_PLOT) { found = x; break; }
      expect(found != null, true);
      expect(found instanceof DG.ScatterPlotViewer, true);
      const vb = v.viewBox;
      expect(vb instanceof DG.Rect, true);
      const info = v.getInfo();
      expect(info != null, true);
      expect('canvas' in info, true);
      expect('overlay' in info, true);
    }
    finally {
      tv.close();
    }
  });
});
