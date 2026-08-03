import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectPropAndLook, look, subscribeAll, withAttachedViewer} from '../helpers';

// ScatterPlot JS API: factory, meta helpers, zoom, events, addViewer geometry.
category('AI: Viewers: ScatterPlot JS API', () => {
  test('df.plot.scatter shorthand returns typed ScatterPlotViewer', async () => {
    const df = demog();
    const v = df.plot.scatter({xColumnName: 'age', yColumnName: 'height', colorColumnName: 'race'});
    expect(v instanceof DG.ScatterPlotViewer, true);
    expect(v.dataFrame === df, true);
    expectPropAndLook(v, {xColumnName: 'age', yColumnName: 'height', colorColumnName: 'race'});
  });

  test('meta.formulaLines helper round-trip', async () => {
    const v = demog().plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    const fl = v.meta.formulaLines;
    expect(fl.items.length, 0);
    fl.add({title: 'AI test line', formula: '${height} = ${age} + 100', type: 'line'});
    expect(fl.items.length, 1);
    expect(fl.items[0].title, 'AI test line');
    const parsed = JSON.parse(v.props['formulaLines']);
    expect(parsed.length, 1);
    expect(parsed[0].title, 'AI test line');
    fl.clear();
    expect(fl.items.length, 0);
  });

  test('meta.annotationRegions helper round-trip', async () => {
    const v = demog().plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    const ar = v.meta.annotationRegions;
    expect(ar.items.length, 0);
    ar.add({header: 'AI test region', description: 'helper round-trip', opacity: 50});
    expect(ar.items.length, 1);
    expect(ar.items[0].header, 'AI test region');
    const parsed = JSON.parse(v.props['annotationRegions']);
    expect(parsed.length, 1);
    expect(parsed[0].header, 'AI test region');
    ar.clear();
    expect(ar.items.length, 0);
  });

  test('zoom() reflects xMin/xMax/yMin/yMax in attached look', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {xColumnName: 'age', yColumnName: 'height'}, (v) => {
        v.zoom(20, 150, 60, 200);
        const l = look(v);
        for (const k of ['xMin', 'xMax', 'yMin', 'yMax'])
          expect(k in l, true);
      });
  });

  test('Scatter event Observables are subscribable', async () => {
    const v = demog(20).plot.scatter({xColumnName: 'age', yColumnName: 'height'});
    subscribeAll([v.onZoomed, v.onResetView, v.onViewportChanged,
      v.onAfterDrawScene, v.onBeforeDrawScene, v.onPointClicked, v.onPointDoubleClicked])();
  });

  test('view.addViewer typed instance + viewBox/getInfo geometry', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {xColumnName: 'age', yColumnName: 'height'}, (v) => {
        expect(v instanceof DG.ScatterPlotViewer, true);
        expect(v.viewBox instanceof DG.Rect, true);
        const info = v.getInfo();
        expect('canvas' in info, true);
        expect('overlay' in info, true);
      });
  });
});
