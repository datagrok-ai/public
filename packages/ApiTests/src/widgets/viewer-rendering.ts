import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Viewer: rendering', () => {
  // the BDD harness arms every viewer the page holds or adds, including the ones a layout
  // re-creates, and only after they are attached: the debounces armed during attach must not be
  // overtaken by the ones armed after the switch
  test('immediateRendering: a layout restores the filter panel', async () => {
    const df = grok.data.demo.demog(1000);
    df.name = 'demog immediate layout';
    const tv = grok.shell.addTableView(df);
    const arm = (v: DG.Viewer) => v.immediateRendering = true;
    const sub = grok.events.onViewerAdded.subscribe((a) => arm(a.args.viewer));
    try {
      tv.addViewer(DG.VIEWER.HISTOGRAM, {valueColumnName: 'age', splitColumnName: 'race'});
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      for (const v of tv.viewers) arm(v);
      await delay(100);
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'weight', min: 100, max: 170}, true);
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'race', selected: ['Black', 'Other']}, true);
      await awaitCheck(() => df.filter.trueCount < df.rowCount, 'the filters did not apply', 3000);
      const filtered = df.filter.trueCount;

      tv.loadLayout(tv.saveLayout());
      await awaitCheck(() => {
        const panel = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.FILTERS) as DG.FilterGroup | undefined;
        return panel?.filters.length === 2 && df.filter.trueCount === filtered;
      }, 'the layout came back with an empty filter panel', 3000);
    } finally {
      sub.unsubscribe();
      tv.close();
      grok.shell.closeTable(df);
    }
  });

  test('immediateRendering', async () => {
    const tv = grok.shell.addTableView(grok.data.demo.demog(100));
    const plot = tv.scatterPlot();
    let renders = 0;
    const sub = plot.onViewerRendered.subscribe(() => renders++);
    try {
      expect(plot.immediateRendering, false);
      await awaitCheck(() => renders > 0, 'plot did not paint', 3000);
      plot.immediateRendering = true;
      expect(plot.immediateRendering, true);
      renders = 0;
      plot.props.markerDefaultSize = 12;
      plot.props.markerOpacity = 50;
      await new Promise((r) => setTimeout(r, 0));
      expect(renders, 1);
    } finally {
      sub.unsubscribe();
      tv.close();
    }
  });
});
