import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Source: core/client/d4/lib/src/viewers/bar_chart/bar_chart_core.dart:298
// (canShowLegend) and :69 (legendCol => stackCol when stackCol is
// categorical/datetime). When stackColumnName == splitColumnName, the
// legend would duplicate the categorical X axis, so the auto legend is
// suppressed. Logic: canShowLegend() returns false iff
//   legendVisibility == Auto && sameLegendColumnDisplayed(splitCol)
// where sameLegendColumnDisplayed(col) checks col == legendCol == stackCol.
// Test relies on real attach + repaint: detached viewers skip legend layout
// (see _refreshLegend early-return on dataFrame == null), so we use a real
// TableView and a 300ms delay to let the deferred legend refresh land.
// Same-column case: .d4-legend should be removed from the DOM
// (see legend.remove() at legend_mixin.dart:530 and legend.dart:168
// 'root?.remove()'). Different-column case: .d4-legend should be present.
category('AI: Viewers: BarChart legend visibility', () => {
  test('legend hidden when stack column equals split column', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART,
        {value: 'age', split: 'race', stack: 'race'}) as DG.BarChartViewer;
      expect(v instanceof DG.BarChartViewer, true);
      const look = v.getOptions(true).look;
      expect(look['splitColumnName'], 'race');
      expect(look['stackColumnName'], 'race');
      await DG.delay(300);
      const legendEl = v.root.querySelector('.d4-legend');
      expect(legendEl == null, true);
    }
    finally {
      tv.close();
    }
  });

  test('legend visible when stack and split columns differ', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART,
        {value: 'age', split: 'race', stack: 'sex'}) as DG.BarChartViewer;
      expect(v instanceof DG.BarChartViewer, true);
      const look = v.getOptions(true).look;
      expect(look['splitColumnName'], 'race');
      expect(look['stackColumnName'], 'sex');
      await DG.delay(300);
      const legendEl = v.root.querySelector('.d4-legend');
      expect(legendEl != null, true);
    }
    finally {
      tv.close();
    }
  });

  test('legend toggles when stack switches to/from split column', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART,
        {value: 'age', split: 'race', stack: 'sex'}) as DG.BarChartViewer;
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') != null, true);

      v.setOptions({stackColumnName: 'race'});
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') == null, true);

      v.setOptions({stackColumnName: 'sex'});
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') != null, true);
    }
    finally {
      tv.close();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
