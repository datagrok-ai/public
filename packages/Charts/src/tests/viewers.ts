import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, settled, test, testViewer} from '@datagrok-libraries/test/src/test';

const viewersToSkip: {[v: string]: string} = {
  // 'Globe': 'GROK-14320',
  'Surface plot':
    'GROK-18695: deterministic EXECUTION TIMEOUT in headless CI (echarts-gl WebGL); reproduces on echarts 5 and 6',
};

const chartsViewers = (): string[] =>
  DG.Func.find({package: 'Charts', meta: {role: DG.FUNC_TYPES.VIEWER}}).map((f) => f.friendlyName);

category('Viewers', () => {
  const df = grok.data.demo.demog(100);

  for (const v of chartsViewers()) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        if (['Sankey', 'Chord'].includes(v)) return (await grok.data.getDemoTable('energy_uk.csv'));
        else if (['Tree', 'Sunburst'].includes(v)) return (await grok.data.getDemoTable('demog.csv'));
        else if (v === 'Globe') return (await grok.data.getDemoTable('geo/earthquakes.csv'));
        return df.clone();
      })(), {detectSemanticTypes: true, arbitraryDfTest: v === 'Globe' ? false : undefined});
    }, v in viewersToSkip ? {skipReason: viewersToSkip[v]} : {});
  }
});

// The core viewers get the same two changes in the Dart "Viewers | Data Manipulation" tests.
category('Viewers: no rows', () => {
  for (const v of chartsViewers()) {
    test(v, async () => {
      const df = ['Sankey', 'Chord'].includes(v) ? await grok.data.getDemoTable('energy_uk.csv') :
        v === 'Globe' ? await grok.data.getDemoTable('geo/earthquakes.csv') : grok.data.demo.demog(100);
      const tv = grok.shell.addTableView(df);
      try {
        const viewer = tv.addViewer(v);
        await settled(viewer);
        grok.shell.clearLastError();
        df.filter.setAll(false);
        await settled(viewer);
        const shown = viewer.getWidgetStatus().values?.['rows shown'];
        if (shown !== undefined)
          expect(shown, 0);
        df.filter.setAll(true);
        df.rows.removeWhere(() => true);
        await settled(viewer);
        expect(await grok.shell.lastError ?? '', '');
      } finally {
        tv.close();
        grok.shell.closeTable(df);
      }
    }, v in viewersToSkip ? {skipReason: viewersToSkip[v]} : {});
  }
});
