import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, testViewer, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from '../gui-utils';

category('Viewers: Scatter Plot', () => {
  test('NX #1613', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('scatter.plot.broken.csv'), {readOnly: true});
  });

  test('NX #1744', async () => {
    const df = await readDataframe('fruits_.csv');
    const col = df.getCol('Col');
    const tv = grok.shell.addTableView(df);
    for (let i = 0; i < 6; i++) col.set(i, null);
    tv.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: 'Col', xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    await awaitCheck(() => document.querySelector('[name=viewer-Scatter-plot i] canvas') !== null,
      'cannot load viewer', 2000);
    expect(df.filter.trueCount, 7);
  });

  test('NX #1764', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('BrokenScatterPlots.csv'), {readOnly: true});
  });

  test('NX #1858', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('scatter.plot.broken2.csv'), {readOnly: true});
  });
});
