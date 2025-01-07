import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, testViewer, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {readDataframe} from '../gui-utils';

category('Viewers: Scatter Plot', () => {
  test('Wrong range #1613', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('#1613.csv'), {readOnly: true});
  });

  test('Table is unexpectedly filtered after empty column creation #1744', async () => {
    const df = await readDataframe('#1744.csv');
    const tv = grok.shell.addTableView(df);
    const col = df.getCol('Col');
    for (let i = 0; i < 6; i++) col.set(i, null);
    tv.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: 'Col', xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    await awaitCheck(() => document.querySelector('[name=viewer-Scatter-plot i] canvas') !== null,
      'cannot load viewer', 2000);
    expect(df.filter.trueCount, 7);
  });

  test('Wrong min/max in log scale if near-zero values present #1764', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('#1764.csv'), {readOnly: true});
  });

  test('Wrong range #1858', async () => {
    await testViewer(DG.VIEWER.SCATTER_PLOT, await readDataframe('#1858.csv'), {readOnly: true});
  });
}, { owner: 'dkovalyov@datagrok.ai' });
