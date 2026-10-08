import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {SurfacePlot} from '../viewers/surface-plot/surface-plot';

category('Surface plot', () => {
  after(async () => grok.shell.closeAll());

  test('Datetime column on an axis', async () => {
    const df = await grok.data.getDemoTable('geo/earthquakes.csv');
    expect(df.columns.byIndex(0).type, DG.TYPE.DATE_TIME);
    const tv = grok.shell.addTableView(df);
    const viewer = new SurfacePlot();
    tv.addViewer(viewer);
    expect(viewer.option.xAxis3D.type, 'time');
    expect(viewer.option.xAxis3D.axisLabel.formatter === undefined);
    viewer.chart.resize();
  });
});
