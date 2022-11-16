import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {before, after, category, expect, test, delay} from '@datagrok-libraries/utils/src/test';
import { DataFrame } from 'datagrok-api/dg';

import './gui-utils';
import '../gis-semtypes';
import { SEMTYPEGIS } from '../gis-semtypes';


category('GIS-viewer', async () => {
  let testDF: DG.DataFrame | null = null;

  before(async () => {
    testDF = grok.data.demo.geo(30);

    // await grok.functions.call('');
    // await delay(9000);
    // const v = grok.shell.getTableView('');
  });

  test('GIS.openMapVewer', async () => {
    const v = grok.shell.getTableView('Map');
    expect(v.name === 'Map', true);
    // expect(allViews.every(item => grok.shell.view(item) !== undefined), true);
    // expect(allTableViews.every(item => grok.shell.view(item) !== undefined), true);
  });

  test('GIS.detectLonLat', async () => {
    if (testDF) {
      let col = testDF.columns.byName('Lon');
      await grok.data.detectSemanticTypes(testDF);
      expect(col.semType, SEMTYPEGIS.LONGITUDE);
      col = testDF.columns.byName('Lat');
      expect(col.semType, SEMTYPEGIS.LATIITUDE);
    }
  });

  after(async () => {
    grok.shell.closeTable(testDF as DataFrame);
    grok.shell.closeAll();
  });
});
