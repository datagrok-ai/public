import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {before, after, category, expect, test, delay} from '@datagrok-libraries/utils/src/test';

import './gui-utils';
import '../gis-semtypes';
import {SEMTYPEGIS} from '../gis-semtypes';


category('GIS: MapViewer', async () => {
  let testDF: DG.DataFrame; // | null = null;

  before(async () => {
    testDF = grok.data.demo.geo(300);

    // await grok.functions.call('');
    // await delay(9000);
    // const v = grok.shell.getTableView('');
  });

  test('GIS:openMapVewer', async () => {
    const viewer = DG.Viewer.fromType('Map', testDF);
    expect(viewer instanceof DG.JsViewer, true);
    expect(viewer.type, 'Map');
    expect(viewer.table.id, testDF.id);

    // const v = grok.shell.getTableView('Map');
    // expect(v.name === 'Map', true);
    // expect(allViews.every(item => grok.shell.view(item) !== undefined), true);
    // expect(allTableViews.every(item => grok.shell.view(item) !== undefined), true);
  });

  test('GIS:detectLonLat', async () => {
    if (!testDF)
      testDF = grok.data.demo.geo(300);
    if (testDF) {
      await grok.data.detectSemanticTypes(testDF);
      let col = testDF.columns.byName('lng');
      expect(col.semType, SEMTYPEGIS.LONGITUDE);
      col = testDF.columns.byName('lat');
      expect(col.semType, SEMTYPEGIS.LATIITUDE);
    }
  });

  after(async () => {
    // grok.shell.closeTable(testDF as DataFrame);
    // grok.shell.closeAll();
  });
});
