import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import {before, after, category, expect, test} from '@datagrok-libraries/utils/src/test';

import * as GUIUTILS from './gui-utils';
import {SEMTYPEGIS} from '../gis-semtypes';
import {GisViewer} from '../gis-viewer';


category('GIS: MapViewer', async () => {
  let testDF: DG.DataFrame; // | null = null;
  let testDF2: DG.DataFrame; // | null = null;

  before(async () => {
    testDF = grok.data.demo.geo(300);
    testDF2 = DG.DataFrame.fromCsv(await GUIUTILS.loadFileAsText('gistesttable.csv'));

    // await grok.functions.call('');
    // await delay(9000);
    // const v = grok.shell.getTableView('');
  });

  test('GIS:openMapViewer', async () => {
    const mapViewer = DG.Viewer.fromType('Map', testDF);
    expect(mapViewer instanceof DG.JsViewer, true);
    expect(mapViewer.type, 'Map');
    expect(mapViewer.table.id, testDF.id);

    // const v = grok.shell.getTableView('Map');
    // expect(v.name === 'Map', true);
    // expect(allViews.every(item => grok.shell.view(item) !== undefined), true);
    // expect(allTableViews.every(item => grok.shell.view(item) !== undefined), true);
  });

  test('GIS:viewerProperties', async () => {
    const mapViewer = (DG.Viewer.fromType('Map', testDF) as GisViewer);
    let options = await GUIUTILS.getOptions(mapViewer);
    expect(options.markerOpacity, 80);
    expect(options.markerDefaultSize, 5);
    expect(options.defaultColor, 0x1f77b4);
    expect(options.selectedColor, 0xff8c00);
    // mapViewer.markerDefaultSize = 13;
    // mapViewer.defaultColor = 0xff00ff;
    mapViewer.setOptions({markerDefaultSize: 13, defaultColor: 0xff00ff});
    //delay(200);
    options = await GUIUTILS.getOptions(mapViewer);
    expect(options.markerDefaultSize, 13);
    expect(options.defaultColor, 0xff00ff);
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

  test('GIS:detectAddress', async () => {
    if (!testDF2)
      testDF2 = DG.DataFrame.fromCsv(await GUIUTILS.loadFileAsText('gistesttable.csv'));
    if (testDF2) {
      await grok.data.detectSemanticTypes(testDF2);
      let col = testDF2.columns.byName('Street Address');
      expect(col.semType, SEMTYPEGIS.GISADDRESS);
      col = testDF2.columns.byName('Store Name');
      expect(false, col.semType === SEMTYPEGIS.GISADDRESS);
    }
  });

  test('GIS:detectCountry', async () => {
    if (!testDF2)
      testDF2 = DG.DataFrame.fromCsv(await GUIUTILS.loadFileAsText('gistesttable.csv'));
    if (testDF2) {
      await grok.data.detectSemanticTypes(testDF2);
      let col = testDF2.columns.byName('Country');
      expect(col.semType, SEMTYPEGIS.GISCOUNTRY);
      col = testDF2.columns.byName('City');
      expect(false, col.semType === SEMTYPEGIS.GISCOUNTRY);
    }
  });

  test('GIS:detectGisZipcode', async () => {
    if (!testDF2)
      testDF2 = DG.DataFrame.fromCsv(await GUIUTILS.loadFileAsText('gistesttable.csv'));
    if (testDF2) {
      await grok.data.detectSemanticTypes(testDF2);
      let col = testDF2.columns.byName('Postcode');
      expect(col.semType, SEMTYPEGIS.GISZIPCODE);
      col = testDF2.columns.byName('Store Number');
      expect(false, col.semType === SEMTYPEGIS.GISZIPCODE);
    }
  });

  after(async () => {
    grok.shell.closeTable(testDF as DG.DataFrame);
    grok.shell.closeTable(testDF2 as DG.DataFrame);
    grok.shell.closeAll();
  });
});
