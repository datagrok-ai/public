/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';
import {gisDemo} from './gis-demo';
import * as GisTypes from '../src/gis-semtypes';
import {OpenLayers, getKMZData} from '../src/gis-openlayer';
import {useGeographic} from 'ol/proj';

import {SEMTYPEGIS} from '../src/gis-semtypes';

export * from './package.g';
export const _package = new DG.Package();

function uiCensusDialog(): DG.Dialog {
  const data = grok.data.demo.demog(5);
  const data2 = grok.data.demo.demog(10);

  const reslutGridStyle = {
    'allowEdit': false,
    'showAddNewRowIcon': true,
    'allowColReordering': false,
    'allowColResizing': false,
    'allowColHeaderResizing': false,
    'allowColSelection': false,
    'allowRowResizing': false,
    'allowRowSelection': false,
    'allowBlockSelection': false,
    'autoScrollColumnIntoView': false,
    'showCurrentRowIndicator': false,
    'showCurrentCellOutline': false,
    'showRowHeader': false,
    'showColumnTooltip': false,
    'showColumnGridlines': false,
    'topLevelDefaultMenu': true,
  };

  const inputGridStyle = {
    'allowEdit': false,
    'showAddNewRowIcon': true,
    'allowColReordering': false,
    'allowColResizing': false,
    'allowColHeaderResizing': false,
    'allowColSelection': false,
    'allowRowResizing': false,
    'allowRowSelection': false,
    'allowBlockSelection': false,
    'autoScrollColumnIntoView': false,
    'showCurrentRowIndicator': true,
    'showCurrentCellOutline': false,
    'showRowHeader': false,
    'showColumnTooltip': false,
    'showColumnGridlines': false,
    'topLevelDefaultMenu': true,
  };

  const inputGrid = DG.Viewer.grid(data2, inputGridStyle);

  const inpSearch: DG.InputBase<string> = ui.input.string('Search', {value: ''});
  const inpLocation: DG.InputBase<string> = ui.input.string('Location', {value: ''});
  const choinpPeriod: DG.InputBase<string | null> = ui.input.choice('Period', {value: '', items: ['1990', '1991']});
  const txtDatasets = ui.divText('Dataset', {style: {color: 'var(--grey-4)', marginTop: '5px'}});
  const arr = [inpSearch, inpLocation, choinpPeriod, txtDatasets, inputGrid.root];

  const inputs = ui.div(arr, 'ui-form-condensed');
  const title = ui.h2('1990 Population Esimates - 1990-2000 Intercensal Esimates: United States, Nevada');

  const description = ui.p('Monthly Intercensal Esimates...');
  description.style.cssText = `
  overflow: hidden;
  display: -webkit-box;
  -webkit-line-clamp: 3;
  -webkit-box-orient: vertical;
  margin: 10px 0;
  `;

  const gridControl = ui.divH([
    ui.link('All', ()=>{grok.shell.info('Select all');}, 'Select all', {style: {marginRight: '10px'}}),
    ui.link('None', ()=>{grok.shell.info('Deselect all');}, 'Deselect all', ''),
    ui.divText('N Cheked'),
  ]);
  (gridControl.children[1] as HTMLElement).style.margin = '0 10px';
  gridControl.style.bottom = '5px';
  gridControl.style.position = 'Absolute';

  const resultGrid = DG.Viewer.grid(data, reslutGridStyle);
  resultGrid.columns.setVisible(['site', 'race']);

  const result = ui.splitV([
    ui.divV([
      title,
      description,
      ui.link('Read more', '#', 'Click to read more about this dataset'),
      gridControl,
    ]),
    ui.box(resultGrid.root),
  ], {style: {marginLeft: '24px'}});

  const root = ui.splitH([
    ui.box(inputs, {style: {maxWidth: '230px'}}),
    result]);

  const dlg = ui.dialog({title: 'Windows', helpUrl: '/help/transform/add-new-column.md'}) //TODO: place my help here
    .add(root)
    .show({width: 700, height: 500, center: true});

  dlg.root.querySelector('.d4-dialog-contents')?.classList.remove('ui-panel');
  dlg.root.querySelector('.d4-dialog-contents')?.classList.add('ui-box');

  return dlg;
}


async function fetchCensus() {
  const url = 'https://api.census.gov/data/';
  let censusRes: any = null;
  try {
    censusRes = await (await grok.dapi.fetchProxy(url)).json();
  } catch (e: any) {
    grok.shell.error(`Census fetch error: ${e.message}`);
    censusRes = '';
  } finally {
    return censusRes;
  }
}

async function handleKMLOrKMZFile(filecontent: string | Uint8Array, isKmz: boolean): Promise<DG.DataFrame[]> {
  let dfFromKML: DG.DataFrame | undefined = undefined;
  const ol = new OpenLayers();

  let kmlData = '';
  if (isKmz) kmlData = await getKMZData(filecontent);
  else kmlData = filecontent as string;
  const newLayer = ol.addKMLLayerFromStream(kmlData);
  const arrFeatures = ol.exportLayerToArray(newLayer);

  if (arrFeatures) {
    if (arrFeatures.length > 0) {
      dfFromKML = DG.DataFrame.fromObjects(arrFeatures);
      if (dfFromKML) {
        const gisCol = dfFromKML.col('gisObject');
        if (gisCol)
          gisCol.semType = SEMTYPEGIS.GISAREA; //SEMTYPEGIS.GISOBJECT;

        dfFromKML.name = newLayer.get('layerName');

        const tableView = grok.shell.addTableView(dfFromKML);
        tableView.name = 'dfFromKML.name' + ' (manual)';

        setTimeout((tableView: DG.TableView, kmlData: string) => {
          const v = tableView;
          const mapViewer = v.addViewer(new GisViewer()) as GisViewer;

          mapViewer.ol.initMap('map-container');
          mapViewer.ol.addKMLLayerFromStream(kmlData);
        }, 1000, tableView, kmlData); //should use it because we need visible and active View for map initializing
      }
    }
  }
  if (dfFromKML) return [dfFromKML];
  return [];
}
function gisGeoJSONFileDetector(strBuf: string): [boolean, boolean] {
  let arrTmp: any[] | null;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?[''|"](?:Multi)?Polygon/ig);
  let cntTypeGeo = arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?\'|\'(?:Multi)?Point/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?\'|\'(?:Multi)?LineString/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?[''|"]Feature/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?[''|"]GeometryCollection/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?[''|"]Topology[''|"]/ig);
  const cntTypeTopo = arrTmp ? arrTmp.length : 0;

  if (cntTypeGeo === 0) return [false, false];
  if (cntTypeTopo > 0) return [true, true];

  return [true, false];
}

function openGISViewer(): void {
  let htmlStyle: DG.ElementOptions = { };
  htmlStyle = {style: {'width': '100%', 'height': '100%'}};
  const boxMap = ui.box(null, htmlStyle);
  boxMap.id = 'map-container'; //boxMap - div that contains map
  grok.shell.newView('OpenLayers preview view', [boxMap]);

  const ol = new OpenLayers();
  ol.initMap('map-container');
}


export class PackageFunctions {
  @grok.decorators.func()
  static async getCensusInfo() {
    let htmlStyle: DG.ElementOptions = { };
    let censusRes: any = null;
    const mapVintages = new Map<string, any[]>();
    let infoDataset: HTMLElement | null = null;

    censusRes = await fetchCensus();

    //TODO: save fetch result into buffer to prevent frequent uploading
    if (!censusRes)
      return 'Fetch error';
    const catalogCensus = censusRes.dataset;
    for (let i = 0; i < catalogCensus.length; i++) {
      if (catalogCensus[i].c_vintage) {
        if (!mapVintages.has(catalogCensus[i].c_vintage)) mapVintages.set(catalogCensus[i].c_vintage, []);
        mapVintages.get(catalogCensus[i].c_vintage)?.push(catalogCensus[i]);
      }
    }

    //create dialog UI>>
    const colVintages = DG.Column.fromList('int', 'Vintage', [...mapVintages.keys()]);
    const dfVintages = DG.DataFrame.fromColumns([colVintages]);
    //@ts-ignore until we fix it
    const gridVintages = DG.Viewer.grid(dfVintages, {'min-width': '70px', 'border': '1px green'});
    gridVintages.autoSize(68, 300);
    const col1 = gridVintages.columns.byName('Vintage');
    if (col1)
      col1.width = 55;
    if (gridVintages.columns.rowHeader)
      gridVintages.columns.rowHeader.visible = false;

    const dfDatasets = DG.DataFrame.create(0);
    //@ts-ignore
    const gridDatasets = DG.Viewer.grid(dfDatasets, {'min-width': '100%', 'border': 'solid 1px red'});
    gridDatasets.autoSize(400, 300, 400, 300);
    const dfVariables = DG.DataFrame.create(0);
    //@ts-ignore
    const gridVariables = DG.Viewer.grid(dfVariables, {'min-width': '170px', 'border': 'solid 1px yellow'});
    gridVariables.autoSize(340, 300, 340, 300);

    //Vintage selection handler>>
    gridVintages.onCurrentCellChanged.subscribe((ev) => {
      const cellVal = (ev.cell.value) as string;
      if (cellVal) {
        const datasets = mapVintages.get(cellVal);
        if (datasets) {
          const df = DG.DataFrame.fromObjects(datasets);
          if (df)
            gridDatasets.dataFrame = df;
          gridDatasets.autoSize(400, 300, 400, 300);
          gridDatasets.columns.setVisible(['title']);
          if (gridDatasets.columns.rowHeader)
            gridDatasets.columns.rowHeader.visible = false;
          const col1 = gridDatasets.columns.byName('title');
          if (col1)
            col1.name = 'Dataset name';
        }
      }
    }); //<<vintage selection handler

    //Dataset selection handler>>
    gridDatasets.onCurrentCellChanged.subscribe(async (ev) => {
      const descr = ev.cell.dataFrame.get('description', ev.cell.rowIndex);
      if (infoDataset && descr)
        infoDataset.innerText = descr;
      //get list of variables for dataset
      const urlVariables = ev.cell.dataFrame.get('c_variablesLink', ev.cell.rowIndex);
      const variables = await(await grok.dapi.fetchProxy(urlVariables)).json();
      if (variables) {
        const varList = [];
        for (const v in variables['variables']) {
          if (variables['variables'].hasOwnProperty(v)) {
            const varObj = variables['variables'][v];
            varObj.varname = v;
            varObj.use = false;
            varList.push(varObj);
          }
        }
        // const df = DG.DataFrame.fromJson(JSON.stringify(variables));
        const df = DG.DataFrame.fromObjects(varList);
        if (df) gridVariables.dataFrame = df;
        if (gridVariables.columns.rowHeader) gridVariables.columns.rowHeader.visible = false;
        gridVariables.columns.setOrder(['use', 'label', 'varname']);
        gridVariables.columns.setVisible(['use', 'label', 'varname']);
        let col1 = gridVariables.columns.byName('label');
        if (col1) col1.width = 190;
        col1 = gridVariables.columns.byName('use');
        if (col1) col1.width = 30;
      }
    }); //<<dataset selection handler

    // htmlStyle = {style: {'width': '200px', 'border-right': 'solid 1px darkgray', 'min-width': '150px',
    //   'margin-right': '10px', 'padding-right': '2px'}};
    htmlStyle = {style: {'border': 'solid 1px lightgray', 'width': '75px'}};
    const panelVintages = ui.box(null, htmlStyle);
    panelVintages.style.minWidth = '75px';
    panelVintages.style.maxWidth = '75px';
    panelVintages.append(gridVintages.root);
    htmlStyle = {style: {'border': 'solid 1px lightgray', 'borderLeft': 'none', 'width': '408px', 'minWidth': '408px'}};
    const panelDatasets = ui.box(null, htmlStyle);
    panelDatasets.style.minWidth = '400px';
    panelDatasets.append(gridDatasets.root);
    htmlStyle = {style: {'border': 'solid 1px lightgray', 'borderLeft': 'none', 'width': '405px', 'minWidth': '402px'}};
    const panelVariables = ui.box(null, htmlStyle);
    panelVariables.style.minWidth = '250';
    panelVariables.append(gridVariables.root);
    htmlStyle = {style: {'width': '100%', 'minHeight': '305px'}};
    const gridsHolder = ui.splitH([
      panelVintages,
      panelDatasets,
      panelVariables,
    ], htmlStyle);

    htmlStyle = {style: {'border': 'none', 'width': '100%', 'minHeight': '95px'}};
    infoDataset = ui.box(null, htmlStyle);
    // @ts-ignore
    infoDataset.style.height = '100%';

    htmlStyle = {style: {'width': '100%'}};
    ui.dialog('US Census bureau data import')
      .add(ui.splitV([
        gridsHolder,
        infoDataset,
      ], htmlStyle))
      .onOK(async (event: any) => {
        const linkStr = (gridDatasets.dataFrame.get('c_valuesLink', gridDatasets.dataFrame.currentRowIdx) as string);
        const endInd = linkStr.lastIndexOf('/');
        let url = linkStr.slice(0, endInd);
        let varList = '?get=';
        for (let i = 0; i < gridVariables.dataFrame.rowCount; i++ ) {
          if (gridVariables.dataFrame.get('use', i)) {
            if ((gridVariables.dataFrame.get('varname', i) != 'for') &&
            (gridVariables.dataFrame.get('varname', i) != 'in')) {
              if (varList === '?get=') varList = varList + (gridVariables.dataFrame.get('varname', i) as string);
              else varList = varList + ',' + (gridVariables.dataFrame.get('varname', i) as string);
            }
          }
        }
        // url = url + varList + '&for=county:*&in=state:02';  //as an example
        url = url + varList + '&for=state:*';
        url = url + '&key=2647d704d8734665d5c417dae1546887c2c90513'; //TODO: we need to hide key in credentials
        try {
          censusRes = await(await grok.dapi.fetchProxy(url)).text();
          if (censusRes.toLowerCase().includes('error'))
            throw new Error(censusRes);

          const df = DG.DataFrame.fromCsv(censusRes);
          grok.shell.addTableView(df);
        } catch (e) {
          grok.shell.error(`Census request error: ${e}`);
        }
      })
      .show({x: 200, y: 80, width: 850, height: 500}); //showModal()
  //<<end of getCensusInfo function
  }

  @grok.decorators.func({
    'name': 'getCensusInfoTest UI',
  })
  static async getCensusInfoT() {
    const dlg = uiCensusDialog();
  }

  @grok.decorators.func()
  static async info() {
    grok.shell.info('GIS Package info: ' +_package.webRoot);
  }


  @grok.decorators.init()
  static init() {
    //Register handlers
    DG.ObjectHandler.register(new GisTypes.GisPointHandler());
    DG.ObjectHandler.register(new GisTypes.GisAreaHandler());
  }


  @grok.decorators.func({
    'name': 'GISGeocoding',
    'description': 'GIS geocoding - receive coordinates for address',
  })
  static async gisGeocoding(
    address: string): Promise<string> {
    let url = 'https://geocoding.geo.census.gov/geocoder/';
    let fetchResult = null;
    // eslint-disable-next-line max-len
    url += `locations/onelineaddress?address=${address}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
    fetchResult = await (await grok.dapi.fetchProxy(url)).json();

    if (!fetchResult) return 'Fetch error';
    const resStr = JSON.stringify(fetchResult);

    return resStr;
  }


  @grok.decorators.func({
    'description': 'GIS geocoding - receive address for coordinates',
  })
  static async gisReverseGeocoding(
    x: number,
    y: number): Promise<string> {
    let url = 'https://geocoding.geo.census.gov/geocoder/';
    let fetchResult = null;
    url += `geographies/coordinates?x=${x}&y=${y}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
    fetchResult = await (await grok.dapi.fetchProxy(url)).json();

    if (!fetchResult)
      return 'Fetch error';

    const resStr = JSON.stringify(fetchResult);

    return resStr;
  }


  @grok.decorators.func({
    'description': 'GIS geocoding - receive coordinates from address',
  })
  static async gisBatchGeocoding(
    address: string): Promise<string> {
    let url = 'https://geocoding.geo.census.gov/geocoder/locations/addressbatch';
    let fetchResult = null;
    // url += `&vintage=ACS2021_Current&benchmark=Public_AR_Current`;
    url += `&benchmark=Public_AR_Current`;
    url += `&addressFile=1,4600 Silver Hill Road,Washington,DC,20233`;
    fetchResult = await (await grok.dapi.fetchProxy(url)).text();

    if (!fetchResult)
      return 'Fetch error';

    return fetchResult;
  }


  @grok.decorators.func({
    'meta': {
      'icon': 'icons/icon.svg',
      'toolbox': 'true',
      'role': 'viewer',
    },
    'name': 'Map',
    'description': 'GIS map viewer',
    'outputs': [{'name': 'result', 'type': 'viewer'}],
  })
  static gisViewer(): GisViewer {
    // setTimeout(() => {grok.shell.windows.showProperties = true;}, 500);
    return new GisViewer();
  }


  @grok.decorators.fileViewer({
    'fileViewer': 'kmz,kml',
  })
  static async gisKMZAndKMLFileViewer(
    file: DG.FileInfo): Promise<DG.View> {
    const viewFile = DG.View.create();
    viewFile.name = 'Preview of: ' + file.name;
    viewFile.root.id = 'map-container'; //boxMap - div that contains map
    viewFile.root.style.padding = '0px';

    let kmlData = '';
    if ((file.extension).toLowerCase() === 'kmz') {
      const strBuf = await file.readAsBytes();
      kmlData = await getKMZData(strBuf);
    } else if ((file.extension).toLowerCase() === 'kml')
      kmlData = await file.readAsString();

    setTimeout(() => {
      const ol = new OpenLayers();
      ol.initMap('map-container');
      useGeographic();
      ol.addKMLLayerFromStream(kmlData);
      ol.setViewOptions({projection: 'EPSG:3857'});
    }, 200); //should use it because we need visible and active View for map initializing

    return viewFile;
  }


  @grok.decorators.fileHandler({
    'ext': 'kml',
  })
  static async gisKMLFileHandler(
    filecontent: string) : Promise<DG.DataFrame[]> {
    return handleKMLOrKMZFile(filecontent, false);
  }


  @grok.decorators.fileHandler({
    'ext': 'kmz',
  })
  static async gisKMZFileHandler(
    @grok.decorators.param({'type': 'list'}) filecontent: Uint8Array): Promise<DG.DataFrame[]> {
    return handleKMLOrKMZFile(filecontent, true);
  }


  @grok.decorators.fileViewer({
    'fileViewer': 'geojson, topojson',
    'outputs': [{'name': 'result', 'type': 'view'}],
  })
  static async gisGeoJSONFileViewer(
    file: DG.FileInfo): Promise<DG.View | null | DG.DataFrame> {
    //read file
    const strBuf = await file.readAsString();
    const isGeoTopo = gisGeoJSONFileDetector(strBuf);

    if (isGeoTopo[0] === false) {
      //if json file is not kind of geoJson or topoJson - show it as a table
      const df = DG.DataFrame.fromJson(strBuf);
      const viewFile = DG.TableView.create(df);
      return viewFile;
      // return df;
    }

    const viewFile = DG.View.create();
    viewFile.name = 'Preview of: ' + file.name;
    viewFile.root.id = 'map-container';
    viewFile.root.style.padding = '0px';

    setTimeout(() => {
      const ol = new OpenLayers();
      ol.initMap('map-container');
      useGeographic();
      ol.setViewOptions({
        projection: 'EPSG:4326',
      });
      if (isGeoTopo[1] === false) ol.addGeoJSONLayerFromStream(strBuf);
      else ol.addTopoJSONLayerFromStream(strBuf);
    }, 200);

    return viewFile;
  }


  @grok.decorators.fileHandler({
    'ext': 'geojson, topojson',
  })
  static gisGeoJSONFileHandler(
    filecontent: string): DG.DataFrame[] {
    //detect the kind of json file
    const isGeoTopo = gisGeoJSONFileDetector(filecontent);
    if (isGeoTopo[0] === false) //if this a simple JSON - just return DF with JSON content
      return [DG.DataFrame.fromJson(filecontent)];

    let dfFromJSON: DG.DataFrame | undefined = undefined;
    const ol = new OpenLayers();
    let newLayer = null;
    if (isGeoTopo[1] === false)
      newLayer = ol.addGeoJSONLayerFromStream(filecontent);
    else
      newLayer = ol.addTopoJSONLayerFromStream(filecontent);

    //create the dataframe
    const arrFeatures = ol.exportLayerToArray(newLayer);
    if (arrFeatures) {
      if (arrFeatures.length > 0) {
        dfFromJSON = DG.DataFrame.fromObjects(arrFeatures);
        if (dfFromJSON) {
          const gisCol = dfFromJSON.col('gisObject');
          if (gisCol)
            gisCol.semType = SEMTYPEGIS.GISAREA; //SEMTYPEGIS.GISOBJECT;

          dfFromJSON.name = newLayer.get('layerName');
          const tableView = grok.shell.addTableView(dfFromJSON);
          tableView.name = 'dfFromJSON.name'; //TODO: add file name here instead of this hardcode

          //show the map viewer
          setTimeout((tableView: DG.TableView) => {
            const v = tableView;
            if (v) {
              const mapViewer = v.addViewer(new GisViewer()) as GisViewer;
              //initialize the map
              mapViewer.ol.initMap('map-container');

              if (isGeoTopo[1] === false)
                mapViewer.ol.addGeoJSONLayerFromStream(filecontent);
              else
                mapViewer.ol.addTopoJSONLayerFromStream(filecontent);
            }
          }, 1000, tableView);
        }
      }
    }
    if (dfFromJSON)
      return [dfFromJSON];

    return [DG.DataFrame.fromJson(filecontent)];
  }


  @grok.decorators.panel({
    'condition': 'true',
    'meta': {role: 'widgets'},
  })
  static gisAreaWidget(
    @grok.decorators.param({'type': 'string', 'options': {'semType': 'gis-area'}}) gisArea: any): DG.Widget | null {
    //this is temporary code - should be filled with usefull functionality
    // if ((!gisArea) || !(gisArea instanceof GisArea)) return null;
    if ((!gisArea)) return null;

    const strToAdd: string = 'test'; //(gisArea as GisArea).semtype;
    let widgetStyle: DG.ElementOptions = { };
    widgetStyle = {style: {'color': '#F55'}};

    return new DG.Widget(ui.divText('gis Area widget ' + strToAdd, widgetStyle));
  }


  @grok.decorators.demo({
    'meta': {
      'demoPath': 'Visualization | Geographical | Map',
    },
    'name': 'mapDemo',
    'description': 'Map viewer shows geospatial data on a map as either markers, or a heat map.',
    'test': {
      'test': '_mapDemo()',
      'wait': '2000',
    },
  })
  static async _mapDemo() {
    await gisDemo({renderType: 'both'});
  }
}
