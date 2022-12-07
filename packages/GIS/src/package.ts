/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';
import * as GisTypes from '../src/gis-semtypes';
//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer';
import {useGeographic} from 'ol/proj';
// import * as cns from '../node_modules/citysdk';

//GIS semantic types import
import {SEMTYPEGIS} from '../src/gis-semtypes';

const census = require('citysdk');

//ZIP utilities
import JSZip from 'jszip';
import {DataFrame, InputBase} from 'datagrok-api/dg';

//USEFUL
// contents = fs.readFileSync(detectorsPath, 'utf8');


export const _package = new DG.Package();

//census SDK async wrapper function
async function censusPromise(args: any) {
  return new Promise(function(resolve, reject) {
    census(args, function(err: any, json: any) {
      if (!err)
        resolve(json);
      else
        reject(err);
    });
  });
}

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

  const inpSearch: DG.InputBase<string> = ui.stringInput('Search', '');
  const inpLocation: DG.InputBase<string> = ui.stringInput('Location', '');
  const choinpPeriod: DG.InputBase<string | null> = ui.choiceInput('Period', '', ['1990', '1991']);
  const txtDatasets = ui.divText('Dataset', {style: {color: 'var(--grey-4)', marginTop: '5px'}});
  const arr = [inpSearch, inpLocation, choinpPeriod, txtDatasets, inputGrid.root];

  const inputs = ui.inputs(arr as Iterable<InputBase<any>>, 'ui-form-condensed');
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

//name: getCensusInfoTest UI
export async function getCensusInfoT() {
  const dlg = uiCensusDialog();
}

//name: getCensusInfo
export async function getCensusInfo() {
  let htmlStyle: DG.ElementOptions = { };
  let censusRes: any = null;
  const mapVintages = new Map<string, any[]>();
  let infoDataset: HTMLElement | null = null;

  //TODO: add try-catch block for fetch
  const url = 'https://api.census.gov/data/';
  censusRes = await (await grok.dapi.fetchProxy(url)).json();

  //TODO: put fetch into separate function
  //TODO: save fetch result into buffer to prevent frequent uploading
  if (!censusRes) return 'Fetch error';
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
  const gridVintages = DG.Viewer.grid(dfVintages, {'min-width': '70px', 'border': '1px green'});
  gridVintages.autoSize(68, 300);
  const col1 = gridVintages.columns.byName('Vintage');
  if (col1) col1.width = 55;
  if (gridVintages.columns.rowHeader) gridVintages.columns.rowHeader.visible = false;
  const dfDatasets = DG.DataFrame.create(0);
  const gridDatasets = DG.Viewer.grid(dfDatasets, {'min-width': '100%', 'border': 'solid 1px red'});
  gridDatasets.autoSize(400, 300, 400, 300);
  const dfVariables = DG.DataFrame.create(0);
  const gridVariables = DG.Viewer.grid(dfVariables, {'min-width': '170px', 'border': 'solid 1px yellow'});
  gridVariables.autoSize(340, 300, 340, 300);

  //Vintage selection handler>>
  gridVintages.onCurrentCellChanged.subscribe((ev) => {
    const cellVal = (ev.cell.value) as string;
    if (cellVal) {
      const datasets = mapVintages.get(cellVal);
      if (datasets) {
        const df = DG.DataFrame.fromObjects(datasets);
        if (df) gridDatasets.dataFrame = df;
        gridDatasets.autoSize(400, 300, 400, 300);
        gridDatasets.columns.setVisible(['title']);
        if (gridDatasets.columns.rowHeader) gridDatasets.columns.rowHeader.visible = false;
        const col1 = gridDatasets.columns.byName('title');
        if (col1) col1.name = 'Dataset name';
      }
    }
  }); //<<vintage selection handler

  //Dataset selection handler>>
  gridDatasets.onCurrentCellChanged.subscribe(async (ev) => {
    const descr = ev.cell.dataFrame.get('description', ev.cell.rowIndex);
    if (infoDataset && descr) infoDataset.innerText = descr;
    //get list of variables for dataset
    const urlVariables = ev.cell.dataFrame.get('c_variablesLink', ev.cell.rowIndex);
    const variables = await(await grok.dapi.fetchProxy(urlVariables)).json();
    if (variables) {
      const varList = [];
      for (let v in variables['variables']) {
        let varObj = variables['variables'][v];
        varObj.varname = v;
        varObj.use = false;
        varList.push(varObj);
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
  htmlStyle = {style: {'border': 'solid 1px lightgray', 'border-left': 'none', 'width': '408px', 'min-width': '408px'}};
  const panelDatasets = ui.box(null, htmlStyle);
  panelDatasets.style.minWidth = '400px';
  panelDatasets.append(gridDatasets.root);
  htmlStyle = {style: {'border': 'solid 1px lightgray', 'border-left': 'none', 'width': '405px', 'min-width': '402px'}};
  const panelVariables = ui.box(null, htmlStyle);
  panelVariables.style.minWidth = '250';
  panelVariables.append(gridVariables.root);
  htmlStyle = {style: {'width': '100%', 'minHeight': '305px'}};
  const gridsHolder = ui.splitH([
    panelVintages,
    panelDatasets,
    panelVariables,
  ], htmlStyle);

  htmlStyle = {style: {'border': 'none', 'width': '100%', 'min-height': '95px'}};
  infoDataset = ui.box(null, htmlStyle);
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
      // url = url + varList + '&for=county:*&in=state:02';
      url = url + varList + '&for=state:*';
      url = url + '&key=2647d704d8734665d5c417dae1546887c2c90513';
      try {
        censusRes = await(await grok.dapi.fetchProxy(url)).text();
        if (censusRes.toLowerCase().includes('error')) throw new Error(censusRes);
        // alert(censusRes);
        // const jsonData = JSON.parse(censusRes);
        // const df = DataFrame.fromJson(JSON.stringify(censusRes));
        // const df = DataFrame.fromJson(jsonData);
        const df = DataFrame.fromCsv(censusRes);
        grok.shell.addTableView(df);
      } catch (e) {
        grok.shell.error(`Census request error: ${e}`);
      }
    })
    .show({x: 200, y: 80, width: 850, height: 500}); //showModal()

//<<end of getCensusInfo function
}

//name: info
export async function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);
}

//tags: init
export function init() {
  //Register handlers
  DG.ObjectHandler.register(new GisTypes.GisPointHandler());
  DG.ObjectHandler.register(new GisTypes.GisAreaHandler());
}

//name: GISGeocoding
//description: GIS geocoding - receive coordinates for address
//input: string address
//output: string result
export async function gisGeocoding(address: string): Promise<string> {
// export async function gisGeocoding(address?: string, x?: number, y?: number): Promise<string> {
  let url = 'https://geocoding.geo.census.gov/geocoder/';
  let fetchResult = null;
  // eslint-disable-next-line max-len
  url += `locations/onelineaddress?address=${address}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
  fetchResult = await (await grok.dapi.fetchProxy(url)).json();
  // if ((address) && (address != '')) {
  //   // eslint-disable-next-line max-len
  //   url += `locations/onelineaddress?address=${address}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
  //   fetchResult = await (await grok.dapi.fetchProxy(url)).json();
  // } else if ((x) && (y)) {
  //   url += `geographies/coordinates?x=${x}&y=${y}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
  //   fetchResult = await (await grok.dapi.fetchProxy(url)).json();
  // }

  if (!fetchResult) return 'Fetch error';
  alert(JSON.stringify(fetchResult));
  // const df = DG.DataFrame.fromJson(fetchResult.addressMatches);
  // const viewFile = DG.TableView.create(df);
  const resStr = JSON.stringify(fetchResult);

  return resStr;
}

//name: gisReverseGeocoding
//description: GIS geocoding - receive address for coordinates
//input: double x
//input: double y
//output: string result
export async function gisReverseGeocoding(x: number, y: number): Promise<string> {
  let url = 'https://geocoding.geo.census.gov/geocoder/';
  let fetchResult = null;
  url += `geographies/coordinates?x=${x}&y=${y}&vintage=Census2020_Current&benchmark=Public_AR_Current&format=json`;
  fetchResult = await (await grok.dapi.fetchProxy(url)).json();

  if (!fetchResult) return 'Fetch error';
  alert(JSON.stringify(fetchResult));
  // const df = DG.DataFrame.fromJson(fetchResult);
  // const viewFile = DG.TableView.create(df);
  const resStr = JSON.stringify(fetchResult);

  return resStr;
}

//name: gisBatchGeocoding
//description: GIS geocoding - receive coordinates from address
//input: string address
//output: string result
export async function gisBatchGeocoding(address: string): Promise<string> {
  let url = 'https://geocoding.geo.census.gov/geocoder/locations/addressbatch';
  let fetchResult = null;
  // url += `&vintage=ACS2021_Current&benchmark=Public_AR_Current`;
  url += `&benchmark=Public_AR_Current`;
  url += `&addressFile=1,4600 Silver Hill Road,Washington,DC,20233`;
  fetchResult = await (await grok.dapi.fetchProxy(url)).text();

  if (!fetchResult) return 'Fetch error';
  alert((fetchResult));
  // const df = DG.DataFrame.fromJson(fetchresult);
  // const viewFile = DG.TableView.create(df);
  // const resStr = JSON.stringify(fetchResult);

  // census(
  //   {
  //     statsKey: '2647d704d8734665d5c417dae1546887c2c90513', //TODO: hide credentials
  //     vintage: '2017',
  //     geoHierarchy: {
  //       state: {
  //         lat: 38.8482,
  //         lng: -76.9312,
  //       },
  //       county: '*',
  //     },
  //     sourcePath: [
  //       'acs',
  //       'acs5',
  //     ],
  //     // values: [
  //     //   'B00001_001E',
  //     // ],
  //     geoResolution: '20m',
  //   },
  //   (err: any, res: any) => {
  //     if (!err)
  //       alert(err);
  //     else
  //       console.log(res);
  //   }
  // );

  // const df = DG.DataFrame.fromJson(fetchresult);
  // const viewFile = DG.TableView.create(df);

  return '';
}

//name: Map
//description: GIS map viewer
//tags: viewer
//meta.icon: icons/package.png
//output: viewer result
export function gisViewer(): GisViewer {
  // setTimeout(() => {grok.shell.windows.showProperties = true;}, 500);
  return new GisViewer();
}

async function getKMZData(buffer: any): Promise<string> {
  const zip = new JSZip();
  let kmlData: string = '';
  await zip.loadAsync(buffer);
  const kmlFile = zip.file(/.kml$/i)[0];
  if (kmlFile)
    kmlData = await kmlFile.async('string');
  return kmlData;
}

//name: gisKMZFileViewer
//tags: fileViewer, fileViewer-kmz, fileViewer-kml
//input: file file
//output: view result
export async function gisKMZFileViewer(file: DG.FileInfo): Promise<DG.View> {
  const viewFile = DG.View.create();
  viewFile.name = 'Preview of: ' + file.name;
  viewFile.root.id = 'map-cont'; //boxMap - div that contains map
  viewFile.root.style.padding = '0px';

  let kmlData = '';
  if ((file.extension).toLowerCase() === 'kmz') {
    const strBuf = await file.readAsBytes();
    kmlData = await getKMZData(strBuf);
  }
  else if ((file.extension).toLowerCase() === 'kml')
    kmlData = await file.readAsString();

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-cont');
    useGeographic();
    ol.addKMLLayerFromStream(kmlData);
    ol.setViewOptions({projection: 'EPSG:3857'});
  }, 200); //should use it because we need visible and active View for map initializing

  return viewFile;
}

//name: gisKMLFileHandler
//tags: file-handler
//meta.ext: kmz
//input: string filecontent
//output: list tables
export function gisKMLFileHandler(filecontent: string): DG.DataFrame[] {
  //TODO: detect the kind of file and join KML/KMZ handler

  let dfFromKML: DG.DataFrame | undefined = undefined;
  const ol = new OpenLayers();

  const kmlData = filecontent;
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
        const tv = grok.shell.addTableView(dfFromKML as DG.DataFrame);
        tv.name = 'dfFromKML.name' + ' (manual)';

        const mapViewer = tv.addViewer(DG.Viewer.fromType('Map', dfFromKML)) as GisViewer;
        if (mapViewer) {
          (mapViewer as GisViewer).init();
          (mapViewer as GisViewer).ol.addLayer(newLayer);
        }
      }
    }
  }
  if (dfFromKML) return [dfFromKML];
  return [];
}

//name: gisKMZFileHandler
//tags: file-handler
//meta.ext: kmz
//input: list bytes
//output: list tables
export async function gisKMZFileHandler(filecontent: Uint8Array): Promise<DG.DataFrame[]> {
  // export function gisKMLFileHandler(filecontent: string): DG.DataFrame[] {
  //detect the kind of file

  let dfFromKML: DG.DataFrame | undefined = undefined;
  const ol = new OpenLayers();

  const kmlData = await getKMZData(filecontent);
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
        const tv = grok.shell.addTableView(dfFromKML as DG.DataFrame);
        tv.name = 'dfFromKML.name' + ' (manual)';

        const mapViewer = tv.addViewer(DG.Viewer.fromType('Map', dfFromKML)) as GisViewer;
        if (mapViewer) {
          (mapViewer as GisViewer).init();
          (mapViewer as GisViewer).ol.addLayer(newLayer);
        }
      }
    }
  }
  if (dfFromKML) return [dfFromKML];
  return [];
}

//gisGeoJSONFileDetector
function gisGeoJSONFileDetector(strBuf: string): [boolean, boolean] {
  let arrTmp: any[] | null;
  //searching for patterns of GeoJSON data
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
  //searching for patterns of TopoJSON data
  arrTmp = strBuf.match(/[''|"]type[''|"]\s?:\s?[''|"]Topology[''|"]/ig);
  const cntTypeTopo = arrTmp ? arrTmp.length : 0;

  if (cntTypeGeo === 0) return [false, false];
  if (cntTypeTopo > 0) return [true, true];

  return [true, false];
}

//name: gisGeoJSONFileViewer
//tags: fileViewer, fileViewer-geojson, fileViewer-topojson, fileViewer-json
//input: file file
//output: view result
export async function gisGeoJSONFileViewer(file: DG.FileInfo): Promise<DG.View | null | DG.DataFrame> {
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
  viewFile.root.id = 'map-cont';
  viewFile.root.style.padding = '0px';

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-cont');
    useGeographic();
    ol.setViewOptions({
      projection: 'EPSG:4326',
    });
    if (isGeoTopo[1] === false) ol.addGeoJSONLayerFromStream(strBuf);
    else ol.addTopoJSONLayerFromStream(strBuf);
  }, 200);

  return viewFile;
}

//name: gisGeoJSONFileHandler
//tags: file-handler
//meta.ext: geojson, topojson, json
//input: string filecontent
//output: list tables
export function gisGeoJSONFileHandler(filecontent: string): DG.DataFrame[] {
  //detect the kind of json file
  const isGeoTopo = gisGeoJSONFileDetector(filecontent);
  if (isGeoTopo[0] === false) return [DG.DataFrame.fromJson(filecontent)];

  let dfFromJSON = null;
  const ol = new OpenLayers();
  // ol.initMap('map-container');
  // useGeographic();
  // ol.setViewOptions({
  //   projection: 'EPSG:4326',
  // });
  let newLayer;
  if (isGeoTopo[1] === false) newLayer = ol.addGeoJSONLayerFromStream(filecontent);
  else newLayer = ol.addTopoJSONLayerFromStream(filecontent);
  // const arrFeatures = ol.getFeaturesFromLayer(newLayer);
  const arrFeatures = ol.exportLayerToArray(newLayer);
  if (arrFeatures) {
    if (arrFeatures.length > 0) {
      dfFromJSON = DG.DataFrame.fromObjects(arrFeatures);
      if (dfFromJSON) {
        const gisCol = dfFromJSON.col('gisObject');
        if (gisCol)
          gisCol.semType = SEMTYPEGIS.GISAREA; //SEMTYPEGIS.GISOBJECT;

        dfFromJSON.name = newLayer.get('layerName');
        const tv = grok.shell.addTableView(dfFromJSON as DG.DataFrame);
        tv.name = 'dfFromJSON.name' + ' (manual)';

        // const mapViewer = tv.addViewer(DG.Viewer.fromType('Map', dfFromJSON)) as GisViewer;
        // if (mapViewer) {
        //   //
        // (mapViewer as GisViewer).init();
        //   (mapViewer as GisViewer).ol.addLayer(newLayer);
        // // if (isGeoTopo[1] === false) newLayer = (mapViewer as GisViewer).ol.addGeoJSONLayerFromStream(filecontent);
        // //   else newLayer = (mapViewer as GisViewer).ol.addTopoJSONLayerFromStream(filecontent);
        // }

        // const mapViewer = DG.Viewer.fromType('Map', (dfFromJSON as DG.DataFrame));
        //TODO: fix issue with GisViewer opening for a table with data
        setTimeout((tv, dfFromJSON, newLayer) => {
          // const v = grok.shell.v;
          const v = tv;
          // if ((v) && (v instanceof DG.TableView))
          //   (v as DG.TableView).addViewer(new GisViewer());
          if ((v) && (v instanceof DG.TableView)) {
            //
            const mapViewer = ((v as DG.TableView).addViewer('Map') as GisViewer);
            // const mapViewer = ((v as DG.TableView).addViewer(new GisViewer()) as GisViewer);
            if (newLayer)
              mapViewer.ol.addLayer(newLayer);
          // const mapViewer = ((v as DG.TableView).addViewer(DG.Viewer.fromType('Map', dfFromJSON)) as GisViewer);
            // if (isGeoTopo[1] === false) newLayer = mapViewer.ol.addGeoJSONLayerFromStream(filecontent);
            // else newLayer = mapViewer.ol.addTopoJSONLayerFromStream(filecontent);
          }
          // (v as DG.TableView).addViewer(DG.Viewer.fromType('Map', (v as DG.TableView).dataFrame));
        }, 4000, tv, dfFromJSON, newLayer);
      }
    }
  }
  if (dfFromJSON) return [dfFromJSON];
  return [DG.DataFrame.fromJson(filecontent)];
}

//openGISViewer
//TODO: if we want to show this in UI we need to annotate it properly and add export
function openGISViewer(): void {
  let htmlStyle: DG.ElementOptions = { };
  htmlStyle = {style: {'width': '100%', 'height': '100%'}};
  const boxMap = ui.box(null, htmlStyle);
  boxMap.id = 'map-container'; //boxMap - div that contains map
  grok.shell.newView('OpenLayers preview view', [boxMap]);

  const ol = new OpenLayers();
  ol.initMap('map-container');
}

//name: gisAreaWidget
//tags: panel, widgets
//input: string gisArea {semType: gis-area}
//output: widget result
//condition: true
export function gisAreaWidget(gisArea: any): DG.Widget | null {
  //this is temporary code - should be filled with usefull functionality
  // if ((!gisArea) || !(gisArea instanceof GisArea)) return null;
  if ((!gisArea)) return null;

  const strToAdd: string = 'test'; //(gisArea as GisArea).semtype;
  let widgetStyle: DG.ElementOptions = { };
  widgetStyle = {style: {'color': '#F55'}};

  return new DG.Widget(ui.divText('gis Area widget ' + strToAdd, widgetStyle));
}

// //name: gisAreaWidget
// //tags: panel, widgets
// //input: string area {semType: gis-area}
// //output: widget result
// //condition:  true
// export function gisAreaWidget(area: string): DG.Widget {
//   return new DG.Widget(ui.divText('gis area widget'));
// }

