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

const census = require('citysdk');

//ZIP utilities
import JSZip from 'jszip';
import { DataFrame } from 'datagrok-api/dg';

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

//name: getCensusInfo
export async function getCensusInfo() {
  let censusRes: any = null;
  let url = 'https://api.census.gov/data/';
  censusRes = await (await grok.dapi.fetchProxy(url)).json();

  if (!censusRes) return 'Fetch error';
  // alert(JSON.stringify(censusRes));
  const df = DataFrame.fromJson(JSON.stringify(censusRes));
  grok.shell.addTableView(df);
}

//name: info
export async function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);

  //TODO: remove this temporary code (was added for testing census feature)
  let censusRes: any = null;
  try {
    censusRes = await censusPromise(
      {
        statsKey: '2647d704d8734665d5c417dae1546887c2c90513', //TODO: hide credentials
        vintage: '2017',
        geoHierarchy: {
          county: {
            lat: 28.2639,
            lng: -80.7214,
          },
        },
      },
    );
  } catch (err: any) {
    grok.shell.error('Census error: ' + err);
  }
  if (censusRes)
    alert(censusRes);
  //<<remove all above
}

//tags: init, autostart
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
  // const df = DG.DataFrame.fromJson(fetchresult);
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
  // const df = DG.DataFrame.fromJson(fetchresult);
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

//name: GISViewer
//description: GIS map viewer
//tags: viewer
//output: viewer result
export function gisViewer(): GisViewer {
  setTimeout(() => {grok.shell.windows.showProperties = true;}, 200);

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

/*
//name: gisKMLFileViewer
//tags: fileViewer, fileViewer-kml
//input: file file
//output: view result
export async function gisKMLFileViewer(file: DG.FileInfo): Promise<DG.View> {
  const viewFile = DG.View.create();
  //alternative way - with embedding additional div into view window>>
  // const boxMap = ui.box(null);
  // boxMap.id = 'map-container'; //boxMap - div that contains map
  // viewFile.append(boxMap);
  viewFile.name = 'Preview of: ' + file.name;
  viewFile.root.id = 'map-cont';
  viewFile.root.style.padding = '0px';

  const strBuf = await file.readAsString();

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-cont');
    useGeographic();
    ol.setViewOptions({
      projection: 'EPSG:4326',
      // center: OLProj.fromLonLat([34.109565, 45.452962]),
      // zoom: 7
    });
    ol.addKMLLayerFromStream(strBuf);
  }, 200);

  return viewFile;
}
*/

//gisGeoJSONFileDetector
function gisGeoJSONFileDetector(strBuf: string): [boolean, boolean] {
  let arrTmp: any[] | null;
  //searching for patterns of GeoJSON data
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"](?:Multi)?Polygon/ig);
  let cntTypeGeo = arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?\'|\"(?:Multi)?Point/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?\'|\"(?:Multi)?LineString/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]Feature/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]GeometryCollection/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  //searching for patterns of TopoJSON data
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]Topology['"]/ig);
  const cntTypeTopo = arrTmp ? arrTmp.length : 0;

  if (cntTypeGeo == 0) return [false, false];
  if (cntTypeTopo > 0) return [true, true];

  return [true, false];
}

//name: gisGeoJSONFileViewer
//tags: fileViewer, fileViewer-geojson, fileViewer-topojson, fileViewer-json
//input: file file
//output: view result
export async function gisGeoJSONFileViewer(file: DG.FileInfo): Promise<DG.View | null> {
  //read file
  const strBuf = await file.readAsString();
  const isGeoTopo = gisGeoJSONFileDetector(strBuf);

  if (isGeoTopo[0] === false) {
    //if json file is not kind of geoJson or topoJson - show it as a table
    const df = DG.DataFrame.fromJson(strBuf);
    const viewFile = DG.TableView.create(df);
    return viewFile;
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

  let dfFormJSON = null;
  const ol = new OpenLayers();
  // ol.initMap('map-container');
  // useGeographic();
  // ol.setViewOptions({
  //   projection: 'EPSG:4326',
  // });
  let newLayer;
  if (isGeoTopo[1] === false) newLayer = ol.addGeoJSONLayerFromStream(filecontent);
  else newLayer = ol.addTopoJSONLayerFromStream(filecontent);
  const arrFeatures = ol.getFeaturesFromLayer(newLayer);
  if (arrFeatures) {
    if (arrFeatures.length > 0) {
      const dfFormJSON = DG.DataFrame.fromObjects(arrFeatures);
      if (dfFormJSON) {
        dfFormJSON.name = newLayer.get('layerName');
        const tv = grok.shell.addTableView(dfFormJSON as DG.DataFrame);
        tv.name = 'dfFormJSON.name' + ' (manual)';
        //TODO: fix issue with GisViewer opening for a table with data
        setTimeout((tv) => {
          // const v = grok.shell.v;
          const v = tv;
          if ((v) && (v instanceof DG.TableView)) (v as DG.TableView).addViewer(new GisViewer());
        }, 500);
      }
    }
  }
  if (dfFormJSON) return [dfFormJSON];
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
