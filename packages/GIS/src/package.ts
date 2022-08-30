/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';
//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer'; //TODO: remove it further (it needs just for preview)
import {useGeographic} from 'ol/proj';

//TODO: remove imports above when we'll be sure we don't need it anymore>>
// import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
// import TileLayer from 'ol/layer/Tile';
// import VectorLayer from 'ol/layer/Vector';
// import * as OLProj from 'ol/proj';
// import OSM from 'ol/source/OSM';
//ZIP utilities
import JSZip from 'jszip';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);
}

// name: GISGeocoding
// description: GIS geocoding - receive coortinates from address
// input: string address
// output: string result
// export function gisGeocoding(address: string): string {
//   // return new GisViewer();
// }

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
//tags: fileViewer, fileViewer-kmz
//input: file file
//output: view result
export async function gisKMZFileViewer(file: DG.FileInfo): Promise<DG.View> {
  const viewFile = DG.View.create();
  viewFile.name = 'Preview of: ' + file.name;
  viewFile.root.id = 'map-container'; //boxMap - div that contains map
  viewFile.root.style.padding = '0px';

  const strBuf = await file.readAsBytes();
  const kmlData = await getKMZData(strBuf);

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-container');
    ol.addKMLLayerFromStream(kmlData);
    ol.setViewOptions({projection: 'EPSG:3857'});
  }, 200); //should use it because we need visible and active View for map initializing

  return viewFile;
}

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
  viewFile.root.id = 'map-container'; //boxMap - div that contains map
  viewFile.root.style.padding = '0px';

  const strBuf = await file.readAsString();

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-container');
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

//name: gisGeoJSONFileViewer
//tags: fileViewer, fileViewer-geojson, fileViewer-topojson, fileViewer-json
//input: file file
//output: view result
export async function gisGeoJSONFileViewer(file: DG.FileInfo): Promise<DG.View | null> {
  //read file
  const strBuf = await file.readAsString();

  let arrTmp: any[] | null;
  //searching for patterns of GeoJSON data
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"](?:Multi)?Polygon/ig);
  let cntTypeGeo = arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/\'|\"type\'|\"\s?:\s?\'|\"(?:Multi)?Point/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/\'|\"type\'|\"\s?:\s?\'|\"(?:Multi)?LineString/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]Feature/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]GeometryCollection/ig);
  cntTypeGeo += arrTmp ? arrTmp.length : 0;
  // alert('type count = ' + cntTypeStr + ' ' + 'geo = ' + cntTypeGeo);
  //searching for patterns of TopoJSON data
  arrTmp = strBuf.match(/['"]type['"]\s?:\s?['"]Topology['"]/ig);
  let cntTypeTopo = arrTmp ? arrTmp.length : 0;
  // alert('geo = ' + cntTypeGeo);

  if (cntTypeGeo == 0) return null;

  const viewFile = DG.View.create();
  viewFile.name = 'Preview of: ' + file.name;
  viewFile.root.id = 'map-container'; //boxMap - div that contains map
  viewFile.root.style.padding = '0px';

  setTimeout(() => {
    const ol = new OpenLayers();
    ol.initMap('map-container');
    useGeographic();
    ol.setViewOptions({
      projection: 'EPSG:4326',
    });
    if (cntTypeTopo == 0) ol.addGeoJSONLayerFromStream(strBuf);
    else ol.addTopoJSONLayerFromStream(strBuf);
  }, 200);

  return viewFile;
}

//name: gisGeoJSONFileHandler
//tags: file-handler
//meta.ext: geojson
//input: file file
//output: list tables
export function gisGeoJSONFileHandler(filecontent: string): DG.DataFrame[] {
  // ... processing files ...
  const ol = new OpenLayers();
  // ol.initMap('map-container');
  // useGeographic();
  // ol.setViewOptions({
  //   projection: 'EPSG:4326',
  // });
  // if (cntTypeTopo == 0) ol.addGeoJSONLayerFromStream(filecontent);
  // else ol.addTopoJSONLayerFromStream(filecontent);
  const newLayer = ol.addTopoJSONLayerFromStream(filecontent, false);
  ol.getFeaturesFromLayer(newLayer);

  // return [new DG.DataFrame(null)];
  return [DG.DataFrame.fromJson(filecontent)];
}


//name: openGISViewer
export function openGISViewer(): void {
  let htmlStyle: DG.ElementOptions = { };
  htmlStyle = {style: {'width': '100%', 'height': '100%'}};
  const boxMap = ui.box(null, htmlStyle);
  boxMap.id = 'map-container'; //boxMap - div that contains map
  grok.shell.newView('OpenLayers preview view', [boxMap]);

  const ol = new OpenLayers();
  ol.initMap('map-container');
}
