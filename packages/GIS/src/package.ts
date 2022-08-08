/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';
//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer'; //TODO: remove it further (it needs just for preview)
import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
import TileLayer from 'ol/layer/Tile';
import VectorLayer from 'ol/layer/Vector';
import * as OLProj from 'ol/proj';
import OSM from 'ol/source/OSM';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);
}

//name: GISGeocoding
//description: GIS geocoding - receive coortinates from address
//input: string address
//output: string result
export function gisGeocoding(address: string): string {
  // return new GisViewer();
}

//name: GISViewer
//description: GIS map viewer
//tags: viewer
//output: viewer result
export function gisViewer(): GisViewer {
  return new GisViewer();
}

//name: gisKMLFileViewer
//tags: fileViewer, fileViewer-kml
//input: file file
//output: view result
export async function gisKMLFileViewer(file: DG.FileInfo): Promise<DG.View> {
// export function gisKMLFileViewer(file: DG.FileInfo): DG.View {
  let htmlStyle: DG.ElementOptions = {style: {'width': '100%', 'height': '100%', 'border': 'solid 1px yellow'}};
  const viewFile = DG.View.create();
  // const boxMap = ui.box(null, htmlStyle);
  //const divMap = ui.div(file.fullPath);
  const divMap = ui.div([], htmlStyle);
  divMap.id = 'map-cnt'; //boxMap - div that contains map
  viewFile.append(divMap);

  // const boxMap = ui.box(divMap, htmlStyle);
  // boxMap.id = 'map-cnt'; //boxMap - div that contains map
  // viewFile.append(boxMap);

  let strBuf = await file.readAsString();
  // const boxMap2 = ui.box(ui.div(strBuf), htmlStyle);
  // viewFile.append(boxMap2);

  setTimeout(() => {
    let ol = new OpenLayers();
    ol.initMap('map-cnt');
    ol.addKMLLayerFromStream(strBuf);
  }, 200);

  // setTimeout(() => {
  //   const newLayer = new TileLayer({
  //     visible: true,
  //     preload: Infinity,
  //     source: new OSM()});

  //   olMap = new OLMap({
  //     target: 'map-cnt',
  //     layers: [newLayer],
  //     // controls: defaultControls({attribution: false, rotate: false}),
  //     view: new OLView({
  //       center: OLProj.fromLonLat([34.109565, 45.452962]),
  //       zoom: 7,
  //     }),
  //   });
  // }, 2000);

  // openGISViewer();

  return viewFile;
}

//name: openGISViewer
export function openGISViewer(): void {
  let htmlStyle: DG.ElementOptions = { };
  htmlStyle = {style: {'width': '100%', 'height': '100%', 'border': 'solid 1px yellow'}};
  const boxMap = ui.box(null, htmlStyle);
  boxMap.id = 'map-container'; //boxMap - div that contains map
  grok.shell.newView('Leaflet preview view', [boxMap]);

  const ol = new OpenLayers();
  ol.initMap('map-container');
}
