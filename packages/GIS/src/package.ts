/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';
//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer'; //TODO: remobve it further (it needs just for preview)

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);
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
//output: view v
export function gisKMLFileViewer(file: DG.FileInfo): DG.View {
  const viewFile = DG.View.create();
  let htmlStyle: DG.ElementOptions = {style: {'width': '100%', 'height': '100%', 'border': 'solid 1px yellow'}};
  // const boxMap = ui.box(null, htmlStyle);
  const boxMap = ui.box(ui.div(file.fullPath), htmlStyle);
  boxMap.id = 'map-container'; //boxMap - div that contains map
  viewFile.append(boxMap);

  // let strBuf = await file.readAsString();
  // const boxMap2 = ui.box(ui.div(strBuf), htmlStyle);

  // viewFile.append(boxMap2);
  // const ol = new OpenLayers();
  // ol.initMap('map-container');

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
