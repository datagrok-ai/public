import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getCensusInfo
export async function getCensusInfo() {
  return PackageFunctions.getCensusInfo();
}

//name: getCensusInfoTest UI
export async function getCensusInfoT() {
  return PackageFunctions.getCensusInfoT();
}

//name: info
export async function info() {
  return PackageFunctions.info();
}

//tags: init
export function init() {
  return PackageFunctions.init();
}

//name: GISGeocoding
//description: GIS geocoding - receive coordinates for address
//input: string address 
//output: string result
export async function gisGeocoding(address: string) {
  return PackageFunctions.gisGeocoding(address);
}

//name: gisReverseGeocoding
//description: GIS geocoding - receive address for coordinates
//input: double x 
//input: double y 
//output: string result
export async function gisReverseGeocoding(x: number, y: number) {
  return PackageFunctions.gisReverseGeocoding(x, y);
}

//name: gisBatchGeocoding
//description: GIS geocoding - receive coordinates from address
//input: string address 
//output: string result
export async function gisBatchGeocoding(address: string) {
  return PackageFunctions.gisBatchGeocoding(address);
}

//name: Map
//description: GIS map viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/icon.svg
//meta.toolbox: true
export function gisViewer() {
  return PackageFunctions.gisViewer();
}

//name: gisKMZAndKMLFileViewer
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: kmz,kml
export async function gisKMZAndKMLFileViewer(file: DG.FileInfo) {
  return PackageFunctions.gisKMZAndKMLFileViewer(file);
}

//name: gisKMLFileHandler
//tags: file-handler
//input: string filecontent 
//output: list<dataframe> result
//meta.ext: kml
export async function gisKMLFileHandler(filecontent: string) {
  return PackageFunctions.gisKMLFileHandler(filecontent);
}

//name: gisKMZFileHandler
//tags: file-handler
//input: list filecontent 
//output: list<dataframe> result
//meta.ext: kmz
export async function gisKMZFileHandler(filecontent: Uint8Array) {
  return PackageFunctions.gisKMZFileHandler(filecontent);
}

//name: gisGeoJSONFileViewer
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: geojson, topojson
export async function gisGeoJSONFileViewer(file: DG.FileInfo) {
  return PackageFunctions.gisGeoJSONFileViewer(file);
}

//name: gisGeoJSONFileHandler
//tags: file-handler
//input: string filecontent 
//output: list<dataframe> result
//meta.ext: geojson, topojson
export function gisGeoJSONFileHandler(filecontent: string) {
  return PackageFunctions.gisGeoJSONFileHandler(filecontent);
}

//name: gisAreaWidget
//tags: panel, widgets
//input: string gisArea { semType: gis-area }
//output: widget result
//condition: true
export function gisAreaWidget(gisArea: any) {
  return PackageFunctions.gisAreaWidget(gisArea);
}

//name: mapDemo
//description: Map viewer shows geospatial data on a map as either markers, or a heat map.
//meta.demoPath: Visualization | Geographical | Map
//test: _mapDemo() //wait: 2000 , skip: GROK-11670 
export async function _mapDemo() {
  return PackageFunctions._mapDemo();
}
