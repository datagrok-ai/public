import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getCensusInfo
export async function getCensusInfo() : Promise<void> {
  await PackageFunctions.getCensusInfo();
}

//name: getCensusInfoTest UI
export async function getCensusInfoT() : Promise<void> {
  await PackageFunctions.getCensusInfoT();
}

//name: info
export async function info() : Promise<void> {
  await PackageFunctions.info();
}

//meta.role: init
export function init() : void {
  PackageFunctions.init();
}

//name: GISGeocoding
//description: GIS geocoding - receive coordinates for address
//input: string address 
//output: string result
export async function gisGeocoding(address: string) : Promise<string> {
  return await PackageFunctions.gisGeocoding(address);
}

//description: GIS geocoding - receive address for coordinates
//input: double x 
//input: double y 
//output: string result
export async function gisReverseGeocoding(x: number, y: number) : Promise<string> {
  return await PackageFunctions.gisReverseGeocoding(x, y);
}

//description: GIS geocoding - receive coordinates from address
//input: string address 
//output: string result
export async function gisBatchGeocoding(address: string) : Promise<string> {
  return await PackageFunctions.gisBatchGeocoding(address);
}

//name: Map
//description: GIS map viewer
//output: viewer result
//meta.icon: icons/icon.svg
//meta.toolbox: true
//meta.role: viewer
export function gisViewer() : any {
  return PackageFunctions.gisViewer();
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: kmz,kml
export async function gisKMZAndKMLFileViewer(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.gisKMZAndKMLFileViewer(file);
}

//input: string filecontent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: kml
export async function gisKMLFileHandler(filecontent: string) : Promise<any> {
  return await PackageFunctions.gisKMLFileHandler(filecontent);
}

//input: list filecontent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: kmz
export async function gisKMZFileHandler(filecontent: Uint8Array) : Promise<any> {
  return await PackageFunctions.gisKMZFileHandler(filecontent);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: geojson, topojson
export async function gisGeoJSONFileViewer(file: DG.FileInfo) {
  return await PackageFunctions.gisGeoJSONFileViewer(file);
}

//input: string filecontent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: geojson, topojson
export function gisGeoJSONFileHandler(filecontent: string) : any {
  return PackageFunctions.gisGeoJSONFileHandler(filecontent);
}

//input: string gisArea { semType: gis-area }
//output: widget result
//meta.role: widgets,panel
//condition: true
export function gisAreaWidget(gisArea: any) : any {
  return PackageFunctions.gisAreaWidget(gisArea);
}

//name: mapDemo
//description: Map viewer shows geospatial data on a map as either markers, or a heat map.
//meta.demoPath: Visualization | Geographical | Map
//test: _mapDemo() //wait: 2000 
export async function _mapDemo() : Promise<void> {
  await PackageFunctions._mapDemo();
}
