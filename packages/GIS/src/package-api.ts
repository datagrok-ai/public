import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function getCensusInfoT(): Promise<any> {
    return await grok.functions.call('Gis:GetCensusInfoT', {});
  }

  export async function getCensusInfo(): Promise<any> {
    return await grok.functions.call('Gis:GetCensusInfo', {});
  }

  export async function info(): Promise<any> {
    return await grok.functions.call('Gis:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('Gis:Init', {});
  }

  //GIS geocoding - receive coordinates for address
  export async function gisGeocoding(address: string): Promise<any> {
    return await grok.functions.call('Gis:GisGeocoding', { address });
  }

  //GIS geocoding - receive address for coordinates
  export async function gisReverseGeocoding(x: number, y: number): Promise<any> {
    return await grok.functions.call('Gis:GisReverseGeocoding', { x, y });
  }

  //GIS geocoding - receive coordinates from address
  export async function gisBatchGeocoding(address: string): Promise<any> {
    return await grok.functions.call('Gis:GisBatchGeocoding', { address });
  }

  //GIS map viewer
  export async function gisViewer(): Promise<any> {
    return await grok.functions.call('Gis:GisViewer', {});
  }

  export async function gisKMZAndKMLFileViewer(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Gis:GisKMZAndKMLFileViewer', { file });
  }

  export async function gisKMLFileHandler(filecontent: string): Promise<any> {
    return await grok.functions.call('Gis:GisKMLFileHandler', { filecontent });
  }

  export async function gisKMZFileHandler(filecontent: any): Promise<any> {
    return await grok.functions.call('Gis:GisKMZFileHandler', { filecontent });
  }

  export async function gisGeoJSONFileViewer(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Gis:GisGeoJSONFileViewer', { file });
  }

  export async function gisGeoJSONFileHandler(filecontent: string): Promise<any> {
    return await grok.functions.call('Gis:GisGeoJSONFileHandler', { filecontent });
  }

  export async function gisAreaWidget(gisArea: string): Promise<any> {
    return await grok.functions.call('Gis:GisAreaWidget', { gisArea });
  }

  //Map viewer shows geospatial data on a map as either markers, or a heat map.
  export async function mapDemo(): Promise<any> {
    return await grok.functions.call('Gis:MapDemo', {});
  }
}
