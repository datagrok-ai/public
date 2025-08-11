import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemspace
//tags: app
//meta.browsePath: Chem
export async function app() {
  return PackageFunctions.app();
}

//name: Databases | Chemspace
//description: Chemspace Samples
//tags: panel, widgets
//input: string smiles { semType: Molecule }
//output: widget result
//condition: true
export async function samplesPanel(smiles: string) {
  return PackageFunctions.samplesPanel(smiles);
}

//name: createSearchPanel
//description: Creates search panel
//input: dynamic searchMode 
//input: string smiles 
//input: dynamic category 
//input: dynamic shipToCountry 
//output: dynamic result
export function createSearchPanel(searchMode: any, smiles: string, category: any, shipToCountry: any) {
  return PackageFunctions.createSearchPanel(searchMode, smiles, category, shipToCountry);
}

//name: Chemspace Prices
//description: Chemspace Prices
//tags: panel, widgets
//input: string id { semType: chemspace-id }
//output: widget result
//condition: true
export async function pricesPanel(id: string) {
  return PackageFunctions.pricesPanel(id);
}

//name: getChemspaceIds
//input: column<string> molColumn { semType: Molecule }
//input: string shipToCountry 
//output: column result
//meta.vectorFunc: true
export async function getChemspaceIds(molColumn: DG.Column, shipToCountry: string) {
  return PackageFunctions.getChemspaceIds(molColumn, shipToCountry);
}

//name: getChemspacePrices
//input: dataframe data 
//input: column<string> idsColumn { semType: chemspace-id }
//input: string shipToCountry 
//output: dataframe res { action: join(data) }
export async function getChemspacePrices(data: DG.DataFrame, idsColumn: DG.Column, shipToCountry: string) {
  return PackageFunctions.getChemspacePrices(data, idsColumn, shipToCountry);
}

//name: queryMultipart
//description: Perform query with multipart form data
//input: string path 
//input: string formParamsStr 
//input: string paramsStr { optional: true }
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function queryMultipart(path: string, formParamsStr: string, paramsStr?: string) {
  return PackageFunctions.queryMultipart(path, formParamsStr, paramsStr);
}
