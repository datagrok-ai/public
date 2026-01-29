import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemspace
//meta.role: app
//meta.browsePath: Chem
export async function app() : Promise<void> {
  await PackageFunctions.app();
}

//name: Databases | Chemspace
//description: Chemspace Samples
//input: string smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function samplesPanel(smiles: string) : Promise<any> {
  return await PackageFunctions.samplesPanel(smiles);
}

//name: Chemspace Prices
//description: Chemspace Prices
//input: string id { semType: chemspace-id }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function pricesPanel(id: string) : Promise<any> {
  return await PackageFunctions.pricesPanel(id);
}

//input: column<string> molColumn { semType: Molecule }
//input: string shipToCountry 
//output: column result
//meta.vectorFunc: true
export async function getChemspaceIds(molColumn: DG.Column, shipToCountry: string) : Promise<any> {
  return await PackageFunctions.getChemspaceIds(molColumn, shipToCountry);
}

//input: dataframe data 
//input: column<string> idsColumn { semType: chemspace-id }
//input: string shipToCountry 
//output: dataframe res { action: join(data) }
export async function getChemspacePrices(data: DG.DataFrame, idsColumn: DG.Column, shipToCountry: string) : Promise<any> {
  return await PackageFunctions.getChemspacePrices(data, idsColumn, shipToCountry);
}

//description: Perform query with multipart form data
//input: string path 
//input: string formParamsStr 
//input: string paramsStr { optional: true }
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function queryMultipart(path: string, formParamsStr: string, paramsStr?: string) : Promise<string> {
  return await PackageFunctions.queryMultipart(path, formParamsStr, paramsStr);
}
