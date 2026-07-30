import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemspace
//meta.role: app
//meta.browsePath: Chem
export async function app() : Promise<void> {
  await PackageFunctions.app();
}

//name: Databases | Chemspace
//description: Finds similar and substructure matches for a molecule in the Chemspace catalog.
//input: string smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function samplesPanel(smiles: string) : Promise<any> {
  return await PackageFunctions.samplesPanel(smiles);
}

//name: Chemspace Prices
//description: Shows vendor prices and pack sizes for a Chemspace compound.
//input: string id { semType: chemspace-id; description: Chemspace compound id (e.g. CSCS00000000000) }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function pricesPanel(id: string) : Promise<any> {
  return await PackageFunctions.pricesPanel(id);
}

//description: Countries Chemspace can ship to, for the shipToCountry parameter
//output: list<string> result
export function getShipToCountries() : string[] {
  return PackageFunctions.getShipToCountries();
}

//name: Get Chemspace Ids
//description: Adds a column of Chemspace compound ids matched by exact structure for each molecule.
//input: column<string> molColumn { semType: Molecule; nullable: false }
//input: string shipToCountry { caption: Ship to; nullable: false; choices: Chemspace:getShipToCountries(); description: Destination country for pricing and availability }
//output: column result
//meta.vectorFunc: true
export async function getChemspaceIds(molColumn: DG.Column, shipToCountry: string) : Promise<any> {
  return await PackageFunctions.getChemspaceIds(molColumn, shipToCountry);
}

//name: Get Chemspace Prices
//description: Looks up vendor, pack size and price for each Chemspace id and joins the result back into the table.
//input: dataframe data { nullable: false }
//input: column<string> idsColumn { semType: chemspace-id; nullable: false; description: Column of Chemspace compound ids to price }
//input: string shipToCountry { caption: Ship to; nullable: false; choices: Chemspace:getShipToCountries(); description: Destination country for pricing and availability }
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
