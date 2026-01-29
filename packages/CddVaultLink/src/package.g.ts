import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: CDD Vault
//input: string path { meta.url: true; optional: true }
//input: string filter { optional: true }
//output: view result
//meta.role: app
//meta.icon: images/cdd-icon-small.png
//meta.browsePath: Chem
export async function cddVaultApp(path: string, filter: string) : Promise<any> {
  return await PackageFunctions.cddVaultApp(path, filter);
}

//input: dynamic treeNode 
export async function cddVaultAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.cddVaultAppTreeBrowser(treeNode);
}

//name: Databases | CDD Vault
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: panel
export function molColumnPropertyPanel(molecule: string) : any {
  return PackageFunctions.molColumnPropertyPanel(molecule);
}

//input: funccall call 
//meta.role: editor
export async function CDDVaultSearchEditor(call: DG.FuncCall) : Promise<void> {
  await PackageFunctions.CDDVaultSearchEditor(call);
}

//name: Get Vault Stats
//input: int vaultId 
//input: string vaultName 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getVaultStats(vaultId: number, vaultName: string) : Promise<string> {
  return await PackageFunctions.getVaultStats(vaultId, vaultName);
}

//name: Get Vaults
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getVaults() : Promise<string> {
  return await PackageFunctions.getVaults();
}

//name: Get Molecules
//input: int vaultId { nullable: true }
//input: string moleculesIds 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getMolecules(vaultId: number, moleculesIds: string) : Promise<any> {
  return await PackageFunctions.getMolecules(vaultId, moleculesIds);
}

//name: Get Molecules Async
//input: int vaultId { nullable: true }
//input: string moleculesIds 
//input: int timeoutMinutes 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getMoleculesAsync(vaultId: number, moleculesIds: string, timeoutMinutes: number) : Promise<any> {
  return await PackageFunctions.getMoleculesAsync(vaultId, moleculesIds, timeoutMinutes);
}

//name: Get Protocols Async
//input: int vaultId { nullable: true }
//input: int timeoutMinutes 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getProtocolsAsync(vaultId: number, timeoutMinutes: number) : Promise<string> {
  return await PackageFunctions.getProtocolsAsync(vaultId, timeoutMinutes);
}

//name: Get Collections Async
//input: int vaultId { nullable: true }
//input: int timeoutMinutes 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getCollectionsAsync(vaultId: number, timeoutMinutes: number) : Promise<string> {
  return await PackageFunctions.getCollectionsAsync(vaultId, timeoutMinutes);
}

//name: Get Saved Searches
//input: int vaultId { nullable: true }
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getSavedSearches(vaultId: number) : Promise<string> {
  return await PackageFunctions.getSavedSearches(vaultId);
}

//name: Get Saved Search Results
//input: int vaultId 
//input: int searchId 
//input: int timeoutMinutes 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getSavedSearchResults(vaultId: number, searchId: number, timeoutMinutes: number) : Promise<any> {
  return await PackageFunctions.getSavedSearchResults(vaultId, searchId, timeoutMinutes);
}

//name: CDD Vault Search Async
//input: int vaultId { nullable: true }
//input: string structure { category: Structure; nullable: true; semType: Molecule; description: SMILES; cxsmiles or mol string }
//input: string structure_search_type { category: Structure; nullable: true; choices: ["exact","similarity","substructure"]; description: SMILES; cxsmiles or mol string }
//input: double structure_similarity_threshold { category: Structure; nullable: true; description: A number between 0 and 1 }
//input: int protocol { category: Protocol; nullable: true; description: Protocol id }
//input: int run { category: Protocol; nullable: true; description: Specific run id }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//editor: Cddvaultlink:CDDVaultSearchEditor
export async function cDDVaultSearchAsync(vaultId: number, structure?: string, structure_search_type?: any, structure_similarity_threshold?: number, protocol?: number, run?: number) : Promise<any> {
  return await PackageFunctions.cDDVaultSearchAsync(vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run);
}

//name: CDD Vault search
//input: int vaultId { nullable: true }
//input: string molecules { category: General; nullable: true; description: Comma separated list of ids }
//input: string names { category: General; nullable: true; description: Comma separated list of names/synonyms }
//input: bool include_original_structures { category: General; nullable: true; description: If true,include the original user defined structure for each molecule }
//input: bool only_ids { category: General; nullable: true; description: If true,only the Molecule IDs are returned,allowing for a smaller and faster response }
//input: bool only_batch_ids { category: General; nullable: true; description: If true,the full Molecule details are still returned but the Batch-level information is left out of the JSON results. (Only the IDs of the Batches belonging to the Molecules are still included.) }
//input: string created_before { category: General; nullable: true; description: ISO 8601 date }
//input: string created_after { category: General; nullable: true; description: ISO 8601 date }
//input: string modified_before { category: General; nullable: true; description: ISO 8601 date }
//input: string modified_after { category: General; nullable: true; description: ISO 8601 date }
//input: string batch_created_before { category: Batch fields; nullable: true; description: ISO 8601 date. A molecule with any batch that has a creation date on or before the parameter will be included }
//input: string batch_created_after { category: Batch fields; nullable: true; description: ISO 8601 date. A molecule with any batch that has a creation date on or after the parameter will be included }
//input: string batch_field_before_name { category: Batch fields; nullable: true; description: Specifes a user-defined batch field for batch_field_before_date }
//input: string batch_field_before_date { category: Batch fields; nullable: true; description: ISO 8601 date. A molecule with any batch that has a batch_field_before_name value date on or before the parameter will be included }
//input: string batch_field_after_name { category: Batch fields; nullable: true; description: Specifes a user-defined batch field for batch_field_after_date }
//input: string batch_field_after_date { category: Batch fields; nullable: true; description: ISO 8601 date. A molecule with any batch that has a batch_field_after_name value date on or after the parameter will be included }
//input: string projects { category: Projects; nullable: true; description: Comma separated list of project ids }
//input: string data_sets { category: Datasets; nullable: true; description: Comma separated list of dataset ids }
//input: string structure { category: Structure; nullable: true; semType: Molecule; description: SMILES,cxsmiles or mol string }
//input: string structure_search_type { category: Structure; nullable: true; choices: ["exact","similarity","substructure"]; description: SMILES,cxsmiles or mol string }
//input: double structure_similarity_threshold { category: Structure; nullable: true; description: A number between 0 and 1 }
//input: string inchikey { category: Structure; nullable: true; description: Use this parameter instead of the 'structure' and 'structure_search_type' parameters }
//input: list<string> molecule_fields { category: Filelds; nullable: true; description: Use this parameter to limit the number of Molecule UDF Fields to return }
//input: list<string> batch_fields { category: Filelds; nullable: true; description: Use this parameter to limit the number of Batch UDF Fields to return }
//input: list<string> fields_search { category: Molecules filelds search; nullable: true; description: This parameter is used for searching across the custom user-defined Molecule fields created by your Vault Administrator }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function cDDVaultSearch(vaultId: number, molecules: string, names: string, include_original_structures: boolean, only_ids: boolean, only_batch_ids: boolean, created_before: string, created_after: string, modified_before: string, modified_after: string, batch_created_before: string, batch_created_after: string, batch_field_before_name: string, batch_field_before_date: string, batch_field_after_name: string, batch_field_after_date: string, projects: string, data_sets: string, structure: string, structure_search_type: any, structure_similarity_threshold: number, inchikey: string, molecule_fields: string[], batch_fields: string[], fields_search: string[]) : Promise<any> {
  return await PackageFunctions.cDDVaultSearch(vaultId, molecules, names, include_original_structures, only_ids, only_batch_ids, created_before, created_after, modified_before, modified_after, batch_created_before, batch_created_after, batch_field_before_name, batch_field_before_date, batch_field_after_name, batch_field_after_date, projects, data_sets, structure, structure_search_type, structure_similarity_threshold, inchikey, molecule_fields, batch_fields, fields_search);
}

//name: CDD Vault search 2
//input: int vaultId { nullable: true }
//input: string structure { category: Structure; nullable: true; semType: Molecule; description: SMILES,cxsmiles or mol string }
//input: string structure_search_type { category: Structure; nullable: true; choices: ["exact","similarity","substructure"]; description: SMILES,cxsmiles or mol string }
//input: double structure_similarity_threshold { category: Structure; nullable: true; description: A number between 0 and 1 }
//input: int protocol { category: Protocol; nullable: true; description: Protocol id }
//input: int run { category: Protocol; nullable: true; description: Specific run id }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//editor: Cddvaultlink:CDDVaultSearchEditor
export async function cDDVaultSearch2(vaultId: number, structure?: string, structure_search_type?: any, structure_similarity_threshold?: number, protocol?: number, run?: number) : Promise<any> {
  return await PackageFunctions.cDDVaultSearch2(vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run);
}
