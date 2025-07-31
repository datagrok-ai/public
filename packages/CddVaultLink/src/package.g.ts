import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: CDD Vault
//tags: app
//input: string path { meta.url: true; optional: true }
//input: string filter { optional: true }
//output: view result
//meta.icon: images/cdd-icon-small.png
//meta.browsePath: Chem
export async function cddVaultApp(path: string, filter: string) {
  return PackageFunctions.cddVaultApp(path, filter);
}

//name: cddVaultAppTreeBrowser
//input: dynamic treeNode 
//output: dynamic result
export async function cddVaultAppTreeBrowser(treeNode: any) {
  return PackageFunctions.cddVaultAppTreeBrowser(treeNode);
}

//name: Databases | CDD Vault
//tags: panel
//input: string molecule { semType: Molecule }
//output: widget result
export function molColumnPropertyPanel(molecule: string) {
  return PackageFunctions.molColumnPropertyPanel(molecule);
}

//name: CDDVaultSearchEditor
//tags: editor
//input: funccall call 
export async function CDDVaultSearchEditor(call: DG.FuncCall) {
  return PackageFunctions.CDDVaultSearchEditor(call);
}

//name: Get Vault Stats
//input: int vaultId 
//input: string vaultName 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getVaultStats(vaultId: number, vaultName: string) {
  return PackageFunctions.getVaultStats(vaultId, vaultName);
}

//name: Get Vaults
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getVaults() {
  return PackageFunctions.getVaults();
}

//name: Get Molecules
//input: int vaultId { nullable: true }
//input: string moleculesIds 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getMolecules(vaultId: number, moleculesIds: string) {
  return PackageFunctions.getMolecules(vaultId, moleculesIds);
}

//name: Get Molecules Async
//input: int vaultId { nullable: true }
//input: string moleculesIds 
//input: int timeoutMinutes 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getMoleculesAsync(vaultId: number, moleculesIds: string, timeoutMinutes: number) {
  return PackageFunctions.getMoleculesAsync(vaultId, moleculesIds, timeoutMinutes);
}

//name: Get Protocols Async
//input: int vaultId { nullable: true }
//input: int timeoutMinutes 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getProtocolsAsync(vaultId: number, timeoutMinutes: number) {
  return PackageFunctions.getProtocolsAsync(vaultId, timeoutMinutes);
}

//name: Get Collections Async
//input: int vaultId { nullable: true }
//input: int timeoutMinutes 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getCollectionsAsync(vaultId: number, timeoutMinutes: number) {
  return PackageFunctions.getCollectionsAsync(vaultId, timeoutMinutes);
}

//name: Get Saved Searches
//input: int vaultId { nullable: true }
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getSavedSearches(vaultId: number) {
  return PackageFunctions.getSavedSearches(vaultId);
}

//name: Get Saved Search Results
//input: int vaultId 
//input: int searchId 
//input: int timeoutMinutes 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function getSavedSearchResults(vaultId: number, searchId: number, timeoutMinutes: number) {
  return PackageFunctions.getSavedSearchResults(vaultId, searchId, timeoutMinutes);
}

//name: CDD Vault Search Async
//input: int vaultId { nullable: true }
//input: string structure { category: Structure; nullable: true; semType: Molecule }
//input: string structure_search_type { category: Structure; nullable: true; choices: ['exact','similarity','substructure'] }
//input: double structure_similarity_threshold { category: Structure; nullable: true }
//input: int protocol { category: Protocol; nullable: true }
//input: int run { category: Protocol; nullable: true }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//editor: Cddvaultlink
export async function cDDVaultSearchAsync(vaultId: number, structure: string, structure_search_type: any, structure_similarity_threshold: number, protocol: number, run: number) {
  return PackageFunctions.cDDVaultSearchAsync(vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run);
}

//name: CDD Vault search
//input: int vaultId { nullable: true }
//input: string molecules { category: General; nullable: true }
//input: string names { category: General; nullable: true }
//input: bool include_original_structures { category: General; nullable: true }
//input: bool only_ids { category: General; nullable: true }
//input: bool only_batch_ids { category: General; nullable: true }
//input: string created_before { category: General; nullable: true }
//input: string created_after { category: General; nullable: true }
//input: string modified_before { category: General; nullable: true }
//input: string modified_after { category: General; nullable: true }
//input: string batch_created_before { category: Batch fields; nullable: true }
//input: string batch_created_after { category: Batch fields; nullable: true }
//input: string batch_field_before_name { category: Batch fields; nullable: true }
//input: string batch_field_before_date { category: Batch fields; nullable: true }
//input: string batch_field_after_name { category: Batch fields; nullable: true }
//input: string batch_field_after_date { category: Batch fields; nullable: true }
//input: string projects { category: Projects; nullable: true }
//input: string data_sets { category: Datasets; nullable: true }
//input: string structure { category: Structure; nullable: true; semType: Molecule }
//input: string structure_search_type { category: Structure; nullable: true; choices: ['exact','similarity','substructure'] }
//input: double structure_similarity_threshold { category: Structure; nullable: true }
//input: string inchikey { category: Structure; nullable: true }
//input: list<string> molecule_fields { category: Filelds; nullable: true }
//input: list<string> batch_fields { category: Filelds; nullable: true }
//input: list<string> fields_search { category: Molecules filelds search; nullable: true }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function cDDVaultSearch(vaultId: number, molecules: string, names: string, include_original_structures: boolean, only_ids: boolean, only_batch_ids: boolean, created_before: string, created_after: string, modified_before: string, modified_after: string, batch_created_before: string, batch_created_after: string, batch_field_before_name: string, batch_field_before_date: string, batch_field_after_name: string, batch_field_after_date: string, projects: string, data_sets: string, structure: string, structure_search_type: any, structure_similarity_threshold: number, inchikey: string, molecule_fields: string[], batch_fields: string[], fields_search: string[]) {
  return PackageFunctions.cDDVaultSearch(vaultId, molecules, names, include_original_structures, only_ids, only_batch_ids, created_before, created_after, modified_before, modified_after, batch_created_before, batch_created_after, batch_field_before_name, batch_field_before_date, batch_field_after_name, batch_field_after_date, projects, data_sets, structure, structure_search_type, structure_similarity_threshold, inchikey, molecule_fields, batch_fields, fields_search);
}

//name: CDD Vault search 2
//input: int vaultId { nullable: true }
//input: string structure { category: Structure; nullable: true; semType: Molecule }
//input: string structure_search_type { category: Structure; nullable: true; choices: ['exact','similarity','substructure'] }
//input: double structure_similarity_threshold { category: Structure; nullable: true }
//input: int protocol { category: Protocol; nullable: true }
//input: int run { category: Protocol; nullable: true }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//editor: Cddvaultlink
export async function cDDVaultSearch2(vaultId: number, structure: string, structure_search_type: any, structure_similarity_threshold: number, protocol: number, run: number) {
  return PackageFunctions.cDDVaultSearch2(vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run);
}
