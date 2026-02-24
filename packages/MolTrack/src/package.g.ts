import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: initDB
export async function initDB() : Promise<void> {
  await PackageFunctions.initDB();
}

//name: MolTrack
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.browsePath: Chem
//meta.icon: images/moltrack.png
//meta.role: app
export async function molTrackApp(path: string) : Promise<any> {
  return await PackageFunctions.molTrackApp(path);
}

//input: dynamic appNode 
//input: view browseView 
//meta.role: appTreeBrowser
//meta.app: MolTrack
export async function molTrackAppTreeBrowser(appNode: any, browseView: any) : Promise<void> {
  await PackageFunctions.molTrackAppTreeBrowser(appNode, browseView);
}

//description: Retrieves all properties defined for the 'compound' scope
//output: string result
//cache: all
//cache.invalidateOn: 0 0 1 * *
export async function fetchCompoundProperties() : Promise<string> {
  return await PackageFunctions.fetchCompoundProperties();
}

//description: Retrieves all properties defined for the 'batch' scope
//output: string result
//cache: all
//cache.invalidateOn: 0 0 1 * *
export async function fetchBatchProperties() : Promise<string> {
  return await PackageFunctions.fetchBatchProperties();
}

//description: Retrieves all dynamic fields
//output: string result
export async function fetchSchema() : Promise<string> {
  return await PackageFunctions.fetchSchema();
}

//description: Retrieves all static fields
//output: string result
export async function fetchDirectSchema() : Promise<string> {
  return await PackageFunctions.fetchDirectSchema();
}

//input: string corporateValue 
//output: object compound
export async function getCompoundByCorporateId(corporateValue: string) {
  return await PackageFunctions.getCompoundByCorporateId(corporateValue);
}

//input: string corporateValue 
//output: object batch
export async function getBatchByCorporateId(corporateValue: string) {
  return await PackageFunctions.getBatchByCorporateId(corporateValue);
}

//description: Registers compound properties in the MolTrack service using the provided JSON payload
//input: string jsonPayload 
//output: string result
export async function registerMolTrackProperties(jsonPayload: string) : Promise<string> {
  return await PackageFunctions.registerMolTrackProperties(jsonPayload);
}

//input: string assayPayload 
//output: string result
export async function registerAssays(assayPayload: string) : Promise<string> {
  return await PackageFunctions.registerAssays(assayPayload);
}

//input: file csvFile 
//input: string scope 
//input: string mapping 
//input: string errorHandling 
//output: dataframe result
export async function registerBulk(csvFile: DG.FileInfo, scope: string, mapping: string, errorHandling: string) : Promise<any> {
  return await PackageFunctions.registerBulk(csvFile, scope, mapping, errorHandling);
}

//input: string query 
//input: string entityEndpoint 
//output: dataframe df
export async function search(query: string, entityEndpoint: string) {
  return await PackageFunctions.search(query, entityEndpoint);
}

//description: Performs a structured MolTrack compound search. The caller must provide the "output": a list of fully-qualified field names to return; "filter": a structured filter tree describing search conditions.
//input: list outputFields { description: List of fields to return. All fields must use MolTrack notation. Valid formats: "<table>.<field>" (direct database column) or "<table>.details.<property>" (dynamic detail property). Valid direct compound fields include (non-exhaustive): "compounds.created_at","compounds.updated_at","compounds.created_by","compounds.updated_by","compounds.canonical_smiles","compounds.original_molfile","compounds.molregno","compounds.inchi","compounds.inchikey","compounds.formula","compounds.hash_mol","compounds.hash_tautomer","compounds.hash_canonical_smiles","compounds.hash_no_stereo_smiles","compounds.hash_no_stereo_tautomer","compounds.is_archived". }
//input: object filter { description: A structured boolean filter supporting simple and nested conditions. Simple format: {"field":"<field_name>","operator":"<operator>","value":<value>,"threshold":<number|null>} where threshold is required only for "IS SIMILAR". String operators: "=","!=","IN","STARTS WITH","ENDS WITH","LIKE","CONTAINS". Numeric operators: "<",">","<=",">=","RANGE" (expects {value:[min,max]}). Datetime operators: "BEFORE","AFTER","ON" (ISO 8601). Molecular operators (only for compounds.structure): "IS SIMILAR" (requires SMILES + numeric similarity threshold), "IS SUBSTRUCTURE OF", "HAS SUBSTRUCTURE". Complex nested conditions: {"operator":"AND"|"OR","conditions":[...]}. Notes: nested conditions may be arbitrarily deep, each branch must be simple or complex. Examples: simple {"field":"compounds.formula","operator":"=","value":"C6H6"}, complex {"operator":"AND","conditions":[{"field":"compounds.created_at","operator":"AFTER","value":"2025-01-01"},{"operator":"OR","conditions":[{"field":"compounds.formula","operator":"=","value":"C6H6"},{"field":"compounds.structure","operator":"IS SIMILAR","value":"c1ccccc1","threshold":0.8}]}]}. }
//output: dataframe searchResult
export async function advancedSearch(outputFields: string[], filter: any) : Promise<any> {
  return await PackageFunctions.advancedSearch(outputFields, filter);
}

//input: string scope 
//output: dataframe result
export async function retrieveEntity(scope: string) : Promise<any> {
  return await PackageFunctions.retrieveEntity(scope);
}

//name: Databases | MolTrack
//input: semantic_value id { semType: Grok ID }
//output: widget res
//meta.role: panel
export async function getMoltrackPropPanelById(id: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.getMoltrackPropPanelById(id);
}
