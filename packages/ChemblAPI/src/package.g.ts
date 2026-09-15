import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Databases | ChEMBL | Substructure Search API
//description: Finds ChEMBL molecules that contain the query structure as a substructure.
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function chemblSubstructureSearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.chemblSubstructureSearchPanel(mol);
}

//name: Databases | ChEMBL | Similarity Search API
//description: Finds ChEMBL molecules most similar to the query structure.
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function chemblSimilaritySearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.chemblSimilaritySearchPanel(mol);
}

//name: GetCompoundsIds
//input: string inchiKey 
//output: object result
export async function getCompoundsIds(inchiKey: string) : Promise<any> {
  return await PackageFunctions.getCompoundsIds(inchiKey);
}

//name: Chembl Get by Id
//description: Fetches a molecule record from ChEMBL by its ChEMBL ID.
//input: string id { description: ChEMBL ID (e.g. CHEMBL25) a bare number is prefixed with CHEMBL }
//output: dataframe result
export async function getById(id: string) : Promise<any> {
  return await PackageFunctions.getById(id);
}
