import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: ChEMBL Search Widget
//input: string mol { semType: Molecule }
//input: string searchType 
//output: widget result
//meta.role: widgets
export async function chemblSearchWidget(mol: string, substructure: boolean) : Promise<any> {
  return await PackageFunctions.chemblSearchWidget(mol, substructure);
}

//name: Databases | ChEMBL | Substructure Search API
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function chemblSubstructureSearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.chemblSubstructureSearchPanel(mol);
}

//name: Databases | ChEMBL | Similarity Search API
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
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
//input: string id 
//output: dataframe result
export async function getById(id: string) : Promise<any> {
  return await PackageFunctions.getById(id);
}
