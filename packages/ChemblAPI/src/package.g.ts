import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: ChEMBL Search Widget
//tags: widgets
//input: string mol { semType: Molecule }
//input: string searchType 
//output: widget result
export async function chemblSearchWidget(mol: string, substructure: boolean) {
  return PackageFunctions.chemblSearchWidget(mol, substructure);
}

//name: Databases | ChEMBL | Substructure Search API
//tags: panel, widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function chemblSubstructureSearchPanel(mol: string) {
  return PackageFunctions.chemblSubstructureSearchPanel(mol);
}

//name: Databases | ChEMBL | Similarity Search API
//tags: panel, widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function chemblSimilaritySearchPanel(mol: string) {
  return PackageFunctions.chemblSimilaritySearchPanel(mol);
}

//name: getCompoundsIds
//input: string inchiKey 
//output: dynamic result
export async function getCompoundsIds(inchiKey: string) {
  return PackageFunctions.getCompoundsIds(inchiKey);
}
