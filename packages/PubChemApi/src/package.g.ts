import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Databases | PubChem | Substructure Search
//tags: panel, widgets
//input: string molString { semType: Molecule }
//output: widget result
export async function pubChemSubstructureSearchPanel(molString: string) {
  return PackageFunctions.pubChemSubstructureSearchPanel(molString);
}

//name: Databases | PubChem | Similarity Search
//tags: panel, widgets
//input: string molString { semType: Molecule }
//output: widget result
export async function pubChemSimilaritySearchPanel(molString: string) {
  return PackageFunctions.pubChemSimilaritySearchPanel(molString);
}

//name: Databases | PubChem | Identity Search
//tags: panel, widgets
//input: string molString { semType: Molecule }
//output: widget result
export async function pubChemIdentitySearch(molString: string) {
  return PackageFunctions.pubChemIdentitySearch(molString);
}

//name: pubChemToSmiles
//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: (^\s*[Pp][Uu][Bb][Cc][Hh][Ee][Mm]\s*\:\s*[0-9]+\s*$)
//meta.connection: PubChemApi
export async function pubChemToSmiles(id: string) {
  return PackageFunctions.pubChemToSmiles(id);
}

//name: inchiKeysToSmiles
//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: ([A-Z]{14}-[A-Z]{10}-N)
//meta.connection: PubChemApi
export async function inchiKeysToSmiles(id: string) {
  return PackageFunctions.inchiKeysToSmiles(id);
}

//name: GetIupacName
//input: string smiles 
//output: string result
//meta.connection: PubChemApi
export async function GetIupacName(smiles: string) {
  return PackageFunctions.GetIupacName(smiles);
}
