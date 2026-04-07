import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Databases | PubChem | Substructure Search
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function pubChemSubstructureSearchPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemSubstructureSearchPanel(molString);
}

//name: Databases | PubChem | Similarity Search
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function pubChemSimilaritySearchPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemSimilaritySearchPanel(molString);
}

//name: Databases | PubChem | Identity Search
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function pubChemIdentitySearch(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemIdentitySearch(molString);
}

//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: (^\s*[Pp][Uu][Bb][Cc][Hh][Ee][Mm]\s*\:\s*[0-9]+\s*$)
//connection: PubChemApi
export async function pubChemToSmiles(id: string) : Promise<string> {
  return await PackageFunctions.pubChemToSmiles(id);
}

//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: ([A-Z]{14}-[A-Z]{10}-N)
//connection: PubChemApi
export async function inchiKeysToSmiles(id: string) : Promise<string> {
  return await PackageFunctions.inchiKeysToSmiles(id);
}

//input: string smiles 
//output: string result
//connection: PubChemApi
export async function GetIupacName(smiles: string) : Promise<string> {
  return await PackageFunctions.GetIupacName(smiles);
}
