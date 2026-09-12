import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Databases | PubChem | Info
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function pubChemPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemPanel(molString);
}

//name: Databases | PubChem | Substructure Search
//description: Finds PubChem compounds that contain the query structure as a substructure.
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function pubChemSubstructureSearchPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemSubstructureSearchPanel(molString);
}

//name: Databases | PubChem | Similarity Search
//description: Finds PubChem compounds most similar to the query structure.
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function pubChemSimilaritySearchPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemSimilaritySearchPanel(molString);
}

//name: Databases | PubChem | Identity Search
//description: Finds the PubChem compound identical to the query structure.
//input: string molString { semType: Molecule }
//output: widget result
//meta.role: widgets,Panel
export async function pubChemIdentitySearchPanel(molString: string) : Promise<any> {
  return await PackageFunctions.pubChemIdentitySearchPanel(molString);
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
