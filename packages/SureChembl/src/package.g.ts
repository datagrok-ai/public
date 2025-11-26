import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: string molecule { semType: Molecule }
//input: int limit 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSubstructureSearch(molecule: string, limit: number) : Promise<any> {
  return await PackageFunctions.sureChemblSubstructureSearch(molecule, limit);
}

//input: string molecule { semType: Molecule }
//input: int limit 
//input: double similarityThreshold 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold?: number) : Promise<any> {
  return await PackageFunctions.sureChemblSimilaritySearch(molecule, limit, similarityThreshold);
}

//name: Databases | SureChEMBL | Substructure Search
//tags: panel, widgets
//input: string molecule { semType: Molecule }
//output: widget result
//condition: true
export function sureChemblSubstructureSearchWidget(molecule: string) : any {
  return PackageFunctions.sureChemblSubstructureSearchWidget(molecule);
}

//name: Databases | SureChEMBL | Similarity Search
//tags: panel, widgets
//input: string molecule { semType: Molecule }
//output: widget result
//condition: true
export function sureChemblSimilaritySearchWidget(molecule: string) : any {
  return PackageFunctions.sureChemblSimilaritySearchWidget(molecule);
}
