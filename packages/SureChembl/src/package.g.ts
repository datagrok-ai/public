import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: sureChemblSubstructureSearch
//input: string molecule { semType: Molecule }
//input: int limit 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSubstructureSearch(molecule: string, limit: number) {
  return PackageFunctions.sureChemblSubstructureSearch(molecule, limit);
}

//name: sureChemblSimilaritySearch
//input: string molecule { semType: Molecule }
//input: int limit 
//input: double similarityThreshold 
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold: number) {
  return PackageFunctions.sureChemblSimilaritySearch(molecule, limit, similarityThreshold);
}

//name: Databases | SureChEMBL | Substructure Search
//tags: panel, widgets
//input: string molecule { semType: Molecule }
//output: widget result
//condition: true
export function sureChemblSubstructureSearchWidget(molecule: string) {
  return PackageFunctions.sureChemblSubstructureSearchWidget(molecule);
}

//name: Databases | SureChEMBL | Similarity Search
//tags: panel, widgets
//input: string molecule { semType: Molecule }
//output: widget result
//condition: true
export function sureChemblSimilaritySearchWidget(molecule: string) {
  return PackageFunctions.sureChemblSimilaritySearchWidget(molecule);
}
