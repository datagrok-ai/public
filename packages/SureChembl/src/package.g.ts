import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: SureChEMBL Substructure Search
//description: Searches SureCHEMBL patents for molecules containing the query structure as a substructure.
//input: string molecule { semType: Molecule }
//input: int limit { description: Maximum number of matching molecules to return }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSubstructureSearch(molecule: string, limit: number) : Promise<any> {
  return await PackageFunctions.sureChemblSubstructureSearch(molecule, limit);
}

//name: SureChEMBL Similarity Search
//description: Searches SureCHEMBL patents for molecules similar to the query structure.
//input: string molecule { semType: Molecule }
//input: int limit { description: Maximum number of matching molecules to return }
//input: double similarityThreshold { description: Minimum Tanimoto similarity, 0-1 (default 0.6) }
//output: dataframe result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold?: number) : Promise<any> {
  return await PackageFunctions.sureChemblSimilaritySearch(molecule, limit, similarityThreshold);
}

//name: Databases | SureChEMBL | Substructure Search
//input: string molecule { semType: Molecule }
//output: widget result
//meta.role: panel,widgets
//condition: true
export function sureChemblSubstructureSearchWidget(molecule: string) : any {
  return PackageFunctions.sureChemblSubstructureSearchWidget(molecule);
}

//name: Databases | SureChEMBL | Similarity Search
//input: string molecule { semType: Molecule }
//output: widget result
//meta.role: panel,widgets
//condition: true
export function sureChemblSimilaritySearchWidget(molecule: string) : any {
  return PackageFunctions.sureChemblSimilaritySearchWidget(molecule);
}
