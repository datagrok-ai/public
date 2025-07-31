import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: init
//tags: autostart
//output: dynamic result
export function init() {
  return PackageFunctions.init();
}

//name: chemblSearchWidgetLocalDb
//tags: widgets
//input: string mol { semType: Molecule }
//input: bool substructure 
//output: widget result
export async function chemblSearchWidgetLocalDb(mol: string, substructure: boolean) {
  return PackageFunctions.chemblSearchWidgetLocalDb(mol, substructure);
}

//name: Databases | ChEMBL | Substructure Search (Internal)
//tags: panel, widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function chemblSubstructureSearchPanel(mol: string) {
  return PackageFunctions.chemblSubstructureSearchPanel(mol);
}

//name: Databases | ChEMBL | Similarity Search (Internal)
//tags: widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function chemblSimilaritySearchPanel(mol: string) {
  return PackageFunctions.chemblSimilaritySearchPanel(mol);
}

//name: Chembl targets by organism
//tags: HitTriageDataSource
//input: int maxNumberOfMolecules { default: 1000 }
//input: string organism { default: 'Shigella' }
//output: dataframe result
export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string) {
  return PackageFunctions.getChemblCompoundsByOrganism(maxNumberOfMolecules, organism);
}

//name: Chembl Compounds
//tags: HitTriageDataSource
//input: int maxNumberOfMolecules { default: 1000 }
//output: dataframe result
export async function getChemblCompounds(maxNumberOfMolecules: number) {
  return PackageFunctions.getChemblCompounds(maxNumberOfMolecules);
}

//name: Chembl molregno
//tags: HitTriageFunction
//input: dataframe table { caption: Table }
//input: column molecules { caption: Molecules; semType: Molecule }
//output: dataframe result
export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column) {
  return PackageFunctions.chemblMolregno(table, molecules);
}

//name: chemblIdToSmilesTs
//input: string id { semType: CHEMBL_ID; default: 'CHEMBL1185' }
//output: string result
//meta.role: converter
//meta.inputRegexp: (CHEMBL[0-9]+)
export async function chemblIdToSmilesTs(id: string) {
  return PackageFunctions.chemblIdToSmilesTs(id);
}

//name: Database Queries
//description: Running various queries to chemical databases using convenient input forms
//meta.demoPath: Cheminformatics | Database Queries
export async function demoDatabasesChembl() {
  return PackageFunctions.demoDatabasesChembl();
}
