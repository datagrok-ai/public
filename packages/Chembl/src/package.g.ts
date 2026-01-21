import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: autostart
export function init() : void {
  PackageFunctions.init();
}

//input: string mol { semType: Molecule }
//input: bool substructure 
//output: widget result
//meta.role: widgets
export async function chemblSearchWidgetLocalDb(mol: string, substructure: boolean) : Promise<any> {
  return await PackageFunctions.chemblSearchWidgetLocalDb(mol, substructure);
}

//name: Databases | ChEMBL | Substructure Search (Internal)
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function chemblSubstructureSearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.chemblSubstructureSearchPanel(mol);
}

//name: Databases | ChEMBL | Similarity Search (Internal)
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function chemblSimilaritySearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.chemblSimilaritySearchPanel(mol);
}

//name: Chembl targets by organism
//input: int maxNumberOfMolecules = 1000 { description: Maximum number of rows to return }
//input: string organism = 'Shigella' { description: Organism name }
//output: dataframe result
//meta.role: hitTriageDataSource
export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string) : Promise<any> {
  return await PackageFunctions.getChemblCompoundsByOrganism(maxNumberOfMolecules, organism);
}

//name: Chembl Compounds
//input: int maxNumberOfMolecules = 1000 { description: Maximum number of rows to return }
//output: dataframe result
//meta.role: hitTriageDataSource
export async function getChemblCompounds(maxNumberOfMolecules: number) : Promise<any> {
  return await PackageFunctions.getChemblCompounds(maxNumberOfMolecules);
}

//name: Chembl molregno
//input: dataframe table { caption: Table; description: Input data table }
//input: column molecules { caption: Molecules; semType: Molecule }
//output: dataframe result
//meta.role: hitTriageFunction
export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column) : Promise<any> {
  return await PackageFunctions.chemblMolregno(table, molecules);
}

//input: string id = 'CHEMBL1185' { semType: CHEMBL_ID }
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: (CHEMBL[0-9]+)
export async function chemblIdToSmilesTs(id: string) : Promise<string> {
  return await PackageFunctions.chemblIdToSmilesTs(id);
}

//name: Database Queries
//description: Running various queries to chemical databases using convenient input forms
//meta.demoPath: Cheminformatics | Database Queries
export async function demoDatabasesChembl() : Promise<void> {
  await PackageFunctions.demoDatabasesChembl();
}
