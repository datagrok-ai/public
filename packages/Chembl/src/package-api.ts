import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Queries {
  export async function compoundActivityDetailsForTarget(target: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CompoundActivityDetailsForTarget', { target });
  }

  export async function bioactivityDataForBacterialTargetsForOrganism(organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:BioactivityDataForBacterialTargetsForOrganism', { organism });
  }

  export async function cbAllChemblStructures(): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CbAllChemblStructures', {});
  }

  export async function cbFindByMolregno(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CbFindByMolregno', { molregno });
  }

  export async function cbChemblBrowserQuery(substructure: string, molecule_type: string, subname: string, max_phase: number, num_ro5_violations: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CbChemblBrowserQuery', { substructure, molecule_type, subname, max_phase, num_ro5_violations });
  }

  export async function patternSimilaritySearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:PatternSimilaritySearch', { pattern, maxRows });
  }

  export async function patternSimilaritySearchWithThreshold(pattern: string, threshold: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:PatternSimilaritySearchWithThreshold', { pattern, threshold });
  }

  export async function patternSubstructureSearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:PatternSubstructureSearch', { pattern, maxRows });
  }

  export async function chemblNumberOfStructures(maxNumberOfMolecules: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblNumberOfStructures', { maxNumberOfMolecules });
  }

  export async function chemblMolregNoBySmiles(smiles: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblMolregNoBySmiles', { smiles });
  }

  export async function structuresByOrganism(maxNumberOfMolecules: number, organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:StructuresByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function chemblIdToSmiles(id: string): Promise<string> {
    return await grok.data.query('Chembl:ChemblIdToSmiles', { id });
  }

  export async function molregnoToSmiles(molregno: number): Promise<string> {
    return await grok.data.query('Chembl:MolregnoToSmiles', { molregno });
  }

  export async function nameToSmiles(compoundName: string): Promise<string> {
    return await grok.data.query('Chembl:NameToSmiles', { compoundName });
  }

  export async function namesToSmiles(): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:NamesToSmiles', {});
  }

  export async function inchiKeyToChembl(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:InchiKeyToChembl', { ids });
  }

  export async function inchiKeyToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:InchiKeyToSmiles', { ids });
  }

  export async function inchiKeyToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:InchiKeyToInchi', { ids });
  }

  export async function chemblToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblToSmiles', { ids });
  }

  export async function chemblToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblToInchi', { ids });
  }

  export async function chemblToInchiKey(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblToInchiKey', { ids });
  }

  export async function pkDataFromCuratedDrugPharmacokineticDataSourceForDrug(drug: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:PKDataFromCuratedDrugPharmacokineticDataSourceForDrug', { drug });
  }

  export async function proteinClassification(): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ProteinClassification', {});
  }

  export async function compoundsWhichAreSelectiveToOneTargetOverASecondTarget(selectiveFor: string, over: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CompoundsWhichAreSelectiveToOneTargetOverASecondTarget', { selectiveFor, over });
  }

  export async function compoundActivityDetailsForAllTargetsContainingProtein(protein: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CompoundActivityDetailsForAllTargetsContainingProtein', { protein });
  }

  export async function unichemUnitTestQuery(): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:UnichemUnitTestQuery', {});
  }

  export async function fracClassification(level1: string, level2: string, level3: string, level4: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:FracClassification', { level1, level2, level3, level4 });
  }

  export async function queryBySubstructure(substructure: string, threshold: string, actionType: string, mechanismOfAction: string, country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:QueryBySubstructure', { substructure, threshold, actionType, mechanismOfAction, country });
  }

  export async function byChemblIds(): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ByChemblIds', {});
  }

  export async function molregnoInfo(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:MolregnoInfo', { molregno });
  }

  export async function chemblInfo(chemblId: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ChemblInfo', { chemblId });
  }

  export async function fracClassificationWithSubstructure(level1: string, level2: string, level3: string, level4: string, substructure: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:FracClassificationWithSubstructure', { level1, level2, level3, level4, substructure });
  }

  export async function compoundNames(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:CompoundNames', { sub });
  }

  export async function organisms(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:Organisms', { sub });
  }

  export async function proteinTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:ProteinTypes', { sub });
  }

  export async function targetTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:TargetTypes', { sub });
  }

  export async function assayTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:AssayTypes', { sub });
  }

  export async function relationshipTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('Chembl:RelationshipTypes', { sub });
  }
}

export namespace Funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('Chembl:Init', {});
  }

  export async function chemblSearchWidgetLocalDb(mol: string, substructure: boolean): Promise<any> {
    return await grok.functions.call('Chembl:ChemblSearchWidgetLocalDb', { mol, substructure });
  }

  export async function chemblSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('Chembl:ChemblSubstructureSearchPanel', { mol });
  }

  export async function chemblSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('Chembl:ChemblSimilaritySearchPanel', { mol });
  }

  export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string): Promise<any> {
    return await grok.functions.call('Chembl:GetChemblCompoundsByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function getChemblCompounds(maxNumberOfMolecules: number): Promise<any> {
    return await grok.functions.call('Chembl:GetChemblCompounds', { maxNumberOfMolecules });
  }

  export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('Chembl:ChemblMolregno', { table, molecules });
  }

  export async function chemblIdToSmilesTs(id: string): Promise<any> {
    return await grok.functions.call('Chembl:ChemblIdToSmilesTs', { id });
  }

  //Running various queries to chemical databases using convenient input forms
  export async function demoDatabasesChembl(): Promise<any> {
    return await grok.functions.call('Chembl:DemoDatabasesChembl', {});
  }
}
