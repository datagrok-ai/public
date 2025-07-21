import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function compoundActivityDetailsForTarget(target: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CompoundActivityDetailsForTarget', { target });
  }

  export async function bioactivityDataForBacterialTargetsForOrganism(organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:BioactivityDataForBacterialTargetsForOrganism', { organism });
  }

  export async function cbAllChemblStructures(): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CbAllChemblStructures', {});
  }

  export async function cbFindByMolregno(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CbFindByMolregno', { molregno });
  }

  export async function cbChemblBrowserQuery(substructure: string, molecule_type: string, subname: string, max_phase: number, num_ro5_violations: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CbChemblBrowserQuery', { substructure, molecule_type, subname, max_phase, num_ro5_violations });
  }

  export async function patternSimilaritySearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:PatternSimilaritySearch', { pattern, maxRows });
  }

  export async function patternSimilaritySearchWithThreshold(pattern: string, threshold: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:PatternSimilaritySearchWithThreshold', { pattern, threshold });
  }

  export async function patternSubstructureSearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:PatternSubstructureSearch', { pattern, maxRows });
  }

  export async function chemblNumberOfStructures(maxNumberOfMolecules: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblNumberOfStructures', { maxNumberOfMolecules });
  }

  export async function chemblMolregNoBySmiles(smiles: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblMolregNoBySmiles', { smiles });
  }

  export async function structuresByOrganism(maxNumberOfMolecules: number, organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:StructuresByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function chemblIdToSmiles(id: string): Promise<string> {
    return await grok.data.query('ChEMBL:ChemblIdToSmiles', { id });
  }

  export async function molregnoToSmiles(molregno: number): Promise<string> {
    return await grok.data.query('ChEMBL:MolregnoToSmiles', { molregno });
  }

  export async function nameToSmiles(compoundName: string): Promise<string> {
    return await grok.data.query('ChEMBL:NameToSmiles', { compoundName });
  }

  export async function namesToSmiles(): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:NamesToSmiles', {});
  }

  export async function inchiKeyToChembl(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:InchiKeyToChembl', { ids });
  }

  export async function inchiKeyToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:InchiKeyToSmiles', { ids });
  }

  export async function inchiKeyToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:InchiKeyToInchi', { ids });
  }

  export async function chemblToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblToSmiles', { ids });
  }

  export async function chemblToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblToInchi', { ids });
  }

  export async function chemblToInchiKey(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblToInchiKey', { ids });
  }

  export async function pkDataFromCuratedDrugPharmacokineticDataSourceForDrug(drug: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:PKDataFromCuratedDrugPharmacokineticDataSourceForDrug', { drug });
  }

  export async function proteinClassification(): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ProteinClassification', {});
  }

  export async function compoundsWhichAreSelectiveToOneTargetOverASecondTarget(selectiveFor: string, over: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CompoundsWhichAreSelectiveToOneTargetOverASecondTarget', { selectiveFor, over });
  }

  export async function compoundActivityDetailsForAllTargetsContainingProtein(protein: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CompoundActivityDetailsForAllTargetsContainingProtein', { protein });
  }

  export async function unichemUnitTestQuery(): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:UnichemUnitTestQuery', {});
  }

  export async function fracClassification(level1: string, level2: string, level3: string, level4: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:FracClassification', { level1, level2, level3, level4 });
  }

  export async function queryBySubstructure(substructure: string, threshold: string, actionType: string, mechanismOfAction: string, country: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:QueryBySubstructure', { substructure, threshold, actionType, mechanismOfAction, country });
  }

  export async function byChemblIds(): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ByChemblIds', {});
  }

  export async function molregnoInfo(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:MolregnoInfo', { molregno });
  }

  export async function chemblInfo(chemblId: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ChemblInfo', { chemblId });
  }

  export async function fracClassificationWithSubstructure(level1: string, level2: string, level3: string, level4: string, substructure: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:FracClassificationWithSubstructure', { level1, level2, level3, level4, substructure });
  }

  export async function compoundNames(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:CompoundNames', { sub });
  }

  export async function organisms(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:Organisms', { sub });
  }

  export async function proteinTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:ProteinTypes', { sub });
  }

  export async function targetTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:TargetTypes', { sub });
  }

  export async function assayTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:AssayTypes', { sub });
  }

  export async function relationshipTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('ChEMBL:RelationshipTypes', { sub });
  }
}

export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('ChEMBL:Init', {});
  }

  export async function chemblSearchWidgetLocalDb(mol: string, substructure: boolean): Promise<any> {
    return await grok.functions.call('ChEMBL:ChemblSearchWidgetLocalDb', { mol, substructure });
  }

  export async function chemblSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChEMBL:ChemblSubstructureSearchPanel', { mol });
  }

  export async function chemblSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChEMBL:ChemblSimilaritySearchPanel', { mol });
  }

  export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string): Promise<any> {
    return await grok.functions.call('ChEMBL:GetChemblCompoundsByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function getChemblCompounds(maxNumberOfMolecules: number): Promise<any> {
    return await grok.functions.call('ChEMBL:GetChemblCompounds', { maxNumberOfMolecules });
  }

  export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('ChEMBL:ChemblMolregno', { table, molecules });
  }

  export async function chemblIdToSmilesTs(id: string): Promise<any> {
    return await grok.functions.call('ChEMBL:ChemblIdToSmilesTs', { id });
  }

  //Running various queries to chemical databases using convenient input forms
  export async function demoDatabasesChembl(): Promise<any> {
    return await grok.functions.call('ChEMBL:DemoDatabasesChembl', {});
  }
}
