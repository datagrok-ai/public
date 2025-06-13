import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function compoundActivityDetailsForTarget(target: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CompoundActivityDetailsForTarget', { target });
  }

  export async function bioactivityDataForBacterialTargetsForOrganism(organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:BioactivityDataForBacterialTargetsForOrganism', { organism });
  }

  export async function cbAllChemblStructures(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CbAllChemblStructures', {});
  }

  export async function cbFindByMolregno(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CbFindByMolregno', { molregno });
  }

  export async function cbChemblBrowserQuery(substructure: string, molecule_type: string, subname: string, max_phase: number, num_ro5_violations: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CbChemblBrowserQuery', { substructure, molecule_type, subname, max_phase, num_ro5_violations });
  }

  export async function patternSimilaritySearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:PatternSimilaritySearch', { pattern, maxRows });
  }

  export async function patternSimilaritySearchWithThreshold(pattern: string, threshold: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:PatternSimilaritySearchWithThreshold', { pattern, threshold });
  }

  export async function patternSubstructureSearch(pattern: string, maxRows: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:PatternSubstructureSearch', { pattern, maxRows });
  }

  export async function chemblNumberOfStructures(maxNumberOfMolecules: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblNumberOfStructures', { maxNumberOfMolecules });
  }

  export async function chemblMolregNoBySmiles(smiles: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblMolregNoBySmiles', { smiles });
  }

  export async function structuresByOrganism(maxNumberOfMolecules: number, organism: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:StructuresByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function chemblIdToSmiles(id: string): Promise<string> {
    return await grok.data.query('@datagrok/chembl:ChemblIdToSmiles', { id });
  }

  export async function molregnoToSmiles(molregno: number): Promise<string> {
    return await grok.data.query('@datagrok/chembl:MolregnoToSmiles', { molregno });
  }

  export async function nameToSmiles(compoundName: string): Promise<string> {
    return await grok.data.query('@datagrok/chembl:NameToSmiles', { compoundName });
  }

  export async function namesToSmiles(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:NamesToSmiles', {});
  }

  export async function inchiKeyToChembl(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:InchiKeyToChembl', { ids });
  }

  export async function inchiKeyToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:InchiKeyToSmiles', { ids });
  }

  export async function inchiKeyToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:InchiKeyToInchi', { ids });
  }

  export async function chemblToSmiles(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblToSmiles', { ids });
  }

  export async function chemblToInchi(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblToInchi', { ids });
  }

  export async function chemblToInchiKey(ids: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblToInchiKey', { ids });
  }

  export async function pkDataFromCuratedDrugPharmacokineticDataSourceForDrug(drug: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:PKDataFromCuratedDrugPharmacokineticDataSourceForDrug', { drug });
  }

  export async function proteinClassification(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ProteinClassification', {});
  }

  export async function compoundsWhichAreSelectiveToOneTargetOverASecondTarget(selectiveFor: string, over: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CompoundsWhichAreSelectiveToOneTargetOverASecondTarget', { selectiveFor, over });
  }

  export async function compoundActivityDetailsForAllTargetsContainingProtein(protein: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CompoundActivityDetailsForAllTargetsContainingProtein', { protein });
  }

  export async function unichemUnitTestQuery(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:UnichemUnitTestQuery', {});
  }

  export async function fracClassification(level1: string, level2: string, level3: string, level4: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:FracClassification', { level1, level2, level3, level4 });
  }

  export async function queryBySubstructure(substructure: string, threshold: string, actionType: string, mechanismOfAction: string, country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:QueryBySubstructure', { substructure, threshold, actionType, mechanismOfAction, country });
  }

  export async function byChemblIds(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ByChemblIds', {});
  }

  export async function molregnoInfo(molregno: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:MolregnoInfo', { molregno });
  }

  export async function chemblInfo(chemblId: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ChemblInfo', { chemblId });
  }

  export async function fracClassificationWithSubstructure(level1: string, level2: string, level3: string, level4: string, substructure: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:FracClassificationWithSubstructure', { level1, level2, level3, level4, substructure });
  }

  export async function compoundNames(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:CompoundNames', { sub });
  }

  export async function organisms(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:Organisms', { sub });
  }

  export async function proteinTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:ProteinTypes', { sub });
  }

  export async function targetTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:TargetTypes', { sub });
  }

  export async function assayTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:AssayTypes', { sub });
  }

  export async function relationshipTypes(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/chembl:RelationshipTypes', { sub });
  }
}

export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:Init', {});
  }

  export async function chemblSearchWidgetLocalDb(mol: string, substructure: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:ChemblSearchWidgetLocalDb', { mol, substructure });
  }

  export async function chemblSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:ChemblSubstructureSearchPanel', { mol });
  }

  export async function chemblSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:ChemblSimilaritySearchPanel', { mol });
  }

  export async function getChemblCompoundsByOrganism(maxNumberOfMolecules: number, organism: string): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:GetChemblCompoundsByOrganism', { maxNumberOfMolecules, organism });
  }

  export async function getChemblCompounds(maxNumberOfMolecules: number): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:GetChemblCompounds', { maxNumberOfMolecules });
  }

  export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:ChemblMolregno', { table, molecules });
  }

  export async function chemblIdToSmilesTs(id: string): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:ChemblIdToSmilesTs', { id });
  }

  //Running various queries to chemical databases using convenient input forms
  export async function demoDatabasesChembl(): Promise<any> {
    return await grok.functions.call('@datagrok/chembl:DemoDatabasesChembl', {});
  }
}
