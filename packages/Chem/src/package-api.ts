import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  //Standardizes the dataset
  export async function curate(data: DG.DataFrame, molecules: DG.Column, kekulization: boolean, normalization: boolean, reionization: boolean, neutralization: boolean, tautomerization: boolean, mainFragment: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:Curate', { data, molecules, kekulization, normalization, reionization, neutralization, tautomerization, mainFragment });
  }

  export async function desc(smiles: string, df1: DG.DataFrame, selected: string, df2: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:Desc', { smiles, df1, selected, df2 });
  }

  export async function findMCS(molecules: string, df: DG.DataFrame, exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
    return await grok.functions.call('@datagrok/chem:FindMCS', { molecules, df, exactAtomSearch, exactBondSearch });
  }

  export async function findRGroupsWithCore(molecules: string, df: DG.DataFrame, core: string, onlyMatchAtRGroups: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:FindRGroupsWithCore', { molecules, df, core, onlyMatchAtRGroups });
  }

  export async function findRGroups(molecules: string, df: DG.DataFrame, core: string, prefix: string): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:FindRGroups', { molecules, df, core, prefix });
  }

  //RDKit-based script.
  export async function chemistryGasteigerPartialCharges(mol: string, contours: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemistryGasteigerPartialCharges', { mol, contours });
  }

  export async function inchiToMol(id: string): Promise<string> {
    return await grok.functions.call('@datagrok/chem:InchiToMol', { id });
  }

  //Generates the new dataset based on the given structure
  export async function mutate(molecule: string, steps: number, randomize: boolean, maxRandomResults: number): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:Mutate', { molecule, steps, randomize, maxRandomResults });
  }

  //generation scaffold tree from dataset
  export async function generateScaffoldTree(data: DG.DataFrame, smilesColumn: string, ringCutoff: number, dischargeAndDeradicalize: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GenerateScaffoldTree', { data, smilesColumn, ringCutoff, dischargeAndDeradicalize });
  }

  export async function smilesTo3DCoordinates(molecule: string): Promise<string> {
    return await grok.functions.call('@datagrok/chem:SmilesTo3DCoordinates', { molecule });
  }

  export async function amideReaction(amines: DG.DataFrame, amine_molecules: DG.Column, acids: DG.DataFrame, acid_molecules: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:AmideReaction', { amines, amine_molecules, acids, acid_molecules });
  }

  //Implementation of the clustering algorithm published in: Butina JCICS 39 747-750 (1999)
  export async function butinaMoleculesClustering(data: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:ButinaMoleculesClustering', { data, molecules });
  }

  //USRCAT - real-time ultrafast shape recognition with pharmacophoric constraints
  export async function usrcat(data: DG.DataFrame, smiles: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:USRCAT', { data, smiles });
  }

  //Finds undesireable molecules based on various criteria
  export async function filterByCatalogs(data: DG.DataFrame, smiles: DG.Column, catalog: string): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:FilterByCatalogs', { data, smiles, catalog });
  }

  //Generation of Murcko scaffolds from a molecule
  export async function murckoScaffolds(data: DG.DataFrame, smiles: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:MurckoScaffolds', { data, smiles });
  }

  //Similarity Maps Using Fingerprints, RDKit based
  export async function similarityMapsUsingFingerprints(mol: string, refmol: string, radius: number): Promise<number> {
    return await grok.functions.call('@datagrok/chem:SimilarityMapsUsingFingerprints', { mol, refmol, radius });
  }

  //Chemical space using t-distributed Stochastic Neighbor Embedding
  export async function chemicalSpaceUsingTSNE(data: DG.DataFrame, smiles: DG.Column, components: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemicalSpaceUsingTSNE', { data, smiles, components, minClusterSize });
  }

  //Two component reaction
  export async function twoComponentReaction(data1: DG.DataFrame, reactants1: DG.Column, data2: DG.DataFrame, reactants2: DG.Column, reaction: string, matrixExpansion: boolean, randomize: boolean, seed: number, maxRandomReactions: number): Promise<DG.DataFrame> {
    return await grok.functions.call('@datagrok/chem:TwoComponentReaction', { data1, reactants1, data2, reactants2, reaction, matrixExpansion, randomize, seed, maxRandomReactions });
  }

  //Chemical space using Uniform Manifold Approximation and Projection
  export async function chemicalSpaceUsingUMAP(data: DG.DataFrame, smiles: DG.Column, neighbors: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemicalSpaceUsingUMAP', { data, smiles, neighbors, minClusterSize });
  }

  //to be used in tests to ensure JKG is up and running
  export async function testPythonRunning(x: number, y: number): Promise<number> {
    return await grok.functions.call('@datagrok/chem:TestPythonRunning', { x, y });
  }
}

export namespace funcs {
  export async function getRdKitModule(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetRdKitModule', {});
  }

  export async function getMolFileHandler(molString: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMolFileHandler', { molString });
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Init', {});
  }

  export async function initChemAutostart(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:InitChemAutostart', {});
  }

  export async function chemTooltip(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemTooltip', { col });
  }

  export async function scaffoldTreeViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ScaffoldTreeViewer', {});
  }

  //RDKit-based substructure filter
  export async function substructureFilter(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SubstructureFilter', {});
  }

  export async function canvasMol(x: number, y: number, w: number, h: number, canvas: any, molString: string, scaffoldMolString: string, options: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CanvasMol', { x, y, w, h, canvas, molString, scaffoldMolString, options });
  }

  export async function drawMolecule(molStr: string, w: number, h: number, popupMenu: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DrawMolecule', { molStr, w, h, popupMenu });
  }

  export async function getCLogP(smiles: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetCLogP', { smiles });
  }

  export async function rdKitCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RdKitCellRenderer', {});
  }

  export async function rdKitReactionRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RdKitReactionRenderer', {});
  }

  export async function chemCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemCellRenderer', {});
  }

  export async function getMorganFingerprints(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMorganFingerprints', {});
  }

  export async function getMorganFingerprint(molString: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMorganFingerprint', { molString });
  }

  export async function getSimilarities(molStringsColumn: DG.Column, molString: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetSimilarities', { molStringsColumn, molString });
  }

  export async function getDiversities(molStringsColumn: DG.Column, limit: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetDiversities', { molStringsColumn, limit });
  }

  export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number, cutoff: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:FindSimilar', { molStringsColumn, molString, limit, cutoff });
  }

  export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, molBlockFailover: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SearchSubstructure', { molStringsColumn, molString, molBlockFailover });
  }

  //As SDF
  export async function saveAsSdf(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SaveAsSdf', {});
  }

  export async function similaritySearchViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SimilaritySearchViewer', {});
  }

  export async function similaritySearchTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SimilaritySearchTopMenu', {});
  }

  export async function diversitySearchViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DiversitySearchViewer', {});
  }

  export async function diversitySearchTopMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DiversitySearchTopMenu', {});
  }

  export async function descriptorsDocker(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DescriptorsDocker', {});
  }

  export async function chemDescriptorsTree(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemDescriptorsTree', {});
  }

  export async function getMapIdentifiers(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMapIdentifiers', {});
  }

  export async function freeTextToSmiles(molfile: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:FreeTextToSmiles', { molfile });
  }

  export async function chemDescriptors(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemDescriptors', { table, molecules });
  }

  export async function chemDescriptor(descriptor: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemDescriptor', { descriptor });
  }

  export async function searchSubstructureEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SearchSubstructureEditor', { call });
  }

  export async function substructureSearchTopMenu(molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SubstructureSearchTopMenu', { molecules });
  }

  //Calculates most common substructures for each cluster
  export async function clusterMCSTopMenu(table: DG.DataFrame, molCol: DG.Column, clusterCol: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ClusterMCSTopMenu', { table, molCol, clusterCol });
  }

  //Calculates most common substructures for each cluster
  export async function performClusterMCS(molCol: DG.Column, clusterCol: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:PerformClusterMCS', { molCol, clusterCol });
  }

  export async function chemSpaceEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemSpaceEditor', { call });
  }

  export async function getFingerprints(col: DG.Column, _metric: string, fingerprintType: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetFingerprints', { col, _metric, fingerprintType });
  }

  //Maps the dataset to 2D plot based on similarity
  export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: string, similarityMetric: string, plotEmbeddings: boolean, options: any, preprocessingFunction: any, clusterEmbeddings: boolean, clusterMCS: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemSpaceTopMenu', { table, molecules, methodName, similarityMetric, plotEmbeddings, options, preprocessingFunction, clusterEmbeddings, clusterMCS });
  }

  export async function chemSpaceTransform(table: DG.DataFrame, molecules: DG.Column, methodName: string, similarityMetric: string, plotEmbeddings: boolean, options: string, clusterEmbeddings: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ChemSpaceTransform', { table, molecules, methodName, similarityMetric, plotEmbeddings, options, clusterEmbeddings });
  }

  export async function getChemSpaceEmbeddings(col: string, methodName: string, similarityMetric: string, xAxis: string, yAxis: string, options: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetChemSpaceEmbeddings', { col, methodName, similarityMetric, xAxis, yAxis, options });
  }

  export async function getChemSimilaritiesMatrix(dim: number, col: DG.Column, df: DG.DataFrame, colName: string, simArr: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetChemSimilaritiesMatrix', { dim, col, df, colName, simArr });
  }

  export async function elementalAnalysis(table: DG.DataFrame, molecules: DG.Column, radarViewer: boolean, radarGrid: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ElementalAnalysis', { table, molecules, radarViewer, radarGrid });
  }

  export async function runElementalAnalysis(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RunElementalAnalysis', { table, molecules });
  }

  export async function rGroupsAnalysisMenu(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RGroupsAnalysisMenu', {});
  }

  export async function rGroupDecomposition(df: DG.DataFrame, molColName: string, core: string, rGroupName: string, rGroupMatchingStrategy: string, onlyMatchAtRGroups: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RGroupDecomposition', { df, molColName, core, rGroupName, rGroupMatchingStrategy, onlyMatchAtRGroups });
  }

  export async function activityCliffsEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ActivityCliffsEditor', { call });
  }

  //Detects pairs of molecules with similar structure and significant difference in any given property
  export async function activityCliffs(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: string, similarityMetric: string, preprocessingFunction: any, options: any, isDemo: boolean, isTest: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ActivityCliffs', { table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, isDemo, isTest });
  }

  export async function activityCliffsInitFunction(v: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ActivityCliffsInitFunction', { v });
  }

  export async function activityCliffsTransform(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: string, similarityMetric: string, options: string, isDemo: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ActivityCliffsTransform', { table, molecules, activities, similarity, methodName, similarityMetric, options, isDemo });
  }

  export async function addInchisTopMenu(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:AddInchisTopMenu', { table, molecules });
  }

  export async function getInchis(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetInchis', {});
  }

  export async function addInchisKeysTopMenu(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:AddInchisKeysTopMenu', { table, molecules });
  }

  export async function getInchiKeys(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetInchiKeys', {});
  }

  //Highlights the fragments that could lead to potential chemical hazards
  export async function structuralAlertsTopMenu(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:StructuralAlertsTopMenu', { table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo });
  }

  export async function runStructuralAlerts(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RunStructuralAlerts', { table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo });
  }

  export async function runStructuralAlert(alert: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RunStructuralAlert', { alert });
  }

  export async function molColumnPropertyPanel(molColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MolColumnPropertyPanel', { molColumn });
  }

  export async function molColumnHighlights(molColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MolColumnHighlights', { molColumn });
  }

  export async function descriptorsWidget(smiles: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DescriptorsWidget', { smiles });
  }

  //Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
  export async function drugLikeness(smiles: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DrugLikeness', { smiles });
  }

  //Basic molecule properties
  export async function properties(smiles: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Properties', { smiles });
  }

  //Return chem property function
  export async function getChemPropertyFunction(name: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetChemPropertyFunction', { name });
  }

  //Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
  export async function structuralAlerts(smiles: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:StructuralAlerts', { smiles });
  }

  export async function identifiers(smiles: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Identifiers', { smiles });
  }

  //3D molecule representation
  export async function structure3D(molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Structure3D', { molecule });
  }

  //2D molecule representation
  export async function structure2d(molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Structure2d', { molecule });
  }

  //Toxicity prediction. Calculated by openchemlib
  export async function toxicity(smiles: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Toxicity', { smiles });
  }

  export async function convertMoleculeNotation(targetNotation: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ConvertMoleculeNotation', { targetNotation });
  }

  //RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000
  export async function convertMolNotation(molecule: string, sourceNotation: string, targetNotation: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ConvertMolNotation', { molecule, sourceNotation, targetNotation });
  }

  export async function convertNotation(data: DG.DataFrame, molecules: DG.Column, targetNotation: string, overwrite: boolean, join: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ConvertNotation', { data, molecules, targetNotation, overwrite, join });
  }

  export async function convertMolNotationAction(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ConvertMolNotationAction', { col });
  }

  //Molecule
  export async function editMoleculeCell(cell: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:EditMoleculeCell', { cell });
  }

  export async function openChemLibSketcher(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:OpenChemLibSketcher', {});
  }

  //Opens SDF file
  export async function importSdf(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ImportSdf', { bytes });
  }

  //Opens smi file
  export async function importSmi(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ImportSmi', { bytes });
  }

  //Opens smi file
  export async function importMol2(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ImportMol2', { bytes });
  }

  //Opens MOL file
  export async function importMol(content: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ImportMol', { content });
  }

  export async function oclCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:OclCellRenderer', {});
  }

  //Sorts a molecular column by similarity
  export async function sortBySimilarity(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:SortBySimilarity', { value });
  }

  //Adds this structure as a substructure filter
  export async function useAsSubstructureFilter(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:UseAsSubstructureFilter', { value });
  }

  //Copies structure in different formats
  export async function copyAsAction(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsAction', { value });
  }

  //Copies structure as smiles
  export async function copyAsSmiles(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsSmiles', { value });
  }

  //Copies structure as molfile V2000
  export async function copyAsMolfileV2000(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsMolfileV2000', { value });
  }

  //Copies structure as molfile V3000
  export async function copyAsMolfileV3000(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsMolfileV3000', { value });
  }

  //Copies structure as smarts
  export async function copyAsSmarts(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsSmarts', { value });
  }

  //Copies structure as Image
  export async function copyAsImage(value: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CopyAsImage', { value });
  }

  export async function isSmiles(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:IsSmiles', { s });
  }

  export async function isSmarts(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:IsSmarts', { s });
  }

  export async function detectSmiles(col: DG.Column, min: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DetectSmiles', { col, min });
  }

  export async function callChemSimilaritySearch(df: DG.DataFrame, col: DG.Column, molecule: string, metricName: string, limit: number, minScore: number, fingerprint: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CallChemSimilaritySearch', { df, col, molecule, metricName, limit, minScore, fingerprint });
  }

  export async function callChemDiversitySearch(col: DG.Column, metricName: string, limit: number, fingerprint: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CallChemDiversitySearch', { col, metricName, limit, fingerprint });
  }

  export async function addChemPropertiesColumns(table: DG.DataFrame, molecules: DG.Column, MW: boolean, HBA: boolean, HBD: boolean, logP: boolean, logS: boolean, PSA: boolean, rotatableBonds: boolean, stereoCenters: boolean, moleculeCharge: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:AddChemPropertiesColumns', { table, molecules, MW, HBA, HBD, logP, logS, PSA, rotatableBonds, stereoCenters, moleculeCharge });
  }

  export async function getMolProperty(property: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMolProperty', { property });
  }

  export async function addChemRisksColumns(table: DG.DataFrame, molecules: DG.Column, mutagenicity: boolean, tumorigenicity: boolean, irritatingEffects: boolean, reproductiveEffects: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:AddChemRisksColumns', { table, molecules, mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects });
  }

  //Generates a hierarchical tree based on the scaffolds presented in dataset
  export async function addScaffoldTree(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:AddScaffoldTree', {});
  }

  export async function mmpViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MmpViewer', {});
  }

  export async function mmpeditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MMPEditor', { call });
  }

  export async function mmpAnalysis(table: DG.DataFrame, molecules: DG.Column, activities: string[], diffTypes: any, fragmentCutoff: number): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MmpAnalysis', { table, molecules, activities, diffTypes, fragmentCutoff });
  }

  //Scaffold Tree filter
  export async function scaffoldTreeFilter(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ScaffoldTreeFilter', {});
  }

  export async function getScaffoldTree(data: DG.DataFrame, ringCutoff: number, dischargeAndDeradicalize: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetScaffoldTree', { data, ringCutoff, dischargeAndDeradicalize });
  }

  export async function removeDuplicates(molecules: any, molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:RemoveDuplicates', { molecules, molecule });
  }

  //Overview of Cheminformatics functionality
  export async function demoChemOverview(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoChemOverview', {});
  }

  //Searching for most similar or diverse molecules in dataset
  export async function demoSimilarityDiversitySearch(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoSimilarityDiversitySearch', {});
  }

  //Detect matched molecule pairs calculate the difference in activity values between them
  export async function demoMMPA(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoMMPA', {});
  }

  //R Group Analysis including R-group decomposition and  visual analysis of the obtained R-groups
  export async function demoRgroupAnalysis(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoRgroupAnalysis', {});
  }

  //Searching similar structures with significant activity difference
  export async function demoMoleculeActivityCliffs(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoMoleculeActivityCliffs', {});
  }

  //Maps the dataset to 2D plot based on similarity
  export async function demoChemicalSpace(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoChemicalSpace', {});
  }

  //Running scaffold analysis with hierarchical tree
  export async function demoScaffold(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:DemoScaffold', {});
  }

  export async function namesToSmiles(data: DG.DataFrame, names: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:NamesToSmiles', { data, names });
  }

  export async function canonicalize(molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Canonicalize', { molecule });
  }

  export async function getMolecularFormula(molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:GetMolecularFormula', { molecule });
  }

  export async function validateMolecule(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ValidateMolecule', { s });
  }

  //To be added
  export async function trainChemprop(df: DG.DataFrame, predictColumn: DG.Column, dataset_type: string, metric: string, multiclass_num_classes: number, num_folds: number, data_seed: number, split_sizes: any, split_type: string, activation: string, atom_messages: boolean, message_bias: boolean, ensemble_size: number, message_hidden_dim: number, depth: number, dropout: number, ffn_hidden_dim: number, ffn_num_layers: number, epochs: number, batch_size: number, warmup_epochs: number, init_lr: number, max_lr: number, final_lr: number, no_descriptor_scaling: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/chem:TrainChemprop', { df, predictColumn, dataset_type, metric, multiclass_num_classes, num_folds, data_seed, split_sizes, split_type, activation, atom_messages, message_bias, ensemble_size, message_hidden_dim, depth, dropout, ffn_hidden_dim, ffn_num_layers, epochs, batch_size, warmup_epochs, init_lr, max_lr, final_lr, no_descriptor_scaling });
  }

  export async function applyChemprop(df: DG.DataFrame, model: any): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ApplyChemprop', { df, model });
  }

  export async function isApplicableNN(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:IsApplicableNN', { df, predictColumn });
  }

  export async function isInteractiveNN(df: DG.DataFrame, predictColumn: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/chem:IsInteractiveNN', { df, predictColumn });
  }

  //Generates the new dataset based on the given structure
  export async function deprotect(table: DG.DataFrame, molecules: DG.Column, fragment: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Deprotect', { table, molecules, fragment });
  }

  //Beautifies the list of molecules and returns the list of beautified molecules
  export async function beautifyMols(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:BeautifyMols', {});
  }

  //Converts the list of molecules to V3K format using OCL
  export async function convertToV3KViaOCL(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:ConvertToV3KViaOCL', {});
  }

  //Calculates the MPO score for the column of molecules
  export async function mpo(): Promise<any> {
    return await grok.functions.call('@datagrok/chem:Mpo', {});
  }

  export async function mpoProfileEditor(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/chem:MpoProfileEditor', { file });
  }

  export async function checkJsonMpoProfile(content: string): Promise<any> {
    return await grok.functions.call('@datagrok/chem:CheckJsonMpoProfile', { content });
  }
}
