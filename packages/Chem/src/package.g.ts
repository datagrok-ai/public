import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//output: object module
export function getRdKitModule() : any {
  return PackageFunctions.getRdKitModule();
}

//input: string molString 
//output: object handler
export function getMolFileHandler(molString: string) : any {
  return PackageFunctions.getMolFileHandler(molString);
}

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//meta.role: autostart
export async function initChemAutostart() : Promise<void> {
  await PackageFunctions.initChemAutostart();
}

//name: Recalculate Coordinates
//description: Recalculates 2D coordinates for molecules in the column using Open Chem Lib
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: bool join = true 
//output: column result
//meta.role: transform
//top-menu: Chem | Transform | Recalculate Coordinates...
export async function recalculateCoordsViaOCL(table: DG.DataFrame, molecules: DG.Column, join: boolean) : Promise<any> {
  return await PackageFunctions.recalculateCoordsViaOCL(table, molecules, join);
}

//name: Chemistry | Most Diverse Structures
//input: column col { semType: Molecule }
//output: widget result
//meta.role: tooltip
export async function chemTooltip(col: DG.Column) : Promise<any> {
  return await PackageFunctions.chemTooltip(col);
}

//name: Scaffold Tree
//output: viewer result
//meta.icon: files/icons/scaffold-tree-icon.svg
//meta.role: viewer
export function scaffoldTreeViewer() : any {
  return PackageFunctions.scaffoldTreeViewer();
}

//name: Substructure Filter
//description: RDKit-based substructure filter
//output: filter result
//meta.semType: Molecule
//meta.primaryFilter: true
//meta.allowMultipleFiltersForColumn: false
//meta.role: filter
export function substructureFilter() : any {
  return PackageFunctions.substructureFilter();
}

//input: int x 
//input: int y 
//input: int w 
//input: int h 
//input: object canvas 
//input: string molString 
//input: string scaffoldMolString 
//input: object options { optional: true }
export function canvasMol(x: number, y: number, w: number, h: number, canvas: any, molString: string, scaffoldMolString: any, options: any) : void {
  PackageFunctions.canvasMol(x, y, w, h, canvas, molString, scaffoldMolString, options);
}

//input: string molStr 
//input: int w { optional: true }
//input: int h { optional: true }
//input: bool popupMenu { optional: true }
//output: object canvas
export function drawMolecule(molStr: string, w?: number, h?: number, popupMenu?: boolean) : any {
  return PackageFunctions.drawMolecule(molStr, w, h, popupMenu);
}

//input: string smiles { semType: Molecule }
//output: double result
export function getCLogP(smiles: string) : number {
  return PackageFunctions.getCLogP(smiles);
}

//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdKitCellRenderer() : Promise<any> {
  return await PackageFunctions.rdKitCellRenderer();
}

//name: chemCellRenderer
//output: grid_cell_renderer result
//meta.cellType: ChemicalReaction
//meta.role: cellRenderer
export async function rdKitReactionRenderer() : Promise<any> {
  return await PackageFunctions.rdKitReactionRenderer();
}

//name: chemMixtureRenderer
//output: grid_cell_renderer result
//meta.cellType: ChemicalMixture
//meta.role: cellRenderer
export async function rdKitMixtureRenderer() : Promise<any> {
  return await PackageFunctions.rdKitMixtureRenderer();
}

//output: grid_cell_renderer result
//meta.cellType: Molecule
//meta.role: cellRenderer
export async function chemCellRenderer() : Promise<any> {
  return await PackageFunctions.chemCellRenderer();
}

//input: column molColumn { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export async function getMorganFingerprints(molColumn: DG.Column) : Promise<any> {
  return await PackageFunctions.getMorganFingerprints(molColumn);
}

//input: string molString { semType: Molecule }
//output: object fingerprintBitset { description: Fingerprints }
export function getMorganFingerprint(molString: string) : any {
  return PackageFunctions.getMorganFingerprint(molString);
}

//input: column molStringsColumn 
//input: string molString 
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) : Promise<any> {
  return await PackageFunctions.getSimilarities(molStringsColumn, molString);
}

//input: column molStringsColumn 
//input: int limit 
//output: dataframe result
export async function getDiversities(molStringsColumn: DG.Column, limit: number) : Promise<any> {
  return await PackageFunctions.getDiversities(molStringsColumn, limit);
}

//input: column molStringsColumn 
//input: string molString 
//input: int limit 
//input: int cutoff 
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number, cutoff: number) : Promise<any> {
  return await PackageFunctions.findSimilar(molStringsColumn, molString, limit, cutoff);
}

//input: column molStringsColumn 
//input: string molString 
//input: string molBlockFailover 
//output: column result
export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, molBlockFailover: string) : Promise<any> {
  return await PackageFunctions.searchSubstructure(molStringsColumn, molString, molBlockFailover);
}

//description: As SDF...
//meta.role: fileExporter
export async function saveAsSdf() : Promise<void> {
  await PackageFunctions.saveAsSdf();
}

//name: Chem Similarity Search
//output: viewer result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
//meta.role: viewer
export function similaritySearchViewer() : any {
  return PackageFunctions.similaritySearchViewer();
}

//name: Similarity Search
//top-menu: Chem | Search | Similarity Search...
export function similaritySearchTopMenu() : void {
  PackageFunctions.similaritySearchTopMenu();
}

//name: Chem Diversity Search
//output: viewer result
//meta.icon: files/icons/chem-diversity-search-viewer.svg
//meta.role: viewer
export function diversitySearchViewer() : any {
  return PackageFunctions.diversitySearchViewer();
}

//name: Diversity Search
//top-menu: Chem | Search | Diversity Search...
export function diversitySearchTopMenu() : void {
  PackageFunctions.diversitySearchTopMenu();
}

//top-menu: Chem | Calculate | Descriptors...
export async function descriptorsDocker() : Promise<void> {
  await PackageFunctions.descriptorsDocker();
}

//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: list<string> selected 
//meta.role: transform
export async function calculateDescriptorsTransform(table: DG.DataFrame, molecules: DG.Column, selected: string[]) : Promise<void> {
  await PackageFunctions.calculateDescriptorsTransform(table, molecules, selected);
}

//input: column molecules { semType: Molecule }
//input: list<string> selected { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function getDescriptors(molecules: DG.Column, selected?: string[]) : Promise<any> {
  return await PackageFunctions.getDescriptors(molecules, selected);
}

//output: object descriptors
export async function chemDescriptorsTree() : Promise<any> {
  return await PackageFunctions.chemDescriptorsTree();
}

//name: Map Identifiers
//top-menu: Chem | Calculate | Map Identifiers...
export async function getMapIdentifiers() : Promise<void> {
  await PackageFunctions.getMapIdentifiers();
}

//input: string molfile 
//output: string result
export async function freeTextToSmiles(molfile: string) : Promise<any> {
  return await PackageFunctions.freeTextToSmiles(molfile);
}

//input: dataframe table 
//input: column molecules 
//input: list<string> descriptors 
export async function chemDescriptors(table: DG.DataFrame, molecules: DG.Column, descriptors: string[]) : Promise<void> {
  await PackageFunctions.chemDescriptors(table, molecules, descriptors);
}

//name: SearchSubstructureEditor
//input: funccall call 
//meta.role: editor
export function searchSubstructureEditor(call: DG.FuncCall) : void {
  PackageFunctions.searchSubstructureEditor(call);
}

//name: Substructure Search
//input: column molecules { semType: Molecule }
//top-menu: Chem | Search | Substructure Search...
//editor: Chem:SearchSubstructureEditor
export function SubstructureSearchTopMenu(molecules: DG.Column) : void {
  PackageFunctions.SubstructureSearchTopMenu(molecules);
}

//name: Cluster MCS
//description: Calculates most common substructures for each cluster
//input: dataframe table 
//input: column molCol { semType: Molecule }
//input: column clusterCol { type: categorical }
//top-menu: Chem | Calculate | Cluster MCS...
export async function clusterMCSTopMenu(table: DG.DataFrame, molCol: DG.Column, clusterCol: DG.Column) : Promise<void> {
  await PackageFunctions.clusterMCSTopMenu(table, molCol, clusterCol);
}

//name: clusterMCS
//description: Calculates most common substructures for each cluster
//input: column molCol { semType: Molecule }
//input: string clusterCol { type: string }
//output: column result { semType: Molecule }
//meta.vectorFunc: true
//friendlyName: Cluster MCS
export async function performClusterMCS(molCol: DG.Column, clusterCol: DG.Column) : Promise<any> {
  return await PackageFunctions.performClusterMCS(molCol, clusterCol);
}

//input: funccall call 
//meta.role: editor
export function ChemSpaceEditor(call: DG.FuncCall) : void {
  PackageFunctions.ChemSpaceEditor(call);
}

//name: Fingerprints
//input: column col { semType: Molecule }
//input: string _metric { optional: true }
//input: string fingerprintType = 'Morgan' { caption: Fingerprint type; optional: true; choices: ["Morgan","RDKit","Pattern","AtomPair","MACCS","TopologicalTorsion"] }
//output: object result
//meta.supportedSemTypes: Molecule
//meta.supportedDistanceFunctions: Tanimoto,Asymmetric,Cosine,Sokal
//meta.role: dimRedPreprocessingFunction
export async function getFingerprints(col: DG.Column, _metric: string | undefined, fingerprintType: any) {
  return await PackageFunctions.getFingerprints(col, _metric, fingerprintType);
}

//name: Chem Space
//description: Maps the dataset to 2D plot based on similarity
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string methodName { choices: ["UMAP","t-SNE"] }
//input: string similarityMetric { choices: ["Tanimoto","Asymmetric","Cosine","Sokal"] }
//input: bool plotEmbeddings = true 
//input: object options { optional: true }
//input: func preprocessingFunction { optional: true }
//input: bool clusterEmbeddings { optional: true }
//input: bool clusterMCS { optional: true }
//output: viewer result
//top-menu: Chem | Analyze | Chemical Space...
//editor: Chem:ChemSpaceEditor
export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options?: any, preprocessingFunction?: any, clusterEmbeddings?: boolean, clusterMCS?: boolean) : Promise<any> {
  return await PackageFunctions.chemSpaceTopMenu(table, molecules, methodName, similarityMetric, plotEmbeddings, options, preprocessingFunction, clusterEmbeddings, clusterMCS);
}

//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string methodName 
//input: string similarityMetric 
//input: bool plotEmbeddings = true 
//input: string options { optional: true }
//input: bool clusterEmbeddings { optional: true }
//output: viewer result
//meta.role: transform
export async function chemSpaceTransform(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options?: string, clusterEmbeddings?: boolean) : Promise<any> {
  return await PackageFunctions.chemSpaceTransform(table, molecules, methodName, similarityMetric, plotEmbeddings, options, clusterEmbeddings);
}

//name: Chem Space Embeddings
//input: string col 
//input: string methodName 
//input: string similarityMetric 
//input: string xAxis 
//input: string yAxis 
//input: object options { optional: true }
//output: object result
export async function getChemSpaceEmbeddings(col: DG.Column, methodName: any, similarityMetric: any, xAxis: string, yAxis: string, options?: any) : Promise<any> {
  return await PackageFunctions.getChemSpaceEmbeddings(col, methodName, similarityMetric, xAxis, yAxis, options);
}

//name: Chem Similarities Matrix
//input: int dim 
//input: column col 
//input: dataframe df 
//input: string colName 
//input: object simArr 
//output: object result
export async function getChemSimilaritiesMatrix(dim: number, col: DG.Column, df: DG.DataFrame, colName: string, simArr: DG.Column[]) : Promise<any> {
  return await PackageFunctions.getChemSimilaritiesMatrix(dim, col, df, colName, simArr);
}

//name: Elemental Analysis
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: bool radarViewer = false { description: Add a standalone radar viewer }
//input: bool radarGrid = false { description: Show radar in grid cells }
//top-menu: Chem | Analyze | Elemental Analysis...
export async function elementalAnalysis(table: DG.DataFrame, molecules: DG.Column, radarViewer: boolean, radarGrid: boolean) : Promise<void> {
  await PackageFunctions.elementalAnalysis(table, molecules, radarViewer, radarGrid);
}

//input: dataframe table 
//input: column molecules { semType: Molecule }
//output: list res
//meta.role: transform
export function runElementalAnalysis(table: DG.DataFrame, molecules: DG.Column) : string[] {
  return PackageFunctions.runElementalAnalysis(table, molecules);
}

//name: R-Groups Analysis
//top-menu: Chem | Analyze | R-Groups Analysis...
export function rGroupsAnalysisMenu() : void {
  PackageFunctions.rGroupsAnalysisMenu();
}

//input: dataframe df 
//input: string molColName 
//input: string core 
//input: string rGroupName 
//input: string rGroupMatchingStrategy 
//input: bool onlyMatchAtRGroups = false { optional: true }
//output: object result
//meta.role: transform
export async function rGroupDecomposition(df: DG.DataFrame, molColName: string, core: string, rGroupName: string, rGroupMatchingStrategy: string, onlyMatchAtRGroups: boolean) : Promise<any> {
  return await PackageFunctions.rGroupDecomposition(df, molColName, core, rGroupName, rGroupMatchingStrategy, onlyMatchAtRGroups);
}

//input: funccall call 
//meta.role: editor
export function ActivityCliffsEditor(call: DG.FuncCall) : void {
  PackageFunctions.ActivityCliffsEditor(call);
}

//name: Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table { description: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity = 80 { description: Similarity cutoff }
//input: string methodName { choices: ["UMAP","t-SNE"] }
//input: string similarityMetric { choices: ["Tanimoto","Asymmetric","Cosine","Sokal"] }
//input: func preprocessingFunction { optional: true }
//input: object options { optional: true }
//input: bool isDemo { optional: true }
//input: bool isTest { optional: true }
//top-menu: Chem | Analyze | Activity Cliffs...
//editor: Chem:ActivityCliffsEditor
export async function activityCliffs(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, preprocessingFunction: any, options?: any, isDemo?: boolean, isTest?: boolean) : Promise<void> {
  await PackageFunctions.activityCliffs(table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, isDemo, isTest);
}

//input: viewer sp 
export async function activityCliffsInitFunction(sp: any) : Promise<void> {
  await PackageFunctions.activityCliffsInitFunction(sp);
}

//input: dataframe table { description: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity = 80 { description: Similarity cutoff }
//input: string methodName { choices: ["UMAP","t-SNE"] }
//input: string similarityMetric { choices: ["Tanimoto","Asymmetric","Cosine","Sokal"] }
//input: string options { optional: true }
//input: bool isDemo { optional: true }
//meta.role: transform
export async function activityCliffsTransform(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, options?: string, isDemo?: boolean) : Promise<void> {
  await PackageFunctions.activityCliffsTransform(table, molecules, activities, similarity, methodName, similarityMetric, options, isDemo);
}

//name: To InchI
//input: dataframe table 
//input: column molecules { semType: Molecule }
//meta.role: transform
//top-menu: Chem | Calculate | To InchI...
export function addInchisTopMenu(table: DG.DataFrame, col: DG.Column) : void {
  PackageFunctions.addInchisTopMenu(table, col);
}

//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchis(molecules: DG.Column) : any {
  return PackageFunctions.getInchis(molecules);
}

//name: To InchI Keys
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//meta.role: transform
//top-menu: Chem | Calculate | To InchI Keys...
export function addInchisKeysTopMenu(table: DG.DataFrame, col: DG.Column) : void {
  PackageFunctions.addInchisKeysTopMenu(table, col);
}

//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchiKeys(molecules: DG.Column) : any {
  return PackageFunctions.getInchiKeys(molecules);
}

//name: Structural Alerts
//description: Highlights the fragments that could lead to potential chemical hazards
//input: dataframe table { description: Input data table; caption: Table }
//input: column molecules { caption: Molecules; semType: Molecule; type: categorical }
//input: bool pains = true { caption: PAINS; description: "Pan Assay Interference Compounds filters" }
//input: bool bms = false { caption: BMS; description: "Bristol-Myers Squibb HTS Deck filters" }
//input: bool sureChembl = false { caption: SureChEMBL; description: "MedChem unfriendly compounds from SureChEMBL" }
//input: bool mlsmr = false { caption: MLSMR; description: "NIH MLSMR Excluded Functionality filters" }
//input: bool dundee = false { caption: Dundee; description: "University of Dundee NTD Screening Library filters" }
//input: bool inpharmatica = false { caption: Inpharmatica; description: "Inpharmatica filters" }
//input: bool lint = false { caption: LINT; description: "Pfizer LINT filters" }
//input: bool glaxo = false { caption: Glaxo; description: "Glaxo Wellcome Hard filters" }
//output: dataframe result
//meta.role: hitTriageFunction
//top-menu: Chem | Analyze | Structural Alerts...
export async function structuralAlertsTopMenu(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) : Promise<any> {
  return await PackageFunctions.structuralAlertsTopMenu(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//input: dataframe table { caption: Table; description: Input data table }
//input: column molecules { caption: Molecules; type: categorical; semType: Molecule }
//input: bool pains = true { caption: PAINS; description: "Pan Assay Interference Compounds filters" }
//input: bool bms = false { caption: BMS; description: "Bristol-Myers Squibb HTS Deck filters" }
//input: bool sureChembl = false { caption: SureChEMBL; description: "MedChem unfriendly compounds from SureChEMBL" }
//input: bool mlsmr = false { caption: MLSMR; description: "NIH MLSMR Excluded Functionality filters" }
//input: bool dundee = false { caption: Dundee; description: "University of Dundee NTD Screening Library filters" }
//input: bool inpharmatica = false { caption: Inpharmatica; description: "Inpharmatica filters" }
//input: bool lint = false { caption: LINT; description: "Pfizer LINT filters" }
//input: bool glaxo = false { caption: Glaxo; description: "Glaxo Wellcome Hard filters" }
//output: dataframe result
//meta.role: transform
export async function runStructuralAlerts(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) {
  return await PackageFunctions.runStructuralAlerts(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//input: column<string> molecules { semType: Molecule }
//input: list<string> alerts { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function getStructuralAlerts(molecules: DG.Column, alerts?: string[]) : Promise<any> {
  return await PackageFunctions.getStructuralAlerts(molecules, alerts);
}

//name: Chemistry | Rendering
//input: column molColumn { semType: Molecule }
//output: widget result
//meta.exclude-actions-panel: true
//meta.role: panel
export function molColumnPropertyPanel(molColumn: DG.Column) : any {
  return PackageFunctions.molColumnPropertyPanel(molColumn);
}

//name: Chemistry | Highlight
//input: column molColumn { semType: Molecule }
//output: widget result
//meta.exclude-actions-panel: true
//meta.role: panel
export function molColumnHighlights(molColumn: DG.Column) : any {
  return PackageFunctions.molColumnHighlights(molColumn);
}

//name: Chemistry | Descriptors
//input: string smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export function descriptorsWidget(smiles: string) : any {
  return PackageFunctions.descriptorsWidget(smiles);
}

//name: Biology | Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//help-url: /help/domains/chem/info-panels/drug-likeness.md
export function drugLikeness(smiles: DG.SemanticValue) : any {
  return PackageFunctions.drugLikeness(smiles);
}

//name: Chemistry | Properties
//description: Basic molecule properties
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export function properties(smiles: DG.SemanticValue) : any {
  return PackageFunctions.properties(smiles);
}

//description: Return chem property function
//input: string name 
//output: object result
export function getChemPropertyFunction(name: string) : any {
  return PackageFunctions.getChemPropertyFunction(name);
}

//name: Biology | Structural Alerts
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
//input: string smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//help-url: /help/domains/chem/info-panels/structural-alerts.md
export async function structuralAlerts(smiles: string) : Promise<any> {
  return await PackageFunctions.structuralAlerts(smiles);
}

//name: Structure | Identifiers
//input: string smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function identifiers(smiles: string) : Promise<any> {
  return await PackageFunctions.identifiers(smiles);
}

//name: Structure | 3D Structure
//description: 3D molecule representation
//input: string molecule { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function structure3D(molecule: string) : Promise<any> {
  return await PackageFunctions.structure3D(molecule);
}

//name: Structure | 2D Structure
//description: 2D molecule representation
//input: string molecule { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export function structure2d(molecule: string) : any {
  return PackageFunctions.structure2d(molecule);
}

//name: Biology | Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
export function toxicity(smiles: DG.SemanticValue) : any {
  return PackageFunctions.toxicity(smiles);
}

//input: column molecule { semType: Molecule }
//input: string targetNotation 
//input: bool kekulize = false { optional: true; nullable: true }
//output: column result
//meta.vectorFunc: true
export async function convertMoleculeNotation(molecule: DG.Column, targetNotation: any, kekulize?: boolean) : Promise<any> {
  return await PackageFunctions.convertMoleculeNotation(molecule, targetNotation, kekulize);
}

//description: RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000
//input: string molecule { semType: Molecule }
//input: string sourceNotation { choices: ["smiles","cxsmiles","smarts","cxsmarts","molblock","v3Kmolblock"] }
//input: string targetNotation { choices: ["smiles","cxsmiles","smarts","cxsmarts","molblock","v3Kmolblock"] }
//output: string result { semType: Molecule }
//meta.role: unitConverter
export function convertMolNotation(molecule: string, sourceNotation: any, targetNotation: any) : string {
  return PackageFunctions.convertMolNotation(molecule, sourceNotation, targetNotation);
}

//name: Convert Notation
//input: dataframe data 
//input: column molecules { semType: Molecule }
//input: string targetNotation = 'smiles' { choices: ["smiles","cxsmiles","smarts","cxsmarts","molblock","v3Kmolblock"] }
//input: bool overwrite = false 
//input: bool join = true 
//input: bool kekulize = false { optional: true; nullable: true }
//output: column result
//meta.role: transform
//top-menu: Chem | Transform | Convert Notation...
export async function convertNotation(data: DG.DataFrame, molecules: DG.Column<any>, targetNotation: any, overwrite: boolean, join: boolean, kekulize?: boolean) : Promise<any> {
  return await PackageFunctions.convertNotation(data, molecules, targetNotation, overwrite, join, kekulize);
}

//name: Convert Notation...
//input: column col { semType: Molecule }
//meta.action: Convert Notation...
export function convertMolNotationAction(col: DG.Column) : void {
  PackageFunctions.convertMolNotationAction(col);
}

//name: Convert Mixture To Smiles...
//input: column col { semType: ChemicalMixture }
//meta.action: Convert mixture to smiles...
export function convertMixtureToSmiles(col: DG.Column) : void {
  PackageFunctions.convertMixtureToSmiles(col);
}

//description: Molecule
//input: grid_cell cell 
//meta.role: cellEditor
export async function editMoleculeCell(cell: any) : Promise<void> {
  await PackageFunctions.editMoleculeCell(cell);
}

//name: OpenChemLib
//output: widget sketcher
//meta.role: moleculeSketcher
export function openChemLibSketcher() : any {
  return PackageFunctions.openChemLibSketcher();
}

//description: Opens SDF file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: sdf,mol
export function importSdf(bytes: Uint8Array) : any {
  return PackageFunctions.importSdf(bytes);
}

//description: Opens smi file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: smi
export function importSmi(bytes: Uint8Array) : any {
  return PackageFunctions.importSmi(bytes);
}

//description: Opens smi file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mol2
export function importMol2(bytes: Uint8Array) : any {
  return PackageFunctions.importMol2(bytes);
}

//description: Opens MOL file
//input: string content 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mol
export function importMol(content: string) : any {
  return PackageFunctions.importMol(content);
}

//output: grid_cell_renderer result
//meta.chemRendererName: OpenChemLib
export async function oclCellRenderer() : Promise<any> {
  return await PackageFunctions.oclCellRenderer();
}

//name: Sort by similarity
//description: Sorts a molecular column by similarity
//input: semantic_value value { semType: Molecule }
//meta.action: Sort by similarity
export async function sortBySimilarity(value: DG.SemanticValue) : Promise<void> {
  await PackageFunctions.sortBySimilarity(value);
}

//name: Use as filter
//description: Adds this structure as a substructure filter
//input: semantic_value value { semType: Molecule }
//meta.action: Use as filter
export function useAsSubstructureFilter(value: DG.SemanticValue) : void {
  PackageFunctions.useAsSubstructureFilter(value);
}

//name: Copy as...
//description: Copies structure in different formats
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as...
//meta.exclude-current-value-menu: true
export function copyAsAction(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsAction(value);
}

//name: Copy as SMILES
//description: Copies structure as smiles
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as SMILES
//meta.exclude-actions-panel: true
export function copyAsSmiles(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsSmiles(value);
}

//name: Copy as MOLFILE V2000
//description: Copies structure as molfile V2000
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as MOLFILE V2000
//meta.exclude-actions-panel: true
export function copyAsMolfileV2000(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsMolfileV2000(value);
}

//name: Copy as MOLFILE V3000
//description: Copies structure as molfile V3000
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as MOLFILE V3000
//meta.exclude-actions-panel: true
export function copyAsMolfileV3000(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsMolfileV3000(value);
}

//name: Copy as SMARTS
//description: Copies structure as smarts
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as SMARTS
//meta.exclude-actions-panel: true
export function copyAsSmarts(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsSmarts(value);
}

//name: Copy as IMAGE
//description: Copies structure as Image
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as Image
//meta.exclude-actions-panel: true
export function copyAsImage(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsImage(value);
}

//input: string s 
//output: bool result
export function isSmiles(s: string) : boolean {
  return PackageFunctions.isSmiles(s);
}

//input: string s 
//output: bool result
export function isSmarts(s: string) : boolean {
  return PackageFunctions.isSmarts(s);
}

//input: column col 
//input: int min 
export function detectSmiles(col: DG.Column, min: number) : void {
  PackageFunctions.detectSmiles(col, min);
}

//name: chemSimilaritySearch
//input: dataframe df 
//input: column col 
//input: string molecule 
//input: string metricName 
//input: string fingerprint 
//input: int limit 
//input: double minScore 
//output: dataframe result
export async function callChemSimilaritySearch(df: DG.DataFrame, col: DG.Column, molecule: string, metricName: any, fingerprint: string, limit: number, minScore: number) : Promise<any> {
  return await PackageFunctions.callChemSimilaritySearch(df, col, molecule, metricName, fingerprint, limit, minScore);
}

//name: chemDiversitySearch
//input: column col 
//input: string metricName 
//input: string fingerprint 
//input: int limit 
//output: dataframe result
export async function callChemDiversitySearch(col: DG.Column, metricName: any, fingerprint: string, limit: number) : Promise<number[]> {
  return await PackageFunctions.callChemDiversitySearch(col, metricName, fingerprint, limit);
}

//name: Chemical Properties
//description: Calculates chemical properties and adds them as columns to the input table. properties include Molecular Weight (MW), Hydrogen Bond Acceptors (HBA), Hydrogen Bond Donors (HBD), LogP (Partition), LogS (Solubility), Polar Surface Area (PSA), Rotatable Bonds, Stereo Centers, Molecule Charge.
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: bool MW = true 
//input: bool HBA = false 
//input: bool HBD = false 
//input: bool logP = false 
//input: bool logS = false 
//input: bool PSA = false 
//input: bool rotatableBonds = false 
//input: bool stereoCenters = false 
//input: bool moleculeCharge = false 
//meta.function_family: biochem-calculator
//meta.method_info.author: Open Chem Lib Team
//meta.method_info.year: 2024
//meta.method_info.github: https://github.com/actelion/openchemlib
//meta.role: hitTriageFunction,transform
//top-menu: Chem | Calculate | Chemical Properties...
export async function addChemPropertiesColumns(table: DG.DataFrame, molecules: DG.Column, MW?: boolean, HBA?: boolean, HBD?: boolean, logP?: boolean, logS?: boolean, PSA?: boolean, rotatableBonds?: boolean, stereoCenters?: boolean, moleculeCharge?: boolean) : Promise<void> {
  await PackageFunctions.addChemPropertiesColumns(table, molecules, MW, HBA, HBD, logP, logS, PSA, rotatableBonds, stereoCenters, moleculeCharge);
}

//input: column molecules { semType: Molecule }
//input: list<string> selected { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function getProperties(molecules: DG.Column, selected?: string[]) : Promise<any> {
  return await PackageFunctions.getProperties(molecules, selected);
}

//name: Toxicity Risks
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: bool mutagenicity = true 
//input: bool tumorigenicity = false 
//input: bool irritatingEffects = false 
//input: bool reproductiveEffects = false 
//meta.role: hitTriageFunction,transform
//top-menu: Chem | Calculate | Toxicity Risks...
export async function addChemRisksColumns(table: DG.DataFrame, molecules: DG.Column, mutagenicity?: boolean, tumorigenicity?: boolean, irritatingEffects?: boolean, reproductiveEffects?: boolean) : Promise<void> {
  await PackageFunctions.addChemRisksColumns(table, molecules, mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects);
}

//input: column molecules { semType: Molecule }
//input: list<string> risks { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function getToxicityRisks(molecules: DG.Column, risks?: string[]) : Promise<any> {
  return await PackageFunctions.getToxicityRisks(molecules, risks);
}

//description: Generates a hierarchical tree based on the scaffolds presented in dataset
//top-menu: Chem | Analyze | Scaffold Tree
export function addScaffoldTree() : void {
  PackageFunctions.addScaffoldTree();
}

//name: Matched Molecular Pairs Analysis
//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function mmpViewer() : any {
  return PackageFunctions.mmpViewer();
}

//input: funccall call 
//meta.role: editor
export function MMPEditor(call: DG.FuncCall) : void {
  PackageFunctions.MMPEditor(call);
}

//name: Matched Molecular Pairs
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: column_list activities { type: numerical }
//input: string_list diffTypes 
//input: string_list scalings 
//input: double fragmentCutoff = 0.4 { description: Maximum fragment size relative to core }
//editor: Chem:MMPEditor
//top-menu: Chem | Analyze | Matched Molecular Pairs...
export async function mmpAnalysis(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column[], diffTypes: any, scalings: any, fragmentCutoff: number) : Promise<void> {
  await PackageFunctions.mmpAnalysis(table, molecules, activities, diffTypes, scalings, fragmentCutoff);
}

//name: Scaffold Tree Filter
//description: Scaffold Tree filter
//output: filter result
//meta.semType: Molecule
//meta.allowMultipleFiltersForColumn: false
//meta.role: filter
export function scaffoldTreeFilter() : any {
  return PackageFunctions.scaffoldTreeFilter();
}

//input: dataframe data 
//input: int ringCutoff = 10 { description: Ignore molecules with # rings > N }
//input: bool dischargeAndDeradicalize = false { description: Remove charges and radicals from scaffolds }
//output: string result
export async function getScaffoldTree(data: DG.DataFrame, ringCutoff: number, dischargeAndDeradicalize: boolean) : Promise<string> {
  return await PackageFunctions.getScaffoldTree(data, ringCutoff, dischargeAndDeradicalize);
}

//name: filterMoleculeDuplicates
//input: list<string> molecules 
//input: string molecule 
//output: list<string> result
export function removeDuplicates(molecules: string[], molecule: string) : string[] {
  return PackageFunctions.removeDuplicates(molecules, molecule);
}

//name: Demo Chem Overview
//description: Overview of Cheminformatics functionality
//meta.isDemoScript: true
//meta.demoSkip: GROK-14320
//meta.demoPath: Cheminformatics | Overview
export async function demoChemOverview() : Promise<void> {
  await PackageFunctions.demoChemOverview();
}

//name: Demo Similarity Search
//description: Searching for most similar or diverse molecules in dataset
//meta.demoPath: Cheminformatics | Similarity & Diversity Search
export async function demoSimilarityDiversitySearch() : Promise<void> {
  await PackageFunctions.demoSimilarityDiversitySearch();
}

//name: Demo Matched Molecular Pairs
//description: Detect matched molecule pairs calculate the difference in activity values between them
//meta.demoPath: Cheminformatics | Matched Molecular Pairs
export async function demoMMPA() : Promise<void> {
  await PackageFunctions.demoMMPA();
}

//name: Demo R Group Analysis
//description: R Group Analysis including R-group decomposition and  visual analysis of the obtained R-groups
//meta.demoPath: Cheminformatics | R-Group Analysis
//meta.demoSkip: GROK-14320
//meta.isDemoDashboard: true
export async function demoRgroupAnalysis() : Promise<void> {
  await PackageFunctions.demoRgroupAnalysis();
}

//name: Demo Activity Cliffs
//description: Searching similar structures with significant activity difference
//meta.demoPath: Cheminformatics | Molecule Activity Cliffs
//meta.demoSkip: GROK-14320
//meta.isDemoDashboard: true
export async function demoMoleculeActivityCliffs() : Promise<void> {
  await PackageFunctions.demoMoleculeActivityCliffs();
}

//name: Demo Chemical Space
//description: Maps the dataset to 2D plot based on similarity
//meta.demoPath: Cheminformatics | Chemical Space
//meta.demoSkip: GROK-14320
export async function demoChemicalSpace() : Promise<void> {
  await PackageFunctions.demoChemicalSpace();
}

//name: Demo Scaffold Tree
//description: Running scaffold analysis with hierarchical tree
//meta.demoPath: Cheminformatics | Scaffold Tree
export async function demoScaffold() : Promise<void> {
  await PackageFunctions.demoScaffold();
}

//name: Names To Smiles
//input: dataframe data 
//input: column names 
//meta.role: transform
//top-menu: Chem | Transform | Names To Smiles...
export async function namesToSmiles(data: DG.DataFrame, names: DG.Column<any>) : Promise<void> {
  await PackageFunctions.namesToSmiles(data, names);
}

//input: string molecule { semType: Molecule }
//output: string smiles { semType: Molecule }
//meta.role: canonicalizer
export function canonicalize(molecule: string) : string {
  return PackageFunctions.canonicalize(molecule);
}

//input: string molecule { semType: Molecule }
//output: string molecularFormula
export function getMolecularFormula(molecule: string) : string {
  return PackageFunctions.getMolecularFormula(molecule);
}

//input: string s 
//output: object result
export function validateMolecule(s: string) : any {
  return PackageFunctions.validateMolecule(s);
}

//description: To be added
//input: dataframe df 
//input: column predictColumn 
//input: string dataset_type = 'regression' { category: General; choices: ["regression","classification"]; description: Type of dataset,e.g. classification or regression. This determines the loss function used during training. }
//input: string metric = 'rmse' { category: General; choices: ["mse","mae","rmse","bounded-mse","bounded-mae","bounded-rmse","r2","binary-mcc","multiclass-mcc","roc","prc","accuracy","f1"]; description: Metric to use during evaluation. Note:Does NOT affect loss function used during training (loss is determined by the `dataset_type` argument). }
//input: int multiclass_num_classes = 3 { category: General; description: Number of classes when running multiclass classification }
//input: int num_folds = 1 { category: General; description: Number of folds when performing cross validation }
//input: int data_seed = 0 { category: General; description: Random seed to use when splitting data into train/val/test sets. When `num_folds` > 1,the first fold uses this seed and all subsequent folds add 1 to the seed. }
//input: list split_sizes = [0.8,0.1,0.1] { category: General; description: Split proportions for train/validation/test sets }
//input: string split_type = 'random' { category: General; choices: ["random","scaffold_balanced","cv","cv_no_val","kennard_stone","kmeans","random_with_repeated_smiles"]; description: Method of splitting the data into train/val/test }
//input: string activation = 'ReLU' { category: Model; choices: ["ReLU","LeakyReLU","PReLU","tanh","SELU","ELU"]; description: Activation function }
//input: bool atom_messages = false { category: Model; description: Use messages on atoms instead of messages on bonds }
//input: bool message_bias = false { category: Model; description: Whether to add bias to linear layers }
//input: int ensemble_size = 1 { category: Model; description: Number of models in ensemble }
//input: int message_hidden_dim = 300 { category: Model; description: Dimensionality of hidden layers in MPN }
//input: int depth = 3 { category: Model; description: Number of message passing step }
//input: double dropout = 0.0 { category: Model; description: Dropout probability }
//input: int ffn_hidden_dim = 300 { category: Model; description: Hidden dim for higher-capacity FFN (defaults to hidden_size) }
//input: int ffn_num_layers = 2 { category: Model; description: Number of layers in FFN after MPN encoding }
//input: int epochs = 50 { category: Training; description: Number of epochs to run }
//input: int batch_size = 64 { category: Training; description: Batch size }
//input: double warmup_epochs = 2.0 { category: Training; description: Number of epochs during which learning rate increases linearly from init_lr to max_lr. Afterwards,learning rate decreases exponentially from max_lr to final_lr. }
//input: double init_lr = 0.001 { category: Training; description: Initial learning rate }
//input: double max_lr = 0.001 { category: Training; description: Maximum learning rate }
//input: double final_lr = 0.0001 { category: Training; description: Final learning rate }
//input: bool no_descriptor_scaling = false { category: Training; description: Turn off scaling of features }
//output: dynamic model
//meta.mlname: Chemprop
//meta.mlrole: train
export async function trainChemprop(df: DG.DataFrame, predictColumn: DG.Column, dataset_type: string, metric: string, multiclass_num_classes: number, num_folds: number, data_seed: number, split_sizes: any, split_type: string, activation: string, atom_messages: boolean, message_bias: boolean, ensemble_size: number, message_hidden_dim: number, depth: number, dropout: number, ffn_hidden_dim: number, ffn_num_layers: number, epochs: number, batch_size: number, warmup_epochs: number, init_lr: number, max_lr: number, final_lr: number, no_descriptor_scaling: boolean) : Promise<any> {
  return await PackageFunctions.trainChemprop(df, predictColumn, dataset_type, metric, multiclass_num_classes, num_folds, data_seed, split_sizes, split_type, activation, atom_messages, message_bias, ensemble_size, message_hidden_dim, depth, dropout, ffn_hidden_dim, ffn_num_layers, epochs, batch_size, warmup_epochs, init_lr, max_lr, final_lr, no_descriptor_scaling);
}

//input: dataframe df 
//input: dynamic model 
//output: dataframe data_out
//meta.mlname: Chemprop
//meta.mlrole: apply
export async function applyChemprop(df: DG.DataFrame, model: Uint8Array) {
  return await PackageFunctions.applyChemprop(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Chemprop
//meta.mlrole: isApplicable
export async function isApplicableNN(df: DG.DataFrame, predictColumn: DG.Column) {
  return await PackageFunctions.isApplicableNN(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Chemprop
//meta.mlrole: isInteractive
//meta.mlupdate: false
export async function isInteractiveNN(df: DG.DataFrame, predictColumn: DG.Column) {
  return await PackageFunctions.isInteractiveNN(df, predictColumn);
}

//name: Deprotect
//description: Removes drawn protecting groups / fragments from molecules
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: string fragment = 'O=C([N:1])OCC1c2ccccc2-c2ccccc21' { semType: Molecule }
//meta.role: transform
//top-menu: Chem | Transform | Deprotect...
//editor: Chem:DeprotectEditor
export async function deprotect(table: DG.DataFrame, molecules: DG.Column, fragment: string) : Promise<void> {
  await PackageFunctions.deprotect(table, molecules, fragment);
}

//name: Deprotect Editor
//input: funccall call 
//output: widget result
//meta.role: editor
export function deprotectEditor(call: DG.FuncCall) : any {
  return PackageFunctions.deprotectEditor(call);
}

//description: Beautifies the list of molecules and returns the list of beautified molecules
//input: list<string> mols 
//output: list<string> result
export async function beautifyMols(mols: string[]) : Promise<string[]> {
  return await PackageFunctions.beautifyMols(mols);
}

//description: Converts the list of molecules to V3K format using OCL
//input: list<string> mols 
//output: list<string> result
export async function convertToV3KViaOCL(mols: string[]) : Promise<string[]> {
  return await PackageFunctions.convertToV3KViaOCL(mols);
}

//name: mpo
//description: Calculates the MPO score for the column of molecules
//top-menu: Chem | Calculate | MPO Score...
export async function _mpo() : Promise<void> {
  await PackageFunctions._mpo();
}

//input: dataframe df 
//input: string profileName 
//input: dynamic aggregation 
//input: object currentProperties 
//output: list res
//meta.role: transform
export async function mpoTransformFunction(df: DG.DataFrame, profileName: string, aggregation: any, currentProperties: any) : Promise<string[]> {
  return await PackageFunctions.mpoTransformFunction(df, profileName, aggregation, currentProperties);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: Chem:checkJsonMpoProfile
export function mpoProfileEditor(file: DG.FileInfo) : any {
  return PackageFunctions.mpoProfileEditor(file);
}

//input: string content 
//output: bool result
export function checkJsonMpoProfile(content: string) {
  return PackageFunctions.checkJsonMpoProfile(content);
}

//name: Chemistry | Mixture
//input: string mixture { semType: ChemicalMixture }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function mixtureWidget(mixture: string) : Promise<any> {
  return await PackageFunctions.mixtureWidget(mixture);
}

//name: Chemistry | MixtureTree
//input: string mixture { semType: ChemicalMixture }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function mixtureTreeWidget(mixture: string) : Promise<any> {
  return await PackageFunctions.mixtureTreeWidget(mixture);
}

//name: Biochemical Properties
//description: Dynamically discovers and executes tagged biochemical calculators
//top-menu: Chem | Calculate | Biochemical Properties
export async function biochemPropsWidget() : Promise<void> {
  await PackageFunctions.biochemPropsWidget();
}

//name: MPO profiles
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.browsePath: Chem
//meta.icon: images/mpo.png
//meta.role: app
export async function mpoProfilesApp(path?: string) : Promise<any> {
  return await PackageFunctions.mpoProfilesApp(path);
}

//input: dynamic treeNode 
//input: view browseView 
export async function mpoProfilesAppTreeBrowser(treeNode: any, browseView: any) : Promise<void> {
  await PackageFunctions.mpoProfilesAppTreeBrowser(treeNode, browseView);
}

//name: Chemistry | MPO
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function mpoWidget(semValue: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.mpoWidget(semValue);
}
