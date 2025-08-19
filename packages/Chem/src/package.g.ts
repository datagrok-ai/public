import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getRdKitModule
//output: object module
export function getRdKitModule() : any {
  return PackageFunctions.getRdKitModule();
}

//name: getMolFileHandler
//input: string molString 
//output: object handler
export function getMolFileHandler(molString: string) : any {
  return PackageFunctions.getMolFileHandler(molString);
}

//tags: init
export async function init() : Promise<void> {
  PackageFunctions.init();
}

//name: initChemAutostart
//tags: autostart
export async function initChemAutostart() : Promise<void> {
  PackageFunctions.initChemAutostart();
}

//name: Chemistry | Most Diverse Structures
//tags: tooltip
//input: column col { semType: Molecule }
//output: widget result
export async function chemTooltip(col: DG.Column) : Promise<any> {
  return PackageFunctions.chemTooltip(col);
}

//name: Scaffold Tree
//tags: viewer
//output: viewer result
//meta.icon: files/icons/scaffold-tree-icon.svg
export function scaffoldTreeViewer() : any {
  return PackageFunctions.scaffoldTreeViewer();
}

//name: Substructure Filter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
//meta.semType: Molecule
//meta.primaryFilter: true
//meta.allowMultipleFiltersForColumn: false
export function substructureFilter() : any {
  return PackageFunctions.substructureFilter();
}

//name: canvasMol
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

//name: drawMolecule
//input: string molStr 
//input: int w { optional: true }
//input: int h { optional: true }
//input: bool popupMenu { optional: true }
//output: object canvas
export function drawMolecule(molStr: string, w?: number, h?: number, popupMenu?: boolean) : any {
  return PackageFunctions.drawMolecule(molStr, w, h, popupMenu);
}

//name: getCLogP
//input: string smiles { semType: Molecule }
//output: double result
export function getCLogP(smiles: string) : number {
  return PackageFunctions.getCLogP(smiles);
}

//name: rdKitCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdKitCellRenderer() : Promise<any> {
  return PackageFunctions.rdKitCellRenderer();
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-ChemicalReaction
//output: grid_cell_renderer result
//meta.cellType: ChemicalReaction
export async function rdKitReactionRenderer() : Promise<any> {
  return PackageFunctions.rdKitReactionRenderer();
}

//name: chemMixtureRenderer
//tags: cellRenderer, cellRenderer-ChemicalMixture
//output: grid_cell_renderer result
//meta.cellType: ChemicalMixture
export async function rdKitMixtureRenderer() : Promise<any> {
  return PackageFunctions.rdKitMixtureRenderer();
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//output: grid_cell_renderer result
//meta.cellType: Molecule
export async function chemCellRenderer() : Promise<any> {
  return PackageFunctions.chemCellRenderer();
}

//name: getMorganFingerprints
//input: column molColumn { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export async function getMorganFingerprints(molColumn: DG.Column) : Promise<any> {
  return PackageFunctions.getMorganFingerprints(molColumn);
}

//name: getMorganFingerprint
//input: string molString { semType: Molecule }
//output: object fingerprintBitset {  }
export function getMorganFingerprint(molString: string) : any {
  return PackageFunctions.getMorganFingerprint(molString);
}

//name: getSimilarities
//input: column molStringsColumn 
//input: string molString 
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) : Promise<any> {
  return PackageFunctions.getSimilarities(molStringsColumn, molString);
}

//name: getDiversities
//input: column molStringsColumn 
//input: int limit 
//output: dataframe result
export async function getDiversities(molStringsColumn: DG.Column, limit: number) : Promise<any> {
  return PackageFunctions.getDiversities(molStringsColumn, limit);
}

//name: findSimilar
//input: column molStringsColumn 
//input: string molString 
//input: int limit 
//input: int cutoff 
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number, cutoff: number) : Promise<any> {
  return PackageFunctions.findSimilar(molStringsColumn, molString, limit, cutoff);
}

//name: searchSubstructure
//input: column molStringsColumn 
//input: string molString 
//input: string molBlockFailover 
//output: column result
export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, molBlockFailover: string) : Promise<any> {
  return PackageFunctions.searchSubstructure(molStringsColumn, molString, molBlockFailover);
}

//name: saveAsSdf
//description: As SDF...
//tags: fileExporter
export async function saveAsSdf() : Promise<void> {
  PackageFunctions.saveAsSdf();
}

//name: Chem Similarity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
export function similaritySearchViewer() : any {
  return PackageFunctions.similaritySearchViewer();
}

//name: Similarity Search
//top-menu: Chem | Search | Similarity Search...
export function similaritySearchTopMenu() : void {
  PackageFunctions.similaritySearchTopMenu();
}

//name: Chem Diversity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-diversity-search-viewer.svg
export function diversitySearchViewer() : any {
  return PackageFunctions.diversitySearchViewer();
}

//name: Diversity Search
//top-menu: Chem | Search | Diversity Search...
export function diversitySearchTopMenu() : void {
  PackageFunctions.diversitySearchTopMenu();
}

//name: descriptorsDocker
//top-menu: Chem | Calculate | Descriptors...
export async function descriptorsDocker() : Promise<void> {
  PackageFunctions.descriptorsDocker();
}

//name: calculateDescriptorsTransform
//tags: Transform
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: list<string> selected 
export async function calculateDescriptorsTransform(table: DG.DataFrame, molecules: DG.Column, selected: string[]) : Promise<void> {
  PackageFunctions.calculateDescriptorsTransform(table, molecules, selected);
}

//name: chemDescriptorsTree
//output: object descriptors
export async function chemDescriptorsTree() : Promise<any> {
  return PackageFunctions.chemDescriptorsTree();
}

//name: Map Identifiers
//top-menu: Chem | Calculate | Map Identifiers...
export async function getMapIdentifiers() : Promise<void> {
  PackageFunctions.getMapIdentifiers();
}

//name: freeTextToSmiles
//input: string molfile 
//output: string result
export async function freeTextToSmiles(molfile: string) : Promise<string> {
  return PackageFunctions.freeTextToSmiles(molfile);
}

//name: chemDescriptors
//input: dataframe table 
//input: column molecules 
//input: list<string> descriptors 
export async function chemDescriptors(table: DG.DataFrame, molecules: DG.Column, descriptors: string[]) : Promise<void> {
  PackageFunctions.chemDescriptors(table, molecules, descriptors);
}

//name: chemDescriptor
//input: column molecules { semType: Molecule }
//input: string descriptor 
//output: column result
//meta.vectorFunc: true
export async function chemDescriptor(molecules: DG.Column, descriptor: string) : Promise<any> {
  return PackageFunctions.chemDescriptor(molecules, descriptor);
}

//name: SearchSubstructureEditor
//tags: editor
//input: funccall call 
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

//name: clusterMCSTopMenu
//description: Calculates most common substructures for each cluster
//input: dataframe table 
//input: column molCol { semType: Molecule }
//input: column clusterCol { type: string }
//friendlyName: Cluster MCS
//top-menu: Chem | Calculate | Cluster MCS...
export async function clusterMCSTopMenu(table: DG.DataFrame, molCol: DG.Column, clusterCol: DG.Column) : Promise<void> {
  PackageFunctions.clusterMCSTopMenu(table, molCol, clusterCol);
}

//name: clusterMCS
//description: Calculates most common substructures for each cluster
//input: column molCol { semType: Molecule }
//input: column clusterCol { type: string }
//output: column result {  }
//meta.vectorFunc: true
//friendlyName: Cluster MCS
export async function performClusterMCS(molCol: DG.Column, clusterCol: DG.Column) : Promise<any> {
  return PackageFunctions.performClusterMCS(molCol, clusterCol);
}

//name: ChemSpaceEditor
//tags: editor
//input: funccall call 
export function ChemSpaceEditor(call: DG.FuncCall) : void {
  PackageFunctions.ChemSpaceEditor(call);
}

//name: Fingerprints
//tags: dim-red-preprocessing-function
//input: column col { semType: Molecule }
//input: string fingerprintType { caption: Fingerprint type; optional: true; choices: ['Morgan','RDKit','Pattern','AtomPair','MACCS','TopologicalTorsion']; default: Morgan }
//input: string _metric { optional: true }
//output: object result
//meta.supportedSemTypes: Molecule
//meta.supportedDistanceFunctions: Tanimoto, Asymmetric, Cosine, Sokal
export async function getFingerprints(col: DG.Column, fingerprintType: any, _metric?: string) {
  return PackageFunctions.getFingerprints(col, fingerprintType, _metric);
}

//name: Chem Space
//description: Maps the dataset to 2D plot based on similarity
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: bool plotEmbeddings { default: true }
//input: object options { optional: true }
//input: func preprocessingFunction { optional: true }
//input: bool clusterEmbeddings { optional: true }
//input: bool clusterMCS { optional: true }
//output: viewer result
//top-menu: Chem | Analyze | Chemical Space...
//editor: Chem:ChemSpaceEditor
export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options?: any, preprocessingFunction?: any, clusterEmbeddings?: boolean, clusterMCS?: boolean) : Promise<any> {
  return PackageFunctions.chemSpaceTopMenu(table, molecules, methodName, similarityMetric, plotEmbeddings, options, preprocessingFunction, clusterEmbeddings, clusterMCS);
}

//name: chemSpaceTransform
//tags: Transform
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string methodName 
//input: string similarityMetric 
//input: bool plotEmbeddings { default: true }
//input: string options { optional: true }
//input: bool clusterEmbeddings { optional: true }
//output: viewer result
export async function chemSpaceTransform(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options?: string, clusterEmbeddings?: boolean) : Promise<any> {
  return PackageFunctions.chemSpaceTransform(table, molecules, methodName, similarityMetric, plotEmbeddings, options, clusterEmbeddings);
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
  return PackageFunctions.getChemSpaceEmbeddings(col, methodName, similarityMetric, xAxis, yAxis, options);
}

//name: Chem Similarities Matrix
//input: int dim 
//input: column col 
//input: dataframe df 
//input: string colName 
//input: object simArr 
//output: object result
export async function getChemSimilaritiesMatrix(dim: number, col: DG.Column, df: DG.DataFrame, colName: string, simArr: DG.Column[]) : Promise<any> {
  return PackageFunctions.getChemSimilaritiesMatrix(dim, col, df, colName, simArr);
}

//name: Elemental Analysis
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: bool radarViewer { description: Add a standalone radar viewer; default: false }
//input: bool radarGrid { description: Show radar in grid cells; default: false }
//top-menu: Chem | Analyze | Elemental Analysis...
export async function elementalAnalysis(table: DG.DataFrame, molecules: DG.Column, radarViewer: boolean, radarGrid: boolean) : Promise<void> {
  PackageFunctions.elementalAnalysis(table, molecules, radarViewer, radarGrid);
}

//name: runElementalAnalysis
//tags: Transform
//input: dataframe table 
//input: column molecules { semType: Molecule }
//output: list res
export function runElementalAnalysis(table: DG.DataFrame, molecules: DG.Column) : string[] {
  return PackageFunctions.runElementalAnalysis(table, molecules);
}

//name: R-Groups Analysis
//top-menu: Chem | Analyze | R-Groups Analysis...
export function rGroupsAnalysisMenu() : void {
  PackageFunctions.rGroupsAnalysisMenu();
}

//name: rGroupDecomposition
//tags: Transform
//input: dataframe df 
//input: string molColName 
//input: string core 
//input: string rGroupName 
//input: string rGroupMatchingStrategy 
//input: bool onlyMatchAtRGroups { optional: true; default: false }
//output: object result
export async function rGroupDecomposition(df: DG.DataFrame, molColName: string, core: string, rGroupName: string, rGroupMatchingStrategy: string, onlyMatchAtRGroups: boolean) : Promise<any> {
  return PackageFunctions.rGroupDecomposition(df, molColName, core, rGroupName, rGroupMatchingStrategy, onlyMatchAtRGroups);
}

//name: ActivityCliffsEditor
//tags: editor
//input: funccall call 
export function ActivityCliffsEditor(call: DG.FuncCall) : void {
  PackageFunctions.ActivityCliffsEditor(call);
}

//name: Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table { description: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity { description: Similarity cutoff; default: 80 }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: func preprocessingFunction { optional: true }
//input: object options { optional: true }
//input: bool isDemo { optional: true }
//input: bool isTest { optional: true }
//top-menu: Chem | Analyze | Activity Cliffs...
//editor: Chem:ActivityCliffsEditor
export async function activityCliffs(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, preprocessingFunction: any, options?: any, isDemo?: boolean, isTest?: boolean) : Promise<void> {
  PackageFunctions.activityCliffs(table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, isDemo, isTest);
}

//name: activityCliffsInitFunction
//input: viewer sp 
export async function activityCliffsInitFunction(sp: any) : Promise<void> {
  PackageFunctions.activityCliffsInitFunction(sp);
}

//name: activityCliffsTransform
//tags: Transform
//input: dataframe table { description: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity { description: Similarity cutoff; default: 80 }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: string options { optional: true }
//input: bool isDemo { optional: true }
export async function activityCliffsTransform(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, options?: string, isDemo?: boolean) : Promise<void> {
  PackageFunctions.activityCliffsTransform(table, molecules, activities, similarity, methodName, similarityMetric, options, isDemo);
}

//name: To InchI
//tags: Transform
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//top-menu: Chem | Calculate | To InchI...
export function addInchisTopMenu(table: DG.DataFrame, col: DG.Column) : void {
  PackageFunctions.addInchisTopMenu(table, col);
}

//name: getInchis
//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchis(molecules: DG.Column) : any {
  return PackageFunctions.getInchis(molecules);
}

//name: To InchI Keys
//tags: Transform
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//top-menu: Chem | Calculate | To InchI Keys...
export function addInchisKeysTopMenu(table: DG.DataFrame, col: DG.Column) : void {
  PackageFunctions.addInchisKeysTopMenu(table, col);
}

//name: getInchiKeys
//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchiKeys(molecules: DG.Column) : any {
  return PackageFunctions.getInchiKeys(molecules);
}

//name: Structural Alerts
//description: Highlights the fragments that could lead to potential chemical hazards
//tags: HitTriageFunction
//input: dataframe table { description: Input data table; caption: Table }
//input: column molecules { caption: Molecules; semType: Molecule; type: categorical }
//input: bool pains { caption: PAINS; default: true; description: 'Pan Assay Interference Compounds filters' }
//input: bool bms { caption: BMS; default: false; description: 'Bristol-Myers Squibb HTS Deck filters' }
//input: bool sureChembl { caption: SureChEMBL; default: false; description: 'MedChem unfriendly compounds from SureChEMBL' }
//input: bool mlsmr { caption: MLSMR; default: false; description: 'NIH MLSMR Excluded Functionality filters' }
//input: bool dundee { caption: Dundee; default: false; description: 'University of Dundee NTD Screening Library filters' }
//input: bool inpharmatica { caption: Inpharmatica; default: false; description: 'Inpharmatica filters' }
//input: bool lint { caption: LINT; default: false; description: 'Pfizer LINT filters' }
//input: bool glaxo { caption: Glaxo; default: false; description: 'Glaxo Wellcome Hard filters' }
//output: dataframe result
//top-menu: Chem | Analyze | Structural Alerts...
export async function structuralAlertsTopMenu(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) : Promise<any> {
  return PackageFunctions.structuralAlertsTopMenu(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//name: runStructuralAlerts
//tags: Transform
//input: dataframe table { caption: Table; description: Input data table }
//input: column molecules { caption: Molecules; type: categorical; semType: Molecule }
//input: bool pains { caption: PAINS; default: true; description: 'Pan Assay Interference Compounds filters' }
//input: bool bms { caption: BMS; default: false; description: 'Bristol-Myers Squibb HTS Deck filters' }
//input: bool sureChembl { caption: SureChEMBL; default: false; description: 'MedChem unfriendly compounds from SureChEMBL' }
//input: bool mlsmr { caption: MLSMR; default: false; description: 'NIH MLSMR Excluded Functionality filters' }
//input: bool dundee { caption: Dundee; default: false; description: 'University of Dundee NTD Screening Library filters' }
//input: bool inpharmatica { caption: Inpharmatica; default: false; description: 'Inpharmatica filters' }
//input: bool lint { caption: LINT; default: false; description: 'Pfizer LINT filters' }
//input: bool glaxo { caption: Glaxo; default: false; description: 'Glaxo Wellcome Hard filters' }
//output: dataframe result
export async function runStructuralAlerts(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) : Promise<any> {
  return PackageFunctions.runStructuralAlerts(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//name: runStructuralAlert
//input: column<string> molecules { semType: Molecule }
//input: string alert 
//output: column res
//meta.vectorFunc: true
export async function runStructuralAlert(molecules: DG.Column, alert: any) : Promise<any> {
  return PackageFunctions.runStructuralAlert(molecules, alert);
}

//name: Chemistry | Rendering
//tags: panel, exclude-actions-panel
//input: column molColumn { semType: Molecule }
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column) : any {
  return PackageFunctions.molColumnPropertyPanel(molColumn);
}

//name: Chemistry | Highlight
//tags: panel, exclude-actions-panel
//input: column molColumn { semType: Molecule }
//output: widget result
export function molColumnHighlights(molColumn: DG.Column) : any {
  return PackageFunctions.molColumnHighlights(molColumn);
}

//name: Chemistry | Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string) : any {
  return PackageFunctions.descriptorsWidget(smiles);
}

//name: Biology | Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//help-url: /help/domains/chem/info-panels/drug-likeness.md
export function drugLikeness(smiles: DG.SemanticValue) : any {
  return PackageFunctions.drugLikeness(smiles);
}

//name: Chemistry | Properties
//description: Basic molecule properties
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export function properties(smiles: DG.SemanticValue) : any {
  return PackageFunctions.properties(smiles);
}

//name: getChemPropertyFunction
//description: Return chem property function
//input: string name 
//output: object result
export function getChemPropertyFunction(name: string) : any {
  return PackageFunctions.getChemPropertyFunction(name);
}

//name: Biology | Structural Alerts
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
//help-url: /help/domains/chem/info-panels/structural-alerts.md
export async function structuralAlerts(smiles: string) : Promise<any> {
  return PackageFunctions.structuralAlerts(smiles);
}

//name: Structure | Identifiers
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function identifiers(smiles: string) : Promise<any> {
  return PackageFunctions.identifiers(smiles);
}

//name: Structure | 3D Structure
//description: 3D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export async function structure3D(molecule: string) : Promise<any> {
  return PackageFunctions.structure3D(molecule);
}

//name: Structure | 2D Structure
//description: 2D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export function structure2d(molecule: string) : any {
  return PackageFunctions.structure2d(molecule);
}

//name: Biology | Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//help-url: /help/domains/chem/info-panels/toxicity-risks.md
export function toxicity(smiles: DG.SemanticValue) : any {
  return PackageFunctions.toxicity(smiles);
}

//name: convertMoleculeNotation
//input: column molecule { semType: Molecule }
//input: string targetNotation 
//output: column result
//meta.vectorFunc: true
export async function convertMoleculeNotation(molecule: DG.Column, targetNotation: any) : Promise<any> {
  return PackageFunctions.convertMoleculeNotation(molecule, targetNotation);
}

//name: convertMolNotation
//description: RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000
//tags: unitConverter
//input: string molecule { semType: Molecule }
//input: string sourceNotation { choices: ['smiles','smarts','molblock','v3Kmolblock'] }
//input: string targetNotation { choices: ['smiles','smarts','molblock','v3Kmolblock'] }
//output: string result {  }
export function convertMolNotation(molecule: string, sourceNotation: any, targetNotation: any) : string {
  return PackageFunctions.convertMolNotation(molecule, sourceNotation, targetNotation);
}

//name: Convert Notation
//tags: Transform
//input: dataframe data 
//input: column molecules { semType: Molecule }
//input: string targetNotation { choices: ['smiles','smarts','molblock','v3Kmolblock']; default: smiles }
//input: bool overwrite { default: false }
//input: bool join { default: true }
//output: column result
//top-menu: Chem | Transform | Convert Notation...
export async function convertNotation(data: DG.DataFrame, molecules: DG.Column<any>, targetNotation: any, overwrite: boolean, join: boolean) : Promise<any> {
  return PackageFunctions.convertNotation(data, molecules, targetNotation, overwrite, join);
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

//name: editMoleculeCell
//description: Molecule
//tags: cellEditor
//input: grid_cell cell 
export async function editMoleculeCell(cell: any) : Promise<void> {
  PackageFunctions.editMoleculeCell(cell);
}

//name: OpenChemLib
//tags: moleculeSketcher
//output: widget sketcher
export function openChemLibSketcher() : any {
  return PackageFunctions.openChemLibSketcher();
}

//name: importSdf
//description: Opens SDF file
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: sdf,mol
export function importSdf(bytes: Uint8Array) : any {
  return PackageFunctions.importSdf(bytes);
}

//name: importSmi
//description: Opens smi file
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: smi
export function importSmi(bytes: Uint8Array) : any {
  return PackageFunctions.importSmi(bytes);
}

//name: importMol2
//description: Opens smi file
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: mol2
export function importMol2(bytes: Uint8Array) : any {
  return PackageFunctions.importMol2(bytes);
}

//name: importMol
//description: Opens MOL file
//tags: file-handler
//input: string content 
//output: list<dataframe> result
//meta.ext: mol
export function importMol(content: string) : any {
  return PackageFunctions.importMol(content);
}

//name: oclCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: OpenChemLib
export async function oclCellRenderer() : Promise<any> {
  return PackageFunctions.oclCellRenderer();
}

//name: Sort by similarity
//description: Sorts a molecular column by similarity
//input: semantic_value value { semType: Molecule }
//meta.action: Sort by similarity
export async function sortBySimilarity(value: DG.SemanticValue) : Promise<void> {
  PackageFunctions.sortBySimilarity(value);
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
//tags: exclude-current-value-menu
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as...
export function copyAsAction(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsAction(value);
}

//name: Copy as SMILES
//description: Copies structure as smiles
//tags: exclude-actions-panel
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as SMILES
export function copyAsSmiles(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsSmiles(value);
}

//name: Copy as MOLFILE V2000
//description: Copies structure as molfile V2000
//tags: exclude-actions-panel
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as MOLFILE V2000
export function copyAsMolfileV2000(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsMolfileV2000(value);
}

//name: Copy as MOLFILE V3000
//description: Copies structure as molfile V3000
//tags: exclude-actions-panel
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as MOLFILE V3000
export function copyAsMolfileV3000(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsMolfileV3000(value);
}

//name: Copy as SMARTS
//description: Copies structure as smarts
//tags: exclude-actions-panel
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as SMARTS
export function copyAsSmarts(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsSmarts(value);
}

//name: Copy as IMAGE
//description: Copies structure as Image
//tags: exclude-actions-panel
//input: semantic_value value { semType: Molecule }
//meta.action: Copy as Image
export function copyAsImage(value: DG.SemanticValue) : void {
  PackageFunctions.copyAsImage(value);
}

//name: isSmiles
//input: string s 
//output: bool result
export function isSmiles(s: string) : boolean {
  return PackageFunctions.isSmiles(s);
}

//name: isSmarts
//input: string s 
//output: bool result
export function isSmarts(s: string) : boolean {
  return PackageFunctions.isSmarts(s);
}

//name: detectSmiles
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
  return PackageFunctions.callChemSimilaritySearch(df, col, molecule, metricName, fingerprint, limit, minScore);
}

//name: chemDiversitySearch
//input: column col 
//input: string metricName 
//input: string fingerprint 
//input: int limit 
//output: dataframe result
export async function callChemDiversitySearch(col: DG.Column, metricName: any, fingerprint: string, limit: number) : Promise<number[]> {
  return PackageFunctions.callChemDiversitySearch(col, metricName, fingerprint, limit);
}

//name: Chemical Properties
//tags: HitTriageFunction, Transform
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: bool MW { default: true }
//input: bool HBA { default: false }
//input: bool HBD { default: false }
//input: bool logP { default: false }
//input: bool logS { default: false }
//input: bool PSA { default: false }
//input: bool rotatableBonds { default: false }
//input: bool stereoCenters { default: false }
//input: bool moleculeCharge { default: false }
//output: dataframe result
//meta.function_family: biochem-calculator
//meta.method_info.author: Open Chem Lib Team
//meta.method_info.year: 2024
//meta.method_info.github: https://github.com/actelion/openchemlib
//top-menu: Chem | Calculate | Properties...
export async function addChemPropertiesColumns(table: DG.DataFrame, molecules: DG.Column, MW?: boolean, HBA?: boolean, HBD?: boolean, logP?: boolean, logS?: boolean, PSA?: boolean, rotatableBonds?: boolean, stereoCenters?: boolean, moleculeCharge?: boolean) : Promise<any> {
  return PackageFunctions.addChemPropertiesColumns(table, molecules, MW, HBA, HBD, logP, logS, PSA, rotatableBonds, stereoCenters, moleculeCharge);
}

//name: getMolProperty
//input: column molecules { semType: Molecule }
//input: string property { choices: ['MW','HBA','HBD','LogP','LogS','PSA','Rotatable bonds','Stereo centers','Molecule charge'] }
//output: column result
//meta.vectorFunc: true
export async function getMolProperty(molecules: DG.Column, property: string) : Promise<any> {
  return PackageFunctions.getMolProperty(molecules, property);
}

//name: Toxicity Risks
//tags: HitTriageFunction, Transform
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: bool mutagenicity { default: true }
//input: bool tumorigenicity { default: false }
//input: bool irritatingEffects { default: false }
//input: bool reproductiveEffects { default: false }
//output: dataframe result
//top-menu: Chem | Calculate | Toxicity Risks...
export async function addChemRisksColumns(table: DG.DataFrame, molecules: DG.Column, mutagenicity?: boolean, tumorigenicity?: boolean, irritatingEffects?: boolean, reproductiveEffects?: boolean) : Promise<any> {
  return PackageFunctions.addChemRisksColumns(table, molecules, mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects);
}

//name: addScaffoldTree
//description: Generates a hierarchical tree based on the scaffolds presented in dataset
//top-menu: Chem | Analyze | Scaffold Tree
export function addScaffoldTree() : void {
  PackageFunctions.addScaffoldTree();
}

//name: Matched Molecular Pairs Analysis
//tags: viewer
//output: viewer result
//meta.showInGallery: false
export function mmpViewer() : any {
  return PackageFunctions.mmpViewer();
}

//name: MMPEditor
//tags: editor
//input: funccall call 
export function MMPEditor(call: DG.FuncCall) : void {
  PackageFunctions.MMPEditor(call);
}

//name: Matched Molecular Pairs
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: column_list activities { type: numerical }
//input: string_list diffTypes 
//input: string_list scalings 
//input: double fragmentCutoff { description: Maximum fragment size relative to core; default: 0.4 }
//editor: Chem:MMPEditor
//top-menu: Chem | Analyze | Matched Molecular Pairs...
export async function mmpAnalysis(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column[], diffTypes: any, scalings: any, fragmentCutoff: number) : Promise<void> {
  PackageFunctions.mmpAnalysis(table, molecules, activities, diffTypes, scalings, fragmentCutoff);
}

//name: Scaffold Tree Filter
//description: Scaffold Tree filter
//tags: filter
//output: filter result
//meta.semType: Molecule
export function scaffoldTreeFilter() : any {
  return PackageFunctions.scaffoldTreeFilter();
}

//name: getScaffoldTree
//input: dataframe data 
//input: int ringCutoff { description: Ignore molecules with # rings > N; default: 10 }
//input: bool dischargeAndDeradicalize { description: Remove charges and radicals from scaffolds; default: false }
//output: string result
export async function getScaffoldTree(data: DG.DataFrame, ringCutoff: number, dischargeAndDeradicalize: boolean) : Promise<string> {
  return PackageFunctions.getScaffoldTree(data, ringCutoff, dischargeAndDeradicalize);
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
  PackageFunctions.demoChemOverview();
}

//name: Demo Similarity Search
//description: Searching for most similar or diverse molecules in dataset
//meta.demoPath: Cheminformatics | Similarity & Diversity Search
export async function demoSimilarityDiversitySearch() : Promise<void> {
  PackageFunctions.demoSimilarityDiversitySearch();
}

//name: Demo Matched Molecular Pairs
//description: Detect matched molecule pairs calculate the difference in activity values between them
//meta.demoPath: Cheminformatics | Matched Molecular Pairs
export async function demoMMPA() : Promise<void> {
  PackageFunctions.demoMMPA();
}

//name: Demo R Group Analysis
//description: R Group Analysis including R-group decomposition and  visual analysis of the obtained R-groups
//meta.demoPath: Cheminformatics | R-Group Analysis
//meta.demoSkip: GROK-14320
export async function demoRgroupAnalysis() : Promise<void> {
  PackageFunctions.demoRgroupAnalysis();
}

//name: Demo Activity Cliffs
//description: Searching similar structures with significant activity difference
//meta.demoPath: Cheminformatics | Molecule Activity Cliffs
//meta.demoSkip: GROK-14320
export async function demoMoleculeActivityCliffs() : Promise<void> {
  PackageFunctions.demoMoleculeActivityCliffs();
}

//name: Demo Chemical Space
//description: Maps the dataset to 2D plot based on similarity
//meta.demoPath: Cheminformatics | Chemical Space
//meta.demoSkip: GROK-14320
export async function demoChemicalSpace() : Promise<void> {
  PackageFunctions.demoChemicalSpace();
}

//name: Demo Scaffold Tree
//description: Running scaffold analysis with hierarchical tree
//meta.demoPath: Cheminformatics | Scaffold Tree
export async function demoScaffold() : Promise<void> {
  PackageFunctions.demoScaffold();
}

//name: Names To Smiles
//tags: Transform
//input: dataframe data 
//input: column names 
//top-menu: Chem | Transform | Names To Smiles...
export async function namesToSmiles(data: DG.DataFrame, names: DG.Column<any>) : Promise<void> {
  PackageFunctions.namesToSmiles(data, names);
}

//name: canonicalize
//input: string molecule { semType: molecule }
//output: string smiles {  }
//meta.role: canonicalizer
export function canonicalize(molecule: string) : string {
  return PackageFunctions.canonicalize(molecule);
}

//name: getMolecularFormula
//input: string molecule { semType: molecule }
//output: string molecularFormula
export function getMolecularFormula(molecule: string) : string {
  return PackageFunctions.getMolecularFormula(molecule);
}

//name: validateMolecule
//input: string s 
//output: object result
export function validateMolecule(s: string) : string {
  return PackageFunctions.validateMolecule(s);
}

//name: trainChemprop
//description: To be added
//input: dataframe df 
//input: column predictColumn 
//input: string dataset_type { category: General; choices: ['regression','classification']; default: regression; description: Type of dataset,e.g. classification or regression. This determines the loss function used during training. }
//input: string metric { category: General; choices: ['mse','mae','rmse','bounded-mse','bounded-mae','bounded-rmse','r2','binary-mcc','multiclass-mcc','roc','prc','accuracy','f1']; default: rmse; description: Metric to use during evaluation. Note:Does NOT affect loss function used during training (loss is determined by the `dataset_type` argument). }
//input: int multiclass_num_classes { category: General; default: 3; description: Number of classes when running multiclass classification }
//input: int num_folds { category: General; default: 1; description: Number of folds when performing cross validation }
//input: int data_seed { category: General; default: 0; description: Random seed to use when splitting data into train/val/test sets. When `num_folds` > 1,the first fold uses this seed and all subsequent folds add 1 to the seed. }
//input: list split_sizes { category: General; default: [0.8,0.1,0.1]; description: Split proportions for train/validation/test sets }
//input: string split_type { category: General; choices: ['random','scaffold_balanced','cv','cv_no_val','kennard_stone','kmeans','random_with_repeated_smiles']; default: random; description: Method of splitting the data into train/val/test }
//input: string activation { category: Model; choices: ['ReLU','LeakyReLU','PReLU','tanh','SELU','ELU']; default: ReLU; description: Activation function }
//input: bool atom_messages { category: Model; default: false; description: Use messages on atoms instead of messages on bonds }
//input: bool message_bias { category: Model; default: false; description: Whether to add bias to linear layers }
//input: int ensemble_size { category: Model; default: 1; description: Number of models in ensemble }
//input: int message_hidden_dim { category: Model; default: 300; description: Dimensionality of hidden layers in MPN }
//input: int depth { category: Model; default: 3; description: Number of message passing step }
//input: double dropout { category: Model; default: 0.0; description: Dropout probability }
//input: int ffn_hidden_dim { category: Model; default: 300; description: Hidden dim for higher-capacity FFN (defaults to hidden_size) }
//input: int ffn_num_layers { category: Model; default: 2; description: Number of layers in FFN after MPN encoding }
//input: int epochs { category: Training; default: 50; description: Number of epochs to run }
//input: int batch_size { category: Training; default: 64; description: Batch size }
//input: double warmup_epochs { category: Training; default: 2.0; description: Number of epochs during which learning rate increases linearly from init_lr to max_lr. Afterwards,learning rate decreases exponentially from max_lr to final_lr. }
//input: double init_lr { category: Training; default: 0.001; description: Initial learning rate }
//input: double max_lr { category: Training; default: 0.001; description: Maximum learning rate }
//input: double final_lr { category: Training; default: 0.0001; description: Final learning rate }
//input: bool no_descriptor_scaling { category: Training; default: false; description: Turn off scaling of features }
//output: dynamic model
//meta.mlname: Chemprop
//meta.mlrole: train
export async function trainChemprop(df: DG.DataFrame, predictColumn: DG.Column, dataset_type: string, metric: string, multiclass_num_classes: number, num_folds: number, data_seed: number, split_sizes: any, split_type: string, activation: string, atom_messages: boolean, message_bias: boolean, ensemble_size: number, message_hidden_dim: number, depth: number, dropout: number, ffn_hidden_dim: number, ffn_num_layers: number, epochs: number, batch_size: number, warmup_epochs: number, init_lr: number, max_lr: number, final_lr: number, no_descriptor_scaling: boolean) : Promise<Uint8Array> {
  return PackageFunctions.trainChemprop(df, predictColumn, dataset_type, metric, multiclass_num_classes, num_folds, data_seed, split_sizes, split_type, activation, atom_messages, message_bias, ensemble_size, message_hidden_dim, depth, dropout, ffn_hidden_dim, ffn_num_layers, epochs, batch_size, warmup_epochs, init_lr, max_lr, final_lr, no_descriptor_scaling);
}

//name: applyChemprop
//input: dataframe df 
//input: dynamic model 
//output: dataframe data_out
//meta.mlname: Chemprop
//meta.mlrole: apply
export async function applyChemprop(df: DG.DataFrame, model: Uint8Array) {
  return PackageFunctions.applyChemprop(df, model);
}

//name: isApplicableNN
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Chemprop
//meta.mlrole: isApplicable
export async function isApplicableNN(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableNN(df, predictColumn);
}

//name: isInteractiveNN
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Chemprop
//meta.mlrole: isInteractive
//meta.mlupdate: false
export async function isInteractiveNN(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveNN(df, predictColumn);
}

//name: Deprotect
//description: Generates the new dataset based on the given structure
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: string fragment { semType: Molecule; default: O=C([N:1])OCC1c2ccccc2-c2ccccc21 }
//top-menu: Chem | Transform | Deprotect...
export async function deprotect(table: DG.DataFrame, molecules: DG.Column, fragment: string) : Promise<void> {
  PackageFunctions.deprotect(table, molecules, fragment);
}

//name: beautifyMols
//description: Beautifies the list of molecules and returns the list of beautified molecules
//input: list<string> mols 
//output: list<string> result
export async function beautifyMols(mols: string[]) : Promise<string[]> {
  return PackageFunctions.beautifyMols(mols);
}

//name: convertToV3KViaOCL
//description: Converts the list of molecules to V3K format using OCL
//input: list<string> mols 
//output: list<string> result
export async function convertToV3KViaOCL(mols: string[]) : Promise<string[]> {
  return PackageFunctions.convertToV3KViaOCL(mols);
}

//name: mpo
//description: Calculates the MPO score for the column of molecules
//top-menu: Chem | Calculate | MPO Score...
export async function _mpo() : Promise<void> {
  PackageFunctions._mpo();
}

//name: mpoProfileEditor
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: json
//meta.fileViewerCheck: Chem:checkJsonMpoProfile
export function mpoProfileEditor(file: DG.FileInfo) : any {
  return PackageFunctions.mpoProfileEditor(file);
}

//name: checkJsonMpoProfile
//input: string content 
//output: bool result
export function checkJsonMpoProfile(content: string) {
  return PackageFunctions.checkJsonMpoProfile(content);
}

//name: Chemistry | Mixture
//tags: panel, chem, widgets
//input: string mixture { semType: ChemicalMixture }
//output: widget result
export async function mixtureWidget(mixture: string) : Promise<any> {
  return PackageFunctions.mixtureWidget(mixture);
}

//name: Chemistry | MixtureTree
//tags: panel, chem, widgets
//input: string mixture { semType: ChemicalMixture }
//output: widget result
export async function mixtureTreeWidget(mixture: string) : Promise<any> {
  return PackageFunctions.mixtureTreeWidget(mixture);
}

//name: Biochemical Properties
//description: Dynamically discovers and executes tagged biochemical calculators
//top-menu: Chem | Calculate | Biochemical Properties
export async function biochemPropsWidget() : Promise<void> {
  PackageFunctions.biochemPropsWidget();
}
