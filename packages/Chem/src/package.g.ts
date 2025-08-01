import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getRdKitModule
//output: dynamic result
export function getRdKitModule() {
  return PackageFunctions.getRdKitModule();
}

//name: getMolFileHandler
//input: string molString 
//output: dynamic result
export function getMolFileHandler(molString: string) {
  return PackageFunctions.getMolFileHandler(molString);
}

//name: init
//tags: init
export async function init() {
  return PackageFunctions.init();
}

//name: initChemAutostart
//tags: autostart
export async function initChemAutostart() {
  return PackageFunctions.initChemAutostart();
}

//name: Chemistry | Most Diverse Structures
//tags: tooltip
//input: column col { semType: Molecule }
//output: dynamic result
export async function chemTooltip(col: DG.Column) {
  return PackageFunctions.chemTooltip(col);
}

//name: Scaffold Tree
//tags: viewer
//output: viewer result
//meta.icon: files/icons/scaffold-tree-icon.svg
export function scaffoldTreeViewer() {
  return PackageFunctions.scaffoldTreeViewer();
}

//name: Substructure Filter
//description: RDKit-based substructure filter
//tags: filter
//output: filter result
//meta.semType: Molecule
//meta.primaryFilter: true
//meta.allowMultipleFiltersForColumn: false
export function substructureFilter() {
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
export function canvasMol(x: number, y: number, w: number, h: number, canvas: any, molString: string, scaffoldMolString: any, options: any) {
  return PackageFunctions.canvasMol(x, y, w, h, canvas, molString, scaffoldMolString, options);
}

//name: drawMolecule
//input: string molStr 
//input: int w { optional: true }
//input: int h { optional: true }
//input: bool popupMenu { optional: true }
//output: dynamic result
export function drawMolecule(molStr: string, w: number, h: number, popupMenu: boolean) {
  return PackageFunctions.drawMolecule(molStr, w, h, popupMenu);
}

//name: getCLogP
//input: string smiles { semType: Molecule }
//output: double result
export function getCLogP(smiles: string) {
  return PackageFunctions.getCLogP(smiles);
}

//name: rdKitCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: RDKit
export async function rdKitCellRenderer() {
  return PackageFunctions.rdKitCellRenderer();
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-ChemicalReaction
//output: grid_cell_renderer result
//meta.cellType: ChemicalReaction
//meta.cell-renderer-sem-type: ChemicalReaction
export async function rdKitReactionRenderer() {
  return PackageFunctions.rdKitReactionRenderer();
}

//name: chemMixtureRenderer
//tags: cellRenderer, cellRenderer-ChemicalMixture
//output: grid_cell_renderer result
//meta.cellType: ChemicalMixture
//meta.cell-renderer-sem-type: ChemicalMixture
export async function rdKitMixtureRenderer() {
  return PackageFunctions.rdKitMixtureRenderer();
}

//name: chemCellRenderer
//tags: cellRenderer, cellRenderer-Molecule
//output: grid_cell_renderer result
//meta.cellType: Molecule
//meta.cell-renderer-sem-type: Molecule
export async function chemCellRenderer() {
  return PackageFunctions.chemCellRenderer();
}

//name: getMorganFingerprints
//input: column molColumn { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export async function getMorganFingerprints(molColumn: DG.Column) {
  return PackageFunctions.getMorganFingerprints(molColumn);
}

//name: getMorganFingerprint
//input: string molString { semType: Molecule }
//output: dynamic result
export function getMorganFingerprint(molString: string) {
  return PackageFunctions.getMorganFingerprint(molString);
}

//name: getSimilarities
//input: column molStringsColumn 
//input: string molString 
//output: dataframe result
export async function getSimilarities(molStringsColumn: DG.Column, molString: string) {
  return PackageFunctions.getSimilarities(molStringsColumn, molString);
}

//name: getDiversities
//input: column molStringsColumn 
//input: int limit { optional: true }
//output: dataframe result
export async function getDiversities(molStringsColumn: DG.Column, limit: number) {
  return PackageFunctions.getDiversities(molStringsColumn, limit);
}

//name: findSimilar
//input: column molStringsColumn 
//input: string molString 
//input: int limit { optional: true }
//input: double cutoff { optional: true }
//output: dataframe result
export async function findSimilar(molStringsColumn: DG.Column, molString: string, limit: number, cutoff: number) {
  return PackageFunctions.findSimilar(molStringsColumn, molString, limit, cutoff);
}

//name: searchSubstructure
//input: column molStringsColumn 
//input: string molString 
//input: string molBlockFailover 
//output: column result
export async function searchSubstructure(molStringsColumn: DG.Column, molString: string, molBlockFailover: string) {
  return PackageFunctions.searchSubstructure(molStringsColumn, molString, molBlockFailover);
}

//name: saveAsSdf
//description: As SDF...
//tags: fileExporter
export async function saveAsSdf() {
  return PackageFunctions.saveAsSdf();
}

//name: Chem Similarity Search
//tags: viewer
//output: dynamic result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
export function similaritySearchViewer() {
  return PackageFunctions.similaritySearchViewer();
}

//name: Similarity Search
//top-menu: Chem | Search | Similarity Search...
export function similaritySearchTopMenu() {
  return PackageFunctions.similaritySearchTopMenu();
}

//name: Chem Diversity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-diversity-search-viewer.svg
export function diversitySearchViewer() {
  return PackageFunctions.diversitySearchViewer();
}

//name: Diversity Search
//top-menu: Chem | Search | Diversity Search...
export function diversitySearchTopMenu() {
  return PackageFunctions.diversitySearchTopMenu();
}

//name: descriptorsDocker
//top-menu: Chem | Calculate | Descriptors...
export async function descriptorsDocker() {
  return PackageFunctions.descriptorsDocker();
}

//name: chemDescriptorsTree
//output: dynamic result
export async function chemDescriptorsTree() {
  return PackageFunctions.chemDescriptorsTree();
}

//name: getMapIdentifiers
//output: dynamic result
//top-menu: Chem | Search | Map Identifiers...
export async function getMapIdentifiers() {
  return PackageFunctions.getMapIdentifiers();
}

//name: freeTextToSmiles
//input: string molfile 
//output: string result
export async function freeTextToSmiles(molfile: string) {
  return PackageFunctions.freeTextToSmiles(molfile);
}

//name: chemDescriptors
//input: dataframe table 
//input: column molecules 
//input: list<string> descriptors 
export async function chemDescriptors(table: DG.DataFrame, molecules: DG.Column, descriptors: string[]) {
  return PackageFunctions.chemDescriptors(table, molecules, descriptors);
}

//name: chemDescriptor
//input: column molecules { semType: Molecule }
//input: string descriptor 
//output: column result
//meta.vectorFunc: true
export async function chemDescriptor(molecules: DG.Column, descriptor: string) {
  return PackageFunctions.chemDescriptor(molecules, descriptor);
}

//name: searchSubstructureEditor
//tags: editor
//input: funccall call 
export function searchSubstructureEditor(call: DG.FuncCall) {
  return PackageFunctions.searchSubstructureEditor(call);
}

//name: Substructure Search
//input: column molecules { semType: Molecule }
//top-menu: Chem | Search | Substructure Search...
//editor: Chem:SearchSubstructureEditor
export function SubstructureSearchTopMenu(molecules: DG.Column) {
  return PackageFunctions.SubstructureSearchTopMenu(molecules);
}

//name: clusterMCSTopMenu
//description: Calculates most common substructures for each cluster
//input: dataframe table 
//input: column molCol { semType: Molecule }
//input: column clusterCol { type: string }
//friendlyName: Cluster MCS
//top-menu: Chem | Calculate | Cluster MCS...
export async function clusterMCSTopMenu(table: DG.DataFrame, molCol: DG.Column, clusterCol: DG.Column) {
  return PackageFunctions.clusterMCSTopMenu(table, molCol, clusterCol);
}

//name: performClusterMCS
//description: Calculates most common substructures for each cluster
//input: column molCol { semType: Molecule }
//input: column clusterCol { type: string }
//output: column result { semType: Molecule }
//meta.vectorFunc: true
//friendlyName: Cluster MCS
export async function performClusterMCS(molCol: DG.Column, clusterCol: DG.Column) {
  return PackageFunctions.performClusterMCS(molCol, clusterCol);
}

//name: ChemSpaceEditor
//tags: editor
//input: funccall call 
export function ChemSpaceEditor(call: DG.FuncCall) {
  return PackageFunctions.ChemSpaceEditor(call);
}

//name: Fingerprints
//tags: dim-red-preprocessing-function
//input: column col { semType: Molecule }
//input: string _metric { optional: true }
//input: dynamic fingerprintType { caption: Fingerprint type; optional: true; choices: ['Morgan','RDKit','Pattern','AtomPair','MACCS','TopologicalTorsion']; default: Morgan }
//output: dynamic result
//meta.supportedSemTypes: Molecule
//meta.supportedDistanceFunctions: Tanimoto, Asymmetric, Cosine, Sokal
export async function getFingerprints(col: DG.Column, _metric: string, fingerprintType: any) {
  return PackageFunctions.getFingerprints(col, _metric, fingerprintType);
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
//output: dynamic result
//top-menu: Chem | Analyze | Chemical Space...
//editor: Chem:ChemSpaceEditor
export async function chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options: any, preprocessingFunction: any, clusterEmbeddings: boolean, clusterMCS: boolean) {
  return PackageFunctions.chemSpaceTopMenu(table, molecules, methodName, similarityMetric, plotEmbeddings, options, preprocessingFunction, clusterEmbeddings, clusterMCS);
}

//name: chemSpaceTransform
//tags: Transform
//input: dataframe table 
//input: column molecules 
//input: dynamic methodName 
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: bool plotEmbeddings { default: true }
//input: string options { optional: true }
//input: bool clusterEmbeddings { optional: true }
//output: dynamic result
export async function chemSpaceTransform(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, options: string, clusterEmbeddings: boolean) {
  return PackageFunctions.chemSpaceTransform(table, molecules, methodName, similarityMetric, plotEmbeddings, options, clusterEmbeddings);
}

//name: Chem Space Embeddings
//input: column col 
//input: dynamic methodName 
//input: dynamic similarityMetric 
//input: string xAxis 
//input: string yAxis 
//input: dynamic options 
//output: dynamic result
export async function getChemSpaceEmbeddings(col: DG.Column, methodName: any, similarityMetric: any, xAxis: string, yAxis: string, options: any) {
  return PackageFunctions.getChemSpaceEmbeddings(col, methodName, similarityMetric, xAxis, yAxis, options);
}

//name: Chem Similarities Matrix
//input: double dim 
//input: column col 
//input: dataframe df 
//input: string colName 
//input: object simArr 
//output: dynamic result
export async function getChemSimilaritiesMatrix(dim: number, col: DG.Column, df: DG.DataFrame, colName: string, simArr: DG.Column[]) {
  return PackageFunctions.getChemSimilaritiesMatrix(dim, col, df, colName, simArr);
}

//name: Elemental Analysis
//description: Add a standalone radar viewer
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: bool radarViewer { description: Add a standalone radar viewer; default: false }
//input: bool radarGrid { description: Show radar in grid cells; default: false }
//top-menu: Chem | Analyze | Elemental Analysis...
export async function elementalAnalysis(table: DG.DataFrame, molecules: DG.Column, radarViewer: boolean, radarGrid: boolean) {
  return PackageFunctions.elementalAnalysis(table, molecules, radarViewer, radarGrid);
}

//name: runElementalAnalysis
//tags: Transform
//input: dataframe table 
//input: column molecules 
//output: list res
export function runElementalAnalysis(table: DG.DataFrame, molecules: DG.Column) {
  return PackageFunctions.runElementalAnalysis(table, molecules);
}

//name: R-Groups Analysis
//top-menu: Chem | Analyze | R-Groups Analysis...
export function rGroupsAnalysisMenu() {
  return PackageFunctions.rGroupsAnalysisMenu();
}

//name: rGroupDecomposition
//tags: Transform
//input: dataframe df 
//input: string molColName 
//input: string core 
//input: string rGroupName 
//input: string rGroupMatchingStrategy 
//input: bool onlyMatchAtRGroups { optional: true; default: false }
//output: dynamic result
export async function rGroupDecomposition(df: DG.DataFrame, molColName: string, core: string, rGroupName: string, rGroupMatchingStrategy: string, onlyMatchAtRGroups: boolean) {
  return PackageFunctions.rGroupDecomposition(df, molColName, core, rGroupName, rGroupMatchingStrategy, onlyMatchAtRGroups);
}

//name: ActivityCliffsEditor
//tags: editor
//input: funccall call 
export function ActivityCliffsEditor(call: DG.FuncCall) {
  return PackageFunctions.ActivityCliffsEditor(call);
}

//name: Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table { caption: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity { caption: Similarity cutoff; default: 80 }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: func preprocessingFunction { optional: true }
//input: object options { optional: true }
//input: bool isDemo { optional: true }
//input: bool isTest { optional: true }
//top-menu: Chem | Analyze | Activity Cliffs...
//editor: Chem:ActivityCliffsEditor
export async function activityCliffs(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, preprocessingFunction: any, options: any, isDemo: boolean, isTest: boolean) {
  return PackageFunctions.activityCliffs(table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, isDemo, isTest);
}

//name: activityCliffsInitFunction
//input: viewer sp 
export async function activityCliffsInitFunction(sp: any) {
  return PackageFunctions.activityCliffsInitFunction(sp);
}

//name: activityCliffsTransform
//tags: Transform
//input: dataframe table { caption: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: column activities { type: numerical }
//input: double similarity { caption: Similarity cutoff; default: 80 }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: string options { optional: true }
//input: bool isDemo { optional: true }
export async function activityCliffsTransform(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, options: string, isDemo: boolean) {
  return PackageFunctions.activityCliffsTransform(table, molecules, activities, similarity, methodName, similarityMetric, options, isDemo);
}

//name: To InchI
//tags: Transform
//input: dataframe table { caption: Input data table }
//input: column col { semType: Molecule }
//top-menu: Chem | Calculate | To InchI...
export function addInchisTopMenu(table: DG.DataFrame, col: DG.Column) {
  return PackageFunctions.addInchisTopMenu(table, col);
}

//name: getInchis
//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchis(molecules: DG.Column) {
  return PackageFunctions.getInchis(molecules);
}

//name: To InchI Keys
//tags: Transform
//input: dataframe table { caption: Input data table }
//input: column col { semType: Molecule }
//top-menu: Chem | Calculate | To InchI Keys...
export function addInchisKeysTopMenu(table: DG.DataFrame, col: DG.Column) {
  return PackageFunctions.addInchisKeysTopMenu(table, col);
}

//name: getInchiKeys
//input: column<string> molecules { semType: Molecule }
//output: column result
//meta.vectorFunc: true
export function getInchiKeys(molecules: DG.Column) {
  return PackageFunctions.getInchiKeys(molecules);
}

//name: Structural Alerts
//description: Highlights the fragments that could lead to potential chemical hazards
//input: dataframe table { caption: Input data table }
//input: column molecules { caption: Molecules; semType: Molecule }
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
export async function structuralAlertsTopMenu(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) {
  return PackageFunctions.structuralAlertsTopMenu(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//name: runStructuralAlerts
//tags: Transform
//input: dataframe table { caption: Input data table }
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
export async function runStructuralAlerts(table: DG.DataFrame, molecules: DG.Column, pains: boolean, bms: boolean, sureChembl: boolean, mlsmr: boolean, dundee: boolean, inpharmatica: boolean, lint: boolean, glaxo: boolean) {
  return PackageFunctions.runStructuralAlerts(table, molecules, pains, bms, sureChembl, mlsmr, dundee, inpharmatica, lint, glaxo);
}

//name: runStructuralAlert
//input: column<string> molecules { semType: Molecule }
//input: string alert 
//output: column res
//meta.vectorFunc: true
export async function runStructuralAlert(molecules: DG.Column, alert: any) {
  return PackageFunctions.runStructuralAlert(molecules, alert);
}

//name: Chemistry | Rendering
//tags: panel, exclude-actions-panel
//input: column molColumn { semType: Molecule }
//output: widget result
export function molColumnPropertyPanel(molColumn: DG.Column) {
  return PackageFunctions.molColumnPropertyPanel(molColumn);
}

//name: Chemistry | Highlight
//tags: panel, exclude-actions-panel
//input: column molColumn { semType: Molecule }
//output: widget result
export function molColumnHighlights(molColumn: DG.Column) {
  return PackageFunctions.molColumnHighlights(molColumn);
}

//name: Chemistry | Descriptors
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function descriptorsWidget(smiles: string) {
  return PackageFunctions.descriptorsWidget(smiles);
}

//name: Biology | Drug Likeness
//description: Drug Likeness score, with explanations on molecule fragments contributing to the score. OCL.
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//helpUrl: /help/domains/chem/info-panels/drug-likeness.md
export function drugLikeness(smiles: DG.SemanticValue) {
  return PackageFunctions.drugLikeness(smiles);
}

//name: Chemistry | Properties
//description: Basic molecule properties
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export function properties(smiles: DG.SemanticValue) {
  return PackageFunctions.properties(smiles);
}

//name: getChemPropertyFunction
//description: Return chem property function
//input: string name 
//output: dynamic result
export function getChemPropertyFunction(name: string) {
  return PackageFunctions.getChemPropertyFunction(name);
}

//name: Biology | Structural Alerts
//description: Screening drug candidates against structural alerts i.e. fragments associated to a toxicological response
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
//helpUrl: /help/domains/chem/info-panels/structural-alerts.md
export async function structuralAlerts(smiles: string) {
  return PackageFunctions.structuralAlerts(smiles);
}

//name: Structure | Identifiers
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export async function identifiers(smiles: string) {
  return PackageFunctions.identifiers(smiles);
}

//name: Structure | 3D Structure
//description: 3D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export async function structure3D(molecule: string) {
  return PackageFunctions.structure3D(molecule);
}

//name: Structure | 2D Structure
//description: 2D molecule representation
//tags: panel, chem, widgets
//input: string molecule { semType: Molecule }
//output: widget result
export function structure2d(molecule: string) {
  return PackageFunctions.structure2d(molecule);
}

//name: Biology | Toxicity
//description: Toxicity prediction. Calculated by openchemlib
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//helpUrl: /help/domains/chem/info-panels/toxicity-risks.md
export function toxicity(smiles: DG.SemanticValue) {
  return PackageFunctions.toxicity(smiles);
}

//name: convertMoleculeNotation
//input: column molecule { semType: Molecule }
//input: string targetNotation { choices: ['smiles','smarts','molblock','v3Kmolblock'] }
//output: column result
//meta.vectorFunc: true
export async function convertMoleculeNotation(molecule: DG.Column, targetNotation: any) {
  return PackageFunctions.convertMoleculeNotation(molecule, targetNotation);
}

//name: convertMolNotation
//description: RDKit-based conversion for SMILES, SMARTS, InChi, Molfile V2000 and Molfile V3000
//tags: unitConverter
//input: string molecule { semType: Molecule }
//input: string sourceNotation { choices: ['smiles','smarts','molblock','v3Kmolblock'] }
//input: string targetNotation { choices: ['smiles','smarts','molblock','v3Kmolblock'] }
//output: string result { semType: Molecule }
export function convertMolNotation(molecule: string, sourceNotation: any, targetNotation: any) {
  return PackageFunctions.convertMolNotation(molecule, sourceNotation, targetNotation);
}

//name: Convert Notation
//tags: Transform
//input: dataframe data 
//input: dynamic molecules 
//input: dynamic targetNotation { choices: ['smiles','smarts','molblock','v3Kmolblock']; default: smiles }
//input: bool overwrite { optional: true; default: false }
//input: bool join { optional: true; default: true }
//output: column result
//top-menu: Chem | Transform | Convert Notation...
export async function convertNotation(data: DG.DataFrame, molecules: any, targetNotation: any, overwrite: boolean, join: boolean) {
  return PackageFunctions.convertNotation(data, molecules, targetNotation, overwrite, join);
}

//name: Convert Notation...
//input: column col { semType: Molecule }
//output: dynamic result
//meta.action: Convert Notation...
export function convertMolNotationAction(col: DG.Column) {
  return PackageFunctions.convertMolNotationAction(col);
}

//name: Convert Mixture To Smiles...
//input: column col { semType: ChemicalMixture }
//meta.action: Convert mixture to smiles...
export function convertMixtureToSmiles(col: DG.Column) {
  return PackageFunctions.convertMixtureToSmiles(col);
}

//name: Molecule
//tags: cellEditor
//input: grid_cell cell 
export async function editMoleculeCell(cell: any) {
  return PackageFunctions.editMoleculeCell(cell);
}

//name: OpenChemLib
//tags: moleculeSketcher
//output: widget sketcher
export function openChemLibSketcher() {
  return PackageFunctions.openChemLibSketcher();
}

//name: importSdf
//description: Opens SDF file
//tags: file-handler
//input: list bytes 
//output: list result
//meta.ext: sdf,mol
export function importSdf(bytes: Uint8Array) {
  return PackageFunctions.importSdf(bytes);
}

//name: importSmi
//description: Opens smi file
//tags: file-handler
//input: list bytes 
//output: list result
//meta.ext: smi
export function importSmi(bytes: Uint8Array) {
  return PackageFunctions.importSmi(bytes);
}

//name: importMol2
//description: Opens smi file
//tags: file-handler
//input: blob bytes 
//output: list result
//meta.ext: mol2
export function importMol2(bytes: Uint8Array) {
  return PackageFunctions.importMol2(bytes);
}

//name: importMol
//description: Opens MOL file
//tags: file-handler
//input: string content 
//output: list result
//meta.ext: mol
export function importMol(content: string) {
  return PackageFunctions.importMol(content);
}

//name: oclCellRenderer
//output: grid_cell_renderer result
//meta.chemRendererName: OpenChemLib
export async function oclCellRenderer() {
  return PackageFunctions.oclCellRenderer();
}

//name: Sort by similarity
//description: Sorts a molecular column by similarity
//input: semantic_value value { semType: Molecule }
//meta.action: Sort by similarity
export async function sortBySimilarity(value: DG.SemanticValue) {
  return PackageFunctions.sortBySimilarity(value);
}

//name: Use as filter
//description: Adds this structure as a substructure filter
//input: semantic_value value { semType: Molecule }
//meta.action: Use as filter
export function useAsSubstructureFilter(value: DG.SemanticValue) {
  return PackageFunctions.useAsSubstructureFilter(value);
}

//name: Copy as...
//description: Copies structure in different formats
//tags: exclude-current-value-menu
//input: semantic_value value 
//output: dynamic result
//meta.action: Copy as...
export function copyAsAction(value: DG.SemanticValue) {
  return PackageFunctions.copyAsAction(value);
}

//name: Copy as SMILES
//description: Copies structure as smiles
//tags: exclude-actions-panel
//input: semantic_value value 
//meta.action: Copy as SMILES
export function copyAsSmiles(value: DG.SemanticValue) {
  return PackageFunctions.copyAsSmiles(value);
}

//name: Copy as MOLFILE V2000
//description: Copies structure as molfile V2000
//tags: exclude-actions-panel
//input: semantic_value value 
//meta.action: Copy as MOLFILE V2000
export function copyAsMolfileV2000(value: DG.SemanticValue) {
  return PackageFunctions.copyAsMolfileV2000(value);
}

//name: Copy as MOLFILE V3000
//description: Copies structure as molfile V3000
//tags: exclude-actions-panel
//input: semantic_value value 
//meta.action: Copy as MOLFILE V3000
export function copyAsMolfileV3000(value: DG.SemanticValue) {
  return PackageFunctions.copyAsMolfileV3000(value);
}

//name: Copy as SMARTS
//description: Copies structure as smarts
//tags: exclude-actions-panel
//input: semantic_value value 
//meta.action: Copy as SMARTS
export function copyAsSmarts(value: DG.SemanticValue) {
  return PackageFunctions.copyAsSmarts(value);
}

//name: Copy as IMAGE
//description: Copies structure as Image
//tags: exclude-actions-panel
//input: semantic_value value 
//meta.action: Copy as Image
export function copyAsImage(value: DG.SemanticValue) {
  return PackageFunctions.copyAsImage(value);
}

//name: isSmiles
//input: string s 
//output: bool result
export function isSmiles(s: string) {
  return PackageFunctions.isSmiles(s);
}

//name: isSmarts
//input: string s 
//output: bool result
export function isSmarts(s: string) {
  return PackageFunctions.isSmarts(s);
}

//name: detectSmiles
//input: column col 
//input: int min 
export function detectSmiles(col: DG.Column, min: number) {
  return PackageFunctions.detectSmiles(col, min);
}

//name: chemSimilaritySearch
//input: dataframe df 
//input: column col 
//input: string molecule 
//input: string metricName { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: int limit { optional: true }
//input: double minScore { optional: true }
//input: string fingerprint 
//output: dataframe result
export async function callChemSimilaritySearch(df: DG.DataFrame, col: DG.Column, molecule: string, metricName: any, limit: number, minScore: number, fingerprint: string) {
  return PackageFunctions.callChemSimilaritySearch(df, col, molecule, metricName, limit, minScore, fingerprint);
}

//name: chemDiversitySearch
//input: column col 
//input: string metricName { choices: ['Tanimoto','Asymmetric','Cosine','Sokal'] }
//input: int limit { optional: true }
//input: string fingerprint 
//output: dataframe result
export async function callChemDiversitySearch(col: DG.Column, metricName: any, limit: number, fingerprint: string) {
  return PackageFunctions.callChemDiversitySearch(col, metricName, limit, fingerprint);
}

//name: Chemical Properties
//tags: HitTriageFunction, Transform
//input: dataframe table { caption: Input data table }
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
//top-menu: Chem | Calculate | Chemical Properties...
export async function addChemPropertiesColumns(table: DG.DataFrame, molecules: DG.Column, MW: boolean, HBA: boolean, HBD: boolean, logP: boolean, logS: boolean, PSA: boolean, rotatableBonds: boolean, stereoCenters: boolean, moleculeCharge: boolean) {
  return PackageFunctions.addChemPropertiesColumns(table, molecules, MW, HBA, HBD, logP, logS, PSA, rotatableBonds, stereoCenters, moleculeCharge);
}

//name: getMolProperty
//input: column molecules { semType: Molecule }
//input: string property { choices: ['MW','HBA','HBD','LogP','LogS','PSA','Rotatable bonds','Stereo centers','Molecule charge'] }
//output: column result
//meta.vectorFunc: true
export async function getMolProperty(molecules: DG.Column, property: string) {
  return PackageFunctions.getMolProperty(molecules, property);
}

//name: Toxicity Risks
//tags: HitTriageFunction, Transform
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: bool mutagenicity { default: true }
//input: bool tumorigenicity { default: false }
//input: bool irritatingEffects { default: false }
//input: bool reproductiveEffects { default: false }
//output: dataframe result
//top-menu: Chem | Calculate | Toxicity Risks...
export async function addChemRisksColumns(table: DG.DataFrame, molecules: DG.Column, mutagenicity: boolean, tumorigenicity: boolean, irritatingEffects: boolean, reproductiveEffects: boolean) {
  return PackageFunctions.addChemRisksColumns(table, molecules, mutagenicity, tumorigenicity, irritatingEffects, reproductiveEffects);
}

//name: addScaffoldTree
//description: Generates a hierarchical tree based on the scaffolds presented in dataset
//top-menu: Chem | Analyze | Scaffold Tree
export function addScaffoldTree() {
  return PackageFunctions.addScaffoldTree();
}

//name: Matched Molecular Pairs Analysis
//tags: viewer
//output: viewer result
//meta.showInGallery: false
export function mmpViewer() {
  return PackageFunctions.mmpViewer();
}

//name: MMPEditor
//tags: editor
//input: funccall call 
export function MMPEditor(call: DG.FuncCall) {
  return PackageFunctions.MMPEditor(call);
}

//name: Matched Molecular Pairs
//tags: top-menu
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: column_list activities { type: numerical }
//input: string_list diffTypes 
//input: string_list scalings 
//input: double fragmentCutoff { description: Maximum fragment size relative to core; default: 0.4 }
//input: bool demo 
//editor: Chem:MMPEditor
//top-menu: Chem | Analyze | Matched Molecular Pairs...
export async function mmpAnalysis(table: DG.DataFrame, molecules: DG.Column, activities: DG.Column[], diffTypes: any, scalings: any, fragmentCutoff: number, demo: boolean) {
  return PackageFunctions.mmpAnalysis(table, molecules, activities, diffTypes, scalings, fragmentCutoff, demo);
}

//name: Scaffold Tree Filter
//description: Scaffold Tree filter
//tags: filter
//output: filter result
//meta.semType: Molecule
export function scaffoldTreeFilter() {
  return PackageFunctions.scaffoldTreeFilter();
}

//name: getScaffoldTree
//input: dataframe data 
//input: int ringCutoff { description: Ignore molecules with # rings > N; default: 10 }
//input: bool dischargeAndDeradicalize { description: Remove charges and radicals from scaffolds; default: false }
//output: string result
export async function getScaffoldTree(data: DG.DataFrame, ringCutoff: number, dischargeAndDeradicalize: boolean) {
  return PackageFunctions.getScaffoldTree(data, ringCutoff, dischargeAndDeradicalize);
}

//name: filterMoleculeDuplicates
//input: list<string> molecules 
//input: string molecule 
//output: list result
export function removeDuplicates(molecules: string[], molecule: string) {
  return PackageFunctions.removeDuplicates(molecules, molecule);
}

//name: Demo Chem Overview
//description: Overview of Cheminformatics functionality
//meta.isDemoScript: true
//meta.demoSkip: GROK-14320
//meta.demoPath: Cheminformatics | Overview
export async function demoChemOverview() {
  return PackageFunctions.demoChemOverview();
}

//name: Demo Similarity Search
//description: Searching for most similar or diverse molecules in dataset
//meta.demoPath: Cheminformatics | Similarity & Diversity Search
export async function demoSimilarityDiversitySearch() {
  return PackageFunctions.demoSimilarityDiversitySearch();
}

//name: Demo Matched Molecular Pairs
//description: Detect matched molecule pairs calculate the difference in activity values between them
//meta.demoPath: Cheminformatics | Matched Molecular Pairs
export async function demoMMPA() {
  return PackageFunctions.demoMMPA();
}

//name: Demo R Group Analysis
//description: R Group Analysis including R-group decomposition and  visual analysis of the obtained R-groups
//meta.demoPath: Cheminformatics | R-Group Analysis
//meta.demoSkip: GROK-14320
export async function demoRgroupAnalysis() {
  return PackageFunctions.demoRgroupAnalysis();
}

//name: Demo Activity Cliffs
//description: Searching similar structures with significant activity difference
//meta.demoPath: Cheminformatics | Molecule Activity Cliffs
//meta.demoSkip: GROK-14320
export async function demoMoleculeActivityCliffs() {
  return PackageFunctions.demoMoleculeActivityCliffs();
}

//name: Demo Chemical Space
//description: Maps the dataset to 2D plot based on similarity
//meta.demoPath: Cheminformatics | Chemical Space
//meta.demoSkip: GROK-14320
export async function demoChemicalSpace() {
  return PackageFunctions.demoChemicalSpace();
}

//name: Demo Scaffold Tree
//description: Running scaffold analysis with hierarchical tree
//meta.demoPath: Cheminformatics | Scaffold Tree
export async function demoScaffold() {
  return PackageFunctions.demoScaffold();
}

//name: Names To Smiles
//tags: Transform
//input: dataframe data 
//input: dynamic names 
//top-menu: Chem | Transform | Names To Smiles...
export async function namesToSmiles(data: DG.DataFrame, names: any) {
  return PackageFunctions.namesToSmiles(data, names);
}

//name: canonicalize
//input: string molecule { semType: molecule }
//output: string smiles { semType: molecule }
//meta.role: canonicalizer
export function canonicalize(molecule: string) {
  return PackageFunctions.canonicalize(molecule);
}

//name: getMolecularFormula
//input: string molecule { semType: molecule }
//output: string molecularFormula
export function getMolecularFormula(molecule: string) {
  return PackageFunctions.getMolecularFormula(molecule);
}

//name: validateMolecule
//input: string s 
//output: string result
export function validateMolecule(s: string) {
  return PackageFunctions.validateMolecule(s);
}

//name: trainChemprop
//description: To be added
//input: dataframe df 
//input: column predictColumn 
//input: string dataset_type { category: General; choices: ['regression','classification']; default: regression }
//input: string metric { category: General; choices: ['mse','mae','rmse','bounded-mse','bounded-mae','bounded-rmse','r2','binary-mcc','multiclass-mcc','roc','prc','accuracy','f1']; default: rmse }
//input: int multiclass_num_classes { category: General; default: 3 }
//input: int num_folds { category: General; default: 1 }
//input: int data_seed { category: General; default: 0 }
//input: list split_sizes { category: General }
//input: string split_type { category: General; choices: ['random','scaffold_balanced','cv','cv_no_val','kennard_stone','kmeans','random_with_repeated_smiles']; default: random }
//input: string activation { category: Model; choices: ['ReLU','LeakyReLU','PReLU','tanh','SELU','ELU']; default: ReLU }
//input: bool atom_messages { category: Model; default: false }
//input: bool message_bias { category: Model; default: false }
//input: int ensemble_size { category: Model; default: 1 }
//input: int message_hidden_dim { category: Model; default: 300 }
//input: int depth { category: Model; default: 3 }
//input: double dropout { category: Model; default: 0.0 }
//input: int ffn_hidden_dim { category: Model; default: 300 }
//input: int ffn_num_layers { category: Model; default: 2 }
//input: int epochs { category: Training; default: 50 }
//input: int batch_size { category: Training; default: 64 }
//input: double warmup_epochs { category: Training; default: 2.0 }
//input: double init_lr { category: Training; default: 0.001 }
//input: double max_lr { category: Training; default: 0.001 }
//input: double final_lr { category: Training; default: 0.0001 }
//input: bool no_descriptor_scaling { category: Training; default: false }
//output: dynamic model modelBlob
//meta.mlname: Chemprop
//meta.mlrole: train
export async function trainChemprop(df: DG.DataFrame, predictColumn: DG.Column, dataset_type: string, metric: string, multiclass_num_classes: number, num_folds: number, data_seed: number, split_sizes: any, split_type: string, activation: string, atom_messages: boolean, message_bias: boolean, ensemble_size: number, message_hidden_dim: number, depth: number, dropout: number, ffn_hidden_dim: number, ffn_num_layers: number, epochs: number, batch_size: number, warmup_epochs: number, init_lr: number, max_lr: number, final_lr: number, no_descriptor_scaling: boolean) {
  return PackageFunctions.trainChemprop(df, predictColumn, dataset_type, metric, multiclass_num_classes, num_folds, data_seed, split_sizes, split_type, activation, atom_messages, message_bias, ensemble_size, message_hidden_dim, depth, dropout, ffn_hidden_dim, ffn_num_layers, epochs, batch_size, warmup_epochs, init_lr, max_lr, final_lr, no_descriptor_scaling);
}

//name: applyChemprop
//input: dataframe df 
//input: blob model 
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
//input: dataframe table { caption: Input data table }
//input: column molecules { semType: Molecule }
//input: string fragment { semType: Molecule; default: O=C([N:1])OCC1c2ccccc2-c2ccccc21 }
//top-menu: Chem | Transform | Deprotect...
export async function deprotect(table: DG.DataFrame, molecules: DG.Column, fragment: string) {
  return PackageFunctions.deprotect(table, molecules, fragment);
}

//name: beautifyMols
//description: Beautifies the list of molecules and returns the list of beautified molecules
//input: list<string> mols 
//output: list result
export async function beautifyMols(mols: string[]) {
  return PackageFunctions.beautifyMols(mols);
}

//name: convertToV3KViaOCL
//description: Converts the list of molecules to V3K format using OCL
//input: list<string> mols 
//output: list result
export async function convertToV3KViaOCL(mols: string[]) {
  return PackageFunctions.convertToV3KViaOCL(mols);
}

//name: mpo
//description: Calculates the MPO score for the column of molecules
//top-menu: Chem | Calculate | MPO Score...
export async function _mpo() {
  return PackageFunctions._mpo();
}

//name: mpoProfileEditor
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: json
//meta.fileViewerCheck: Chem:checkJsonMpoProfile
export function mpoProfileEditor(file: DG.FileInfo) {
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
export async function mixtureWidget(mixture: string) {
  return PackageFunctions.mixtureWidget(mixture);
}

//name: Chemistry | MixtureTree
//tags: panel, chem, widgets
//input: string mixture { semType: ChemicalMixture }
//output: widget result
export async function mixtureTreeWidget(mixture: string) {
  return PackageFunctions.mixtureTreeWidget(mixture);
}

//name: Biochemical Properties
//description: Dynamically discovers and executes tagged biochemical calculators
//top-menu: Chem | Calculate | Biochemical Properties
export async function biochemPropsWidget() {
  return PackageFunctions.biochemPropsWidget();
}
