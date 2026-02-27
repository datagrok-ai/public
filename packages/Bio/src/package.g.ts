import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Returns an instance of the monomer library helper
//output: object result
export async function getMonomerLibHelper() : Promise<any> {
  return await PackageFunctions.getMonomerLibHelper();
}

//tags: init
//meta.role: init
export async function initBio() : Promise<void> {
  await PackageFunctions.initBio();
}

//tags: tooltip
//input: column col { semType: Macromolecule }
//output: widget result
//meta.role: tooltip
export function sequenceTooltip(col: DG.Column) : any {
  return PackageFunctions.sequenceTooltip(col);
}

//input: string library 
//output: string result
export async function standardiseMonomerLibrary(library: string) : Promise<string> {
  return await PackageFunctions.standardiseMonomerLibrary(library);
}

//description: Matches molecules in a column with monomers from the selected library(s)
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string polymerType = 'PEPTIDE' { choices: ["PEPTIDE","RNA","CHEM"]; caption: Polymer Type }
//top-menu: Bio | Manage | Match with Monomer Library...
export async function matchWithMonomerLibrary(table: DG.DataFrame, molecules: DG.Column, polymerType: any) : Promise<void> {
  await PackageFunctions.matchWithMonomerLibrary(table, molecules, polymerType);
}

//output: object monomerLib
export function getBioLib() : any {
  return PackageFunctions.getBioLib();
}

//input: column sequence { semType: Macromolecule }
//output: object result
export function getSeqHandler(sequence: DG.Column<any>) : any {
  return PackageFunctions.getSeqHandler(sequence);
}

//name: Bioinformatics | Get Region
//description: Creates a new column with sequences of the region between start and end
//tags: panel
//input: column seqCol { semType: Macromolecule }
//output: widget result
//meta.role: panel
export function getRegionPanel(seqCol: DG.Column<any>) : any {
  return PackageFunctions.getRegionPanel(seqCol);
}

//name: Bioinformatics | Manage Monomer Libraries
//tags: exclude-actions-panel
//input: column seqColumn { semType: Macromolecule }
//output: widget result
//meta.exclude-actions-panel: true
//meta.role: panel
export async function libraryPanel(_seqColumn: DG.Column) : Promise<any> {
  return await PackageFunctions.libraryPanel(_seqColumn);
}

//tags: editor
//input: funccall call 
//meta.role: editor
export function GetRegionEditor(call: DG.FuncCall) : void {
  PackageFunctions.GetRegionEditor(call);
}

//tags: editor
//input: funccall call 
//meta.role: editor
export function SplitToMonomersEditor(call: DG.FuncCall) : void {
  PackageFunctions.SplitToMonomersEditor(call);
}

//tags: editor
//input: funccall call 
//meta.role: editor
export function SequenceSpaceEditor(call: DG.FuncCall) : void {
  PackageFunctions.SequenceSpaceEditor(call);
}

//tags: editor
//input: funccall call 
//meta.role: editor
export function SeqActivityCliffsEditor(call: DG.FuncCall) : void {
  PackageFunctions.SeqActivityCliffsEditor(call);
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=custom
//meta.role: cellRenderer
export function customSequenceCellRenderer() : any {
  return PackageFunctions.customSequenceCellRenderer();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=fasta
//meta.role: cellRenderer
export function fastaSequenceCellRenderer() : any {
  return PackageFunctions.fastaSequenceCellRenderer();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=separator
//meta.role: cellRenderer
export function separatorSequenceCellRenderer() : any {
  return PackageFunctions.separatorSequenceCellRenderer();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=biln
//meta.role: cellRenderer
export function bilnSequenceCellRenderer() : any {
  return PackageFunctions.bilnSequenceCellRenderer();
}

//tags: notationRefiner
//input: column col 
//input: object stats 
//input: string separator { nullable: true; optional: true }
//output: bool result
//meta.role: notationRefiner
export function refineNotationProviderForBiln(col: DG.Column<any>, stats: any, separator: any) : boolean {
  return PackageFunctions.refineNotationProviderForBiln(col, stats, separator);
}

//name: Bioinformatics | Sequence Renderer
//tags: panel
//input: column molColumn { semType: Macromolecule }
//output: widget result
//meta.role: panel
export function macroMolColumnPropertyPanel(molColumn: DG.Column) : any {
  return PackageFunctions.macroMolColumnPropertyPanel(molColumn);
}

//name: Composition analysis
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export function compositionAnalysisWidget(sequence: DG.SemanticValue) : any {
  return PackageFunctions.compositionAnalysisWidget(sequence);
}

//name: MacromoleculeDifferenceCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: MacromoleculeDifference
//meta.columnTags: quality=MacromoleculeDifference
//meta.role: cellRenderer
export function macromoleculeDifferenceCellRenderer() : any {
  return PackageFunctions.macromoleculeDifferenceCellRenderer();
}

//input: string alignType { choices: ["Local alignment","Global alignment"] }
//input: string alignTable { choices: ["AUTO","NUCLEOTIDES","BLOSUM45","BLOSUM50","BLOSUM62","BLOSUM80","BLOSUM90","PAM30","PAM70","PAM250","SCHNEIDER","TRANS"] }
//input: double gap 
//input: string seq1 
//input: string seq2 
//output: object result
export function sequenceAlignment(alignType: string, alignTable: string, gap: number, seq1: string, seq2: string) {
  return PackageFunctions.sequenceAlignment(alignType, alignTable, gap, seq1, seq2);
}

//name: WebLogo
//description: WebLogo
//tags: viewer, panel
//output: viewer result
//meta.icon: files/icons/weblogo-viewer.svg
//meta.role: viewer,panel
export function webLogoViewer() {
  return PackageFunctions.webLogoViewer();
}

//name: VdRegions
//description: V-Domain regions viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/vdregions-viewer.svg
//meta.role: viewer,panel
export function vdRegionsViewer() {
  return PackageFunctions.vdRegionsViewer();
}

//description: Gets a new column with sequences of the region between start and end
//input: column sequence 
//input: string start { optional: true }
//input: string end { optional: true }
//input: string name { optional: true; description: Name of the column to be created }
//output: column result
export function getRegion(sequence: DG.Column<any>, start?: string, end?: string, name?: string) : any {
  return PackageFunctions.getRegion(sequence, start, end, name);
}

//name: Get Region Top Menu
//description: Get sequences for a region specified from a Macromolecule
//input: dataframe table { description: Input data table }
//input: column sequence { semType: Macromolecule; description: Sequence column }
//input: string start { optional: true; description: Region start position name }
//input: string end { optional: true; description: Region end position name }
//input: string name { optional: true; description: Region column name }
//top-menu: Bio | Calculate | Get Region...
//editor: Bio:GetRegionEditor
export async function getRegionTopMenu(table: DG.DataFrame, sequence: DG.Column, start?: string, end?: string, name?: string) : Promise<void> {
  await PackageFunctions.getRegionTopMenu(table, sequence, start, end, name);
}

//name: Sequence Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table { description: Input data table }
//input: string molecules { semType: Macromolecule; description: Input data table }
//input: column activities 
//input: double similarity = 80 { description: Similarity cutoff }
//input: string methodName { choices: ["UMAP","t-SNE"] }
//input: string similarityMetric { choices: ["Hamming","Levenshtein","Monomer chemical distance"] }
//input: func preprocessingFunction 
//input: object options { optional: true }
//input: bool demo { optional: true }
//top-menu: Bio | Analyze | Activity Cliffs...
//editor: Bio:SeqActivityCliffsEditor
export async function activityCliffs(table: DG.DataFrame, molecules: DG.Column<any>, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, preprocessingFunction: any, options?: any, demo?: boolean) : Promise<any> {
  return await PackageFunctions.activityCliffs(table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, demo);
}

//name: Encode Sequences
//tags: dim-red-preprocessing-function
//input: column col { semType: Macromolecule }
//input: string metric 
//input: double gapOpen = 1 { caption: Gap open penalty; optional: true }
//input: double gapExtend = 0.6 { caption: Gap extension penalty; optional: true }
//input: string fingerprintType = 'Morgan' { caption: Fingerprint type; choices: ["Morgan","RDKit","Pattern","AtomPair","MACCS","TopologicalTorsion"]; optional: true }
//output: object result
//meta.supportedSemTypes: Macromolecule
//meta.supportedTypes: string
//meta.supportedDistanceFunctions: Hamming,Levenshtein,Monomer chemical distance,Needlemann-Wunsch
//meta.role: dim-red-preprocessing-function
export async function macromoleculePreprocessingFunction(col: DG.Column, metric: any, gapOpen: number, gapExtend: number, fingerprintType: string) : Promise<any> {
  return await PackageFunctions.macromoleculePreprocessingFunction(col, metric, gapOpen, gapExtend, fingerprintType);
}

//name: Helm Fingerprints
//input: column col { semType: Macromolecule }
//input: string _metric 
//output: object result
//meta.supportedSemTypes: Macromolecule
//meta.supportedTypes: string
//meta.supportedUnits: helm
//meta.supportedDistanceFunctions: Tanimoto,Asymmetric,Cosine,Sokal
export async function helmPreprocessingFunction(col: DG.Column<any>, _metric: any) : Promise<any> {
  return await PackageFunctions.helmPreprocessingFunction(col, _metric);
}

//name: Sequence Space
//description: Creates 2D sequence space with projected sequences by pairwise distance
//input: dataframe table 
//input: column molecules { semType: Macromolecule }
//input: string methodName { choices: ["UMAP","t-SNE"] }
//input: string similarityMetric { choices: ["Hamming","Levenshtein","Monomer chemical distance"] }
//input: bool plotEmbeddings = true 
//input: func preprocessingFunction { optional: true }
//input: object options { optional: true }
//input: bool clusterEmbeddings = true { optional: true }
//input: bool isDemo { optional: true }
//top-menu: Bio | Analyze | Sequence Space...
//editor: Bio:SequenceSpaceEditor
export async function sequenceSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, preprocessingFunction?: any, options?: any, clusterEmbeddings?: boolean, isDemo?: boolean) : Promise<any> {
  return await PackageFunctions.sequenceSpaceTopMenu(table, molecules, methodName, similarityMetric, plotEmbeddings, preprocessingFunction, options, clusterEmbeddings, isDemo);
}

//name: Molecules to HELM
//description: Converts Peptide molecules to HELM notation by matching with monomer library
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule; description: Molecule column }
//top-menu: Bio | Transform | Molecules to HELM...
export async function moleculesToHelmTopMenu(table: DG.DataFrame, molecules: DG.Column) : Promise<void> {
  await PackageFunctions.moleculesToHelmTopMenu(table, molecules);
}

//name: Molecule to HELM Single
//description: Converts a single molecule to HELM notation without requiring a table or column
//input: string molecule { semType: Molecule; description: Input molecule }
//output: string result { semType: Macromolecule; units: helm }
export async function moleculeToHelmSingle(molecule: string) : Promise<string> {
  return await PackageFunctions.moleculeToHelmSingle(molecule);
}

//name: To Atomic Level
//description: Converts sequences to molblocks
//input: dataframe table { description: Input data table }
//input: column seqCol { semType: Macromolecule; caption: Sequence }
//input: bool nonlinear = true { caption: Non-linear; description: Slower mode for cycling/branching HELM structures }
//input: bool highlight = false { caption: Highlight monomers; description: Highlight monomers' substructures of the molecule }
//top-menu: Bio | Transform | To Atomic Level...
export async function toAtomicLevel(table: DG.DataFrame, seqCol: DG.Column, nonlinear: boolean, highlight: boolean) : Promise<void> {
  await PackageFunctions.toAtomicLevel(table, seqCol, nonlinear, highlight);
}

//name: To Atomic Level...
//input: column seqCol { semType: Macromolecule }
//meta.action: to atomic level
export async function toAtomicLevelAction(seqCol: DG.Column) : Promise<void> {
  await PackageFunctions.toAtomicLevelAction(seqCol);
}

//name: Molecular Structure
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export async function toAtomicLevelPanel(sequence: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.toAtomicLevelPanel(sequence);
}

//name: To Atomic Level Single sequence
//description: Converts a single sequence to molblock
//input: string sequence { semType: Macromolecule }
//output: string molfile { semType: Molecule }
export async function toAtomicLevelSingleSeq(sequence: string) : Promise<string> {
  return await PackageFunctions.toAtomicLevelSingleSeq(sequence);
}

//name: Molecular 3D Structure
//tags: bio, widgets, panel
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: bio
export async function sequence3dStructureWidget(sequence: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.sequence3dStructureWidget(sequence);
}

//name: MSA
//description: Performs multiple sequence alignment
//tags: bio, panel
//meta.domain: bio
//meta.role: panel
//top-menu: Bio | Analyze | MSA...
export function multipleSequenceAlignmentDialog() : void {
  PackageFunctions.multipleSequenceAlignmentDialog();
}

//name: Multiple Sequence Alignment
//description: Multiple sequence alignment
//input: column sequenceCol { semType: Macromolecule }
//input: column clustersCol 
//input: object options { optional: true }
//output: column result
//meta.domain: bio
export async function alignSequences(sequenceCol: any, clustersCol: any, options?: any) : Promise<any> {
  return await PackageFunctions.alignSequences(sequenceCol, clustersCol, options);
}

//name: Composition Analysis
//description: Visualizes sequence composition on a WebLogo plot
//output: viewer result
//meta.icon: files/icons/composition-analysis.svg
//top-menu: Bio | Analyze | Composition
export async function compositionAnalysis() {
  return await PackageFunctions.compositionAnalysis();
}

//description: Opens FASTA file
//tags: fileHandler
//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: fasta, fna, ffn, faa, frn, fa, fst
export function importFasta(fileContent: string) : any {
  return PackageFunctions.importFasta(fileContent);
}

//description: Opens Bam file
//tags: fileHandler
//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: bam, bai
export function importBam(fileContent: string) : any {
  return PackageFunctions.importBam(fileContent);
}

//top-menu: Bio | Transform | Convert Sequence Notation...
export function convertDialog() : void {
  PackageFunctions.convertDialog();
}

//name: Convert Notation...
//input: column col { semType: Macromolecule }
//meta.action: Convert Notation...
export function convertColumnAction(col: DG.Column) : void {
  PackageFunctions.convertColumnAction(col);
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: Monomer
//meta.columnTags: quality=Monomer
//meta.role: cellRenderer
export function monomerCellRenderer() : any {
  return PackageFunctions.monomerCellRenderer();
}

//input: string path { choices: ["Demo:Files/","System:AppData/"] }
//output: dataframe result
export async function testDetectMacromolecule(path: string) : Promise<any> {
  return await PackageFunctions.testDetectMacromolecule(path);
}

//name: Split to Monomers
//input: dataframe table 
//input: column sequence { semType: Macromolecule }
//output: dataframe result
//top-menu: Bio | Transform | Split to Monomers...
//editor: Bio:SplitToMonomersEditor
export async function splitToMonomersTopMenu(table: DG.DataFrame, sequence: DG.Column) : Promise<any> {
  return await PackageFunctions.splitToMonomersTopMenu(table, sequence);
}

//name: Bio: getHelmMonomers
//input: column sequence { semType: Macromolecule }
//output: object result
export function getHelmMonomers(sequence: DG.Column<any>) : string[] {
  return PackageFunctions.getHelmMonomers(sequence);
}

//name: Sequence Similarity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/sequence-similarity-viewer.svg
//meta.role: viewer
export function similaritySearchViewer() : any {
  return PackageFunctions.similaritySearchViewer();
}

//name: similaritySearch
//description: Finds similar sequences
//output: viewer result
//top-menu: Bio | Search | Similarity Search
export function similaritySearchTopMenu() {
  return PackageFunctions.similaritySearchTopMenu();
}

//name: Sequence Diversity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/sequence-diversity-viewer.svg
//meta.role: viewer
export function diversitySearchViewer() : any {
  return PackageFunctions.diversitySearchViewer();
}

//name: diversitySearch
//description: Finds the most diverse sequences
//output: viewer result
//top-menu: Bio | Search | Diversity Search
export function diversitySearchTopMenu() {
  return PackageFunctions.diversitySearchTopMenu();
}

//name: SearchSubsequenceEditor
//tags: editor
//input: funccall call 
//meta.role: editor
export function searchSubsequenceEditor(call: DG.FuncCall) : void {
  PackageFunctions.searchSubsequenceEditor(call);
}

//name: Subsequence Search
//input: column macromolecules 
//top-menu: Bio | Search | Subsequence Search ...
//editor: Bio:SearchSubsequenceEditor
export function SubsequenceSearchTopMenu(macromolecules: DG.Column) : void {
  PackageFunctions.SubsequenceSearchTopMenu(macromolecules);
}

//name: Identity
//description: Adds a column with fraction of matching monomers
//input: dataframe table { description: Table containing Macromolecule column }
//input: column macromolecule { semType: Macromolecule; description: Sequences to score }
//input: string reference { description: Sequence,matching column format }
//output: column result
//top-menu: Bio | Calculate | Identity...
export async function sequenceIdentityScoring(table: DG.DataFrame, macromolecule: DG.Column, reference: string) : Promise<any> {
  return await PackageFunctions.sequenceIdentityScoring(table, macromolecule, reference);
}

//name: Similarity
//description: Adds a column with similarity scores, calculated as sum of monomer fingerprint similarities
//input: dataframe table { description: Table containing Macromolecule column }
//input: column macromolecule { semType: Macromolecule; description: Sequences to score }
//input: string reference { description: Sequence,matching column format }
//output: column result
//top-menu: Bio | Calculate | Similarity...
export async function sequenceSimilarityScoring(table: DG.DataFrame, macromolecule: DG.Column, reference: string) : Promise<any> {
  return await PackageFunctions.sequenceSimilarityScoring(table, macromolecule, reference);
}

//name: Manage Monomer Libraries
//description: Manage HELM monomer libraries
export async function manageMonomerLibraries() : Promise<void> {
  await PackageFunctions.manageMonomerLibraries();
}

//name: Manage Monomer Libraries View
//top-menu: Bio | Manage | Monomer Libraries
export async function manageLibrariesView() : Promise<void> {
  await PackageFunctions.manageLibrariesView();
}

//description: Edit and create monomers
//top-menu: Bio | Manage | Monomers
export async function manageMonomersView() : Promise<void> {
  await PackageFunctions.manageMonomersView();
}

//name: Manage Monomer Libraries
//tags: app
//output: view result
//meta.role: app
//meta.browsePath: Peptides
//meta.icon: files/icons/monomers.png
export async function manageMonomerLibrariesView() : Promise<any> {
  return await PackageFunctions.manageMonomerLibrariesView();
}

//name: Monomer Manager Tree Browser
//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Manage Monomer Libraries
export async function manageMonomerLibrariesViewTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.manageMonomerLibrariesViewTreeBrowser(treeNode);
}

//description: As FASTA...
//meta.role: fileExporter
export function saveAsFasta() : void {
  PackageFunctions.saveAsFasta();
}

//name: Bio Substructure Filter
//description: Substructure filter for macromolecules
//tags: filter
//output: filter result
//meta.semType: Macromolecule
//meta.role: filter
export function bioSubstructureFilter() : any {
  return PackageFunctions.bioSubstructureFilter();
}

//name: Bio Substructure Filter Test
//description: Substructure filter for Helm package tests
//output: object result
export function bioSubstructureFilterTest() : any {
  return PackageFunctions.bioSubstructureFilterTest();
}

//name: webLogoLargeApp
export async function webLogoLargeApp() : Promise<void> {
  await PackageFunctions.webLogoLargeApp();
}

//name: webLogoAggApp
export async function webLogoAggApp() : Promise<void> {
  await PackageFunctions.webLogoAggApp();
}

//name: getRegionApp
export async function getRegionApp() : Promise<void> {
  await PackageFunctions.getRegionApp();
}

//name: getRegionHelmApp
export async function getRegionHelmApp() : Promise<void> {
  await PackageFunctions.getRegionHelmApp();
}

//name: longSeqTableSeparator
export function longSeqTableSeparator() : void {
  PackageFunctions.longSeqTableSeparator();
}

//name: longSeqTableFasta
export function longSeqTableFasta() : void {
  PackageFunctions.longSeqTableFasta();
}

//name: longSeqTableHelm
export function longSeqTableHelm() : void {
  PackageFunctions.longSeqTableHelm();
}

//input: object cell 
//input: object menu 
export function addCopyMenu(cell: any, menu: any) : void {
  PackageFunctions.addCopyMenu(cell, menu);
}

//description: Sequence similarity tracking and evaluation dataset diversity
//meta.demoPath: Bioinformatics | Similarity, Diversity
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Similarity,%20Diversity
export async function demoBioSimilarityDiversity() : Promise<void> {
  await PackageFunctions.demoBioSimilarityDiversity();
}

//description: Exploring sequence space of Macromolecules, comparison with hierarchical clustering results
//meta.isDemoDashboard: true
//meta.demoPath: Bioinformatics | Sequence Space
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Sequence%20Space
export async function demoBioSequenceSpace() : Promise<void> {
  await PackageFunctions.demoBioSequenceSpace();
}

//description: Activity Cliffs analysis on Macromolecules data
//meta.demoPath: Bioinformatics | Activity Cliffs
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Activity%20Cliffs
export async function demoBioActivityCliffs() : Promise<void> {
  await PackageFunctions.demoBioActivityCliffs();
}

//description: Atomic level structure of Macromolecules
//meta.demoPath: Bioinformatics | Atomic Level
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Atomic%20Level
export async function demoBioAtomicLevel() : Promise<void> {
  await PackageFunctions.demoBioAtomicLevel();
}

//name: SDF to JSON Library
//input: dataframe table 
export async function sdfToJsonLib(table: DG.DataFrame) : Promise<void> {
  await PackageFunctions.sdfToJsonLib(table);
}

//description: Converts a `Macromolecule` sequence to its atomic level `Molecule` representation
//input: string seq { semType: Macromolecule }
//input: bool nonlinear 
//output: string molfile { semType: Molecule }
//friendlyName: seq2atomic
export async function seq2atomic(seq: string, nonlinear: boolean) : Promise<any> {
  return await PackageFunctions.seq2atomic(seq, nonlinear);
}

//description: Gets identity to a reference sequence
//input: string seq { semType: Macromolecule }
//input: string ref { semType: Macromolecule }
//output: double result
//friendlyName: seqIdentity
export async function seqIdentity(seq: string, ref: string) : Promise<any> {
  return await PackageFunctions.seqIdentity(seq, ref);
}

//input: file file 
//input: string colName 
//input: double probeCount = 100 
export async function detectMacromoleculeProbe(file: DG.FileInfo, colName: string, probeCount: number) : Promise<void> {
  await PackageFunctions.detectMacromoleculeProbe(file, colName, probeCount);
}

//output: object result
export async function getSeqHelper() : Promise<any> {
  return await PackageFunctions.getSeqHelper();
}

//input: dataframe df 
//input: column helmCol 
//input: bool chiralityEngine = true 
//output: column result
export async function getMolFromHelm(df: DG.DataFrame, helmCol: DG.Column<any>, chiralityEngine: boolean) : Promise<any> {
  return await PackageFunctions.getMolFromHelm(df, helmCol, chiralityEngine);
}
