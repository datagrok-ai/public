import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getMonomerLibHelper
//description: Returns an instance of the monomer library helper
//output: object result
export async function getMonomerLibHelper() {
  return PackageFunctions.getMonomerLibHelper();
}

//name: initBio
//tags: init
export async function initBio() {
  return PackageFunctions.initBio();
}

//name: sequenceTooltip
//tags: tooltip
//input: column col { semType: Macromolecule }
//output: widget result
export function sequenceTooltip(col: DG.Column) {
  return PackageFunctions.sequenceTooltip(col);
}

//name: getBioLib
//output: object monomerLib
export function getBioLib() {
  return PackageFunctions.getBioLib();
}

//name: getSeqHandler
//input: column sequence { semType: Macromolecule }
//output: object result
export function getSeqHandler(sequence: any) {
  return PackageFunctions.getSeqHandler(sequence);
}

//name: Bioinformatics | Get Region
//description: Creates a new column with sequences of the region between start and end
//tags: panel
//input: column seqCol { semType: Macromolecule }
//output: widget result
export function getRegionPanel(seqCol: any) {
  return PackageFunctions.getRegionPanel(seqCol);
}

//name: Bioinformatics | Manage Monomer Libraries
//tags: panel, exclude-actions-panel
//input: column seqColumn { semType: Macromolecule }
//output: widget result
export async function libraryPanel(_seqColumn: DG.Column) {
  return PackageFunctions.libraryPanel(_seqColumn);
}

//name: GetRegionEditor
//tags: editor
//input: funccall call 
export function GetRegionEditor(call: DG.FuncCall) {
  return PackageFunctions.GetRegionEditor(call);
}

//name: SplitToMonomersEditor
//tags: editor
//input: funccall call 
export function SplitToMonomersEditor(call: DG.FuncCall) {
  return PackageFunctions.SplitToMonomersEditor(call);
}

//name: SequenceSpaceEditor
//tags: editor
//input: funccall call 
export function SequenceSpaceEditor(call: DG.FuncCall) {
  return PackageFunctions.SequenceSpaceEditor(call);
}

//name: SeqActivityCliffsEditor
//tags: editor
//input: funccall call 
export function SeqActivityCliffsEditor(call: DG.FuncCall) {
  return PackageFunctions.SeqActivityCliffsEditor(call);
}

//name: customSequenceCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=custom
export function customSequenceCellRenderer() {
  return PackageFunctions.customSequenceCellRenderer();
}

//name: fastaSequenceCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=fasta
export function fastaSequenceCellRenderer() {
  return PackageFunctions.fastaSequenceCellRenderer();
}

//name: separatorSequenceCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sequence
//meta.columnTags: quality=Macromolecule, units=separator
export function separatorSequenceCellRenderer() {
  return PackageFunctions.separatorSequenceCellRenderer();
}

//name: Bioinformatics | Sequence Renderer
//tags: panel
//input: column molColumn { semType: Macromolecule }
//output: widget result
export function macroMolColumnPropertyPanel(molColumn: DG.Column) {
  return PackageFunctions.macroMolColumnPropertyPanel(molColumn);
}

//name: Composition analysis
//tags: panel, bio, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export function compositionAnalysisWidget(sequence: DG.SemanticValue) {
  return PackageFunctions.compositionAnalysisWidget(sequence);
}

//name: MacromoleculeDifferenceCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: MacromoleculeDifference
//meta.columnTags: quality=MacromoleculeDifference
export function macromoleculeDifferenceCellRenderer() {
  return PackageFunctions.macromoleculeDifferenceCellRenderer();
}

//name: sequenceAlignment
//input: string alignType { choices: ['Local alignment','Global alignment'] }
//input: string alignTable { choices: ['AUTO','NUCLEOTIDES','BLOSUM45','BLOSUM50','BLOSUM62','BLOSUM80','BLOSUM90','PAM30','PAM70','PAM250','SCHNEIDER','TRANS'] }
//input: double gap 
//input: string seq1 
//input: string seq2 
//output: object result
export function sequenceAlignment(alignType: string, alignTable: string, gap: number, seq1: string, seq2: string) {
  return PackageFunctions.sequenceAlignment(alignType, alignTable, gap, seq1, seq2);
}

//name: WebLogo
//description: WebLogo
//tags: panel, viewer
//output: viewer result
//meta.icon: files/icons/weblogo-viewer.svg
export function webLogoViewer() {
  return PackageFunctions.webLogoViewer();
}

//name: VdRegions
//description: V-Domain regions viewer
//tags: panel, viewer
//output: viewer result
//meta.icon: files/icons/vdregions-viewer.svg
export function vdRegionsViewer() {
  return PackageFunctions.vdRegionsViewer();
}

//name: getRegion
//description: Gets a new column with sequences of the region between start and end
//input: column sequence 
//input: string start { optional: true }
//input: string end { optional: true }
//input: string name { optional: true; description: Name of the column to be created }
//output: column result
export function getRegion(sequence: any, start: any, end: any, name: any) {
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
export async function getRegionTopMenu(table: DG.DataFrame, sequence: DG.Column, start: any, end: any, name: any) {
  return PackageFunctions.getRegionTopMenu(table, sequence, start, end, name);
}

//name: Sequence Activity Cliffs
//description: Detects pairs of molecules with similar structure and significant difference in any given property
//input: dataframe table 
//input: string molecules { semType: Macromolecule; description: Input data table }
//input: column activities 
//input: double similarity { default: 80; description: Similarity cutoff }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Hamming','Levenshtein','Monomer chemical distance'] }
//input: func preprocessingFunction 
//input: object options { optional: true }
//input: bool demo { optional: true }
//top-menu: Bio | Analyze | Activity Cliffs...
//editor: Bio:SeqActivityCliffsEditor
export async function activityCliffs(table: DG.DataFrame, molecules: any, activities: DG.Column, similarity: number, methodName: any, similarityMetric: any, preprocessingFunction: any, options?: any, demo?: boolean) {
  return PackageFunctions.activityCliffs(table, molecules, activities, similarity, methodName, similarityMetric, preprocessingFunction, options, demo);
}

//name: Encode Sequences
//tags: dim-red-preprocessing-function
//input: column col { semType: Macromolecule }
//input: string metric 
//input: double gapOpen { default: 1; caption: Gap open penalty; optional: true }
//input: double gapExtend { default: 0.6; caption: Gap extension penalty; optional: true }
//input: string fingerprintType { caption: Fingerprint type; default: Morgan; choices: ['Morgan','RDKit','Pattern','AtomPair','MACCS','TopologicalTorsion']; optional: true }
//output: object result
//meta.supportedSemTypes: Macromolecule
//meta.supportedTypes: string
//meta.supportedDistanceFunctions: Hamming,Levenshtein,Monomer chemical distance,Needlemann-Wunsch
export async function macromoleculePreprocessingFunction(col: DG.Column, metric: any, gapOpen: number, gapExtend: number, fingerprintType: string) {
  return PackageFunctions.macromoleculePreprocessingFunction(col, metric, gapOpen, gapExtend, fingerprintType);
}

//name: Helm Fingerprints
//input: column col { semType: Macromolecule }
//input: string _metric 
//output: object result
//meta.supportedSemTypes: Macromolecule
//meta.supportedTypes: string
//meta.supportedUnits: helm
//meta.supportedDistanceFunctions: Tanimoto,Asymmetric,Cosine,Sokal
export async function helmPreprocessingFunction(col: any, _metric: any) {
  return PackageFunctions.helmPreprocessingFunction(col, _metric);
}

//name: Sequence Space
//description: Creates 2D sequence space with projected sequences by pairwise distance
//input: dataframe table 
//input: column molecules { semType: Macromolecule }
//input: string methodName { choices: ['UMAP','t-SNE'] }
//input: string similarityMetric { choices: ['Hamming','Levenshtein','Monomer chemical distance'] }
//input: bool plotEmbeddings { default: true }
//input: func preprocessingFunction { optional: true }
//input: object options { optional: true }
//input: bool clusterEmbeddings { optional: true; default: true }
//input: bool isDemo { optional: true }
//top-menu: Bio | Analyze | Sequence Space...
//editor: Bio:SequenceSpaceEditor
export async function sequenceSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: any, similarityMetric: any, plotEmbeddings: boolean, preprocessingFunction?: any, options?: any, clusterEmbeddings?: boolean, isDemo?: boolean) {
  return PackageFunctions.sequenceSpaceTopMenu(table, molecules, methodName, similarityMetric, plotEmbeddings, preprocessingFunction, options, clusterEmbeddings, isDemo);
}

//name: To Atomic Level
//description: Converts sequences to molblocks
//input: dataframe table { description: Input data table }
//input: column seqCol { semType: Macromolecule; caption: Sequence }
//input: bool nonlinear { default: false; caption: Non-linear; description: Slower mode for cycling/branching HELM structures }
//input: bool highlight { default: false; caption: Highlight monomers; description: Highlight monomers' substructures of the molecule }
//top-menu: Bio | Transform | To Atomic Level...
export async function toAtomicLevel(table: DG.DataFrame, seqCol: DG.Column, nonlinear: boolean, highlight: boolean) {
  return PackageFunctions.toAtomicLevel(table, seqCol, nonlinear, highlight);
}

//name: To Atomic Level...
//input: column seqCol { semType: Macromolecule }
//meta.action: to atomic level
export async function toAtomicLevelAction(seqCol: DG.Column) {
  return PackageFunctions.toAtomicLevelAction(seqCol);
}

//name: Molecular Structure
//tags: panel, bio, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export async function toAtomicLevelPanel(sequence: DG.SemanticValue) {
  return PackageFunctions.toAtomicLevelPanel(sequence);
}

//name: Molecular 3D Structure
//tags: panel, bio, widgets
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export async function sequence3dStructureWidget(sequence: DG.SemanticValue) {
  return PackageFunctions.sequence3dStructureWidget(sequence);
}

//name: MSA
//description: Performs multiple sequence alignment
//tags: panel, bio
//top-menu: Bio | Analyze | MSA...
export function multipleSequenceAlignmentDialog() {
  return PackageFunctions.multipleSequenceAlignmentDialog();
}

//name: Multiple Sequence Alignment
//description: Multiple sequence alignment
//tags: bio
//input: column sequenceCol { semType: Macromolecule }
//input: column clustersCol 
//input: object options { optional: true; default: undefined }
//output: column result
export async function alignSequences(sequenceCol: any, clustersCol: any, options?: any) {
  return PackageFunctions.alignSequences(sequenceCol, clustersCol, options);
}

//name: Composition Analysis
//description: Visualizes sequence composition on a WebLogo plot
//output: viewer result
//meta.icon: files/icons/composition-analysis.svg
//top-menu: Bio | Analyze | Composition
export async function compositionAnalysis() {
  return PackageFunctions.compositionAnalysis();
}

//name: importFasta
//description: Opens FASTA file
//tags: file-handler
//input: string fileContent 
//output: list<dataframe> result
//meta.ext: fasta, fna, ffn, faa, frn, fa, fst
export function importFasta(fileContent: string) {
  return PackageFunctions.importFasta(fileContent);
}

//name: importBam
//description: Opens Bam file
//tags: file-handler
//input: string fileContent 
//output: list<dataframe> result
//meta.ext: bam, bai
export function importBam(fileContent: string) {
  return PackageFunctions.importBam(fileContent);
}

//name: convertDialog
//top-menu: Bio | Transform | Convert Notation...
export function convertDialog() {
  return PackageFunctions.convertDialog();
}

//name: Convert Notation...
//input: column col { semType: Macromolecule }
//meta.action: Convert Notation...
export function convertColumnAction(col: DG.Column) {
  return PackageFunctions.convertColumnAction(col);
}

//name: monomerCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: Monomer
//meta.columnTags: quality=Monomer
export function monomerCellRenderer() {
  return PackageFunctions.monomerCellRenderer();
}

//name: testDetectMacromolecule
//input: string path { choices: ['Demo:Files/','System:AppData/'] }
//output: dataframe result
export async function testDetectMacromolecule(path: string) {
  return PackageFunctions.testDetectMacromolecule(path);
}

//name: Split to Monomers
//input: dataframe table 
//input: column sequence { semType: Macromolecule }
//output: dataframe result
//top-menu: Bio | Transform | Split to Monomers...
//editor: Bio:SplitToMonomersEditor
export async function splitToMonomersTopMenu(table: DG.DataFrame, sequence: DG.Column) {
  return PackageFunctions.splitToMonomersTopMenu(table, sequence);
}

//name: Bio: getHelmMonomers
//input: column sequence { semType: Macromolecule }
//output: object result
export function getHelmMonomers(sequence: any) {
  return PackageFunctions.getHelmMonomers(sequence);
}

//name: Sequence Similarity Search
//tags: viewer
//output: viewer result
//meta.icon: files/icons/sequence-similarity-viewer.svg
export function similaritySearchViewer() {
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
export function diversitySearchViewer() {
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
export function searchSubsequenceEditor(call: DG.FuncCall) {
  return PackageFunctions.searchSubsequenceEditor(call);
}

//name: Subsequence Search
//input: column macromolecules 
//top-menu: Bio | Search | Subsequence Search ...
//editor: Bio:SearchSubsequenceEditor
export function SubsequenceSearchTopMenu(macromolecules: DG.Column) {
  return PackageFunctions.SubsequenceSearchTopMenu(macromolecules);
}

//name: Identity Scoring
//description: Adds a column with fraction of matching monomers
//input: dataframe table { description: Table containing Macromolecule column }
//input: column macromolecule { semType: Macromolecule; description: Sequences to score }
//input: string reference { description: Sequence,matching column format }
//output: column result
//top-menu: Bio | Calculate | Identity...
export async function sequenceIdentityScoring(table: DG.DataFrame, macromolecule: DG.Column, reference: string) {
  return PackageFunctions.sequenceIdentityScoring(table, macromolecule, reference);
}

//name: Similarity Scoring
//description: Adds a column with similarity scores, calculated as sum of monomer fingerprint similarities
//input: dataframe table { description: Table containing Macromolecule column }
//input: column macromolecule { semType: Macromolecule; description: Sequences to score }
//input: string reference { description: Sequence,matching column format }
//output: column result
//top-menu: Bio | Calculate | Similarity...
export async function sequenceSimilarityScoring(table: DG.DataFrame, macromolecule: DG.Column, reference: string) {
  return PackageFunctions.sequenceSimilarityScoring(table, macromolecule, reference);
}

//name: Manage Monomer Libraries
//description: Manage HELM monomer libraries
export async function manageMonomerLibraries() {
  return PackageFunctions.manageMonomerLibraries();
}

//name: Manage Monomer Libraries View
//top-menu: Bio | Manage | Monomer Libraries
export async function manageLibrariesView() {
  return PackageFunctions.manageLibrariesView();
}

//name: manageMonomersView
//description: Edit and create monomers
//top-menu: Bio | Manage | Monomers
export async function manageMonomersView() {
  return PackageFunctions.manageMonomersView();
}

//name: Manage Monomer Libraries
//tags: app
//output: view result
//meta.browsePath: Peptides
//meta.icon: files/icons/monomers.png
export async function manageMonomerLibrariesView() {
  return PackageFunctions.manageMonomerLibrariesView();
}

//name: Monomer Manager Tree Browser
//input: dynamic treeNode 
//input: dynamic browsePanel 
export async function manageMonomerLibrariesViewTreeBrowser(treeNode: any, browsePanel: any) {
  return PackageFunctions.manageMonomerLibrariesViewTreeBrowser(treeNode, browsePanel);
}

//name: saveAsFasta
//description: As FASTA...
//tags: fileExporter
export function saveAsFasta() {
  return PackageFunctions.saveAsFasta();
}

//name: Bio Substructure Filter
//description: Substructure filter for macromolecules
//tags: filter
//output: filter result
//meta.semType: Macromolecule
export function bioSubstructureFilter() {
  return PackageFunctions.bioSubstructureFilter();
}

//name: Bio Substructure Filter Test
//description: Substructure filter for Helm package tests
//output: object result
export function bioSubstructureFilterTest() {
  return PackageFunctions.bioSubstructureFilterTest();
}

//name: webLogoLargeApp
export async function webLogoLargeApp() {
  return PackageFunctions.webLogoLargeApp();
}

//name: webLogoAggApp
export async function webLogoAggApp() {
  return PackageFunctions.webLogoAggApp();
}

//name: getRegionApp
export async function getRegionApp() {
  return PackageFunctions.getRegionApp();
}

//name: getRegionHelmApp
export async function getRegionHelmApp() {
  return PackageFunctions.getRegionHelmApp();
}

//name: longSeqTableSeparator
export function longSeqTableSeparator() {
  return PackageFunctions.longSeqTableSeparator();
}

//name: longSeqTableFasta
export function longSeqTableFasta() {
  return PackageFunctions.longSeqTableFasta();
}

//name: longSeqTableHelm
export function longSeqTableHelm() {
  return PackageFunctions.longSeqTableHelm();
}

//name: addCopyMenu
//input: object cell 
//input: object menu 
export function addCopyMenu(cell: any, menu: any) {
  return PackageFunctions.addCopyMenu(cell, menu);
}

//name: demoBioSimilarityDiversity
//description: Sequence similarity tracking and evaluation dataset diversity
//meta.demoPath: Bioinformatics | Similarity, Diversity
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Similarity,%20Diversity
//meta.demoSkip: GROK-14320
export async function demoBioSimilarityDiversity() {
  return PackageFunctions.demoBioSimilarityDiversity();
}

//name: demoBioSequenceSpace
//description: Exploring sequence space of Macromolecules, comparison with hierarchical clustering results
//meta.demoPath: Bioinformatics | Sequence Space
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Sequence%20Space
//meta.demoSkip: GROK-14320
export async function demoBioSequenceSpace() {
  return PackageFunctions.demoBioSequenceSpace();
}

//name: demoBioActivityCliffs
//description: Activity Cliffs analysis on Macromolecules data
//meta.demoPath: Bioinformatics | Activity Cliffs
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Activity%20Cliffs
//meta.demoSkip: GROK-14320
export async function demoBioActivityCliffs() {
  return PackageFunctions.demoBioActivityCliffs();
}

//name: demoBioAtomicLevel
//description: Atomic level structure of Macromolecules
//meta.demoPath: Bioinformatics | Atomic Level
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Atomic%20Level
//meta.demoSkip: GROK-14320
export async function demoBioAtomicLevel() {
  return PackageFunctions.demoBioAtomicLevel();
}

//name: SDF to JSON Library
//input: dataframe table 
export async function sdfToJsonLib(table: DG.DataFrame) {
  return PackageFunctions.sdfToJsonLib(table);
}

//name: seq2atomic
//description: Converts a `Macromolecule` sequence to its atomic level `Molecule` representation
//input: string seq { semType: Macromolecule }
//input: bool nonlinear 
//output: string molfile { semType: Molecule }
//friendlyName: seq2atomic
export async function seq2atomic(seq: string, nonlinear: boolean) {
  return PackageFunctions.seq2atomic(seq, nonlinear);
}

//name: seqIdentity
//description: Gets identity to a reference sequence
//input: string seq { semType: Macromolecule }
//input: string ref { semType: Macromolecule }
//output: double result
//friendlyName: seqIdentity
export async function seqIdentity(seq: string, ref: string) {
  return PackageFunctions.seqIdentity(seq, ref);
}

//name: detectMacromoleculeProbe
//input: file file 
//input: string colName { default:  }
//input: double probeCount { default: 100 }
export async function detectMacromoleculeProbe(file: DG.FileInfo, colName: string, probeCount: number) {
  return PackageFunctions.detectMacromoleculeProbe(file, colName, probeCount);
}

//name: getSeqHelper
//output: object result
export async function getSeqHelper() {
  return PackageFunctions.getSeqHelper();
}

//name: getMolFromHelm
//input: dataframe df 
//input: column helmCol 
//input: bool chiralityEngine { default: true }
//output: column result
export async function getMolFromHelm(df: DG.DataFrame, helmCol: any, chiralityEngine: boolean) {
  return PackageFunctions.getMolFromHelm(df, helmCol, chiralityEngine);
}
