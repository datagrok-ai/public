import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Dendrogram
//description: Dendrogram tree visualization
//output: viewer result
//meta.icon: files/icons/dendrogram-viewer.svg
//meta.role: viewer
export function dendrogram() : any {
  return PackageFunctions.dendrogram();
}

//output: object result
export function getTreeHelper() : any {
  return PackageFunctions.getTreeHelper();
}

//output: object result
export function getDendrogramService() : any {
  return PackageFunctions.getDendrogramService();
}

//name: generateTreeDialog
export function generateTreeDialog() : void {
  PackageFunctions.generateTreeDialog();
}

//description: Test/demo app for Dendrogram
export async function dendrogramApp() : Promise<void> {
  await PackageFunctions.dendrogramApp();
}

//description: Test/demo app for Dendrogram Large
export async function dendrogramLargeApp() : Promise<void> {
  await PackageFunctions.dendrogramLargeApp();
}

//description: Test/demo app for TreeForGrid (custom renderer)
export async function treeForGridApp() : Promise<void> {
  await PackageFunctions.treeForGridApp();
}

//description: Test/demo app for TreeForGridFilter (custom renderer)
export async function treeForGridFilterApp() : Promise<void> {
  await PackageFunctions.treeForGridFilterApp();
}

//description: Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)
export async function treeForGridCutApp() : Promise<void> {
  await PackageFunctions.treeForGridCutApp();
}

//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringApp() : Promise<void> {
  await PackageFunctions.hierarchicalClusteringApp();
}

//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringSequencesApp() : Promise<void> {
  await PackageFunctions.hierarchicalClusteringSequencesApp();
}

//name: Hierarchical Clustering
//description: Calculates hierarchical clustering on features and injects tree to grid
//input: dataframe df 
//input: column_list colNameList 
//input: string distance = 'euclidean' { choices: ["euclidean","manhattan"] }
//input: string linkage = 'ward' { choices: ["single","complete","average","weighted","centroid","median","ward"] }
export async function hierarchicalClustering(df: DG.DataFrame, colNameList: DG.ColumnList, distance: any, linkage: string) : Promise<void> {
  await PackageFunctions.hierarchicalClustering(df, colNameList, distance, linkage);
}

//name: importNwk
//description: Opens Newick file
//input: string fileContent 
//output: list tables
//meta.role: fileHandler
//meta.ext: nwk, newick
export async function importNewick(fileContent: string) : Promise<any> {
  return await PackageFunctions.importNewick(fileContent);
}

//input: file file 
//output: view preview
//meta.role: fileViewer
//meta.fileViewer: nwk,newick
export async function previewNewick(file: DG.FileInfo) {
  return await PackageFunctions.previewNewick(file);
}

//name: Hierarchical Clustering (Sequences)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: Bio | Analyze | Hierarchical Clustering...
export async function hierarchicalClusteringSequences() : Promise<void> {
  await PackageFunctions.hierarchicalClusteringSequences();
}

//name: Hierarchical Clustering (Molecules)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: Chem | Analyze | Hierarchical Clustering...
export async function hierarchicalClusteringMolecules() : Promise<void> {
  await PackageFunctions.hierarchicalClusteringMolecules();
}

//name: Hierarchical Clustering (All)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: ML | Cluster | Hierarchical Clustering...
export async function hierarchicalClustering2() : Promise<void> {
  await PackageFunctions.hierarchicalClustering2();
}

//name: heatMapDemo
//description: Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
//meta.demoPath: Visualization | General | Heatmap
export async function _heatMapDemo() : Promise<void> {
  await PackageFunctions._heatMapDemo();
}
