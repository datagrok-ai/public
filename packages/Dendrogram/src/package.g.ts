import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() {
  return PackageFunctions.info();
}

//name: Dendrogram
//description: Dendrogram tree visualization
//tags: viewer
//output: viewer result
//meta.icon: files/icons/dendrogram-viewer.svg
export function dendrogram() {
  return PackageFunctions.dendrogram();
}

//name: getTreeHelper
//output: object result
export function getTreeHelper() {
  return PackageFunctions.getTreeHelper();
}

//name: getDendrogramService
//output: object result
export function getDendrogramService() {
  return PackageFunctions.getDendrogramService();
}

//name: generateTreeDialog
export function generateTreeDialog() {
  return PackageFunctions.generateTreeDialog();
}

//name: dendrogramApp
//description: Test/demo app for Dendrogram
export async function dendrogramApp() {
  return PackageFunctions.dendrogramApp();
}

//name: dendrogramLargeApp
//description: Test/demo app for Dendrogram Large
export async function dendrogramLargeApp() {
  return PackageFunctions.dendrogramLargeApp();
}

//name: treeForGridApp
//description: Test/demo app for TreeForGrid (custom renderer)
export async function treeForGridApp() {
  return PackageFunctions.treeForGridApp();
}

//name: treeForGridFilterApp
//description: Test/demo app for TreeForGridFilter (custom renderer)
export async function treeForGridFilterApp() {
  return PackageFunctions.treeForGridFilterApp();
}

//name: treeForGridCutApp
//description: Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)
export async function treeForGridCutApp() {
  return PackageFunctions.treeForGridCutApp();
}

//name: hierarchicalClusteringApp
//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringApp() {
  return PackageFunctions.hierarchicalClusteringApp();
}

//name: hierarchicalClusteringSequencesApp
//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringSequencesApp() {
  return PackageFunctions.hierarchicalClusteringSequencesApp();
}

//name: Hierarchical Clustering
//description: Calculates hierarchical clustering on features and injects tree to grid
//input: dataframe df 
//input: column_list colNameList 
//input: string distance { default: euclidean; choices: ['euclidean','manhattan'] }
//input: string linkage { default: ward; choices: ['single','complete','average','weighted','centroid','median','ward'] }
export async function hierarchicalClustering(df: DG.DataFrame, colNameList: DG.ColumnList, distance: any, linkage: string) {
  return PackageFunctions.hierarchicalClustering(df, colNameList, distance, linkage);
}

//name: importNwk
//description: Opens Newick file
//tags: file-handler
//input: string fileContent 
//output: list tables
//meta.ext: nwk, newick
export async function importNewick(fileContent: string) {
  return PackageFunctions.importNewick(fileContent);
}

//name: previewNewick
//tags: fileViewer
//input: file file 
//output: view preview
//meta.fileViewer: nwk,newick
export async function previewNewick(file: DG.FileInfo) {
  return PackageFunctions.previewNewick(file);
}

//name: Hierarchical Clustering (Sequences)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: Bio | Analyze | Hierarchical Clustering...
export async function hierarchicalClusteringSequences() {
  return PackageFunctions.hierarchicalClusteringSequences();
}

//name: Hierarchical Clustering (Molecules)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: Chem | Analyze | Hierarchical Clustering...
export async function hierarchicalClusteringMolecules() {
  return PackageFunctions.hierarchicalClusteringMolecules();
}

//name: Hierarchical Clustering (All)
//description: Calculates hierarchical clustering on features and injects tree to grid
//top-menu: ML | Cluster | Hierarchical Clustering...
export async function hierarchicalClustering2() {
  return PackageFunctions.hierarchicalClustering2();
}

//name: heatMapDemo
//description: Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
//meta.demoPath: Visualization | General | Heatmap
export async function _heatMapDemo() {
  return PackageFunctions._heatMapDemo();
}
