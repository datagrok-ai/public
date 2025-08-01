import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';
//name: info
//output: dynamic result
//meta.undefined: info
export function info() {
  return PackageFunctions.info();
}

//name: dendrogram
//output: dynamic result
//meta.undefined: Dendrogram tree visualization
export function dendrogram() {
  return PackageFunctions.dendrogram();
}

//name: getTreeHelper
//output: dynamic result
//meta.undefined: getTreeHelper
export function getTreeHelper() {
  return PackageFunctions.getTreeHelper();
}

//name: getDendrogramService
//output: dynamic result
//meta.undefined: getDendrogramService
export function getDendrogramService() {
  return PackageFunctions.getDendrogramService();
}

//name: generateTreeDialog
//output: dynamic result
//meta.undefined: generateTreeDialog
export function generateTreeDialog() {
  return PackageFunctions.generateTreeDialog();
}

//name: dendrogramApp
//meta.undefined: Test/demo app for Dendrogram
export async function dendrogramApp() {
  return PackageFunctions.dendrogramApp();
}

//name: dendrogramLargeApp
//meta.undefined: Test/demo app for Dendrogram Large
export async function dendrogramLargeApp() {
  return PackageFunctions.dendrogramLargeApp();
}

//name: treeForGridApp
//meta.undefined: Test/demo app for TreeForGrid (custom renderer)
export async function treeForGridApp() {
  return PackageFunctions.treeForGridApp();
}

//name: treeForGridFilterApp
//meta.undefined: Test/demo app for TreeForGridFilter (custom renderer)
export async function treeForGridFilterApp() {
  return PackageFunctions.treeForGridFilterApp();
}

//name: treeForGridCutApp
//meta.undefined: Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)
export async function treeForGridCutApp() {
  return PackageFunctions.treeForGridCutApp();
}

//name: hierarchicalClusteringApp
//meta.undefined: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringApp() {
  return PackageFunctions.hierarchicalClusteringApp();
}

//name: hierarchicalClusteringSequencesApp
//meta.undefined: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringSequencesApp() {
  return PackageFunctions.hierarchicalClusteringSequencesApp();
}

//name: hierarchicalClustering
//input: dataframe df 
//input: column_list colNameList 
//input: dynamic distance 
//input: string linkage 
//meta.undefined: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClustering(df: DG.DataFrame, colNameList: DG.ColumnList, distance: any, linkage: string) {
  return PackageFunctions.hierarchicalClustering(df, colNameList, distance, linkage);
}

//name: importNewick
//input: string fileContent 
//output: list result
//meta.undefined: Opens Newick file
export async function importNewick(fileContent: string) {
  return PackageFunctions.importNewick(fileContent);
}

//name: previewNewick
//input: file file 
//output: dynamic result
//meta.undefined: fileViewer
export async function previewNewick(file: DG.FileInfo) {
  return PackageFunctions.previewNewick(file);
}

//name: hierarchicalClusteringSequences
//meta.undefined: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClusteringSequences() {
  return PackageFunctions.hierarchicalClusteringSequences();
}

//name: hierarchicalClusteringMolecules
//meta.undefined: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClusteringMolecules() {
  return PackageFunctions.hierarchicalClusteringMolecules();
}

//name: hierarchicalClustering2
//meta.undefined: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClustering2() {
  return PackageFunctions.hierarchicalClustering2();
}

//name: _heatMapDemo
//output: dynamic result
//meta.undefined: Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
export async function _heatMapDemo() {
  return PackageFunctions._heatMapDemo();
}

