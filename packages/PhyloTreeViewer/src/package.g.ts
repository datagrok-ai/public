import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: PhylocanvasGL
//description: Phylogenetic tree visualization
//tags: viewer
//output: dynamic result
//meta.icon: files/icons/phylocanvasgl-viewer.svg
export function phylocanvasGlViewer() {
  return PackageFunctions.phylocanvasGlViewer();
}

//name: phylocanvasGlViewerApp
//description: Test/demo app for PhylocanvasGlViewer
export async function phylocanvasGlViewerApp() {
  return PackageFunctions.phylocanvasGlViewerApp();
}

//name: TreeToGrid
//description: Test/demo app for TreeToGrid (PhylocanvasGL based)
//output: dynamic result
export async function treeToGridApp() {
  return PackageFunctions.treeToGridApp();
}

//name: TreeCutAsTree
//description: Test/demo app for TreeCutAsTree
export async function treeCutAsTreeApp() {
  return PackageFunctions.treeCutAsTreeApp();
}

//name: TreeInGridCell
//description: Test/demo app for TreeInGridCell
export async function treeInGridCellApp() {
  return PackageFunctions.treeInGridCellApp();
}

//name: injectTree
//description: Opens Newick file
//input: viewer grid 
//input: string newickText 
//input: string leafColName { optional: true }
//output: dynamic result
export async function injectTreeToGrid(grid: any, newickText: string, leafColName: string) {
  return PackageFunctions.injectTreeToGrid(grid, newickText, leafColName);
}

//name: getPhylocanvasGlService
//output: dynamic result
export function getPhylocanvasGlService() {
  return PackageFunctions.getPhylocanvasGlService();
}
