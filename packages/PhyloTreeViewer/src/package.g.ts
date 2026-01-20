import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: PhylocanvasGL
//description: Phylogenetic tree visualization
//output: viewer result
//meta.icon: files/icons/phylocanvasgl-viewer.svg
//meta.role: viewer
export function phylocanvasGlViewer() : any {
  return PackageFunctions.phylocanvasGlViewer();
}

//description: Test/demo app for PhylocanvasGlViewer
export async function phylocanvasGlViewerApp() : Promise<void> {
  await PackageFunctions.phylocanvasGlViewerApp();
}

//name: TreeToGrid
//description: Test/demo app for TreeToGrid (PhylocanvasGL based)
//output: object result
export async function treeToGridApp() : Promise<any> {
  return await PackageFunctions.treeToGridApp();
}

//name: TreeCutAsTree
//description: Test/demo app for TreeCutAsTree
export async function treeCutAsTreeApp() : Promise<void> {
  await PackageFunctions.treeCutAsTreeApp();
}

//name: TreeInGridCell
//description: Test/demo app for TreeInGridCell
export async function treeInGridCellApp() : Promise<void> {
  await PackageFunctions.treeInGridCellApp();
}

//name: injectTree
//description: Opens Newick file
//input: viewer grid 
//input: string newickText 
//input: string leafColName { optional: true }
export async function injectTreeToGrid(grid: any, newickText: string, leafColName?: string) : Promise<void> {
  await PackageFunctions.injectTreeToGrid(grid, newickText, leafColName);
}

//output: object result
export function getPhylocanvasGlService() : any {
  return PackageFunctions.getPhylocanvasGlService();
}
