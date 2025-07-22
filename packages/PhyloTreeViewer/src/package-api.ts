import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  //Phylogenetic tree visualization
  export async function phylocanvasGlViewer(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:PhylocanvasGlViewer', {});
  }

  //Test/demo app for PhylocanvasGlViewer
  export async function phylocanvasGlViewerApp(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:PhylocanvasGlViewerApp', {});
  }

  //Test/demo app for TreeToGrid (PhylocanvasGL based)
  export async function treeToGridApp(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:TreeToGridApp', {});
  }

  //Test/demo app for TreeCutAsTree
  export async function treeCutAsTreeApp(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:TreeCutAsTreeApp', {});
  }

  //Test/demo app for TreeInGridCell
  export async function treeInGridCellApp(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:TreeInGridCellApp', {});
  }

  //Opens Newick file
  export async function injectTreeToGrid(grid: any, newickText: string, leafColName: string): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:InjectTreeToGrid', { grid, newickText, leafColName });
  }

  export async function getPhylocanvasGlService(): Promise<any> {
    return await grok.functions.call('PhyloTreeViewer:GetPhylocanvasGlService', {});
  }
}
