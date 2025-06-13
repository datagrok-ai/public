import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  //Phylogenetic tree visualization
  export async function phylocanvasGlViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:PhylocanvasGlViewer', {});
  }

  //Test/demo app for PhylocanvasGlViewer
  export async function phylocanvasGlViewerApp(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:PhylocanvasGlViewerApp', {});
  }

  //Test/demo app for TreeToGrid (PhylocanvasGL based)
  export async function treeToGridApp(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:TreeToGridApp', {});
  }

  //Test/demo app for TreeCutAsTree
  export async function treeCutAsTreeApp(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:TreeCutAsTreeApp', {});
  }

  //Test/demo app for TreeInGridCell
  export async function treeInGridCellApp(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:TreeInGridCellApp', {});
  }

  //Opens Newick file
  export async function injectTreeToGrid(grid: any, newickText: string, leafColName: string): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:InjectTreeToGrid', { grid, newickText, leafColName });
  }

  export async function getPhylocanvasGlService(): Promise<any> {
    return await grok.functions.call('@datagrok/phylo-tree-viewer:GetPhylocanvasGlService', {});
  }
}
