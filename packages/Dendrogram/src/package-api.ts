import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  //Returns distance matrix condensed
  export async function distanceScript(data: DG.DataFrame, distance_name: string): Promise<DG.DataFrame> {
    return await grok.functions.call('Dendrogram:DistanceScript', { data, distance_name });
  }

  //Returns the newick representation of the tree for given dataset
  export async function hierarchicalClusteringScript(data: DG.DataFrame, distance_name: string, linkage_name: string): Promise<string> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringScript', { data, distance_name, linkage_name });
  }

  //Returns the newick representation of the tree for given distance matrix
  export async function hierarchicalClusteringByDistanceScript(data: DG.DataFrame, size: number, linkage_name: string): Promise<string> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringByDistanceScript', { data, size, linkage_name });
  }
}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Dendrogram:Info', {});
  }

  //Dendrogram tree visualization
  export async function dendrogram(): Promise<any> {
    return await grok.functions.call('Dendrogram:Dendrogram', {});
  }

  export async function getTreeHelper(): Promise<any> {
    return await grok.functions.call('Dendrogram:GetTreeHelper', {});
  }

  export async function getDendrogramService(): Promise<any> {
    return await grok.functions.call('Dendrogram:GetDendrogramService', {});
  }

  export async function generateTreeDialog(): Promise<any> {
    return await grok.functions.call('Dendrogram:GenerateTreeDialog', {});
  }

  //Test/demo app for Dendrogram
  export async function dendrogramApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:DendrogramApp', {});
  }

  //Test/demo app for Dendrogram Large
  export async function dendrogramLargeApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:DendrogramLargeApp', {});
  }

  //Test/demo app for TreeForGrid (custom renderer)
  export async function treeForGridApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridApp', {});
  }

  //Test/demo app for TreeForGridFilter (custom renderer)
  export async function treeForGridFilterApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridFilterApp', {});
  }

  //Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)
  export async function treeForGridCutApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridCutApp', {});
  }

  //Test/demo app for hierarchical clustering (inject tree to grid)
  export async function hierarchicalClusteringApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringApp', {});
  }

  //Test/demo app for hierarchical clustering (inject tree to grid)
  export async function hierarchicalClusteringSequencesApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringSequencesApp', {});
  }

  //Calculates hierarchical clustering on features and injects tree to grid
  export async function hierarchicalClustering(df: DG.DataFrame, colNameList: string[], distance: string, linkage: string): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClustering', { df, colNameList, distance, linkage });
  }

  //Opens Newick file
  export async function importNewick(fileContent: string): Promise<any> {
    return await grok.functions.call('Dendrogram:ImportNewick', { fileContent });
  }

  export async function previewNewick(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Dendrogram:PreviewNewick', { file });
  }

  //Calculates hierarchical clustering on features and injects tree to grid
  export async function hierarchicalClusteringSequences(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringSequences', {});
  }

  //Calculates hierarchical clustering on features and injects tree to grid
  export async function hierarchicalClusteringMolecules(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringMolecules', {});
  }

  //Calculates hierarchical clustering on features and injects tree to grid
  export async function hierarchicalClustering2(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClustering2', {});
  }

  //Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
  export async function heatMapDemo(): Promise<any> {
    return await grok.functions.call('Dendrogram:HeatMapDemo', {});
  }
}
