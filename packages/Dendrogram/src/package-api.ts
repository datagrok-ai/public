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

  export async function dendrogramApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:DendrogramApp', {});
  }

  export async function dendrogramLargeApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:DendrogramLargeApp', {});
  }

  export async function treeForGridApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridApp', {});
  }

  export async function treeForGridFilterApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridFilterApp', {});
  }

  export async function treeForGridCutApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:TreeForGridCutApp', {});
  }

  export async function hierarchicalClusteringApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringApp', {});
  }

  export async function hierarchicalClusteringSequencesApp(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringSequencesApp', {});
  }

  export async function hierarchicalClustering(df: DG.DataFrame, colNameList: string[], distance: any, linkage: string): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClustering', { df, colNameList, distance, linkage });
  }

  export async function importNewick(fileContent: string): Promise<any> {
    return await grok.functions.call('Dendrogram:ImportNewick', { fileContent });
  }

  export async function previewNewick(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Dendrogram:PreviewNewick', { file });
  }

  export async function hierarchicalClusteringSequences(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringSequences', {});
  }

  export async function hierarchicalClusteringMolecules(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClusteringMolecules', {});
  }

  export async function hierarchicalClustering2(): Promise<any> {
    return await grok.functions.call('Dendrogram:HierarchicalClustering2', {});
  }

  export async function heatMapDemo(): Promise<any> {
    return await grok.functions.call('Dendrogram:HeatMapDemo', {});
  }
}
