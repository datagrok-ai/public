import {SparseMatrix} from '@datagrok-libraries/math';
import {MCLOptions} from './types';
import {runMarkovClustering} from './clustering-steps';

import {getWebColaLayot} from './webCola';

export const defaultMCLOptions: MCLOptions = {
  expandFactor: 2,
  maxIterations: 5,
  inflateFactor: 2,
  multFactor: 1,
  pruneValue: 1e-17,
};

export class MCLSparseReducer {
    private _options: MCLOptions;

    constructor(opts: Partial<MCLOptions> = {}) {
      this._options = {...defaultMCLOptions, ...opts};
    }

    public async transform(sparseMatrix: SparseMatrix, nRows: number) {
      // get the superClusters before the MCL
      // this will be used for layouting
      const superClusters = this.assignClusters(sparseMatrix, nRows);
      this.correctClusters(superClusters.clusters);
      // before beggining the mcl, we need to save the connections between original points,
      // as MCL might mutate the matrix
      const clusterConnectionMap = this.splitConnectionsIntoClusters(sparseMatrix, superClusters.clusters);
      // perform the MCL
      const mclMatrix = await runMarkovClustering(
        sparseMatrix, nRows, this._options.inflateFactor, this._options.maxIterations, this._options.pruneValue);
      const {clusters} = this.assignClusters(mclMatrix, nRows);
      this.correctClusters(clusters);
      const embeddings = this.layout(superClusters.clusters, clusterConnectionMap, nRows, clusters);
      return {clusters, embedX: embeddings.embedX, embedY: embeddings.embedY, is: sparseMatrix.i, js: sparseMatrix.j};
    }

    public async transformWebGPU(sparseMatrix: SparseMatrix, nRows: number) {
      // TODO: implement correct webGPU version
      return this.transform(sparseMatrix, nRows);
    }

    // here as we operate on original sparse matrix, we know for sure that there are no duplicates or self loops
    private splitConnectionsIntoClusters(sparseMatrix: SparseMatrix, clusters: number[]) {
      const clusterConnections = new Map<number, {i: number[], j: number[], v: number[]}>();
      for (let i = 0; i < clusters.length; i++) {
        if (!clusterConnections.has(clusters[i]))
          clusterConnections.set(clusters[i], {i: [], j: [], v: []});
      }
      for (let i = 0; i < sparseMatrix.i.length; i++) {
        const cluster = clusters[sparseMatrix.i[i]];
        const cc = clusterConnections.get(cluster)!;
        cc.i.push(sparseMatrix.i[i]);
        cc.j.push(sparseMatrix.j[i]);
        cc.v.push(1 - sparseMatrix.distance[i]);
      }
      return clusterConnections;
    }

    private assignClusters(sparseMatrix: SparseMatrix, nRows: number) {
      let clusterNum = 0;
      const is: number[] = [];
      const js: number[] = [];
      // array containing for each point which cluster it belongs to
      const clusters: number[] = new Array(nRows).fill(-1);
      for (let it = 0; it < sparseMatrix.i.length; it++) {
        const i = sparseMatrix.i[it];
        const j = sparseMatrix.j[it];
        if (i === j)
          continue;
        is.push(i);
        js.push(j);
        if (clusters[i] !== -1 && clusters[j] !== -1) {
          if (clusters[i] !== clusters[j])
            this.mergeClusters(clusters, i, j);
        } else if (clusters[i] !== -1) {
          clusters[j] = clusters[i];
        } else if (clusters[j] !== -1) {
          clusters[i] = clusters[j];
        } else {
          clusterNum++;
          clusters[i] = clusterNum;
          clusters[j] = clusterNum;
        }
      }

      // final step for the clusters that are not connected to anything
      for (let i=0; i < clusters.length; i++) {
        if (clusters[i] === -1) {
          clusterNum ++;
          clusters[i] = clusterNum;
        }
      }
      return {clusters, is: new Uint32Array(is), js: new Uint32Array(js)};
    }

    private mergeClusters(clusters: number[], i: number, j: number) {
      const iCluster = clusters[i];
      const jCluster = clusters[j];
      for (let k = 0; k < clusters.length; k++) {
        if (clusters[k] === jCluster)
          clusters[k] = iCluster;
      }
    }

    /** After assigning and merging, clusters will need reordering according to size */
    private correctClusters(clusters: number[]) {
      const clusterSizeMap: {[_: number]: number} = {};
      for (const cluster of clusters) {
        if (!clusterSizeMap[cluster])
          clusterSizeMap[cluster] = 0;
        clusterSizeMap[cluster]++;
      }
      const sortedIndexes =
        Object.keys(clusterSizeMap).map(Number).sort((a, b) => clusterSizeMap[b] - clusterSizeMap[a]);
      const clusterMap: {[_: number]: number} = {};
      sortedIndexes.forEach((clusterIdx, i) => clusterMap[clusterIdx] = i + 1);
      for (let i = 0; i < clusters.length; i++)
        clusters[i] = clusterMap[clusters[i]];
    }


    /** notice that here, first argument is the superclusters and last is the subClusters
     * the second argument is the original sparse matrix, and the third is the number of rows
     */
    private layout(clusters: number[], clusterConnectionMap: Map<number, {
        i: number[];
        j: number[];
        v: number[];}>, nRows: number, subCluster: number[]) {
      const embedX = new Float32Array(nRows).fill(0);
      const embedY = new Float32Array(nRows).fill(0);
      const clusterMap: {[_: number]: number[]} = {};
      clusters.forEach((cluster, i) => {
        if (!clusterMap[cluster])
          clusterMap[cluster] = [];
        clusterMap[cluster].push(i);
      });

      // split sparse matrix connections into super-clusters, save only indexes
      let clusterNum = 0;
      const sortedClusterNames = Object.keys(clusterMap);
      sortedClusterNames.sort((a, b) => clusterMap[b as any].length - clusterMap[a as any].length);
      let perRow = 1;

      let yOffset = 0;
      const layoutSize = 5;
      for (const clusterName of sortedClusterNames) {
        const cluster = clusterMap[clusterName as any]!;
        const clusterConnections = clusterConnectionMap.get(Number(clusterName))!;
        const embeddings = getWebColaLayot(cluster, clusterConnections, subCluster);
        if (clusterNum === perRow) {
          clusterNum = 0;
          yOffset += layoutSize / perRow;
          perRow = Math.min(Math.ceil(perRow * 3), 45);
        }
        //const clustersPerRow = Math.ceil(perRow / 1.5);
        const offsetX = ((clusterNum) * layoutSize / perRow + layoutSize / perRow * (1 / 1.2 / 4));
        // const offsetY = Math.floor(clusterNum / perRow) * 2;

        for (let i = 0; i < embeddings.embedX.length; i++) {
          embedX[cluster[i]] = embeddings.embedX[i] * layoutSize / perRow / 1.2 + offsetX;
          embedY[cluster[i]] = embeddings.embedY[i] * layoutSize / perRow / 1.2 + yOffset;
        }
        clusterNum++;
      }
      return {embedX, embedY};
    }
}


