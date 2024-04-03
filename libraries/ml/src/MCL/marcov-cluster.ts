import {SparseMatrixResult} from '../distance-matrix/sparse-matrix-service';
import {bioLayout} from './bio-layout';
import {MCLOptions, SparseMatrixObject} from './types';

export const defaultMCLOptions: MCLOptions = {
  expandFactor: 2,
  maxIterations: 5,
  inflateFactor: 2,
  multFactor: 1,
};

export class MCLSparseReducer {
  private _options: MCLOptions;

  constructor(opts: Partial<MCLOptions> = {}) {
    this._options = {...defaultMCLOptions, ...opts};
  }

  public async transform(sparseMatrix: SparseMatrixResult, nRows: number) {
    // testWorkerMultiply();
    // return new Int32Array(nRows);
    let sparseObject = this.toObjectForm(sparseMatrix);
    if (this._options.maxIterations > 0) {
      this.addLoops(sparseObject, nRows);
      this.normalize(sparseObject);
      for (let i = 0; i < this._options.maxIterations; i++) {
        sparseObject = this.expand(sparseObject, nRows);
        this.inflate(sparseObject);
        this.normalize(sparseObject);
      }
    }
    const {clusters, is, js} = this.assignClusters(sparseObject, nRows);
    this.correctClusters(clusters);
    const embeddings = await this.layout(clusters, sparseObject, nRows);
    return {clusters, embedX: embeddings.embedX, embedY: embeddings.embedY, is, js};
  }

  private correctClusters(clusters: number[]) {
    const clusterSizeMap: {[_: number]: number} = {};
    for (const cluster of clusters) {
      if (!clusterSizeMap[cluster])
        clusterSizeMap[cluster] = 0;
      clusterSizeMap[cluster]++;
    }
    const sortedIndexes = Object.keys(clusterSizeMap).map(Number).sort((a, b) => clusterSizeMap[a] - clusterSizeMap[b]);
    const clusterMap: {[_: number]: number} = {};
    sortedIndexes.forEach((clusterIdx, i) => clusterMap[clusterIdx] = i + 1);
    for (let i = 0; i < clusters.length; i++)
      clusters[i] = clusterMap[clusters[i]];
    // let curCluster = 1;
    // for (let i = 0; i < clusters.length; i++) {
    //   if (!clusterMap[clusters[i]]) {
    //     clusterMap[clusters[i]] = curCluster;
    //     clusters[i] = curCluster;
    //     curCluster++;
    //   } else {
    //     clusters[i] = clusterMap[clusters[i]];
    //   }
    // }
  }

  private async layout(clusters: number[], sparseMatrix: SparseMatrixObject, nRows: number) {
    const embedX = new Float32Array(nRows).fill(0);
    const embedY = new Float32Array(nRows).fill(0);
    const clusterMap: {[_: number]: number[]} = {};
    clusters.forEach((cluster, i) => {
      if (!clusterMap[cluster])
        clusterMap[cluster] = [];
      clusterMap[cluster].push(i);
    });
    // const nClusters = Object.keys(clusterMap).length;
    // const perRow = Math.floor(Math.sqrt(nClusters));
    let clusterNum = 0;
    const sortedClusterNames = Object.keys(clusterMap);
    sortedClusterNames.sort((a, b) => clusterMap[b as any].length - clusterMap[a as any].length);
    let perRow = 6;

    let yOffset = 0;
    const layoutSize = 5;
    for (const clusterName of sortedClusterNames) {
      const cluster = clusterMap[clusterName as any]!;
      const embeddings = await bioLayout(cluster, sparseMatrix, 0.001);
      if (clusterNum === Math.ceil(perRow / 1.5)) {
        clusterNum = 0;
        yOffset += layoutSize / perRow;
        perRow = Math.ceil(perRow * 1.5);
      }
      const offsetX = ((clusterNum % perRow) * layoutSize / perRow) * 1.5;
      // const offsetY = Math.floor(clusterNum / perRow) * 2;

      for (let i = 0; i < embeddings.embedX.length; i++) {
        embedX[cluster[i]] = embeddings.embedX[i] * layoutSize / perRow + offsetX;
        embedY[cluster[i]] = embeddings.embedY[i] * layoutSize / perRow + yOffset;
      }
      clusterNum++;
    }
    return {embedX, embedY};
  }

  private mergeClusters(clusters: number[], i: number, j: number) {
    const iCluster = clusters[i];
    const jCluster = clusters[j];
    for (let k = 0; k < clusters.length; k++) {
      if (clusters[k] === jCluster)
        clusters[k] = iCluster;
    }
  }
  public assignClusters(sparseMatrix: SparseMatrixObject, nRows: number) {
    let clusterNum = 0;
    const is: number[] = [];
    const js: number[] = [];
    const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 1;
    const minOrder = Math.pow(10, order);
    const clusters: number[] = new Array(nRows).fill(-1);
    for (const i of Object.keys(sparseMatrix)) {
      for (const j of Object.keys(sparseMatrix[i as any])) {
        if (Math.round(sparseMatrix[i as any][j as any] * minOrder) / minOrder > 0 &&
          sparseMatrix[i as any][j as any] !== Number(i) && Number(j) > Number(i)) {
          is.push(Number(i));
          js.push(Number(j));
          if (clusters[Number(i)] !== -1 && clusters[Number(j)] !== -1) {
            if (clusters[Number(i)] !== clusters[Number(j)])
              this.mergeClusters(clusters, Number(i), Number(j));
          } else if (clusters[Number(i)] !== -1) {
            clusters[Number(j)] = clusters[Number(i)];
          } else if (clusters[Number(j)] !== -1) {
            clusters[Number(i)] = clusters[Number(j)];
          } else {
            clusterNum++;
            clusters[Number(i)] = clusterNum;
            clusters[Number(j)] = clusterNum;
          }
        }
      }
    }
    for (let i=0; i < clusters.length; i++) {
      if (clusters[i] === -1) {
        clusterNum ++;
        clusters[i] = clusterNum;
      }
    }
    return {clusters, is: new Uint32Array(is), js: new Uint32Array(js)};
  }

  public toObjectForm(sparseMatrix: SparseMatrixResult): SparseMatrixObject {
    const sparseObject: {[_: number]: {[_: number]: number}} = {};
    for (let i = 0; i < sparseMatrix.i.length; i++) {
      if (!sparseObject[sparseMatrix.i[i]])
        sparseObject[sparseMatrix.i[i]] = {};
      sparseObject[sparseMatrix.i[i]][sparseMatrix.j[i]] = 1 - sparseMatrix.distance[i];
      if (!sparseObject[sparseMatrix.j[i]])
        sparseObject[sparseMatrix.j[i]] = {};
      sparseObject[sparseMatrix.j[i]][sparseMatrix.i[i]] = 1 - sparseMatrix.distance[i];
    }
    return sparseObject;
  }

  private addLoops(sparseObject: SparseMatrixObject, nRows: number) {
    for (let i = 0; i < nRows; i++) {
      if (!sparseObject[i])
        sparseObject[i] = {};
      sparseObject[i][i] = this._options.multFactor;
    }
  }

  private normalize(sparseObject: SparseMatrixObject) {
    for (const i of Object.keys(sparseObject)) {
      const row = sparseObject[i as any];
      let sum = 0;
      for (const j of Object.keys(row))
        sum += row[j as any];
      if (sum === 0) continue;
      for (const j of Object.keys(row))
        sparseObject[i as any][j as any] /= sum;
    }
  }

  private expand(sparseObject: SparseMatrixObject, nRows: number) {
    const expandedObject: SparseMatrixObject = {};
    const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 1;
    const minOrder = Math.pow(10, order);
    for (let i = 0; i < nRows; i++) {
      if (!sparseObject[i])
        continue;
      // const row = sparseObject[i];
      expandedObject[i] = {};
      for (let j = i; j < nRows; j++) {
        if (!sparseObject[i]?.[j])
          continue;
        const val = this.getExpandValue(sparseObject, i, j); //pruning step
        if (Math.round(val * minOrder) / minOrder > 0) {
          expandedObject[i][j] = val;
          if (!expandedObject[j])
            expandedObject[j] = {};
          expandedObject[j][i] = val;
        }
      }
    }
    return expandedObject;
  }

  // private prune(row: SparseMatrixObject[number]) {

  // }

  private inflate(sparseObject: SparseMatrixObject) {
    for (const i of Object.keys(sparseObject)) {
      const row = sparseObject[i as any];
      for (const j of Object.keys(row))
        sparseObject[i as any][j as any] = Math.pow(sparseObject[i as any][j as any], this._options.inflateFactor);
    }
  }

  private getExpandValue(sparseObject: SparseMatrixObject, i: any, j: any) {
    let val = 0;
    const currentIndexes = Object.keys(sparseObject[i] ?? {});
    const otherIndexes = Object.keys(sparseObject[j] ?? {});
    for (const k of currentIndexes) {
      if (otherIndexes.includes(k))
        val += sparseObject[i][k as any] * sparseObject[j][k as any];
    }
    return val;
  }
}

