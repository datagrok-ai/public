/* eslint-disable max-len */
// import {SparseMatrixResult} from '../distance-matrix/sparse-matrix-service';
// import {bioLayout} from './bio-layout';
// import {MCLOptions, SparseMatrixObject} from './types';
// import {MCLOpReturnType} from '@datagrok-libraries/math/src/webGPU/MCL/types';
// import {markovClusterWebGPU} from '@datagrok-libraries/math/src/webGPU/MCL/MCL-webGPU';
// export const defaultMCLOptions: MCLOptions = {
//   expandFactor: 2,
//   maxIterations: 5,
//   inflateFactor: 2,
//   multFactor: 1,
//   pruneValue: 1e-15,
// };

// // P.S. This code is legacy!!!!!!!!!!!!! It is not used in the current version of the library.
// // It is kept here for reference purposes only.

// export class MCLSparseReducer {
//   private _options: MCLOptions;

//   constructor(opts: Partial<MCLOptions> = {}) {
//     this._options = {...defaultMCLOptions, ...opts};
//   }

//   public async transform(sparseMatrix: SparseMatrixResult, nRows: number) {
//     // testWorkerMultiply();
//     // return new Int32Array(nRows);
//     const originalSparseObject = this.toObjectForm(sparseMatrix);
//     const sparseClusters = this.assignClusters(originalSparseObject, nRows);
//     this.correctClusters(sparseClusters.clusters);
//     let sparseObject = originalSparseObject;
//     if (this._options.maxIterations > 0) {
//       this.addLoops(sparseObject, nRows);
//       this.normalize(sparseObject);
//       for (let i = 0; i < this._options.maxIterations; i++) {
//         sparseObject = this.expand(sparseObject, nRows);
//         this.normalize(sparseObject);
//         this.inflate(sparseObject);
//         this.normalize(sparseObject);
//       }
//     }
//     const {clusters} = this.assignClusters(sparseObject, nRows);
//     this.correctClusters(clusters);
//     //const embeddings = this.layoutInterConnect(clusters, nRows, sparseObject, sparseMatrix);
//     const embeddings = this.layout(sparseClusters.clusters, originalSparseObject, nRows, clusters);
//     return {clusters, embedX: embeddings.embedX, embedY: embeddings.embedY, is: sparseMatrix.i, js: sparseMatrix.j};
//   }

//   public async transformWebGPU(sparseMatrix: SparseMatrixResult, nRows: number) {
//     if (this._options.maxIterations === 0)
//       return this.transform(sparseMatrix, nRows);

//     const mclRes = await markovClusterWebGPU(
//       sparseMatrix, nRows, this._options.maxIterations, this._options.inflateFactor
//     );
//     const originalSparseObject = this.toObjectForm(sparseMatrix);
//     const sparseClusters = this.assignClusters(originalSparseObject, nRows);
//     this.correctClusters(sparseClusters.clusters);
//     const sparseObject = this.csrToSparseObject(mclRes, nRows);
//     const {clusters} = this.assignClusters(sparseObject, nRows);
//     this.correctClusters(clusters);
//     const embeddings = this.layout(sparseClusters.clusters, originalSparseObject, nRows, clusters);
//     return {clusters, embedX: embeddings.embedX, embedY: embeddings.embedY, is: sparseMatrix.i, js: sparseMatrix.j};
//   }

//   private correctClusters(clusters: number[]) {
//     const clusterSizeMap: {[_: number]: number} = {};
//     for (const cluster of clusters) {
//       if (!clusterSizeMap[cluster])
//         clusterSizeMap[cluster] = 0;
//       clusterSizeMap[cluster]++;
//     }
//     const sortedIndexes = Object.keys(clusterSizeMap).map(Number).sort((a, b) => clusterSizeMap[b] - clusterSizeMap[a]);
//     const clusterMap: {[_: number]: number} = {};
//     sortedIndexes.forEach((clusterIdx, i) => clusterMap[clusterIdx] = i + 1);
//     for (let i = 0; i < clusters.length; i++)
//       clusters[i] = clusterMap[clusters[i]];
//   }

//   private csrToSparseObject(mclRes: MCLOpReturnType, nRows: number) {
//     const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 1;
//     const minOrder = 1 / Math.pow(10, order);
//     const sparseObject: SparseMatrixObject = {};
//     for (let i = 0; i < nRows; i++) {
//       sparseObject[i] = {};
//       for (let k = mclRes.indexOffsets[i]; k < mclRes.indexOffsets[i + 1]; k++) {
//         const j = mclRes.KNNIndexes[k];
//         if (j <= i || mclRes.KNNSimilarities[k] < minOrder)
//           continue;
//         sparseObject[i][j] = mclRes.KNNSimilarities[k];
//       }
//     }
//     return sparseObject;
//   }

//   // private layoutInterConnect(
//   //   clusters: number[], nRows: number, sparseMatrix: SparseMatrixObject, originalSparse: SparseMatrixResult
//   // ) {
//   //   const embedX = new Float32Array(nRows).fill(0);
//   //   const embedY = new Float32Array(nRows).fill(0);
//   //   const clusterMap: {[_: number]: number[]} = {};
//   //   clusters.forEach((cluster, i) => {
//   //     if (!clusterMap[cluster])
//   //       clusterMap[cluster] = [];
//   //     clusterMap[cluster].push(i);
//   //   });
//   //   // first we layout the individual clusters
//   //   const sortedClusterNames = Object.keys(clusterMap);
//   //   sortedClusterNames.sort((a, b) => clusterMap[b as any].length - clusterMap[a as any].length);
//   //   const clusterEmbeddings: {[clusterNum: number]: {embedX: Float32Array, embedY: Float32Array}} = {};
//   //   for (const clusterName of sortedClusterNames) {
//   //     const cluster = clusterMap[clusterName as any]!;
//   //     // these embeddings will be in range 0-1
//   //     clusterEmbeddings[clusterName as any] = bioLayout(cluster, sparseMatrix, 0.001);
//   //     // lets say layout size is 10x10
//   //     // normilize the size of embeddings according to the size of the cluster
//   //     const clusterSize = Math.sqrt(cluster.length);
//   //     const embedX = clusterEmbeddings[clusterName as any].embedX;
//   //     const embedY = clusterEmbeddings[clusterName as any].embedY;
//   //     let minX = embedX[0];
//   //     let maxX = embedX[0];
//   //     let minY = embedY[0];
//   //     let maxY = embedY[0];
//   //     for (let i = 1; i < cluster.length; i++) {
//   //       minX = Math.min(minX, embedX[i]);
//   //       maxX = Math.max(maxX, embedX[i]);
//   //       minY = Math.min(minY, embedY[i]);
//   //       maxY = Math.max(maxY, embedY[i]);
//   //     }
//   //     const rangeX = maxX == minX ? 1 : maxX - minX;
//   //     const rangeY = maxY == minY ? 1: maxY - minY;
//   //     for (let i = 0; i < cluster.length; i++) {
//   //       embedX[i] = (embedX[i] - minX) / rangeX / 1.5 + 0.16;
//   //       embedY[i] = (embedY[i] - minY) / rangeY / 1.5 + 0.16;
//   //     }
//   //   }

//   //   // numbers of interCluster connections
//   //   const residualClusterConnections: {[clusterNum: number]: {[cluster2Num: number]: number}} = {};
//   //   const superClusters: {[_: number]: number} = {}; // maps original cluster number to super cluster number
//   //   let superClusterNum = 1;
//   //   // initialize superclusters to be the same as clusters
//   //   for (const clusterNum of Object.keys(clusterMap))
//   //     superClusters[Number(clusterNum)] = Number(clusterNum);
//   //   for (let con = 0; con < originalSparse.i.length; con ++) {
//   //     const i = originalSparse.i[con];
//   //     const j = originalSparse.j[con];
//   //     if (clusters[i] === clusters[j])
//   //       continue;
//   //     if (!residualClusterConnections[clusters[i]])
//   //       residualClusterConnections[clusters[i]] = {};
//   //     if (!residualClusterConnections[clusters[j]])
//   //       residualClusterConnections[clusters[j]] = {};
//   //     if (!residualClusterConnections[clusters[i]][clusters[j]])
//   //       residualClusterConnections[clusters[i]][clusters[j]] = 0;
//   //     if (!residualClusterConnections[clusters[j]][clusters[i]])
//   //       residualClusterConnections[clusters[j]][clusters[i]] = 0;
//   //     residualClusterConnections[clusters[i]][clusters[j]]++;
//   //     residualClusterConnections[clusters[j]][clusters[i]]++;
//   //     if (!superClusters[clusters[i]] && !superClusters[clusters[j]]) {
//   //       superClusters[clusters[i]] = superClusterNum;
//   //       superClusters[clusters[j]] = superClusterNum;
//   //       superClusterNum++;
//   //     } else if (superClusters[clusters[i]] && superClusters[clusters[j]]) {
//   //       const compVal = Number(superClusters[clusters[j]]);
//   //       for (const clusterNum of (Object.keys(superClusters) as unknown as number[])) {
//   //         if (Number(superClusters[clusterNum]) === compVal)
//   //           superClusters[clusterNum] = superClusters[clusters[i]];
//   //       }
//   //       superClusters[clusters[j]] = superClusters[clusters[i]];
//   //     } else {
//   //       const superCluster = superClusters[clusters[i]] ?? superClusters[clusters[j]];
//   //       superClusters[clusters[i]] = superCluster;
//   //       superClusters[clusters[j]] = superCluster;
//   //     }
//   //   }

//   //   // now we can get the reverse mapping of super clusters to clusters
//   //   const superClusterMap: {[_: number]: number[]} = {};
//   //   for (const clusterNum of Object.keys(superClusters)) {
//   //     const superClusterNum = superClusters[Number(clusterNum)];
//   //     if (!superClusterMap[superClusterNum])
//   //       superClusterMap[superClusterNum] = [];
//   //     superClusterMap[superClusterNum].push(Number(clusterNum));
//   //   }

//   //   // we need to have superClusters sorted by total size of points in them
//   //   const superClusterNums = Object.keys(superClusterMap).map(Number);
//   //   const superClusterSizes: {[supClustNum: number]: number} = {};
//   //   superClusterNums.forEach((superClusterNum) => {
//   //     let size = 0;
//   //     for (const clusterNum of superClusterMap[superClusterNum])
//   //       size += clusterMap[clusterNum].length;
//   //     superClusterSizes[superClusterNum] = size;
//   //   });
//   //   const sortedSuperClusterNums = superClusterNums.sort((a, b) => superClusterSizes[b] - superClusterSizes[a]);

//   //   // need to also sort clusters in each superCluster by size
//   //   for (const superClusterNum of sortedSuperClusterNums) {
//   //     const clusterNums = superClusterMap[superClusterNum];
//   //     clusterNums.sort((a, b) => clusterMap[b].length - clusterMap[a].length);
//   //   }

//   //   // now we need to layout superClusters


//   //   superClusterNum = 0;
//   //   let supPerRow = 3;
//   //   let supYOffset = 0;
//   //   const supLayoutSize = 15;
//   //   for (const superClusterName of sortedSuperClusterNums) {
//   //     // cs contains ids of clusters which are in this superCluster (superClusterName)
//   //     const cs = superClusterMap[superClusterName];
//   //     // space for 3 clusters in row(start with that), that actually fits 2 and rest is space between and arround
//   //     let perRow = 3;

//   //     let yOffset = 0;
//   //     let clusterNum = 0;
//   //     if (superClusterNum === Math.ceil(supPerRow / 1.5)) {
//   //       superClusterNum = 0;
//   //       supYOffset += supLayoutSize / supPerRow;
//   //       supPerRow = Math.ceil(supPerRow * 1.5);
//   //     }
//   //     const layoutSize = supLayoutSize / supPerRow;
//   //     const supOffsetX = ((superClusterNum % supPerRow + 0.25) * supLayoutSize / supPerRow) * 1.5;
//   //     let minX = Number.MAX_VALUE;
//   //     let maxX = Number.MIN_VALUE;
//   //     let minY = Number.MAX_VALUE;
//   //     let maxY = Number.MIN_VALUE;
//   //     for (const cn of cs) {
//   //       // cn is the id of subclaster
//   //       const clusterIs = clusterMap[cn];
//   //       // clusterIs contains indexes of points in this cluster
//   //       const embeddings = clusterEmbeddings[cn];
//   //       if (clusterNum === Math.ceil(perRow / 1.5)) {
//   //         clusterNum = 0;
//   //         yOffset += layoutSize / perRow;
//   //         perRow = Math.ceil(perRow * 1.5);
//   //       }
//   //       const offsetX = ((clusterNum % perRow) * layoutSize / perRow) * 1.5;
//   //       for (let i = 0; i < embeddings.embedX.length; i++) {
//   //         embedX[clusterIs[i]] = embeddings.embedX[i] * layoutSize / perRow + offsetX + supOffsetX;
//   //         embedY[clusterIs[i]] = embeddings.embedY[i] * layoutSize / perRow + yOffset + supYOffset;
//   //         minX = Math.min(minX, embedX[clusterIs[i]]);
//   //         maxX = Math.max(maxX, embedX[clusterIs[i]]);
//   //         minY = Math.min(minY, embedY[clusterIs[i]]);
//   //         maxY = Math.max(maxY, embedY[clusterIs[i]]);
//   //       }
//   //       clusterNum++;
//   //     }
//   //     // rescale all superclusters to their bounding box
//   //     const rangeX = maxX - minX == 0 ? 1 : maxX - minX;
//   //     const rangeY = maxY - minY == 0 ? 1 : maxY - minY;
//   //     const boundXMin = supOffsetX;
//   //     const boundYMin = supYOffset;
//   //     const boundXMax = supOffsetX + layoutSize;
//   //     const boundYMax = supYOffset + layoutSize;
//   //     for (const cn of cs) {
//   //       const clusterIs = clusterMap[cn];
//   //       for (let i = 0; i < clusterIs.length; i++) {
//   //         embedX[clusterIs[i]] = (embedX[clusterIs[i]] - minX) / rangeX / 1.2 * (boundXMax - boundXMin) + boundXMin;
//   //         embedY[clusterIs[i]] = (embedY[clusterIs[i]] - minY) / rangeY / 1.2 * (boundYMax - boundYMin) + boundYMin;
//   //       }
//   //     }

//   //     superClusterNum++;
//   //   }

//   //   return {embedX, embedY};
//   // }

//   private layout(clusters: number[], sparseMatrix: SparseMatrixObject, nRows: number, subCluster: number[]) {
//     const embedX = new Float32Array(nRows).fill(0);
//     const embedY = new Float32Array(nRows).fill(0);
//     const clusterMap: {[_: number]: number[]} = {};
//     clusters.forEach((cluster, i) => {
//       if (!clusterMap[cluster])
//         clusterMap[cluster] = [];
//       clusterMap[cluster].push(i);
//     });
//     // const nClusters = Object.keys(clusterMap).length;
//     // const perRow = Math.floor(Math.sqrt(nClusters));
//     let clusterNum = 0;
//     const sortedClusterNames = Object.keys(clusterMap);
//     sortedClusterNames.sort((a, b) => clusterMap[b as any].length - clusterMap[a as any].length);
//     let perRow = 1;

//     let yOffset = 0;
//     const layoutSize = 5;
//     for (const clusterName of sortedClusterNames) {
//       const cluster = clusterMap[clusterName as any]!;
//       const embeddings = bioLayout(cluster, sparseMatrix, 0.001, subCluster);
//       if (clusterNum === perRow) {
//         clusterNum = 0;
//         yOffset += layoutSize / perRow;
//         perRow = Math.min(Math.ceil(perRow * 3), 45);
//       }
//       //const clustersPerRow = Math.ceil(perRow / 1.5);
//       const offsetX = ((clusterNum) * layoutSize / perRow + layoutSize / perRow * (1 / 1.2 / 4));
//       // const offsetY = Math.floor(clusterNum / perRow) * 2;

//       for (let i = 0; i < embeddings.embedX.length; i++) {
//         embedX[cluster[i]] = embeddings.embedX[i] * layoutSize / perRow / 1.2 + offsetX;
//         embedY[cluster[i]] = embeddings.embedY[i] * layoutSize / perRow / 1.2 + yOffset;
//       }
//       clusterNum++;
//     }
//     return {embedX, embedY};
//   }

//   private mergeClusters(clusters: number[], i: number, j: number) {
//     const iCluster = clusters[i];
//     const jCluster = clusters[j];
//     for (let k = 0; k < clusters.length; k++) {
//       if (clusters[k] === jCluster)
//         clusters[k] = iCluster;
//     }
//   }
//   public assignClusters(sparseMatrix: SparseMatrixObject, nRows: number) {
//     let clusterNum = 0;
//     const is: number[] = [];
//     const js: number[] = [];
//     //const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 7;
//     const minOrder = Math.pow(10, -14);
//     const clusters: number[] = new Array(nRows).fill(-1);
//     for (const i of Object.keys(sparseMatrix)) {
//       for (const j of Object.keys(sparseMatrix[i as any])) {
//         if (sparseMatrix[i as any][j as any] > minOrder &&
//           sparseMatrix[i as any][j as any] !== Number(i) && Number(j) > Number(i)) {
//           is.push(Number(i));
//           js.push(Number(j));
//           if (clusters[Number(i)] !== -1 && clusters[Number(j)] !== -1) {
//             if (clusters[Number(i)] !== clusters[Number(j)])
//               this.mergeClusters(clusters, Number(i), Number(j));
//           } else if (clusters[Number(i)] !== -1) {
//             clusters[Number(j)] = clusters[Number(i)];
//           } else if (clusters[Number(j)] !== -1) {
//             clusters[Number(i)] = clusters[Number(j)];
//           } else {
//             clusterNum++;
//             clusters[Number(i)] = clusterNum;
//             clusters[Number(j)] = clusterNum;
//           }
//         }
//       }
//     }
//     for (let i=0; i < clusters.length; i++) {
//       if (clusters[i] === -1) {
//         clusterNum ++;
//         clusters[i] = clusterNum;
//       }
//     }
//     return {clusters, is: new Uint32Array(is), js: new Uint32Array(js)};
//   }

//   /** same as assign clusters but working on the level of CSR format returned from webGPU */
//   private assignClustersCSR(mclRes: MCLOpReturnType, nRows: number) {
//     let clusterNum = 0;
//     const is: number[] = [];
//     const js: number[] = [];
//     const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 1;
//     const minOrder = 1 / Math.pow(10, order);
//     const clusters: number[] = new Array(nRows).fill(-1);
//     const correctedOffsets = new Uint32Array(nRows + 1);
//     let offsetCounter = 0;
//     correctedOffsets[0] = 0;
//     for (let i = 0; i < nRows; i++) {
//       for (let k = mclRes.indexOffsets[i]; k < mclRes.indexOffsets[i + 1]; k++) {
//         const j = mclRes.KNNIndexes[k];
//         if (j <= i || mclRes.KNNSimilarities[k] <= minOrder)
//           continue;
//         is.push(i);
//         js.push(j);
//         offsetCounter++;
//         if (clusters[i] !== -1 && clusters[j] !== -1) {
//           if (clusters[i] !== clusters[j])
//             this.mergeClusters(clusters, i, j);
//         } else if (clusters[i] !== -1) {
//           clusters[j] = clusters[i];
//         } else if (clusters[j] !== -1) {
//           clusters[i] = clusters[j];
//         } else {
//           clusterNum++;
//           clusters[i] = clusterNum;
//           clusters[j] = clusterNum;
//         }
//       }
//       correctedOffsets[i + 1] = offsetCounter;
//     }
//     for (let i=0; i < clusters.length; i++) {
//       if (clusters[i] === -1) {
//         clusterNum ++;
//         clusters[i] = clusterNum;
//       }
//     }
//     return {clusters, is: new Uint32Array(is), js: new Uint32Array(js), correctedOffsets};
//   }

//   public toObjectForm(sparseMatrix: SparseMatrixResult): SparseMatrixObject {
//     const sparseObject: {[_: number]: {[_: number]: number}} = {};
//     for (let i = 0; i < sparseMatrix.i.length; i++) {
//       if (!sparseObject[sparseMatrix.i[i]])
//         sparseObject[sparseMatrix.i[i]] = {};
//       sparseObject[sparseMatrix.i[i]][sparseMatrix.j[i]] = 1 - sparseMatrix.distance[i];
//       if (!sparseObject[sparseMatrix.j[i]])
//         sparseObject[sparseMatrix.j[i]] = {};
//       sparseObject[sparseMatrix.j[i]][sparseMatrix.i[i]] = 1 - sparseMatrix.distance[i];
//     }
//     return sparseObject;
//   }

//   private addLoops(sparseObject: SparseMatrixObject, nRows: number) {
//     for (let i = 0; i < nRows; i++) {
//       if (!sparseObject[i])
//         sparseObject[i] = {};
//       sparseObject[i][i] = this._options.multFactor;
//     }
//   }

//   private normalize(sparseObject: SparseMatrixObject) {
//     for (const i of Object.keys(sparseObject)) {
//       const row = sparseObject[i as any];
//       let sum = 0;
//       for (const j of Object.keys(row))
//         sum += row[j as any];
//       if (sum === 0) continue;
//       for (const j of Object.keys(row))
//         sparseObject[i as any][j as any] /= sum;
//     }
//   }

//   private expand(sparseObject: SparseMatrixObject, nRows: number) {
//     const expandedObject: SparseMatrixObject = {};
//     const minOrder = 1000000000;
//     if (nRows < 12000) {
//       // if rows are less than 12000, translate to matrix form and use the dense matrix multiplication
//       // this is faster than the sparse matrix multiplication
//       const denseMatrix = new Array(nRows).fill(0).map(() => new Float32Array(nRows).fill(0));
//       for (const i of Object.keys(sparseObject)) {
//         const row = sparseObject[i as any];
//         for (const j of Object.keys(row))
//           denseMatrix[i as any][j as any] = row[j as any];
//       }
//       // expansion step
//       for (let i = 0; i < nRows; i++) {
//         for (let j = i + 1; j < nRows; j++) {
//           let val = 0;
//           for (let k = 0; k < nRows; k++)
//             val += denseMatrix[i][k] * denseMatrix[j][k];
//           if (Math.round(val * minOrder) > 0) {
//             expandedObject[i] ??= {};
//             expandedObject[i][j] = val;
//             expandedObject[j] ??= {};
//             expandedObject[j][i] = val;
//           }
//         }
//       }

//       return sparseObject;
//     }
//     for (let i = 0; i < nRows; i++) {
//       if (!sparseObject[i])
//         continue;
//       // const row = sparseObject[i];
//       expandedObject[i] ??= {};
//       for (let j = i; j < nRows; j++) {
//         if (!sparseObject[i]?.[j])
//           continue;
//         const val = this.getExpandValue(sparseObject, i, j); //pruning step
//         if (Math.round(val * minOrder) > 0) {
//           expandedObject[i][j] = val;
//           if (!expandedObject[j])
//             expandedObject[j] = {};
//           expandedObject[j][i] = val;
//         }
//       }
//     }
//     return expandedObject;
//   }

//   // private prune(row: SparseMatrixObject[number]) {

//   // }

//   private inflate(sparseObject: SparseMatrixObject) {
//     for (const i of Object.keys(sparseObject)) {
//       const row = sparseObject[i as any];
//       for (const j of Object.keys(row))
//         sparseObject[i as any][j as any] = Math.pow(sparseObject[i as any][j as any], this._options.inflateFactor);
//     }
//   }

//   private getExpandValue(sparseObject: SparseMatrixObject, i: any, j: any) {
//     let val = 0;
//     const currentIndexes = Object.keys(sparseObject[i] ?? {});
//     const otherIndexes = Object.keys(sparseObject[j] ?? {});
//     for (const k of currentIndexes) {
//       if (otherIndexes.includes(k))
//         val += sparseObject[i][k as any] * sparseObject[j][k as any];
//     }
//     return val;
//   }
// }

// // legacy
// function _getLayoutEmbeddings(is: number[], js: number[], cluster: number[], subCluster: number[]) {
//   // count number of unique clusters
//   const clusterSet = new Map<number, number>();
//   let counter = 0;
//   for (const c of cluster) {
//     if (!clusterSet.has(subCluster[c]))
//       clusterSet.set(subCluster[c], counter++);
//   }
//   const localConnectivity: {[key: number]: number} = {};
//   const interConnectivity: {[key: number]: number} = {};
//   // how many connections between clusters
//   const interClusterConnectivity: {[key: number]: {[key2: number]: number}} = {};
//   for (let i = 0; i < is.length; i++) {
//     if (subCluster[is[i]] === subCluster[js[i]]) {
//       localConnectivity[is[i]] = (localConnectivity[is[i]] ?? 0) + 1;
//       localConnectivity[js[i]] = (localConnectivity[js[i]] ?? 0) + 1;
//     } else {
//       interConnectivity[subCluster[is[i]]] = (interConnectivity[subCluster[is[i]]] ?? 0) + 1;
//       interConnectivity[subCluster[js[i]]] = (interConnectivity[subCluster[js[i]]] ?? 0) + 1;
//       interClusterConnectivity[subCluster[is[i]]] = interClusterConnectivity[subCluster[is[i]]] ?? {};
//       interClusterConnectivity[subCluster[is[i]]][subCluster[js[i]]] = (interClusterConnectivity[subCluster[is[i]]][subCluster[js[i]]] ?? 0) + 1;
//       interClusterConnectivity[subCluster[js[i]]] = interClusterConnectivity[subCluster[js[i]]] ?? {};
//       interClusterConnectivity[subCluster[js[i]]][subCluster[is[i]]] = (interClusterConnectivity[subCluster[js[i]]][subCluster[is[i]]] ?? 0) + 1;
//     }
//   }


//   const numSubClusters = clusterSet.size;
//   // we need to start with embeddings such that different subclusters are separated
//   // we will start with grid that can accommodate all subclusters on its perimeter
//   // perimeter of the grid is 2 * numCols + 2 * numRows - 4
//   // lets say we have square grid, then the perimeter is 4 * numCols - 4
//   // so numCols = (perimeter + 4) / 4

//   // const numCols = Math.ceil((numSubClusters + 4) / 4);
//   // const numRows = numCols;
//   const embedSize = 2;
//   // const emptySpaceMultiplier = 2;
//   const embedX = new Float32Array(cluster.length).fill(0).map((_, _i) => {
//     // const clusterNum = clusterSet.get(subCluster[cluster[i]])!;
//     // if (clusterNum < numCols)
//     //   return (clusterNum % numCols) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
//     // if (clusterNum < numCols + (numRows - 2) * 2)
//     //   return ((clusterNum - numCols) % 2) * (numCols - 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
//     // return ((clusterNum - numCols + (numRows - 2) * 2) % numCols) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
//     return Math.random() * embedSize;
//   });
//   const embedY = new Float32Array(cluster.length).fill(0).map((_, _i) => {
//     // const clusterNum = clusterSet.get(subCluster[cluster[i]])!;
//     // if (clusterNum < numCols)
//     //   return Math.random() * embedSize;
//     // if (clusterNum < numCols + (numRows - 2) * 2)
//     //   return Math.floor((clusterNum - numCols) / 2 + 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
//     // return (numRows - 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
//     return Math.random() * embedSize;
//   });
//   const iterations = 100;
//   const xVelocities = new Float32Array(cluster.length).fill(0);
//   const yVelocities = new Float32Array(cluster.length).fill(0);
//   // loop 2 for subClusters
//   for (let i = 0; i < iterations; i++) {
//     const temperature = (1 - i / iterations);
//     //xVelocities.fill(0);
//     //yVelocities.fill(0);
//     // loop one for inner cluster layout
//     for (let idx = 0; idx < is.length; idx++) {
//       const a = is[idx];
//       const b = js[idx];
//       if (subCluster[a] !== subCluster[b])
//         continue;
//       const diffX = embedX[a] - embedX[b];
//       const diffY = embedY[a] - embedY[b];
//       const dxa = diffX / localConnectivity[a];
//       const dya = diffY / localConnectivity[a];
//       const dxb = diffX / localConnectivity[b];
//       const dyb = diffY / localConnectivity[b];
//       const inverseForceXa = (0.2 - Math.abs(diffX)) / localConnectivity[a] * sign(diffX);
//       const inverseForceYa = (0.2 - Math.abs(diffY)) / localConnectivity[a] * sign(diffY);
//       const inverseForceXb = (0.2 - Math.abs(diffX)) / localConnectivity[b] * sign(diffX);
//       const inverseForceYb = (0.2 - Math.abs(diffY)) / localConnectivity[b] * sign(diffY);

//       if (Math.abs(diffX) > 0.01) {
//         xVelocities[a] -= dxa;
//         xVelocities[b] += dxb;
//       } else {
//         xVelocities[a] += inverseForceXa;
//         xVelocities[b] -= inverseForceXb;
//       }
//       if (Math.abs(diffY) > 0.01) {
//         yVelocities[a] -= dya;
//         yVelocities[b] += dyb;
//       } else {
//         yVelocities[a] += inverseForceYa;
//         yVelocities[b] -= inverseForceYb;
//       }
//       // yVelocities[a] -= dy;
//       // yVelocities[b] += dy;
//     }
//     for (let idx = 0; idx < cluster.length; idx++) {
//       xVelocities[idx] *= temperature;
//       yVelocities[idx] *= temperature;

//       embedX[idx] += xVelocities[idx];
//       embedY[idx] += yVelocities[idx];
//     }
//   }

//   // rescale each subClaster

//   // first fit everything in [0, 1] to just retain the scales
//   let minX = embedX[0];
//   let minY = embedY[0];
//   let maxX = embedX[0];
//   let maxY = embedY[0];
//   for (let idx = 1; idx < cluster.length; idx++) {
//     minX = Math.min(minX, embedX[idx]);
//     minY = Math.min(minY, embedY[idx]);
//     maxX = Math.max(maxX, embedX[idx]);
//     maxY = Math.max(maxY, embedY[idx]);
//   }
//   let rangeX = maxX - minX;
//   let rangeY = maxY - minY;
//   if (rangeX === 0)
//     rangeX = maxX === 0 ? 1 : maxX;
//   for (let idx = 0; idx < cluster.length; idx++)
//     embedX[idx] = (embedX[idx] - minX) / rangeX;
//   if (rangeY === 0)
//     rangeY = maxY === 0 ? 1 : maxY;
//   for (let idx = 0; idx < cluster.length; idx++)
//     embedY[idx] = (embedY[idx] - minY) / rangeY;


//   const minXs = new Float32Array(numSubClusters).fill(Number.MAX_VALUE);
//   const minYs = new Float32Array(numSubClusters).fill(Number.MAX_VALUE);
//   const maxXs = new Float32Array(numSubClusters).fill(Number.MIN_VALUE);
//   const maxYs = new Float32Array(numSubClusters).fill(Number.MIN_VALUE);
//   for (let idx = 0; idx < cluster.length; idx++) {
//     const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
//     minXs[clusterNum] = Math.min(minXs[clusterNum], embedX[idx]);
//     minYs[clusterNum] = Math.min(minYs[clusterNum], embedY[idx]);
//     maxXs[clusterNum] = Math.max(maxXs[clusterNum], embedX[idx]);
//     maxYs[clusterNum] = Math.max(maxYs[clusterNum], embedY[idx]);
//   }
//   const rangeXs = new Float32Array(numSubClusters).fill(0);
//   const rangeYs = new Float32Array(numSubClusters).fill(0);
//   for (let idx = 0; idx < numSubClusters; idx++) {
//     rangeXs[idx] = maxXs[idx] - minXs[idx];
//     rangeYs[idx] = maxYs[idx] - minYs[idx];
//     rangeXs[idx] = rangeXs[idx] === 0 ? 1 : rangeXs[idx];
//     rangeYs[idx] = rangeYs[idx] === 0 ? 1 : rangeYs[idx];
//   }

//   const centerXs = new Float32Array(numSubClusters).fill(0);
//   const centerYs = new Float32Array(numSubClusters).fill(0);
//   for (let idx = 0; idx < numSubClusters; idx++) {
//     centerXs[idx] = (minXs[idx] + maxXs[idx]) / 2;
//     centerYs[idx] = (minYs[idx] + maxYs[idx]) / 2;
//   }

//   // get cluster sizes
//   const clusterSizes = new Float32Array(numSubClusters).fill(0);
//   for (let idx = 0; idx < cluster.length; idx++) {
//     const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
//     clusterSizes[clusterNum]++;
//   }

//   // normalize cluster sizes and take square root as the size is in area
//   const maxClusterSizeSQ = Math.sqrt(clusterSizes.reduce((a, b) => Math.max(a, b), 0));
//   for (let idx = 0; idx < numSubClusters; idx++)
//     clusterSizes[idx] = Math.sqrt(clusterSizes[idx]) / maxClusterSizeSQ;

//   // position each centroid based on interconnectivity
//   xVelocities.fill(0);
//   yVelocities.fill(0);
//   for (let i = 0; i < iterations; i++) {
//     const temperature = (1 - i / iterations);
//     for (const clusterNum1 in interClusterConnectivity) {
//       for (const clusterNum2 in interClusterConnectivity[clusterNum1 as unknown as number]) {
//         const a = clusterSet.get(clusterNum1 as unknown as number)!;
//         const b = clusterSet.get(clusterNum2 as unknown as number)!;
//         const diffX = centerXs[a] - centerXs[b];
//         const diffY = centerYs[a] - centerYs[b];
//         const dxa = diffX * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum1];
//         const dya = diffY * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum1];
//         const dxb = diffX * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum2];
//         const dyb = diffY * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum2];

//         xVelocities[a] -= dxa;
//         xVelocities[b] += dxb;
//         yVelocities[a] -= dya;
//         yVelocities[b] += dyb;
//       }

//       for (let idx = 0; idx < numSubClusters; idx++) {
//         xVelocities[idx] *= temperature;
//         yVelocities[idx] *= temperature;

//         centerXs[idx] += xVelocities[idx];
//         centerYs[idx] += yVelocities[idx];
//       }
//     }
//   }

//   // normilize centroids
//   // minX = centerXs[0];
//   // minY = centerYs[0];
//   // maxX = centerXs[0];
//   // maxY = centerYs[0];

//   // for (let idx = 1; idx < numSubClusters; idx++) {
//   //   minX = Math.min(minX, centerXs[idx]);
//   //   minY = Math.min(minY, centerYs[idx]);
//   //   maxX = Math.max(maxX, centerXs[idx]);
//   //   maxY = Math.max(maxY, centerYs[idx]);
//   // }

//   // rangeX = maxX - minX;
//   // rangeY = maxY - minY;
//   // if (rangeX === 0)
//   //   rangeX = maxX === 0 ? 1 : maxX;
//   // for (let idx = 0; idx < numSubClusters; idx++)
//   //   centerXs[idx] = (centerXs[idx] - minX) / rangeX;
//   // if (rangeY === 0)
//   //   rangeY = maxY === 0 ? 1 : maxY;
//   // for (let idx = 0; idx < numSubClusters; idx++)
//   //   centerYs[idx] = (centerYs[idx] - minY) / rangeY;

//   // inflate centroids
//   inflateCentroids(centerXs, centerYs, clusterSizes);

//   // map the points to the center of the cluster
//   for (let idx = 0; idx < cluster.length; idx++) {
//     const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
//     const size = clusterSizes[clusterNum];
//     const centerX = centerXs[clusterNum];
//     const centerY = centerYs[clusterNum];
//     const clusterMidX = (minXs[clusterNum] + maxXs[clusterNum]) / 2;
//     const clusterMidY = (minYs[clusterNum] + maxYs[clusterNum]) / 2;
//     const dx = (embedX[idx] - clusterMidX) / rangeXs[clusterNum];
//     const dy = (embedY[idx] - clusterMidY) / rangeYs[clusterNum];
//     // make each cluster more circular
//     const angle = Math.atan2(dy, dx);
//     const dist = Math.sqrt(dx * dx + dy * dy);
//     const maxDist = size / 2;

//     embedX[idx] = centerX + Math.cos(angle) * Math.min(maxDist, dist);
//     embedY[idx] = centerY + Math.sin(angle) * Math.min(maxDist, dist);

//     //    embedX[idx] = centerX + dx * size;
//     //    embedY[idx] = centerY + dy * size;
//   }

//   // put all points in each cluster in center
//   // for (let idx = 0; idx < cluster.length; idx++) {
//   //   const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
//   //   embedX[idx] = centerXs[clusterNum];
//   //   embedY[idx] = centerYs[clusterNum];
//   // }


//   // normalize
//   minX = embedX[0];
//   minY = embedY[0];
//   maxX = embedX[0];
//   maxY = embedY[0];
//   for (let idx = 1; idx < cluster.length; idx++) {
//     minX = Math.min(minX, embedX[idx]);
//     minY = Math.min(minY, embedY[idx]);
//     maxX = Math.max(maxX, embedX[idx]);
//     maxY = Math.max(maxY, embedY[idx]);
//   }
//   rangeX = maxX - minX;
//   rangeY = maxY - minY;
//   if (rangeX === 0)
//     rangeX = maxX === 0 ? 1 : maxX;
//   for (let idx = 0; idx < cluster.length; idx++)
//     embedX[idx] = (embedX[idx] - minX) / rangeX;
//   if (rangeY === 0)
//     rangeY = maxY === 0 ? 1 : maxY;
//   for (let idx = 0; idx < cluster.length; idx++)
//     embedY[idx] = (embedY[idx] - minY) / rangeY;

//   return {embedX, embedY};
// }

// function sign(x: number) {
//   return x < 0 ? -1 : 1;
// }

// // this function will inflate the centroids of the clusters such that their size is proportional to the number of points in the cluster
// // sizes should be normalized to [0, 1]
// function inflateCentroids(centerXs: Float32Array, centerYs: Float32Array, clusterSizes: Float32Array) {
//   const sortedIndexes = Array.from(clusterSizes.keys()).sort((a, b) => clusterSizes[b] - clusterSizes[a]);
//   const numClusters = centerXs.length;
//   for (let uidx = 0; uidx < numClusters; uidx++) {
//     const idx = sortedIndexes[uidx];
//     const size = clusterSizes[idx];
//     const centerX = centerXs[idx];
//     const centerY = centerYs[idx];
//     for (let ui = 0; ui < numClusters; ui++) {
//       const i = sortedIndexes[ui];
//       const size2 = clusterSizes[i];
//       const minDist = (size + size2);
//       if (i === idx)
//         continue;
//       const dx = centerX - centerXs[i];
//       const dy = centerY - centerYs[i];
//       const dist = Math.sqrt(dx * dx + dy * dy);
//       if (dist < minDist) {
//         const diff = size / 2;
//         const angle = Math.atan2(dy, dx);
//         centerXs[i] -= diff * Math.cos(angle);
//         centerYs[i] -= diff * Math.sin(angle);
//       }
//     }
//   }
// }
