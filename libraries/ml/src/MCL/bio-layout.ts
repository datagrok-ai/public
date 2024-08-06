/* eslint-disable max-len */
import {SparseMatrixObject} from './types';
import {getWebColaLayot} from './webCola';

export function bioLayout(cluster: number[], sparseObject: SparseMatrixObject,
  _threshold: number, subCluster: number[]) {
  const n = cluster.length;
  const is: number[] = [];
  const js: number[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a = cluster[i];
      const b = cluster[j];
      if ((sparseObject[a]?.[b] ?? 0) > 0) {
        is.push(i);
        js.push(j);
      }
    }
  }

  return getWebColaLayot(is, js, cluster, subCluster);
}

// legacy
function _getLayoutEmbeddings(is: number[], js: number[], cluster: number[], subCluster: number[]) {
  // count number of unique clusters
  const clusterSet = new Map<number, number>();
  let counter = 0;
  for (const c of cluster) {
    if (!clusterSet.has(subCluster[c]))
      clusterSet.set(subCluster[c], counter++);
  }
  const localConnectivity: {[key: number]: number} = {};
  const interConnectivity: {[key: number]: number} = {};
  // how many connections between clusters
  const interClusterConnectivity: {[key: number]: {[key2: number]: number}} = {};
  for (let i = 0; i < is.length; i++) {
    if (subCluster[is[i]] === subCluster[js[i]]) {
      localConnectivity[is[i]] = (localConnectivity[is[i]] ?? 0) + 1;
      localConnectivity[js[i]] = (localConnectivity[js[i]] ?? 0) + 1;
    } else {
      interConnectivity[subCluster[is[i]]] = (interConnectivity[subCluster[is[i]]] ?? 0) + 1;
      interConnectivity[subCluster[js[i]]] = (interConnectivity[subCluster[js[i]]] ?? 0) + 1;
      interClusterConnectivity[subCluster[is[i]]] = interClusterConnectivity[subCluster[is[i]]] ?? {};
      interClusterConnectivity[subCluster[is[i]]][subCluster[js[i]]] = (interClusterConnectivity[subCluster[is[i]]][subCluster[js[i]]] ?? 0) + 1;
      interClusterConnectivity[subCluster[js[i]]] = interClusterConnectivity[subCluster[js[i]]] ?? {};
      interClusterConnectivity[subCluster[js[i]]][subCluster[is[i]]] = (interClusterConnectivity[subCluster[js[i]]][subCluster[is[i]]] ?? 0) + 1;
    }
  }


  const numSubClusters = clusterSet.size;
  // we need to start with embeddings such that different subclusters are separated
  // we will start with grid that can accommodate all subclusters on its perimeter
  // perimeter of the grid is 2 * numCols + 2 * numRows - 4
  // lets say we have square grid, then the perimeter is 4 * numCols - 4
  // so numCols = (perimeter + 4) / 4

  // const numCols = Math.ceil((numSubClusters + 4) / 4);
  // const numRows = numCols;
  const embedSize = 2;
  // const emptySpaceMultiplier = 2;
  const embedX = new Float32Array(cluster.length).fill(0).map((_, _i) => {
    // const clusterNum = clusterSet.get(subCluster[cluster[i]])!;
    // if (clusterNum < numCols)
    //   return (clusterNum % numCols) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
    // if (clusterNum < numCols + (numRows - 2) * 2)
    //   return ((clusterNum - numCols) % 2) * (numCols - 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
    // return ((clusterNum - numCols + (numRows - 2) * 2) % numCols) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
    return Math.random() * embedSize;
  });
  const embedY = new Float32Array(cluster.length).fill(0).map((_, _i) => {
    // const clusterNum = clusterSet.get(subCluster[cluster[i]])!;
    // if (clusterNum < numCols)
    //   return Math.random() * embedSize;
    // if (clusterNum < numCols + (numRows - 2) * 2)
    //   return Math.floor((clusterNum - numCols) / 2 + 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
    // return (numRows - 1) * embedSize * emptySpaceMultiplier + Math.random() * embedSize;
    return Math.random() * embedSize;
  });
  const iterations = 100;
  const xVelocities = new Float32Array(cluster.length).fill(0);
  const yVelocities = new Float32Array(cluster.length).fill(0);
  // loop 2 for subClusters
  for (let i = 0; i < iterations; i++) {
    const temperature = (1 - i / iterations);
    //xVelocities.fill(0);
    //yVelocities.fill(0);
    // loop one for inner cluster layout
    for (let idx = 0; idx < is.length; idx++) {
      const a = is[idx];
      const b = js[idx];
      if (subCluster[a] !== subCluster[b])
        continue;
      const diffX = embedX[a] - embedX[b];
      const diffY = embedY[a] - embedY[b];
      const dxa = diffX / localConnectivity[a];
      const dya = diffY / localConnectivity[a];
      const dxb = diffX / localConnectivity[b];
      const dyb = diffY / localConnectivity[b];
      const inverseForceXa = (0.2 - Math.abs(diffX)) / localConnectivity[a] * sign(diffX);
      const inverseForceYa = (0.2 - Math.abs(diffY)) / localConnectivity[a] * sign(diffY);
      const inverseForceXb = (0.2 - Math.abs(diffX)) / localConnectivity[b] * sign(diffX);
      const inverseForceYb = (0.2 - Math.abs(diffY)) / localConnectivity[b] * sign(diffY);

      if (Math.abs(diffX) > 0.01) {
        xVelocities[a] -= dxa;
        xVelocities[b] += dxb;
      } else {
        xVelocities[a] += inverseForceXa;
        xVelocities[b] -= inverseForceXb;
      }
      if (Math.abs(diffY) > 0.01) {
        yVelocities[a] -= dya;
        yVelocities[b] += dyb;
      } else {
        yVelocities[a] += inverseForceYa;
        yVelocities[b] -= inverseForceYb;
      }
      // yVelocities[a] -= dy;
      // yVelocities[b] += dy;
    }
    for (let idx = 0; idx < cluster.length; idx++) {
      xVelocities[idx] *= temperature;
      yVelocities[idx] *= temperature;

      embedX[idx] += xVelocities[idx];
      embedY[idx] += yVelocities[idx];
    }
  }

  // rescale each subClaster

  // first fit everything in [0, 1] to just retain the scales
  let minX = embedX[0];
  let minY = embedY[0];
  let maxX = embedX[0];
  let maxY = embedY[0];
  for (let idx = 1; idx < cluster.length; idx++) {
    minX = Math.min(minX, embedX[idx]);
    minY = Math.min(minY, embedY[idx]);
    maxX = Math.max(maxX, embedX[idx]);
    maxY = Math.max(maxY, embedY[idx]);
  }
  let rangeX = maxX - minX;
  let rangeY = maxY - minY;
  if (rangeX === 0)
    rangeX = maxX === 0 ? 1 : maxX;
  for (let idx = 0; idx < cluster.length; idx++)
    embedX[idx] = (embedX[idx] - minX) / rangeX;
  if (rangeY === 0)
    rangeY = maxY === 0 ? 1 : maxY;
  for (let idx = 0; idx < cluster.length; idx++)
    embedY[idx] = (embedY[idx] - minY) / rangeY;


  const minXs = new Float32Array(numSubClusters).fill(Number.MAX_VALUE);
  const minYs = new Float32Array(numSubClusters).fill(Number.MAX_VALUE);
  const maxXs = new Float32Array(numSubClusters).fill(Number.MIN_VALUE);
  const maxYs = new Float32Array(numSubClusters).fill(Number.MIN_VALUE);
  for (let idx = 0; idx < cluster.length; idx++) {
    const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
    minXs[clusterNum] = Math.min(minXs[clusterNum], embedX[idx]);
    minYs[clusterNum] = Math.min(minYs[clusterNum], embedY[idx]);
    maxXs[clusterNum] = Math.max(maxXs[clusterNum], embedX[idx]);
    maxYs[clusterNum] = Math.max(maxYs[clusterNum], embedY[idx]);
  }
  const rangeXs = new Float32Array(numSubClusters).fill(0);
  const rangeYs = new Float32Array(numSubClusters).fill(0);
  for (let idx = 0; idx < numSubClusters; idx++) {
    rangeXs[idx] = maxXs[idx] - minXs[idx];
    rangeYs[idx] = maxYs[idx] - minYs[idx];
    rangeXs[idx] = rangeXs[idx] === 0 ? 1 : rangeXs[idx];
    rangeYs[idx] = rangeYs[idx] === 0 ? 1 : rangeYs[idx];
  }

  const centerXs = new Float32Array(numSubClusters).fill(0);
  const centerYs = new Float32Array(numSubClusters).fill(0);
  for (let idx = 0; idx < numSubClusters; idx++) {
    centerXs[idx] = (minXs[idx] + maxXs[idx]) / 2;
    centerYs[idx] = (minYs[idx] + maxYs[idx]) / 2;
  }

  // get cluster sizes
  const clusterSizes = new Float32Array(numSubClusters).fill(0);
  for (let idx = 0; idx < cluster.length; idx++) {
    const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
    clusterSizes[clusterNum]++;
  }

  // normalize cluster sizes and take square root as the size is in area
  const maxClusterSizeSQ = Math.sqrt(clusterSizes.reduce((a, b) => Math.max(a, b), 0));
  for (let idx = 0; idx < numSubClusters; idx++)
    clusterSizes[idx] = Math.sqrt(clusterSizes[idx]) / maxClusterSizeSQ;

  // position each centroid based on interconnectivity
  xVelocities.fill(0);
  yVelocities.fill(0);
  for (let i = 0; i < iterations; i++) {
    const temperature = (1 - i / iterations);
    for (const clusterNum1 in interClusterConnectivity) {
      for (const clusterNum2 in interClusterConnectivity[clusterNum1 as unknown as number]) {
        const a = clusterSet.get(clusterNum1 as unknown as number)!;
        const b = clusterSet.get(clusterNum2 as unknown as number)!;
        const diffX = centerXs[a] - centerXs[b];
        const diffY = centerYs[a] - centerYs[b];
        const dxa = diffX * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum1];
        const dya = diffY * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum1];
        const dxb = diffX * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum2];
        const dyb = diffY * interClusterConnectivity[clusterNum1][clusterNum2] / interConnectivity[clusterNum2];

        xVelocities[a] -= dxa;
        xVelocities[b] += dxb;
        yVelocities[a] -= dya;
        yVelocities[b] += dyb;
      }

      for (let idx = 0; idx < numSubClusters; idx++) {
        xVelocities[idx] *= temperature;
        yVelocities[idx] *= temperature;

        centerXs[idx] += xVelocities[idx];
        centerYs[idx] += yVelocities[idx];
      }
    }
  }

  // normilize centroids
  // minX = centerXs[0];
  // minY = centerYs[0];
  // maxX = centerXs[0];
  // maxY = centerYs[0];

  // for (let idx = 1; idx < numSubClusters; idx++) {
  //   minX = Math.min(minX, centerXs[idx]);
  //   minY = Math.min(minY, centerYs[idx]);
  //   maxX = Math.max(maxX, centerXs[idx]);
  //   maxY = Math.max(maxY, centerYs[idx]);
  // }

  // rangeX = maxX - minX;
  // rangeY = maxY - minY;
  // if (rangeX === 0)
  //   rangeX = maxX === 0 ? 1 : maxX;
  // for (let idx = 0; idx < numSubClusters; idx++)
  //   centerXs[idx] = (centerXs[idx] - minX) / rangeX;
  // if (rangeY === 0)
  //   rangeY = maxY === 0 ? 1 : maxY;
  // for (let idx = 0; idx < numSubClusters; idx++)
  //   centerYs[idx] = (centerYs[idx] - minY) / rangeY;

  // inflate centroids
  inflateCentroids(centerXs, centerYs, clusterSizes);

  // map the points to the center of the cluster
  for (let idx = 0; idx < cluster.length; idx++) {
    const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
    const size = clusterSizes[clusterNum];
    const centerX = centerXs[clusterNum];
    const centerY = centerYs[clusterNum];
    const clusterMidX = (minXs[clusterNum] + maxXs[clusterNum]) / 2;
    const clusterMidY = (minYs[clusterNum] + maxYs[clusterNum]) / 2;
    const dx = (embedX[idx] - clusterMidX) / rangeXs[clusterNum];
    const dy = (embedY[idx] - clusterMidY) / rangeYs[clusterNum];
    // make each cluster more circular
    const angle = Math.atan2(dy, dx);
    const dist = Math.sqrt(dx * dx + dy * dy);
    const maxDist = size / 2;

    embedX[idx] = centerX + Math.cos(angle) * Math.min(maxDist, dist);
    embedY[idx] = centerY + Math.sin(angle) * Math.min(maxDist, dist);

    //    embedX[idx] = centerX + dx * size;
    //    embedY[idx] = centerY + dy * size;
  }

  // put all points in each cluster in center
  // for (let idx = 0; idx < cluster.length; idx++) {
  //   const clusterNum = clusterSet.get(subCluster[cluster[idx]])!;
  //   embedX[idx] = centerXs[clusterNum];
  //   embedY[idx] = centerYs[clusterNum];
  // }


  // normalize
  minX = embedX[0];
  minY = embedY[0];
  maxX = embedX[0];
  maxY = embedY[0];
  for (let idx = 1; idx < cluster.length; idx++) {
    minX = Math.min(minX, embedX[idx]);
    minY = Math.min(minY, embedY[idx]);
    maxX = Math.max(maxX, embedX[idx]);
    maxY = Math.max(maxY, embedY[idx]);
  }
  rangeX = maxX - minX;
  rangeY = maxY - minY;
  if (rangeX === 0)
    rangeX = maxX === 0 ? 1 : maxX;
  for (let idx = 0; idx < cluster.length; idx++)
    embedX[idx] = (embedX[idx] - minX) / rangeX;
  if (rangeY === 0)
    rangeY = maxY === 0 ? 1 : maxY;
  for (let idx = 0; idx < cluster.length; idx++)
    embedY[idx] = (embedY[idx] - minY) / rangeY;

  return {embedX, embedY};
}

function sign(x: number) {
  return x < 0 ? -1 : 1;
}

// this function will inflate the centroids of the clusters such that their size is proportional to the number of points in the cluster
// sizes should be normalized to [0, 1]
function inflateCentroids(centerXs: Float32Array, centerYs: Float32Array, clusterSizes: Float32Array) {
  const sortedIndexes = Array.from(clusterSizes.keys()).sort((a, b) => clusterSizes[b] - clusterSizes[a]);
  const numClusters = centerXs.length;
  for (let uidx = 0; uidx < numClusters; uidx++) {
    const idx = sortedIndexes[uidx];
    const size = clusterSizes[idx];
    const centerX = centerXs[idx];
    const centerY = centerYs[idx];
    for (let ui = 0; ui < numClusters; ui++) {
      const i = sortedIndexes[ui];
      const size2 = clusterSizes[i];
      const minDist = (size + size2);
      if (i === idx)
        continue;
      const dx = centerX - centerXs[i];
      const dy = centerY - centerYs[i];
      const dist = Math.sqrt(dx * dx + dy * dy);
      if (dist < minDist) {
        const diff = size / 2;
        const angle = Math.atan2(dy, dx);
        centerXs[i] -= diff * Math.cos(angle);
        centerYs[i] -= diff * Math.sin(angle);
      }
    }
  }
}
