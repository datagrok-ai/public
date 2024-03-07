import {SparseMatrixObject} from './types';

export async function bioLayout(cluster: number[], sparseObject: SparseMatrixObject, threshold: number) {
  const n = cluster.length;
  const is: number[] = [];
  const js: number[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a = cluster[i];
      const b = cluster[j];
      if (sparseObject[a]?.[b] >= threshold) {
        is.push(i);
        js.push(j);
      }
    }
  }
  return getLayoutEmbeddings(is, js, cluster);
}

function getLayoutEmbeddings(is: number[], js: number[], cluster: number[]) {
  const embedX = new Float32Array(cluster.length).fill(0).map(() => Math.random() * 10);
  const embedY = new Float32Array(cluster.length).fill(0).map(() => Math.random() * 10);
  const iterations = 100;
  const xVelocities = new Float32Array(cluster.length).fill(0);
  const yVelocities = new Float32Array(cluster.length).fill(0);
  for (let i = 0; i < iterations; i++) {
    const temperature = (1 - i / iterations);
    xVelocities.fill(0);
    yVelocities.fill(0);
    for (let idx = 0; idx < is.length; idx++) {
      const a = is[idx];
      const b = js[idx];
      const dx = embedX[a] - embedX[b];
      const dy = embedY[a] - embedY[b];
      //   const distance = Math.sqrt(dx * dx + dy * dy);
      //   const factor = (distance - 1) / distance;
      //   const offsetX = dx * factor;
      //   const offsetY = dy * factor;
      if (Math.abs(dx) >= 1) {
        xVelocities[a] -= temperature * dx;
        xVelocities[b] += temperature * dx;
      }
      if (Math.abs(dy) >= 1) {
        yVelocities[a] -= temperature * dy;
        yVelocities[b] += temperature * dy;
      }
    //   embedY[a] += offsetY * temperature;
    //   embedY[b] -= offsetY * temperature;
    }
    for (let idx = 0; idx < cluster.length; idx++) {
      const vecSize = Math.sqrt(xVelocities[idx] * xVelocities[idx] + yVelocities[idx] * yVelocities[idx]);
      if (vecSize > 0) {
        embedX[idx] += xVelocities[idx] / vecSize * temperature;
        embedY[idx] += yVelocities[idx] / vecSize * temperature;
      }
    }
  }

  // normalize
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
    rangeX = maxX;
  for (let idx = 0; idx < cluster.length; idx++)
    embedX[idx] = (embedX[idx] - minX) / rangeX / 2 + 0.5;
  if (rangeY === 0)
    rangeY = maxY;
  for (let idx = 0; idx < cluster.length; idx++)
    embedY[idx] = (embedY[idx] - minY) / rangeY / 2 + 0.5;

  return {embedX, embedY};
}
