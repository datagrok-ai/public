export function smoothKNNDistanceWGSL(
  threadsPerWorkgroupDim: number = 10, threadsPerYDimension: number,
  knnDistances: number[][] | Float32Array[], nNeighbors: number = 15,
  localConnectivity: number = 1.0,
  nIter: number = 64, bandwidth: number = 1.0
) {
  const numOfEntries = knnDistances.length;
  let meanDistancesSum: number = 0.0;
  for (let i = 0; i < numOfEntries; i++) {
    let currentMean: number = 0.0;
    for (let j = 0; j < nNeighbors; j++)
      currentMean += knnDistances[i][j];
    meanDistancesSum += currentMean / nNeighbors;
  }
  const meanDistances = meanDistancesSum / numOfEntries;
  const target = Math.log(nNeighbors) / Math.log(2) * bandwidth;
  return `  
        var<private> SMOOTH_K_TOLERANCE: f32 = 1e-5;
        var<private> MIN_K_DIST_SCALE: f32 = 1e-3;
        var<private> nNeighbors: u32 = ${nNeighbors};
        var<private> localConnectivity: f32 = ${localConnectivity}; 
        var<private> nIter: u32 = ${nIter};
        var<private> bandwidth: u32 = ${bandwidth};
        var<private> meanDistances: f32 = ${meanDistances};
        var<private> targetValue: f32 = ${target};
        @group(0) @binding(0) var<storage, read_write> sigmas: array<f32, ${numOfEntries}>;
        @group(0) @binding(1) var<storage, read_write> rhos: array<f32, ${numOfEntries}>;
        @group(0) @binding(2) var<storage, read> distances: array<array<f32, ${nNeighbors}>, ${numOfEntries}>;
        @compute @workgroup_size(${threadsPerWorkgroupDim}, ${threadsPerWorkgroupDim}) fn smoothKNNDistance(
            @builtin(global_invocation_id) id: vec3<u32>,
        ) {
            let col: u32 = id.x;
            let row: u32 = id.y;
            let workingIndex: u32 = row * ${threadsPerYDimension} + col;

            if (workingIndex >= ${numOfEntries}) {
                return;
            }

            var lo: f32 = 0.0;
            var hi: f32 = 1.0e+30; // treated as Infinity
            var mid: f32 = 1.0;

            var nonZeroDists: array<f32, ${nNeighbors}>;
            var nonZeroDistsSize: u32 = 0;
            for (var j = 0u; j < nNeighbors; j = j + 1u) {
                if (distances[workingIndex][j] > 0.0) {
                    nonZeroDists[nonZeroDistsSize] = distances[workingIndex][j];
                    nonZeroDistsSize = nonZeroDistsSize + 1u;
                }
            }

            if (f32(nonZeroDistsSize) >= localConnectivity) {
                let index: u32 = u32(floor(localConnectivity));
                let interpolation: f32 = localConnectivity - f32(index);
                if (index > 0) {
                    rhos[workingIndex] = nonZeroDists[index - 1];
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rhos[workingIndex] += interpolation * (nonZeroDists[index] - nonZeroDists[index - 1]);
                    }
                }
                else {
                    rhos[workingIndex] = interpolation * nonZeroDists[0];
                }
            }
            else if (nonZeroDistsSize > 0u) {
                var maxDist: f32 = 0.0;
                for (var j = 0u; j < nonZeroDistsSize; j = j + 1u) {
                    maxDist = max(nonZeroDists[j], maxDist);
                }
                rhos[workingIndex] = maxDist;
            }

            for (var n = 0u; n < nIter; n = n + 1u) {
                var psum: f32 = 0.0;
                for (var j = 1u; j < ${nNeighbors}; j = j + 1u) {
                    let d: f32 = distances[workingIndex][j] - rhos[workingIndex];
                    if (d > 0.0) {
                        psum += exp(0.0 - (d / mid));
                    }
                    else {
                        psum += 1.0;
                    }
                }

                if (abs(psum - targetValue) < SMOOTH_K_TOLERANCE) {
                    break;
                }

                if (psum > targetValue) {
                    hi = mid;
                    mid = (lo + hi) / 2.0;
                }
                else {
                    lo = mid;
                    if (hi == 1.0e+30) {
                        mid *= 2.0;
                    }
                    else {
                        mid = (lo + hi) / 2.0;
                    }
                }
            }

            sigmas[workingIndex] = mid;

            if (rhos[workingIndex] > 0.0) {
                var sum: f32 = 0.0;
                for (var j = 0u; j < ${nNeighbors}; j = j + 1u) {
                    sum += distances[workingIndex][j];
                }
                let meanIthDistances: f32 = sum / ${nNeighbors}.0;
                if (sigmas[workingIndex] < MIN_K_DIST_SCALE * meanIthDistances) {
                    sigmas[workingIndex] = MIN_K_DIST_SCALE * meanIthDistances;
                }
            }
            else {
                if (sigmas[workingIndex] < MIN_K_DIST_SCALE * meanDistances) {
                    sigmas[workingIndex] = MIN_K_DIST_SCALE * meanDistances;
                }
            }
        }
    `;
}
