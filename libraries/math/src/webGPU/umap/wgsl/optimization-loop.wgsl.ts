/* eslint-disable max-len */
export function optimizeLayoutWGSL(
  headSize: number, embedSize: number, moveOther: boolean, initialAlpha: number, gamma: number,
  a: number, b: number, dim: number, nEpochs: number, nVertices: number, workgroupDim: number, threadsPerRow: number,
  divisionFactor = 100000, stepsPerEpoch = 1, threadsPerStep = headSize
) {
  return `
// first define some functions

fn rDist(ar1: array<f32, ${dim}>, ar2: array<f32, ${dim}>) -> f32 {
    var res = 0.0;
    for (var i1 = 0u; i1 < ${dim}; i1++) {
        let diff = ar1[i1] - ar2[i1];
        res += pow(diff, 2);
    }
    return res;
}

fn clip(x: f32, clipValue: f32) -> f32 {
    if (x > clipValue) {
        return clipValue;
    }
    else if (x < 0.0 - clipValue) {
        return 0.0 - clipValue;
    }
    return x;
}

var<private> rngState: u32 = 0;
var<private> nVertices: u32 = ${nVertices};
var<private> initialAlpha: f32 = ${initialAlpha};
var<private> gamma: f32 = ${gamma};
var<private> a: f32 = ${a};
var<private> b: f32 = ${b};
var<private> dim: u32 = ${dim};
var<private> nEpochs: f32 = ${nEpochs};
var<private> moveOther: u32  = ${moveOther === true ? 1 : 0};

// PCG hash (Jarzynski and Olano, Hash Functions for GPU Rendering, 2020)
fn pcgHash(v: u32) -> u32 {
    let state = v * 747796405u + 2891336453u;
    let word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

fn randInt(maxVal: u32) -> u32 {
    rngState = pcgHash(rngState);
    return rngState % maxVal;
}

struct MatrixStorage {
    head: array<u32, ${headSize}>,
    tail: array<u32, ${headSize}>,
}

struct EpochStorage {
    epochsPerSample: array<f32, ${headSize}>,
    epochOfNextSample: array<f32, ${headSize}>,
    epochsPerNegativeSample: array<f32, ${headSize}>,
    epochOfNextNegativeSample: array<f32, ${headSize}>
}

struct EmbeddingStorage {
    // embeddings will be in range of -10 to 10, as atomics only can store i32 or u32, we will store them as floor(f32 * 100_000)
    headEmbeddings: array<array<atomic<i32>, ${dim}>, ${embedSize}>,
    tailEmbeddings: array<array<atomic<i32>, ${dim}>, ${embedSize}>
}

// the embeddings as they were when the step started. Threads read only these and add their moves to EmbeddingStorage:
// read from EmbeddingStorage, a thread would see whichever moves of other threads happened to land first, and the
// same seed would give other coordinates on every run. Integer sums do not depend on the order of the adds.
struct EmbeddingSnapshot {
    headEmbeddings: array<array<i32, ${dim}>, ${embedSize}>,
    tailEmbeddings: array<array<i32, ${dim}>, ${embedSize}>
}

struct ComputeInfo {
    n: f32,
    alpha: f32,
    seed: u32,
}
    @group(0) @binding(0) var<storage, read_write> matrixStorage: MatrixStorage;
    @group(0) @binding(1) var<storage, read_write> epochStorage: EpochStorage;
    @group(0) @binding(2) var<storage, read_write> embedStorage: EmbeddingStorage;
    @group(0) @binding(3) var<storage, read_write> computeInfo: ComputeInfo;
    @group(0) @binding(4) var<storage, read> embedSnapshot: EmbeddingSnapshot;
    @group(0) @binding(5) var<storage, read> step: u32;
    @compute @workgroup_size(${workgroupDim}, ${workgroupDim}) fn optimizeStep(
        @builtin(global_invocation_id) id: vec3<u32>,
     ) {
        let divisionFactor: f32 = ${divisionFactor};
        let col: u32 = id.x;
        let row: u32 = id.y;
        let threadIndex: u32 = row * ${threadsPerRow} + col;
        let clipValue: f32 = 4.0;

        if (threadIndex >= ${threadsPerStep}u) {
            return;
        }

        let workingIndex = threadIndex * ${stepsPerEpoch}u + step;
        if (workingIndex >= ${headSize}u) {
            return;
        }
        //let startAtOffset = computeInfo.tailOffsets[threadIndex];
        //let endAtOffset = computeInfo.tailOffsets[threadIndex + 1];
        //for (var workingIndex = startAtOffset; workingIndex < endAtOffset; workingIndex++) {
            // a stream of its own for every edge and epoch, all from the one seed
            rngState = pcgHash(workingIndex + pcgHash(u32(computeInfo.n) + pcgHash(computeInfo.seed)));

            let epochOfNextSample = epochStorage.epochOfNextSample[workingIndex];
            if (epochOfNextSample > computeInfo.n) {
                return;
            }

            let j = matrixStorage.head[workingIndex];
            let k = matrixStorage.tail[workingIndex];

            // as said before, we will store embeddings as floor(f32 * 100_000)
            var current: array<f32, ${dim}>;
            var other: array<f32, ${dim}>;
            for (var i = 0u; i < ${dim}u; i++) {
                current[i] = f32(embedSnapshot.headEmbeddings[j][i]) / divisionFactor;
                other[i] = f32(embedSnapshot.tailEmbeddings[k][i]) / divisionFactor;
            }

            let distSquared = rDist(current, other);
            
            var gradCoeff: f32 = 0.0;
            if (distSquared > 0.0) {
                gradCoeff = (0.0 - 2.0) * ${a} * ${b} * pow(distSquared, b - 1.0);
                gradCoeff /= (a * pow(distSquared, b) + 1.0);
            }

            var gradBuff: array<i32, ${dim}>;
            for (var d = 0u; d < ${dim}u; d++) {
                let gradD = clip((gradCoeff * (current[d] - other[d])), clipValue);
                let toAdd = gradD * computeInfo.alpha;
                gradBuff[d] = i32(floor(toAdd * divisionFactor));
                current[d] += toAdd;
                other[d] -= toAdd;
            }

            for (var d = 0u; d < ${dim}u; d++) {
                atomicAdd(&embedStorage.headEmbeddings[j][d], gradBuff[d]);
                atomicSub(&embedStorage.tailEmbeddings[k][d], gradBuff[d]);
            }

            epochStorage.epochOfNextSample[workingIndex] += epochStorage.epochsPerSample[workingIndex];
            var nNegSamples: i32 = i32(floor((computeInfo.n - epochStorage.epochOfNextNegativeSample[workingIndex]) / epochStorage.epochsPerNegativeSample[workingIndex]));

            for (var p = 0i; p < nNegSamples; p++) {
                let k1 = randInt(nVertices);
                var other1: array<f32, ${dim}>;
                for (var i = 0u; i < ${dim}u; i++) {
                    other1[i] = f32(embedSnapshot.tailEmbeddings[k1][i]) / divisionFactor;
                }

                let distSquared1 = rDist(current, other1);

                var gradCoeff1: f32 = 0.0;
                if (distSquared1 > 0.0) {
                    gradCoeff1 = 2.0 * gamma * b;
                    gradCoeff1 /= (0.001 + distSquared1) * (a * pow(distSquared1, b) + 1);
                }
                else if (j == k1) {
                    continue;
                }

                var gradDBuff: array<i32, ${dim}>;
                for (var d = 0u; d < ${dim}u; d++) {
                    var gradD: f32 = 4.0;
                    if (gradCoeff1 > 0.0) {
                        gradD = clip(gradCoeff1 * (current[d] - other1[d]), clipValue);
                    }

                    gradDBuff[d] = i32(floor(gradD * computeInfo.alpha * divisionFactor));
                    
                }
                for (var d = 0u; d < ${dim}u; d++) {
                    atomicAdd(&embedStorage.headEmbeddings[j][d], gradDBuff[d]);
                }
            }
            epochStorage.epochOfNextNegativeSample[workingIndex] += f32(nNegSamples) * epochStorage.epochsPerNegativeSample[workingIndex];
    //}
    }
`;
}
