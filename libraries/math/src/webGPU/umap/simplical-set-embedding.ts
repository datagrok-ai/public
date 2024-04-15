import {performMatrixOps} from './fuzzy-simplical-set';
import {findABParams, getNEpochs, toOffsetForm} from './utils';

export async function initializeSimplicialSetEmbedding(
  graph: NonNullable<Awaited<ReturnType<typeof performMatrixOps>>>,
  entryLen: number, spread: number, minDist: number,
  negativeSampleRate: number = 5
) {
  if (!graph!.unionSizes)
    return;

  const nEpochs = getNEpochs(entryLen);
  const nComponents = 2; // todo: support mult dims
  // get max value:
  const maxVal = graph.res!.resultKnnDistances.reduce((prev, cur) => Math.max(prev, cur));

  //let offsetAccum = 0;
  // for (let i = 0; i < graph.unionSizes.length; i++) {
  //     for (let j = 0; j < graph.unionSizes)

  // }
  const minGraphVal = maxVal / nEpochs;
  // this will hold the sizes of arrays with all zeros removed.
  const unionSizesCopy = new Uint32Array(graph.unionSizes.length);
  unionSizesCopy.set(graph.unionSizes);

  let offsetAccum = 0;
  for (let i = 0; i < graph.unionSizes.length; i++) {
    for (let j = 0; j < graph.unionSizes[i]; j++) {
      const index = offsetAccum + j;
      if (graph.res!.resultKnnDistances[index] < minGraphVal) {
                graph.res!.resultKnnDistances[index] = 0;
                unionSizesCopy[i]--;
      }
    }
    offsetAccum += graph.unionSizes[i];
  }


  const nonZeroGraphOffsets = toOffsetForm(unionSizesCopy);
  const nonZeroGraphSize = nonZeroGraphOffsets[nonZeroGraphOffsets.length - 1];
  const weights = new Float32Array(nonZeroGraphSize);
  const head = new Uint32Array(nonZeroGraphSize);
  const tail = new Uint32Array(nonZeroGraphSize);

  offsetAccum = 0;
  for (let i = 0; i < graph.unionSizes.length; i++) {
    for (let j = 0; j < graph.unionSizes[i]; j++) {
      const index = graph.unionMatrixOffsets[i] + j;
      if (graph.res!.resultKnnDistances[index] != 0) {
        weights[offsetAccum] = graph.res!.resultKnnDistances[index];
        head[offsetAccum] = graph.res!.resultKnnIndexes[index];
        tail[offsetAccum] = i;
        offsetAccum +=1;
      }
    }
  }

  const epochsPerSample = new Float32Array(weights.length).fill(-1);
  const nSamples = weights.map((w) => (w / maxVal) * nEpochs);
  nSamples.forEach((n, i) => {
    if (n > 0) epochsPerSample[i] = nEpochs / nSamples[i];
  });

  const {a, b} = findABParams(spread, minDist);
  const epochsPerNegativeSample = epochsPerSample.map(
    (e) => e / negativeSampleRate
  );

  const gamma = 1.0;
  const initialAlpha = 1.0;
  return {epochsPerSample, epochsPerNegativeSample, a, b,
    head, tail, entryLen, nEpochs, nComponents, initialAlpha, gamma};
}
