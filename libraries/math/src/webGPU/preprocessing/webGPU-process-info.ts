/* eslint-disable max-len */
import {BioDistanceFnOptions, distanceFunctionComplexity,
  SupportedEntryTypes, TypeSupportedDistances, WEBGPUDISTANCE, WGPUENTRYTYPE} from '../multi-col-distances/webGPU-multicol-distances';

export function webGPUProcessInfo(
  entryList: Array<SupportedEntryTypes>,
  distanceMetric: WEBGPUDISTANCE = WEBGPUDISTANCE.HAMMING,
  entryIndex: number, // index of the entries in the list of lists that we want to process
  options: BioDistanceFnOptions & {[key: string]: any} = {gapOpenPenalty: 1.0, gapExtensionPenalty: 0.6}
) {
  let entryType: WGPUENTRYTYPE | null = null;
  const encodedList = (() => {
    if (entryList.some((e) => typeof e === 'string')) {
      entryType = WGPUENTRYTYPE.STRING;
      return entryList.map((entry) => new Uint32Array((entry as string).split('').map((c) => c.charCodeAt(0))));
    }
    if (entryList.some((e) => typeof e === 'number')) {
      entryType = WGPUENTRYTYPE.NUMBER;
      return entryList.map((entry) => new Float32Array([entry as number]));
    }
    if (typeof entryList[0] == 'object' && entryList.some((e) => '_data' in (e as Object) && '_length' in (e as Object))) {
      entryType = WGPUENTRYTYPE.BITARRAY;
      return entryList.map((entry) => (entry as {_data: Uint32Array, _length: number})._data);
    }
    if (entryList.some((e) => e instanceof Float32Array)) {
      entryType = WGPUENTRYTYPE.FLOAT32ARRAY;
      return entryList as Float32Array[];
    }
    if (entryList.some((e) => e instanceof Uint32Array)) {
      entryType = WGPUENTRYTYPE.UINT32ARRAY;
      return entryList as Uint32Array[];
    }
    if (entryList.some((e) => e instanceof Int32Array)) {
      entryType = WGPUENTRYTYPE.INT32ARRAY;
      return entryList as Int32Array[];
    }
    //return entryList as Uint32Array[];
  })();

  if (!encodedList || !entryType)
    throw new Error('Invalid entry type, could not determine entry type from input list');

  const encodedListType = encodedList[0] instanceof Int32Array ? WGPUENTRYTYPE.INT32ARRAY :
    encodedList[0] instanceof Float32Array ? WGPUENTRYTYPE.FLOAT32ARRAY : WGPUENTRYTYPE.UINT32ARRAY;
  // sizes of each entries might differ, so we need to keep track of that for some distance metrics, like hamming for example
  const arraySizes = new Uint32Array(encodedList.map((arr) => arr.length));

  if (!TypeSupportedDistances[entryType] || !TypeSupportedDistances[entryType].has(distanceMetric))
    throw new Error(`Distance metric '${distanceMetric}' not supported for entry type '${entryType}'`);


  const maxEntryLen = arraySizes.reduce((a, b) => Math.max(a, b), 0);
  // get the complexity of used algorithm
  const complexity = distanceFunctionComplexity[distanceMetric](maxEntryLen);

  const EncodedArrayConstructor = encodedListType === WGPUENTRYTYPE.INT32ARRAY ? Int32Array :
    (encodedListType === WGPUENTRYTYPE.FLOAT32ARRAY ? Float32Array : Uint32Array);

  const flatSourceArray = new EncodedArrayConstructor(encodedList.length * maxEntryLen);

  // when setting, we need to set each array at a specific offset, which is controlled by maxArrayLen because each array might have different sizes.
  // this way we will get correct matrix representation in the compute shader
  encodedList.forEach((seq, i) => {
    flatSourceArray.set(seq, i * maxEntryLen);
  });

  // NB! all this before the line was generic, now we need to calculate some specific things for some specific distance metrics
  // initialize supp info line that will be included in the final shader;
  let suppInfoStructWgsl: string = ''; // the code that will be included in the struct of suppInfo
  let suppInfoSize: number = 0;
  let suppInfoType: WGPUENTRYTYPE.FLOAT32ARRAY | WGPUENTRYTYPE.UINT32ARRAY = WGPUENTRYTYPE.FLOAT32ARRAY;
  let suppInfoBuffer: Float32Array | Uint32Array | null = null;
  if (distanceMetric === WEBGPUDISTANCE.NEEDLEMAN_WUNSCH || distanceMetric === WEBGPUDISTANCE.MONOMER_CHEMICAL_DISTANCE) {
    let maxMonomerIndex = options.scoringMatrix && options.alphabetIndexes ?
      Object.keys(options.alphabetIndexes).reduce((prev, n) => Math.max(prev, n.charCodeAt(0)), 0) : -1;
    // generate default similarity matrix if it is not provided
    if (!options.alphabetIndexes || !options.scoringMatrix) {
      for (let i = 0; i < flatSourceArray.length; i++) {
        if (flatSourceArray[i] > maxMonomerIndex)
          maxMonomerIndex = flatSourceArray[i];
      }
      options.scoringMatrix =
                new Array(maxMonomerIndex + 1).fill(null).map(() => new Array(maxMonomerIndex + 1).fill(0));
      options.alphabetIndexes = {};

      for (let i = 0; i < options.scoringMatrix.length; i++) {
        options.scoringMatrix[i][i] = 1;
        options.alphabetIndexes[String.fromCharCode(i)] = i;
      }
    }
    const similarityMatrixSize = (maxMonomerIndex + 1) * (maxMonomerIndex + 1);
    const transferedSimilarityMatrix = new Array(maxMonomerIndex + 1).fill(null).map(() => new Float32Array(maxMonomerIndex + 1));
    // set diagonal to 1
    for (let i = 0; i < maxMonomerIndex + 1; i++)
      transferedSimilarityMatrix[i][i] = 1;

    const alphabetIndexes = options.alphabetIndexes!;
    for (const key of Object.keys(alphabetIndexes)) {
      for (const key2 of Object.keys(alphabetIndexes)) {
        if (key === key2)
          continue;

        transferedSimilarityMatrix[key.charCodeAt(0)][key2.charCodeAt(0)] =
                options.scoringMatrix![alphabetIndexes[key]][alphabetIndexes[key2]];
      }
    }

    // in memory layout, we will have 2 float32 s for gapOpen and gapExtension penalties, and then f32 array<array<f32>> for similarity matrix.
    // because of primitives, there will be no padding, so we can calculate the size directly
    suppInfoSize = 2 + similarityMatrixSize;
    suppInfoType = WGPUENTRYTYPE.FLOAT32ARRAY;
    suppInfoBuffer = new Float32Array(suppInfoSize);
    suppInfoBuffer[0] = options.gapOpenPenalty ?? 1.0;
    suppInfoBuffer[1] = options.gapExtensionPenalty ?? 0.6;
    let offset = 2;
    for (let i = 0; i < transferedSimilarityMatrix.length; i++) {
      suppInfoBuffer.set(transferedSimilarityMatrix[i], offset);
      offset += transferedSimilarityMatrix[i].length;
    }
    suppInfoStructWgsl = `
            gapOpenPenalty${entryIndex}: f32,
            gapExtensionPenalty${entryIndex}: f32,
            similarityMatrix${entryIndex}: array<array<f32, ${maxMonomerIndex + 1}>, ${maxMonomerIndex + 1}>`;
  } else if (distanceMetric === WEBGPUDISTANCE.Difference) {
    // for difference, we need range of values for normalization of the difference
    if (!options.range || typeof options.range !== 'number' || options.range <= 0) {
      const min = (flatSourceArray as Float32Array).reduce((a, b) => Math.min(a, b), flatSourceArray[0]);
      const max = (flatSourceArray as Float32Array).reduce((a, b) => Math.max(a, b), flatSourceArray[0]);
      options.range = max - min;
    }
    if (options.range <= 0)
      throw new Error('Invalid range for difference distance metric: ' + options.range);

    suppInfoSize = 1;
    suppInfoType = WGPUENTRYTYPE.FLOAT32ARRAY;
    suppInfoBuffer = new Float32Array([options.range]);
    suppInfoStructWgsl = `
            range${entryIndex}: f32`;
  }

  const dataTypeWGSL = flatSourceArray instanceof Int32Array ? 'i32' : (flatSourceArray instanceof Float32Array ? 'f32' : 'u32');
  const dataStructWgsl = `data${entryIndex}: array<array<${dataTypeWGSL}, ${maxEntryLen}>, ${encodedList.length}>`;

  // for now, other distances do not require any additional information, so we can skip that
  return {
    flatSourceArray,
    sourceArraySize: flatSourceArray.length,
    maxEntryLen,
    arraySizes,
    complexity,
    suppInfoBuffer,
    suppInfoSize,
    suppInfoType: suppInfoType as WGPUENTRYTYPE.FLOAT32ARRAY | WGPUENTRYTYPE.UINT32ARRAY,
    suppInfoStructWgsl,
    entryType,
    dataTypeWGSL,
    dataStructWgsl,
    EncodedArrayConstructor
  };
}

export type WebGPUProcessInfo = ReturnType<typeof webGPUProcessInfo>;
