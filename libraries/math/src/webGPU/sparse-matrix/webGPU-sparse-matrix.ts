/* eslint-disable max-len */
import {SparseMatrix} from '../types';
import {
  WEBGSLAGGREGATION,
  WEBGSLAGGREGATIONFUNCTIONS,
} from '../multi-col-distances/webGPU-aggregation';
import {
  SupportedEntryTypes,
  WEBGPUDISTANCE,
  webGPUFunctions,
  WGPUENTRYTYPE,
} from '../multi-col-distances/webGPU-multicol-distances';
import {webGPUProcessInfo} from '../preprocessing/webGPU-process-info';
import {getGPUDevice} from '../getGPUDevice';

/** generate sparse matrix based on list of lists of entries.
 *  these entries are each encoded as Uint32Array or FLOAT32Array (depending on their type).
 * for example, sequences would be encoded as Uint32Array based on char code of the letter at each position.
 * [65, 66, 67, 68, 69] would be a sequence of 5 letters.
 * for chemical fingerprints, it would be a binary array of 0s and 1s,
 * represented as Uint32Array(_data property of DG bitarray).
 *
 * Be ware that size of entryList, distanceMetrics, weights and options must be the same.
 * if there are no options for entries i, pass an empty object.
 * for now options are needed for
 * needleman-wunsch and monomer chemical distances: see {@link BioDistanceFnOptions} as for how it should be passed
 * numeric distances (Difference): {range: number} where range is the range of the values in the column (max - min).
 * in both cases, if options are not provided, they will be calculated automatically.
 */
export async function multiColWebGPUSparseMatrix(
  entryList: Array<Array<SupportedEntryTypes>>, // list of lists of entries, for multiple columns
  threshold: number = 0.8, // similarity threshold, be ware that if you use too small threshold, there might be memory overflow...
  distanceMetrics: WEBGPUDISTANCE[], // distance metrics for each column
  aggregationFunction: WEBGSLAGGREGATION, // aggregation function for the distances
  weights: number[], // weights for each column
  options: { [key: string]: any }[] // supplementary options for each column
): Promise<SparseMatrix | null> {
  const device = await getGPUDevice();
  if (!device)
    return null; // if no device, return null, as we cannot do anything without it.

  const availableDistanceMetrics = Object.values(WEBGPUDISTANCE);
  if (distanceMetrics.some((metric) => !availableDistanceMetrics.includes(metric)))
    throw new Error('Invalid distance metrics provided: ' + distanceMetrics.join(', '));

  const availableAggregationFunctions = Object.values(WEBGSLAGGREGATION);
  if (!availableAggregationFunctions.includes(aggregationFunction))
    throw new Error('Invalid aggregation function provided: ' + aggregationFunction);

  const maxDistance = 1 - threshold; // maximum distance

  // first, check that all the supplementary options are provided and are the same length:
  if (
    options.length !== entryList.length ||
    options.length !== distanceMetrics.length ||
    options.length !== weights.length
  ) {
    throw new Error(
      'Options, weigths and distance functions must be provided for each column'
    );
  }

  // check that all the entry lists are the same length
  if (entryList.some((list) => list.length !== entryList[0].length))
    throw new Error('All entry lists must be the same length');

  const numOfColumns = entryList.length; // number of columns
  const listSize = entryList[0].length; // size of each list (or column)
  const processInfo = entryList.map((entry, i) => {
    return webGPUProcessInfo(entry, distanceMetrics[i], i, options[i]);
  });

  if (numOfColumns === 0) {
    throw new Error(
      'No columns provided. Please provide at least one column of data.'
    );
  }

  if (numOfColumns === 1) aggregationFunction = WEBGSLAGGREGATION.MANHATTAN; // save a bit of time

  // combine all struct types into one to put into the suppInfo struct.
  let suppInfoWgsl = processInfo
    .map((info) => info.suppInfoStructWgsl)
    .filter((wgsl) => !!wgsl && wgsl != '')
    .join(',\n');
  // structures in wgsl must have at least one member, so if we have no structures, we need to add a dummy one
  let needsDummy = false;
  if (!suppInfoWgsl || suppInfoWgsl.trim() == '') {
    needsDummy = true;
    suppInfoWgsl = '\ndummy: f32\n';
  }

  // combine all data wgsl struct code into one
  const dataWgsl = processInfo.map((info) => info.dataStructWgsl).filter((wgsl) => !!wgsl && wgsl != '').join(',\n');
  // combine all array sizes into one array (easier for setting)
  const arraySizes = new Uint32Array(numOfColumns * listSize);
  processInfo.forEach((info, i) => {
    arraySizes.set(info.arraySizes, i * listSize);
  }); // array.flat is not as optimized as this

  // if we try to map large arrays directly from GPU, sometimes, device disconnects. so we need to do it in chunks, a good number
  // we found is 10000. So we will perform computations in chunks of 10000. meaning that we will dispatch 10000 threads at a time.
  const numOfThreads = 10000;

  // in this case we do not need to worry about complexity of the algorithm, as the 100 is low enaugh number, which is limited by memory usage.
  const sparseResultSizePerThread = 100; // number of iterations per thread (number of pair comparisons)

  const combinedComplexity = processInfo.reduce((a, b) => a + b.complexity, 0); // combined complexity of all the columns
  const maxIterationsPerThread = Math.ceil(6000 / combinedComplexity); // maximum number of iterations per thread
  const workGroupDivision = 10; // how many threads inside of one workgroup dimension (in this case 10 * 10 threads per workgroup)
  const threadsPerWorkgroup = workGroupDivision * workGroupDivision;
  const workgroupsDim = Math.ceil(
    Math.sqrt(Math.ceil(numOfThreads / threadsPerWorkgroup))
  ); // how many workgroups per 2d dimension
  const globalThreadDimSize = workgroupsDim * workGroupDivision; // how many threads per 2d dimension


  const condensedDistanceMatrixSize = listSize * (listSize - 1) / 2; // size of the condensed distance matrix, this many comparisons will be made.
  const dmChunkSizePerThread = Math.ceil(condensedDistanceMatrixSize / numOfThreads); // how many comparisons per thread
  const module = device.createShaderModule({
    label: 'Sparse matrix compute shader',
    code: `
    // each thread will perform ${sparseResultSizePerThread} iterations at one time, comparing ${sparseResultSizePerThread} pairs of entries.
    // in total, each thread will perform at most ${dmChunkSizePerThread} comparisons.
    // first is the result struct, containing is, js, and distances. each array with length of ${sparseResultSizePerThread},
    // and also integer for how many pairs were found to be below threshold.
    struct SparseResult {
      i: array<array<u32, ${sparseResultSizePerThread}>, ${numOfThreads}>,
      j: array<array<u32, ${sparseResultSizePerThread}>, ${numOfThreads}>,
      distances: array<array<f32, ${sparseResultSizePerThread}>, ${numOfThreads}>,
      found: array<u32, ${numOfThreads}>,
      done: array<u32, ${numOfThreads}>
    }
    // struct for the data
    struct ComputeInfo {
      // start at cols and rows, and end at cols and rows for each thread, these will be calculated on cpu and passed to gpu.
      startAtCols: array<u32, ${numOfThreads}>,
      startAtRows: array<u32, ${numOfThreads}>,
      endAtCols: array<u32, ${numOfThreads}>,
      endAtRows: array<u32, ${numOfThreads}>,

      // the ACTUALLY sizes of each entry
      entrySizes: array<array<u32, ${listSize}>, ${numOfColumns}>,
      // the weights for each entry
      weights: array<f32, ${numOfColumns}>,
      // the data for each entry
      ${dataWgsl} // an example of the dataWgsl would be:
      //data0: array<array<u32,20>,100>,
      //data1: array<array<u32,20>,100>
    }

    // struct for the supplementary information
    struct SuppInfo {
      // struct containing all the supplementary info, like scoring matrix, alphabet indexes, range, etc.
      ${suppInfoWgsl}
    };

    @group(0) @binding(0) var<storage, read_write> computeInfo: ComputeInfo;
    @group(0) @binding(1) var<storage, read_write> suppInfo: SuppInfo;
    @group(0) @binding(2) var<storage, read_write> results: SparseResult;
    @compute @workgroup_size(${workGroupDivision}, ${workGroupDivision}) fn calcSparseMatrix(
      @builtin(global_invocation_id) id: vec3<u32>
    ) {
      ${needsDummy ? `let otherDummy = suppInfo.dummy * 2;` : ''} // just to make sure that the suppInfo is not optimized out
      let threadCol = id.x;
      let threadRow = id.y;
      let linearIndex = threadRow * ${globalThreadDimSize} + threadCol;
      if (linearIndex >= ${numOfThreads}) {
        return; // if we are out of bounds, return
      }      
      var startAtCol: u32 = computeInfo.startAtCols[linearIndex];
      var startAtRow: u32 = computeInfo.startAtRows[linearIndex];
      let endAtCol: u32 = min(computeInfo.endAtCols[linearIndex], ${listSize}u);
      let endAtRow: u32 = min(computeInfo.endAtRows[linearIndex], ${listSize}u);
      let is = &results.i[linearIndex];
      let js = &results.j[linearIndex];
      let distances = &results.distances[linearIndex];
      results.found[linearIndex] = 0; // initialize the found counter
      var found: u32 = 0;
      if (results.done[linearIndex] > 0) {
        return; // if we are done, return
      }
      for (var i = 0; i < ${maxIterationsPerThread}; i++) {
        if (startAtCol >= endAtCol && startAtRow >= endAtRow) {
          results.done[linearIndex] = 1;
          break;
        }
        if (found >= ${sparseResultSizePerThread}) {
          break;
        }
        let dist = combinedDistance(startAtCol, startAtRow);
        if (dist <= ${maxDistance}) {
          (*is)[found] = startAtCol;
          (*js)[found] = startAtRow;
          (*distances)[found] = dist;
          found = found + 1;
        }
        startAtCol = startAtCol + 1;
        if (startAtCol >= ${listSize}u) {
          startAtRow += 1;
          startAtCol = startAtRow + 1;
        }
      }
      results.found[linearIndex] = found;
      // update the startAtCols and startAtRows
      computeInfo.startAtCols[linearIndex] = startAtCol;
      computeInfo.startAtRows[linearIndex] = startAtRow;

    }

    // this will generate the distance script for each distance metric and then combine them into one
    ${getCombinedDistanceScript(distanceMetrics, processInfo.map((info) => info.maxEntryLen), maxDistance, aggregationFunction)}


    `});

  const pipeline = device.createComputePipeline({
    label: 'sparse matrix compute pipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'calcSparseMatrix',
    },
  });

  // generate startAtCols, startAtRows, endAtCols, endAtRows
  const startAtCols = new Uint32Array(numOfThreads);
  const startAtRows = new Uint32Array(numOfThreads);
  const endAtCols = new Uint32Array(numOfThreads);
  const endAtRows = new Uint32Array(numOfThreads);
  const chunkSize = Math.floor(condensedDistanceMatrixSize / numOfThreads); // size of the chunk per thread (in total)
  let startRow = 0;
  let startCol = 1;
  console.time('GPUthreadStarts');
  for (let i = 0; i < numOfThreads; i++) {
    const endIdx = i === numOfThreads - 1 ? condensedDistanceMatrixSize - 1 : (i + 1) * chunkSize;
    // fancy formulas to calculate the start and end indices for the condensed distance matrix for each thread start
    const endRow =
      listSize - 2 - Math.floor(Math.sqrt(-8 * endIdx + 4 * listSize * (listSize - 1) - 7) / 2 - 0.5);
    const endCol =
        endIdx - listSize * endRow + Math.floor((endRow + 1) * (endRow + 2) / 2);

    startAtCols[i] = startCol;
    startAtRows[i] = startRow;
    endAtCols[i] = endCol;
    endAtRows[i] = endRow;
    startRow = endRow;
    startCol = endCol;
    // const startRow = values[0].length - 2 - Math.floor(
    //   Math.sqrt(-8 * startIdx + 4 * values[0].length * (values[0].length - 1) - 7) / 2 - 0.5);
    // const startCol = startIdx - values[0].length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
  }
  console.timeEnd('GPUthreadStarts');

  // size of the computeInfo buffer
  const computeInfoBuffer32Size = numOfThreads * 4 + // startAtCols, startAtRows, endAtCols, endAtRows
    listSize * numOfColumns + // entrySizes
    numOfColumns + // weights
    processInfo.reduce((a, b) => a + b.sourceArraySize, 0);

  // size of the suppInfo buffer
  const suppInfoBuffer32Size = processInfo.reduce((a, b) => a + b.suppInfoSize, 0);

  // size of the results buffer
  const sparseMatrixEachArray32Size = sparseResultSizePerThread * numOfThreads;
  const resultsBuffer32Size = 3 * sparseMatrixEachArray32Size + numOfThreads + numOfThreads; // i, j, distances, found, done
  // create a buffer on the GPU to hold computeInfo
  // beware that struct must be padded to 16 bytes, so we need to calculate the size of the struct in 32bit values
  const computeInfoBufferSize = computeInfoBuffer32Size * Uint32Array.BYTES_PER_ELEMENT;
  let paddedComputeInfoBufferSize = computeInfoBufferSize;
  const remainder = computeInfoBufferSize & 15; // check if the size is a multiple of 16
  if (remainder !== 0)
    paddedComputeInfoBufferSize += 16 - remainder; // pad the size accordingly

  const computeInfoBuffer = device.createBuffer({
    label: 'compute info buffer',
    size: paddedComputeInfoBufferSize,
    usage:
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_SRC |
      GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const mappedComputeInfoArrayBuffer = computeInfoBuffer.getMappedRange(); // get full buffer

  // dynamic offset for the computeInfo buffer
  let computeInfoOffSet = 0;
  // first write the startAtCols, startAtRows, endAtCols, endAtRows
  const startAtColsBufferView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfThreads);
  startAtColsBufferView.set(startAtCols);
  computeInfoOffSet += numOfThreads * Uint32Array.BYTES_PER_ELEMENT; // array of 32bit values
  const startAtRowsBufferView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfThreads);
  startAtRowsBufferView.set(startAtRows);
  computeInfoOffSet += numOfThreads * Uint32Array.BYTES_PER_ELEMENT; // array of 32bit values
  const endAtColsBufferView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfThreads);
  endAtColsBufferView.set(endAtCols);
  computeInfoOffSet += numOfThreads * Uint32Array.BYTES_PER_ELEMENT; // array of 32bit values
  const endAtRowsBufferView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfThreads);
  endAtRowsBufferView.set(endAtRows);
  computeInfoOffSet += numOfThreads * Uint32Array.BYTES_PER_ELEMENT; // array of 32bit values

  // then write the entrySizes
  const entrySizesView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, arraySizes.length);
  entrySizesView.set(arraySizes);
  computeInfoOffSet += arraySizes.length * Uint32Array.BYTES_PER_ELEMENT; // array of 32bit values
  // then write the weights
  const weightsView = new Float32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfColumns);
  weightsView.set(weights);
  computeInfoOffSet += numOfColumns * Float32Array.BYTES_PER_ELEMENT;

  // finally, write the data itself
  for (const info of processInfo) {
    const ArrayConstructor = info.EncodedArrayConstructor;
    const chunkSize = info.sourceArraySize;
    const dataView = new ArrayConstructor(mappedComputeInfoArrayBuffer, computeInfoOffSet, chunkSize);//new ArrayConstructor(computeInfoBuffer.getMappedRange(computeInfoOffSet, chunkByteSize));
    dataView.set(info.flatSourceArray);
    computeInfoOffSet += chunkSize * ArrayConstructor.BYTES_PER_ELEMENT;
  }
  // we are done at this point.
  computeInfoBuffer.unmap();

  // create a buffer on the GPU to hold suppInfo
  // same here, we need to pad the size of the struct to 16 bytes
  const suppInfoBufferSize = suppInfoBuffer32Size * Uint32Array.BYTES_PER_ELEMENT;
  let paddedSuppInfoBufferSize = suppInfoBufferSize;
  const suppInfoRemainder = suppInfoBufferSize & 15; // check if the size is a multiple of 16
  if (suppInfoRemainder !== 0)
    paddedSuppInfoBufferSize += 16 - suppInfoRemainder; // pad the size accordingly

  paddedSuppInfoBufferSize = Math.max(paddedSuppInfoBufferSize, 16);
  const suppInfoBuffer = device.createBuffer({
    label: 'supp info buffer',
    size: paddedSuppInfoBufferSize,
    usage:
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_SRC |
      GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const mappedSuppInfoArrayBuffer = suppInfoBuffer.getMappedRange(); // get full buffer
  let suppInfoOffSet = 0;

  for (const info of processInfo) {
    if (info.suppInfoBuffer && info.suppInfoBuffer.byteLength > 0 && info.suppInfoSize > 0) {
      const ArrayConstructor = info.suppInfoType === WGPUENTRYTYPE.UINT32ARRAY ? Uint32Array : Float32Array;
      const suppInfoView = new ArrayConstructor(mappedSuppInfoArrayBuffer, suppInfoOffSet, info.suppInfoBuffer.length); //new ArrayConstructor(suppInfoBuffer.getMappedRange(suppInfoOffSet, info.suppInfoBuffer.byteLength));
      suppInfoView.set(info.suppInfoBuffer);
      suppInfoOffSet += info.suppInfoBuffer.byteLength; // info.suppInfoBuffer.length * ArrayConstructor.BYTES_PER_ELEMENT;
    }
  }

  if (suppInfoOffSet === 0) {
    const dummyView = new Uint32Array(mappedSuppInfoArrayBuffer, 0, 4);//new Uint32Array(suppInfoBuffer.getMappedRange(0, 16));
    dummyView.set([1, 1, 1, 1]);
  }
  suppInfoBuffer.unmap();

  // create a buffer for the results
  const resultsBufferSize = resultsBuffer32Size * Uint32Array.BYTES_PER_ELEMENT;
  let paddedResultsBufferSize = resultsBufferSize;
  const resultsRemainder = resultsBufferSize & 15; // check if the size is a multiple of 16
  if (resultsRemainder !== 0)
    paddedResultsBufferSize += 16 - resultsRemainder; // pad the size accordingly
  const resultsBuffer = device.createBuffer({
    label: 'results buffer',
    size: paddedResultsBufferSize,
    usage:
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_SRC
  });

  // Setup a bindGroup to tell the shader which
  // buffer to use for the computation
  const bindGroup = device.createBindGroup({
    label: 'bindGroup for sparse matrix buffer',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: computeInfoBuffer}},
      {binding: 1, resource: {buffer: suppInfoBuffer}},
      {binding: 2, resource: {buffer: resultsBuffer}},
    ],
  });

  //const pairComparisonsPerPass = maxIterationsPerThread * numOfThreads;
  //const passes = Math.ceil(condensedDistanceMatrixSize / pairComparisonsPerPass);
  // we will distpatch this many passes to the GPU, and it will handle indexes all by itself.
  // we already copied the start/end information to it, so it will know where to start and end on each pass.
  const resultsOutBuffer = device.createBuffer({
    label: 'results out buffer',
    size: resultsBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  const resultIs: Array<Uint32Array> = [];
  const resultJs: Array<Uint32Array> = [];
  const resultDistances: Array<Float32Array> = [];
  //let combinedFound = 0;
  let isAllDone = false;
  while (!isAllDone) {
    // Encode commands to do the computation
    const encoder = device.createCommandEncoder({
      label: 'distance encoder',
    });
    const pass = encoder.beginComputePass({
      label: 'distance compute pass',
    });
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(
      workgroupsDim,
      workgroupsDim
    );
    pass.end();

    encoder.copyBufferToBuffer(
      resultsBuffer,
      0,
      resultsOutBuffer,
      0,
      resultsOutBuffer.size
    );

    // Finish encoding and submit the commands
    const commandBuffer = encoder.finish();
    device.queue.submit([commandBuffer]);

    // Read the results
    await device.queue.onSubmittedWorkDone();

    await resultsOutBuffer.mapAsync(GPUMapMode.READ);

    const resultsOutArrayBuffer = resultsOutBuffer.getMappedRange();


    // read the results
    let resultOffset = 0;
    const resultsI = new Uint32Array(resultsOutArrayBuffer, resultOffset, sparseMatrixEachArray32Size);
    resultOffset += sparseMatrixEachArray32Size * Uint32Array.BYTES_PER_ELEMENT;
    const resultsJ = new Uint32Array(resultsOutArrayBuffer, resultOffset, sparseMatrixEachArray32Size);
    resultOffset += sparseMatrixEachArray32Size * Uint32Array.BYTES_PER_ELEMENT;
    const resultsDistances = new Float32Array(resultsOutArrayBuffer, resultOffset, sparseMatrixEachArray32Size);
    resultOffset += sparseMatrixEachArray32Size * Float32Array.BYTES_PER_ELEMENT;
    const resultsFound = new Uint32Array(resultsOutArrayBuffer, resultOffset, numOfThreads);
    resultOffset += numOfThreads * Uint32Array.BYTES_PER_ELEMENT;
    const resultsDone = new Uint32Array(resultsOutArrayBuffer, resultOffset, numOfThreads);
    isAllDone = resultsDone.every((d) => d === 1);

    const totalResults = resultsFound.reduce((a, b) => a + b, 0);

    const combinedI = new Uint32Array(totalResults);
    const combinedJ = new Uint32Array(totalResults);
    const combinedDistances = new Float32Array(totalResults);
    let combinedOffset = 0;
    for (let resI = 0; resI < resultsFound.length; resI++) {
      const found = resultsFound[resI];
      if (found === 0) continue;
      combinedI.set(resultsI.subarray(resI * sparseResultSizePerThread, resI * sparseResultSizePerThread + found), combinedOffset);
      combinedJ.set(resultsJ.subarray(resI * sparseResultSizePerThread, resI * sparseResultSizePerThread + found), combinedOffset);
      combinedDistances.set(resultsDistances.subarray(resI * sparseResultSizePerThread, resI * sparseResultSizePerThread + found), combinedOffset);
      combinedOffset += found;
    }
    resultIs.push(combinedI);
    resultJs.push(combinedJ);
    resultDistances.push(combinedDistances);

    resultsOutBuffer.unmap();
  }

  const totalSize = resultIs.reduce((a, b) => a + b.length, 0);
  const finalI = new Uint32Array(totalSize);
  const finalJ = new Uint32Array(totalSize);
  const finalDistances = new Float32Array(totalSize);
  let finalOffset = 0;
  for (let i = 0; i < resultIs.length; i++) {
    finalI.set(resultIs[i], finalOffset);
    finalJ.set(resultJs[i], finalOffset);
    finalDistances.set(resultDistances[i], finalOffset);
    finalOffset += resultIs[i].length;
  }

  // as rule mandates, destroy all buffers.
  computeInfoBuffer.destroy();
  suppInfoBuffer.destroy();
  resultsBuffer.destroy();
  resultsOutBuffer.destroy();

  return {i: finalI, j: finalJ, distance: finalDistances};
}


function getCombinedDistanceScript(distanceMetrics: WEBGPUDISTANCE[], maxEntryLens: number[], maxDistance: number, aggregation: WEBGSLAGGREGATION) {
  const distanceWgsls = distanceMetrics.map((metric, i) => {
    return `
        fn distanceScript${i}(aIndex: u32, bIndex: u32) -> f32 {
          let a = computeInfo.data${i}[aIndex];
          let b = computeInfo.data${i}[bIndex];
          let maxDistance: f32 = ${maxDistance};
          ${webGPUFunctions[metric](maxEntryLens[i], i)}
        }
      `;
  });

  const allDistanceScripts = distanceWgsls.join('\n');

  const combineDistancesScript = `
      fn combinedDistance(aIndex: u32, bIndex: u32) -> f32 {
        var distances: array<f32, ${distanceMetrics.length}>;
        ${distanceMetrics.map((_, i) => `distances[${i}] = distanceScript${i}(aIndex, bIndex);`).join('\n')}
        ${WEBGSLAGGREGATIONFUNCTIONS[aggregation](distanceMetrics.length)}
      }
    
    `;

  return allDistanceScripts + '\n' + combineDistancesScript;
}
