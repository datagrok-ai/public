/* eslint-disable max-len */
import {getGPUDevice} from '../getGPUDevice';
import {WEBGSLAGGREGATION, WEBGSLAGGREGATIONFUNCTIONS} from '../multi-col-distances/webGPU-aggregation';
import {SupportedEntryTypes, WEBGPUDISTANCE, webGPUFunctions,
  WGPUENTRYTYPE} from '../multi-col-distances/webGPU-multicol-distances';
import {webGPUProcessInfo} from '../preprocessing/webGPU-process-info';


/** generate KNN based on list of lists of entries.
 *  these entries are each encoded as Uint32Array or FLOAT32Array (depending on their type).
 * for example, sequences would be encoded as Uint32Array based on char code of the letter at each position.
 * [65, 66, 67, 68, 69] would be a sequence of 5 letters.
 * for chemical fingerprints, it would be a binary array of 0s and 1s,
 * represented as Uint32Array(_data property of DG bitarray).
 *
 * Be ware that size of entryList, distanceMetrics, weights and options must be the same.
 * if there are no options for entries i, pass an empty object.
 * for now options are needed for:
 * needleman-wunsch and monomer chemical distances: see {@link BioDistanceFnOptions} as for how it should be passed
 * numeric distances (Difference): {range: number} where range is the range of the values in the column (max - min).
 * in both cases, if options are not provided, they will be calculated automatically.
 */
export async function multiColWebGPUKNN(
  entryList: Array<Array<SupportedEntryTypes>>, // list of lists of entries, for multiple columns
  knnSize: number = 15, // size of the k-nearest neighbors
  distanceMetrics: WEBGPUDISTANCE[], // distance metrics for each column
  aggregationFunction: WEBGSLAGGREGATION, // aggregation function for the distances
  weights: number[], // weights for each column
  options: {[key: string]: any}[] // supplementary options for each column
) {
  // first, check that all the supplementary options are provided and are the same length:
  if (options.length !== entryList.length ||
    options.length !== distanceMetrics.length || options.length !== weights.length)
    throw new Error('Options, weigths and distance functions must be provided for each column');

  // check that all the entry lists are the same length
  if (entryList.some((list) => list.length !== entryList[0].length))
    throw new Error('All entry lists must be the same length');


  const availableDistanceMetrics = Object.values(WEBGPUDISTANCE);
  if (distanceMetrics.some((metric) => !availableDistanceMetrics.includes(metric)))
    throw new Error('Invalid distance metrics provided: ' + distanceMetrics.join(', '));

  const availableAggregationFunctions = Object.values(WEBGSLAGGREGATION);
  if (!availableAggregationFunctions.includes(aggregationFunction))
    throw new Error('Invalid aggregation function provided: ' + aggregationFunction);

  const numOfColumns = entryList.length; // number of columns
  if (numOfColumns === 0)
    throw new Error('No columns provided. Please provide at least one column of data.');

  const device = await getGPUDevice();
  if (!device)
    return null;

  // device may be lost
  let deviceLost = false;
  device.lost.then(() => {
    deviceLost = true;
  });

  const listSize = entryList[0].length; // size of each list (or column)
  const processInfo = entryList.map((entry, i) => {
    return webGPUProcessInfo(entry, distanceMetrics[i], i, options[i]);
  });


  if (numOfColumns === 1)
    aggregationFunction = WEBGSLAGGREGATION.MANHATTAN; // save a bit of time


  // combine all struct types into one to put into the suppInfo struct.
  let suppInfoWgsl = processInfo.map((info) => info.suppInfoStructWgsl)
    .filter((wgsl) => !!wgsl && wgsl != '').join(',\n');
  // structures in wgsl must have at least one member, so if we have no structures, we need to add a dummy one
  let needsDummy = false;
  if (!suppInfoWgsl || suppInfoWgsl.trim() == '') {
    needsDummy = true;
    suppInfoWgsl = '\ndummy: f32\n';
  }


  // combine all complexities into one
  const combinedComplexity = processInfo.reduce((a, b) => a + b.complexity, 0);
  // combine all data wgsl struct code into one
  const dataWgsl = processInfo.map((info) => info.dataStructWgsl).filter((wgsl) => !!wgsl && wgsl != '').join(',\n');
  // combine all array sizes into one array (easier for setting)
  const arraySizes = new Uint32Array(numOfColumns * listSize);
  processInfo.forEach((info, i) => {
    arraySizes.set(info.arraySizes, i * listSize);
  }); // array.flat is not as optimized as this

  // if we try to map large knn directly from GPU, sometimes, device disconnects. so we need to do it in chunks, a good number
  // we found is 10000. So we will perform computations in chunks of 10000.
  const computationsPerPass = 10000;

  // also, as the computation per thread takes some time, we also need to divide the work into smaller chunks. meaning
  // that we will divide nummper of pair computations per pass. start and end of the pair comparisons will be stored in startAtEndAt buffer (vec4<u32>) as z and w coordinates.
  const pairComparisonsPerPass = Math.ceil(10000 / combinedComplexity);


  const workGroupDivision = 10; // how many threads inside of one workgroup dimension (in this case 10 * 10 threads per workgroup)

  const threadsPerWorkgroup = workGroupDivision * workGroupDivision;

  const workgroupsDim = Math.ceil(
    Math.sqrt(Math.ceil(computationsPerPass / threadsPerWorkgroup))
  ); // how many workgroups per 2d dimension

  const globalThreadDimSize = workgroupsDim * workGroupDivision; // how many threads per 2d dimension

  // console.log(getCombinedDistanceScript(distanceMetrics, processInfo.map((info) => info.maxEntryLen), knnSize, aggregationFunction));
  // return;
  const resultIndexes = new Array(listSize)
    .fill(null)
    .map(() => new Uint32Array(knnSize));
  const resultDistances = new Array(listSize)
    .fill(null)
    .map(() => new Float32Array(knnSize));

  const module = device.createShaderModule({
    label: 'KNN compute shader',
    code: `
          // this struct will contain all the info about the computation, startAtEndAt will contain the start and end of the knnDistances and knnIndexes.
          // array of sizes for each entries, and the main data as arrays of arrays called data0, data1 and so on. good thing is that because the first entry is vec4<u32>,
          // there will be no paddings, so no need to worry about padding data. also, arrays and matrices get stucked together, so no padding there as well.
          // what we need to worry about is the padding of overall struct. the size of overall struct will be factor of 16 bytes, so keep that in mind.
          struct ComputeInfo {
            // the x coordinate will contain the start index of the knnDistances and knnIndexes, while y will contain the end index   
            // the z coordinate will contain the start of the pair comparisons, while w will contain the end of the pair comparisons
            // just keep in mind that this vec4 will be in first 16 bytes of corresponding buffer.
            startAtEndAt: vec4<u32>,
            // the ACTUALLY sizes of each entry
            entrySizes: array<array<u32, ${listSize}>, ${numOfColumns}>,
            // the weights for each entry
            weights: array<f32, ${numOfColumns}>,
            // the data for each entry
            ${dataWgsl} // an example of the dataWgsl would be:
            //data0: array<array<u32,20>,100>,
            //data1: array<array<u32,20>,100>
          };

          struct SuppInfo {
            // struct containing all the supplementary info, like scoring matrix, alphabet indexes, range, etc.
            ${suppInfoWgsl}
          };
    
          @group(0) @binding(0) var<storage, read_write> knnIndexes: array<array<u32, ${knnSize}>, ${computationsPerPass}>;
          @group(0) @binding(1) var<storage, read_write> knnDistances: array<array<f32, ${knnSize}>, ${computationsPerPass}>; // each time just compute for a subset of the list
          @group(0) @binding(2) var<storage, read_write> computeInfo: ComputeInfo;
          @group(0) @binding(3) var<storage, read_write> suppInfo: SuppInfo;
          
          @compute @workgroup_size(${workGroupDivision}, ${workGroupDivision}) fn calcKNN(
            @builtin(global_invocation_id) id: vec3<u32>
          ) {
            ${needsDummy ? `let otherDummy = suppInfo.dummy * 2;` : ''} // just to make sure that the suppInfo is not optimized out
            let col = id.x; //* ${workGroupDivision} + localId.x;
            let row = id.y; //* ${workGroupDivision} + localId.y;
            let graphIndex = row * ${globalThreadDimSize} + col;
            let index = graphIndex + computeInfo.startAtEndAt.x; // add the starting index of the knnDistances and knnIndexes
    
            if (index >= min(${listSize}u, computeInfo.startAtEndAt.y)) {return;}
    
            let pairComparisonStartAt = computeInfo.startAtEndAt.z;
            let pairComparisonEndAt = min(computeInfo.startAtEndAt.w, ${listSize}u);
    
    
            // only clear the knnDistances and knnIndexes if we are at the start of the pair comparison
            if (pairComparisonStartAt == 0u) {
              for (var i = 0u; i < ${knnSize}; i = i + 1u) {
                knnDistances[graphIndex][i] = 99999.0;
                knnIndexes[graphIndex][i] = 0u;
              }
            }
    
            for (var i: u32 = pairComparisonStartAt; i < pairComparisonEndAt; i = i + 1u) {
              if (i == index) {continue;}
              let dist = combinedDistance(index, i);
              insertKnn(graphIndex, dist, i);
            }
            
          }
          // this will generate the distance script for each distance metric and then combine them into one
          ${getCombinedDistanceScript(distanceMetrics, processInfo.map((info) => info.maxEntryLen), knnSize, aggregationFunction)}
    
          fn insertKnn(knnIndex: u32, dist: f32, index: u32) {
            // small optimization, if the distance is larger than the last element in the knnDistances, we can skip
            if (dist >= knnDistances[knnIndex][${knnSize} - 1u]) {return;}
            for (var i = 0u; i < ${knnSize}; i = i + 1u) {
              if (dist < knnDistances[knnIndex][i]) {
                for (var j = ${knnSize} - 1u; j > i; j = j - 1u) {
                  knnDistances[knnIndex][j] = knnDistances[knnIndex][j - 1u];
                  knnIndexes[knnIndex][j] = knnIndexes[knnIndex][j - 1u];
                }
                knnDistances[knnIndex][i] = dist;
                knnIndexes[knnIndex][i] = index;
                return;
              }
            }
          }
        `,
  });

  const pipeline = device.createComputePipeline({
    label: 'hamming compute pipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'calcKNN',
    },
  });

  // calculate the size of computeInfo struct buffer in terms of 32bit values
  //           startAtEndAt=4,          entrySizes,        weights,          data0 + data1 + ...
  const computeInfo32Size = 4 + numOfColumns * listSize + numOfColumns + processInfo.reduce((a, b) => a + b.sourceArraySize, 0);

  // calculate the size of suppInfo struct buffer in terms of 32bit values
  const suppInfo32Size = processInfo.reduce((a, b) => a + b.suppInfoSize, 0);


  // create a buffer on the GPU to hold computeInfo
  // beware that struct must be padded to 16 bytes, so we need to calculate the size of the struct in 32bit values
  const computeInfoBufferSize = computeInfo32Size * Uint32Array.BYTES_PER_ELEMENT;
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
  // create vars to hold startAtEndAt
  let startAt = 0;
  let endAt = computationsPerPass;
  let pairComparisonStartAt = 0;
  let pairComparisonEndAt = pairComparisonsPerPass;
  const startAtEndAtSize = 4;
  // copy the data into the buffer at correct places
  let computeInfoOffSet = 0;
  const startAtEndAtView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, startAtEndAtSize); //new Uint32Array(computeInfoBuffer.getMappedRange(computeInfoOffSet, 4 * Uint32Array.BYTES_PER_ELEMENT));
  startAtEndAtView.set([startAt, endAt, pairComparisonStartAt, pairComparisonEndAt]);
  // write entry sizes
  computeInfoOffSet += startAtEndAtSize * Uint32Array.BYTES_PER_ELEMENT;
  const entrySizesView = new Uint32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, arraySizes.length); //new Uint32Array(computeInfoBuffer.getMappedRange(computeInfoOffSet, arraySizes.byteLength));
  entrySizesView.set(arraySizes);
  // write weights
  computeInfoOffSet += arraySizes.byteLength; // arraySizes.length * Uint32Array.BYTES_PER_ELEMENT;
  const weightsView = new Float32Array(mappedComputeInfoArrayBuffer, computeInfoOffSet, numOfColumns); //new Float32Array(computeInfoBuffer.getMappedRange(computeInfoOffSet, numOfColumns * Float32Array.BYTES_PER_ELEMENT));
  weightsView.set(weights);
  // write data
  computeInfoOffSet += numOfColumns * Float32Array.BYTES_PER_ELEMENT;
  for (const info of processInfo) {
    //device.queue.writeBuffer(computeInfoBuffer, computeInfoOffSet, info.flatSourceArray, 0, chunkByteSize);
    const ArrayConstructor = info.EncodedArrayConstructor;
    const chunkSize = info.sourceArraySize;
    //@ts-ignore reason: new typescript shinanigens
    const dataView = new ArrayConstructor(mappedComputeInfoArrayBuffer, computeInfoOffSet, chunkSize);//new ArrayConstructor(computeInfoBuffer.getMappedRange(computeInfoOffSet, chunkByteSize));
    dataView.set(info.flatSourceArray);
    computeInfoOffSet += chunkSize * ArrayConstructor.BYTES_PER_ELEMENT;
  }
  // we are done at this point.
  computeInfoBuffer.unmap();

  // create a buffer on the GPU to hold suppInfo
  // same here, we need to pad the size of the struct to 16 bytes
  const suppInfoBufferSize = suppInfo32Size * Uint32Array.BYTES_PER_ELEMENT;
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
      //@ts-ignore reason: new typescript shinanigens
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

  const outKnnArraySize =
  knnSize * Uint32Array.BYTES_PER_ELEMENT * computationsPerPass; // size of the slice of the knnIndexes and knnDistances

  // create a buffer on the GPU to hold knnIndexes and knnDistances
  const bufferDistances = device.createBuffer({
    label: 'buffer distances',
    size: outKnnArraySize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
  });

  const bufferIndexes = device.createBuffer({
    label: 'buffer indexes',
    size: outKnnArraySize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
  });
  // no need to copy anything to the knn buffers.

  const bytesPerOutArrayRow = knnSize * Uint32Array.BYTES_PER_ELEMENT; // how many bytes one row of the knnDistances and knnIndexes will take


  // Setup a bindGroup to tell the shader which
  // buffer to use for the computation
  const bindGroup = device.createBindGroup({
    label: 'bindGroup for knn buffer',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: bufferIndexes}},
      {binding: 1, resource: {buffer: bufferDistances}},
      {binding: 2, resource: {buffer: computeInfoBuffer}},
      {binding: 3, resource: {buffer: suppInfoBuffer}},
      // { binding: 4, resource: { buffer: arraySizesBuffer } },
    ],
  });

  //create a buffer on the GPU to get a copy of the results
  const resultBufferDistances = device.createBuffer({
    label: 'result buffer distances',
    size: bufferDistances.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  const resultBufferIndexes = device.createBuffer({
    label: 'result buffer indexes',
    size: bufferIndexes.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  for (let iter = 0; iter < Math.ceil(listSize / computationsPerPass); iter++) {
    startAt = iter * computationsPerPass;
    endAt = Math.min((iter + 1) * computationsPerPass, listSize);
    const pairComparisonPasses = Math.ceil(listSize / pairComparisonsPerPass);
    for (let pairIter = 0; pairIter < pairComparisonPasses; pairIter++) {
      pairComparisonStartAt = pairIter * pairComparisonsPerPass;
      pairComparisonEndAt = Math.min(
        (pairIter + 1) * pairComparisonsPerPass,
        listSize
      );
      // write only four values to the buffer, startAt, endAt, pairComparisonStartAt, pairComparisonEndAt, corresponding to vec4<u32> at the start of info struct
      device.queue.writeBuffer(
        computeInfoBuffer,
        0,
        new Uint32Array([
          startAt,
          endAt,
          pairComparisonStartAt,
          pairComparisonEndAt,
        ])
      );

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
        //workgroupsDim,
        Math.ceil(computationsPerPass / workgroupsDim / threadsPerWorkgroup)
      );
      pass.end();

      // if we are at the last pair comparison pass, we need to read the results
      if (pairIter === pairComparisonPasses - 1) {
        // Encode a command to copy the results to a mappable buffer.
        encoder.copyBufferToBuffer(
          bufferDistances,
          0,
          resultBufferDistances,
          0,
          resultBufferDistances.size
        );

        encoder.copyBufferToBuffer(
          bufferIndexes,
          0,
          resultBufferIndexes,
          0,
          resultBufferIndexes.size
        );

        // Finish encoding and submit the commands
        const commandBuffer = encoder.finish();
        device.queue.submit([commandBuffer]);

        // Read the results
        //console.time('pass end');
        await device.queue.onSubmittedWorkDone();
        // console.timeEnd('pass end');

        // console.time("read");
        //const offsetBytes = startAt * bytesPerOutArrayRow;
        //const sizeBytes = (endAt - startAt) * bytesPerOutArrayRow;
        await resultBufferDistances.mapAsync(GPUMapMode.READ);
        await resultBufferIndexes.mapAsync(GPUMapMode.READ);
        // console.timeEnd("read");
        const indexes = resultBufferIndexes.getMappedRange();
        const distances = resultBufferDistances.getMappedRange();

        // console.time("decode");

        for (let i = 0; i < endAt - startAt; i++) {
          const index = new Uint32Array(
            indexes,
            i * bytesPerOutArrayRow,
            knnSize
          );
          const distance = new Float32Array(
            distances,
            i * bytesPerOutArrayRow,
            knnSize
          );
          resultIndexes[startAt + i].set(index);
          resultDistances[startAt + i].set(distance);
        }

        resultBufferIndexes.unmap();
        resultBufferDistances.unmap();
      } else {
        const commandBuffer = encoder.finish();
        device.queue.submit([commandBuffer]);

        // Read the results
        // console.time('pass between');
        await device.queue.onSubmittedWorkDone();
        // console.timeEnd('pass between');
      }
      // device may get lost during opperation, so in this case, we need to return null
      if (deviceLost)
        return null;
    }
  }

  bufferDistances.destroy();
  bufferIndexes.destroy();
  computeInfoBuffer.destroy();
  suppInfoBuffer.destroy();
  resultBufferDistances.destroy();
  resultBufferIndexes.destroy();

  if (deviceLost)
    return null;
  return {knnIndexes: resultIndexes, knnDistances: resultDistances};
}

function getCombinedDistanceScript(distanceMetrics: WEBGPUDISTANCE[], maxEntryLens: number[], knnSize: number, aggregation: WEBGSLAGGREGATION) {
  const distanceWgsls = distanceMetrics.map((metric, i) => {
    return `
        fn distanceScript${i}(aIndex: u32, bIndex: u32) -> f32 {
          let a = computeInfo.data${i}[aIndex];
          let b = computeInfo.data${i}[bIndex];
          let maxDistance: f32 = knnDistances[aIndex - computeInfo.startAtEndAt.x][${knnSize} - 1u];
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
