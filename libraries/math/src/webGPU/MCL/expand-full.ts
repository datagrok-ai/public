/* eslint-disable max-len */

// this file is currently not used but contains two very good and efficient implementations of sparse matrix FULL multiplication

import {toOffsetForm} from '../umap/utils';
import {MCLOpReturnType} from './types';
import {knnToSparseForm} from './utils';


export async function expand2(
  device: GPUDevice, knnSimilarities: Float32Array, knnIndexes: Uint32Array, offsets: Uint32Array, nRows: number
): Promise<MCLOpReturnType> {
  const blockSizeDim = 4000;
  const neededThreads = blockSizeDim * blockSizeDim;
  const workGroupThreadsPerDim = 10;
  const totalWorkgroupThreads = workGroupThreadsPerDim * workGroupThreadsPerDim;
  const neededWorkGroups = Math.ceil(neededThreads / totalWorkgroupThreads);
  const workGroupDim = Math.ceil(Math.sqrt(neededWorkGroups));

  const outKNNIndexes = new Array(nRows).fill(0).map(() => [] as number[]);
  const outKNNSimilarities = new Array(nRows).fill(0).map(() => [] as number[]);
  const module = device.createShaderModule({
    label: 'expand',
    code: `
          struct SparseKNN {
              knnSimilarities: array<f32, ${knnSimilarities.length}>,
              knnIndexes: array<u32, ${knnIndexes.length}>,
              offsets: array<u32, ${offsets.length}>,
          }
  
          @group(0) @binding(0) var<storage, read_write> sparseKNN: SparseKNN;
          @group(0) @binding(1) var<storage, read> startEnd: vec4<u32>; // xy will contain top left corner of the block, zw will contain bottom right corner of the block
          @group(0) @binding(2) var<storage, read_write> resultBlock: array<array<f32, ${blockSizeDim}>, ${blockSizeDim}>;
          @compute @workgroup_size(${workGroupThreadsPerDim}, ${workGroupThreadsPerDim}) fn expand(
              @builtin(global_invocation_id) id: vec3<u32>,
            ) {
              let col = id.x + startEnd.x;
              let row = id.y + startEnd.y;
              if (row >= min(${nRows}, startEnd.w) || col >= min(${nRows},startEnd.z) || id.x >= ${blockSizeDim} || id.y >= ${blockSizeDim}) {
                  return;
              }
              // we take the row subbarray and col subbarray and multiply them elementwise, if we find common indexes
              // keep in mind that in this representation of the sparse knn matrix, values are duplicated along the diagonals
              let offsetBeginRow = sparseKNN.offsets[row];
              let offsetEndRow = sparseKNN.offsets[row + 1];
              let offsetBeginCol = sparseKNN.offsets[col];
              let offsetEndCol = sparseKNN.offsets[col + 1];
              var sum = 0.0;
              if (row == col) { // in case of diagonal values, we just square them, saves some time
                  for (var i = offsetBeginRow; i < offsetEndRow; i = i + 1) {
                      sum = sum + sparseKNN.knnSimilarities[i] * sparseKNN.knnSimilarities[i];
                  }
                  resultBlock[id.y][id.x] = sum;
                  return;
              }
              for (var i = offsetBeginRow; i < offsetEndRow; i = i + 1) {
                  for(var j = offsetBeginCol; j < offsetEndCol; j = j + 1) {
                      if (sparseKNN.knnIndexes[i] == sparseKNN.knnIndexes[j]) {
                          sum = sum + sparseKNN.knnSimilarities[i] * sparseKNN.knnSimilarities[j];
                          break;
                      }
                  }
              }
              resultBlock[id.y][id.x] = sum;
          }
      `});

  const pipeline = device.createComputePipeline({
    label: 'expand compute pipeline',
    layout: 'auto',
    compute: {
      module: module,
      entryPoint: 'expand',
    },
  });

  const sparseKNNBuffer32Size = knnSimilarities.length + knnIndexes.length + offsets.length;
  let sparseKNNBufferByteSize = sparseKNNBuffer32Size * 4;
  const remainder = sparseKNNBufferByteSize & 15;
  if (remainder !== 0)
    sparseKNNBufferByteSize += 16 - remainder;
  const sparseKNNBuffer = device.createBuffer({
    label: 'sparse knn buffer',
    size: sparseKNNBufferByteSize,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const sparseKNNArrayBuffer = sparseKNNBuffer.getMappedRange();
  // set similarities
  new Float32Array(sparseKNNArrayBuffer, 0, knnSimilarities.length).set(knnSimilarities);
  // set indexes
  new Uint32Array(sparseKNNArrayBuffer, knnSimilarities.length * 4, knnIndexes.length).set(knnIndexes);
  // set offsets
  new Uint32Array(sparseKNNArrayBuffer, (knnSimilarities.length + knnIndexes.length) * 4, offsets.length).set(offsets);
  sparseKNNBuffer.unmap();

  const startEndBuffer = device.createBuffer({
    label: 'start end buffer',
    size: 16,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  new Uint32Array(startEndBuffer.getMappedRange()).set([0, 0, blockSizeDim, blockSizeDim]);
  startEndBuffer.unmap();

  const resultBlockBuffer = device.createBuffer({
    label: 'result block buffer',
    size: blockSizeDim * blockSizeDim * 4,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
  });

  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: sparseKNNBuffer}},
      {binding: 1, resource: {buffer: startEndBuffer}},
      {binding: 2, resource: {buffer: resultBlockBuffer}},
    ],
  });

  const order = Math.floor(Math.max(Math.log10(nRows), 2)) + 2;
  // minimum value after expansion.
  const pruneValue = Math.pow(10, -order);

  const outBlockBuffer = device.createBuffer({
    label: 'out block buffer',
    size: blockSizeDim * blockSizeDim * 4,
    usage:
          GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
  });


  for (let i = 0; i < Math.ceil(nRows / blockSizeDim); i++) {
    const startY = i * blockSizeDim;
    const endY = Math.min(nRows, (i + 1) * blockSizeDim);
    for (let j = 0; j < Math.ceil(nRows / blockSizeDim); j++) {
      const startX = j * blockSizeDim;
      const endX = Math.min(nRows, (j + 1) * blockSizeDim);
      device.queue.writeBuffer(startEndBuffer, 0, new Uint32Array([startX, startY, endX, endY]));
      const encoder = device.createCommandEncoder({
        label: 'expand encoder',
      });
      const pass = encoder.beginComputePass({
        label: 'expand compute pass',
      });
      pass.setPipeline(pipeline);
      pass.setBindGroup(0, bindGroup);
      pass.dispatchWorkgroups(
        workGroupDim,
        workGroupDim
      );
      pass.end();

      encoder.copyBufferToBuffer(resultBlockBuffer, 0, outBlockBuffer, 0, blockSizeDim * blockSizeDim * 4);
      device.queue.submit([encoder.finish()]);
      //console.time('onSubmittedWorkDone');
      await device.queue.onSubmittedWorkDone();
      //console.timeEnd('onSubmittedWorkDone');
      await outBlockBuffer.mapAsync(GPUMapMode.READ);
      const outBlock = new Float32Array(outBlockBuffer.getMappedRange());
      //console.time('fillout');
      for (let row = 0; row < blockSizeDim; row++) {
        const rowIdx = row + startY;
        if (rowIdx >= Math.min(nRows, endY))
          break;
        const rowMultiplier = row * blockSizeDim;
        for (let col = 0; col < blockSizeDim; col++) {
          const colIdx = col + startX;
          if (colIdx >= Math.min(nRows, endX))
            break;
          const value = outBlock[rowMultiplier + col];
          if (value > pruneValue) {
            outKNNIndexes[colIdx].push(rowIdx);
            outKNNSimilarities[colIdx].push(value);
          }
        }
      }
      //console.timeEnd('fillout');
      outBlockBuffer.unmap();
    }
  }

  // destroy
  sparseKNNBuffer.destroy();
  startEndBuffer.destroy();
  resultBlockBuffer.destroy();
  outBlockBuffer.destroy();
  return knnToSparseForm(outKNNIndexes, outKNNSimilarities);
}


export async function expand(
  device: GPUDevice, knnSimilarities: Float32Array, knnIndexes: Uint32Array, offsets: Uint32Array, nRows: number
): Promise<MCLOpReturnType> {
  const blockSizeDim = 4000;
  const neededThreads = blockSizeDim * blockSizeDim;
  const workGroupThreadsPerDim = 10;
  const totalWorkgroupThreads = workGroupThreadsPerDim * workGroupThreadsPerDim;
  const neededWorkGroups = Math.ceil(neededThreads / totalWorkgroupThreads);
  const workGroupDim = Math.ceil(Math.sqrt(neededWorkGroups));

  const outSparseSimilarities: Float32Array[] = [];
  const outSparseIindexes: Uint32Array[] = [];
  const outSparseJindexes: Uint32Array[] = [];
  const order = Math.floor(Math.max(Math.log10(nRows), 2));
  // minimum value after expansion.
  const pruneValue = Math.pow(10, -order);
  const module = device.createShaderModule({
    label: 'expand',
    code: `
          struct SparseKNN {
              knnSimilarities: array<f32, ${knnSimilarities.length}>,
              knnIndexes: array<u32, ${knnIndexes.length}>,
              offsets: array<u32, ${offsets.length}>,
          }
  
          struct SparseResult {
            is: array<u32, ${neededThreads}>,
            js: array<u32, ${neededThreads}>,
            similarities: array<f32, ${neededThreads}>,
          }
  
          @group(0) @binding(0) var<storage, read_write> sparseKNN: SparseKNN;
          @group(0) @binding(1) var<storage, read> startEnd: vec4<u32>; // xy will contain top left corner of the block, zw will contain bottom right corner of the block
          @group(0) @binding(2) var<storage, read_write> sparseResult: SparseResult;
          @group(0) @binding(3) var<storage, read_write> foundIndex: atomic<u32>; // this will denote the index of the next element to be written
          @compute @workgroup_size(${workGroupThreadsPerDim}, ${workGroupThreadsPerDim}) fn expand(
              @builtin(global_invocation_id) id: vec3<u32>,
            ) {
              let col = id.x + startEnd.x;
              let row = id.y + startEnd.y;
              if (row >= min(${nRows}, startEnd.w) || col >= min(${nRows},startEnd.z) || id.x >= ${blockSizeDim} || id.y >= ${blockSizeDim}) {
                  return;
              }
              // we take the row subbarray and col subbarray and multiply them elementwise, if we find common indexes
              // keep in mind that in this representation of the sparse knn matrix, values are duplicated along the diagonals
              let offsetBeginRow = sparseKNN.offsets[row];
              let offsetEndRow = sparseKNN.offsets[row + 1];
              let offsetBeginCol = sparseKNN.offsets[col];
              let offsetEndCol = sparseKNN.offsets[col + 1];
              var sum = 0.0;
              if (row == col) { // in case of diagonal values, we just square them, saves some time
                  for (var i = offsetBeginRow; i < offsetEndRow; i = i + 1) {
                      sum = sum + sparseKNN.knnSimilarities[i] * sparseKNN.knnSimilarities[i];
                  }
                  if (sum > ${pruneValue}) {
                    let currentSparseIndex = atomicAdd(&foundIndex, 1); // atomically increment found index
                    sparseResult.is[currentSparseIndex] = row;
                    sparseResult.js[currentSparseIndex] = col;
                    sparseResult.similarities[currentSparseIndex] = sum;
                  }
                  return;
              }
              for (var i = offsetBeginRow; i < offsetEndRow; i = i + 1) {
                  for(var j = offsetBeginCol; j < offsetEndCol; j = j + 1) {
                      if (sparseKNN.knnIndexes[i] == sparseKNN.knnIndexes[j]) {
                          sum = sum + sparseKNN.knnSimilarities[i] * sparseKNN.knnSimilarities[j];
                          break;
                      }
                  }
              }
              if (sum > ${pruneValue}) {
                let currentSparseIndex = atomicAdd(&foundIndex, 1); // atomically increment found index
                sparseResult.is[currentSparseIndex] = row;
                sparseResult.js[currentSparseIndex] = col;
                sparseResult.similarities[currentSparseIndex] = sum;
              }
          }
      `});

  const pipeline = device.createComputePipeline({
    label: 'expand compute pipeline',
    layout: 'auto',
    compute: {
      module: module,
      entryPoint: 'expand',
    },
  });

  const sparseKNNBuffer32Size = knnSimilarities.length + knnIndexes.length + offsets.length;
  let sparseKNNBufferByteSize = sparseKNNBuffer32Size * 4;
  const remainder = sparseKNNBufferByteSize & 15;
  if (remainder !== 0)
    sparseKNNBufferByteSize += 16 - remainder;
  const sparseKNNBuffer = device.createBuffer({
    label: 'sparse knn buffer',
    size: sparseKNNBufferByteSize,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const sparseKNNArrayBuffer = sparseKNNBuffer.getMappedRange();
  // set similarities
  new Float32Array(sparseKNNArrayBuffer, 0, knnSimilarities.length).set(knnSimilarities);
  // set indexes
  new Uint32Array(sparseKNNArrayBuffer, knnSimilarities.length * 4, knnIndexes.length).set(knnIndexes);
  // set offsets
  new Uint32Array(sparseKNNArrayBuffer, (knnSimilarities.length + knnIndexes.length) * 4, offsets.length).set(offsets);
  sparseKNNBuffer.unmap();

  const startEndBuffer = device.createBuffer({
    label: 'start end buffer',
    size: 16,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  new Uint32Array(startEndBuffer.getMappedRange()).set([0, 0, blockSizeDim, blockSizeDim]);
  startEndBuffer.unmap();

  const resultBlock32Size = neededThreads * 3;
  let resultBlockByteSize = resultBlock32Size * 4;
  const resultBlockRemainder = resultBlockByteSize & 15;
  if (resultBlockRemainder !== 0)
    resultBlockByteSize += 16 - resultBlockRemainder;
  const resultBlockBuffer = device.createBuffer({
    label: 'result block buffer',
    size: resultBlockByteSize,
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
  });

  const fouundResultsBuffer = device.createBuffer({
    label: 'found results buffer',
    size: 4, // just a single u32
    usage:
              GPUBufferUsage.STORAGE |
              GPUBufferUsage.COPY_SRC |
              GPUBufferUsage.COPY_DST,
  });

  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: sparseKNNBuffer}},
      {binding: 1, resource: {buffer: startEndBuffer}},
      {binding: 2, resource: {buffer: resultBlockBuffer}},
      {binding: 3, resource: {buffer: fouundResultsBuffer}},
    ],
  });

  const outBlockBuffer = device.createBuffer({
    label: 'out block buffer',
    size: resultBlockBuffer.size,
    usage:
          GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
  });

  const outFoundResultsBuffer = device.createBuffer({
    label: 'out found results buffer',
    size: fouundResultsBuffer.size,
    usage:
          GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
  });


  for (let i = 0; i < Math.ceil(nRows / blockSizeDim); i++) {
    const startY = i * blockSizeDim;
    const endY = Math.min(nRows, (i + 1) * blockSizeDim);
    for (let j = 0; j < Math.ceil(nRows / blockSizeDim); j++) {
      const startX = j * blockSizeDim;
      const endX = Math.min(nRows, (j + 1) * blockSizeDim);
      device.queue.writeBuffer(startEndBuffer, 0, new Uint32Array([startX, startY, endX, endY]));
      device.queue.writeBuffer(fouundResultsBuffer, 0, new Uint32Array([0])); // initialize found results to 0
      const encoder = device.createCommandEncoder({
        label: 'expand encoder',
      });
      const pass = encoder.beginComputePass({
        label: 'expand compute pass',
      });
      pass.setPipeline(pipeline);
      pass.setBindGroup(0, bindGroup);
      pass.dispatchWorkgroups(
        workGroupDim,
        workGroupDim
      );
      pass.end();

      encoder.copyBufferToBuffer(resultBlockBuffer, 0, outBlockBuffer, 0, outBlockBuffer.size);
      encoder.copyBufferToBuffer(fouundResultsBuffer, 0, outFoundResultsBuffer, 0, outFoundResultsBuffer.size);
      device.queue.submit([encoder.finish()]);
      //console.time('onSubmittedWorkDone');
      await device.queue.onSubmittedWorkDone();
      //console.timeEnd('onSubmittedWorkDone');
      await outBlockBuffer.mapAsync(GPUMapMode.READ);
      await outFoundResultsBuffer.mapAsync(GPUMapMode.READ);

      //console.time('fillout');
      const foundArrayBuffer = outFoundResultsBuffer.getMappedRange();
      const foundArray = new Uint32Array(foundArrayBuffer);
      const numOfFound = foundArray[0];
      const outIs = new Uint32Array(numOfFound);
      const outJs = new Uint32Array(numOfFound);
      const outSims = new Float32Array(numOfFound);
      const blockArrayBuffer = outBlockBuffer.getMappedRange();
      const isArray = new Uint32Array(blockArrayBuffer, 0, numOfFound);
      const jsArray = new Uint32Array(blockArrayBuffer, neededThreads * 4, numOfFound);
      const simsArray = new Float32Array(blockArrayBuffer, neededThreads * 8, numOfFound);
      outIs.set(isArray);
      outJs.set(jsArray);
      outSims.set(simsArray);

      outSparseIindexes.push(outIs);
      outSparseJindexes.push(outJs);
      outSparseSimilarities.push(outSims);

      //console.timeEnd('fillout');
      outFoundResultsBuffer.unmap();
      outBlockBuffer.unmap();
    }
  }

  // destroy
  outFoundResultsBuffer.destroy();
  sparseKNNBuffer.destroy();
  startEndBuffer.destroy();
  resultBlockBuffer.destroy();
  outBlockBuffer.destroy();
  fouundResultsBuffer.destroy();
  return expandResToSparseForm(outSparseIindexes, outSparseJindexes, outSparseSimilarities, nRows);
}


function expandResToSparseForm(is: Uint32Array[], js: Uint32Array[], similarities: Float32Array[], nRows: number): MCLOpReturnType {
  console.time('count');
  const occurances = new Uint32Array(nRows);
  for (let i = 0; i < is.length; i++) {
    for (let j = 0; j < is[i].length; j++) {
      // only increment the occurance of the row, since we are expanding the rows
      occurances[is[i][j]]++;
    }
  }
  const offsetForm = toOffsetForm(occurances);

  const totalLength = offsetForm[nRows];
  const outIndexes = new Uint32Array(totalLength);
  const outSims = new Float32Array(totalLength);
  const insertCounter = new Uint32Array(nRows).fill(0); // for counting how many things we inserted into the knn arrays.

  for (let i = 0; i < is.length; i++) {
    for (let j = 0; j < is[i].length; j++) {
      const row = is[i][j];
      const col = js[i][j];
      const similarity = similarities[i][j];
      outIndexes[offsetForm[row] + insertCounter[row]] = col;
      outSims[offsetForm[row] + insertCounter[row]] = similarity;
      insertCounter[row]++;
    }
  }
  console.timeEnd('count');
  return {KNNIndexes: outIndexes, KNNSimilarities: outSims, indexOffsets: offsetForm};
}

