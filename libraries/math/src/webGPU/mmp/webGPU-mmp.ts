/* eslint-disable max-len */
import {getGPUDevice} from '../getGPUDevice';
import {toOffsetForm} from '../umap/utils';

export async function webGPUMMP(frags: [number, number][][], fragSizes: Uint32Array, threshold: number = 0.4) {
  const device = await getGPUDevice();
  if (!device)
    return null; // if no device, return null, as we cannot do anything without it.

  const threadsPerWorkgroupDim = 10;
  const threadsInWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const neededThreads = 10000; //TBD
  const neededWorkgroups = Math.ceil(neededThreads / threadsInWorkgroup);
  const workgroupsPerDim = Math.ceil(Math.sqrt(neededWorkgroups));
  const totalThreadsPerDim = workgroupsPerDim * threadsPerWorkgroupDim; // threads per x axis

  const maxComputations = 5000;
  const maxFound = 100;

  const {frags1, frags2, offsets} = fragsToOffsetForm(frags);

  const listSize = offsets.length - 1;
  const module = device.createShaderModule({
    label: 'mmp compute shader',
    code: `
      struct MMPData {
        frags1: array<u32, ${frags1.length}>,
        frags2: array<u32, ${frags2.length}>,
        sizes: array<u32, ${fragSizes.length}>,
        offsets: array<u32, ${offsets.length}>
      }

      struct Result {
        molecules1: array<array<u32, ${maxFound}>,${neededThreads}>,
        molecules2: array<array<u32, ${maxFound}>,${neededThreads}>,
        cores: array<array<u32, ${maxFound}>,${neededThreads}>,
        frag1: array<array<u32, ${maxFound}>,${neededThreads}>,
        frag2: array<array<u32, ${maxFound}>,${neededThreads}>,
        found: array<u32, ${neededThreads}>,
        done: array<u32, ${neededThreads}>
      }

      struct ComputeInfo {
        // start at cols and rows, and end at cols and rows for each thread, these will be calculated on cpu and passed to gpu.
        startAtCols: array<u32, ${neededThreads}>,
        startAtRows: array<u32, ${neededThreads}>,
        endAtCols: array<u32, ${neededThreads}>,
        endAtRows: array<u32, ${neededThreads}>
      }

      @group(0) @binding(0) var<storage, read_write> mmpInfo: MMPData;
      @group(0) @binding(1) var<storage, read_write> res: Result;
      @group(0) @binding(2) var<storage, read_write> computeInfo: ComputeInfo;
      @compute @workgroup_size(${threadsPerWorkgroupDim}, ${threadsPerWorkgroupDim}) fn calcMMP(
        @builtin(global_invocation_id) id: vec3<u32>
      ) {
        
        let col = id.x;
        let row = id.y;
        let workingIndex = row * ${totalThreadsPerDim} + col;
        if (workingIndex >= ${neededThreads}) {
          return; // if we are out of bounds, return
        }

      var startAtCol: u32 = computeInfo.startAtCols[workingIndex];
      var startAtRow: u32 = computeInfo.startAtRows[workingIndex];
      let endAtCol: u32 = min(computeInfo.endAtCols[workingIndex], ${listSize}u);
      let endAtRow: u32 = min(computeInfo.endAtRows[workingIndex], ${listSize}u);
      res.found[workingIndex] = 0; // initialize the found counter
      var found: u32 = 0;
      if (res.done[workingIndex] > 0) {
        return; // if we are done, return
      }

      for (var i = 0; i < ${maxComputations}; i++) {
        if (startAtCol >= endAtCol && startAtRow >= endAtRow) {
          res.done[workingIndex] = 1;
          break;
        }

        if (found >= ${maxFound}) {
          break;
        }

        let molIdx1 = startAtCol;
        let molIdx2 = startAtRow;
        var core = -1i;
        var r1 = -1i;
        var r2 = -1i;
        var coreSize = 0i;
        for (var j = mmpInfo.offsets[molIdx1]; j < mmpInfo.offsets[molIdx1 + 1]; j++) {
          let frag1Core = mmpInfo.frags1[j];
          let frag1Frag = mmpInfo.frags2[j];
          for (var k = mmpInfo.offsets[molIdx2]; k < mmpInfo.offsets[molIdx2 + 1]; k++) {
            let frag2Core = mmpInfo.frags1[k];
            let frag2Frag = mmpInfo.frags2[k];
            if (frag1Core == frag2Core && i32(mmpInfo.sizes[frag1Core]) > coreSize) {
              coreSize = i32(mmpInfo.sizes[frag1Core]);
              r1 = i32(frag1Frag);
              r2 = i32(frag2Frag);
              core = i32(frag2Core);
            }
          }
        }
        if (coreSize > 0 && core != -1i &&
          f32(mmpInfo.sizes[r1]) / f32(coreSize) <= ${threshold} &&
          f32(mmpInfo.sizes[r2]) / f32(coreSize) <= ${threshold}) {
            res.molecules1[workingIndex][found] = molIdx1;
            res.molecules2[workingIndex][found] = molIdx2;
            res.cores[workingIndex][found] = u32(core);
            res.frag1[workingIndex][found] = u32(r1);
            res.frag2[workingIndex][found] = u32(r2);
            found = found + 1;
        }

        
        startAtCol = startAtCol + 1;
        if (startAtCol >= ${listSize}u) {
          startAtRow += 1;
          startAtCol = startAtRow + 1;
        }

      }

      res.found[workingIndex] = found;
      // update the startAtCols and startAtRows
      computeInfo.startAtCols[workingIndex] = startAtCol;
      computeInfo.startAtRows[workingIndex] = startAtRow;
    }
    `
  });

  const pipeline = device.createComputePipeline({
    label: 'mmp compute pipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'calcMMP',
    },
  });

  const startAtCols = new Uint32Array(neededThreads);
  const startAtRows = new Uint32Array(neededThreads);
  const endAtCols = new Uint32Array(neededThreads);
  const endAtRows = new Uint32Array(neededThreads);

  const condensedDistanceMatrixSize = (listSize) * (listSize - 1) / 2;
  const chunkSize = Math.floor(condensedDistanceMatrixSize / neededThreads); // size of the chunk per thread (in total)
  let startRow = 0;
  let startCol = 1;
  for (let i = 0; i < neededThreads; i++) {
    const endIdx = i === neededThreads - 1 ? condensedDistanceMatrixSize - 1 : (i + 1) * chunkSize;
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
  }

  const mppDataStruct32Size = frags1.length + frags2.length + fragSizes.length + offsets.length;

  let mppDataStructBiteSize = mppDataStruct32Size * Uint32Array.BYTES_PER_ELEMENT;
  const remainder = mppDataStructBiteSize & 15;
  if (remainder !== 0)
    mppDataStructBiteSize += 16 - remainder;

  //                 mol1, mol2, cores, frag1, frag2        // found, done
  const resultStruct32Size = maxFound * neededThreads * 5 + neededThreads *2;

  let resultStructBiteSize = resultStruct32Size * Uint32Array.BYTES_PER_ELEMENT;
  const rem1 = resultStructBiteSize & 15;
  if (rem1 !== 0)
    resultStructBiteSize += 16 - rem1;

  const computeInfoStruct32Size = neededThreads * 4;
  let computeInfoStructBiteSize = computeInfoStruct32Size * Uint32Array.BYTES_PER_ELEMENT;
  const rem2 = computeInfoStructBiteSize & 15;
  if (rem2 !== 0)
    computeInfoStructBiteSize += 16 - rem2;

  const mmpDataBuffer = device.createBuffer({
    label: 'mmp data buffer',
    size: mppDataStructBiteSize,
    usage: GPUBufferUsage.STORAGE |
    GPUBufferUsage.COPY_SRC |
    GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });


  const mmpDataArrayBuffer = mmpDataBuffer.getMappedRange();
  const mmpDataU32View = new Uint32Array(mmpDataArrayBuffer);
  mmpDataU32View.set(frags1, 0);
  mmpDataU32View.set(frags2, frags1.length);
  mmpDataU32View.set(fragSizes, frags1.length + frags2.length);
  mmpDataU32View.set(offsets, frags1.length + frags2.length + fragSizes.length);
  mmpDataBuffer.unmap();

  const resultBuffer = device.createBuffer({
    label: 'mmp result buffer',
    size: resultStructBiteSize,
    usage: GPUBufferUsage.STORAGE |
    GPUBufferUsage.COPY_SRC |
    GPUBufferUsage.COPY_DST,
  });

  const computeInfoBuffer = device.createBuffer({
    label: 'compute info buffer',
    size: computeInfoStructBiteSize,
    usage: GPUBufferUsage.STORAGE |
    GPUBufferUsage.COPY_SRC |
    GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const computInfoArrayBuffer = computeInfoBuffer.getMappedRange();
  const computeInfo32View = new Uint32Array(computInfoArrayBuffer);

  computeInfo32View.set(startAtCols, 0);
  computeInfo32View.set(startAtRows, neededThreads);
  computeInfo32View.set(endAtCols, neededThreads * 2);
  computeInfo32View.set(endAtRows, neededThreads * 3);
  computeInfoBuffer.unmap();

  const resultsOutBuffer = device.createBuffer({
    label: 'results out buffer',
    size: resultBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  const bindGroup = device.createBindGroup({
    label: 'bindGroup for mmp buffers',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: mmpDataBuffer}},
      {binding: 1, resource: {buffer: resultBuffer}},
      {binding: 2, resource: {buffer: computeInfoBuffer}},
    ],
  });

  const outMols1: Array<Uint32Array> = [];
  const outMols2: Array<Uint32Array> = [];
  const outFrag1: Array<Uint32Array> = [];
  const outFrag2: Array<Uint32Array> = [];
  const outCores: Array<Uint32Array> = [];

  let isAllDone = false;
  while (!isAllDone) {
    const encoder = device.createCommandEncoder({
      label: 'distance encoder',
    });
    const pass = encoder.beginComputePass({
      label: 'distance compute pass',
    });
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(
      workgroupsPerDim,
      workgroupsPerDim
    );
    pass.end();
    encoder.copyBufferToBuffer(
      resultBuffer,
      0,
      resultsOutBuffer,
      0,
      resultsOutBuffer.size
    );

    const commandBuffer = encoder.finish();
    device.queue.submit([commandBuffer]);

    // Read the results
    await device.queue.onSubmittedWorkDone();

    await resultsOutBuffer.mapAsync(GPUMapMode.READ);

    const resultsOutArrayBuffer = resultsOutBuffer.getMappedRange();
    let resultOffset = 0;
    const resMols1 = new Uint32Array(resultsOutArrayBuffer, resultOffset, maxFound * neededThreads);
    resultOffset += maxFound * neededThreads * 4;
    const resMols2 = new Uint32Array(resultsOutArrayBuffer, resultOffset, maxFound * neededThreads);
    resultOffset += maxFound * neededThreads * 4;
    const resCores = new Uint32Array(resultsOutArrayBuffer, resultOffset, maxFound * neededThreads);
    resultOffset += maxFound * neededThreads * 4;
    const resFrag1 = new Uint32Array(resultsOutArrayBuffer, resultOffset, maxFound * neededThreads);
    resultOffset += maxFound * neededThreads * 4;
    const resFrag2 = new Uint32Array(resultsOutArrayBuffer, resultOffset, maxFound * neededThreads);
    resultOffset += maxFound * neededThreads * 4;
    const resFound = new Uint32Array(resultsOutArrayBuffer, resultOffset, neededThreads);
    resultOffset += neededThreads * 4;
    const resDone = new Uint32Array(resultsOutArrayBuffer, resultOffset, neededThreads);
    resultOffset += neededThreads * 4;

    isAllDone = resDone.every((i) => i == 1);
    const totalResults = resFound.reduce((a, b) => a + b, 0);
    const combinedMols1 = new Uint32Array(totalResults);
    const combinedMols2 = new Uint32Array(totalResults);
    const combinedCores = new Uint32Array(totalResults);
    const combinedFrags1 = new Uint32Array(totalResults);
    const combinedFrags2 = new Uint32Array(totalResults);
    let combinedOffset = 0;
    for (let resI = 0; resI < resFound.length; resI++) {
      const found = resFound[resI];
      if (found === 0) continue;
      combinedMols1.set(resMols1.subarray(resI * maxFound, resI* maxFound + found), combinedOffset);
      combinedMols2.set(resMols2.subarray(resI * maxFound, resI* maxFound + found), combinedOffset);
      combinedCores.set(resCores.subarray(resI * maxFound, resI* maxFound + found), combinedOffset);
      combinedFrags1.set(resFrag1.subarray(resI * maxFound, resI* maxFound + found), combinedOffset);
      combinedFrags2.set(resFrag2.subarray(resI * maxFound, resI* maxFound + found), combinedOffset);
      combinedOffset += found;
    }
    outMols1.push(combinedMols1);
    outMols2.push(combinedMols2);
    outCores.push(combinedCores);
    outFrag1.push(combinedFrags1);
    outFrag2.push(combinedFrags2);
    resultsOutBuffer.unmap();
  }


  const totalSize = outMols1.reduce((a, b) => a + b.length, 0);
  const finalMols1 = new Uint32Array(totalSize);
  const finalMols2 = new Uint32Array(totalSize);
  const finalFrags1 = new Uint32Array(totalSize);
  const finalFrags2 = new Uint32Array(totalSize);
  const finalCores = new Uint32Array(totalSize);
  let finalOffset = 0;
  for (let i = 0; i < outMols1.length; i++) {
    finalMols1.set(outMols1[i], finalOffset);
    finalMols2.set(outMols2[i], finalOffset);
    finalFrags1.set(outFrag1[i], finalOffset);
    finalFrags2.set(outFrag2[i], finalOffset);
    finalCores.set(outCores[i], finalOffset);
    finalOffset += outMols1[i].length;
  }

  resultsOutBuffer.destroy();
  computeInfoBuffer.destroy();
  resultBuffer.destroy();
  mmpDataBuffer.destroy();

  return {finalMols1, finalMols2, finalFrags1, finalFrags2, finalCores};
}


function fragsToOffsetForm(frags: [number, number][][]) {
  const fragsTotalLength = frags.reduce((prev, cur) => prev + cur.length, 0);
  const frags1 = new Uint32Array(fragsTotalLength);
  const frags2 = new Uint32Array(fragsTotalLength);
  let curOffset = 0;
  for (let i = 0; i < frags.length; i++) {
    for (let j = 0; j < frags[i].length; j++) {
      frags1[curOffset] = frags[i][j][0];
      frags2[curOffset] = frags[i][j][1];
      curOffset ++;
    }
  }
  const offsets = toOffsetForm(frags.map((i) => i.length));
  return {offsets, frags1, frags2};
}
