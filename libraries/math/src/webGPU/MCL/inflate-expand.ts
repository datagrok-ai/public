/* eslint-disable max-len */
import {MCLOpReturnType} from './types';

/**
 * inflation operator in place
 */
export function inflate(knnSimilarities: Float32Array, factor: number = 2) {
  for (let i = 0; i < knnSimilarities.length; i++)
    knnSimilarities[i] = Math.pow(knnSimilarities[i], factor);
}

function getRowIndexes(offsets: Uint32Array) {
  const res = new Uint32Array(offsets[offsets.length - 1]);
  for (let i = 0; i < offsets.length - 1; i++) {
    for (let j = offsets[i]; j < offsets[i + 1]; j++)
      res[j] = i;
  }
  return res;
}


// this implementation is not looking at already zeroed out values, and only looks at live cells.
export async function expandNoRevive(
  device: GPUDevice, knnSimilarities: Float32Array, knnIndexes: Uint32Array, offsets: Uint32Array, _nRows: number
): Promise<MCLOpReturnType> {
  const neededThreads = 90000;
  const workGroupThreadsPerDim = 10;
  const totalWorkgroupThreads = workGroupThreadsPerDim * workGroupThreadsPerDim;
  const neededWorkGroups = Math.ceil(neededThreads / totalWorkgroupThreads);
  const workGroupDim = Math.ceil(Math.sqrt(neededWorkGroups));
  const threadsPerDim = workGroupThreadsPerDim * workGroupDim;

  //const order = Math.floor(Math.max(Math.log(nRows), 2));
  // minimum value after expansion.
  const pruneValue = 0.000000001;//Math.pow(10, -order);

  const outKNNSimilarities = new Float32Array(knnSimilarities.length);
  const module = device.createShaderModule({
    label: 'expand',
    code: `
        struct SparseKNN {
            knnSimilarities: array<f32, ${knnSimilarities.length}>,
            knnIndexes: array<u32, ${knnIndexes.length}>,
            offsets: array<u32, ${offsets.length}>,
            rowIndexes: array<u32, ${knnIndexes.length}>,
        }

        @group(0) @binding(0) var<storage, read_write> sparseKNN: SparseKNN;
        @group(0) @binding(1) var<storage, read_write> resultSimBlock: array<f32, ${neededThreads}>;
        @group(0) @binding(2) var<storage, read_write> startAt: u32;
        @compute @workgroup_size(${workGroupThreadsPerDim}, ${workGroupThreadsPerDim}) fn expand(
            @builtin(global_invocation_id) id: vec3<u32>,
          ) {
            let col = id.x;
            let row = id.y;
            let index = row * ${threadsPerDim} + col;
            if (index >= ${neededThreads}) {
                return;
            }
            let workingIndex = index + startAt;
            if (workingIndex >= ${knnIndexes.length}) {
                return;
            }
            
            let rowIdx = sparseKNN.rowIndexes[workingIndex];
            let colIdx = sparseKNN.knnIndexes[workingIndex];
            let offsetBeginRow = sparseKNN.offsets[rowIdx];
            let offsetEndRow = sparseKNN.offsets[rowIdx + 1];
            let offsetBeginCol = sparseKNN.offsets[colIdx];
            let offsetEndCol = sparseKNN.offsets[colIdx + 1];
            var sum = 0.0;
            for (var i = offsetBeginRow; i < offsetEndRow; i = i + 1) {
                for(var j = offsetBeginCol; j < offsetEndCol; j = j + 1) {
                    if (sparseKNN.knnIndexes[i] == sparseKNN.knnIndexes[j]) {
                        sum = sum + sparseKNN.knnSimilarities[i] * sparseKNN.knnSimilarities[j];
                        break;
                    }
                }
            }
            if (sum > ${pruneValue}) {
                resultSimBlock[index] = sum;
            } else {
                resultSimBlock[index] = 0.0;
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

  const sparseKNNBuffer32Size = knnSimilarities.length + knnIndexes.length + offsets.length + knnIndexes.length;
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
  // set row indexes
  const rowIndexes = getRowIndexes(offsets);
  new Uint32Array(sparseKNNArrayBuffer, (knnSimilarities.length + knnIndexes.length + offsets.length) * 4, rowIndexes.length).set(rowIndexes);

  sparseKNNBuffer.unmap();

  const startBuffer = device.createBuffer({
    label: 'start end buffer',
    size: 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  new Uint32Array(startBuffer.getMappedRange()).set([0]);
  startBuffer.unmap();

  const resultBlockBuffer = device.createBuffer({
    label: 'result block buffer',
    size: neededThreads * 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
  });

  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: sparseKNNBuffer}},
      {binding: 1, resource: {buffer: resultBlockBuffer}},
      {binding: 2, resource: {buffer: startBuffer}},
    ],
  });

  const outBlockBuffer = device.createBuffer({
    label: 'out block buffer',
    size: resultBlockBuffer.size,
    usage:
        GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
  });

  for (let i = 0; i < Math.ceil(knnIndexes.length / neededThreads); i++) {
    const start = i * neededThreads;
    const curLength = Math.min(neededThreads, knnIndexes.length - start);
    device.queue.writeBuffer(startBuffer, 0, new Uint32Array([start]));
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
    device.queue.submit([encoder.finish()]);
    //console.time('onSubmittedWorkDone');
    await device.queue.onSubmittedWorkDone();
    //console.timeEnd('onSubmittedWorkDone');
    await outBlockBuffer.mapAsync(GPUMapMode.READ);
    const outBlock = new Float32Array(outBlockBuffer.getMappedRange(), 0, curLength);
    //console.time('fillout');
    outKNNSimilarities.set(outBlock, start);
    //console.timeEnd('fillout');
    outBlockBuffer.unmap();
  }

  // destroy
  sparseKNNBuffer.destroy();
  startBuffer.destroy();
  resultBlockBuffer.destroy();
  outBlockBuffer.destroy();
  return {KNNIndexes: knnIndexes, KNNSimilarities: outKNNSimilarities, indexOffsets: offsets};
}
