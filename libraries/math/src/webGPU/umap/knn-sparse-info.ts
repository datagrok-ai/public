import {getGPUDevice} from '../getGPUDevice';
import {SparseKNNInfo} from './types';
import {knnMatrixOpInfoWGSL} from './wgsl/knn-sparse-info.wgsl';

/**
 * this will get all the info needed for subsequent operation in umap
 * like sizes of each array of transposed knn and sizes of union of knn and its transpose rows.
 * after this point, all knn structures will be in single array form, with corresponding offsets.
 * @param knnIndexes
 * @returns
 */
export async function getKnnSparseOpInfo(
  knnIndexes: number[][] | Int32Array[] | Uint32Array[]
): Promise<SparseKNNInfo | undefined> {
  const threadsPerWorkgroupDim = 10;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const numOfEntries = knnIndexes.length;
  const requiredWorkgroups = Math.ceil(numOfEntries / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));
  const knnSize = knnIndexes[0].length;
  const processWgsl = knnMatrixOpInfoWGSL(
    threadsPerWorkgroupDim, threadsPerWorkgroupDim * numOfWorkgroupsPerDim, numOfEntries, knnSize
  );

  const indexesBuffer32Size = numOfEntries * knnIndexes[0].length;

  const device = await getGPUDevice();
  if (device == null)
    return;


  const module = device.createShaderModule({
    label: 'transposedSizesCalulation',
    code: processWgsl
  });

  const pipeline = device.createComputePipeline({
    label: 'transposedSizesPipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'countTransposedCols',
    },
  });

  const indexBuffer = device.createBuffer({
    label: 'indexes buffer',
    size: indexesBuffer32Size * 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedIndexesBuffer = indexBuffer.getMappedRange();
  const mappedIndexesArray = new Int32Array(mappedIndexesBuffer);
  for (let i = 0; i < knnIndexes.length; i ++)
    mappedIndexesArray.set(knnIndexes[i], i * knnIndexes[i].length);

  indexBuffer.unmap();

  const sizesBuffer = device.createBuffer({
    label: 'transpose sizes buffer',
    size: knnIndexes.length * 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
  });
  const unionSizesBuffer = device.createBuffer({
    label: 'union sizes buffer',
    size: knnIndexes.length * 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const mappedUnionSizesBuffer = unionSizesBuffer.getMappedRange();
  const mappedUnionSizesArray = new Uint32Array(mappedUnionSizesBuffer);
  mappedUnionSizesArray.fill(knnSize);
  unionSizesBuffer.unmap();
  const bindGroup = device.createBindGroup({
    label: 'bindGroup for count ops',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: indexBuffer}},
      {binding: 1, resource: {buffer: sizesBuffer}},
      {binding: 2, resource: {buffer: unionSizesBuffer}}
      // { binding: 4, resource: { buffer: arraySizesBuffer } },
    ],
  });

  const encoder = device.createCommandEncoder({
    label: 'matrix ops encoder',
  });

  const pass = encoder.beginComputePass({
    label: 'matrix ops compute pass',
  });
  pass.setPipeline(pipeline);
  pass.setBindGroup(0, bindGroup);
  pass.dispatchWorkgroups(
    numOfWorkgroupsPerDim,
    numOfWorkgroupsPerDim
  );
  pass.end();


  const resultTransposeSizesBuffer = device.createBuffer({
    label: 'result transpose sizes',
    size: sizesBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  const resultUnionSizesBuffer = device.createBuffer({
    label: 'result union sizes',
    size: unionSizesBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  encoder.copyBufferToBuffer(
    sizesBuffer,
    0,
    resultTransposeSizesBuffer,
    0,
    resultTransposeSizesBuffer.size
  );

  encoder.copyBufferToBuffer(
    unionSizesBuffer,
    0,
    resultUnionSizesBuffer,
    0,
    resultUnionSizesBuffer.size
  );

  const commandBuffer = encoder.finish();
  device.queue.submit([commandBuffer]);
  await device.queue.onSubmittedWorkDone();

  await resultTransposeSizesBuffer.mapAsync(GPUMapMode.READ);
  await resultUnionSizesBuffer.mapAsync(GPUMapMode.READ);
  const resultTransposeSizesArrayBuffer = resultTransposeSizesBuffer.getMappedRange();
  const resultUnionSizesArrayBuffer = resultUnionSizesBuffer.getMappedRange();

  const resultTransposeSizesArray = new Uint32Array(knnIndexes.length);
  resultTransposeSizesArray.set(new Uint32Array(resultTransposeSizesArrayBuffer, 0, knnIndexes.length));
  resultTransposeSizesBuffer.unmap();

  const resultUnionSizesArray = new Uint32Array(knnIndexes.length);
  resultUnionSizesArray.set(new Uint32Array(resultUnionSizesArrayBuffer, 0, knnIndexes.length));
  resultUnionSizesBuffer.unmap();

  const unionMatrixSize = resultUnionSizesArray.reduce((old, val) => old + val, 0);


  resultUnionSizesBuffer.destroy();
  resultTransposeSizesBuffer.destroy();
  sizesBuffer.destroy();
  indexBuffer.destroy();
  unionSizesBuffer.destroy();

  // sizes of each transposed knn row, sizes of each union of knn and its transpose row, total size of union matrix
  return {resultTransposeSizesArray, resultUnionSizesArray, unionMatrixSize};

  //counteTransposedCols
}
