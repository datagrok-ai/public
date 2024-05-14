import {getGPUDevice} from '../getGPUDevice';
import {optimizeLayoutWGSL} from './wgsl/optimization-loop.wgsl';

export async function optimizationLoop(
  head: Uint32Array, tail: Uint32Array, embeddingSize: number,
  epochsPerSample: Float32Array, epochsPerNegativeSample: Float32Array, initialAlpha: number, gamma: number,
  a: number, b: number, dim: number, nEpochs: number, nVertices: number
) {
  const threadsPerWorkgroupDim = 10;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const requiredWorkgroups = Math.ceil(head.length / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));
  const threadsPerYDimension = threadsPerWorkgroupDim * numOfWorkgroupsPerDim;
  const device = await getGPUDevice();
  if (!device)
    return;

  const divisionFactor = 10_000_000;
  const headSize = head.length;
  const processWGSL = optimizeLayoutWGSL(headSize, embeddingSize, true, initialAlpha,
    gamma, a, b, dim, nEpochs, nVertices, threadsPerWorkgroupDim, threadsPerYDimension, divisionFactor
  );


  //          head/tail
  const matrixStorage32Size = headSize * 2;
  let paddedMatrixStorageSize = matrixStorage32Size * 4;
  const remainder1 = paddedMatrixStorageSize & 15;
  if (remainder1 !== 0)
    paddedMatrixStorageSize += 16 - remainder1;


  const epochStorage32Size = headSize * 4;
  let paddedepochStorageSize = epochStorage32Size * 4;
  const remainder2 = paddedepochStorageSize & 15;
  if (remainder2 !== 0)
    paddedepochStorageSize += 16 - remainder2;


  const embedStorage32Size = embeddingSize * 2 * 2;
  let paddedembedStorageSize = embedStorage32Size * 4;
  const remainder3 = paddedembedStorageSize & 15;
  if (remainder3 !== 0)
    paddedembedStorageSize += 16 - remainder3;


  const comuteInfo32Size = 2 + headSize;
  let paddedComputeInfoSize = comuteInfo32Size * 4;
  const remainder4 = paddedComputeInfoSize & 15;
  if (remainder4 !== 0)
    paddedComputeInfoSize += 16 - remainder4;


  const module = device.createShaderModule({
    label: 'optimizeShader',
    code: processWGSL
  });


  const pipeline = device.createComputePipeline({
    label: 'optimizePipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'optimizeStep',
    },
  });


  const matrixStorageBuffer = device.createBuffer({
    label: 'matrix storage buffer',
    size: paddedMatrixStorageSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const matrixStorageArrayBuffer = matrixStorageBuffer.getMappedRange();

  const headView = new Uint32Array(matrixStorageArrayBuffer, 0, headSize);
  headView.set(head);
  const tailView = new Uint32Array(matrixStorageArrayBuffer, headSize * 4, headSize);
  tailView.set(tail);

  matrixStorageBuffer.unmap();


  const epochStorageBuffer = device.createBuffer({
    label: 'epoch storage buffer',
    size: paddedepochStorageSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const epochStorageArrayBuffer = epochStorageBuffer.getMappedRange();

  // all are float32 so we get the whole view
  const epochStorageArrayView = new Float32Array(epochStorageArrayBuffer);

  // first epochsPerSample
  epochStorageArrayView.set(epochsPerSample, 0);
  //epochOfNextSample, set to same
  epochStorageArrayView.set(epochsPerSample, headSize);
  //epochsPerNegativeSample
  epochStorageArrayView.set(epochsPerNegativeSample, headSize*2);
  //epochOfNextNegativeSample
  epochStorageArrayView.set(epochsPerNegativeSample, headSize*3);
  epochStorageBuffer.unmap();


  const headEmbeddings =
    new Int32Array(embeddingSize * dim).map(() => Math.floor((Math.random() * 2 -1) * divisionFactor));
  const tailEmbeddings = new Int32Array(embeddingSize * dim);
  tailEmbeddings.set(headEmbeddings);
  const embeddingStorageBuffer = device.createBuffer({
    label: 'embedding storage buffer',
    size: paddedembedStorageSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedEmbeddingStorageBuffer = embeddingStorageBuffer.getMappedRange();
  const mappedEmbeddingStorageArray = new Int32Array(mappedEmbeddingStorageBuffer);
  mappedEmbeddingStorageArray.set(headEmbeddings, 0);
  mappedEmbeddingStorageArray.set(headEmbeddings, headEmbeddings.length);
  embeddingStorageBuffer.unmap();


  const computeInfoBuffer = device.createBuffer({
    label: 'compute info buffer',
    size: paddedComputeInfoSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const randomNumbers = new Uint32Array(headSize).map(() => Math.floor(Math.random() * divisionFactor));
  const computeInfoArrayBuffer = computeInfoBuffer.getMappedRange();
  const computeInfoMappedArray = new Float32Array(computeInfoArrayBuffer, 0, 2);
  computeInfoMappedArray.set([0.0, initialAlpha]);
  const randomNumbersView = new Uint32Array(computeInfoArrayBuffer, 8, headSize);
  randomNumbersView.set(randomNumbers);

  computeInfoBuffer.unmap();

  const bindGroup = device.createBindGroup({
    label: 'bindGroup for optimize',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: matrixStorageBuffer}},
      {binding: 1, resource: {buffer: epochStorageBuffer}},
      {binding: 2, resource: {buffer: embeddingStorageBuffer}},
      {binding: 3, resource: {buffer: computeInfoBuffer}}
    ],
  });


  let alpha = initialAlpha;


  const outEmbedViewBuffer = device.createBuffer({
    label: 'result buffer',
    size: embeddingStorageBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  for (let iter = 0; iter < nEpochs; iter++) {
    // write n and alpha


    device.queue.writeBuffer(computeInfoBuffer, 0, new Float32Array([iter, alpha]), 0);

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

    // here, we copy the head/tail views to their correct places in embeddings.

    if (iter === nEpochs - 1)
      encoder.copyBufferToBuffer(embeddingStorageBuffer, 0, outEmbedViewBuffer, 0, outEmbedViewBuffer.size);

    const commandBuffer = encoder.finish();
    device.queue.submit([commandBuffer]);
    await device.queue.onSubmittedWorkDone();


    alpha = initialAlpha * (1.0 - iter / nEpochs);
  }
  await outEmbedViewBuffer.mapAsync(GPUMapMode.READ);

  const resArrayBuffer = outEmbedViewBuffer.getMappedRange();

  const headEmbeddingsView = new Int32Array(resArrayBuffer, 0, headEmbeddings.length);
  headEmbeddings.set(headEmbeddingsView);

  outEmbedViewBuffer.unmap();

  const result = new Array(dim).fill(null).map((_ ) => new Float32Array(embeddingSize));
  for (let i = 0; i < embeddingSize; i++) {
    for (let j = 0; j < dim; j++)
      result[j][i] = headEmbeddings[i * dim + j] / divisionFactor;
  }

  matrixStorageBuffer.destroy();
  epochStorageBuffer.destroy();
  embeddingStorageBuffer.destroy();
  computeInfoBuffer.destroy();
  outEmbedViewBuffer.destroy();

  return result;
}
