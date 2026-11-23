import {getGPUDevice} from '../getGPUDevice';
import {optimizeLayoutWGSL} from './wgsl/optimization-loop.wgsl';
import {WORKGROUP_DIM} from '../workgroup';

/** An epoch runs its edges in this many steps, edge i in step i % STEPS_PER_EPOCH, each step reading the embeddings
 * the previous one left. With one step, a vertex took all its moves of an epoch from the same position, and a layout of
 * 30K rows kept 18% fewer of its nearest neighbors than the CPU one; with 8 it keeps as many. */
const STEPS_PER_EPOCH = 8;

export async function optimizationLoop(
  head: Uint32Array, tail: Uint32Array, embeddingSize: number,
  epochsPerSample: Float32Array, epochsPerNegativeSample: Float32Array, initialAlpha: number, gamma: number,
  a: number, b: number, dim: number, nEpochs: number, nVertices: number, random: () => number
) {
  const threadsPerWorkgroupDim = WORKGROUP_DIM;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const threadsPerStep = Math.ceil(head.length / STEPS_PER_EPOCH);
  const requiredWorkgroups = Math.ceil(threadsPerStep / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));
  const threadsPerYDimension = threadsPerWorkgroupDim * numOfWorkgroupsPerDim;
  const device = await getGPUDevice();
  if (!device)
    return;

  const divisionFactor = 10_000_000;
  const headSize = head.length;
  const processWGSL = optimizeLayoutWGSL(headSize, embeddingSize, true, initialAlpha,
    gamma, a, b, dim, nEpochs, nVertices, threadsPerWorkgroupDim, threadsPerYDimension, divisionFactor,
    STEPS_PER_EPOCH, threadsPerStep
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


  const comuteInfo32Size = 3;
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
    new Int32Array(embeddingSize * dim).map(() => Math.floor((random() * 2 -1) * divisionFactor));
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

  const embeddingSnapshotBuffer = device.createBuffer({
    label: 'embedding snapshot buffer',
    size: paddedembedStorageSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
  });


  const computeInfoBuffer = device.createBuffer({
    label: 'compute info buffer',
    size: paddedComputeInfoSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const computeInfoArrayBuffer = computeInfoBuffer.getMappedRange();
  const computeInfoMappedArray = new Float32Array(computeInfoArrayBuffer, 0, 2);
  computeInfoMappedArray.set([0.0, initialAlpha]);
  new Uint32Array(computeInfoArrayBuffer, 8, 1)[0] = Math.floor(random() * 0x100000000);

  computeInfoBuffer.unmap();

  // the step indexes, one per bind group, each at an offset the device can bind
  const stepStride = device.limits.minStorageBufferOffsetAlignment;
  const stepBuffer = device.createBuffer({
    label: 'step buffer',
    size: STEPS_PER_EPOCH * stepStride,
    usage: GPUBufferUsage.STORAGE,
    mappedAtCreation: true,
  });
  const stepArray = new Uint32Array(stepBuffer.getMappedRange());
  for (let step = 0; step < STEPS_PER_EPOCH; step++)
    stepArray[step * stepStride / 4] = step;
  stepBuffer.unmap();

  const bindGroups = Array.from({length: STEPS_PER_EPOCH}, (_, step) => device.createBindGroup({
    label: `bindGroup for optimize step ${step}`,
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: matrixStorageBuffer}},
      {binding: 1, resource: {buffer: epochStorageBuffer}},
      {binding: 2, resource: {buffer: embeddingStorageBuffer}},
      {binding: 3, resource: {buffer: computeInfoBuffer}},
      {binding: 4, resource: {buffer: embeddingSnapshotBuffer}},
      {binding: 5, resource: {buffer: stepBuffer, offset: step * stepStride, size: 4}},
    ],
  }));


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
    for (let step = 0; step < STEPS_PER_EPOCH; step++) {
      encoder.copyBufferToBuffer(embeddingStorageBuffer, 0, embeddingSnapshotBuffer, 0, embeddingSnapshotBuffer.size);
      const pass = encoder.beginComputePass({
        label: 'matrix ops compute pass',
      });
      pass.setPipeline(pipeline);
      pass.setBindGroup(0, bindGroups[step]);
      pass.dispatchWorkgroups(
        numOfWorkgroupsPerDim,
        numOfWorkgroupsPerDim
      );
      pass.end();
    }

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
  embeddingSnapshotBuffer.destroy();
  stepBuffer.destroy();
  computeInfoBuffer.destroy();
  outEmbedViewBuffer.destroy();

  return result;
}
