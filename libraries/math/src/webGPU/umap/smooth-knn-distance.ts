import {getGPUDevice} from '../getGPUDevice';
import {smoothKNNDistanceWGSL} from './wgsl/smooth-knn-distance.wgsl';

export async function smoothKNNDistance(
  knnDistances: number[][] | Float32Array[], nNeighbors: number = 15, localConnectivity: number = 1.0,
  nIter: number = 64, bandwidth: number = 1.0
) {
  const threadsPerWorkgroupDim = 10;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const numOfEntries = knnDistances.length;
  const requiredWorkgroups = Math.ceil(numOfEntries / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));

  const processWgsl = smoothKNNDistanceWGSL(
    threadsPerWorkgroupDim, threadsPerWorkgroupDim * numOfWorkgroupsPerDim,
    knnDistances, nNeighbors, localConnectivity, nIter, bandwidth);

  const distanceBuffer32Size = numOfEntries * knnDistances[0].length;

  // const flatDistanceArray = new Float32Array(distanceBuffer32Size);
  // for (let i = 0; i < knnDistances.length; i ++) {
  //     flatDistanceArray.set(knnDistances[i], i * knnDistances[i].length);
  // }

  const device = await getGPUDevice();
  if (device == null)
    return;


  const resultSigmas = new Float32Array(numOfEntries);
  const resultRhos = new Float32Array(numOfEntries);

  const module = device.createShaderModule({
    label: 'sigmaRhoCalulation',
    code: processWgsl
  });

  const pipeline = device.createComputePipeline({
    label: 'sigmaRhoPipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'smoothKNNDistance',
    },
  });

  const distancesBuffer = device.createBuffer({
    label: 'distance buffer',
    size: distanceBuffer32Size * 4,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedDistanceBuffer = distancesBuffer.getMappedRange();
  const mappedDistanceArray = new Float32Array(mappedDistanceBuffer);
  for (let i = 0; i < knnDistances.length; i ++)
    mappedDistanceArray.set(knnDistances[i], i * knnDistances[i].length);

  distancesBuffer.unmap();

  const bufferSigmas = device.createBuffer({
    label: 'buffer sigmas',
    size: numOfEntries * 4,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
  });

  const bufferRhos = device.createBuffer({
    label: 'buffer Rhos',
    size: numOfEntries * 4,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
  });
  const bindGroup = device.createBindGroup({
    label: 'bindGroup for sigmarho buffers',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: bufferSigmas}},
      {binding: 1, resource: {buffer: bufferRhos}},
      {binding: 2, resource: {buffer: distancesBuffer}},
      // { binding: 4, resource: { buffer: arraySizesBuffer } },
    ],
  });

  const encoder = device.createCommandEncoder({
    label: 'sigmarho encoder',
  });

  const pass = encoder.beginComputePass({
    label: 'sigmarho compute pass',
  });
  pass.setPipeline(pipeline);
  pass.setBindGroup(0, bindGroup);
  pass.dispatchWorkgroups(
    numOfWorkgroupsPerDim,
    numOfWorkgroupsPerDim
  );
  pass.end();

  const resultBufferRhos = device.createBuffer({
    label: 'result buffer rhos',
    size: bufferRhos.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });
  const resultBufferSigmas = device.createBuffer({
    label: 'result buffer sigmas',
    size: bufferSigmas.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  encoder.copyBufferToBuffer(
    bufferRhos,
    0,
    resultBufferRhos,
    0,
    resultBufferRhos.size
  );
  encoder.copyBufferToBuffer(
    bufferSigmas,
    0,
    resultBufferSigmas,
    0,
    resultBufferSigmas.size
  );
  const commandBuffer = encoder.finish();
  device.queue.submit([commandBuffer]);
  await device.queue.onSubmittedWorkDone();

  await resultBufferRhos.mapAsync(GPUMapMode.READ);
  await resultBufferSigmas.mapAsync(GPUMapMode.READ);

  resultRhos.set(new Float32Array(resultBufferRhos.getMappedRange()));
  resultSigmas.set(new Float32Array(resultBufferSigmas.getMappedRange()));

  resultBufferRhos.unmap();
  resultBufferSigmas.unmap();

  device.destroy();

  return {resultRhos, resultSigmas};
}
