import {getGPUDevice} from '../getGPUDevice';
import {computeMembershipStrengthsWGSL} from './wgsl/membership-strengths.wgsl';

export async function computeMembershipStrengths(
  knnDistances: number[][] | Float32Array[], sigmas: number[] | Float32Array, rhos: number[] | Float32Array
) {
  const threadsPerWorkgroupDim = 10;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const numOfEntries = knnDistances.length;
  const requiredWorkgroups = Math.ceil(numOfEntries / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));

  const processWgsl = computeMembershipStrengthsWGSL(
    threadsPerWorkgroupDim, threadsPerWorkgroupDim * numOfWorkgroupsPerDim,
    knnDistances, sigmas, rhos);

  const distanceBuffer32Size = numOfEntries * knnDistances[0].length;


  const device = await getGPUDevice();
  if (device == null)
    return;


  const module = device.createShaderModule({
    label: 'rowsColsValsCalulation',
    code: processWgsl
  });

  const pipeline = device.createComputePipeline({
    label: 'rowsColsValsPipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'computeMembershipStrengths',
    },
  });

  const rhos32Size = rhos.length;
  const sigmas32Size = sigmas.length;

  const combinedMembershipStrengthsInfo32Size = rhos32Size + sigmas32Size + distanceBuffer32Size;
  let paddedCombinedMembershipStructByteSize = combinedMembershipStrengthsInfo32Size * 4;
  const remainder = paddedCombinedMembershipStructByteSize & 15;
  if (remainder != 0)
    paddedCombinedMembershipStructByteSize += 16 - remainder;


  // input params
  const membershipStrengthsInfoBuffer = device.createBuffer({
    label: 'info buffer',
    size: paddedCombinedMembershipStructByteSize,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedInfoBuffer = membershipStrengthsInfoBuffer.getMappedRange();
  const mappedInfoArray = new Float32Array(mappedInfoBuffer);
  // first set knn distances
  for (let i = 0; i < knnDistances.length; i ++)
    mappedInfoArray.set(knnDistances[i], i * knnDistances[i].length);

  let offset = knnDistances[0].length * knnDistances.length;
  // now set sigmas
  mappedInfoArray.set(sigmas, offset);
  offset += sigmas.length;
  // now set rhos
  mappedInfoArray.set(rhos, offset);
  membershipStrengthsInfoBuffer.unmap();

  const knnDistancesBuffer = device.createBuffer({
    label: 'knn distance buffer',
    size: knnDistances[0].length * knnDistances.length * Float32Array.BYTES_PER_ELEMENT,
    usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC |
            GPUBufferUsage.COPY_DST,
  });


  const bindGroup = device.createBindGroup({
    label: 'bindGroup for membership strhegths buffers',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: membershipStrengthsInfoBuffer}},
      {binding: 1, resource: {buffer: knnDistancesBuffer}},
      // { binding: 4, resource: { buffer: arraySizesBuffer } },
    ],
  });

  const encoder = device.createCommandEncoder({
    label: 'membership strengths encoder',
  });

  const pass = encoder.beginComputePass({
    label: 'membership strengths compute pass',
  });
  pass.setPipeline(pipeline);
  pass.setBindGroup(0, bindGroup);
  pass.dispatchWorkgroups(
    numOfWorkgroupsPerDim,
    numOfWorkgroupsPerDim
  );
  pass.end();

  const resultBufferDistances = device.createBuffer({
    label: 'result buffer distances',
    size: knnDistancesBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  encoder.copyBufferToBuffer(
    knnDistancesBuffer,
    0,
    resultBufferDistances,
    0,
    resultBufferDistances.size
  );

  const commandBuffer = encoder.finish();
  device.queue.submit([commandBuffer]);
  await device.queue.onSubmittedWorkDone();

  await resultBufferDistances.mapAsync(GPUMapMode.READ);

  const resultKnnDistancesArrayBuffer = resultBufferDistances.getMappedRange();

  const bytesPerKnnRow = knnDistances[0].length * Float32Array.BYTES_PER_ELEMENT;
  const resKnnDistances = new Array(knnDistances.length).fill(null).map((_, i) => {
    const dataView = new Float32Array(resultKnnDistancesArrayBuffer, bytesPerKnnRow * i, knnDistances[0].length);
    const out = new Float32Array(knnDistances[0].length);
    out.set(dataView);
    return out;
  });

  resultBufferDistances.unmap();

  resultBufferDistances.destroy();
  knnDistancesBuffer.destroy();
  membershipStrengthsInfoBuffer.destroy();

  return resKnnDistances;
}
