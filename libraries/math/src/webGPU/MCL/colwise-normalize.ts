
/**
 * colwise normalize the knn similarities
 * BE WARE, modifies the input array in place
 * @param device
 * @param knnSimilarities
 * @param offsets
 * @param nRows
 */
export async function sparseKNNNoralizeColwise(
  device: GPUDevice,
  knnSimilarities: Float32Array,
  offsets: Uint32Array,
  nRows: number
) {
  const neededThreads = nRows;
  const workGroupThreadsPerDim = 10;
  const totalWorkgroupThreads = workGroupThreadsPerDim * workGroupThreadsPerDim;
  const neededWorkGroups = Math.ceil(neededThreads / totalWorkgroupThreads);
  const workGroupDim = Math.ceil(Math.sqrt(neededWorkGroups));
  const threadsPerDim = workGroupThreadsPerDim * workGroupDim;

  const module = device.createShaderModule({
    label: 'colwise-normalize',
    code: `
        @group(0) @binding(0) var<storage, read_write> knnSimilarities : array<f32>;
        @group(0) @binding(1) var<storage, read> offsets : array<u32>;
        @compute @workgroup_size(${workGroupThreadsPerDim}, ${workGroupThreadsPerDim}) fn normalize(
            @builtin(global_invocation_id) id: vec3<u32>,
          ) {
            let row = id.x;
            let col = id.y;
            let index = row * ${threadsPerDim} + col;
            if (index >= ${nRows}) {
                return;
            }
            let offsetBegin = offsets[index];
            let offsetEnd = offsets[index + 1];
            var sum = 0.0;
            if (offsetEnd - offsetBegin == 0) {
                return;
            }
            for (var i = offsetBegin; i < offsetEnd; i = i + 1) {
                sum = sum + knnSimilarities[i];
            }
            if (sum > 0.0) {
              for (var i = offsetBegin; i < offsetEnd; i = i + 1) {
                  knnSimilarities[i] = knnSimilarities[i] / sum;
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
      entryPoint: 'normalize',
    },
  });

  const simmilaritiesBuffer = device.createBuffer({
    label: 'simmilarities info buffer',
    size: knnSimilarities.byteLength,
    usage:
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_SRC |
      GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  new Float32Array(simmilaritiesBuffer.getMappedRange()).set(knnSimilarities);
  simmilaritiesBuffer.unmap();

  const offsetsBuffer = device.createBuffer({
    label: 'offsets info buffer',
    size: offsets.byteLength,
    usage:
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_SRC |
      GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  new Uint32Array(offsetsBuffer.getMappedRange()).set(offsets);
  offsetsBuffer.unmap();

  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {
        binding: 0,
        resource: {
          buffer: simmilaritiesBuffer,
        },
      },
      {
        binding: 1,
        resource: {
          buffer: offsetsBuffer,
        },
      },
    ],
  });

  const commandEncoder = device.createCommandEncoder();
  const passEncoder = commandEncoder.beginComputePass();
  passEncoder.setPipeline(pipeline);
  passEncoder.setBindGroup(0, bindGroup);
  passEncoder.dispatchWorkgroups(workGroupDim, workGroupDim);
  passEncoder.end();

  const outSimilarityBuffer = device.createBuffer({
    label: 'out similarity buffer',
    size: knnSimilarities.byteLength,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  commandEncoder.copyBufferToBuffer(
    simmilaritiesBuffer,
    0,
    outSimilarityBuffer,
    0,
    knnSimilarities.byteLength
  );
  device.queue.submit([commandEncoder.finish()]);
  await device.queue.onSubmittedWorkDone();
  await outSimilarityBuffer.mapAsync(GPUMapMode.READ);
  const arrayBuffer = new Float32Array(outSimilarityBuffer.getMappedRange());
  knnSimilarities.set(arrayBuffer);
  outSimilarityBuffer.unmap();
  simmilaritiesBuffer.destroy();
  offsetsBuffer.destroy();
  outSimilarityBuffer.destroy();
}
