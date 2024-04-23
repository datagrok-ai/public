export async function testAtomics() {
  const addapter = await navigator.gpu.requestAdapter();
  const device = await addapter?.requestDevice();
  if (!device)
    return;

  const module = device.createShaderModule({
    label: 'KNN compute shader',
    code: `
    
    // WGSL shader code snippet to use atomic operations
@group(0) @binding(0) var<storage, read_write> atomicArray : array<atomic<u32>>;
@group(0) @binding(1) var<storage, read_write> addsArray: array<u32>;
@compute @workgroup_size(10)
fn main(@builtin(global_invocation_id) GlobalInvocationID : vec3<u32>) {
    let index = GlobalInvocationID.x;
    let oldVal = atomicAdd(&atomicArray[0], 1);
    addsArray[index] = oldVal;
    }
    `});

  const bufferSize = 4 * 10000; // For example, an array of 100 uint32 atomic counters
  const atomicBuffer = device.createBuffer({
    size: bufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC,
  });
  const addsBuffer = device.createBuffer({
    size: bufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC,
  });

  const pipeline = device.createComputePipeline({
    label: 'hamming compute pipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'main',
    },
  });
  const bindGroup = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [{
      binding: 0,
      resource: {
        buffer: atomicBuffer,
      },

    }, {
      binding: 1, resource: {buffer: addsBuffer}}],
  });

  const encoder = device.createCommandEncoder({
    label: 'distance encoder',
  });
  const pass = encoder.beginComputePass({
    label: 'distance compute pass',
  });
  pass.setPipeline(pipeline);
  pass.setBindGroup(0, bindGroup);
  pass.dispatchWorkgroups(1000, 1);
  pass.end();

  const bufferDistances = device.createBuffer({
    label: 'buffer distances',
    size: atomicBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  const bufferAdds = device.createBuffer({
    label: 'buffer adds',
    size: atomicBuffer.size,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });
  encoder.copyBufferToBuffer(
    atomicBuffer,
    0,
    bufferDistances,
    0,
    bufferDistances.size
  );
  encoder.copyBufferToBuffer(
    addsBuffer,
    0,
    bufferAdds,
    0,
    bufferAdds.size
  );

  const commandBuffer = encoder.finish();
  device.queue.submit([commandBuffer]);

  // Read the results
  //console.time('pass end');
  await device.queue.onSubmittedWorkDone();
  await bufferDistances.mapAsync(GPUMapMode.READ);
  await bufferAdds.mapAsync(GPUMapMode.READ);
  const distances = bufferDistances.getMappedRange();
  const adds = bufferAdds.getMappedRange();
  const res = new Uint32Array(10000);
  const resAdds = new Uint32Array(10000);
  res.set(new Uint32Array(distances));
  resAdds.set(new Uint32Array(adds));
  bufferDistances.unmap();
  bufferAdds.unmap();
  atomicBuffer.destroy();
  addsBuffer.destroy();
  bufferDistances.destroy();
  bufferAdds.destroy();
  console.log(res);
  console.log('adds', resAdds);
}
