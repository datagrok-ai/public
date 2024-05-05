import {getGPUDevice} from '../getGPUDevice';

/* eslint-disable max-len */

/**
 * wgsl code for pairwise operation on two sparse matrices.
 * @param indexes1 - indexes of first sparse matrix
 * @param distances1 - values1 of first sparse matrix
 * @param offsets1 - 32 bit offsets of first sparse matrix
 * @param indexes2 - indexes of second sparse matrix
 * @param distances2 - values of second sparse matrix
 * @param offsets2 - 32 bit offsets of second sparse matrix
 * @param numOfEntries - number of entries in original data
 * @param unionMatrixOffsets - offsets of union matrix, used to store result
 * @param op - operation to perform
 * @returns
 */
export async function pairwiseOpSparse(
  indexes1: number[] | Int32Array | Uint32Array,
  distances1: number[] | Float32Array | Uint32Array,
  offsets1: number[] | Int32Array | Uint32Array,
  indexes2: number[] | Int32Array | Uint32Array,
  distances2: number[] | Float32Array | Uint32Array,
  offsets2: number[] | Int32Array | Uint32Array,
  numOfEntries: number,
  unionMatrixOffsets: number[] | Int32Array | Uint32Array, // array of size numOfEntries + 1
  op: '+' | '*' | '-' | '/'
) {
  // this function is only meant for regular knn matrix, not for sparse matrix for of knn.
  const threadsPerWorkgroupDim = 10;
  const threadsPerWorkgroup = threadsPerWorkgroupDim * threadsPerWorkgroupDim;
  const requiredWorkgroups = Math.ceil(numOfEntries / threadsPerWorkgroup);
  const numOfWorkgroupsPerDim = Math.ceil(Math.sqrt(requiredWorkgroups));
  const threadsPerYDimension = threadsPerWorkgroupDim * numOfWorkgroupsPerDim;
  const sparseSize1 = offsets1[offsets1.length - 1];
  const sparseSize2 = offsets2[offsets2.length - 1];

  // this array will have offsets of each row in union matrix with length of numOfEntries + 1

  const unionMatrixSize = unionMatrixOffsets[unionMatrixOffsets.length - 1]; // last index will be last offset i.e. full size
  const processWgsl = `
        struct SparseKNNStorage1 {
            indexes:array<i32, ${sparseSize1}>,
            knnDistances: array<f32, ${sparseSize1}>,
            offsets: array<u32, ${numOfEntries + 1}>
        }

        struct SparseKNNStorage2 {
            indexes:array<i32, ${sparseSize2}>,
            knnDistances: array<f32, ${sparseSize2}>,
            offsets: array<u32, ${numOfEntries + 1}>
        }

        struct ResKnn {
            knnIndexes: array<i32, ${unionMatrixSize}>,
            knnDistances: array<f32, ${unionMatrixSize}>
        }

        @group(0) @binding(0) var<storage, read_write> source1: SparseKNNStorage1;
        @group(0) @binding(1) var<storage, read_write> source2: SparseKNNStorage2;
        @group(0) @binding(2) var<storage, read_write> result: ResKnn;
        @group(0) @binding(3) var<storage, read_write> unionMatrixOffsets: array<u32, ${numOfEntries + 1}>;
        @compute @workgroup_size(${threadsPerWorkgroupDim}, ${threadsPerWorkgroupDim}) fn pairwiseOp(
            @builtin(global_invocation_id) id: vec3<u32>,
        ) {

            let col: u32 = id.x;
            let row: u32 = id.y;
            let workingIndex: u32 = row * ${threadsPerYDimension} + col;

            if (workingIndex >= ${numOfEntries}) {
                return;
            }
            var curUnionOffset: u32 = unionMatrixOffsets[workingIndex]; // offset at which to start writing in union matrix

            let start1 = source1.offsets[workingIndex];
            let end1 = source1.offsets[workingIndex + 1];
            let start2 = source2.offsets[workingIndex];
            let end2 = source2.offsets[workingIndex + 1];
            // TODO: sort a copy of these arrays, will make operation faster
            for (var i: u32 = start1; i < end1; i++) {
                let val1: f32 = source1.knnDistances[i];
                var val2: f32 = 0.0;
                let curIndex: i32 = source1.indexes[i];
                for (var j: u32 = start2; j < end2; j++) {
                    if (source2.indexes[j] == curIndex) {
                        val2 = source2.knnDistances[j];
                        break;
                    }
                }
                
                result.knnIndexes[curUnionOffset] = curIndex;
                result.knnDistances[curUnionOffset] = val1 ${op} val2;
                curUnionOffset += 1;
            }

            if (curUnionOffset >= unionMatrixOffsets[workingIndex + 1]) {
                // small optimization applicable only to our case, but saves a good amount of time
                return;
            }

            // do the same for source2 but skip the indexes where they are already present in source1
            for (var i: u32 = start2; i < end2; i++) {
                let val1: f32 = source2.knnDistances[i];
                let curIndex: i32 = source2.indexes[i];
                var found = false;
                for (var j: u32 = start1; j < end1; j++) {
                    if (source1.indexes[j] == curIndex) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.knnIndexes[curUnionOffset] = curIndex;
                    result.knnDistances[curUnionOffset] = val1 ${op} 0.0;
                    curUnionOffset += 1;
                }
            }
        }
`;


  const storage32Size1 = sparseSize1 * 2 + numOfEntries + 1;
  const storage32Size2 = sparseSize2 * 2 + numOfEntries + 1;
  const resStorage32Size = unionMatrixSize * 2;
  // standard stuff to pad the storage size to multiple of 16 bytes
  let paddedStorage1ByteSize = storage32Size1 * 4;
  const remainder1 = paddedStorage1ByteSize & 15;
  if (remainder1 !== 0)
    paddedStorage1ByteSize += 16 - remainder1;


  let paddedStorage2ByteSize = storage32Size2 * 4;
  const remainder2 = paddedStorage2ByteSize & 15;
  if (remainder2 !== 0)
    paddedStorage2ByteSize += 16 - remainder2;


  let paddedResStorageByteSize = resStorage32Size * 4;
  const remainder3 = paddedResStorageByteSize & 15;
  if (remainder3 !== 0)
    paddedResStorageByteSize += 16 - remainder3;


  const device = await getGPUDevice();
  if (device == null)
    return;


  const module = device.createShaderModule({
    label: 'pairwiseOpShader',
    code: processWgsl
  });


  const pipeline = device.createComputePipeline({
    label: 'pairwiseOpPipeline',
    layout: 'auto',
    compute: {
      module,
      entryPoint: 'pairwiseOp',
    },
  });


  // input params
  const source1Buffer = device.createBuffer({
    label: 'source 1 buffer',
    size: paddedStorage1ByteSize,
    usage:
        GPUBufferUsage.STORAGE |
        GPUBufferUsage.COPY_SRC |
        GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedSource1Buffer = source1Buffer.getMappedRange();

  const source1IndexesArray = new Int32Array(mappedSource1Buffer, 0, sparseSize1);
  source1IndexesArray.set(indexes1);
  const source1DistancesArray = new Float32Array(mappedSource1Buffer, sparseSize1 * 4, sparseSize1);
  source1DistancesArray.set(distances1);
  const source1OffsetsArray = new Uint32Array(mappedSource1Buffer, sparseSize1 * 8, numOfEntries + 1);
  source1OffsetsArray.set(offsets1);
  source1Buffer.unmap();

  const source2Buffer = device.createBuffer({
    label: 'source 2 buffer',
    size: paddedStorage2ByteSize,
    usage:
        GPUBufferUsage.STORAGE |
        GPUBufferUsage.COPY_SRC |
        GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });
  const mappedSource2Buffer = source2Buffer.getMappedRange();

  const source2IndexesArray = new Int32Array(mappedSource2Buffer, 0, sparseSize2);
  source2IndexesArray.set(indexes2);
  const source2DistancesArray = new Float32Array(mappedSource2Buffer, sparseSize2 * 4, sparseSize2);
  source2DistancesArray.set(distances2);
  const source2OffsetsArray = new Uint32Array(mappedSource2Buffer, sparseSize2 * 8, numOfEntries + 1);
  source2OffsetsArray.set(offsets2);
  source2Buffer.unmap();

  const resBuffer = device.createBuffer({
    label: 'res buffer',
    size: paddedResStorageByteSize,
    usage:
        GPUBufferUsage.STORAGE |
        GPUBufferUsage.COPY_SRC |
        GPUBufferUsage.COPY_DST,
  });

  const unionMatrixOffsetsBuffer = device.createBuffer({
    label: 'union matrix offsets buffer',
    size: (numOfEntries + 1) * 4,
    usage:
        GPUBufferUsage.STORAGE |
        GPUBufferUsage.COPY_SRC |
        GPUBufferUsage.COPY_DST,
    mappedAtCreation: true,
  });

  const mappedUnionMatrixOffsetsBuffer = unionMatrixOffsetsBuffer.getMappedRange();
  const unionMatrixOffsetsArray = new Uint32Array(mappedUnionMatrixOffsetsBuffer);
  unionMatrixOffsetsArray.set(unionMatrixOffsets);
  unionMatrixOffsetsBuffer.unmap();

  const bindGroup = device.createBindGroup({
    label: 'bindGroup for pairwise ops',
    layout: pipeline.getBindGroupLayout(0),
    entries: [
      {binding: 0, resource: {buffer: source1Buffer}},
      {binding: 1, resource: {buffer: source2Buffer}},
      {binding: 2, resource: {buffer: resBuffer}},
      {binding: 3, resource: {buffer: unionMatrixOffsetsBuffer}}
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

  const resultBuffer = device.createBuffer({
    label: 'result buffer',
    size: paddedResStorageByteSize,
    usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
  });

  encoder.copyBufferToBuffer(
    resBuffer,
    0,
    resultBuffer,
    0,
    resultBuffer.size
  );

  const commandBuffer = encoder.finish();
  device.queue.submit([commandBuffer]);
  await device.queue.onSubmittedWorkDone();

  await resultBuffer.mapAsync(GPUMapMode.READ);
  const resultArrayBuffer = resultBuffer.getMappedRange();

  const resultKnnIndexes = new Int32Array(unionMatrixSize);
  const resultKnnDistances = new Float32Array(unionMatrixSize);

  // its actually a sparse matrix.
  resultKnnIndexes.set(new Int32Array(resultArrayBuffer, 0, unionMatrixSize));
  resultKnnDistances.set(new Float32Array(resultArrayBuffer, unionMatrixSize * 4, unionMatrixSize));
  resultBuffer.unmap();

  resultBuffer.destroy();
  source1Buffer.destroy();
  source2Buffer.destroy();
  resBuffer.destroy();
  unionMatrixOffsetsBuffer.destroy();

  return {resultKnnIndexes, resultKnnDistances};
}
