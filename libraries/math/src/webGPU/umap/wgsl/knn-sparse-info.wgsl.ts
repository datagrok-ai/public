/* eslint-disable max-len */
export function knnMatrixOpInfoWGSL(threadsPerWorkgroupDim: number, threadsPerYDimension: number, entryLen: number, knnSize = 15) {
  return `

        @group(0) @binding(0) var<storage, read_write> knnIndexes: array<array<i32, ${knnSize}>, ${entryLen}>;
        @group(0) @binding(1) var<storage, read_write> resTransposedSizes: array<atomic<u32>, ${entryLen}>;
        // this array will store the union of knn and its transpose sizes per index. should be initialized to knnSize each.
        @group(0) @binding(2) var<storage, read_write> resUnionMatrixSizes: array<atomic<u32>, ${entryLen}>;

        @compute @workgroup_size(${threadsPerWorkgroupDim}, ${threadsPerWorkgroupDim}) fn countTransposedCols(
            @builtin(global_invocation_id) id: vec3<u32>,
         ) {
            let col: u32 = id.x;
            let row: u32 = id.y;
            let workingIndex: u32 = row * ${threadsPerYDimension} + col;

            if (workingIndex >= ${entryLen}) {
                return;
            }
            
            for (var i = 0u; i < ${knnSize}u; i++) {
                let otherIndex: i32 = knnIndexes[workingIndex][i];
                if (otherIndex != -1 && otherIndex < ${entryLen}) {
                    atomicAdd(&resTransposedSizes[otherIndex], 1u);
                    atomicAdd(&resUnionMatrixSizes[otherIndex], 1u);
                    let otherIndexes = &knnIndexes[otherIndex];
                    for(var j = 0u; j < ${knnSize}; j++) {
                        if ((*otherIndexes)[j] == i32(workingIndex)) {
                            atomicSub(&resUnionMatrixSizes[workingIndex], 1u);
                            // if same is found in other array, decrement by one;
                            break;
                        }
                    }
                }

            }
         }
    `;
}
