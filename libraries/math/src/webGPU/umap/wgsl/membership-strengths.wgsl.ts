/* eslint-disable max-len */
export function computeMembershipStrengthsWGSL(
  threadsPerWorkgroupDim: number = 10, threadsPerYDimension: number,
  knnDistances: number[][] | Float32Array[],
  sigmas: number[] | Float32Array, rhos: number[] | Float32Array
) {
  const numOfEntries = knnDistances.length;
  const sigmasLength = sigmas.length;
  const rhosLength = rhos.length;
  // TODO: do it with structs
  return `
           struct MembershipStrengthsInfo {
               knnDistances: array<array<f32, ${knnDistances[0].length}>, ${numOfEntries}>,
               sigmas: array<f32, ${sigmasLength}>,
               rhos: array<f32, ${rhosLength}>
           };
           var<private> nNeighbors: u32 = ${knnDistances[0].length};
           @group(0) @binding(0) var<storage, read> membershipStrengthsInfo: MembershipStrengthsInfo;
           @group(0) @binding(1) var<storage, read_write> resultKnnDistances: array<array<f32, ${knnDistances[0].length}>, ${numOfEntries}>;
           @compute @workgroup_size(${threadsPerWorkgroupDim}, ${threadsPerWorkgroupDim}) fn computeMembershipStrengths(
               @builtin(global_invocation_id) id: vec3<u32>,
           ) {
               let col: u32 = id.x;
               let row: u32 = id.y;
               let workingIndex: u32 = row * ${threadsPerYDimension} + col;

               if (workingIndex >= ${numOfEntries}) {
                   return;
               }

               let workingResKnnDistances = &resultKnnDistances[workingIndex];
               let knnDistances = &membershipStrengthsInfo.knnDistances[workingIndex];
               for (var j = 0u; j < nNeighbors; j = j + 1u) {
                   var val: f32 = 0.0;

                   if ((*knnDistances)[j] - membershipStrengthsInfo.rhos[workingIndex] <= 0.0) {
                       val = 1.0;
                   } else {
                       val = exp(0.0 - ((knnDistances[j] - membershipStrengthsInfo.rhos[workingIndex]) / membershipStrengthsInfo.sigmas[workingIndex]));   
                   }

                   (*workingResKnnDistances)[j] = val;
               }
           }
   `;
}
