import {exportCppLib} from './wasmCluster.js';

/**
 calls wasm function to get clustering merges and heights from distance matrix
 distmat: distance matrix as an array of floats comprising the top right triangle of the matrix
 n: number of observables
 method: 0 = single linkage, 1 = complete linkage, 2 = average linkage, 3 = median linkage
 result: { mergeResult: Int32Array, heightsResult: Float32Array }
 mergeResult has length of (n-1)*2 and heightResult has length of n-1
 together they describe the clustering as a dendrogram as a matrix of (n-1) X 3
 each row in the given matrix comprises the merge event at step i and last column is the height from the leaves.
*/
export async function getClustersFromDistMatWasm(distmat, n, method) {
  const wasmUrl = new URL('./wasmCluster.wasm', import.meta.url).href;
  const wasmPath = wasmUrl.substring(0, wasmUrl.lastIndexOf('/') + 1) + 'wasmCluster.wasm';
  let wasmInstance;
  try {
      wasmInstance = await exportCppLib({locateFile: () => wasmPath, printErr: (_) => {}});
  } catch (e) {
      try {
          wasmInstance = await exportCppLib({locateFile: () => wasmUrl, printErr: (_) => {}});
      }
      catch (e) {
          console.error(e);
          throw new Error('Unable to load wasm file for hierarchical clustering');
      }
  }

  const getDendrogramWasm =
    wasmInstance.cwrap('getDendrogram', 'null', ['number', 'number', 'number', 'number', 'number']);
  // allocate a memory block on the wasm heap for the distance matrix
  const bytesForDistMat = Float32Array.BYTES_PER_ELEMENT * distmat.length;
  const dataPtr = wasmInstance._malloc(bytesForDistMat);
  const matHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, dataPtr, bytesForDistMat);
  matHeap.set(new Uint8Array(distmat.buffer));

  //allocate a memory block on the wasm heap for the merges array
  const bytesForMergeArray = Int32Array.BYTES_PER_ELEMENT * (n - 1) * 2;
  const mergePtr = wasmInstance._malloc(bytesForMergeArray);
  const mergeHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, mergePtr, bytesForMergeArray);

  //allocate a memory block on the wasm heap for the heigt array
  const bytesForHeightsArray = Float32Array.BYTES_PER_ELEMENT * (n - 1);
  const heightsPtr = wasmInstance._malloc(bytesForHeightsArray);
  const heightsHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, heightsPtr, bytesForHeightsArray);

  // call the wasm function
  getDendrogramWasm(dataPtr, n, method, mergePtr, heightsPtr);
  // get the merges and heights from the wasm heap
  const mergeResult = new Int32Array(mergeHeap.buffer, mergeHeap.byteOffset, (n-1) * 2);
  const heightsResult = (new Float32Array(heightsHeap.buffer, heightsHeap.byteOffset, n-1)).slice(0, n-1);

  const mergeRow1 = mergeResult.slice(0, n-1);
  const mergeRow2 = mergeResult.slice(n-1, 2*n - 2);
  // free the memory blocks on the wasm heap
  wasmInstance._free(matHeap.byteOffset);
  wasmInstance._free(mergeHeap.byteOffset);
  wasmInstance._free(heightsHeap.byteOffset);
  return {mergeRow1, mergeRow2, heightsResult};
}
