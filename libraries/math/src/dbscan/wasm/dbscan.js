import { exportCppDbscanLib } from "./wasmDbscan.js";

export async function dbscan(embedX, embedY, epsilon, minPts) {
    const wasmUrl = new URL('./wasmDbscan.wasm', import.meta.url).href;
    const wasmPath = wasmUrl.substring(0, wasmUrl.lastIndexOf('/') + 1) + 'wasmDbscan.wasm';
    let wasmInstance;
    try {
        wasmInstance = await exportCppDbscanLib({locateFile: () => wasmUrl, printErr: (_) => {}});
    } catch (e) {
        try {
            wasmInstance = await exportCppDbscanLib({locateFile: () => wasmPath, printErr: (_) => {}});
        }
        catch (e) {
            console.error(e);
            throw new Error('Unable to load wasm file for dbscan');
        }
    }
    const getDbscanWasm =
    wasmInstance.cwrap('dbscan', 'null', ['number', 'number', 'number', 'number', 'number', 'number']);
    const bytesForEmbedX = Float32Array.BYTES_PER_ELEMENT * embedX.length;
    const dataPtrX = wasmInstance._malloc(bytesForEmbedX); // allocate a memory block on the wasm heap for the embedX
    const dataPtrY = wasmInstance._malloc(bytesForEmbedX); // allocate a memory block on the wasm heap for the embedY
    const xHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, dataPtrX, bytesForEmbedX);
    xHeap.set(new Uint8Array(embedX.buffer)); // copy embedX to the wasm heap
    const yHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, dataPtrY, bytesForEmbedX);
    yHeap.set(new Uint8Array(embedY.buffer)); // copy embedY to the wasm heap
    
    const clusterPtr = wasmInstance._malloc(bytesForEmbedX);
    const clusterHeap = new Uint8Array(wasmInstance.HEAPU8.buffer, clusterPtr, bytesForEmbedX);

    getDbscanWasm(dataPtrX, dataPtrY, embedX.length, epsilon, minPts, clusterPtr);
    const clusterResult = new Int32Array(clusterHeap.buffer, clusterHeap.byteOffset, embedX.length);
    wasmInstance._free(xHeap.byteOffset);
    wasmInstance._free(yHeap.byteOffset);
    wasmInstance._free(clusterHeap.byteOffset);
    return clusterResult;
}