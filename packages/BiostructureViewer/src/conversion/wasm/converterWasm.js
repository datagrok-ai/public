import { convertCifToPdb } from './wasmConvert.js';

export async function addWasm(inputString) {
    // Convert input string to Uint8Array
    const inputArray = new TextEncoder().encode(inputString);

    const wasmUrl = new URL('./wasmConvert.wasm', import.meta.url).href;
    const wasmPath = wasmUrl.substring(0, wasmUrl.lastIndexOf('/') + 1) + 'wasmConvert.wasm';
    let wasmInstance;
    try {
      wasmInstance = await convertCifToPdb({ locateFile: () => wasmPath, printErr: (_) => {} });
    } catch (e) {
      try {
        wasmInstance = await convertCifToPdb({ locateFile: () => wasmUrl, printErr: (_) => {} });
      } catch (e) {
        console.error(e);
        throw new Error('Unable to load wasm file');
      }
    }

    const convertFunction = wasmInstance.cwrap('convert', 'number', ['number', 'number']);

    // Allocate memory for input array in the WASM heap
    const inputPointer = wasmInstance._malloc(inputArray.length);
    wasmInstance.HEAPU8.set(inputArray, inputPointer);

    // Call the WASM function
    const resultPointer = convertFunction(inputPointer, inputArray.length);

    // Create a copy of the result array from the WASM heap
    const resultArray = new Uint8Array(inputArray.length);
    resultArray.set(wasmInstance.HEAPU8.subarray(resultPointer, resultPointer + inputArray.length));

    // Free the allocated memory
    wasmInstance._free(inputPointer);
    wasmInstance._free(resultPointer);

    // Convert result Uint8Array back to string
    const resultString = new TextDecoder().decode(resultArray);

    return resultString;
}
