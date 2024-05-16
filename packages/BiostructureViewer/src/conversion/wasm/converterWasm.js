import { convertCifToPdb } from './wasmConvert.js';

export async function addWasm(a) {
    const wasmUrl = new URL('./wasmCluster.wasm', import.meta.url).href;
    const wasmPath = wasmUrl.substring(0, wasmUrl.lastIndexOf('/') + 1) + 'wasmCluster.wasm';
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
  
    // Wrap the cwrap call in a non-strict mode function
    const wrapFunction = () => {
      const convertFunction = wasmInstance.cwrap('convert', 'string', ['string']);
      const result = convertFunction(a);
      return result;
    };
  
    // Call the wrapper function
    const result = wrapFunction();
  
    return result;
  }
  
