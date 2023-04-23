import exportCppLib from './CppLib.js';
import {_package} from '../src/package.ts';
export async function testWasm (){

    const a = await exportCppLib({locateFile: () =>_package.webRoot + '/wasm/CppLib.wasm'});

    let js_wrapped_sma = a.cwrap("add", "null", ["number", "number", "number", "number"]);
    let ar1 = new Int32Array(3);
    ar1[0]=1;
    ar1[1]=2;
    ar1[2]=3;
    let ar2 = new Int32Array(3);
    ar2[0]=4;
    ar2[1]=5;
    ar2[2]=6;
    let dataPtr1 = a._malloc(12);
    let dataHeap1 = new Uint8Array(a.HEAPU8.buffer, dataPtr1, 12);
    dataHeap1.set(new Uint8Array(ar1.buffer));

    let dataPtr2 = a._malloc(12);
    let dataHeap2 = new Uint8Array(a.HEAPU8.buffer, dataPtr2, 12);
    dataHeap2.set(new Uint8Array(ar2.buffer));

    let resPointer = a._malloc(12);
    let resHeap = new Uint8Array(a.HEAPU8.buffer, resPointer, 12);

    js_wrapped_sma(dataPtr1,dataPtr2,3,resPointer);
    let result = new Int32Array(resHeap.buffer, resHeap.byteOffset, 3);
    console.log(result);


}
