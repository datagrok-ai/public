import {exportFAEexplicit} from '../wasm/FAEexplicitForWebWorker';
import {allocateMemoryForBuffer, clearMemoryForBuffer,
    putDataToBuffer, getDataFromBuffer,
    getArrOfWasmParams, getArrOfWasmTypes} from '../wasm/wasmWebWorkerUtils';

onmessage = async function (evt) {
    exportFAEexplicit().then(module => 
    {   
        let data = evt.data;          
        let cFuncName = data.cFuncName;
        let inputs = data.inputs;
        
        console.log(cFuncName);
        console.log(inputs);

        console.log('Allocating memory for buffers!');

        allocateMemoryForBuffer(module, inputs);

        putDataToBuffer(module, inputs);

        console.log(inputs);

        let params = getArrOfWasmParams(inputs);

        let types = getArrOfWasmTypes(inputs);

        console.log(params);
        console.log(types);

        // call wasm-function
        let result = module.ccall(cFuncName, 'number', types, params);        

        getDataFromBuffer(module, inputs);

        console.log(inputs);

        clearMemoryForBuffer(module, inputs);

        postMessage(result);
    } )
}