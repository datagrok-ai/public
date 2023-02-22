// Runtime-system for launch wasm-functions via webworkers
// REMARK! Currently, several parts are hard-coded. This is a basis.

import {getCppInput, extractNewlyCreatedData, getOutput} from './wasmUtils';

// The main tool that combines all together
export async function callWasmWebWorker(module, cFuncName, inputs) {        

    // get specification of exported C/C++-function
    let funcSpecification = module[cFuncName];

    // array of arguments that further are used by cpp-wrapper
    let cppFuncInput = getCppInput(funcSpecification.arguments, inputs);        

    // name of js-worker is hard-coded 
    // TODO: provide automatic naming!
    var worker = new Worker(new URL('../wasm/solveFAEexplicitWorker.js', import.meta.url)); 
        
    // run wasm-function in webworker
    worker.postMessage({'cFuncName': cFuncName, 'cppFuncInputs': cppFuncInput});
  
    // process results of wasm-function execution
    worker.onmessage = function(e) {   
        
        // get result, which is returned by wasm-function
        funcSpecification.arguments._callResult = e.data.callResult;

        // extract data (for example, new columns), created in wasm-function
        extractNewlyCreatedData(funcSpecification.arguments, e.data.args);

        // get output, i.e. a goal of wasm-computation (for example, new columns)
        let output = getOutput(funcSpecification);

        // THE FOLLOWING IS HARD-CODED
        // TODO: provide automatization or RichFunctionView
        output.name = 'Solution';
        let view = grok.shell.addTableView(output);
        view.lineChart({ markerType: 'dot', sharex: 'true', multiAxis: 'true'});        
    } 
  } // callWasmWebWorkers
