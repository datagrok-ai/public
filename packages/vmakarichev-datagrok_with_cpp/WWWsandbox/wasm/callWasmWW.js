// callWasmWW.js

// Runtime-system for launch wasm-functions via webworkers

import {getCppInput, foo} from './wasmWebWorkerUtils';

// The main tool that combines all together
export async function callWasmWebWorker(module, cFuncName, inputs) {        

    // get specification of exported C/C++-function
    let funcSpecification = module[cFuncName];

    // get argumnets
    //let args = funcSpecification.arguments; 

    // array of arguments that further are used by cpp-wrapper
    let cppFuncInput = getCppInput(funcSpecification.arguments, inputs); 
    
    //grok.shell.addTableView(foo());    
    
    var worker = new Worker(new URL('../wasm/workerFAE.js', import.meta.url));    
    
    worker.postMessage({'cFuncName': cFuncName, 'inputs': cppFuncInput});
  
    worker.onmessage = function(e) {  
      alert(e.data);
    } 
  } // callWasmWebWorkers
