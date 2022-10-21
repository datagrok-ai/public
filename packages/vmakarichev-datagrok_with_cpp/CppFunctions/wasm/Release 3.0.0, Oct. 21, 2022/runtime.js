
// Wrapper for exported C-functions
function cFuncWrapper(module, cFuncName, returnType, argTypes, args, argsToUpdate){    

    // type-to-heap correspondence
    let heapMap = {"Int32Array": "HEAP32",
                   "Float32Array": "HEAPF32"}; // TODO: add other types
  
    let params = []; // parameters that are further substituted to C-function
    let paramTypes = []; // types of parameters that are further substituted to C-function
    let bufs = []; // buffers that are used for arrays    
    let buffer; // the current buffer
  
    // creation of parameters that are substituted to the exported C-function
    for(let i = 0; i < args.length; i++){      
  
        if(argTypes[i] == 'number') {
            params.push(args[i]);
            paramTypes.push('number');
            bufs.push(null); // i.e. no allocated memory is used
        }
  
        else { // the current argument is an array, so some buffer/heap routine shoude be performed
            buffer = module._malloc(args[i].BYTES_PER_ELEMENT * args[i].length); // allocate memory
  
            let heap = heapMap[argTypes[i]]; // choose the required heap            
            let shift = 0; // current shift
            
            // compute shift
            switch(heap) {
                case "HEAP32": case "HEAPF32": 
                    shift = 2; // this shift is specified by EMSCRIPTEN 
                    break;
                //TODO: add processing with respect to other types
            }
  
            // assign the data to the heap
            module[heap].set(args[i], buffer >> shift);
  
            // append arrays that are passed to C-function
            params.push(buffer);
            paramTypes.push('number'); // in the case of array, 'number' is also used
            bufs.push(buffer);
        }
    }
  
    // call C-function
    let result;   
    if(returnType)
        result = module.ccall(cFuncName, returnType, paramTypes, params);
    else
        module.ccall(cFuncName, returnType, paramTypes, params);
  
    // process arrays that are modified, when executing C-function
    for(let j = 0; j < argsToUpdate.length; j++){
        let m = argsToUpdate[j]; 
  
        if(argTypes[m] != 'number') {            
            let heap = heapMap[argTypes[m]]; 
            let buffer = bufs[m];
            let bytes = args[m].BYTES_PER_ELEMENT;
            
            // update the current array
            for(let n = 0; n < args[m].length; n++)
                args[m][n] = module[heap][buffer/bytes + n];                
        } 
    }
    
    // clear allocated memory
    for(let k = 0; k < bufs.length; k++)
        if(bufs[k])
            module._free(bufs[k]);
    
    if(returnType)
        return result;
  }