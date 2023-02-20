// Utilities for calling wasm-functions via webworker.

// Constants for wasm-functions in webworkers runtime system
export const TYPE = 'type';
export const NUM_TYPE = 'num';
export const NEW_FLOAT_COLUMNS_TYPE = 'newFloatColumns';
export const CALL_RESULT = '_callResult';
export const NUM_OF_ROWS = 'numOfRows';
export const NUM_OF_COLUMNS = 'numOfColumns';
export const REF = 'ref';
export const VALUE = 'value';

// type-to-heap correspondence
const heapMap = { 
                 'newFloatColumns': "HEAPF32" // Float32Array
                };

// type signature to typed array ,app
const typeMap = { 
                 'newFloatColumns': Float32Array // Float32Array
                }; 
 

// Get input for C++-function
export function getCppInput(argsSpecification, inputVals)
{
    let cppFuncInput = [];

    // complete an input for cpp
    let i = 0;
    for(let key in argsSpecification){
        let arg = argsSpecification[key]; 
        
        // skip auxiliry element
        if(key == CALL_RESULT)             
            continue;

        // create an input
        switch(arg.type){
            case NUM_TYPE:
                arg.data = inputVals[i];
                i++;
                break;
            case NEW_FLOAT_COLUMNS_TYPE:
                let val1 = argsSpecification[ arg[NUM_OF_ROWS][REF] ].data;
                let val2 = argsSpecification[ arg[NUM_OF_COLUMNS][REF] ].data;                
                arg.data = {'numOfRows': val1,
                  'numOfColumns': val2};
                i++;
                break;
            // TODO: add other cases
            default: 
                return; // TODO: specify behaviour           
            } // switch        

        cppFuncInput.push(arg);
    } // for key

    return cppFuncInput;
}

// Allocate memory for buffers for array data
export function allocateMemoryForBuffer(module, inputs) {
    for(let arg of inputs) 
        switch(arg.type) {
            case NUM_TYPE:
                break;
            case NEW_FLOAT_COLUMNS_TYPE:
                arg.data.buf = module._malloc(arg.data.numOfRows * arg.data.numOfColumns 
                    * Float32Array.BYTES_PER_ELEMENT);
                break;
            default:
                break; // TODO: process other cases and mistakes
        }    
}

// Get array of values that are put to wasm-function
export function getArrOfWasmParams(inputs) {
    let params = [];

    for(let arg of inputs) {
        switch (arg.type) {
            case NUM_TYPE:
                params.push(arg.data);
                break;
            
            case NEW_FLOAT_COLUMNS_TYPE:
                params.push(arg.data.buf);
                params.push(arg.data.numOfRows);
                params.push(arg.data.numOfColumns);
                break;
        
            default:
                break;
        }
    }

    return params;
}

// Get array of types that are put to wasm-function
export function getArrOfWasmTypes(inputs) {
    let types = [];

    for(let arg of inputs) {
        switch (arg.type) {
            case NUM_TYPE:
                types.push('number');
                break;
            
            case NEW_FLOAT_COLUMNS_TYPE:
                types.push('number');
                types.push('number');
                types.push('number');
                break;
        
            default:
                break;
        }
    }

    return types;
}

// Put array data to buffer
export function putDataToBuffer(module, inputs) {
}

// Get array data from buffer
export function getDataFromBuffer(module, inputs) {

    for(let arg of inputs) {
        let type = arg.type;
        switch (type) {          
            case NEW_FLOAT_COLUMNS_TYPE:
                let heap = module[heapMap[type]];
                let numOfRows = arg.data.numOfRows;
                let numOfCols = arg.data.numOfColumns;
                let numOfBytes = typeMap[type].BYTES_PER_ELEMENT;    
                let buf = arg.data.buf;
                let arrays = [];

                for(let i = 0; i < numOfCols; i++) {
                   let arr = new typeMap[type](numOfRows);

                   for(let j = 0; j < numOfRows; j++)
                      arr[j] = heap[buf / numOfBytes + j + i * numOfRows];               
                
                   arrays.push(arr);
                }

                arg.arrays = arrays;
                
                break;
        
            default:
                break;
        }
    }            
}

// Clear memory allocated for array data
export function clearMemoryForBuffer(module, inputs) {
    for(let arg of inputs) 
        switch(arg.type) {
            case NUM_TYPE:
                break;
            case NEW_FLOAT_COLUMNS_TYPE:
                module._free(arg.buf);
                break;
            default:
                break; // TODO: process other cases and mistakes
        }    
}

export function foo() {
    let t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
    t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
    return t;
}