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
export const TABLE_OF_COLUMNS = 'tableFromColumns';

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

            case NUM_TYPE: // no memory allocation for numbers
                break;

            case NEW_FLOAT_COLUMNS_TYPE: // allocation memory for columns that are created
                arg.data.buf = module._malloc(arg.data.numOfRows * arg.data.numOfColumns 
                    * typeMap[arg.type].BYTES_PER_ELEMENT);
                break;

            // TODO: process other cases and mistakes
            default:
                break; 
        }    
}

// Get array of values that are put to wasm-function
export function getArrOfWasmParams(inputs) {
    let params = [];

    for(let arg of inputs) {
        switch (arg.type) {

            case NUM_TYPE: // put number
                params.push(arg.data);
                break;
            
            case NEW_FLOAT_COLUMNS_TYPE: // put buffer and sizes
                params.push(arg.data.buf);
                params.push(arg.data.numOfRows);
                params.push(arg.data.numOfColumns);
                break;
        
            // TODO: process other cases and mistakes
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

            // TODO: process other cases and mistakes
            default:
                break;
        }
    }

    return types;
}

// Put array data to buffer
export function putDataToBuffer(module, inputs) {
    // TODO: implement this
}

// Get array data from buffer
export function getDataFromBuffer(module, inputs) {

    for(let arg of inputs) {
        let type = arg.type;
        switch (type) { 
            
            case NUM_TYPE:
                break;

            // get arrays that will be further used for newly created columns
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
        
            // TODO: process other cases and mistakes
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

            // TODO: process other cases and mistakes
            default:
                break; 
        }    
}

// Extract newly created data: new column(s) are created
export function extractNewlyCreatedData(funcSpecificationArgs, argsAfterWasmCall) {
    // type-to-column_creator map
    const typeToColumnCreatorMap = {'newFloatColumns': DG.Column.fromFloat32Array};

    let i = 0;

    for(let key in funcSpecificationArgs) 
    {
        let arg = funcSpecificationArgs[key];

        switch(arg.type){

            case NUM_TYPE:
                break;

            case NEW_FLOAT_COLUMNS_TYPE:
                let columns = [];

                for(let j = 0; j < argsAfterWasmCall[i].arrays.length; j++)
                    columns.push(typeToColumnCreatorMap[arg.type](arg.names[j],
                        argsAfterWasmCall[i].arrays[j]));

                arg.columns = columns;
                
                break;

            // TODO: process other cases and mistakes
            default:
                break;
        }

        i++;
    }
}

// Get output data: overall output is created
export function getOutput(funcSpecification) {  
    let output = funcSpecification.output;    

    switch(output.type) {
        case TABLE_OF_COLUMNS:
            return DG.DataFrame.fromColumns(funcSpecification.arguments[output.source].columns);
            break;

        // TODO: process other cases and mistakes
        default:
            break;
    }

}