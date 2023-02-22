// Utilities for calling wasm-functions via webworker.

// Constants for wasm-functions in webworkers runtime system
const TYPE = 'type';
const NUM_TYPE = 'num';
const FLOAT_COLUMN_TYPE = 'floatColumn';
const INT_COLUMN_TYPE = 'intColumn';
const NEW_FLOAT_COLUMNS_TYPE = 'newFloatColumns';
const CALL_RESULT = '_callResult';
const NUM_OF_ROWS = 'numOfRows';
const NUM_OF_COLUMNS = 'numOfColumns';
const REF = 'ref';
const VALUE = 'value';
const TABLE_OF_COLUMNS = 'tableFromColumns';
const INT_TYPE = 'int';
const DOUBLE_TYPE = 'double';

// type-to-heap correspondence
const heapMap = {
                 'intColumn': "HEAP32", // Int32Array
                 'floatColumn': "HEAPF32", // Float32Array
                 'newFloatColumns': "HEAPF32" // Float32Array
                };

// type signature to typed array ,app
const typeMap = {
                 'intColumn': Int32Array,
                 'floatColumn': Float32Array, 
                 'newFloatColumns': Float32Array // Float32Array
                }; 

// type-to-shift map
const shiftMap = {'intColumn': 2, // Int32Array
                  'floatColumn': 2 // Float32Array
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
            case INT_TYPE:
            case DOUBLE_TYPE:
                arg.data = inputVals[i];
                i++;
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
                arg.data = { 'array': inputVals[i].getRawData() };
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

    console.log(cppFuncInput);

    return cppFuncInput;
}

// Allocate memory for buffers for array data
function allocateMemoryForBuffer(module, inputs) {
    for(let arg of inputs) 
        switch(arg.type) {

            case NUM_TYPE: // no memory allocation for numbers
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
                arg.data.buf = module._malloc(arg.data.array.length 
                    * arg.data.array.BYTES_PER_ELEMENT);
                break;

            case NEW_FLOAT_COLUMNS_TYPE: // allocation memory for columns that are created
                arg.data.buf = module._malloc(arg.data.numOfRows * arg.data.numOfColumns 
                    * typeMap[arg.type].BYTES_PER_ELEMENT);
                break;

            // TODO: process other cases and mistakes
            default:
                break; 
        }  
    
    console.log(inputs);
}

// Get array of values that are put to wasm-function
function getArrOfWasmParams(inputs) {
    let params = [];

    for(let arg of inputs) {
        switch (arg.type) {

            case NUM_TYPE: // put number            
            case INT_TYPE:
            case DOUBLE_TYPE:
                params.push(arg.data);
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
                params.push(arg.data.buf);                
                params.push(arg.data.array.length);
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
function getArrOfWasmTypes(inputs) {
    let types = [];

    for(let arg of inputs) {
        switch (arg.type) {

            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:
                types.push('number');
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
                types.push('number');                
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
function putDataToBuffer(module, inputs) {
    // TODO: implement this
    for(let arg of inputs) {
        let type = arg.type;

        switch (type) {

            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:                
                break;
                
            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:                
                let shift = shiftMap[type];
                let heap = module[heapMap[type]];  
                heap.set(arg.data.array, arg.data.buf >> shift);
                break;
            
            case NEW_FLOAT_COLUMNS_TYPE:                
                break;

            // TODO: process other cases and mistakes
            default:
                break;
        }
    }
}

// Get array data from buffer
function getDataFromBuffer(module, inputs) {

    for(let arg of inputs) {
        let type = arg.type;
        switch (type) { 
            
            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:            
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
function clearMemoryForBuffer(module, inputs) {
    for(let arg of inputs) 
        switch(arg.type) {

            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
                module._free(arg.data.buf);
                break;
                
            case NEW_FLOAT_COLUMNS_TYPE:
                module._free(arg.data.buf);
                break;

            // TODO: process other cases and mistakes
            default:
                break; 
        }    
}

// Extract newly created data: new column(s) are created
function extractNewlyCreatedData(funcSpecificationArgs, argsAfterWasmCall) {
    // type-to-column_creator map
    const typeToColumnCreatorMap = {'newFloatColumns': DG.Column.fromFloat32Array};

    let i = 0;

    for(let key in funcSpecificationArgs) 
    {
        let arg = funcSpecificationArgs[key];

        switch(arg.type){

            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
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
function getOutput(funcSpecification) {  
    let output = funcSpecification.output;    

    switch(output.type) {

        case NUM_TYPE:
        case INT_TYPE:
        case DOUBLE_TYPE:
            return funcSpecification.arguments[output.source];
            break;

        case TABLE_OF_COLUMNS:
            return DG.DataFrame.fromColumns(funcSpecification.arguments[output.source].columns);
            break;

        // TODO: process other cases and mistakes
        default:
            break;
    }

}

export function cppWrapper(module, args, cppFuncName, returnType)
{
    // allocate memory for arrays that are passed to C++-function
    allocateMemoryForBuffer(module, args);

    // put data (just column(s)) to allocated buffers
    putDataToBuffer(module, args);

    // create array of parameters that are passed to C++-function
    let params = getArrOfWasmParams(args);

    console.log(params);

    // create array of parameters' types that are passed to C++-function
    let types = getArrOfWasmTypes(args);
    
    console.log(types);

    // call wasm-function
    let result = module.ccall(cppFuncName, returnType, types, params);
    
    console.log(result);

    // get data from buffers (just column(s))
    getDataFromBuffer(module, args);        

    // clear memory that was previousely allocated
    clearMemoryForBuffer(module, args);

    console.log('done');

    return result;
}

export function getResult(funcSpecification, dataFromWebWorker)
{
    funcSpecification.arguments._callResult = dataFromWebWorker.callResult;  

    extractNewlyCreatedData(funcSpecification.arguments, dataFromWebWorker.args);

    return getOutput(funcSpecification);
}