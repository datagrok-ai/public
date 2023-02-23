// Utilities for calling wasm-functions via webworker.

// Constants for wasm-functions in webworkers runtime system
const TYPE = 'type';
const NUM_TYPE = 'num';
const FLOAT_COLUMN_TYPE = 'floatColumn';
const INT_COLUMN_TYPE = 'intColumn';
const FLOAT_COLUMNS_TYPE = 'floatColumns';
const NEW_FLOAT_COLUMNS_TYPE = 'newFloatColumns';
const INT_COLUMNS_TYPE = 'intColumns';
const NEW_INT_COLUMNS_TYPE = 'newIntColumns';
const NEW_FLOAT_COLUMN_TYPE = 'newFloatColumn';
const NEW_INT_COLUMN_TYPE = 'newIntColumn';
const COLUMN = 'column';
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
                 'intColumn': "HEAP32",
                 'floatColumn': "HEAPF32",
                 'floatColumns': "HEAPF32",
                 'newFloatColumns': "HEAPF32",
                 'intColumns': "HEAP32",
                 'newIntColumns': "HEAP32",
                 'newFloatColumn': "HEAPF32",
                 'newIntColumn': "HEAP32"
                };

// type signature to typed array ,app
const typeMap = {
                 'intColumn': Int32Array,
                 'floatColumn': Float32Array,
                 'floatColumns': Float32Array,
                 'newFloatColumns': Float32Array,
                 'intColumns': Int32Array,  
                 'newIntColumns': Int32Array,
                 'newFloatColumn': Float32Array,
                 'newIntColumn': Int32Array
                }; 

// type-to-shift map
const shiftMap = {'intColumn': 2, 
                  'floatColumn': 2, 
                  'floatColumns': 2,  
                  'newFloatColumns': 2, 
                  'intColumns': 2,  
                  'newIntColumns': 2,
                  'newFloatColumn': 2,
                  'newIntColumn': 2
                 };         

// Get input for C++-function
export function getCppInput(argsSpecification, inputVals)
{
    let cppFuncInput = [];
    let ref;

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
                let array = inputVals[i].getRawData();
                arg.data = { 'array': array,
                             'numOfRows': array.length};
                i++;
                break;

            case NEW_INT_COLUMN_TYPE:
            case NEW_FLOAT_COLUMN_TYPE:
                let val = 0;                

                ref = arg[NUM_OF_ROWS][REF];

                if (ref.type == INT_TYPE)
                    val = argsSpecification[ref].data;
                else
                    val = argsSpecification[ref].data[arg[NUM_OF_ROWS][VALUE]];                                

                arg.data = {'numOfRows': val};

                i++;
                break;

            case INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:                
                let arrays = [];                

                for(let col of inputVals[i].toList())
                  arrays.push(col.getRawData());

                arg.data = { 'arrays': arrays,
                    'numOfRows': arrays[0].length,
                    'numOfColumns': arrays.length};

                i++;  
                break;

            case NEW_INT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMNS_TYPE:
                let val1 = 0;
                let val2 = 0;

                ref = arg[NUM_OF_ROWS][REF];

                if (ref.type == INT_TYPE)
                    val1 = argsSpecification[ref].data;
                else
                    val1 = argsSpecification[ref].data[arg[NUM_OF_ROWS][VALUE]];

                ref = arg[NUM_OF_COLUMNS][REF];

                if (ref.type == INT_TYPE)
                    val2 = argsSpecification[ref].data;
                else
                    val2 = argsSpecification[ref].data[arg[NUM_OF_COLUMNS][VALUE]];                

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

    console.log('cppFuncInput:');
    console.log(cppFuncInput);

    return cppFuncInput;
}

// Allocate memory for buffers for array data
function allocateMemoryForBuffer(module, inputs) {
    for(let arg of inputs) { 
        let type = arg.type;

        switch(type) {

            case NUM_TYPE: // no memory allocation for numbers
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
                arg.data.buf = module._malloc(arg.data.numOfRows * typeMap[type].BYTES_PER_ELEMENT);
                break;

            case INT_COLUMNS_TYPE:             
            case NEW_INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:             
            case NEW_FLOAT_COLUMNS_TYPE: // allocation memory for columns that are created
                arg.data.buf = module._malloc(arg.data.numOfRows * arg.data.numOfColumns 
                    * typeMap[type].BYTES_PER_ELEMENT);
                break;

            // TODO: process other cases and mistakes
            default:
                break; 
        }  
    }
    
    console.log('inputs after memory allocation:');
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
            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
                params.push(arg.data.buf);                
                params.push(arg.data.numOfRows);
                break;
            
            case INT_COLUMNS_TYPE:
            case NEW_INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:
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
            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
                types.push('number');                
                types.push('number');
                break;
            
            case INT_COLUMNS_TYPE:
            case NEW_INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:
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
    let shift;
    let heap;
    
    for(let arg of inputs) {
        let type = arg.type;

        switch (type) {

            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:                
                break;
                
            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:                
                shift = shiftMap[type];
                heap = module[heapMap[type]];  
                heap.set(arg.data.array, arg.data.buf >> shift);
                break;
            
            case INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:
                shift = shiftMap[type];
                heap = module[heapMap[type]];
                let numOfBytes = typeMap[type].BYTES_PER_ELEMENT;
                let buf = arg.data.buf;
                let numOfColumns = arg.data.numOfColumns;
                let numOfRows = arg.data.numOfRows;
                let arrays = arg.data.arrays;

                for(let i = 0; i < numOfColumns; i++)
                    heap.set(arrays[i], (buf + i * numOfRows * numOfBytes) >> shift);
                                
                break;

            case NEW_INT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:                
                break;

            // TODO: process other cases and mistakes
            default:
                break;
        }
    }
}

// Get array data from buffer
function getDataFromBuffer(module, inputs) {

    let heap;
    let numOfRows;
    let numOfCols;
    let numOfBytes;    
    let buf;

    for(let arg of inputs) {
        let type = arg.type;
        switch (type) { 
            
            case NUM_TYPE:
            case INT_TYPE:
            case DOUBLE_TYPE:
                break;

            case INT_COLUMN_TYPE:
            case FLOAT_COLUMN_TYPE:
            case FLOAT_COLUMNS_TYPE:  
            case INT_COLUMNS_TYPE:          
                break;

            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
                heap = module[heapMap[type]];
                numOfRows = arg.data.numOfRows;
                numOfBytes = typeMap[type].BYTES_PER_ELEMENT;    
                buf = arg.data.buf;
                let array = new typeMap[type](numOfRows);

                for(let j = 0; j < numOfRows; j++)
                    array[j] = heap[buf / numOfBytes + j];

                arg.array = array;

                break;

            // get arrays that will be further used for newly created columns
            case NEW_INT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMNS_TYPE: 
                heap = module[heapMap[type]];
                numOfRows = arg.data.numOfRows;
                numOfCols = arg.data.numOfColumns;
                numOfBytes = typeMap[type].BYTES_PER_ELEMENT;    
                buf = arg.data.buf;
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
            case INT_COLUMNS_TYPE:
            case NEW_INT_COLUMNS_TYPE:
            case FLOAT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
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
    const typeToColumnCreatorMap = {'newFloatColumns': DG.Column.fromFloat32Array,
                                    'newIntColumns': DG.Column.fromInt32Array,
                                    'newFloatColumn': DG.Column.fromFloat32Array,
                                    'newIntColumn': DG.Column.fromInt32Array};

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
            case FLOAT_COLUMNS_TYPE:
            case INT_COLUMNS_TYPE:
                break;

            case NEW_FLOAT_COLUMN_TYPE:
            case NEW_INT_COLUMN_TYPE:
                let name;

                if(arg.name == undefined)
                    name = (0).toString();
                else 
                    names = arg.name;

                arg.column = typeToColumnCreatorMap[arg.type](name,
                    argsAfterWasmCall[i].array);
                break;

            case NEW_INT_COLUMNS_TYPE:
            case NEW_FLOAT_COLUMNS_TYPE:
                let columns = [];
                let length = argsAfterWasmCall[i].arrays.length;

                let names = [];
                if(arg.names == undefined)
                    for(let k = 1; k <= length; k++)
                        names.push((k).toString());
                else names = arg.names;

                for(let j = 0; j < length; j++)
                    columns.push(typeToColumnCreatorMap[arg.type](names[j],
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

        case COLUMN:
            return funcSpecification.arguments[output.source].column;
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

    console.log('params:');
    console.log(params);

    // create array of parameters' types that are passed to C++-function
    let types = getArrOfWasmTypes(args);
    
    console.log('types:');
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