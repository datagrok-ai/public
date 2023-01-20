// Runtime system for exported C/C++-functions call

// type-to-heap correspondence
const heapMap = {"i32": "HEAP32", // Int32Array
                 "f32": "HEAPF32" // Float32Array
                };

// type signature to typed array ,app
const typeMap = {"i32": Int32Array, // Int32Array
                 "f32": Float32Array // Float32Array
                }; 
 
// type-to-shift map
const shiftMap = {"i32": 2, // Int32Array
                  "f32": 2 // Float32Array
                 };

// type-to-column_creator map
const typeToColumnCreatorMap = {"i32": DG.Column.fromInt32Array,
                                "f32": DG.Column.fromFloat32Array};
 

// CLASSES THAT ARE USED BY CPP-WRAPPER

// simple argument: a number 
class Arg{

    constructor(data) {
        this.data = data;        
    }

    complementArrOfParams(arrOfParams) {
        arrOfParams.push(this.data);
    }

    complementArrOfTypes(arrOfTypes) {
        arrOfTypes.push('number');        
    }

    allocateMemoryForBuffer(module) {}

    isMemoryForBufferAllocated(){
        return true;
    }

    putDataToBuffer(module) {}

    getDataFromBuffer(module) {}

    freeBuffer(module) {}
}

// column argument
class ArgColumn extends Arg {

    constructor(data, targetType, toUpdate = false) {
        super(data);
        this.type = targetType;
        this.toUpdate = toUpdate;
        this.buf = 0;    
        this.numOfRows = data.length;
    }

    complementArrOfParams(arrOfParams) {
        arrOfParams.push(this.buf);
        arrOfParams.push(this.numOfRows);
    }

    complementArrOfTypes(arrOfTypes) {
        arrOfTypes.push('number');
        arrOfTypes.push('number');        
    }

    allocateMemoryForBuffer(module) {
        this.buf = module._malloc(this.numOfRows * typeMap[this.type].BYTES_PER_ELEMENT);
    }

    isMemoryForBufferAllocated(){
        return (this.buf != 0);
    }

    putDataToBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let shift = shiftMap[type];
            let heap = module[heapMap[type]];            
            let array = null;
            let col = this.data;

            if(((col.type == 'int') && (type == 'i32')) 
               || ((col.type == 'double') && (type == 'f32')))
              array = col.getRawData();
            else
              array = new typeMap[type](col.getRawData());

            if(array) 
              heap.set(array, this.buf >> shift);
        }
    }

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let buffer = this.buf;
            let bytes = typeMap[type].BYTES_PER_ELEMENT;
            let array = this.data.getRawData();

            for(let i = 0; i < this.numOfRows; i++)
                array[i] = heap[buffer / bytes + i];
        }
    }

    freeBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
           module._free(this.buf);
           this.buf = 0;
        }
    }

}

// new column argument: a new column is created
class ArgNewColumn extends ArgColumn {
    constructor(targetType, numOfRows) {
        super([], targetType, true);
        this.numOfRows = numOfRows;
    }

    putDataToBuffer(module) {}

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let buf = this.buf;
            
            let columnCreator = typeToColumnCreatorMap[type];   
            
            let arr = new typeMap[type](this.numOfRows);            

            for(let i = 0; i < arr.length; i++)
                arr[i] = heap[buf / arr.BYTES_PER_ELEMENT + i];

            this.data = columnCreator('name', arr);

            //console.log(arr);
            
            // This makes mistake, when this.numOfRows = 5 (TODO: investigate why)
            //this.data = columnCreator('name', new typeMap[type](heap.buffer, this.buf, this.numOfRows));

            //console.log(this.data.getRawData());
        }
    }
}

// an array of columns argument
class ArgColumns extends Arg {

    constructor(data, targetType, toUpdate = false) {
        super(data);
        this.type = targetType;
        this.toUpdate = toUpdate;
        this.buf = 0;    
        this.numOfColumns = data.length;
        this.numOfRows = data[0].length;
    }

    complementArrOfParams(arrOfParams) {
        arrOfParams.push(this.buf);
        arrOfParams.push(this.numOfRows);
        arrOfParams.push(this.numOfColumns);
    }

    complementArrOfTypes(arrOfTypes) {
        arrOfTypes.push('number');
        arrOfTypes.push('number');
        arrOfTypes.push('number');        
    }

    allocateMemoryForBuffer(module) {
        this.buf = module._malloc(this.numOfRows * this.numOfColumns 
                                  * typeMap[this.type].BYTES_PER_ELEMENT);
    }

    isMemoryForBufferAllocated(){
        return (this.buf != 0);
    }

    putDataToBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let shift = shiftMap[type];
            let heap = module[heapMap[type]];
            let numOfBytes = typeMap[type].BYTES_PER_ELEMENT;
            
            // put columns data to buffer
            for(let i = 0; i < this.numOfColumns; i++) {
                let array = null;
                let col = this.data[i];

                if(((col.type == 'int') && (type == 'i32')) 
                   || ((col.type == 'double') && (type == 'f32')))
                     array = col.getRawData();
                else
                     array = new typeMap[type](col.getRawData());

                // check data array
                if(array != null) 
                    heap.set(array, (this.buf + i * this.numOfRows * numOfBytes) >> shift);                                   
            }
        }
    }

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let numOfRows = this.numOfRows;
            let numOfCols = this.numOfColumns;                 
            let arr = new typeMap[type](heap.buffer, this.buf, numOfRows * numOfCols);
       
            for(let i = 0; i < numOfCols; i++) {
                let colData = this.data[i].getRawData();
                for(let j = 0; j < numOfRows; j++)
                  colData[j] = arr[j + i * numOfRows];                    
            }
        }
    }

    freeBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
           module._free(this.buf);
           this.buf = 0;
        }
    }

}

// an array of new columns: new columns are created
class ArgNewColumns extends ArgColumns {

    constructor(targetType, numOfRows, numOfColumns, names) {
        super([[]], targetType, true);  
        this.data = [];        
        this.numOfColumns = numOfColumns;
        this.numOfRows = numOfRows;
        this.names = names;
    }

    putDataToBuffer(module) { }

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let numOfRows = this.numOfRows;
            let numOfCols = this.numOfColumns;
            let numOfBytes = typeMap[type].BYTES_PER_ELEMENT;                       
            let columnCreator = typeToColumnCreatorMap[type];
            let buf = this.buf;
            let names = this.names;

            for(let i = 0; i < numOfCols; i++) {
                let arr = new typeMap[type](numOfRows);

                for(let j = 0; j < numOfRows; j++)
                    arr[j] = heap[buf / numOfBytes + j + i * numOfRows];
                
                //this.data.push(columnCreator((i + 1).toString(), arr));
                this.data.push(columnCreator(names[i], arr));
            }

            // create columns: here, may be a problem when numOfRows = 5            
            /*for(let i = 0; i < numOfCols; i++)
                this.data.push(columnCreator((i + 1).toString(), new typeMap[type](heap.buffer, 
                    this.buf + i * numOfRows * numOfBytes, numOfRows)));*/
        }
    }
}

// a wrapper for exported C/C++-function call
function cppFuncWrapper(module, cFuncName, returnType, args)
{
    let result;

    // allocate memory for buffers
    for(let arg of args)         
        arg.allocateMemoryForBuffer(module);
    
    let isEnoughOfMemoryAllocated = true;

    // check memory allocation
    for(let arg of args)
        isEnoughOfMemoryAllocated &= arg.isMemoryForBufferAllocated();  
    
    // run exported function if enough of memory is allocated
    if(isEnoughOfMemoryAllocated) {

        let params = []; // arguments that are put to the exported function
        let types = []; // their types

        // prepare data that is put to exported function
        for(let arg of args) {
            arg.complementArrOfParams(params);
            arg.complementArrOfTypes(types);
            arg.putDataToBuffer(module);
        }

        let extendedTypeOfReturn = (returnType == 'num') ? 'number' : null;   
        
        /*// call exported function
        if(extendedTypeOfReturn)
            result = module.ccall(cFuncName, extendedTypeOfReturn, types, params);
        else
            result = module.ccall(cFuncName, extendedTypeOfReturn, types, params);*/

        result = module.ccall(cFuncName, extendedTypeOfReturn, types, params);
                    
        // update and get data from buffers if required
        for(let arg of args)
            arg.getDataFromBuffer(module);       
    }

    // clear buffers
    for(let arg of args)         
        arg.freeBuffer(module);

    if(result != undefined)
        return result;    
} // cppFuncWrapper


// A LAYER BETWEEN JS AND CPP-WRAPPER

// Parameters creator
let Param = {
    floatColumn (column) {
        return new ArgColumn(column, 'f32');
    },

    newFloatColumn (numOfRows) {
        return new ArgNewColumn('f32', numOfRows);
    },

    floatColumns(columns){
        return new ArgColumns(columns.toList(), 'f32');
    },

    newFloatColumns(numOfRows, numOfColumns, names) {
        return new ArgNewColumns('f32', numOfRows, numOfColumns, names);
    },

    int(number) {
        return new Arg(number);
    },

    num(number) {
        return new Arg(number);
    }
};

// Return value creator
let Return = {
    tableFromColumns(argColumns) {
        return DG.DataFrame.fromColumns(argColumns.data);
    },

    num(number) {
        return number;
    },

    int(number) {
        return number;
    },

    double(number) {
        return number;
    },

    column(argColumn) {
        return argColumn.data;
    }
};

// The main tool that combines all together
export async function callWasm(module, funcName, inputs) {

    let start = new Date().getTime();

    // get specification of exported C/C++-function
    let funcSpecification = module[funcName];

    // get argumnets
    let args = funcSpecification.arguments; 

    // array of arguments that further are used by cpp-wrapper
    let cppFuncInput = [];

    // complete an input for cpp
    let i = 0;
    for(let key in args){
        let arg = args[key]; 
        
        // skip auxiliry element
        if(key == '_callResult')             
            continue;

        // create an argument
        switch(arg.type){
            case 'floatColumns':
            case 'int':
            case 'num':
            case 'floatColumn':
                arg.data = Param[arg.type](inputs[i]);
                i++;
                break;
            case 'newFloatColumns':
                let val1 = args[ arg['numOfRows']['ref'] ].data[arg['numOfRows']['value']];
                let val2 = args[ arg['numOfColumns']['ref'] ].data[arg['numOfColumns']['value']];                
                arg.data = Param[arg.type](val1, val2, arg.names);
                break;
            case 'newFloatColumn':
                let val = args[ arg['numOfRows']['ref'] ].data[arg['numOfRows']['value']];
                arg.data = Param[arg.type](val);
                break;
            } // switch

        cppFuncInput.push(args[key].data);

    } // for key    

    // CALL EXPORTED CPP-FUNCTION
    let callResult = cppFuncWrapper(module, funcName, 'num', cppFuncInput);

    console.log('C++-function call result: ' + callResult);

    // store result that is returned by exported cpp-function
    args._callResult = Param.num(callResult);

    // create output
    let output = funcSpecification.output;

    let finish = new Date().getTime();

    console.log(`Time for C/C++-function is ${finish - start} ms.`)

    // if a single object must be returned
    if(output['type'] != 'objects')
      return Return[output['type']](args[output['source']].data);

    let arrayToReturn = [];

    // push data of the required arguments
    for(let name of output['source'])
      arrayToReturn.push(args[name].data.data);

    return arrayToReturn;
  } // callWasm