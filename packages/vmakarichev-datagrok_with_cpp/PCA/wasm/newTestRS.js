var factory = require('./eigenPCA');

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
};

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
            // TODO: check type & add getRawData()
            let array = new typeMap[this.type](this.data); 
            heap.set(array, this.buf >> shift);
        }
    }

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let buffer = this.buf;
            let bytes = typeMap[type].BYTES_PER_ELEMENT;
            //let array = new typeMap[type](heap.buffer, this.buf, this.numOfRows);
            // TODO: columns routine!
            //this.data = array;
            for(let i = 0; i < this.numOfRows; i++)
                this.data[i] = heap[buffer / bytes + i];
        }
    }

    freeBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
           module._free(this.buf);
           this.buf = 0;
        }
    }

};

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
            let array = new typeMap[type](heap.buffer, this.buf, this.numOfRows);
            // TODO: columns routine!
            this.data = array;            
        }
    }
};

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

                // TODO: add data type check
                array = new typeMap[type](this.data[i]);

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

            // TODO: columns routine!            
            for(let i = 0; i < numOfCols; i++)
                for(let j = 0; j < numOfRows; j++)
                    this.data[i][j] = arr[j + i * numOfRows];                    
        }
    }

    freeBuffer(module) {
        if(this.isMemoryForBufferAllocated()) {
           module._free(this.buf);
           this.buf = 0;
        }
    }

};

class ArgNewColumns extends ArgColumns {

    constructor(targetType, numOfRows, numOfColumns) {
        super([[]], targetType, true);  
        this.data = [];        
        this.numOfColumns = numOfColumns;
        this.numOfRows = numOfRows;
    }

    putDataToBuffer(module) { }

    getDataFromBuffer(module) {
        if(this.toUpdate && this.isMemoryForBufferAllocated()) {
            let type = this.type;
            let heap = module[heapMap[type]];
            let numOfRows = this.numOfRows;
            let numOfCols = this.numOfColumns;
            let numOfBytes = typeMap[type].BYTES_PER_ELEMENT;                 
            //let arr = new typeMap[type](heap.buffer, this.buf, numOfRows * numOfCols);

            // TODO: columns routine!            
            for(let i = 0; i < numOfCols; i++)
                this.data.push(new typeMap[type](heap.buffer, 
                    this.buf + i * numOfRows * numOfBytes, numOfRows));
        }
    }
};

function cppFuncWrapper(module, cFuncName, returnType, arguments)
{
    let result;

    // allocate memory for buffers
    for(let arg of arguments)         
        arg.allocateMemoryForBuffer(module);
    
    let isEnoughOfMemoryAllocated = true;

    // check memory allocation
    for(let arg of arguments)
        isEnoughOfMemoryAllocated &= arg.isMemoryForBufferAllocated();  
    
    // run exported function if enough of memory is allocated
    if(isEnoughOfMemoryAllocated) {

        let params = []; // arguments that are put to the exported function
        let types = []; // their types

        // prepare data that is put to exported function
        for(let arg of arguments) {
            arg.complementArrOfParams(params);
            arg.complementArrOfTypes(types);
            arg.putDataToBuffer(module);
        }

        let extendedTypeOfReturn = (returnType == 'num') ? 'number' : null;   
        
        // call exported function
        if(extendedTypeOfReturn)
            result = module.ccall(cFuncName, extendedTypeOfReturn, types, params);
        else
            result = module.ccall(cFuncName, extendedTypeOfReturn, types, params);
                    
        // update and get data from buffers if required
        for(let arg of arguments)
            arg.getDataFromBuffer(module);

        //console.log(params);
        //console.log(types);
    }

    // clear buffers
    for(let arg of arguments)         
        arg.freeBuffer(module);

    if(result != undefined)
        return result;    
}


factory().then((instance) => {

    let param = [];
    let types = [];

    // test for a number
    console.log('\nNUMBER TEST!');
    let a = 10;
    let arg = new Arg(a);
    arg.complementArrOfParams(param);
    arg.complementArrOfTypes(types);

    console.log(arg);

    console.log(param);
    console.log(types);


    // test for a column/array
    console.log('\nARRAY TEST!');    
    let arr = new Int32Array([10,20,30,40,50]);
    let col = new ArgColumn(arr, 'f32', true);

    console.log(col);
    console.log(arr);

    console.log('Is memory allocated: ' + col.isMemoryForBufferAllocated());
    col.getDataFromBuffer(instance);

    col.allocateMemoryForBuffer(instance);
    console.log('Is memory allocated: ' + col.isMemoryForBufferAllocated());

    col.putDataToBuffer(instance);

    col.getDataFromBuffer(instance);

    console.log(col);    

    col.complementArrOfParams(param);
    col.complementArrOfTypes(types);

    console.log(param);
    console.log(types);

    //col.freeBuffer(instance);
    console.log(col);

    // test for a column/array
    console.log('\nNEW ARRAY TEST!');
    let newCol = new ArgNewColumn('f32', 7);
    console.log(newCol);

    console.log('Is memory allocated: ' + newCol.isMemoryForBufferAllocated());

    newCol.allocateMemoryForBuffer(instance);
    console.log(newCol);
    console.log('Is memory allocated: ' + newCol.isMemoryForBufferAllocated());

    newCol.putDataToBuffer(instance);

    newCol.getDataFromBuffer(instance);
    console.log(newCol);

    newCol.complementArrOfParams(param);
    newCol.complementArrOfTypes(types);

    console.log(param);
    console.log(types);

    //newCol.freeBuffer(instance);
    console.log(newCol);

    console.log(param);
    console.log(types);

    // test for runtime system: column, column
    console.log('\nRUNTIME SYSTEM: column, column!');
    let arr1 = new Int32Array([21,32,43,34,15]);
    let arr2 = new Int32Array([11,12,13,14,15]);

    let arguments = [new ArgColumn(arr1, 'f32'),
                     new ArgColumn(arr2, 'f32')];
    
    console.log(arguments);

    console.log(cppFuncWrapper(instance, 'error', 'num', arguments));


    // test for a columns

    param = [];
    types = [];

    console.log('\nCOLUMNS TEST!');
    let col1 = new Int32Array([11, -22, 23, 14, -50, 34, 21]);
    let col2 = new Float32Array([26, 17, 18, 29, -10, -3, 32]);
    let col3 = new Int32Array([11, -12, -13, -14, 15, 11, 32]);
    let col4 = new Float32Array([-16, -17, 18, -19, -20, 14, 71]);

    let data = new ArgColumns([col1, col2, col3, col4], 'f32');

    console.log(data);

    data.allocateMemoryForBuffer(instance);

    console.log(data);

    data.complementArrOfParams(param);
    data.complementArrOfTypes(types);

    console.log(param);
    console.log(types);

    data.putDataToBuffer(instance);

    data.getDataFromBuffer(instance);

    //data.freeBuffer(instance);

    console.log(data);


    // test for a new columns

    param = [];
    types = [];

    console.log('\nNEW COLUMNS TEST!');

    let newCols = new ArgNewColumns('f32', 10, 3);
    console.log(newCols);

    newCols.allocateMemoryForBuffer(instance);
    console.log(newCols);

    newCols.putDataToBuffer(instance);

    newCols.getDataFromBuffer(instance);
    console.log(newCols);


    console.log('\nPCA TEST!');
    let argData = new ArgColumns([ new Int32Array([11, -22, 23, 14, -50, 34, 21]), 
                                   new Float32Array([26, 17, 18, 29, -10, -3, 32]),
                                   new Int32Array([11, -12, -13, -14, 15, 11, 32]),
                                   new Float32Array([-16, -17, 18, -19, -20, 14, 71]) ],
                                 'f32');

    let numOfComponents = 4;
    let argNumComp = new Arg(numOfComponents);

    let argPComp = new ArgNewColumns('f32', 7, numOfComponents);
    let argApprox = new ArgNewColumns('f32', 7, 4);

    arguments = [argData, argNumComp, argPComp, argApprox];
    console.log(cppFuncWrapper(instance, 'pcaWithApproximation', 'num', arguments)); 
    
    console.log(argData);

    console.log(argPComp);

    console.log(argApprox);
    

    /*let tArr1 = new Float32Array([11, -22, 23, 14, -50, 34, 21]);

    let arg1 = new ArgColumns(tArr1, 'f32', true);
    console.log(arg1);
    arg1.allocateMemoryForBuffer(instance);
    arg1.putDataToBuffer();
    arg1.getDataFromBuffer();
    arg1.freeBuffer();

/*
    //let arr1 = [1, 2, 3, 4, 5];
    let arr1 = [11, -22, 23, 14, -50, 34, 21];
    let tArr1 = new Float32Array(arr1);
    //console.log(tArr1);
    
    //let arr2 = [6, 7, 8, 9, 10];
    let arr2 = [26, 17, 18, 29, -10, -3, 32];
    let tArr2 = new Int32Array(arr2);
    //console.log(tArr2);
    
    //let arr3 = [11, 12, 13, 14, 15];
    let arr3 = [11, -12, -13, -14, 15, 11, 32];
    let tArr3 = new Float32Array(arr3);    
    //console.log(tArr3);

    //let arr4 = [16, 17, 18, 19, 20];
    let arr4 = [-16, -17, 18, -19, -20, 14, 71];
    let tArr4 = new Int32Array(arr4);
    //console.log(tArr4);
    
    let data = [tArr1, tArr2, tArr3, tArr4];
    let numOfRows = tArr1.length;
    let numOfCols = data.length;

    let arg1 = new ArgColumns(tArr1, 'f32');
    console.log(arg1);

    let arg3 = new ArgColumns([], 'f32', 4, true);
    console.log(arg3);

    let numOfPrincipalComponents = 4;

    let arg2 = new Arg(numOfPrincipalComponents);
    console.log(arg2);
    let arr = [];
    arg2.complementArrOfTypes(arr);
    console.log(arr);

    let princComp = [];   

    //let arg3 = new Arg(princComp, true, 'f32', numOfRows, numOfPrincipalComponents);
    //console.log(arg3);

    let approx = [];

    //let arg4 = new Arg(approx, true, 'f32', numOfRows, numOfCols);
    //console.log(arg4);

    //let args = [arg1, arg2, arg3, arg4];
    //console.log(args);*/
    
});
