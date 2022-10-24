/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initNewModule();  
}

// RUNTIME SYSTEM FOR EXPORTED C-FUNCTIONS CALL

// Wrapper for exported C-functions
function cFuncWrapper(module, cFuncName, returnType, argTypes, args, argsToUpdate=[]){    

    // type-to-heap correspondence
    let heapMap = {"i32a": "HEAP32", // Int32Array
                   "f32a": "HEAPF32", // Float32Array
                   "i32c": "HEAP32", // column that contains Int32Array
                   "f32c": "HEAPF32"}; // column that contains Float32Array
                   // TODO: add other items
  
    let params = []; // parameters that are further substituted to C-function
    let paramTypes = []; // types of parameters that are further substituted to C-function
    let bufs = []; // buffers that are used for arrays    
    let buffer; // the current buffer 
    let extendedTypeOfReturn = (returnType == 'num') ? 'number' : null; // C-function return type
  
    // creation of parameters that are substituted to the exported C-function
    for(let i = 0; i < args.length; i++){          
        let argType = argTypes[i];
  
        if(argType == 'num') {
            params.push(args[i]);
            paramTypes.push('number');
            bufs.push(null); // i.e. no allocated memory is used
        }
  
        else { // the current argument is an array, so some buffer/heap routine shoude be performed
            let arr; // raw data array
            
            // set arr taking into account the current structure feature: array or column
            switch(argType[argType.length - 1]) { // check the last symbol: 'a' - array, 'c' - column
              case 'a':
                arr = args[i];
                break;
              case 'c':
                arr = args[i].getRawData();
                break;
            }
            
            buffer = module._malloc(arr.BYTES_PER_ELEMENT * arr.length); // allocate memory
  
            let heap = heapMap[argType]; // choose the required heap            
            let shift = 0; // current shift
            
            // compute shift
            switch(heap) {
                case "HEAP32": case "HEAPF32": 
                    shift = 2; // this shift is specified by EMSCRIPTEN 
                    break;
                //TODO: add processing with respect to other types
            }
  
            // assign the data to the heap
            module[heap].set(arr, buffer >> shift);
  
            // append arrays that are passed to C-function
            params.push(buffer);
            paramTypes.push('number'); // in the case of array, 'number' is also used
            bufs.push(buffer);
  
            // append length of array
            params.push(arr.length);
            paramTypes.push('number');
        }
    }
  
    // call C-function
    let result;   
    if(extendedTypeOfReturn)
        result = module.ccall(cFuncName, extendedTypeOfReturn, paramTypes, params);
    else
        module.ccall(cFuncName, extendedTypeOfReturn, paramTypes, params);
  
    // process arrays that are modified, when executing C-function
    for(let j = 0; j < argsToUpdate.length; j++){
        let m = argsToUpdate[j];
        
        let argType = argTypes[m];
  
        if(argType != 'num') {     
            let arr; // raw data array
  
            // set arr taking into account the current structure feature: array or column
            switch(argType[argType.length - 1]) { // check the last symbol: 'a' - array, 'c' - column
              case 'a':
                arr = args[m];
                break;
              case 'c':
                arr = args[m].getRawData();
                break;
            }
  
            let heap = heapMap[argType]; 
            let buffer = bufs[m];
            let bytes = arr.BYTES_PER_ELEMENT;
            
            // update the current array
            for(let n = 0; n < arr.length; n++)
                arr[n] = module[heap][buffer/bytes + n];                
        } 
    }
    
    // clear allocated memory
    for(let k = 0; k < bufs.length; k++)
        if(bufs[k])
            module._free(bufs[k]);
    
    if(returnType)
        return result;
  }


// THE FOLLOWING CODE IS GENERATED AUTOMATICALLY

// EXPORTED C-FUNCTIONS

// Functions from lib1.c

//name: sumOfIntegers
//input: int x
//input: int y
//output: int z
export function sumOfIntegers(x, y) {
  return cFuncWrapper(NewModule, 'sumOfInts',  'num',
                      ['num', 'num'],
                      [x, y]);
}

//name: minOfColumn
//input: dataframe df
//input: column col
//output: int num
export function minOfColumn(df, col) {
  return cFuncWrapper(NewModule, 'minOfArray',  'num',
                      ['i32c'],
                      [col]);
}



// Functions from lib2.c

/*
   name: sumOfColumns (non-export)
   input: column col1
   input: column col2
   output: column result
*/
function sumOfColumns(col1, col2) {
  let sum = new Float32Array(col1.length);
  cFuncWrapper(NewModule, 'sumOfArrays',  null,
               ['f32c', 'f32c', 'f32a'],
               [col1, col2, sum],
               [2]);
  return DG.Column.fromFloat32Array('sum', sum);
}

/*
   name: sumAndProdOfColumns (non-export)
   input: column col1
   input: column col2
   output: columns result
*/
function sumAndProdOfColumns(col1, col2) {
  let sum = new Float32Array(col1.length);
  let prod = new Float32Array(col1.length);
  cFuncWrapper(NewModule, 'sumAndProdOfArrays',  null,
               ['f32c', 'f32c', 'f32a', 'f32a'],
               [col1, col2, sum, prod],
               [2, 3]);
  return [DG.Column.fromFloat32Array('sum', sum),
          DG.Column.fromFloat32Array('prod', prod)];
}

/*
   name: powersOfInts (non-export)
   input: column col
   output: columns result
*/
function powersOfInts(col) {
  let squared = new Int32Array(col.length);
  let cubic = new Int32Array(col.length);
  let fourth = new Int32Array(col.length);
  cFuncWrapper(NewModule, 'powers',  null,
               ['i32c', 'i32a', 'i32a', 'i32a'],
               [col, squared, cubic, fourth],
               [1, 2, 3]);
  return [DG.Column.fromInt32Array('squared', squared),
          DG.Column.fromInt32Array('cubic', cubic),
          DG.Column.fromInt32Array('fourth', fourth)];
}


// THE FOLLOWING FUNCTIONS ARE ADDED MANUALLY

//name: sumOfColumnsTest
//input: dataframe df
//input: column col1
//input: column col2
export function sumOfColumnsTest(df, col1, col2) {
  df.columns.add(sumOfColumns(col1, col2));
}

//name: sumAndProdOfColumnsTest
//input: dataframe df
//input: column col1
//input: column col2
export function sumAndProdOfColumnsTest(df, col1, col2) {
  let result = sumAndProdOfColumns(col1, col2);  
  for(let col of result) 
      df.columns.add(col);
}

//name: powersOfColumn
//input: dataframe df
//input: column col
export function powersOfColumn(df, col) {
  let result = powersOfInts(col);
  for(let element of result) 
    df.columns.add(element);
 }


