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

//name: sumOfInts
//input: int a
//input: int b
//output: int sum
export function sumOfInts(a, b) {

  let x = a;
  let y = b;
  
  let paramTypes = ['number', 'number'];
  let params = [x, y];
  let paramsToUpdate = [];
  

  return cFuncWrapper(NewModule, 'sumOfInts', 'number', paramTypes, params, paramsToUpdate);
}

//name: minOfColumn
//input: dataframe df
//input: column col
//output: int num
export function minOfColumn(df, col) {

  let arrayLength = col.length;

  let array = col.getRawData();

  let paramTypes = ['Int32Array', 'number'];
  let params = [array, arrayLength];
  let paramsToUpdate = [];

  return cFuncWrapper(NewModule, 'minOfArray', 'number', paramTypes, params, paramsToUpdate);  
}


/*
  name: powersOfInts (non-export)
  input: column col
  output: columns result
*/
function powersOfInts(col){
  let arrLength = col.length;
  let squaredLength = col.length;
  let cubicLength = col.length;
  let fourthLength = col.length;

  let arr = col.getRawData();
  let squared = new Int32Array(squaredLength);
  let cubic = new Int32Array(cubicLength);
  let fourth = new Int32Array(fourthLength);

  let params = [arr, arrLength, 
                squared, squaredLength, 
                cubic, cubicLength, 
                fourth, fourthLength];

  let paramTypes = ['Int32Array', 'number',
                    'Int32Array', 'number',
                    'Int32Array', 'number',
                    'Int32Array', 'number'];
  
  let paramsToUpdate = [2, 4, 6];

  cFuncWrapper(NewModule, 'powers', null, paramTypes, params, paramsToUpdate);

  let result = [];

  result.push(DG.Column.fromInt32Array('squared', squared));
  result.push(DG.Column.fromInt32Array('cubic', cubic));
  result.push(DG.Column.fromInt32Array('fourth', fourth));

  return result;
}

//name: powersOfColumn
//input: dataframe df
//input: column col
export function powersOfColumn(df, col) {
  let result = powersOfInts(col);

  for(let element of result) {
    df.columns.add(element);
  }
}

//name: doubleColumn
//input: dataframe df
//input: column col
export function doubleColumn(df, col) {
  let arrayLength = col.length;

  let array = col.getRawData();

  let paramTypes = ['Float32Array', 'number'];
  let params = [array, arrayLength];
  let paramsToUpdate = [0];

  cFuncWrapper(NewModule, 'doubleArray', 'number', paramTypes, params, paramsToUpdate);
}

/*
  name: sumOfColumns (non-export)
  input: column col1
  input: column col2
  output: column result
*/
function sumOfColumns(col1, col2) {

  let arr1Length = col1.length;
  let arr2Length = col2.length;
  let sumLength = col1.length;

  let arr1 = col1.getRawData();
  let arr2 = col2.getRawData();
  let sum = new Float32Array(sumLength);

  let params = [arr1, arr1Length, 
                arr2, arr2Length, 
                sum, sumLength];

  let paramTypes = ['Float32Array', 'number',
                    'Float32Array', 'number',
                    'Float32Array', 'number'];

  let paramsToUpdate = [4];
  
  cFuncWrapper(NewModule, 'sumOfArrays', null, paramTypes, params, paramsToUpdate);

  let result = DG.Column.fromFloat32Array("result", sum);

  return result;
}

//name: sumOfColumnsTest
//input: dataframe df
//input: column col1
//input: column col2
export function sumOfColumnsTest(df, col1, col2) {
  df.columns.add(sumOfColumns(col1, col2));
}

/*
   name: sumAndProdOfColumns (non-export)
   input: column col1
   input: column col2
   output: columns result
*/
function sumAndProdOfColumns(col1, col2) {
  let arr1Length = col1.length;
  let arr2Length = col2.length;
  let sumLength = col1.length;
  let prodLength = col1.length;

  let arr1 = col1.getRawData();
  let arr2 = col2.getRawData();
  let sum = new Float32Array(sumLength);
  let prod = new Float32Array(prodLength);

  let params = [arr1, arr1Length, 
                arr2, arr2Length, 
                sum, sumLength,
                prod, prodLength];

  let paramTypes = ['Float32Array', 'number',
                    'Float32Array', 'number',
                    'Float32Array', 'number',
                    'Float32Array', 'number'];

  let paramsToUpdate = [4, 6];
  
  cFuncWrapper(NewModule, 'sumAndProdOfArrays', null, paramTypes, params, paramsToUpdate);

  let result = [];
  result.push(DG.Column.fromFloat32Array("sum", sum));
  result.push(DG.Column.fromFloat32Array("prod", prod));

  return result;
}

//name: sumAndProdOfColumnsTest
//input: dataframe df
//input: column col1
//input: column col2
export function sumAndProdOfColumnsTest(df, col1, col2) {
  let result = sumAndProdOfColumns(col1, col2);

  for(let col of result) {
    df.columns.add(col);
  }
}






