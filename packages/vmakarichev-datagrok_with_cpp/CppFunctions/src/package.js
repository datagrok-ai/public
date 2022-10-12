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

// EXPORTED C-FUNCTIONS

// Functions from lib1.c

//name: sumOfIntegers
//input: int x
//input: int y
//output: int z
export function sumOfIntegers(x, y) {
  let a = x;
  let b = y;

  // Call exported C-function
  let exportedFunctionOutput = NewModule._sumOfInts(a, b);

  return exportedFunctionOutput;
}

//name: minOfColumn
//input: dataframe df
//input: column col
//output: int num
export function minOfColumn(df, col) {
  let length = col.length;

  // Buffer routine: array
  let rawData_array = col.getRawData();
  let length_array = rawData_array.length;
  let numOfBytes_array = rawData_array.BYTES_PER_ELEMENT * length_array;
  let dataPtr_array = NewModule._malloc(numOfBytes_array);
  let dataHeap_array = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_array, numOfBytes_array);
  dataHeap_array.set(new Uint8Array(rawData_array.buffer));

  // Call exported C-function
  let exportedFunctionOutput = NewModule._minOfArray(dataPtr_array, length);

  // Cleaning memory allocated above
  NewModule._free(dataPtr_array);

  return exportedFunctionOutput;
}

// Functions from lib2.c

//name: doubleColumn
//input: column col
//output: column result
export function doubleColumn(col) {
  let length = col.length;

  // Buffer routine: array
  let rawData_array = col.getRawData();
  let length_array = rawData_array.length;
  let numOfBytes_array = rawData_array.BYTES_PER_ELEMENT * length_array;
  let dataPtr_array = NewModule._malloc(numOfBytes_array);
  let dataHeap_array = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_array, numOfBytes_array);
  dataHeap_array.set(new Uint8Array(rawData_array.buffer));

  // Call exported C-function
  NewModule._doubleArray(dataPtr_array, length);

  // Creating output: the column "result"
  let dataArray_result = new Float32Array(dataHeap_array.buffer, dataHeap_array.byteOffset, length_array);
  let result = DG.Column.fromFloat32Array("result", dataArray_result);

  // Cleaning memory allocated above
  NewModule._free(dataPtr_array);

  return result;
}

//name: sumOfColumns
//input: column col1
//input: column col2
//output: column result
export function sumOfColumns(col1, col2) {
  let length = col1.length;

  // Buffer routine: arr1
  let rawData_arr1 = col1.getRawData();
  let length_arr1 = rawData_arr1.length;
  let numOfBytes_arr1 = rawData_arr1.BYTES_PER_ELEMENT * length_arr1;
  let dataPtr_arr1 = NewModule._malloc(numOfBytes_arr1);
  let dataHeap_arr1 = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_arr1, numOfBytes_arr1);
  dataHeap_arr1.set(new Uint8Array(rawData_arr1.buffer));

  // Buffer routine: arr2
  let rawData_arr2 = col2.getRawData();
  let length_arr2 = rawData_arr2.length;
  let numOfBytes_arr2 = rawData_arr2.BYTES_PER_ELEMENT * length_arr2;
  let dataPtr_arr2 = NewModule._malloc(numOfBytes_arr2);
  let dataHeap_arr2 = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_arr2, numOfBytes_arr2);
  dataHeap_arr2.set(new Uint8Array(rawData_arr2.buffer));

  // Buffer routine: sum
  let length_sum = length;
  let numOfBytes_sum = Float32Array.BYTES_PER_ELEMENT * length_sum;
  let dataPtr_sum = NewModule._malloc(numOfBytes_sum);
  let dataHeap_sum = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_sum, numOfBytes_sum);

  // Call exported C-function
  NewModule._sumOfArrays(dataPtr_arr1, dataPtr_arr2, dataPtr_sum, length);

  // Creating output: the column "result"
  let dataArray_result = new Float32Array(dataHeap_sum.buffer, dataHeap_sum.byteOffset, length_sum);
  let result = DG.Column.fromFloat32Array("result", dataArray_result);

  // Cleaning memory allocated above
  NewModule._free(dataPtr_arr1);
  NewModule._free(dataPtr_arr2);
  NewModule._free(dataPtr_sum);

  return result;
}

//name: doubleColumnTest
//input: dataframe df
//input: column col
export function doubleColumnTest(df, col) {
  df.columns.add(doubleColumn(col));
}

//name: sumOfColumnsTest
//input: dataframe df
//input: column col1
//input: column col2
export function sumOfColumnsTest(df, col1, col2) {
  df.columns.add(sumOfColumns(col1, col2));
}