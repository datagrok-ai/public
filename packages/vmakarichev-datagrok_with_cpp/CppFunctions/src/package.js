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

// C-FUNCTIONS

// Functions from lib3.c

/* Exported C-function
name: doubleArray
input: float * array
input: int length
output: void result  */
function doubleArray(array, length) { // <--- CHECK THIS!

  // Buffer routine: array. Check it!
  let numOfBytes_array = array.length * array.BYTES_PER_ELEMENT; // <--- CHECK THIS!
  let dataPtr_array = NewModule._malloc(numOfBytes_array);
  let dataHeap_array = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_array, numOfBytes_array);
  dataHeap_array.set(new Uint8Array(array.buffer)); // <--- CHECK THIS!

  NewModule._doubleArray(dataPtr_array, length);

  /* OPTIONAL: each array could be modified when executing exported C-function.
     Remove '//' in order to get the modified array.
     The modified array further can be returned.
     MODIFY if it is required! */

  let arrayModified = new Float32Array(dataHeap_array.buffer, dataHeap_array.byteOffset, array.length);

  // Cleaning allocated memory.
  NewModule._free(dataPtr_array);

  return arrayModified;   //  <--- here, smth can be returned
}

//name: doubleColumn
//input: dataframe dataFrame
//input: column column
export function doubleColumn(dataFrame, column) {

  let result = doubleArray( column.getRawData(), column.length );

  dataFrame.columns.add(DG.Column.fromFloat32Array("2 * " + column.name, result));

}



/* Exported C-function
name: minOfArray
input: int * array
input: int length
output: int result  */
function minOfArray(array, length) { // <--- CHECK THIS!

  // Buffer routine: array. Check it!
  let numOfBytes_array = array.length * array.BYTES_PER_ELEMENT; // <--- CHECK THIS!
  let dataPtr_array = NewModule._malloc(numOfBytes_array);
  let dataHeap_array = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_array, numOfBytes_array);
  dataHeap_array.set(new Uint8Array(array.buffer)); // <--- CHECK THIS!

  let result = NewModule._minOfArray(dataPtr_array, length);

  /* OPTIONAL: each array could be modified when executing exported C-function.
     Remove '//' in order to get the modified array.
     The modified array further can be returned.
     MODIFY if it is required! */

  //let arrayModified = new Int32Array(dataHeap_array.buffer, dataHeap_array.byteOffset, array.length);

  // Cleaning allocated memory.
  NewModule._free(dataPtr_array);

  return result;
}

//name: minOfColumn
//input: dataframe dataFrame
//input: column column
export function minOfColumn(dataFrame, column) {
  let result = minOfArray( column.getRawData(), column.length );

  alert( result );
}



/* Exported C-function
name: sumOfArrays
input: float * arr1
input: float * arr2
input: float * sum
input: int length
output: void result  */
function sumOfArrays(arr1, arr2, length) { // <--- CHECK THIS!

  // Buffer routine: arr1. Check it!
  let numOfBytes_arr1 = arr1.length * arr1.BYTES_PER_ELEMENT; // <--- CHECK THIS!
  let dataPtr_arr1 = NewModule._malloc(numOfBytes_arr1);
  let dataHeap_arr1 = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_arr1, numOfBytes_arr1);
  dataHeap_arr1.set(new Uint8Array(arr1.buffer)); // <--- CHECK THIS!

  // Buffer routine: arr2. Check it!
  let numOfBytes_arr2 = arr2.length * arr2.BYTES_PER_ELEMENT; // <--- CHECK THIS!
  let dataPtr_arr2 = NewModule._malloc(numOfBytes_arr2);
  let dataHeap_arr2 = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_arr2, numOfBytes_arr2);
  dataHeap_arr2.set(new Uint8Array(arr2.buffer)); // <--- CHECK THIS!

  // Buffer routine: sum. Check it!
  //let numOfBytes_sum = sum.length * sum.BYTES_PER_ELEMENT; // <--- CHECK THIS!
  let dataPtr_sum = NewModule._malloc(numOfBytes_arr1);
  let dataHeap_sum = new Uint8Array(NewModule.HEAPU8.buffer, dataPtr_sum, numOfBytes_arr1);
  //dataHeap_sum.set(new Uint8Array(sum.buffer)); // <--- CHECK THIS!

  NewModule._sumOfArrays(dataPtr_arr1, dataPtr_arr2, dataPtr_sum, length);

  /* OPTIONAL: each array could be modified when executing exported C-function.
     Remove '//' in order to get the modified array.
     The modified array further can be returned.
     MODIFY if it is required! */

  //let arr1Modified = new Float32Array(dataHeap_arr1.buffer, dataHeap_arr1.byteOffset, arr1.length);
  //let arr2Modified = new Float32Array(dataHeap_arr2.buffer, dataHeap_arr2.byteOffset, arr2.length);
  let sumModified = new Float32Array(dataHeap_sum.buffer, dataHeap_sum.byteOffset, length);

  // Cleaning allocated memory.
  NewModule._free(dataPtr_arr1);
  NewModule._free(dataPtr_arr2);
  NewModule._free(dataPtr_sum);

  return sumModified;   //  <--- here, smth can be returned
}

//name: sumOfColumns
//input: dataframe dataFrame
//input: column column1
//input: column column2
export function sumOfColumns(dataFrame, column1, column2){
  let result = sumOfArrays(column1.getRawData(), column2.getRawData(), column1.length);

  dataFrame.columns.add(DG.Column.fromFloat32Array(column1.name + " + " + column2.name, result));
}