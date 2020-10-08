/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();
fetch('algs.wasm').then(response => response.arrayBuffer())
.then(bytes => instantiate(bytes, importObject))
.then(instance => instance.exports.e());


//name: Moving Avarage Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: int window_size [SMA Window Size]
export function SMA_filter(data,column_to_filter,window_size) {
    let js_wrapped_sma = Module.cwrap("sma", "null", ["number","number","number"]);
    let column_name = column_to_filter.name +' SMA Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length,window_size);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Exponential Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double filter_ratio [Exponential Filter Parameter]
export function Exp_filter(data,column_to_filter,filter_ratio) {
    let js_wrapped_exp = Module.cwrap("exps", "null", ["number","number","number"]);
    let column_name = column_to_filter.name +' Exponentially Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    console.log(dataPtr);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_exp(dataPtr,filter_array.length,filter_ratio);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}


//name: Kalman Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double Q [Covariance of the process noise]
//input: double R [Covariance of the observation noise]
//input: double P [a posteriori estimate covariance]
export function Kalman_filter(data,column_to_filter,Q,R,P) {
    let js_wrapped_kalm = Module.cwrap("kalman", "null", ["number","number","number","number","number"]);
    let column_name = column_to_filter.name +' Kalman Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_kalm(dataPtr,filter_array.length,Q,R,P)
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) { column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Min Max Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export function MinMax_transform(data,column_to_filter) {
    let js_wrapped_sma = Module.cwrap("minmax", "null", ["number","number"]);
    let column_name = column_to_filter.name +' Min Max Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Z-score Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export function Zscore_transform(data,column_to_filter) {
    let js_wrapped_sma = Module.cwrap("zscore", "null", ["number","number"]);
    let column_name = column_to_filter.name +' Z-score Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Box Cox Transform
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lambda
//input: double ofset
export function box_cox_transform(data,column_to_filter,lambda, ofset) {
    let js_wrapped_sma = Module.cwrap("boxcox", "null", ["number","number","number","number"]);
    let column_name = column_to_filter.name +' Box Cox Transformed';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr,filter_array.length,lambda,ofset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}

//name: Fourier Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lowcut
//input: double hicut
export function fourier_filter(data,column_to_filter,lowcut, hicut) {
    let js_wrapped_ff = Module.cwrap("ffilter", "null", ["number","number","number","number"]);
    let column_name = column_to_filter.name +' Fourier Filtered (L: '+lowcut + '; H: ' + hicut + ')';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = Module._malloc(nDataBytes);
    let dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_ff(dataPtr,filter_array.length,lowcut,hicut);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    let i=0;
    for (i = 0; i < filter_array.length; i++) {column[i]=result[i];}
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    Module._free(dataHeap.byteOffset);
}