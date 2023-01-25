/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//top-menu: Tools | DSP Filters | SMA Filter
//name: Moving Average Filter
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
//input: int windowSize [SMA Window Size]
export async function SMA_filter(dataframe, column, windowSize) {
  await initDSPModule();
  let js_wrapped_sma = DSPModule.cwrap("sma", "null", ["number", "number", "number"]);
  //let column_name = column_to_filter.name + ' SMA Filtered';
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_sma(dataPtr, filter_array.length, windowSize);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' SMA Filtered', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Filters | Exponential Filter
//name: Exponential Filter
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
//input: double filterRatio [Exponential Filter Parameter]
export async function Exp_filter(dataframe, column, filterRatio) {
  await initDSPModule();
  let js_wrapped_exp = DSPModule.cwrap("exps", "null", ["number", "number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_exp(dataPtr, filter_array.length, filterRatio);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Exponentially Filtered', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Filters | Kalmam Filter
//name: Kalman Filter
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
//input: double Q [Covariance of the process noise]
//input: double R [Covariance of the observation noise]
//input: double P [a posteriori estimate covariance]
export async function Kalman_filter(dataframe, column, Q, R, P) {
  await initDSPModule();
  let js_wrapped_kalm = DSPModule.cwrap("kalman", "null", ["number", "number", "number", "number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_kalm(dataPtr, filter_array.length, Q, R, P)
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Kalman Filtered', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Min-Max
//name: Min Max Normalization
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
export async function MinMax_transform(dataframe, column) {
  await initDSPModule();
  let js_wrapped_sma = DSPModule.cwrap("minmax", "null", ["number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_sma(dataPtr, filter_array.length);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Min Max Normalized', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Z-score
//name: Z-score Normalization
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
export async function Zscore_transform(dataframe, column) {
  await initDSPModule();
  let js_wrapped_sma = DSPModule.cwrap("zscore", "null", ["number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_sma(dataPtr, filter_array.length);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Z-score Normalized', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Box Cox
//name: Box Cox Transform
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
//input: double lambda
//input: double offset
export async function box_cox_transform(dataframe, column, lambda, offset) {
  await initDSPModule();
  let js_wrapped_sma = DSPModule.cwrap("boxcox", "null", ["number", "number", "number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_sma(dataPtr, filter_array.length, lambda, offset);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Box Cox Transformed', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Get Trend
//name: Get Trend
//input: dataframe dataframe [Input Dataframe]
//input: column column [Column]
export async function get_trend(dataframe, column) {
  await initDSPModule();
  let js_wrapped_gt = DSPModule.cwrap("gettrend", "null", ["number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_gt(dataPtr, filter_array.length);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Trend', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Detrend
//name: Detrend
//input: dataframe dataframe [Input Dataframe]
//input: column column [Column]
export async function remove_trend(dataframe, column) {
  await initDSPModule();
  let js_wrapped_rt = DSPModule.cwrap("removetrend", "null", ["number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_rt(dataPtr, filter_array.length);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Detrended', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Filters | Fourier Filter
//name: Fourier Filter
//input: dataframe dataframe [Input data table]
//input: column column [Signal to Filter]
//input: double lowcut
//input: double hicut
//input: double observationTime[Time taken to record observed signal]
export async function fourier_filter(dataframe, column, lowcut, hicut, observationTime) {
  await initDSPModule();
  let lowcuthz = lowcut * observationTime - 1;
  let hicuthz = hicut * observationTime - 1;
  let js_wrapped_ff = DSPModule.cwrap("ffilter", "null", ["number", "number", "number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_ff(dataPtr, filter_array.length, lowcuthz, hicuthz);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  dataframe.columns.add(DG.Column.fromFloat32Array(column.name + ' Fourier Filtered (L: ' + lowcut + '; H: ' + hicut +
    ')', new Float32Array(result)));
  DSPModule._free(dataHeap.byteOffset);
}


//top-menu: Tools | DSP Preprocess | Spectral Density
//name: Spectral Density
//input: dataframe dataframe [Input Dataframe]
//input: column column [Column]
//input: double observationTime [Time taken to record observed signal]
export async function spectral_density(dataframe, column, observationTime) {
  await initDSPModule();
  let js_wrapped_sd = DSPModule.cwrap("sdensity", "null", ["number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_sd(dataPtr, filter_array.length);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  let col = new Float32Array(result.slice(0, result.length / 2));
  let newName = dataframe.name + " Densities";
  let tableNames = grok.shell.tableNames;
  let doesExist = tableNames.includes(newName);
  if (doesExist) {
    let t = grok.shell.tableByName(newName);
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Density', col));
  } else {
    let col1 = [];
    for (let i = 0; i < col.length; i++) {
      col1[i] = (i + 1) / observationTime;
    }
    let t = DG.DataFrame.create(col.length);
    t.columns.add(DG.Column.fromList('double', 'Frequency', col1));
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Density', col));
    t.name = newName;
    grok.shell.addTableView(t);
  }
  DSPModule._free(dataHeap.byteOffset);
}

//top-menu: Tools | DSP Preprocess | Subsample
//name: Subsample
//input: dataframe dataframe [Input Dataframe]
//input: column column [Column]
//input: double subsampleSize [Desired Subsample Length]
//input: double offset [Subsampling starting point]
export async function subsample(dataframe, column, subsampleSize, offset) {
  await initDSPModule();
  let js_wrapped_ss = DSPModule.cwrap("subsamp", "null", ["number", "number", "number", "number"]);
  //let column_name = col.name + ' Subsample';
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  js_wrapped_ss(dataPtr, filter_array.length, subsampleSize, offset);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  let col = new Float32Array(result.slice(0, subsampleSize));
  DSPModule._free(dataHeap.byteOffset);
  let newName = dataframe.name + " Subsamples, L=" + subsampleSize;
  let tableNames = grok.shell.tableNames;
  let doesExist = tableNames.includes(newName);
  if (doesExist) {
    let t = grok.shell.tableByName(newName);
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Subsample', col));
  } else {
    let col1 = [];
    for (let i = 0; i < col.length; i++) {
      col1[i] = i;
    }
    let t = DG.DataFrame.create(col.length);
    t.columns.add(DG.Column.fromList('int', 'time', col1));
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Subsample', col));
    t.name = newName;
    grok.shell.addTableView(t);
  }
}

//top-menu: Tools | DSP Preprocess | Averaging Downsampling
//name: Averaging Downsampling
//input: dataframe dataframe [Input Dataframe]
//input: column col [Column]
//input: double windowSize [Desired Window Size]
//input: double offset [Subsampling starting point]
export async function asample(dataframe, column, windowSize, offset) {
  await initDSPModule();
  let js_wrapped_as = DSPModule.cwrap("asamp", "null", ["number", "number", "number", "number"]);
  let filter_array = column.getRawData();
  let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
  let dataPtr = DSPModule._malloc(nDataBytes);
  let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(filter_array.buffer));
  let nLength = js_wrapped_as(dataPtr, filter_array.length, windowSize, offset);
  let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
  let col = new Float32Array(result.slice(0, nLength));
  DSPModule._free(dataHeap.byteOffset);
  let newName = dataframe.name + " Avaraged, L=" + nLength;
  let tableNames = grok.shell.tableNames;
  let doesExist = tableNames.includes(newName);
  if (doesExist) {
    let t = grok.shell.tableByName(newName);
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Subsample', col));
  } else {
    let col1 = [];
    for (let i = 0; i < col.length; i++) {
      col1[i] = i;
    }
    let t = DG.DataFrame.create(col.length);
    t.columns.add(DG.Column.fromList('int', 'time', col1));
    t.columns.add(DG.Column.fromFloat32Array(column.name + ' Subsample', col));
    t.name = newName;
    grok.shell.addTableView(t);
  }
}