/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//top-menu: DSP | Filters | SMA Filter
//name: Moving Avarage Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: int window_size [SMA Window Size]
export async function SMA_filter(data, column_to_filter, window_size) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("sma", "null", ["number", "number", "number"]);
    //let column_name = column_to_filter.name + ' SMA Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length, window_size);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' SMA Filtered', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Filters | Exponential Filter
//name: Exponential Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double filter_ratio [Exponential Filter Parameter]
export async function Exp_filter(data, column_to_filter, filter_ratio) {
    await initDSPModule();
    let js_wrapped_exp = DSPModule.cwrap("exps", "null", ["number", "number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_exp(dataPtr, filter_array.length, filter_ratio);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Exponentially Filtered', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Filters | Kalmam Filter
//name: Kalman Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double Q [Covariance of the process noise]
//input: double R [Covariance of the observation noise]
//input: double P [a posteriori estimate covariance]
export async function Kalman_filter(data, column_to_filter, Q, R, P) {
    await initDSPModule();
    let js_wrapped_kalm = DSPModule.cwrap("kalman", "null", ["number", "number", "number", "number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_kalm(dataPtr, filter_array.length, Q, R, P)
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Kalman Filtered', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Min-Max
//name: Min Max Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export async function MinMax_transform(data, column_to_filter) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("minmax", "null", ["number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Min Max Normalized', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Z-score
//name: Z-score Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export async function Zscore_transform(data, column_to_filter) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("zscore", "null", ["number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Z-score Normalized', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Box Cox
//name: Box Cox Transform
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lambda
//input: double ofset
export async function box_cox_transform(data, column_to_filter, lambda, ofset) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("boxcox", "null", ["number", "number", "number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length, lambda, ofset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Box Cox Transformed', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Get Trend
//name: Get Trend
//input: dataframe data [Input Dataframe]
//input: column col [Column]
export async function get_trend(data, col) {
    await initDSPModule();
    let js_wrapped_gt = DSPModule.cwrap("gettrend", "null", ["number", "number"]);
    let filter_array = col.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_gt(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(col.name + ' Trend', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Detrend
//name: Detrend
//input: dataframe data [Input Dataframe]
//input: column Column [Column]
export async function remove_trend(data, col) {
    await initDSPModule();
    let js_wrapped_rt = DSPModule.cwrap("removetrend", "null", ["number", "number"]);
    let filter_array = col.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_rt(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(col.name + ' Detrended', column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Filters | Fourier Filter
//name: Fourier Filter
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
//input: double lowcut
//input: double hicut
//input: double observationTime[Time taken to record observed signal]
export async function fourier_filter(data, column_to_filter, lowcut, hicut, observationTime) {
    await initDSPModule();
    let lowcuthz = lowcut * observationTime - 1;
    let hicuthz = hicut * observationTime - 1;
    let js_wrapped_ff = DSPModule.cwrap("ffilter", "null", ["number", "number", "number", "number"]);
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_ff(dataPtr, filter_array.length, lowcuthz, hicuthz);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result);
    data.columns.add(DG.Column.fromFloat32Array(column_to_filter.name + ' Fourier Filtered (L: ' + lowcut +
        '; H: ' + hicut + ')', column));
    DSPModule._free(dataHeap.byteOffset);
}


//top-menu: DSP | Preprocess | Spectral Density
//name: Spectral Density
//input: dataframe data [Input Dataframe]
//input: column col [Column]
//input: double observationTime[Time taken to record observed signal]
export async function spectral_density(data, col, observationTime) {
    await initDSPModule();
    let js_wrapped_sd = DSPModule.cwrap("sdensity", "null", ["number", "number"]);
    let filter_array = col.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sd(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result.slice(0, result.length / 2));
    let newName = data.name + " Densities";
    let tableNames = grok.shell.tableNames;
    let doesExist = tableNames.includes(newName);
    if (doesExist) {
        let t = grok.shell.tableByName(newName);
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Density', column));
    }
    else {
        let col1 = [];
        for (let i = 0; i < column.length; i++) {
            col1[i] = (i + 1) / observationTime;
        }
        let t = DG.DataFrame.create(column.length);
        t.columns.add(DG.Column.fromList('double', 'Frequency', col1));
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Density', column));
        t.name = newName;
        grok.shell.addTableView(t);
    }
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Subsample
//name: Subsample
//input: dataframe data [Input Dataframe]
//input: column col [Column]
//input: double subsampleSize [Desired Subsample Length]
//input: double offset [Subsampling starting point]
export async function subsample(data, col, subsampleSize, offset) {
    await initDSPModule();
    let js_wrapped_ss = DSPModule.cwrap("subsamp", "null", ["number", "number", "number", "number"]);
    //let column_name = col.name + ' Subsample';
    let filter_array = col.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_ss(dataPtr, filter_array.length, subsampleSize, offset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result.slice(0, subsampleSize));
    DSPModule._free(dataHeap.byteOffset);
    let newName = data.name + " Subsamples, L=" + subsampleSize;
    let tableNames = grok.shell.tableNames;
    let doesExist = tableNames.includes(newName);
    if (doesExist) {
        let t = grok.shell.tableByName(newName);
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Subsample', column));
    }
    else {
        let col1 = [];
        for (let i = 0; i < column.length; i++) {
            col1[i] = i;
        }
        let t = DG.DataFrame.create(column.length);
        t.columns.add(DG.Column.fromList('int', 'time', col1));
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Subsample', column));
        t.name = newName;
        grok.shell.addTableView(t);
    }
}

//top-menu: DSP | Preprocess | Averaging Downsampling
//name: Averaging Downsampling
//input: dataframe data [Input Dataframe]
//input: column col [Column]
//input: double windowSize [Desired Window Size]
//input: double offset [Subsampling starting point]
export async function asample(data, col, windowSize, offset) {
    await initDSPModule();
    let js_wrapped_as = DSPModule.cwrap("asamp", "null", ["number", "number", "number", "number"]);
    let filter_array = col.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    let nLength = js_wrapped_as(dataPtr, filter_array.length, windowSize, offset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(result.slice(0, nLength));
    DSPModule._free(dataHeap.byteOffset);
    let newName = data.name + " Avaraged, L=" + nLength;
    let tableNames = grok.shell.tableNames;
    let doesExist = tableNames.includes(newName);
    if (doesExist) {
        let t = grok.shell.tableByName(newName);
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Subsample', column));
    }
    else {
        let col1 = [];
        for (let i = 0; i < column.length; i++) {
            col1[i] = i;
        }
        let t = DG.DataFrame.create(column.length);
        t.columns.add(DG.Column.fromList('int', 'time', col1));
        t.columns.add(DG.Column.fromFloat32Array(col.name + ' Subsample', column));
        t.name = newName;
        grok.shell.addTableView(t);
    }
}