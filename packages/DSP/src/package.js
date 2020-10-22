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
    let column_name = column_to_filter.name + ' SMA Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length, window_size);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
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
    let column_name = column_to_filter.name + ' Exponentially Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    console.log(dataPtr);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_exp(dataPtr, filter_array.length, filter_ratio);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
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
    let column_name = column_to_filter.name + ' Kalman Filtered';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_kalm(dataPtr, filter_array.length, Q, R, P)
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Min-Max
//name: Min Max Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export async function MinMax_transform(data, column_to_filter) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("minmax", "null", ["number", "number"]);
    let column_name = column_to_filter.name + ' Min Max Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Z-score
//name: Z-score Normalization
//input: dataframe data [Input data table]
//input: column column_to_filter [Signal to Filter]
export async function Zscore_transform(data, column_to_filter) {
    await initDSPModule();
    let js_wrapped_sma = DSPModule.cwrap("zscore", "null", ["number", "number"]);
    let column_name = column_to_filter.name + ' Z-score Normalized';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
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
    let column_name = column_to_filter.name + ' Box Cox Transformed';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sma(dataPtr, filter_array.length, lambda, ofset);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
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
    let column_name = column_to_filter.name + ' Fourier Filtered (L: ' + lowcut + '; H: ' + hicut + ')';
    let filter_array = column_to_filter.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_ff(dataPtr, filter_array.length, lowcuthz, hicuthz);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    DSPModule._free(dataHeap.byteOffset);
}


//top-menu: DSP | Preprocess | Spectral Density
//name: Spectral Density
//input: dataframe data [Input Dataframe]
//input: column Column [Column]
//input: double observationTime[Time taken to record observed signal]
export async function spectral_density(data, Column, observationTime) {
    await initDSPModule();
    let js_wrapped_sd = DSPModule.cwrap("sdensity", "null", ["number", "number"]);
    let column_name = Column.name + ' Density';
    let filter_array = Column.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sd(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(Math.floor(filter_array.length / 2));
    for (let i = 0; i < column.length; i++) { column[i] = result[i]; }
    let newName = data.name + " Densities";
    let tableNames = grok.shell.tableNames;
    let doesExist = tableNames.includes(newName);
    if (doesExist) {
        let t = grok.shell.tableByName(newName);
        t.columns.add(DG.Column.fromFloat32Array(column_name, column));
    }
    else {
        let col1 = [];
        for (let i = 0; i < column.length; i++) {
            col1[i] = (i + 1) / observationTime;
        }
        let t = DG.DataFrame.create(column.length);
        t.columns.add(DG.Column.fromList('double', 'Frequency', col1));
        t.columns.add(DG.Column.fromFloat32Array(column_name, column));
        t.name = newName;
        grok.shell.addTableView(t);
    }
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Spectral Amplitude
//name: Spectral Amplitude
//input: dataframe data [Input Dataframe]
//input: column Column [Column]
//input: double observationTime[Time taken to record observed signal]
export async function spectral_amplitude(data, Column, observationTime) {
    await initDSPModule();
    let js_wrapped_sd = DSPModule.cwrap("samplitude", "null", ["number", "number"]);
    let column_name = Column.name + ' Amplitude';
    let filter_array = Column.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_sd(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(Math.floor(filter_array.length / 2));
    for (let i = 0; i < column.length; i++) { column[i] = result[i]; }
    let newName = data.name + " Amplitudes";
    let tableNames = grok.shell.tableNames;
    let doesExist = tableNames.includes(newName);
    if (doesExist) {
        let t = grok.shell.tableByName(newName);
        t.columns.add(DG.Column.fromFloat32Array(column_name, column));
    }
    else {
        let col1 = [];
        for (let i = 0; i < column.length; i++) {
            col1[i] = (i + 1) / observationTime;
        }
        let t = DG.DataFrame.create(column.length);
        t.columns.add(DG.Column.fromList('double', 'Frequency', col1));
        t.columns.add(DG.Column.fromFloat32Array(column_name, column));
        t.name = newName;
        grok.shell.addTableView(t);
    }
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Get Trend
//name: Get Trend
//input: dataframe data [Input Dataframe]
//input: column Column [Column]
export async function get_trend(data, Column) {
    await initDSPModule();
    let js_wrapped_gt = DSPModule.cwrap("gettrend", "null", ["number", "number"]);
    let column_name = Column.name + ' Trend';
    let filter_array = Column.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_gt(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    DSPModule._free(dataHeap.byteOffset);
}

//top-menu: DSP | Preprocess | Detrend
//name: Detrend
//input: dataframe data [Input Dataframe]
//input: column Column [Column]
export async function remove_trend(data, Column) {
    await initDSPModule();
    let js_wrapped_rt = DSPModule.cwrap("removetrend", "null", ["number", "number"]);
    let column_name = Column.name + ' Detrended';
    let filter_array = Column.getRawData();
    let nDataBytes = filter_array.length * filter_array.BYTES_PER_ELEMENT;
    let dataPtr = DSPModule._malloc(nDataBytes);
    let dataHeap = new Uint8Array(DSPModule.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(filter_array.buffer));
    js_wrapped_rt(dataPtr, filter_array.length);
    let result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, filter_array.length);
    let column = new Float32Array(filter_array.length);
    for (let i = 0; i < filter_array.length; i++) { column[i] = result[i]; }
    data.columns.add(DG.Column.fromFloat32Array(column_name, column));
    DSPModule._free(dataHeap.byteOffset);
}