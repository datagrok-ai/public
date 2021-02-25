/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

async function asample(data, col, windowSize, offset) {
  return await grok.functions.call('Dsp:asample',
      {
        'data': data,
        'col': col,
        'windowSize': windowSize,
        'offset': offset
      });
}

async function subsample(data, col, subsampleSize, offset) {
  return await grok.functions.call('Dsp:subsample',
      {
        'data': data,
        'col': col,
        'subsampleSize': subsampleSize,
        'offset': offset
      });
}

async function spectral_density(data, col, observationTime) {
  return await grok.functions.call('Dsp:spectral_density',
      {
        'data': data,
        'col': col,
        'observationTime': observationTime
      });
}

async function fourier_filter(data, column_to_filter, lowcut, hicut, observationTime) {
  return await grok.functions.call('Dsp:fourier_filter',
      {
        'data': data,
        'column_to_filter': column_to_filter,
        'lowcut': lowcut,
        'hicut': hicut,
        'observationTime': observationTime
      });
}

async function remove_trend(data, col) {
  return await grok.functions.call('Dsp:remove_trend',
      {
        'data': data,
        'col': col
      });
}

async function get_trend(data, col) {
  return await grok.functions.call('Dsp:get_trend',
      {
        'data': data,
        'col': col
      });
}

async function box_cox_transform(data, columnToFilter, lambda, ofset) {
  return await grok.functions.call('Dsp:box_cox_transform',
      {
        'data': data,
        'column_to_filter': columnToFilter,
        'lambda': lambda,
        'ofset': ofset
      });
}

async function Zscore_transform(data, columnToFilter) {
  return await grok.functions.call('Dsp:Zscore_transform',
      {
        'data': data,
        'column_to_filter': columnToFilter
      });
}

async function MinMax_transform(data, columnToFilter) {
  return await grok.functions.call('Dsp:MinMax_transform',
      {
        'data': data,
        'column_to_filter': columnToFilter
      });
}

async function Exp_filter(data, columnToFilter, filterRatio) {
  return await grok.functions.call('Dsp:Exp_filter',
      {
        'data': data,
        'column_to_filter': columnToFilter,
        'filter_ratio': filterRatio
      });
}

async function SMA_filter(data, columnToFilter, windowSize) {
  return await grok.functions.call('Dsp:SMA_filter',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'window_size': windowSize
    });
}

async function typeDetector(table, npeaks, fsamp) {
  let bsType = await grok.functions.call('BioSignals:typeDetector',
    {
      'dat': table,
      'npeaks': npeaks,
      'fsamp': fsamp
    })
  return (bsType);
}

async function applyFilter(data, fsamp, bsType, paramsT) {
  let f = await grok.functions.eval("BioSignals:filtersPyphysio");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'signalType': bsType,
    'paramsT': paramsT
  });
  await call.call();
  return call.getParamValue('newDf');
}

async function extractInfo(data, fsamp, bsType, paramsT, infoType) {
  let f = await grok.functions.eval("BioSignals:infoPyphysio");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'signalType': bsType,
    'paramsT': paramsT,
    'info': infoType.value
  });
  await call.call();
  return call.getParamValue('newDf');
}

async function toIndicators(data, fsamp, bsType, paramsT, infoType, indicator) {
  let f = await grok.functions.eval("BioSignals:indicatorsPyphysio");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'signalType': bsType,
    'paramsT': paramsT,
    'info': infoType.value,
    'preset': indicator.value
  });
  await call.call();
  return call.getParamValue('out');
}

function paramsToTable(filtersLST, allParams) {
  let paramsT = DG.DataFrame.create(filtersLST.length);
  paramsT.columns.addNew('filter', 'string');
  let string_parameters = ['ftype', 'normMethod', 'kind', 'allnan', 'method', 'irftype'];
  for (let j = 0; j < filtersLST.length; j++) {
    paramsT.columns.byName('filter').set(j, filtersLST[j].value);
    Object.keys(allParams[j]).forEach(key => {
      if (!paramsT.columns.names().includes(key)) {
        if (string_parameters.includes(key)) {
          paramsT.columns.addNew(key, 'string');
        } else {
          paramsT.columns.addNew(key, 'double');
        }
        paramsT.columns.byName(key).set(j, allParams[j][key].value);
      }
    })
  }
  return paramsT;
}

function getMethodsAccordingTo(signalType, dspMethods) {
  let commonFilters = dspMethods.concat(['IIR', 'FIR', 'normalize', 'resample', 'KalmanFilter', 'ImputeNAN', 'RemoveSpikes', 'ConvolutionalFilter']);
  let commonEstimators = ['Local energy'];
  let commonIndicators = [];
  switch (signalType) {
    case 'ECG':
      return {
        filters: commonFilters,
        estimators: commonEstimators.concat(['Beat from ECG']),
        indicators: commonIndicators.concat(['HRV time domain', 'HRV frequency domain', 'HRV nonlinear domain'])
      };
    case 'EDA':
      return {
        filters: commonFilters.concat(['DenoiseEDA']),
        estimators: commonEstimators.concat(['Phasic estimation']),
        indicators: commonIndicators
      };
    case 'Accelerometer':
      return {
        filters: commonFilters,
        estimators: commonEstimators,
        indicators: commonIndicators
      };
    case 'EMG':
      return {
        filters: commonFilters,
        estimators: commonEstimators,
        indicators: commonIndicators
      };
    case 'EEG':
      return {
        filters: commonFilters,
        estimators: commonEstimators,
        indicators: commonIndicators
      };
    case 'ABP':
      return {
        filters: commonFilters,
        estimators: commonEstimators.concat(['BeatFromBP']),
        indicators: commonIndicators
      };
    case 'BVP(PPG)':
      return {
        filters: commonFilters,
        estimators: commonEstimators.concat(['BeatFromBP']),
        indicators: commonIndicators
      };
    case 'Respiration':
      return {
        filters: commonFilters,
        estimators: commonEstimators,
        indicators: commonIndicators
      };
  }
}

function showTheRestOfLayout(formView) {
  ui.dialog('Demo Pipeline')
      .add(formView)
      .showModal(true);
}

function getEmptyChart() {
  return DG.Viewer.fromType(
      DG.VIEWER.LINE_CHART,
      DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.TYPE.FLOAT, 'time', [])
      ])
  );
}

function showMainDialog(table, signalType, column, samplingFreq) {

  let accordionFilters = ui.accordion();
  let accordionEstimators = ui.accordion();
  let accordionIndicators = ui.accordion();

  let dspMethods = ['SMA_filter', 'Exp_filter', 'MinMax_transform', 'Zscore_transform', 'box_cox_transform', 'get_trend', 'remove_trend', 'fourier_filter', 'spectral_density', 'subsample', 'asample'];
  let relevantMethods = getMethodsAccordingTo(signalType.stringValue, dspMethods);
  let signalInputs = {[column.value[0]]: column};

  // Filter dialogue
  let paramsT;
  let filtersList = [];
  let inputsList = [];
  let paramsList = [];
  let containerList = [];
  let emptyCharts = [];
  let addFilterButton = ui.div();
  let filterInputsNew = ui.inputs(filtersList);
  let nameOfLastOutput = '';
  let i = 0;
  addFilterButton.appendChild(ui.button('Add Filter', () => {
    let containerFilter = ui.div();
    containerList[i] = containerFilter;
    filtersList[i] = ui.choiceInput('Filter ' + (i + 1), '', relevantMethods.filters);
    let inputPreset = (Object.keys(signalInputs).length === 1) ? column.value[0] : nameOfLastOutput;
    inputsList[i] = ui.choiceInput('Input', inputPreset, Object.keys(signalInputs));
    emptyCharts[i] = getEmptyChart();
    let filterInputsOld = ui.inputs([filtersList[i]]);
    containerList[i].appendChild(filterInputsOld);
    filtersList[i].onChanged(function () {
      let val = filtersList[i - 1].value;
      paramsList[i - 1] = paramSelector(val);
      filterInputsNew = ui.div([
        ui.block25([
          ui.inputs(
              [filtersList[i - 1]]
                  .concat([inputsList[i - 1]])
                  .concat(Object.values(paramsList[i - 1]))
                  .concat(addChartButton)
          )]
        ),
        ui.block75([emptyCharts[i - 1]])
      ]);
      containerList[i - 1].replaceChild(filterInputsNew, filterInputsOld);
      filterInputsOld = filterInputsNew;
    });
    accordionFilters.addPane('Filter ' + (i + 1), () => containerFilter, true)
    i++;
  }));

  let addChartButton = ui.div();
  addChartButton.appendChild(ui.button('Plot', async () => {
    paramsT = paramsToTable(filtersList, paramsList);
    let t = (inputsList.length === 1) ? DG.DataFrame.fromColumns([column.value[0]]) :
        DG.DataFrame.fromColumns([signalInputs[inputsList[i-1].value].columns.byName('sig')]);

    let plotFL = '';
    let currentlyChosenFilterType = filtersList[filtersList.length-1].value;
    switch (currentlyChosenFilterType) {
      case 'SMA_filter':
        await SMA_filter(t, column.value[0], paramsList[i-1].win_len.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' SMA Filtered')]);
        break;
      case 'Exp_filter':
        await Exp_filter(t, column.value[0], paramsList[i-1].filter_ratio.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Exponentially Filtered')]);
        break;
      case 'MinMax_transform':
        await MinMax_transform(t, column.value[0]);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Min Max Normalized')]);
        break;
      case 'Zscore_transform':
        await Zscore_transform(t, column.value[0]);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Z-score Normalized')]);
        break;
      case 'box_cox_transform':
        await box_cox_transform(t, column.value[0], paramsList[i-1].lambda.value, paramsList[i-1].ofset.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Box Cox Transformed')]);
        break;
      case 'get_trend':
        await get_trend(t, column.value[0]);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Trend')]);
        break;
      case 'remove_trend':
        await remove_trend(t, column.value[0]);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Detrended')]);
        break;
      case 'fourier_filter':
        await fourier_filter(t, column.value[0], paramsList[i-1].lowcut.value, paramsList[i-1].hicut.value, paramsList[i-1].observationTime.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Fourier Filtered (L: ' + paramsList[i-1].lowcut.value + '; H: ' + paramsList[i-1].hicut.value + ')')]);
        break;
      case 'spectral_density':
        await spectral_density(t, column.value[0], paramsList[i-1].observationTime.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Density')]);
        break;
      case 'subsample':
        await subsample(t, column.value[0], paramsList[i-1].subsampleSize.value, paramsList[i-1].offset.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Subsample')]);
        break;
      case 'asample':
        await asample(t, column.value[0], paramsList[i-1].windowSize.value, paramsList[i-1].offset.value);
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(column.stringValue + ' Subsample')]);
        break;
      default:
        plotFL = await applyFilter(t, samplingFreq.value, signalType.stringValue, paramsT);
    }
    emptyCharts[i - 1].dataFrame = plotFL;
    nameOfLastOutput = 'Output of Filter ' + i + ' (' + filtersList[i-1].value + ')';
    Object.assign(signalInputs, {[nameOfLastOutput]: plotFL});
  }));


  // Information extraction dialogue
  let containerWithEstimators = ui.div();
  let containerWithPlotsOfEstimators = ui.div();
  let typesOfEstimators = ui.choiceInput('Estimators', '', relevantMethods.estimators);
  let inputsToEstimators = ui.inputs([typesOfEstimators]);
  containerWithEstimators.appendChild(inputsToEstimators);
  containerWithPlotsOfEstimators.appendChild(ui.button('Extract Info', async () => {
    paramsT = paramsToTable(filtersList, paramsList);
    let t = DG.DataFrame.fromColumns([column.value[0]]);
    let plotInfo = await extractInfo(t, samplingFreq.value, signalType.stringValue, paramsT, typesOfEstimators);
    let estimatorChart = getEmptyChart();
    accordionEstimators.addPane(typesOfEstimators.value, () => ui.block(estimatorChart),true);
    estimatorChart.dataFrame = plotInfo;
    nameOfLastOutput = 'Output of Estimator: ' + typesOfEstimators.value;
    Object.assign(signalInputs, {[nameOfLastOutput]: plotInfo});
  }));

  // Indicators dialogue
  let containerIndicator = ui.div();
  let calculateButton = ui.div();
  let indicator = ui.choiceInput('Indicators', '', relevantMethods.indicators);
  let indicatorInputs = ui.inputs([indicator]);
  containerIndicator.appendChild(indicatorInputs);
  calculateButton.appendChild(ui.button('Calculate', async () => {
    paramsT = paramsToTable(filtersList, paramsList);
    let t = DG.DataFrame.fromColumns([column.value[0]]);
    let indicatorDf = await toIndicators(t, samplingFreq.value, signalType.stringValue, paramsT, typesOfEstimators, indicator);
    let indicatorChart = getEmptyChart();
    accordionIndicators.addPane(indicator.value, () => ui.block(indicatorChart),true);
    indicatorChart.dataFrame = indicatorDf;
    nameOfLastOutput = 'Output of Estimator: ' + indicator.value;
    Object.assign(signalInputs, {[nameOfLastOutput]: indicatorDf});
  }));

  let formView = ui.div([
    ui.block25([
      ui.inputs([
        ui.h2('Filtering and Preprocessing'),
        column,
        samplingFreq,
        signalType
      ])],'formview'),
    ui.block75([
      DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})})
    ]),
    accordionFilters,
    addFilterButton,
    ui.h2('Information extraction'),
    ui.divH([containerWithEstimators, containerWithPlotsOfEstimators]),
    accordionEstimators,
    ui.h2('Physiological Indicators'),
    ui.divH([containerIndicator, calculateButton]),
    accordionIndicators
  ],'formview');
  showTheRestOfLayout(formView);
}

function paramSelector(x) {
  if (x === 'IIR') {
    let passFrequency = ui.floatInput('Pass frequency', '');
    let stopFrequency = ui.floatInput('Stop frequency', '');
    let ftype = ui.choiceInput('Filter type', '', ['butter', 'cheby1', 'cheby2', 'ellip']);
    return {'fp': passFrequency, 'fs': stopFrequency, 'ftype': ftype};
  }
  else if (x === 'FIR') {
    let passFrequency = ui.floatInput('Pass frequency', '');
    let stopFrequency = ui.floatInput('Stop frequency', '');
    //let ftype = ui.choiceInput('Window type', '', ['hamming']);
    return {'fp': passFrequency, 'fs': stopFrequency};
  }
  else if (x === 'normalize') {
    let normMethod = ui.choiceInput('norm_method', '', ['mean', 'standard', 'min', 'maxmin', 'custom']);
    return {'normMethod': normMethod};
  }
  else if (x === 'resample') {
    let fout = ui.intInput('Output sampling frequency', '');
    let kind = ui.choiceInput('Interpolation method', '', ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']);
    return {'fout': fout, 'kind': kind};
  }
  else if (x === 'KalmanFilter') {
    let r = ui.floatInput('R', '');
    r.setTooltip("R should be positive");
    let ratio = ui.floatInput('Ratio', '');
    ratio.setTooltip("Ratio should be >1");
    return {'R': r, 'ratio': ratio};
  }
  else if (x === 'ImputeNAN') {
    let winLen = ui.floatInput('Window Length', '');
    winLen.setTooltip('Window Length should be positive');
    let allNan = ui.choiceInput('All NaN', '', ['zeros', 'nan']);
    return {'win_len': winLen, 'allnan': allNan};
  }
  else if (x === 'RemoveSpikes') {
    let K = ui.floatInput('K', 2);
    K.setTooltip("K should be positive");
    let N = ui.intInput('N', 1);
    N.setTooltip('N should be positive integer');
    let dilate = ui.floatInput('Dilate', 0);
    dilate.setTooltip('dilate should be >= 0.0');
    let D = ui.floatInput('D', 0.95);
    D.setTooltip('D should be >= 0.0');
    let method = ui.choiceInput('Method', '', ['linear', 'step']);
    return {'K': K, 'N': N, 'dilate': dilate, 'D': D, 'method': method};
  }
  else if (x === 'DenoiseEDA') {
    let winLen = ui.floatInput('Window Length', 2);
    winLen.setTooltip('Window Length should be positive');
    let threshold = ui.floatInput('Threshold', '');
    threshold.setTooltip('Threshold should be positive');
    return {'win_len': winLen, 'threshold': threshold};
  }
  else if (x === 'ConvolutionalFilter') {
    let irftype = ui.choiceInput('irftype', '', ['gauss', 'rect', 'triang', 'dgauss', 'custom']);
    irftype.setTooltip('Impulse response function (IRF)');
    let winLen = ui.floatInput('Window Length', 2);
    winLen.setTooltip('Duration of the generated IRF in seconds (if irftype is not \'custom\')');
    return {'win_len': winLen, 'irftype': irftype};
  }
  else if (x === 'SMA_filter') {
    let winLen = ui.floatInput('Window Length', '');
    winLen.setTooltip('Order of filtering');
    return {'win_len': winLen};
  }
  else if (x === 'Exp_filter') {
    let filterRatio = ui.floatInput('Filter ratio', '');
    return {'filter_ratio': filterRatio};
  }
  else if (x === 'MinMax_transform') {
    return {};
  }
  else if (x === 'Zscore_transform') {
    return {};
  }
  else if (x === 'box_cox_transform') {
    let lambda = ui.floatInput('lambda', '');
    let ofset = ui.floatInput('ofset', '');
    return {'lambda': lambda, 'ofset': ofset};
  }
  else if (x === 'get_trend') {
    return {};
  }
  else if (x === 'remove_trend') {
    return {};
  }
  else if (x === 'fourier_filter') {
    let lowcut = ui.floatInput('lowcut', '');
    let hicut = ui.floatInput('hicut', '');
    let observationTime = ui.floatInput('Observation time', '');
    return {'lowcut': lowcut, 'hicut': hicut, 'observationTime': observationTime};
  }
  else if (x === 'spectral_density') {
    let observationTime = ui.floatInput('Observation time', '');
    return {'observationTime': observationTime};
  }
  else if (x === 'subsample') {
    let subsampleSize = ui.floatInput('Subsample size', '');
    let offset = ui.floatInput('Offset', '');
    return {'subsampleSize': subsampleSize, 'offset': offset};
  }
  else if (x === 'asample') {
    let windowSize = ui.floatInput('Window size', '');
    let offset = ui.floatInput('Offset', '');
    return {'windowSize': windowSize, 'offset': offset};
  }
}

function getDescription(i, filtersLST, allParams) {
  let j = filtersLST.length - 1;
  let a = filtersLST[j].value;
  Object.keys(allParams[j]).forEach(key => {
    a = a + ', ' + key + ': ' + allParams[j][key].value;
  });
  return 'Output of Filter ' + i + ': ' + a + '.';
}

//name: BioSignals
//tags: panel, widgets
//input: dataframe table
//output: widget result
//condition: analysisCondition(table)
export function Biosensors(table) {

  let tempButton = ui.div();
  tempButton.appendChild(ui.button('launch', () => {

    let column = ui.columnsInput('Columns', table);
    column.setTooltip('Choose columns to plot');

    let samplingFreq = ui.floatInput('Sampling frequency', '');
    samplingFreq.setTooltip('Number of samples taken per second');

    let signalType = ui.choiceInput('Signal type', '', ['ECG', 'EDA', 'Accelerometer', 'EMG', 'EEG', 'ABP', 'BVP(PPG)', 'Respiration']);
    signalType.onChanged(() => {
      if (!IsSignalTypeDetectedAutomatically) {
        showMainDialog(table, signalType, column, samplingFreq);
      }
    });

    let formView = ui.div([
      ui.inputs([
        ui.h2('Filtering and Preprocessing'),
        column
      ]),
    ],'formview');
    showTheRestOfLayout(formView);

    let IsSignalTypeDetectedAutomatically = null;
    column.onChanged(() => {
      if (table.col(column.stringValue).semType) {
        IsSignalTypeDetectedAutomatically = true;
        signalType.stringValue = table.col(column.stringValue).semType.split('-')[1];
        showMainDialog(table, signalType, column, samplingFreq);
      }
      else {
        IsSignalTypeDetectedAutomatically = false;
        let formView = ui.div([
            ui.block25([
                ui.inputs([
                    ui.h2('Filtering and Preprocessing'),
                  column,
                  samplingFreq,
                  signalType
                ])
            ],'formview'),
          ui.block75([
            DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})})
          ])
        ]);
        showTheRestOfLayout(formView);
      }
    });
  }));
  return new DG.Widget(tempButton);
}
