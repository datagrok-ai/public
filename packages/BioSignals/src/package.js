/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

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

function getMethodsAccordingTo(signalType) {
  let commonFilters = ['IIR', 'FIR', 'normalize', 'resample', 'KalmanFilter', 'ImputeNAN', 'RemoveSpikes', 'ConvolutionalFilter'];
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

function showTheRestOfLayout(formView, accordionCharts, accordionFilters) {
  let rightView = ui.div([ui.h2('Charts'),accordionCharts],'chartview');
  let view = ui.splitH([formView,rightView]);
  ui.dialog('Demo Pipeline')
      .add(view)
      .showModal(true);
  $('.chartview').css('width', '100%');
  $(accordionFilters).css('background', '#fefefe');
  $('.chartview').after('<style>.chart-box{width:100%;height:300px;}</style>');
}

//name: BioSignals
//tags: panel, widgets
//input: dataframe table
//output: widget result
//condition: analysisCondition(table)
export function Biosensors(table) {

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
  }

  function getDescription(i, filtersLST, allParams) {
    let j = filtersLST.length - 1;
    let a = filtersLST[j].value;
    Object.keys(allParams[j]).forEach(key => {
      a = a + ', ' + key + ': ' + allParams[j][key].value;
    });
    return 'Output of Filter №' + i + ': ' + a + '.';
  }

  let tempButton = ui.div();
  tempButton.appendChild(ui.button('launch', () => {

    let accordionFilters = ui.accordion();
    let accordionCharts = ui.accordion();
    let signalInputs = {};

    let column = ui.columnsInput('Columns', table);
    column.setTooltip('Choose columns to plot');

    let samplingFreq = ui.floatInput('Sampling frequency', '');
    samplingFreq.setTooltip('Number of samples taken per second');

    column.onChanged(() => {
      signalInputs = {[column.value[0]]: column};
      accordionCharts.addPane('Raw signal', () => ui.divV([
        ui.div([DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})})],
            'chart-box'
        )]),true
      );
    });

    let signalType = ui.choiceInput('Signal type', '', ['ECG', 'EDA', 'Accelerometer', 'EMG', 'EEG', 'ABP', 'BVP(PPG)', 'Respiration']);

    signalType.onChanged(() => {

      let relevantMethods = getMethodsAccordingTo(signalType.stringValue);

      // Filter dialogue
      let paramsT;
      let filtersList = [];
      let inputsList = [];
      let paramsList = [];
      let containerList = [];
      let addFilterButton = ui.div();
      let filterInputsNew = ui.inputs(filtersList);
      let nameOfLastOutput = '';
      let i = 0;
      addFilterButton.appendChild(ui.button('Add Filter', () => {
        let containerFilter = ui.div();
        containerList[i] = containerFilter;
        filtersList[i] = ui.choiceInput('Filter №' + (i + 1), '', relevantMethods.filters);
        let inputPreset = (Object.keys(signalInputs).length === 1) ? column.value[0] : nameOfLastOutput;
        inputsList[i] = ui.choiceInput('Input', inputPreset, Object.keys(signalInputs));
        let filterInputsOld = ui.inputs([filtersList[i]]);
        containerList[i].appendChild(filterInputsOld);
        filtersList[i].onChanged(function () {
          let val = filtersList[i - 1].value;
          paramsList[i - 1] = paramSelector(val);
          filterInputsNew = ui.inputs(
              [filtersList[i - 1]]
              .concat([inputsList[i - 1]])
              .concat(Object.values(paramsList[i - 1]))
              .concat(addChartButton)
          );
          containerList[i - 1].replaceChild(filterInputsNew, filterInputsOld);
          filterInputsOld = filterInputsNew;
        });
        accordionFilters.addPane('Filter №' + (i + 1), () => containerFilter, true)
        i++;
      }));

      let addChartButton = ui.div();
      addChartButton.appendChild(ui.button('Plot', async () => {
        paramsT = paramsToTable(filtersList, paramsList);
        let t = (inputsList.length === 1) ? DG.DataFrame.fromColumns([column.value[0]]) :
                                            DG.DataFrame.fromColumns([signalInputs[inputsList[i-1].value].columns.byName('sig')]);
        let plotFL = await applyFilter(t, samplingFreq.value, signalType.stringValue, paramsT);
        let name = getDescription(i, filtersList, paramsList);
        accordionCharts.addPane(name, () => ui.divV([
          ui.div([DG.Viewer.fromType('Line chart', plotFL).root], 'chart-box')]),true
        );
        nameOfLastOutput = 'Output of Filter №' + i + ' (' + filtersList[i-1].value + ')';
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
        accordionCharts.addPane(typesOfEstimators.value, () => ui.divV([
          ui.div([DG.Viewer.fromType('Line chart', plotInfo).root], 'chart-box')]),true
        );
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
        accordionCharts.addPane(indicator.value, () => ui.divV([
          ui.div([DG.Viewer.fromType('Line chart', indicatorDf).root], 'chart-box')]),true
        );
        nameOfLastOutput = 'Output of Estimator: ' + indicator.value;
        Object.assign(signalInputs, {[nameOfLastOutput]: indicatorDf});
      }));

      let formView = ui.divV([
        ui.inputs([
          ui.h2('Filtering and Preprocessing'),
          ui.divH([column, samplingFreq]),
          signalType
        ]),
        accordionFilters,
        addFilterButton,
        ui.h2('Information extraction'),
        ui.divH([containerWithEstimators, containerWithPlotsOfEstimators]),
        ui.h2('Physiological Indicators'),
        ui.divH([containerIndicator, calculateButton])
      ],'formview');
      showTheRestOfLayout(formView, accordionCharts, accordionFilters);
    });
    let formView = ui.divV([
      ui.inputs([
        ui.h2('Filtering and Preprocessing'),
        ui.divH([column, samplingFreq]),
        signalType
      ]),
    ],'formview');
    showTheRestOfLayout(formView, accordionCharts, accordionFilters);
  }));
  return new DG.Widget(tempButton);
}
