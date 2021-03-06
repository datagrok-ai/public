/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import {getRelevantMethods} from "./getRelevantMethods.js";
import {getFilterParameters} from "./getFilterParameters.js";
import {getExtractorParameters} from "./getExtractorParameters.js";
import {getIndicatorParameters} from "./getIndicatorParameters.js";

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

async function applyFilter(data, fsamp, paramsT) {
  let f = await grok.functions.eval("BioSignals:filters");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'paramsT': paramsT
  });
  await call.call();
  return call.getParamValue('newDf');
}

async function applyExtractor(data, fsamp, paramsT) {
  let f = await grok.functions.eval("BioSignals:extractors");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'paramsT': paramsT
  });
  await call.call();
  return call.getParamValue('newDf');
}

async function getIndicator(data, fsamp, paramsT, infoType, indicator) {
  let f = await grok.functions.eval("BioSignals:indicators");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'paramsT': paramsT,
    'info': infoType[infoType.length-1].value,
    'preset': indicator.value
  });
  await call.call();
  return call.getParamValue('out');
}

function parametersToDataFrame(filtersLST, allParams) {
  let paramsT = DG.DataFrame.create(filtersLST.length);
  paramsT.columns.addNew('type', 'string');
  let string_parameters = ['ftype', 'normMethod', 'kind', 'allnan', 'method', 'irftype'];
  for (let j = 0; j < filtersLST.length; j++) {
    paramsT.columns.byName('type').set(j, filtersLST[j].value);
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
  let accordionExtractors = ui.accordion();
  let accordionIndicators = ui.accordion();
  let extractorOutputsObj = {};

  let relevantMethods = getRelevantMethods(signalType.stringValue);

  // Filter dialogue
  let paramsT;
  let filterTypesList = [];
  let filterInputsList = [];
  let filterParametersList = [];
  let filterContainerList = [];
  let filterChartsList = [];
  let addFilterButton = ui.div();
  let addFilterChartButton = ui.div();
  let filterInputsNew = ui.inputs(filterTypesList);
  let filterOutputsObj = {[column.value[0]]: column};
  let nameOfLastFiltersOutput = '';
  let i = 0;
  addFilterButton.appendChild(ui.button('Add Filter', () => {
    let containerFilter = ui.div();
    filterContainerList[i] = containerFilter;
    let filterInputPreset = (Object.keys(filterOutputsObj).length === 1) ? column.value[0] : nameOfLastFiltersOutput;
    filterInputsList[i] = ui.choiceInput('Input', filterInputPreset, Object.keys(filterOutputsObj));
    filterChartsList[i] = getEmptyChart();
    filterTypesList[i] = ui.choiceInput('Filter ' + (i + 1), '', relevantMethods.filters);
    let filterInputsOld = ui.inputs([filterTypesList[i]]);
    filterContainerList[i].appendChild(filterInputsOld);
    filterTypesList[i].onChanged(function () {
      let filterType = filterTypesList[i - 1].value;
      filterParametersList[i - 1] = getFilterParameters(filterType);
      filterInputsNew = ui.div([
        ui.block25([
          ui.inputs(
              [filterTypesList[i - 1]]
                  .concat([filterInputsList[i - 1]])
                  .concat(Object.values(filterParametersList[i - 1]))
                  .concat(addFilterChartButton)
          )]
        ),
        ui.block75([filterChartsList[i - 1]])
      ]);
      filterContainerList[i - 1].replaceChild(filterInputsNew, filterInputsOld);
      filterInputsOld = filterInputsNew;
    });
    accordionFilters.addPane('Filter ' + (i + 1), () => containerFilter, true)
    i++;
  }));

  addFilterChartButton.appendChild(ui.button('Plot', async () => {

    let t;
    if (filterInputsList.length === 1) {
      t = DG.DataFrame.fromColumns([column.value[0]]);
    }
    else if (filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)) {
      t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)]);
    }
    else {
      t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName('sig')]);
    }

    let plotFL = '';
    let currentlyChosenFilterType = filterTypesList[filterTypesList.length-1].value;
    switch (currentlyChosenFilterType) {
      case 'Moving Average Filter':
        await SMA_filter(t, column.value[0], filterParametersList[i-1].win_len.value);
        nameOfLastFiltersOutput = column.stringValue + ' SMA Filtered';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Exponential Filter':
        await Exp_filter(t, column.value[0], filterParametersList[i-1].filter_ratio.value);
        nameOfLastFiltersOutput = column.stringValue + ' Exponentially Filtered';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Min Max Normalization':
        await MinMax_transform(t, column.value[0]);
        nameOfLastFiltersOutput = column.stringValue + ' Min Max Normalized';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Z-score Normalization':
        await Zscore_transform(t, column.value[0]);
        nameOfLastFiltersOutput = column.stringValue + ' Z-score Normalized';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Box Cox Transform':
        await box_cox_transform(t, column.value[0], filterParametersList[i-1].lambda.value, filterParametersList[i-1].ofset.value);
        nameOfLastFiltersOutput = column.stringValue + ' Box Cox Transformed';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Get Trend':
        await get_trend(t, column.value[0]);
        nameOfLastFiltersOutput = column.stringValue + ' Trend';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Detrend':
        await remove_trend(t, column.value[0]);
        nameOfLastFiltersOutput = column.stringValue + ' Detrended';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Fourier Filter':
        await fourier_filter(t, column.value[0], filterParametersList[i-1].lowcut.value, filterParametersList[i-1].hicut.value, filterParametersList[i-1].observationTime.value);
        nameOfLastFiltersOutput = column.stringValue + ' Fourier Filtered (L: ' + filterParametersList[i-1].lowcut.value + '; H: ' + filterParametersList[i-1].hicut.value + ')';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Spectral Density':
        await spectral_density(t, column.value[0], filterParametersList[i-1].observationTime.value);
        nameOfLastFiltersOutput = column.stringValue + ' Density';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Subsample':
        await subsample(t, column.value[0], filterParametersList[i-1].subsampleSize.value, filterParametersList[i-1].offset.value);
        nameOfLastFiltersOutput = column.stringValue + ' Subsample';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      case 'Averaging Downsampling':
        await asample(t, column.value[0], filterParametersList[i-1].windowSize.value, filterParametersList[i-1].offset.value);
        nameOfLastFiltersOutput = column.stringValue + ' Subsample';
        plotFL = DG.DataFrame.fromColumns([table.columns.byName('time'), t.columns.byName(nameOfLastFiltersOutput)]);
        break;
      default:
        paramsT = parametersToDataFrame(filterTypesList, filterParametersList);
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
        plotFL = await applyFilter(t, samplingFreq.value, paramsT);
    }
    filterChartsList[i - 1].dataFrame = plotFL;
    Object.assign(filterOutputsObj, {[nameOfLastFiltersOutput]: plotFL});
  }));


  // Information extraction dialogue
  let extractorTypesList = [];
  let extractorChartsList = [];
  let extractorInputsList = [];
  let extractorParametersList = [];
  let extractorContainerList = [];
  let addExtractorButton = ui.div();
  let addExtractorChartButton = ui.div();
  let extractorInputsNew = ui.inputs(extractorTypesList);
  let j = 0;
  addExtractorButton.appendChild(ui.button('Add Extractor', () => {
    let containerExtractor = ui.div();
    extractorContainerList[j] = containerExtractor;
    let extractorInputPreset = (Object.keys(filterOutputsObj).length === 1) ? filterInputsList[0].value : nameOfLastFiltersOutput;
    extractorInputsList[j] = ui.choiceInput('Input', extractorInputPreset, Object.keys(filterOutputsObj));
    extractorChartsList[j] = getEmptyChart();
    extractorTypesList[j] = ui.choiceInput('Extractor ' + (j + 1), '', relevantMethods.extractors);
    let extractorInputsOld = ui.inputs([extractorTypesList[j]]);
    extractorContainerList[j].appendChild(extractorInputsOld);
    extractorTypesList[j].onChanged(function () {
      let extractorType = extractorTypesList[j - 1].value;
      extractorParametersList[j - 1] = getExtractorParameters(extractorType);
      extractorInputsNew = ui.div([
        ui.block25([
          ui.inputs(
              [extractorTypesList[j - 1]]
                  .concat([extractorInputsList[j - 1]])
                  .concat(Object.values(extractorParametersList[j - 1]))
                  .concat(addExtractorChartButton)
          )]
        ),
        ui.block75([extractorChartsList[j - 1]])
      ]);
      extractorContainerList[j - 1].replaceChild(extractorInputsNew, extractorInputsOld);
      extractorInputsOld = extractorInputsNew;
    });
    accordionExtractors.addPane('Extractor ' + (j + 1), () => containerExtractor, true)
    j++;
  }));

  let nameOfLastExtractorsOutput = '';
  addExtractorChartButton.appendChild(ui.button('Plot', async () => {
    let extractorParametersDF = parametersToDataFrame(extractorTypesList, extractorParametersList);
    let t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)]);
    let plotInfo = await applyExtractor(t, samplingFreq.value, extractorParametersDF);
    nameOfLastExtractorsOutput = 'Output of Extractor ' + j + ' (' + extractorTypesList[j-1].value + ')';
    extractorChartsList[j-1].dataFrame = plotInfo;
    Object.assign(extractorOutputsObj, {[nameOfLastExtractorsOutput]: plotInfo});
  }));

  // Indicators dialogue
  let indicatorTypesList = [];
  let indicatorChartsList = [];
  let indicatorInputsList = [];
  let indicatorParametersList = [];
  let indicatorContainerList = [];
  let addIndicatorButton = ui.div();
  let addIndicatorChartButton = ui.div();
  let indicatorInputsNew = ui.inputs(indicatorTypesList);
  let k = 0;
  addIndicatorButton.appendChild(ui.button('Add Indicator', () => {
    let containerIndicator = ui.div();
    indicatorContainerList[k] = containerIndicator;
    let indicatorInputPreset = (Object.keys(extractorOutputsObj).length === 1) ? extractorInputsList[0].value : nameOfLastExtractorsOutput;
    indicatorInputsList[k] = ui.choiceInput('Input', indicatorInputPreset, Object.keys(extractorOutputsObj));
    indicatorChartsList[k] = getEmptyChart();
    indicatorTypesList[k] = ui.choiceInput('Indicator ' + (k + 1), '', relevantMethods.indicators);
    let indicatorInputsOld = ui.inputs([indicatorTypesList[k]]);
    indicatorContainerList[k].appendChild(indicatorInputsOld);
    indicatorTypesList[k].onChanged(function () {
      let indicatorType = indicatorTypesList[k - 1].value;
      indicatorParametersList[k - 1] = getIndicatorParameters(indicatorType);
      indicatorInputsNew = ui.div([
        ui.block25([
          ui.inputs(
              [indicatorTypesList[k - 1]]
                  .concat([indicatorInputsList[k - 1]])
                  .concat(Object.values(indicatorParametersList[k - 1]))
                  .concat(addIndicatorChartButton)
          )]
        ),
        ui.block75([indicatorChartsList[k - 1]])
      ]);
      indicatorContainerList[k - 1].replaceChild(indicatorInputsNew, indicatorInputsOld);
      indicatorInputsOld = indicatorInputsNew;
    });
    accordionIndicators.addPane('Indicator ' + (k + 1), () => containerIndicator, true)
    k++;
  }));

  addIndicatorChartButton.appendChild(ui.button('Plot', async () => {
    let indicatorParametersDF = parametersToDataFrame(indicatorTypesList, indicatorParametersList);
    let t = DG.DataFrame.fromColumns([extractorOutputsObj[indicatorInputsList[k-1].value].columns.byName('RR intervals')]);
    let indicatorDf = await getIndicator(t, samplingFreq.value, indicatorParametersDF, indicatorTypesList, indicatorTypesList[k-1]);
    let nameOfLastIndicatorsOutput = 'Output of Indicator ' + k + ' (' + indicatorTypesList[k-1].value + ')';
    indicatorChartsList[k-1].dataFrame = indicatorDf;
    Object.assign(extractorOutputsObj, {[nameOfLastIndicatorsOutput]: indicatorDf});
  }));

  let formView = ui.div([
    ui.block25([
      ui.inputs([
        column,
        samplingFreq,
        signalType
      ])],'formview'),
    ui.block75([
      DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})})
    ]),
    ui.h2('Filtering and Preprocessing'),
    accordionFilters,
    addFilterButton,
    ui.h2('Information extraction'),
    accordionExtractors,
    addExtractorButton,
    ui.h2('Physiological Indicators'),
    accordionIndicators,
    addIndicatorButton
  ],'formview');
  showTheRestOfLayout(formView);
}

function getDescription(i, filtersLST, allParams) {
  let j = filtersLST.length - 1;
  let a = filtersLST[j].value;
  Object.keys(allParams[j]).forEach(key => {
    a = a + ', ' + key + ': ' + allParams[j][key].value;
  });
  return 'Output of Filter ' + i + ': ' + a + '.';
}

async function readPhysionetRecord(file) {
  let f = await grok.functions.eval("BioSignals:readPhysionetRecord");
  let call = f.prepare({
    'file': file,
    'record_name': file.name
  });
  await call.call();
  return call.getParamValue('df');
}

//tags: fileViewer, fileViewer-tar
//input: file file
//output: view view
export async function bioSignalViewer(file) {
  let view = DG.View.create();
  let t = await readPhysionetRecord(file);
  var host = ui.block([DG.Viewer.lineChart(t)], 'd4-ngl-viewer');
  view.append(host);
  return view;
}

//name: BioSignals
//tags: panel, widgets
//input: dataframe table
//output: widget result
//condition: analysisCondition(table)
export function Biosensors(table) {

  let tempButton = ui.div();
  tempButton.appendChild(ui.bigButton('launch', () => {

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

    let formView = ui.dialog('Demo Pipeline')
        .add(ui.inputs([column]))
        .showModal(true);

    let IsSignalTypeDetectedAutomatically = null;
    column.onChanged(() => {

      formView.close();

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