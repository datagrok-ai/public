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
    let a = '';
    let j = filtersLST.length - 1;
    a = a + filtersLST[j].value;
    Object.keys(allParams[j]).forEach(key => {
      a = a + ', ' + key + ': ' + allParams[j][key].value;
    });
    return 'Output of Filter №' + i + ': ' + a + '.';
  }

  let tempButton = ui.div();
  tempButton.appendChild(ui.button('launch', () => {

    let accordionFilters = ui.accordion('keyThatGivesPersistence');
    let accordionCharts = ui.accordion();

    let column = ui.columnsInput('Columns', table);
    column.setTooltip('Choose columns to plot');

    let samplingFreq = ui.floatInput('Sampling frequency', '');
    samplingFreq.setTooltip('Number of samples taken per second');

    let bsColumn;
    let bsType;
    let npeaks = 10;
    column.onChanged(async () => {
      //bsColumn = column.value[0]; //table.columns.byName('ecg_data');
      //bsColumn = bsColumn.getRawData().slice(0, npeaks * samplingFreq.value);
      //let t = DG.DataFrame.fromColumns([DG.Column.fromList('double', 'x', bsColumn)]);
      //bsType = await typeDetector(t, npeaks, samplingFreq.value);
      bsType = 'ecg';
      accordionCharts.addPane('Raw signal', () => ui.divV([
        ui.div([DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})})],
            'chart-box'
        )]),true
      );
    });

    // Filter dialogue
    let filtersLST = [];
    let allParams = [];
    let paramsT;
    let containerFILTER = ui.div();
    let paramsContainer = ui.div();
    let filterInputs = ui.inputs(filtersLST);
    containerFILTER.appendChild(filterInputs);
    let i = 0;
    let addFilterButton = ui.div();
    addFilterButton.appendChild(ui.bigButton('Add filter', async () => {

      filtersLST[i] = ui.choiceInput('Filter №' + (i + 1), '',
          ['IIR', 'FIR', 'normalize', 'resample', 'KalmanFilter', 'ImputeNAN', 'RemoveSpikes', 'DenoiseEDA', 'ConvolutionalFilter']
      );
      let filterInputs1 = ui.inputs(filtersLST);

      containerFILTER.replaceChild(filterInputs1, filterInputs);
      filterInputs = filterInputs1;

      filtersLST[i].onChanged(function () {
        $(paramsContainer).empty();
        let val = filtersLST[i - 1].value;
        allParams[i - 1] = paramSelector(val);
        paramsContainer.appendChild(ui.inputs(Object.values(allParams[i - 1])));
      });

      accordionFilters.addPane('Filter №' + (i + 1), () => ui.inputs(
          [filtersLST[i], paramsContainer, addChartButton]), true
      );
      i++;
    }));

    // Button to plot filtered signal
    let addChartButton = ui.bigButton('Plot', async () => {
      paramsT = paramsToTable(filtersLST, allParams);
      let t = DG.DataFrame.fromColumns([column.value[0]]);
      let plotFL = await applyFilter(t, samplingFreq.value, bsType, paramsT);
      let name = getDescription(i, filtersLST, allParams);
      accordionCharts.addPane(name, () => ui.divV([
        ui.div([DG.Viewer.fromType('Line chart', plotFL).root], 'chart-box')]),true
      );
    });

    // Information extraction dialogue
    let containerINFO = ui.div();
    let containerINFplot = ui.div();
    let infoType = ui.choiceInput('To extract', '', ['Beat from ECG', 'Phasic estimation']);
    let infoInputs = ui.inputs([infoType]);
    containerINFO.appendChild(infoInputs);
    containerINFplot.appendChild(ui.bigButton('Extract Info', async () => {
      paramsT = paramsToTable(filtersLST, allParams);
      let t = DG.DataFrame.fromColumns([column.value[0]]);
      let plotInfo = await extractInfo(t, samplingFreq.value, bsType, paramsT, infoType);
      accordionCharts.addPane(infoType.value, () => ui.divV([
        ui.div([DG.Viewer.fromType('Line chart', plotInfo).root], 'chart-box')]),true
      );
    }));


    // Indicators dialogue
    let containerIndicator = ui.div();
    let calculateButton = ui.div();
    let indicator = ui.choiceInput('Indicator preset', '', ['HRV']);
    let indicatorInputs = ui.inputs([indicator]);
    containerIndicator.appendChild(indicatorInputs);
    calculateButton.appendChild(ui.bigButton('Calculate', async () => {
      paramsT = paramsToTable(filtersLST, allParams);
      let t = DG.DataFrame.fromColumns([column.value[0]]);
      let indicatorDf = await toIndicators(t, samplingFreq.value, bsType, paramsT, infoType, indicator);
      accordionCharts.addPane(indicator.value, () => ui.divV([
        ui.div([DG.Viewer.fromType('Line chart', indicatorDf).root], 'chart-box')]),true
      );
    }));

    //modal main view
    let formView = ui.divV([
      ui.inputs([
        ui.h2('Preprocessing'),
        column,
        samplingFreq,
        ui.h2('Filtering'),
        accordionFilters
      ]),
      addFilterButton,
      containerINFO,
      containerINFplot,
      containerIndicator,
      calculateButton
    ],'formview');

    let rightView = ui.div([ui.h2('Charts'),accordionCharts],'chartview');
    let view = ui.splitH([formView,rightView]);
    ui.dialog('Demo Pipeline')
        .add(view)
        .showModal(true);
    $('.chartview').css('width', '100%');
    $(accordionFilters).css('background', '#FEFEFE');
    $('.chartview').after('<style>.chart-box{width:100%;height:300px;}</style>');
  }));
  return new DG.Widget(tempButton);
}
