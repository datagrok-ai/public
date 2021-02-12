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
  return call.getParamValue('FD_HRV_df');
}

function paramsToTable(filtersLST, allParams) {
  let paramsT = DG.DataFrame.create(filtersLST.length);
  paramsT.columns.addNew('filter', 'string');
  for (let j = 0; j < filtersLST.length; j++) {
    paramsT.columns.byName('filter').set(j, filtersLST[j].value);
    Object.keys(allParams[j]).forEach(key => {
      if (!paramsT.columns.names().includes(key)) {
        // definitely needs reworking
        if (key === 'fp' || key === 'fs') {
          paramsT.columns.addNew(key, 'double');
        } else if (key === 'fout') {
          paramsT.columns.addNew(key, 'int');
        } else {
          paramsT.columns.addNew(key, typeof (allParams[j][key].value));
        }
      }
      paramsT.columns.byName(key).set(j, allParams[j][key].value);
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
    let methodparams;
    if (x === 'IIR') {
      let fp = ui.floatInput('fp', 45);
      let fs = ui.floatInput('fs', 50);
      let ftype = ui.choiceInput('ftype', 'ellip', ['ellip']);
      methodparams = {'fp': fp, 'fs': fs, 'ftype': ftype};
    }
    else if (x === 'normalize') {
      let normMethod = ui.choiceInput('norm_method', 'standard', ['standard']);
      methodparams = {'normMethod': normMethod};
    }
    else if (x === 'resample') {
      let fout = ui.intInput('fout', 4096);
      let kind = ui.choiceInput('kind', 'cubic', ['cubic']);
      methodparams = {'fout': fout, 'kind': kind};
    }
    return methodparams;
  }

  let tempButton = ui.div();
  tempButton.appendChild(ui.button('launch', () => {

    //INPUTS
    let column = ui.columnsInput('Columns', table);
    column.setTooltip('Choose columns to plot');

    let samplingFreq = ui.intInput('Sampling frequency', 2048);
    samplingFreq.setTooltip('Number of samples taken per second');

    let containerImport = ui.div();
    containerImport.appendChild(ui.inputs([column, samplingFreq]));

    let bsColumn;
    let bsType;
    let npeaks = 10;
    let fsamp = samplingFreq.value;
    let node1;
    let tableView = grok.shell.getTableView(table.name);
    column.onChanged(async () => {
      let viewer = DG.Viewer.fromType('Line chart', table, {yColumnNames: column.value.map((c) => {return c.name})});
      //node1 = tableView.dockManager.dock(viewer, 'right', null, 'Original plot');
      //bsColumn = column.value[0]; //table.columns.byName('ecg_data');
      //bsColumn = bsColumn.getRawData().slice(0, npeaks * fsamp);
      //let t = DG.DataFrame.fromColumns([DG.Column.fromList('double', 'x', bsColumn)]);
      //bsType = await typeDetector(t, npeaks, fsamp);
      bsType = 'ecg';
      view.append(viewer.root);
    });

    // Filter dialogue
    let filtersLST = [];
    let allParams = [];
    let paramsT;
    let containerFILTER = ui.div();
    let filterButton = ui.div();
    let accFILTER = ui.accordion();
    let paramsContainer = ui.div();
    let containerFLplot = ui.div();
    let filterInputs = ui.inputs(filtersLST);
    containerFILTER.appendChild(filterInputs);
    let i = 0;
    filterButton.appendChild(ui.button('Add Filter', async () => {

      filtersLST[i] = ui.choiceInput('filter №' + (i + 1), '', ['IIR', 'normalize', 'resample']);
      let filterInputs1 = ui.inputs(filtersLST);

      containerFILTER.replaceChild(filterInputs1, filterInputs);
      filterInputs = filterInputs1;

      filtersLST[i].onChanged(function () {
        $(paramsContainer).empty();
        let val = filtersLST[i - 1].value;
        allParams[i - 1] = paramSelector(val);
        paramsContainer.appendChild(ui.inputs(Object.values(allParams[i - 1])));
      })
      i++;
      accFILTER.addPane('parameters', () => paramsContainer);
    }));

    let node2;
    containerFLplot.appendChild(ui.bigButton('Plot Filtered', async () => {

      paramsT = paramsToTable(filtersLST, allParams);
      let t = DG.DataFrame.fromColumns([column.value[0]]);
      let plotFL = await applyFilter(t, fsamp, bsType, paramsT);
      // node2 = tableView.dockManager.dock(plotFL, 'fill', node1, 'Filtered plot');
      let viewer2 = DG.Viewer.fromType('Line chart', plotFL);
      view.append(viewer2.root);
      //node2 = tableView.dockManager.dock(viewer.root, 'fill', node1, 'Filtered plot');
    }));

    // Information extraction dialogue
    let containerINFO = ui.div();
    let containerINFplot = ui.div();
    let infoType = ui.choiceInput('To extract', 'Beat from ECG', ['Beat from ECG', 'Phasic estimation']);
    let infoInputs = ui.inputs([infoType]);
    containerINFO.appendChild(infoInputs);
    let node3;
    containerINFplot.appendChild(ui.bigButton('Extract Info', async () => {

      paramsT = paramsToTable(filtersLST, allParams);
      let t = DG.DataFrame.fromColumns([column.value[0]]);
      let plotInfo = await extractInfo(t, fsamp, bsType, paramsT, infoType);
      if (infoType.value === 'Beat from ECG') {
        let viewer3 = DG.Viewer.fromType('Line chart', plotInfo);
        view.append(viewer3.root);
        //node3 = tableView.dockManager.dock(plotInfo, 'fill', node2, infoType.value);
      } else if (infoType.value === 'Phasic estimation') {
        let newView = grok.shell.addTableView(plotInfo);
        newView.lineChart();
      }

    }));


    // Indicators dialogue
    let containerIndicator = ui.div();
    let indicator = ui.choiceInput('Indicator preset', '', ['HRV']);
    let indicatorInputs = ui.inputs([indicator]);
    containerIndicator.appendChild(indicatorInputs);
    let view = ui.div([containerImport, containerFILTER, accFILTER, filterButton, containerFLplot, containerINFO, containerINFplot])
    //CREATE DIALOGUE
    ui.dialog('Demo Pipeline')
        .add(view)
        // .add(containerIndicator).onOK(async () => {
        //
        //   paramsT = paramsToTable(filtersLST, allParams);
        //   let t = DG.DataFrame.fromColumns([column.value[0]]);
        //   let indicatorDf = await toIndicators(t, fsamp, bsType, paramsT, infoType, indicator);
        //   grok.shell.addTableView(indicatorDf);
        //
        // }).show()
        .showModal(true);
    $(view).css('height','100%');
    $(view).css('overflow', 'scroll');
  }));
  return new DG.Widget(tempButton);
}
