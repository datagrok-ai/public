/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import {physionetDatabasesDictionary} from "./physionetDatabasesDictionary.js";
import {AnnotatorViewer} from "./annotatorViewer.js";

export let _package = new DG.Package();

//name: AnnotatorViewer
//tags: viewer
//output: viewer result
export function annotator() {
  return new AnnotatorViewer();
}

function getEmptyChart() {
  return DG.Viewer.fromType(
    DG.VIEWER.LINE_CHART,
    DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'time', [])
    ])
  );
}

function showMainDialog(signals, signalsWithAnnotations, signalType, samplingFrequency, chosenDatabase, chosenRecord) {

  // Filter dialogue
  let accordionFilters = ui.accordion();
  let filterTypesList = [];
  let filterContainerList = [];
  let filterChartsList = [];
  let filterInputsNew = ui.inputs(filterTypesList);
  const dspPackageFilters = ['DSP:SMA_filter', 'DSP:Exp_filter', 'DSP:Kalman_filter', 'DSP:MinMax_transform',
    'DSP:Zscore_transform', 'DSP:box_cox_transform', 'DSP:get_trend', 'DSP:remove_trend', 'DSP:fourier_filter',
    'DSP:spectral_density', 'DSP:subsample', 'DSP:asample'];
  let filterScripts = [];
  (async function () {
    for (let filter of dspPackageFilters) filterScripts.push(await grok.functions.eval(filter));
  })();
  let i = 0;
  let addFilterButton = ui.button('Add Filter', async () => {
    filterChartsList[i] = getEmptyChart();
    let tag = await grok.dapi.scripts.filter('#filters').list();
    filterTypesList[i] = ui.choiceInput('Filter ' + (i + 1), '', filterScripts.concat(tag), async function () {
      let call = filterTypesList[i - 1].value.prepare();
      let context = DG.Context.create();
      context.setVariable('table', signals);
      call.context = context;
      filterInputsNew = ui.div([
        ui.block25([
          ui.inputs([
            filterTypesList[i - 1],
            await call.getEditor(),
            ui.button('Plot', async () => {
              let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting filter\'s output...');
              try {
                await call.call();
                let df = call.getOutputParamValue();
                if (df == null) {
                  df = DG.DataFrame.fromColumns([
                    DG.Column.fromList('int', 'time', Array(signals.columns.byIndex(0).length).fill().map((_, idx) => idx)),
                    signals.columns.byIndex(signals.columns.length - 1)
                  ]);
                } else {
                  signals = signals.append(df);
                  df = DG.DataFrame.fromColumns([
                    DG.Column.fromList('int', 'time', Array(df.columns.byIndex(df.columns.length - 1).length).fill().map((_, idx) => idx)),
                    df.columns.byIndex(df.columns.length - 1)
                  ]);
                }
                filterChartsList[i - 1].dataFrame = df;
              } catch (e) {
                grok.shell.error(e);
                throw e;
              } finally {
                pi.close();
              }
            })
          ])
        ]),
        ui.block75([filterChartsList[i - 1]])
      ]);
      filterContainerList[i - 1].replaceChild(filterInputsNew, filterInputsOld);
      filterInputsOld = filterInputsNew;
    });
    let filterInputsOld = ui.inputs([filterTypesList[i]]);
    filterContainerList[i] = ui.div(filterInputsOld);
    accordionFilters.addPane('Filter ' + (i + 1), () => filterContainerList[i], true);
    i++;
  });

  // Information extraction dialogue
  let accordionExtractors = ui.accordion();
  let extracted;
  let extractorTypesList = [];
  let extractorChartsList = [];
  let extractorContainerList = [];
  let extractorInputsNew = ui.inputs(extractorTypesList);
  let j = 0;
  let addExtractorButton = ui.button('Add Extractor', async () => {
    extractorChartsList[j] = getEmptyChart();
    let extractors = await grok.dapi.scripts.filter('#extractors').list();
    extractorTypesList[j] = ui.choiceInput('Extractor ' + (j + 1), '', extractors, async function () {
      let call = extractorTypesList[j - 1].value.prepare();
      let context = DG.Context.create();
      context.setVariable('table', signals);
      call.context = context;
      extractorInputsNew = ui.div([
        ui.block25([
          ui.inputs([
            extractorTypesList[j - 1],
            await call.getEditor(),
            ui.button('Plot', async () => {
              let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting extractor\'s output...');
              try {
                await call.call();
                let df = call.getOutputParamValue();
                if (extracted == null) {
                  extracted = DG.DataFrame.fromColumns([df.columns.byIndex(0)]);
                  df = DG.DataFrame.fromColumns([
                    DG.Column.fromList('int', 'time', Array(df.columns.byIndex(0).length).fill().map((_, idx) => idx)),
                    df.columns.byIndex(0)
                  ]);
                } else {
                  let len = df.rowCount;
                  let floats = new Float32Array(len);
                  let oldColunm = df.columns.byIndex(0);
                  for (let i = 0; i < len; i++)
                    floats[i] = oldColunm.get(i);
                  extracted.columns.addNewFloat(df.columns.byIndex(0).name).init(floats);
                  df = DG.DataFrame.fromColumns([
                    DG.Column.fromList('int', 'time', Array(df.columns.byIndex(0).length).fill().map((_, idx) => idx)),
                    df.columns.byIndex(0)
                  ]);
                }
                extractorChartsList[j - 1].dataFrame = df;
              } catch (e) {
                grok.shell.error(e);
                throw e;
              } finally {
                pi.close();
              }
            })
          ])
        ]),
        ui.block75([extractorChartsList[j - 1]])
      ]);
      extractorContainerList[j - 1].replaceChild(extractorInputsNew, extractorInputsOld);
      extractorInputsOld = extractorInputsNew;
    });
    let extractorInputsOld = ui.inputs([extractorTypesList[j]]);
    extractorContainerList[j] = ui.div(extractorInputsOld);
    accordionExtractors.addPane('Extractor ' + (j + 1), () => extractorContainerList[j], true)
    j++;
  });

  // Indicators dialogue
  let accordionIndicators = ui.accordion();
  let indicatorTypesList = [];
  let indicatorChartsList = [];
  let indicatorContainerList = [];
  let indicatorInputsNew = ui.inputs(indicatorTypesList);
  let k = 0;
  let addIndicatorButton = ui.button('Add Indicator', async () => {
    indicatorChartsList[k] = getEmptyChart();
    let indicatorsNames = await grok.dapi.scripts.filter('#indicators').list();
    indicatorTypesList[k] = ui.choiceInput('Indicator ' + (k + 1), '', indicatorsNames, async function () {
      let call = indicatorTypesList[k - 1].value.prepare();
      let context = DG.Context.create();
      extracted.name = 'Extracted';
      context.setVariable('table', extracted);
      call.context = context;
      indicatorInputsNew = ui.div([
        ui.block25([
          ui.inputs([
            indicatorTypesList[k - 1],
            await call.getEditor(),
            ui.button('Plot', async () => {
              let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting indicator\'s output...');
              try {
                await call.call();
                indicatorChartsList[k - 1].dataFrame = call.getOutputParamValue();
              } catch (e) {
                grok.shell.error(e);
                throw e;
              } finally {
                pi.close();
              }
            })
          ])
        ]),
        ui.block75([indicatorChartsList[k - 1]])
      ]);
      indicatorContainerList[k - 1].replaceChild(indicatorInputsNew, indicatorInputsOld);
      indicatorInputsOld = indicatorInputsNew;
    });
    let indicatorInputsOld = ui.inputs([indicatorTypesList[k]]);
    indicatorContainerList[k] = ui.div(indicatorInputsOld);
    accordionIndicators.addPane('Indicator ' + (k + 1), () => indicatorContainerList[k], true)
    k++;
  });

  grok.shell.newView('BioSignals', [
    chosenDatabase,
    chosenRecord,
    ui.divText('Input sampling frequency: ' + samplingFrequency + ' samples per second (Hz)'),
    ui.divText('Signal type: ' + signalType),
    ui.block([
      DG.Viewer.fromType('AnnotatorViewer', signalsWithAnnotations)
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
  ]);
}

async function readPhysionetRecord(fileInfos, fileNameWithoutExtension) {
  return grok.functions.call("BioSignals:readPhysionetRecord",
    {
      'fileATR': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.atr')),
      'fileDAT': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.dat')),
      'fileHEA': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.hea')),
      'record_name_without_extension': fileNameWithoutExtension
    });
}

//tags: fileViewer, fileViewer-atr, fileViewer-dat, fileViewer-hea
//input: file file
//output: view view
export async function bioSignalViewer(file) {
  let view = DG.View.create();
  const currentFolder = file.fullPath.substring(0, file.fullPath.length - file.name.length);
  const isRecursive = false;
  const fileNameWithoutExtension = file.name.substring(0, file.name.length - 3);
  const fileInfos = await grok.dapi.files.list(currentFolder, isRecursive, fileNameWithoutExtension);
  if (fileInfos.length === 3) {
    let pi = DG.TaskBarProgressIndicator.create('Reading Physionet record...');
    const t = await readPhysionetRecord(fileInfos, file.name.slice(0, -4));
    pi.close();
    view.append(ui.block([DG.Viewer.lineChart(t)], 'd4-ngl-viewer'));
  } else {
    view.append(ui.divText('In order to view this Physionet recording you need to have \'.atr\', \'.dat\', \'.hea\' in the same folder!'));
  }
  return view;
}

export async function loadPhysionetRecordWithAnnotations(chosenDatabase, chosenRecord) {
  try {
    let f = await grok.functions.eval("BioSignals:loadPhysionetRecordWithAnnotations");
    let call = f.prepare({
      'chosenDatabase': chosenDatabase,
      'chosenRecord': chosenRecord.stringValue
    });
    await call.call();
    const tableWithAnnotations = call.getParamValue('annotations_df');
    const tableWithSignals = call.getParamValue('signals_df');
    const samplingFrequency = call.getParamValue('sampling_frequency');
    return [tableWithSignals, tableWithAnnotations, samplingFrequency];
  } catch (e) {
    grok.shell.error(e);
    throw e;
  }
}

//name: BioSignals
//tags: app
export function BioSignals() {

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let chosenDatabase = ui.choiceInput('Physionet database', '', Object.keys(physionetDatabasesDictionary), () => {
    let chosenRecord = ui.choiceInput('Physionet record', '', physionetDatabasesDictionary[chosenDatabase.stringValue].record_names, async () => {
      let pi = DG.TaskBarProgressIndicator.create('Loading record with annotations from Physionet...');
      let chosenDatabaseShortName = physionetDatabasesDictionary[chosenDatabase.stringValue].short_name;
      let [tableWithSignals, tableWithAnnotations, samplingFrequency] = await loadPhysionetRecordWithAnnotations(chosenDatabaseShortName, chosenRecord);
      let tableWithSignalsAndAnnotations = tableWithSignals.append(tableWithAnnotations);
      let signalType = 'ECG';
      await showMainDialog(tableWithSignals, tableWithSignalsAndAnnotations, signalType, samplingFrequency, chosenDatabase, chosenRecord);
      pi.close();
    });

    grok.shell.newView('BioSignals', [ui.divV([chosenDatabase.root, chosenRecord.root])]);
  });

  let folderName = ui.stringInput('Path to folder', '');
  let runPipelineButton = ui.div(ui.button('Run pipeline', async () => {grok.shell.info('Coming soon')}));

  grok.shell.newView('BioSignals', [
    ui.divV([
      chosenDatabase,
      folderName,
      runPipelineButton
    ])
  ]);
}