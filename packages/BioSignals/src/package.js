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

function showMainDialog(view, tableWithSignals, tableWithSignalsAndAnnotations, signalType, column, samplingFrequency, chosenDatabase, chosenRecord) {

  let inputCase = tableWithSignals.columns.byName('testEcg');

  // Filter dialogue
  let accordionFilters = ui.accordion();
  let filterTypesList = [];
  let filterContainerList = [];
  let filterChartsList = [];
  let addFilterButton = ui.div();
  let filterInputsNew = ui.inputs(filterTypesList);
  let filterOutputsObj = {[inputCase.name]: column};
  const dspPackageFilters = ['DSP:SMA_filter', 'DSP:Exp_filter', 'DSP:Kalman_filter', 'DSP:MinMax_transform',
    'DSP:Zscore_transform', 'DSP:box_cox_transform', 'DSP:get_trend', 'DSP:remove_trend', 'DSP:fourier_filter',
    'DSP:spectral_density', 'DSP:subsample', 'DSP:asample'];
  let filterScripts = [];
  (async function () {
    for (let filter of dspPackageFilters) filterScripts.push(await grok.functions.eval(filter));
  })();
  let i = 0;
  addFilterButton.appendChild(ui.button('Add Filter', async () => {
    filterChartsList[i] = getEmptyChart();
    let tag = await grok.dapi.scripts.filter('#filters').list();
    filterTypesList[i] = ui.choiceInput('Filter ' + (i + 1), '', filterScripts.concat(tag), async function () {
      let call = filterTypesList[i - 1].value.prepare();
      let t;
      if (i === 1)
        t = DG.DataFrame.fromColumns([inputCase]);
      else if (filterOutputsObj[Object.keys(filterOutputsObj)[i - 1]].columns.byName('sig'))
        t = DG.DataFrame.fromColumns([filterOutputsObj[Object.keys(filterOutputsObj)[i - 1]].columns.byName('sig')]);
      else
        t = DG.DataFrame.fromColumns([filterOutputsObj[Object.keys(filterOutputsObj)[i - 1]].columns.byIndex(i - 1)]);
      let context = DG.Context.create();
      t.name = tableWithSignals.name;
      context.setVariable('table', t);
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
                if (df == null)
                  df = DG.DataFrame.fromColumns([
                    DG.Column.fromList('int', 'time', Array(t.columns.byIndex(0).length).fill().map((_, idx) => idx)),
                    t.columns.byIndex(t.columns.length - 1)
                  ]);
                filterChartsList[i - 1].dataFrame = df;
                Object.assign(filterOutputsObj, {[Object.keys(filterOutputsObj)[i]]: df});
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
    accordionFilters.addPane('Filter ' + (i + 1), () => filterContainerList[i], true)
    i++;
  }));

  // Information extraction dialogue
  let accordionExtractors = ui.accordion();
  let extractorOutputsObj = {};
  let extractorTypesList = [];
  let extractorChartsList = [];
  let extractorContainerList = [];
  let addExtractorButton = ui.div();
  let extractorInputsNew = ui.inputs(extractorTypesList);
  let j = 0;
  addExtractorButton.appendChild(ui.button('Add Extractor', async () => {
    let containerExtractor = ui.div();
    extractorContainerList[j] = containerExtractor;
    extractorChartsList[j] = getEmptyChart();
    let extractors = await grok.dapi.scripts.filter('#extractors').list();
    extractorTypesList[j] = ui.choiceInput('Extractor ' + (j + 1), '', extractors);
    let extractorInputsOld = ui.inputs([extractorTypesList[j]]);
    extractorContainerList[j].appendChild(extractorInputsOld);
    extractorTypesList[j].onChanged(async function () {
      let call = extractorTypesList[j - 1].value.prepare();
      let context = DG.Context.create();
      let t = DG.DataFrame.fromColumns([filterOutputsObj[Object.keys(filterOutputsObj)[i - 1]].columns.byName('sig')]);
      t.name = 'dataframe';
      context.setVariable('table', t);
      call.context = context;
      extractorInputsNew = ui.div([
        ui.block25([
          await call.getEditor(),
          ui.inputs([
            ui.buttonsInput([
              ui.button('Plot', async () => {
                let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting extractor\'s output...');
                try {
                  await call.call();
                  let df = call.getOutputParamValue();
                  if (df == null) {
                    df = t;
                  }
                  extractorChartsList[j - 1].dataFrame = df;
                  let nameOfLastExtractorsOutput = 'Output of Extractor ' + j + ' (' + extractorTypesList[j - 1].value + ')';
                  Object.assign(extractorOutputsObj, {[nameOfLastExtractorsOutput]: df});
                } catch (e) {
                  grok.shell.error(e);
                  throw e;
                } finally {
                  pi.close();
                }
              })
            ])
          ])
        ]),
        ui.block75([extractorChartsList[j - 1]])
      ]);
      extractorContainerList[j - 1].replaceChild(extractorInputsNew, extractorInputsOld);
      extractorInputsOld = extractorInputsNew;
    });
    accordionExtractors.addPane('Extractor ' + (j + 1), () => containerExtractor, true)
    j++;
  }));

  // Indicators dialogue
  let accordionIndicators = ui.accordion();
  let indicatorTypesList = [];
  let indicatorChartsList = [];
  let indicatorInputsList = [];
  let indicatorContainerList = [];
  let addIndicatorButton = ui.div();
  let indicatorInputsNew = ui.inputs(indicatorTypesList);
  let k = 0;
  addIndicatorButton.appendChild(ui.button('Add Indicator', async () => {
    let containerIndicator = ui.div();
    indicatorContainerList[k] = containerIndicator;
    let indicatorInputPreset = Object.keys(extractorOutputsObj)[Object.keys(extractorOutputsObj).length - 1];
    indicatorInputsList[k] = ui.choiceInput('Input', indicatorInputPreset, Object.keys(extractorOutputsObj));
    indicatorChartsList[k] = getEmptyChart();
    let indicatorsNames = await grok.dapi.scripts.filter('#indicators').list();
    indicatorTypesList[k] = ui.choiceInput('Indicator ' + (k + 1), '', indicatorsNames);
    let indicatorInputsOld = ui.inputs([indicatorTypesList[k]]);
    indicatorContainerList[k].appendChild(indicatorInputsOld);
    indicatorTypesList[k].onChanged(async function () {
      let call = indicatorTypesList[k - 1].value.prepare();
      let context = DG.Context.create();
      let t = DG.DataFrame.fromColumns([extractorOutputsObj[indicatorInputsList[k - 1].value].columns.byName('RR intervals')]);
      t.name = 'dataframe';
      context.setVariable('table', t);
      call.context = context;
      indicatorInputsNew = ui.div([
        ui.block25([
          await call.getEditor(),
          ui.inputs([
            ui.buttonsInput([
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
          ])
        ]),
        ui.block75([indicatorChartsList[k - 1]])
      ]);
      indicatorContainerList[k - 1].replaceChild(indicatorInputsNew, indicatorInputsOld);
      indicatorInputsOld = indicatorInputsNew;
    });
    accordionIndicators.addPane('Indicator ' + (k + 1), () => containerIndicator, true)
    k++;
  }));

  let formView = ui.div([
    chosenDatabase,
    chosenRecord,
    ui.divText('Input sampling frequency: ' + samplingFrequency + ' samples per second (Hz)'),
    ui.divText('Signal type: ' + signalType),
    ui.block([
      DG.Viewer.fromType('AnnotatorViewer', tableWithSignalsAndAnnotations)
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
  view = grok.shell.newView('BioSignals', []);
  view.append(formView);
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

async function readPhysionetAnnotations(fileInfos, fileNameWithoutExtension) {
  let f = await grok.functions.eval("BioSignals:readPhysionetAnnotations");
  let call = f.prepare({
    'fileATR': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.atr')),
    'fileDAT': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.dat')),
    'fileHEA': fileInfos.find((({name}) => name === fileNameWithoutExtension + '.hea')),
    'record_name_without_extension': fileNameWithoutExtension
  });
  await call.call();
  const df = call.getParamValue('df');
  const age = call.getParamValue('age');
  const sex = call.getParamValue('sex');
  const dateOfRecording = call.getParamValue('date_of_recording');
  const samplingFrequency = call.getParamValue('sampling_frequency');
  const heartRate = call.getParamValue('heart_rate');
  const rrStd = call.getParamValue('rr_std');
  return [df, age, sex, dateOfRecording, samplingFrequency, heartRate, rrStd];
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
  let view = grok.shell.newView('BioSignals', []);

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let chosenDatabase = ui.choiceInput('Physionet database', '', Object.keys(physionetDatabasesDictionary), () => {
    let chosenRecord = ui.choiceInput('Physionet record', '', physionetDatabasesDictionary[chosenDatabase.stringValue].record_names);

    view = grok.shell.newView('BioSignals', []);
    view.append(ui.divV([chosenDatabase.root, chosenRecord.root]));

    chosenRecord.onInput(async () => {
      let pi = DG.TaskBarProgressIndicator.create('Loading record with annotations from Physionet...');
      let chosenDatabaseShortName = physionetDatabasesDictionary[chosenDatabase.stringValue].short_name;
      let [tableWithSignals, tableWithAnnotations, samplingFrequency] = await loadPhysionetRecordWithAnnotations(chosenDatabaseShortName, chosenRecord);
      let col = tableWithSignals.columns.byName('testEcg');
      let tableWithSignalsAndAnnotations = tableWithSignals.append(tableWithAnnotations);
      let signalType = 'ECG';
      await showMainDialog(view, tableWithSignals, tableWithSignalsAndAnnotations, signalType, col, samplingFrequency, chosenDatabase, chosenRecord);
      pi.close();
    });
  });

  let folderName = ui.stringInput('Path to folder', '');
  let runPipelineButton = ui.div();
  runPipelineButton.appendChild(ui.button('Run pipeline', async () => {

    grok.dapi.users.current().then(async (user) => {
      let pi = DG.TaskBarProgressIndicator.create('Calculating table...');

      const pathToFolder = user.login + ':Home/' + folderName.value + '/';
      const personalFoldersInfos = await grok.dapi.files.list(pathToFolder, false, '');
      const personalFoldersNames = personalFoldersInfos.map((folder) => folder.name)
      let subjectsTable = DG.DataFrame.create(personalFoldersNames.length);

      subjectsTable.columns.addNewString('Person');
      subjectsTable.columns.addNewString('Record');
      subjectsTable.columns.addNewString('Sex');
      subjectsTable.columns.addNewInt('Age');
      subjectsTable.columns.addNewString('Date');
      subjectsTable.columns.addNewInt('Heart Rate');
      subjectsTable.columns.addNewFloat('rrStd');

      let sex, age, dateOfRecording, samplingFrequency, annotationsDF, heartRate, rrStd;
      let indexCounter = 0;
      for (const personalFolderName of personalFoldersNames) {
        console.log(personalFolderName);
        let filesInPersonalFolder = await grok.dapi.files.list(pathToFolder + personalFolderName, false, '');
        let uniqueFileNamesWithoutExtension = Array.from(new Set(filesInPersonalFolder.map((file) => file.name.slice(0, -4))));
        for (const fileNameWithoutExtension of uniqueFileNamesWithoutExtension) {
          [annotationsDF, age, sex, dateOfRecording, samplingFrequency, heartRate, rrStd] = await readPhysionetAnnotations(filesInPersonalFolder, fileNameWithoutExtension);
          subjectsTable.columns.byName('Person').set(indexCounter, personalFolderName);
          subjectsTable.columns.byName('Record').set(indexCounter, fileNameWithoutExtension);
          subjectsTable.columns.byName('Sex').set(indexCounter, sex);
          subjectsTable.columns.byName('Age').set(indexCounter, age);
          subjectsTable.columns.byName('Date').set(indexCounter, dateOfRecording);
          subjectsTable.columns.byName('Heart Rate').set(indexCounter, heartRate);
          subjectsTable.columns.byName('rrStd').set(indexCounter, rrStd);
          indexCounter++;
          console.log(fileNameWithoutExtension);
        }
      }
      let view = grok.shell.addTableView(subjectsTable);
      view.boxPlot({x: 'Sex', y: 'Heart Rate'});
      pi.close();
    });
  }));

  let formView = ui.divV([chosenDatabase, folderName, runPipelineButton]);
  view = grok.shell.newView('BioSignals', []);
  view.append(formView);
}