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

function createPipelineObject(pipeline, functionCategory, typesList, parametersList) {
  let c = 0;
  const typeList = typesList.map((choiceInput) => choiceInput.value);
  for (const type of typeList) {
    pipeline[functionCategory][type] = {};
    let keys = Object.keys(parametersList[c]).map(function (key) {
      return key;
    });
    for (const key of keys) {
      pipeline[functionCategory][type][key] = parametersList[c][key].value;
    }
    c++;
  }
  return pipeline;
}

function showMainDialog(view, tableWithSignals, tableWithSignalsAndAnnotations, signalType, column, samplingFrequency, chosenDatabase, chosenRecord) {

  let inputCase = tableWithSignals.columns.byName('testEcg');

  let accordionFilters = ui.accordion();
  let accordionExtractors = ui.accordion();
  let accordionIndicators = ui.accordion();
  let extractorOutputsObj = {};

  // Filter dialogue
  let filterTypesList = [];
  let filterInputsList = [];
  let filterParametersList = [];
  let filterContainerList = [];
  let filterChartsList = [];
  let addFilterButton = ui.div();
  let filterInputsNew = ui.inputs(filterTypesList);
  let filterOutputsObj = {[inputCase.name]: column};
  const dspPackageFilters = ['DSP:SMA_filter', 'DSP:Exp_filter', 'DSP:Kalman_filter', 'DSP:MinMax_transform',
    'DSP:Zscore_transform', 'DSP:box_cox_transform', 'DSP:get_trend', 'DSP:remove_trend', 'DSP:fourier_filter',
    'DSP:spectral_density', 'DSP:subsample', 'DSP:asample'];
  let i = 0;
  addFilterButton.appendChild(ui.button('Add Filter', async () => {
    let containerFilter = ui.div();
    filterContainerList[i] = containerFilter;
    let nameOfLastFiltersOutput = Object.keys(filterOutputsObj)[Object.keys(filterOutputsObj).length - 1];
    let filterInputPreset = (Object.keys(filterOutputsObj).length === 1) ? inputCase.name : nameOfLastFiltersOutput;
    filterInputsList[i] = ui.choiceInput('Input', filterInputPreset, Object.keys(filterOutputsObj));
    filterChartsList[i] = getEmptyChart();
    let filters = await grok.dapi.scripts.filter('#filters').list();
    for (let ind = 0; ind < dspPackageFilters.length; ind++) {
      filters.push(await grok.functions.eval(dspPackageFilters[ind]));
    }
    filterTypesList[i] = ui.choiceInput('Filter ' + (i + 1), '', filters);
    let filterInputsOld = ui.inputs([filterTypesList[i]]);
    filterContainerList[i].appendChild(filterInputsOld);
    filterTypesList[i].onChanged(async function () {
      let call = filterTypesList[i - 1].value.prepare();
      let t;
      if (filterInputsList.length === 1) {
        t = DG.DataFrame.fromColumns([inputCase]);
      } else if (filterOutputsObj[filterInputsList[i - 1].value].columns.byName('sig')) {
        t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i - 1].value].columns.byName('sig')]);
      } else {
        t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i - 1].value].columns.byIndex(i - 1)]);
      }
      let context = DG.Context.create();
      t.name = 'dataframe';
      context.setVariable('table', t);
      call.context = context;
      filterInputsNew = ui.div([
        ui.block25([
          await call.getEditor(),
          ui.inputs([
            ui.buttonsInput([
              ui.button('Plot', async () => {
                let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting filter\'s output...');
                try {
                  await call.call();
                  let df = call.getOutputParamValue();
                  if (df == null) {
                    df = DG.DataFrame.fromColumns([
                      DG.Column.fromList('int', 'time', Array(t.columns.byIndex(0).length).fill().map((_, idx) => idx)),
                      t.columns.byIndex(t.columns.length - 1)
                    ]);
                  }
                  filterChartsList[i - 1].dataFrame = df;
                  Object.assign(filterOutputsObj, {[nameOfLastFiltersOutput]: df});
                } catch (e) {
                  grok.shell.info(e);
                  throw e;
                } finally {
                  pi.close();
                }
              })
            ])
          ])
        ]),
        ui.block75([filterChartsList[i - 1]])
      ]);
      filterContainerList[i - 1].replaceChild(filterInputsNew, filterInputsOld);
      filterInputsOld = filterInputsNew;
    });
    accordionFilters.addPane('Filter ' + (i + 1), () => containerFilter, true)
    i++;
  }));

  // Information extraction dialogue
  let extractorTypesList = [];
  let extractorChartsList = [];
  let extractorInputsList = [];
  let extractorParametersList = [];
  let extractorContainerList = [];
  let addExtractorButton = ui.div();
  let extractorInputsNew = ui.inputs(extractorTypesList);
  let j = 0;
  addExtractorButton.appendChild(ui.button('Add Extractor', async () => {
    let containerExtractor = ui.div();
    extractorContainerList[j] = containerExtractor;
    let extractorInputPreset = Object.keys(filterOutputsObj)[Object.keys(filterOutputsObj).length - 1];
    extractorInputsList[j] = ui.choiceInput('Input', extractorInputPreset, Object.keys(filterOutputsObj));
    extractorChartsList[j] = getEmptyChart();
    let extractors = await grok.dapi.scripts.filter('#extractors').list();
    extractorTypesList[j] = ui.choiceInput('Extractor ' + (j + 1), '', extractors);
    let extractorInputsOld = ui.inputs([extractorTypesList[j]]);
    extractorContainerList[j].appendChild(extractorInputsOld);
    extractorTypesList[j].onChanged(async function () {
      let call = extractorTypesList[j - 1].value.prepare();
      let context = DG.Context.create();
      let t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i - 1].value].columns.byName('sig')]);
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
                  grok.shell.info(e);
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
  let indicatorTypesList = [];
  let indicatorChartsList = [];
  let indicatorInputsList = [];
  let indicatorParametersList = [];
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
                  grok.shell.info(e);
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

  let savePipelineButton = ui.div();
  savePipelineButton.appendChild(ui.button('Save pipeline', () => {

    let pipeline = {'Filters': {}, 'Estimators': {}, 'Indicators': {}};

    pipeline = createPipelineObject(pipeline, 'Filters', filterTypesList, filterParametersList);
    pipeline = createPipelineObject(pipeline, 'Estimators', extractorTypesList, extractorParametersList);
    pipeline = createPipelineObject(pipeline, 'Indicators', indicatorTypesList, indicatorParametersList);

    grok.dapi.users.current().then(async (user) => {
      const pathToFolder = user.login + ':Home/';
      grok.dapi.files.writeAsText(pathToFolder + 'pipeline.txt', JSON.stringify(pipeline));
    });
  }));

  let formView = ui.div([
    chosenDatabase,
    chosenRecord,
    ui.divText('Input sampling frequency: ' + samplingFrequency + ' samples per second (Hz)'),
    ui.divText('Signal type: ' + signalType),
    savePipelineButton,
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

function getDescription(i, filtersLST, allParams) {
  let j = filtersLST.length - 1;
  let a = filtersLST[j].value;
  Object.keys(allParams[j]).forEach(key => {
    a = a + ', ' + key + ': ' + allParams[j][key].value;
  });
  return 'Output of Filter ' + i + ': ' + a + '.';
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

export async function loadPhysionetRecord(chosenDatabase, chosenRecord) {
  let f = await grok.functions.eval("BioSignals:loadPhysionetRecord");
  let call = f.prepare({
    'chosenDatabase': chosenDatabase,
    'chosenRecord': chosenRecord.stringValue
  });
  await call.call();
  let df = call.getParamValue('df');
  let samplingFrequency = call.getParamValue('sampling_frequency');
  return [df, samplingFrequency];
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
    grok.shell.info(e);
    throw e;
  }
}

export async function loadPhysionetAnnotations(chosenDatabase, chosenRecord) {
  let f = await grok.functions.eval("BioSignals:loadPhysionetAnnotations");
  let call = f.prepare({
    'chosenDatabase': chosenDatabase,
    'chosenRecord': chosenRecord.stringValue
  });
  await call.call();
  return call.getParamValue('df');
}

//name: BioSignals
//tags: app
export function BioSignals() {
  let view = grok.shell.newView('BioSignals', []);
  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;
  let chosenDatabase = ui.choiceInput('Physionet database', '', Object.keys(physionetDatabasesDictionary));
  chosenDatabase.onInput(() => {
    let chosenRecord = ui.choiceInput('Physionet record', '', physionetDatabasesDictionary[chosenDatabase.stringValue].record_names);

    let formView = ui.divV([chosenDatabase, chosenRecord]);
    view = grok.shell.newView('BioSignals', []);
    view.append(formView);

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