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

async function readPhysionetRecord(fileInfos, fileNameWithoutExtension) {
  return grok.functions.call("BioSignals:readPhysionetRecord",
    {
      'fileATR': fileInfos.find(({name}) => name === fileNameWithoutExtension + '.atr'),
      'fileDAT': fileInfos.find(({name}) => name === fileNameWithoutExtension + '.dat'),
      'fileHEA': fileInfos.find(({name}) => name === fileNameWithoutExtension + '.hea'),
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
    const annotations = call.getParamValue('annotations_df');
    const signals = call.getParamValue('signals_df');
    signals.columns.byIndex(0).setTag('samplingFrequency', call.getParamValue('sampling_frequency'))
    return {signals, annotations};
  } catch (e) {
    grok.shell.error(e);
    throw e;
  }
}

async function getInitValues(isLocalTable, chosenDatabase, localTables, chosenRecord) {
  let signals;
  if (isLocalTable) {
    signals = DG.DataFrame.fromColumns([localTables.find(({name}) => name === chosenDatabase.stringValue).columns.byName(chosenRecord.stringValue)]);
    let annotations = DG.DataFrame.create(1);
    signals.name = chosenDatabase.stringValue;
    return {signals, annotations};
  } else {
    let pi = DG.TaskBarProgressIndicator.create('Loading record with annotations from Physionet...');
    let chosenDatabaseShortName = physionetDatabasesDictionary[chosenDatabase.stringValue].shortName;
    let {signals, annotations} = await loadPhysionetRecordWithAnnotations(chosenDatabaseShortName, chosenRecord);
    signals.name = chosenDatabase.stringValue + '/' + chosenRecord.stringValue;
    pi.close();
    return {signals, annotations};
  }
}

//name: BioSignals
//tags: app
export function BioSignals() {

  let link = ui.element('a');
  link.href = 'https://github.com/datagrok-ai/public/tree/master/packages/BioSignals'
  link.text = 'here';
  link.target = '_blank';

  let appDescription = ui.div([
    ui.h1('Analyze biomedical signals using built-in and custom scripts.'),
    ui.span(['See details ', link]),
    ui.divText('\n How to add your script:', {style: {'font-weight': 'bolder'}}),
    ui.divText('\n 1. Go to Functions | Scripts | Actions | New <yourScriptLanguage> Script'),
    ui.divText('2. Write your script, and test it on files'),
    ui.divText('3. Set tag #filters, #extractors or #indicators. Now it is available in corresponding app section')
  ], 'grok-datajob-publish-alert');

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let chosenRecordDiv = ui.div();
  let enterSamplingFrequencyDiv = ui.div();
  let samplingFrequencyDiv = ui.div();
  let annotationViewerDiv = ui.div();

  let localTables = grok.shell.tables;
  let namesOfLocalTables = localTables.map((df) => df.name)
  let tablesAndPhysionetDBs = namesOfLocalTables.concat(Object.keys(physionetDatabasesDictionary));
  let chosenDatabase = ui.choiceInput('Database', '', tablesAndPhysionetDBs, () => {
    let isLocalTable = (namesOfLocalTables.includes(chosenDatabase.stringValue));
    let items = (isLocalTable) ?
      localTables.find(({name}) => name === chosenDatabase.stringValue).columns.names() :
      physionetDatabasesDictionary[chosenDatabase.stringValue].namesOfRecords;
    let samplingFreq = ui.floatInput('Sampling frequency: ', '');
    if (isLocalTable) enterSamplingFrequencyDiv.append(samplingFreq.root);
    let columnName = (isLocalTable) ? 'Column' : 'Physionet record';
    let chosenRecord = ui.choiceInput(columnName, '', items, async () => {
      let {signals, annotations} = await getInitValues(isLocalTable, chosenDatabase, localTables, chosenRecord);
      let signalsWithAnnotations = signals.append(annotations);

      if (isLocalTable) {
        enterSamplingFrequencyDiv.innerHTML = '';
        signals.columns.byIndex(0).setTag('samplingFrequency', samplingFreq.value);
        signalsWithAnnotations.columns.byIndex(0).setTag('samplingFrequency', samplingFreq.value);
      }
      annotationViewerDiv.innerHTML = '';
      annotationViewerDiv.append(
        ui.block([
          DG.Viewer.fromType('AnnotatorViewer', signalsWithAnnotations)
        ])
      );

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
      await (async function () {
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
          for (let table of grok.shell.tables)
            if (table.name !== chosenDatabase.stringValue)
              context.setVariable(table.name, table);
          call.context = context;
          filterInputsNew = ui.div([
            ui.block25([
              ui.inputs([
                filterTypesList[i - 1],
                await call.getEditor(),
                ui.buttonsInput([
                  ui.divH([
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
                    }),
                    ui.iconFA('trash', () => {
                      accordionFilters.root.removeChild(accordionFilters.root.childNodes[i - 1]);
                      filterContainerList = filterContainerList.filter(e => e !== filterContainerList[i - 1]);
                      filterTypesList = filterTypesList.filter(e => e !== filterTypesList[i - 1]);
                      filterChartsList = filterChartsList.filter(e => e !== filterChartsList[i - 1]);
                      signals.columns.remove(signals.columns.byIndex(i).name);
                      i = i - 1;
                    })
                  ])
                ])
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
          grok.shell.tables.forEach(function (table, index) {
            context.setVariable('table' + index, table);
          });
          call.context = context;
          extractorInputsNew = ui.div([
            ui.block25([
              ui.inputs([
                extractorTypesList[j - 1],
                await call.getEditor(),
                ui.buttonsInput([
                  ui.divH([
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
                          let floats = new Array(extracted.columns.byIndex(extracted.columns.length - 1).length);
                          let oldColumn = df.columns.byIndex(0);
                          for (let i = 0; i < oldColumn.length; i++)
                            floats[i] = oldColumn.get(i);
                          extracted.columns.addNewFloat(df.columns.byIndex(0).name).init((i) => floats[i]);
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
                    }),
                    ui.iconFA('trash', () => {
                      accordionExtractors.root.removeChild(accordionExtractors.root.childNodes[j - 1]);
                      extractorContainerList = extractorContainerList.filter(e => e !== extractorContainerList[j - 1]);
                      extractorTypesList = extractorTypesList.filter(e => e !== extractorTypesList[j - 1]);
                      extractorChartsList = extractorChartsList.filter(e => e !== extractorChartsList[j - 1]);
                      extracted.columns.remove(extracted.columns.byIndex(j).name);
                      j = j - 1;
                    })
                  ])
                ])
              ])
            ]),
            ui.block75([extractorChartsList[j - 1]])
          ]);
          extractorContainerList[j - 1].replaceChild(extractorInputsNew, extractorInputsOld);
          extractorInputsOld = extractorInputsNew;
        });
        let extractorInputsOld = ui.inputs([extractorTypesList[j]]);
        extractorContainerList[j] = ui.div(extractorInputsOld);
        accordionExtractors.addPane('Extractor ' + (j + 1), () => extractorContainerList[j], true);
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
          grok.shell.tables.forEach(function (table, index) {
            context.setVariable('table' + index, table);
          });
          call.context = context;
          indicatorInputsNew = ui.div([
            ui.block25([
              ui.inputs([
                indicatorTypesList[k - 1],
                await call.getEditor(),
                ui.buttonsInput([
                  ui.divH([
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
                    }),
                    ui.iconFA('trash', () => {
                      accordionIndicators.root.removeChild(accordionIndicators.root.childNodes[k - 1]);
                      indicatorContainerList = indicatorContainerList.filter(e => e !== indicatorContainerList[k - 1]);
                      indicatorTypesList = indicatorTypesList.filter(e => e !== indicatorTypesList[k - 1]);
                      indicatorChartsList = indicatorChartsList.filter(e => e !== indicatorChartsList[k - 1]);
                      k = k - 1;
                    })
                  ])
                ])
              ])
            ]),
            ui.block75([indicatorChartsList[k - 1]])
          ]);
          indicatorContainerList[k - 1].replaceChild(indicatorInputsNew, indicatorInputsOld);
          indicatorInputsOld = indicatorInputsNew;
        });
        let indicatorInputsOld = ui.inputs([indicatorTypesList[k]]);
        indicatorContainerList[k] = ui.div(indicatorInputsOld);
        accordionIndicators.addPane('Indicator ' + (k + 1), () => indicatorContainerList[k], true);
        k++;
      });

      grok.shell.newView('BioSignals', [
        appDescription,
        chosenDatabase,
        enterSamplingFrequencyDiv,
        chosenRecordDiv,
        samplingFrequencyDiv,
        annotationViewerDiv,
        ui.h2('Filtering and preprocessing'),
        accordionFilters,
        addFilterButton,
        ui.h2('Information extraction'),
        accordionExtractors,
        addExtractorButton,
        ui.h2('Physiological indicators'),
        accordionIndicators,
        addIndicatorButton
      ]);
    });
    chosenRecordDiv.innerHTML = '';
    chosenRecordDiv.append(chosenRecord.root);
  });

  grok.shell.newView('BioSignals', [
    appDescription,
    chosenDatabase,
    enterSamplingFrequencyDiv,
    chosenRecordDiv
  ]);
}