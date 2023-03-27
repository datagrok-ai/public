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
      DG.Column.fromList(DG.TYPE.FLOAT, 'time', []),
      DG.Column.fromList(DG.TYPE.FLOAT, 'values', [])
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
      'chosenRecord': chosenRecord
    });
    await call.call();
    //const annotations = call.getParamValue('annotations_df');
    const signals = call.getParamValue('signals_df');
    signals.columns.byIndex(0).setTag('samplingFrequency', call.getParamValue('sampling_frequency').toString())
    return signals;
  } catch (e) {
    grok.shell.error(e);
    throw e;
  }
}

async function getInitValues(isLocalTable, chosenDatabase, localTables, chosenRecord) {
  if (isLocalTable) {
    let signals = DG.DataFrame.fromColumns([
      localTables.find(({name}) => name === chosenDatabase).columns.byName(chosenRecord)
    ]);
    signals.name = chosenDatabase;
    return signals;
  } else {
    let pi = DG.TaskBarProgressIndicator.create('Loading record with annotations from Physionet...');
    let chosenDatabaseShortName = physionetDatabasesDictionary[chosenDatabase].shortName;
    let signals = await loadPhysionetRecordWithAnnotations(chosenDatabaseShortName, chosenRecord);
    signals.name = chosenDatabase + '/' + chosenRecord;
    pi.close();
    return signals;
  }
}

async function main(mainDiv, filterScripts, isLocalTable, chosenDatabase, localTables, chosenRecord, samplingFreq,
                    annotationViewerDiv, enterSamplingFrequencyDiv) {

  let signals = await getInitValues(isLocalTable, chosenDatabase, localTables, chosenRecord);
  signals.columns.byIndex(0).setTag('displayTitle', 'true');
  let context = DG.Context.create();
  context.setVariable(signals.name, signals);

  if (isLocalTable) {
    enterSamplingFrequencyDiv.innerHTML = '';
    signals.columns.byIndex(0).setTag('samplingFrequency', samplingFreq.value.toString());
  } else {
    parent.location.hash = chosenDatabase + '/' + chosenRecord;
  }
  annotationViewerDiv.innerHTML = '';
  annotationViewerDiv.append(
    ui.block([
      DG.Viewer.fromType('AnnotatorViewer', signals).root
    ])
  );

  // Filter dialogue
  let accordionFilters = ui.accordion();
  let filterTypesList = [];
  let filterContainerList = [];
  let filterChartsList = [];
  let i = 0;
  let addFilterButton = ui.button('Add Filter', async () => {
    filterChartsList[i] = ui.div();
    let tag = await grok.dapi.scripts.filter('#filters').list();
    filterTypesList[i] = ui.choiceInput('Filter ' + (i + 1), '', filterScripts.concat(tag), async () => {
      let call = filterTypesList[i - 1].value.prepare();
      for (let table of grok.shell.tables)
        if (table.name !== chosenDatabase)
          context.setVariable(table.name, table);
      call.context = context;
      filterContainerList[i - 1].replaceWith(ui.div([
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
                      df = DG.DataFrame.fromColumns([signals.columns.byIndex(signals.columns.length - 1)]);
                    } else {
                      signals = signals.append(df);
                      df = DG.DataFrame.fromColumns([df.columns.byIndex(df.columns.length - 1)]);
                    }
                    df.columns.byIndex(0).setTag('samplingFrequency', signals.columns.byIndex(0).getTag('samplingFrequency'));
                    filterChartsList[i - 1].replaceWith(DG.Viewer.fromType('AnnotatorViewer', df).root);
                  } catch (e) {
                    grok.shell.error(e);
                    throw e;
                  } finally {
                    pi.close();
                  }
                }),
                ui.button(
                  ui.iconFA('trash-alt', (ev) => {
                    let id = parseInt(ev.currentTarget.offsetParent.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.className.slice(41, 44));
                    let lst = filterTypesList.map((e) => e.caption).filter(function (el) {return el != null;});
                    let idx = lst.indexOf('Filter ' + id);
                    accordionFilters.root.removeChild(accordionFilters.root.childNodes[idx]);
                    filterContainerList.splice(idx, 1);
                    filterTypesList.splice(idx, 1);
                    filterChartsList.splice(idx, 1);
                    signals.columns.remove(signals.columns.byIndex(idx + 1).name);
                  })
                )
              ])
            ])
          ])
        ]),
        ui.block75([filterChartsList[i - 1]])
      ]));
    });
    filterContainerList[i] = ui.inputs([filterTypesList[i]]);
    accordionFilters.addPane('Filter ' + (i + 1), () => filterContainerList[i], true);
    i++;
  });

  // Information extraction dialogue
  let accordionExtractors = ui.accordion();
  let extracted;
  let extractorTypesList = [];
  let extractorChartsList = [];
  let extractorContainerList = [];
  let j = 0;
  let addExtractorButton = ui.button('Add Extractor', async () => {
    extractorChartsList[j] = getEmptyChart();
    let extractors = await grok.dapi.scripts.filter('#extractors').list();
    extractorTypesList[j] = ui.choiceInput('Extractor ' + (j + 1), '', extractors, async function () {
      let call = extractorTypesList[j - 1].value.prepare();
      let contextExtractors = DG.Context.create();
      contextExtractors.setVariable('table', signals);
      grok.shell.tables.forEach(function (table, index) {
        contextExtractors.setVariable('table' + index, table);
      });
      call.context = contextExtractors;
      extractorContainerList[j - 1].replaceWith(ui.div([
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
                ui.button(
                  ui.iconFA('trash-alt', (ev) => {
                    let id = parseInt(ev.currentTarget.offsetParent.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.className.slice(44, 47));
                    let lst = extractorTypesList.map((e) => e.caption).filter(function (el) {return el != null;});;
                    let idx = lst.indexOf('Extractor ' + id);
                    accordionExtractors.root.removeChild(accordionExtractors.root.childNodes[idx]);
                    extractorContainerList.splice(idx, 0);
                    extractorTypesList.splice(idx, 0);
                    extractorChartsList.splice(idx, 0);
                    extracted.columns.remove(extracted.columns.byIndex(idx).name);
                  })
                )
              ])
            ])
          ])
        ]),
        ui.block75([extractorChartsList[j - 1]])
      ]));
    });
    extractorContainerList[j] = ui.inputs([extractorTypesList[j]]);
    accordionExtractors.addPane('Extractor ' + (j + 1), () => extractorContainerList[j], true);
    j++;
  });

  // Indicators dialogue
  let accordionIndicators = ui.accordion();
  let indicatorTypesList = [];
  let indicatorChartsList = [];
  let indicatorContainerList = [];
  let k = 0;
  let addIndicatorButton = ui.button('Add Indicator', async () => {
    indicatorChartsList[k] = getEmptyChart();
    let indicatorsNames = await grok.dapi.scripts.filter('#indicators').list();
    indicatorTypesList[k] = ui.choiceInput('Indicator ' + (k + 1), '', indicatorsNames, async function () {
      let call = indicatorTypesList[k - 1].value.prepare();
      let contextIndicators = DG.Context.create();
      extracted.name = 'Extracted';
      contextIndicators.setVariable('table', extracted);
      grok.shell.tables.forEach(function (table, index) {
        contextIndicators.setVariable('table' + index, table);
      });
      call.context = contextIndicators;
      indicatorContainerList[k - 1].replaceWith(ui.div([
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
                ui.button(
                  ui.iconFA('trash-alt', (ev) => {
                    let id = parseInt(ev.currentTarget.offsetParent.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.className.slice(44, 47));
                    let lst = indicatorTypesList.map((e) => e.caption).filter(function (el) {return el != null;});;
                    let idx = lst.indexOf('Indicator ' + id);
                    accordionIndicators.root.removeChild(accordionIndicators.root.childNodes[idx]);
                    indicatorContainerList.splice(idx, 0);
                    indicatorTypesList.splice(idx, 0);
                    indicatorChartsList.splice(idx, 0);
                  })
                )
              ])
            ])
          ])
        ]),
        ui.block75([indicatorChartsList[k - 1]])
      ]));
    });
    indicatorContainerList[k] = ui.inputs([indicatorTypesList[k]]);
    accordionIndicators.addPane('Indicator ' + (k + 1), () => indicatorContainerList[k], true);
    k++;
  });

  mainDiv.replaceWith(
    ui.div([
      ui.h2('Filtering and preprocessing'),
      accordionFilters,
      addFilterButton,
      ui.h2('Information extraction'),
      accordionExtractors,
      addExtractorButton,
      ui.h2('Physiological indicators'),
      accordionIndicators,
      addIndicatorButton
    ])
  );
}

//name: Bio Signals
//tags: app
export function BioSignals() {

  let link = ui.element('a');
  link.href = 'https://github.com/datagrok-ai/public/tree/master/packages/BioSignals#readme'
  link.text = 'here';
  link.target = '_blank';

  let appDescription = ui.info(
    [
      ui.span(['See details ', link]),
      ui.divText('\n How to add your script:', {style: {'font-weight': 'bolder'}}),
      ui.divText('1. Go to Functions | Scripts | Actions | New <yourScriptLanguage> Script'),
      ui.divText('2. Write your script, and test it on files'),
      ui.divText('3. Set tag #filters, #extractors or #indicators'),
      ui.divText('Now it is available in corresponding app section')
    ], 'Analyze biomedical signals using built-in and custom scripts'
  );

  let windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const dspPackageFilters = ['DSP:SMA_filter', 'DSP:Exp_filter', 'DSP:Kalman_filter', 'DSP:MinMax_transform',
    'DSP:Zscore_transform', 'DSP:box_cox_transform', 'DSP:get_trend', 'DSP:remove_trend', 'DSP:fourier_filter',
    'DSP:spectral_density', 'DSP:subsample', 'DSP:asample'];
  let filterScripts = [];
  (async function () {
    for (let filter of dspPackageFilters) filterScripts.push(await grok.functions.eval(filter));
  })();

  let chosenRecordDiv = ui.div();
  let enterSamplingFrequencyDiv = ui.div();
  let annotationViewerDiv = ui.div();

  let localTables = grok.shell.tables;
  let namesOfLocalTables = localTables.map((df) => df.name);
  let tablesAndPhysionetDBs = namesOfLocalTables.concat(Object.keys(physionetDatabasesDictionary));

  let chosenDatabase = ui.choiceInput('Database', '', tablesAndPhysionetDBs, () => {
    let isLocalTable = (namesOfLocalTables.includes(chosenDatabase.stringValue));
    let items = (isLocalTable) ?
      localTables.find(({name}) => name === chosenDatabase.stringValue).columns.names() :
      physionetDatabasesDictionary[chosenDatabase.stringValue].namesOfRecords;
    let samplingFreq = ui.floatInput('Sampling frequency: ', '');
    if (isLocalTable) enterSamplingFrequencyDiv.append(samplingFreq.root);
    let columnName = (isLocalTable) ? 'Column' : 'Record';
    let chosenRecord = ui.choiceInput(columnName, '', items, async () => {
      main(mainDiv, filterScripts, isLocalTable, chosenDatabase.stringValue, localTables, chosenRecord.stringValue, samplingFreq, annotationViewerDiv, enterSamplingFrequencyDiv);
    });
    chosenRecordDiv.innerHTML = '';
    chosenRecordDiv.append(chosenRecord.root);
  });

  let mainDiv = ui.div();
  grok.shell.newView('BioSignals', [
    appDescription,
    ui.h2('Choose Physionet database or your file'),
    ui.div([
      ui.block25([
        ui.inputs([
          chosenDatabase,
          enterSamplingFrequencyDiv,
          chosenRecordDiv
        ])
      ]),
      ui.block75([annotationViewerDiv])
    ]),
    mainDiv
  ]);
  if (document.location.hash) {
    let s = document.location.hash.split('/');
    let db = s[0].replace(/%20/g, ' ').slice(1);
    let rec = (s.length === 3) ? s[1] + '/' + s[2] : s[1];
    chosenDatabase.value = db;
    chosenRecordDiv.innerHTML = '';
    let chosenRecord = ui.choiceInput('Record', '', physionetDatabasesDictionary[db].namesOfRecords, async () => {
      let isLocalTable = (namesOfLocalTables.includes(chosenDatabase.stringValue));
      let samplingFreq = ui.floatInput('Sampling frequency: ', '');
      if (isLocalTable) enterSamplingFrequencyDiv.append(samplingFreq.root);
      main(mainDiv, filterScripts, isLocalTable, chosenDatabase.stringValue, localTables, chosenRecord.stringValue, samplingFreq, annotationViewerDiv, enterSamplingFrequencyDiv);
    });
    chosenRecordDiv.append(chosenRecord.root);
    main(mainDiv, filterScripts, false, db, [''], rec, '', annotationViewerDiv, '');
  }
}
