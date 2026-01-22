import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package, SUPPORTED_VIEWS, VIEW_CREATE_FUNC} from '../package';
import {CDISC_STANDARD, ClinCaseTableView, ClinStudyConfig} from './types';
import {defineXmlFileName, STENDTC, STSTDTC, StudyConfigFileName,
  studyConfigJsonFileName} from '../constants/constants';
import X2JS from 'x2js';
import {hideValidationColumns, readClinicalFile,
  removeExtension, studyConfigToMap} from './utils';
import {SUBJECT_ID, TSPARM, TSPARMCD, TSVAL} from '../constants/columns-constants';
import {ClinicalStudy} from '../clinical-study';
import {SUMMARY_VIEW_NAME, VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import {TABLE_VIEWS_META} from './views-creation-utils';
import {Subject, Subscription} from 'rxjs';
import {createInitialSatistics} from './initial-statistics-widget';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {ValidationHelper} from '../helpers/validation-helper';
import {getRequiredColumnsByView, handleMouseMoveOverErrorCell, setupValidationErrorColumns,
  setupValidationErrorIndicators} from './views-validation-utils';
import {DOMAINS_CATEGORIES_LIST, DOMAINS_DESCRIPTIONS, SUPP_DOMAIN_CATEGORY} from '../constants/domains-constants';
import {setupTableViewLayout} from './layout-utils';

export const validationNodes: {[key: string]: DG.TreeViewNode} = {};
export let currentOpenedView: DG.ViewBase | null = null;
export const CLINICAL_CASE_APP_PATH: string = '/apps/ClinicalCase/ClinicalCase';
export const PRECLINICAL_CASE_APP_PATH: string = '/apps/ClinicalCase/PreclinicalCase';
export const studyLoadedSubject = new Subject<{name: string, loaded: boolean, errorDomains?: string[]}>();
export const studies: {[key: string]: ClinicalStudy} = {};

export async function cdiscAppTB(treeNode: DG.TreeViewGroup, standard: string,
  currentStudy: string, currentViewName: string) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  //create import study view
  const importStudyItem = treeNode.item('Import study');
  importStudyItem.onSelected.subscribe(() => createImportStudyView(standard));
  //this creates studies objects and tree view nodes
  await createStudiesFromAppData(treeNode, standard as CDISC_STANDARD);
  //opens exact study and view
  openStudy(treeNode, standard, currentStudy, currentViewName);
  loaderItem.remove();
}

export function createImportStudyView(standard: string) {
  currentOpenedView = DG.View.create();
  currentOpenedView.root.classList.add('clinical-case-study-import-view');
  currentOpenedView.name = `Import Study - ${standard === CDISC_STANDARD.SDTM ? 'Clinical Case' : 'Preclinical Case'}`;

  let studyConfig: ClinStudyConfig | null = null;
  // Get the tree node for creating study config
  const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').
    getOrCreateGroup(standard === CDISC_STANDARD.SDTM ? 'Clinical Case' : 'Preclinical Case');

  const fileNamesDiv = ui.div();
  const errorDiv = ui.div();
  const filesSection = ui.divV([
    fileNamesDiv,
    errorDiv,
  ], {style: {marginTop: '10px'}});
  const statisticsDiv = ui.div('', {style: {marginLeft: '10px', marginTop: '10px', position: 'relative'}});
  const importStudyStatusDiv = ui.div('', {style: {position: 'relative', marginLeft: '20px'}});

  const filesInput = ui.input.files('', {
    onValueChanged: async () => {
      // Clear previous content
      ui.empty(fileNamesDiv);
      ui.empty(errorDiv);
      ui.empty(statisticsDiv);
      ui.empty(importStudyStatusDiv);

      // Show list of loaded file names as tags in horizontal flex container
      if (filesInput.value && filesInput.value.length > 0) {
        const fileTagsContainer = ui.div(
          filesInput.value.map((file) => {
            const tag = ui.div(file.name, {
              style: {
                display: 'inline-block',
                padding: '4px 8px',
                marginRight: '8px',
                marginBottom: '4px',
                backgroundColor: '#E7F0F3',
                color: 'var(--blue-6)',
                borderRadius: '4px',
                fontSize: '12px',
              },
            });
            tag.classList.add('clinical-case-import-study-file-tag');
            return tag;
          }),
          {style: {display: 'flex', flexWrap: 'wrap', marginTop: '10px', marginLeft: '20px'}},
        );
        fileNamesDiv.append(fileTagsContainer);

        try {
          ui.setUpdateIndicator(statisticsDiv, true, 'Collecting study summary...');
          importStudyButton.disabled = true;

          // Create study config from files
          studyConfig = await createStudyWithConfig(filesInput.value, clinicalCaseNode, true);

          // Show basic study statistics from configuration
          if (studyConfig) {
            const configMap = studyConfigToMap(studyConfig);
            const statsTable = ui.tableFromMap(configMap);
            statisticsDiv.append(ui.divText('Summary', {style: {marginLeft: '10px'}}));
            statisticsDiv.append(statsTable);
            ui.setUpdateIndicator(statisticsDiv, false);
            importStudyButton.disabled = false;
          }
        } catch (e: any) {
          const errorMessage = e?.message ?? String(e);
          errorDiv.append(ui.divText(`Error: ${errorMessage}`,
            {style: {color: 'red', marginLeft: '20px', marginTop: '10px'}}));
          ui.empty(statisticsDiv);
          importStudyButton.disabled = true;
          grok.shell.error(e);
        }
      }
    },
    acceptExtensions: ['.xpt', '.csv', '.xml', '.json'],
  });
  filesInput.root.style.paddingTop = '6px';

  const importStudyButton = ui.button('Import', async () => {
    importStudyStatusDiv.classList.add('clinical-case-study-import-in-progress-div');
    ui.setUpdateIndicator(importStudyStatusDiv, true, `Loading data for study ${studyConfig.name}`);
    const sub = studyLoadedSubject.subscribe((data) => {
      if (data.name === studyConfig.name) {
        sub.unsubscribe();
        if (data.loaded) {
          importStudyButton.disabled = true;
          importStudyStatusDiv.append(ui.divText(`Study ${studyConfig.name} loaded successfully`,
            {style: {color: 'green'}}));
          const tags = Array.from(fileNamesDiv.querySelectorAll('.clinical-case-import-study-file-tag'));
          if (data.errorDomains) {
            for (const errorDomain of data.errorDomains) {
              const tag = tags.filter((it) => it.textContent.toLocaleLowerCase() === errorDomain.toLowerCase());
              if (tag.length)
                (tag[0] as HTMLElement).style.color = 'red';
            }
          }
        } else {
          importStudyStatusDiv.append(ui.divText(`Error loading study ${studyConfig.name}`,
            {style: {color: 'red'}}));
        }
        importStudyStatusDiv.classList.remove('clinical-case-study-import-in-progress-div');
        ui.setUpdateIndicator(importStudyStatusDiv, false);
      }
    });
    //save all selected files into study directory since all .xpt files are needed to perform validation
    for (const file of filesInput.value) {
      await grok.dapi.files.write(
        `System:AppData/ClinicalCase/${studyConfig.standard}/${studyConfig.name}/${file.name}`,
        Array.from(file.data));
    }
    // eslint-disable-next-line max-len
    grok.dapi.files.writeAsText(`System:AppData/ClinicalCase/${studyConfig.standard}/${studyConfig.name}/${studyConfigJsonFileName}`,
      JSON.stringify(studyConfig));
    addStudyToBrowseTree(studies[studyConfig.name], clinicalCaseNode, filesInput.value);
    openStudy(clinicalCaseNode, studyConfig.standard, studyConfig.name, SUMMARY_VIEW_NAME);
  });

  importStudyButton.disabled = true;

  currentOpenedView.root.append(ui.divV([
    ui.divH([ui.h3('Import study files'),
      filesInput.root, importStudyButton], {style: {gap: '10px', marginLeft: '20px'}}),
    importStudyStatusDiv,
    filesSection,
    statisticsDiv,
  ]));

  grok.shell.addView(currentOpenedView);
}

export async function openApp(standard: CDISC_STANDARD): Promise<DG.ViewBase | void> {
  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + `/img/${standard === CDISC_STANDARD.SDTM ?
      'clin_case_icon.png' : 'preclinical_case_icon.png'}`,
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/ClinicalCase/README.md',
    description:
      `-  Visualize and explore your ${standard} data\n` +
      '-  Find patterns and trends in your data\n' +
      '-  Explore data on patient specific and trial specific levels\n' +
      '-  Browse through AEs and related data\n' +
      `-  Validate your ${standard} data`,
    iconSize: 120,
  });

  const studiesHeader = ui.h1('Studies');
  const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').
    getOrCreateGroup(standard === CDISC_STANDARD.SDTM ? 'Clinical Case' : 'Preclinical Case');
  clinicalCaseNode.currentItem = null;
  //this creates studies objects and tree view nodes
  await createStudiesFromAppData(clinicalCaseNode, standard);
  const studiesDiv = ui.div();
  studiesDiv.append(createInitialSatistics(clinicalCaseNode,
    Object.values(studies).filter((it) => it.config.standard === standard).map((it) => it.config)));
  const importStudyDiv = ui.div('', {style: {position: 'relative'}});
  const importStudyButton = ui.button('Import study...', () => {
    currentOpenedView?.close();
    const item = clinicalCaseNode.items.find((it) => it.text === 'Import study');
    clinicalCaseNode.currentItem = item;
  });
  importStudyButton.classList.add('clinical-case-import-study-button');
  importStudyDiv.append(importStudyButton);
  currentOpenedView?.close();
  currentOpenedView = DG.View.create();
  currentOpenedView.name = standard === CDISC_STANDARD.SDTM ? 'Clinical Case' : 'Preclinical Case';
  currentOpenedView.path = standard === CDISC_STANDARD.SEND ? PRECLINICAL_CASE_APP_PATH :
    CLINICAL_CASE_APP_PATH;
  currentOpenedView.root.append(ui.divV([
    appHeader,
    studiesHeader,
    studiesDiv,
    importStudyDiv,
  ]));
  return currentOpenedView;
}


export async function createStudiesFromAppData(treeNode: DG.TreeViewGroup, standard: CDISC_STANDARD) {
  const folders = await _package.files.list(`${standard}`);
  for (const folder of folders) {
    try {
      const filesList = await _package.files.list(folder);
      //const datasetsFolder = _package.files.list(`${standard}/studies/${folder.name}/datasets`);
      await createStudyWithConfig(filesList, treeNode);
    } catch (e: any) {
      grok.shell.error(`Error reading study config for study ${folder}: ${e?.message ?? e}`);
      continue;
    }
  }
}


async function createStudyWithConfig(files: DG.FileInfo[], treeNode: DG.TreeViewGroup,
  doNotAddToTree?: boolean): Promise<ClinStudyConfig> {
  let config: ClinStudyConfig | null = null;
  let dmDf: DG.DataFrame | null = null;
  let tsDf: DG.DataFrame | null = null;
  try {
    //if study has been opened previously, the config file will be saved within study folder
    const studySavedConfig = files.filter((it) => it.name === studyConfigJsonFileName);
    if (studySavedConfig.length)
      config = JSON.parse(await studySavedConfig[0].readAsString());

    //if study is loaded for the first time and no config has been saved previously
    //firsts extract all key fields from define.xml
    if (!config) {
      config = {standard: CDISC_STANDARD.SDTM};
      //look for define.xml
      const defineXml = files.filter((it) => it.name === defineXmlFileName);
      if (defineXml.length) {
        const parser = new X2JS();
        const defineJson = parser.xml2js(await defineXml[0].readAsString()) as any;
        config.name = defineJson?.ODM?.Study?.GlobalVariables?.StudyName;
        config.protocol = defineJson?.ODM?.Study?.GlobalVariables?.ProtocolName;
        config.description = defineJson?.ODM?.Study?.GlobalVariables?.StudyDescription;
        // eslint-disable-next-line max-len
        if (defineJson?.ODM?.Study?.MetaDataVersion?.['_def:StandardName'] && defineJson?.ODM?.Study.MetaDataVersion?.['_def:StandardName'].toLowerCase().includes('send'))
          config.standard = CDISC_STANDARD.SEND;
        //create fields descriptions
        if (defineJson?.ODM?.Study?.MetaDataVersion?.ItemDef &&
          Array.isArray(defineJson?.ODM?.Study?.MetaDataVersion?.ItemDef)) {
          config.fieldsDefinitions = {};
          for (const field of defineJson?.ODM?.Study?.MetaDataVersion?.ItemDef) {
            if (field._Name) {
              if (field.Description?.__text)
                config.fieldsDefinitions[field._Name] = field.Description.__text;
              else if (field.Description?.TranslatedText?.__text)
                config.fieldsDefinitions[field._Name] = field.Description?.TranslatedText?.__text;
            }
          }
        }
      } else {
      //if define.xml not found or there is no id in it - look for study.json
        const configFile = files.filter((it) => it.name === StudyConfigFileName);
        if (configFile.length) {
          const configJson = JSON.parse(await grok.dapi.files.readAsText(configFile[0]));
          config.name = configJson.name;
        }
      }
      if (!config.name)
        throw new Error(`Invalid or missing define.xml/study.json`);

      //second step - look for dm domain
      const dm = files.filter((it) => removeExtension(it.name) === 'dm');
      if (!dm.length)
        throw new Error(`No dm domain found for study ${config.name}`);
      dmDf = await readClinicalFile(dm[0]);
      if (dmDf) {
        dmDf.name = 'dm';
        config.totalSubjects = dmDf.col(SUBJECT_ID)?.categories?.length;
      }
      //now look for ts (trial summary) domain to extract other statistics
      const ts = files.filter((it) => removeExtension(it.name) === 'ts');
      if (ts.length) {
        try {
          tsDf = await readClinicalFile(ts[0]);
          tsDf.name = 'ts';
          const termCodeCol = tsDf.col(TSPARMCD);
          const termDescCol = tsDf.col(TSPARM);
          const valCol = tsDf.col(TSVAL);
          if (termCodeCol && valCol) {
            for (let i = 0; i < tsDf.rowCount; i++) {
              if (termCodeCol.get(i) === STSTDTC)
                config.startDate = valCol.get(i);
              else if (termCodeCol.get(i) === STENDTC)
                config.endDate = valCol.get(i);
              else {
                if (!config.other)
                  config.other = {};
                config.other[termDescCol ? termDescCol.get(i) : termCodeCol.get(i)] = valCol.get(i);
              }
            }
          }
        } catch (e: any) {
          grok.shell.warning(`ts domain could not be read: ${e?.message ?? e}`);
        }
      }
    }

    if (!studies[config.name]) {
      studies[config.name] = new ClinicalStudy(config);
      studies[config.name].domains.dm = dmDf;
      studies[config.name].domains.ts = tsDf;
      //write config file into folder
      if (!doNotAddToTree) {
        grok.dapi.files.writeAsText(
          `System:AppData/ClinicalCase/${config.standard}/${config.name}/${studyConfigJsonFileName}`,
          JSON.stringify(config));
        addStudyToBrowseTree(studies[config.name], treeNode, files);
      }
    }
    return config;
  } catch (e: any) {
    throw e;
  }
}

function addStudyToBrowseTree(study: ClinicalStudy, treeNode: DG.TreeViewGroup, studyFiles?: DG.FileInfo[]) {
  const node = treeNode.getOrCreateGroup(study.config.name, null, false);
  node.onSelected.subscribe(async (_) => {
    if (studies[study.studyId].changeViewToSummary)
      studies[study.studyId].currentViewName = SUMMARY_VIEW_NAME;
    else
      studies[study.studyId].changeViewToSummary = true;
    if (!node.expanded)
      node.expanded = true;
    else
      await openStudyNode(study, node, SUMMARY_VIEW_NAME);
  });

  node.onNodeExpanding.subscribe(async (_) => {
    if (studies[study.studyId].loadingStudyData === true)
      return;
    if (!studies[study.studyId].initCompleted) {
      const sub = studyLoadedSubject.subscribe((data) => {
        if (data.name === study.studyId) {
          sub.unsubscribe();
          if (data.loaded) {
            //adding views to tree only after data is loaded since we need all domains to perform validation
            for (const viewName of SUPPORTED_VIEWS[study.config.standard]) {
              //do not add view in case validation not passed
              const validator = new ValidationHelper(getRequiredColumnsByView(study.studyId)[viewName], study.studyId);
              if (!validator.validate())
                continue;
              const viewNode = node.item(viewName);
              if (viewName === VALIDATION_VIEW_NAME)
                validationNodes[study.studyId] = viewNode;
              viewNode.onSelected.subscribe(() => {
                if (studies[study.studyId].loadingStudyData === true) {
                  grok.shell.warning(`Loading data for study ${study.studyId}`);
                  treeNode.currentItem = node;
                  return;
                }
                studies[study.studyId].currentViewName = viewName;
                loadView(study, viewName, node);
                const browseTreeRoot = node.root.closest('.d4-tree-view-root');
                (browseTreeRoot as HTMLElement).focus();
              });
            }
            addDomainsToTree(study, node);
            openStudyNode(studies[study.studyId], node, studies[study.studyId].currentViewName);
          }
        }
      });
      initClinicalStudyData(studies[study.studyId], studyFiles);
    } else
      addDomainsToTree(study, node);
  });
}

function addDomainsToTree(study: ClinicalStudy, treeNode: DG.TreeViewGroup) {
  const domains = study.domains.all();
  const domainsNode = treeNode.group('Domains', null, false);
  domainsNode.onSelected.subscribe(() => {
    if (!domainsNode.expanded)
      domainsNode.expanded = true;
  });

  // Group domains by category
  const domainsByCategory: {[category: string]: typeof domains} = {};
  for (const domain of domains) {
    const domainInfo = DOMAINS_DESCRIPTIONS[domain.name];
    let category = domainInfo?.category;
    if (!category)
      category = domain.name.startsWith('supp') ? SUPP_DOMAIN_CATEGORY : 'Other';
    if (!domainsByCategory[category])
      domainsByCategory[category] = [];
    domainsByCategory[category].push(domain);
  }

  // Add domains grouped by category
  for (const category of DOMAINS_CATEGORIES_LIST) {
    if (domainsByCategory[category]) {
      const categoryNode = domainsNode.group(category, null, false);
      categoryNode.onSelected.subscribe(() => {
        if (!categoryNode.expanded)
          categoryNode.expanded = true;
      });
      for (const domain of domainsByCategory[category]) {
        const domainItem = categoryNode.item(domain.name);
        const domainInfo = DOMAINS_DESCRIPTIONS[domain.name];
        if (domainInfo?.description)
          ui.tooltip.bind(domainItem.root, domainInfo.description);
        domainItem.onSelected.subscribe(() => {
          addDomainAsTableView(study.studyId, domain);
          const browseTreeRoot = treeNode.root.closest('.d4-tree-view-root');
          (browseTreeRoot as HTMLElement).focus();
        });
      }
    }
  }
}

export function addDomainAsTableView(studyId: string, df: DG.DataFrame, closeCurrentPreview = true) {
  if (closeCurrentPreview)
    currentOpenedView?.close();
  currentOpenedView = grok.shell.addTablePreview(df);
  setupValidationErrorColumns(df);
  hideValidationColumns(currentOpenedView as DG.TableView);
  setupTableViewLayout(currentOpenedView as DG.TableView, studyId, df.name);

  let errorSubs: Subscription[] = [];
  const ribbons = currentOpenedView.getRibbonPanels();
  const showErrors = ui.input.bool('Show validation errors', {
    value: false,
    onValueChanged: () => {
      if (!showErrors.value) {
        errorSubs.forEach((sub) => sub.unsubscribe());
        (currentOpenedView as DG.TableView).grid.overlay.removeEventListener('mousemove', handleMouseMoveOverErrorCell);
      } else
        errorSubs = setupValidationErrorIndicators((currentOpenedView as DG.TableView).grid, df);
      (currentOpenedView as DG.TableView).grid.invalidate();
    },
  });
  ribbons.push([showErrors.root]);
  currentOpenedView.setRibbonPanels(ribbons);
}

async function openStudyNode(study: ClinicalStudy, node: DG.TreeViewGroup, currentViewName: string) {
  if (!currentViewName)
    currentViewName = SUMMARY_VIEW_NAME;
  if (studies[study.studyId].loadingStudyData) {
    if (studies[study.studyId].loadingStudyData === true) {
      grok.shell.warning(`Loading data for study ${study.studyId}`);
      return;
    }
  }
  loadView(study, currentViewName, node);
}

async function loadView(study: ClinicalStudy, viewName: string, parentNode: DG.TreeViewGroup) {
  let view = studies[study.studyId].views[viewName];
  let helper: any;
  if (!view) {
    if (Object.keys(TABLE_VIEWS_META).includes(viewName)) { //load table view
      const clinCaseTableView = VIEW_CREATE_FUNC[viewName](study.studyId, viewName) as ClinCaseTableView;
      view = clinCaseTableView.view;
      helper = clinCaseTableView.helper;
    } else { // load view
      if (viewName === SUMMARY_VIEW_NAME) {
        const errorLinkHandler = () => {
          parentNode.expanded = true;
          parentNode.currentItem = validationNodes[study.studyId];
        };
        view = VIEW_CREATE_FUNC[viewName](study.studyId, errorLinkHandler) as DG.ViewBase;
      } else
        view = VIEW_CREATE_FUNC[viewName](study.studyId) as DG.ViewBase;
    }
  }
  currentOpenedView?.close();
  currentOpenedView = grok.shell.addPreview(view);
  //return focus to browse tree to enable keyboard navigation
  const browseTreeRoot = parentNode.root.closest('.d4-tree-view-root');
  (browseTreeRoot as HTMLElement).focus();

  //in case of table view search for saved layout
  if (currentOpenedView.type === DG.TYPE.TABLE_VIEW) {
    const tableView = currentOpenedView as DG.TableView;
    await setupTableViewLayout(tableView, study.studyId, viewName);
  }
  if (view.hasOwnProperty('loaded') && !(view as ClinicalCaseViewBase).loaded)
    (view as ClinicalCaseViewBase).load();
  else
    helper?.propertyPanel();
  view.path =
        `${study.config.standard === CDISC_STANDARD.SEND ? PRECLINICAL_CASE_APP_PATH :
          CLINICAL_CASE_APP_PATH}/${encodeURI(study.studyId)}/${encodeURI(viewName)}`;
};

export function openStudy(treeNode: DG.TreeViewGroup, standard: string,
  currentStudyName: string, currentViewName: string) {
  if (currentStudyName && !Object.keys(studies).includes(currentStudyName))
    grok.shell.error(`Study ${currentStudyName} doesn't exist`);
  else if (currentStudyName) {
    const studyNode = treeNode.getOrCreateGroup(currentStudyName);
    console.log(Object.keys(VIEW_CREATE_FUNC));
    if (currentViewName && !Object.keys(VIEW_CREATE_FUNC).includes(currentViewName)) {
      grok.shell.warning(`${currentViewName} view doesn't exist, opening summary view`);
      currentViewName = SUMMARY_VIEW_NAME;
    } else if (!currentViewName)
      currentViewName = SUMMARY_VIEW_NAME;
    studies[currentStudyName].currentViewName = currentViewName;
    studies[currentStudyName].changeViewToSummary = false;
    //this will trigger onNodeExpanding
    if (treeNode.currentItem !== studyNode)
      treeNode.currentItem = studyNode;
    if (studyNode.expanded)
      openStudyNode(studies[currentStudyName], studyNode, currentViewName);
  }
}


export async function initClinicalStudyData(study: ClinicalStudy, studyFiles?: DG.FileInfo[]) {
  //closing the current view in case it was opened
  currentOpenedView?.close();
  const view = DG.View.create();
  currentOpenedView = grok.shell.addPreview(view);
  ui.setUpdateIndicator(view.root, true, `Loading data for study ${study.config.name}...`);
  if (!studies[study.studyId].initCompleted) {
    try {
      studies[study.studyId].loadingStudyData = true;
      const progressBar = DG.TaskBarProgressIndicator.create(`Reading data for study ${study.studyId}`);
      const dataLoaded = await readClinicalData(study, studyFiles);
      if (!dataLoaded.loaded) {
        studies[study.studyId].loadingStudyData = false;
        studyLoadedSubject.next({name: study.studyId, loaded: false});
      }
      studies[study.studyId].init();
      progressBar.close();
      grok.shell.info(`Data for study ${study.studyId} is ready`);
      studies[study.studyId].loadingStudyData = false;
      studyLoadedSubject.next({name: study.studyId, loaded: true, errorDomains: dataLoaded.errorDomains});
    } catch (e: any) {
      studies[study.studyId].loadingStudyData = false;
      studyLoadedSubject.next({name: study.studyId, loaded: false});
      throw e;
    }
  }
}


export async function readClinicalData(study: ClinicalStudy, importedFiles?: DG.FileInfo[]):
  Promise<{loaded: boolean, errorDomains: string[]}> {
  const notLoadedDomains: string[] = [];
  const pb = DG.TaskBarProgressIndicator
    .create(`Reading data for ${study.config.name}...`);
  try {
    const studyFiles = importedFiles ?? await _package.files.list(`${study.config.standard}/${study.studyId}`);
    const domainsNames = Object.keys(study.domains).filter((it) => it !== 'supp');

    //look for d42 file and read it in case it exists
    const d42DataFrames = studyFiles
      .filter((it) => removeExtension(it.fileName) === study.studyId && it.extension === 'd42');
    if (d42DataFrames.length) {
      const dfs = await grok.dapi.files.readBinaryDataFrames(d42DataFrames[0]);
      for (const df of dfs) {
        if (df.name.startsWith('supp') &&
              !studies[study.studyId].domains.supp.filter((it) => it.name == df.name).length)
          studies[study.studyId].domains.supp.push(df);
        else if (!studies[study.studyId].domains[df.name])
          studies[study.studyId].domains[df.name] = df;
      }
    } else { //if there is no .d42 file with dfs list, reading file by file, looking for xpt or csv
      for (let i = 0; i < studyFiles.length; i++) {
        const domainNameWithExt = studyFiles[i].fileName.toLowerCase();
        const domainNameWithoutExt = removeExtension(domainNameWithExt);
        pb.update(i / studyFiles.length * 100, `Reading ${domainNameWithExt}...`);
        if (!studies[study.studyId].domains[domainNameWithoutExt] &&
          (domainsNames.includes(domainNameWithoutExt) || domainNameWithoutExt.startsWith('supp'))) {
          const df = await readClinicalFile(studyFiles[i]);
          if (df) {
            df.name = domainNameWithoutExt;
            if (domainNameWithoutExt.startsWith('supp') &&
              !studies[study.studyId].domains.supp.filter((it) => it.name == domainNameWithoutExt).length)
              studies[study.studyId].domains.supp.push(df);
            else if (!studies[study.studyId].domains[domainNameWithoutExt])
              studies[study.studyId].domains[domainNameWithoutExt] = df;
          } else
            notLoadedDomains.push(domainNameWithExt);
        }
      }

      //saving .d42 format for further fast reading
      const d42FileDirectory =
        `System:AppData/ClinicalCase/${study.config.standard}/${study.studyId}/${study.studyId}.d42`;
      const dfsList = studies[study.studyId].domains.all();
      grok.dapi.files.writeBinaryDataFrames(d42FileDirectory, dfsList)
        .then(() => grok.shell.info(`.d42 file has been saved for study ${study.studyId}`));

      //in case we import study - also save config file
      if (importedFiles) {
        //save define.xml, if we have it. Otherwise, save study.json
        const defineXml = importedFiles.filter((it) => it.name === defineXmlFileName);
        if (defineXml.length) {
          grok.dapi.files.writeAsText(
            `System:AppData/ClinicalCase/${study.config.standard}/${study.studyId}/define.xml`,
            await defineXml[0].readAsString());
        } else {
          grok.dapi.files.writeAsText(
            `System:AppData/ClinicalCase/${study.config.standard}/${study.studyId}/study.json`,
            JSON.stringify(study.config));
        }
      }
    }

    // Ensure validation columns are added to domains
    // If validation is already completed, adds columns immediately
    // If validation is not completed yet, subscribes to add columns when validation completes
    study.ensureValidationColumnsAdded();

    return {loaded: true, errorDomains: notLoadedDomains};
  } finally {
    pb.close();
  }
}
