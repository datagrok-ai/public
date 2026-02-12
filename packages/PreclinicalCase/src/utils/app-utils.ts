import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package, SUPPORTED_VIEWS, VIEW_CREATE_FUNC} from '../package';
import {StudyConfig, TableView} from '../types/types';
import {defineXmlFileName, STENDTC, STSTDTC, StudyConfigFileName,
  studyConfigJsonFileName} from '../constants/constants';
import X2JS from 'x2js';
import {hideValidationColumns, readClinicalFile,
  removeExtension, studyConfigToMap} from './utils';
import {SUBJECT_ID, TSPARM, TSPARMCD, TSVAL} from '../constants/columns-constants';
import {PreclinicalStudy} from '../preclinical-study';
import {SUMMARY_VIEW_NAME, VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import {Subject, Subscription} from 'rxjs';
import {createInitialSatistics} from './initial-statistics-widget';
import {handleMouseMoveOverErrorCell, setupValidationErrorColumns,
  setupValidationErrorIndicators} from './views-validation-utils';
import {DOMAINS_CATEGORIES_LIST, DOMAINS_DESCRIPTIONS, SUPP_DOMAIN_CATEGORY} from '../constants/domains-constants';
import {setupTableViewLayout} from './layout-utils';
import { createImportStudyView } from '../views/import-study-view';

export const validationNodes: {[key: string]: DG.TreeViewNode} = {};
export let currentOpenedView: DG.ViewBase | null = null;
export const PRECLINICAL_CASE_APP_PATH: string = '/apps/Preclinicalcase';
export const studyLoadedSubject = new Subject<{name: string, loaded: boolean, errorDomains?: string[]}>();
export const studies: {[key: string]: PreclinicalStudy} = {};

export async function cdiscAppTB(treeNode: DG.TreeViewGroup,
  currentStudy: string, currentViewName: string) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  const importStudyItem = treeNode.item('Import study');
  importStudyItem.onSelected.subscribe(async () => {
    currentOpenedView?.close();
    currentOpenedView = await createImportStudyView();
    grok.shell.addPreview(currentOpenedView);
  });
  await createStudiesFromAppData(treeNode);
  openStudy(treeNode, currentStudy, currentViewName);
  loaderItem.remove();
}


export async function openApp(): Promise<DG.ViewBase | void> {
  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/img/preclinical_case_icon.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/PreclinicalCase/README.md',
    description:
      '-  Visualize and explore your SEND data\n' +
      '-  Find patterns and trends in your data\n' +
      '-  Explore data on subject specific and study specific levels\n' +
      '-  Browse through findings and related data\n' +
      '-  Validate your SEND data',
    iconSize: 120,
  });

  const studiesHeader = ui.h1('Studies');
  const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').
    getOrCreateGroup('Preclinical Case');
  await createStudiesFromAppData(clinicalCaseNode);
  const studiesDiv = ui.div();
  studiesDiv.append(createInitialSatistics(clinicalCaseNode,
    Object.values(studies).map((it) => it.config)));
  const importStudyDiv = ui.div('', {style: {position: 'relative'}});
  const importStudyButton = ui.button('Import study...', () => {
    currentOpenedView?.close();
    const item = clinicalCaseNode.items.find((it) => it.text === 'Import study');
    clinicalCaseNode.currentItem = item!;
  });
  importStudyButton.classList.add('preclinical-case-import-study-button');
  importStudyDiv.append(importStudyButton);
  currentOpenedView?.close();
  currentOpenedView = DG.View.create();
  currentOpenedView.name = 'Preclinical Case';
  currentOpenedView.path = PRECLINICAL_CASE_APP_PATH;
  currentOpenedView.root.append(ui.divV([
    appHeader,
    studiesHeader,
    studiesDiv,
    importStudyDiv,
  ]));
  return currentOpenedView;
}


export async function createStudiesFromAppData(treeNode: DG.TreeViewGroup) {
  const folders = await _package.files.list('SEND');
  for (const folder of folders) {
    try {
      const filesList = await _package.files.list(folder);
      await createStudyWithConfig(filesList, treeNode);
    } catch (e: any) {
      grok.shell.error(`Error reading study config for study ${folder}: ${e?.message ?? e}`);
      continue;
    }
  }
}


export async function createStudyWithConfig(files: DG.FileInfo[], treeNode: DG.TreeViewGroup,
  doNotAddToTree?: boolean): Promise<StudyConfig> {
  let config: StudyConfig | null = null;
  let dmDf: DG.DataFrame | null = null;
  let tsDf: DG.DataFrame | null = null;
  try {
    const studySavedConfig = files.filter((it) => it.name === studyConfigJsonFileName);
    if (studySavedConfig.length)
      config = JSON.parse(await studySavedConfig[0].readAsString());

    if (!config) {
      config = {} as StudyConfig;
      const defineXml = files.filter((it) => it.name === defineXmlFileName);
      if (defineXml.length) {
        const parser = new X2JS();
        const defineJson = parser.xml2js(await defineXml[0].readAsString()) as any;
        config.name = defineJson?.ODM?.Study?.GlobalVariables?.StudyName;
        config.protocol = defineJson?.ODM?.Study?.GlobalVariables?.ProtocolName;
        config.description = defineJson?.ODM?.Study?.GlobalVariables?.StudyDescription;
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
        const configFile = files.filter((it) => it.name === StudyConfigFileName);
        if (configFile.length) {
          const configJson = JSON.parse(await grok.dapi.files.readAsText(configFile[0]));
          config.name = configJson.name;
        }
      }
      if (!config.name)
        throw new Error('Invalid or missing define.xml/study.json');

      const dm = files.filter((it) => removeExtension(it.name) === 'dm');
      if (!dm.length)
        throw new Error(`No dm domain found for study ${config.name}`);
      dmDf = await readClinicalFile(dm[0]);
      if (dmDf) {
        dmDf.name = 'dm';
        config.totalSubjects = dmDf.col(SUBJECT_ID)?.categories?.length;
      }
      const ts = files.filter((it) => removeExtension(it.name) === 'ts');
      if (ts.length) {
        try {
          tsDf = await readClinicalFile(ts[0]);
          if (tsDf) {
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

          }
        } catch (e: any) {
          grok.shell.warning(`ts domain could not be read: ${e?.message ?? e}`);
        }
      }
    }

    if (!studies[config.name!]) {
      studies[config.name!] = new PreclinicalStudy(config);
      studies[config.name!].domains.dm = dmDf;
      studies[config.name!].domains.ts = tsDf;
      if (!doNotAddToTree) {
        grok.dapi.files.writeAsText(
          `System:AppData/Preclinicalcase/SEND/${config.name}/${studyConfigJsonFileName}`,
          JSON.stringify(config));
        addStudyToBrowseTree(studies[config.name!], treeNode, files);
      }
    }
    return config;
  } catch (e: any) {
    throw e;
  }
}

export function addStudyToBrowseTree(study: PreclinicalStudy, treeNode: DG.TreeViewGroup, studyFiles?: DG.FileInfo[]) {
  const node = treeNode.getOrCreateGroup(study.config.name!, null, false);
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
            for (const viewName of SUPPORTED_VIEWS) {
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
      initStudyData(studies[study.studyId], studyFiles);
    } else
      addDomainsToTree(study, node);
  });
}

function addDomainsToTree(study: PreclinicalStudy, treeNode: DG.TreeViewGroup) {
  const domains = study.domains.all();
  const domainsNode = treeNode.group('Domains', null, false);
  domainsNode.onSelected.subscribe(() => {
    if (!domainsNode.expanded)
      domainsNode.expanded = true;
  });

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

async function openStudyNode(study: PreclinicalStudy, node: DG.TreeViewGroup, currentViewName: string) {
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

async function loadView(study: PreclinicalStudy, viewName: string, parentNode: DG.TreeViewGroup) {
  currentOpenedView?.close();
  let view = studies[study.studyId].views[viewName];
  let helper: any;
  if (!view) {
    view = VIEW_CREATE_FUNC[viewName](study.studyId) as DG.ViewBase;
  }
  currentOpenedView = grok.shell.addPreview(view);
  const browseTreeRoot = parentNode.root.closest('.d4-tree-view-root');
  (browseTreeRoot as HTMLElement).focus();

  if (currentOpenedView.type === DG.TYPE.TABLE_VIEW) {
    const tableView = currentOpenedView as DG.TableView;
    await setupTableViewLayout(tableView, study.studyId, viewName);
  }
  if ((view as any).loaded === false && typeof (view as any).load === 'function')
    (view as any).load();
  else
    helper?.propertyPanel();
  view.path = `${PRECLINICAL_CASE_APP_PATH}/${encodeURI(study.studyId)}/${encodeURI(viewName)}`;
};

export function openStudy(treeNode: DG.TreeViewGroup,
  currentStudyName: string, currentViewName: string) {
  if (currentStudyName && !Object.keys(studies).includes(currentStudyName))
    grok.shell.error(`Study ${currentStudyName} doesn't exist`);
  else if (currentStudyName) {
    const studyNode = treeNode.getOrCreateGroup(currentStudyName);
    if (currentViewName && !Object.keys(VIEW_CREATE_FUNC).includes(currentViewName)) {
      grok.shell.warning(`${currentViewName} view doesn't exist, opening summary view`);
      currentViewName = SUMMARY_VIEW_NAME;
    } else if (!currentViewName)
      currentViewName = SUMMARY_VIEW_NAME;
    studies[currentStudyName].currentViewName = currentViewName;
    studies[currentStudyName].changeViewToSummary = false;
    if (treeNode.currentItem !== studyNode)
      treeNode.currentItem = studyNode;
    if (studyNode.expanded)
      openStudyNode(studies[currentStudyName], studyNode, currentViewName);
  }
}


export async function initStudyData(study: PreclinicalStudy, studyFiles?: DG.FileInfo[]) {
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


export async function readClinicalData(study: PreclinicalStudy, importedFiles?: DG.FileInfo[]):
  Promise<{loaded: boolean, errorDomains: string[]}> {
  const notLoadedDomains: string[] = [];
  const pb = DG.TaskBarProgressIndicator
    .create(`Reading data for ${study.config.name}...`);
  try {
    const studyFiles = importedFiles ?? await _package.files.list(`SEND/${study.studyId}`);
    const domainsNames = Object.keys(study.domains).filter((it) => it !== 'supp');

    const d42DataFrames = studyFiles
      .filter((it) => removeExtension(it.fileName) === study.studyId && it.extension === 'd42');
    if (d42DataFrames.length) {
      const dfs = await grok.dapi.files.readBinaryDataFrames(d42DataFrames[0]);
      for (const df of dfs) {
        if (df.name.startsWith('supp') &&
              !studies[study.studyId].domains.supp.filter((it) => it.name == df.name).length)
          studies[study.studyId].domains.supp.push(df);
        else if (!(studies[study.studyId].domains as any)[df.name])
          (studies[study.studyId].domains as any)[df.name] = df;
      }
    } else {
      for (let i = 0; i < studyFiles.length; i++) {
        const domainNameWithExt = studyFiles[i].fileName.toLowerCase();
        const domainNameWithoutExt = removeExtension(domainNameWithExt);
        pb.update(i / studyFiles.length * 100, `Reading ${domainNameWithExt}...`);
        if (!(studies[study.studyId].domains as any)[domainNameWithoutExt] &&
          (domainsNames.includes(domainNameWithoutExt) || domainNameWithoutExt.startsWith('supp'))) {
          const df = await readClinicalFile(studyFiles[i]);
          if (df) {
            df.name = domainNameWithoutExt;
            if (domainNameWithoutExt.startsWith('supp') &&
              !studies[study.studyId].domains.supp.filter((it) => it.name == domainNameWithoutExt).length)
              studies[study.studyId].domains.supp.push(df);
            else if (!(studies[study.studyId].domains as any)[domainNameWithoutExt])
              (studies[study.studyId].domains as any)[domainNameWithoutExt] = df;
          } else
            notLoadedDomains.push(domainNameWithExt);
        }
      }

      const d42FileDirectory =
        `System:AppData/Preclinicalcase/SEND/${study.studyId}/${study.studyId}.d42`;
      const dfsList = studies[study.studyId].domains.all();
      grok.dapi.files.writeBinaryDataFrames(d42FileDirectory, dfsList)
        .then(() => grok.shell.info(`.d42 file has been saved for study ${study.studyId}`));

      if (importedFiles) {
        const defineXml = importedFiles.filter((it) => it.name === defineXmlFileName);
        if (defineXml.length) {
          grok.dapi.files.writeAsText(
            `System:AppData/Preclinicalcase/SEND/${study.studyId}/define.xml`,
            await defineXml[0].readAsString());
        } else {
          grok.dapi.files.writeAsText(
            `System:AppData/Preclinicalcase/SEND/${study.studyId}/study.json`,
            JSON.stringify(study.config));
        }
      }
    }

    study.ensureValidationColumnsAdded();

    return {loaded: true, errorDomains: notLoadedDomains};
  } finally {
    pb.close();
  }
}
