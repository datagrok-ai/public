/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StudySummaryView} from './views/study-summary-view';
import {TimelinesView} from './views/timelines-view';
import {PatientProfileView} from './views/patient-profile-view';
import {AdverseEventsView} from './views/adverse-events-view';
import {ValidationView} from './views/validation-view';
import {LaboratoryView} from './views/laboratory-view';
import {AERiskAssessmentView} from './views/ae-risk-assessment-view';
import {SurvivalAnalysisView} from './views/survival-analysis-view';
import {BoxPlotsView} from './views/boxplots-view';
import {MatrixesView} from './views/matrixes-view';
import {TimeProfileView} from './views/time-profile-view';
import {STUDY_ID, SUBJECT_ID, TSPARM, TSPARMCD, TSVAL} from './constants/columns-constants';
import {TreeMapView} from './views/tree-map-view';
import {MedicalHistoryView} from './views/medical-history-view';
import {VisitsView} from './views/visits-view';
import {StudyConfigurationView} from './views/study-config-view';
import {ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME,
  ANIMAL_PROFILE_VIEW_NAME,
  CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME,
  PATIENT_PROFILE_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, SUMMARY_VIEW_NAME,
  SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME,
  VALIDATION_VIEW_NAME, VISITS_VIEW_NAME} from './constants/view-names-constants';
import {createClinCaseTableView, TABLE_VIEWS_META} from './utils/views-creation-utils';
//import {CohortView} from './views/cohort-view';
import {QuestionnaiesView} from './views/questionnaires-view';
import {CDISC_STANDARD, ClinCaseTableView, ClinStudyConfig} from './utils/types';
import {defineXmlFileName, STENDTC, STSTDTC, StudyConfigFileName, studyConfigJsonFileName} from './constants/constants';
import {ClinicalStudy} from './clinical-study';
import {ClinicalCaseViewBase} from './model/ClinicalCaseViewBase';
import '../css/clinical-case.css';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {scripts} from './package-api';
import {Subject} from 'rxjs';
import X2JS from 'x2js';
import {readClinicalFile, removeExtension} from './utils/utils';
import dayjs from 'dayjs';
import {createInitialSatistics} from './utils/initial-statistics-widget';

export * from './package.g';

export const _package = new DG.Package();

export let validationRulesList = null;

const CLINICAL_CASE_APP_PATH: string = '/apps/ClinicalCase';

let currentOpenedView: DG.ViewBase | null = null;

const studyLoadedSubject = new Subject<{name: string, loaded: boolean}>();

//export let study: ClinicalStudy = new ClinicalStudy();
export const studies: {[key: string]: ClinicalStudy} = {};

type CurrentStudyAndView = {
  study: string;
  standard: string;
  viewName: string;
}

const validationNodes: {[key: string]: DG.TreeViewNode} = {};

const domains = (studyId: string, exactDomains?: string[]) =>
  (exactDomains ?? Object.keys(studies[studyId].domains)).map((it) => it.toLocaleLowerCase());
export let c: DG.FuncCall;

function getCurrentStudyAndView(path: string): CurrentStudyAndView {
  let currentStandard = '';
  let currentStudy = '';
  let currentViewName = '';
  if (path && !path.startsWith('/browse')) {
    const pathSegments = path.split('/');
    currentStandard = pathSegments[1];
    currentStudy = pathSegments.length > 1 ? pathSegments[1] : '';
    currentViewName = pathSegments.length > 2 ? pathSegments[2] : SUMMARY_VIEW_NAME;
  }
  return {study: currentStudy, viewName: currentViewName, standard: currentStandard};
}

async function clinicalCaseAppTB(treeNode: DG.TreeViewGroup, standard: string,
  currentStudy: string, currentViewName: string) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  //this creates studies objects and tree view nodes
  await createStudiesFromAppData(treeNode);
  //opens exact study and view
  openStudy(treeNode, standard, currentStudy, currentViewName);
  loaderItem.remove();
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
};


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
  if (view.hasOwnProperty('loaded') && !(view as ClinicalCaseViewBase).loaded)
    (view as ClinicalCaseViewBase).load();
  else
    helper?.propertyPanel();
  view.path =
        `${CLINICAL_CASE_APP_PATH}/${study.config.standard}/${study.studyId}/${viewName.replaceAll(' ', '')}`;
};

function addStudyToBrowseTree(study: ClinicalStudy, treeNode: DG.TreeViewGroup, studyFiles?: DG.FileInfo[]) {
  const cdiscStandardNode = treeNode.getOrCreateGroup(study.config.standard, null, false);
  const node = cdiscStandardNode.getOrCreateGroup(study.config.name, null, false);

  node.onSelected.subscribe(async (_) => {
    studies[study.studyId].currentViewName = SUMMARY_VIEW_NAME;
    if (!node.expanded && !studies[study.studyId].initCompleted)
      node.expanded = true;
    else
      await openStudyNode(study, node, SUMMARY_VIEW_NAME);
  });

  node.onNodeExpanding.subscribe(async (_) => {
    if (studies[study.studyId].loadingStudyData === true)
      return;

    for (const viewName of SUPPORTED_VIEWS[study.config.standard]) {
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
      });
    }
    if (!studies[study.studyId].initCompleted) {
      const sub = studyLoadedSubject.subscribe((data) => {
        if (data.name === study.studyId) {
          sub.unsubscribe();
          if (data.loaded)
            openStudyNode(studies[study.studyId], node, studies[study.studyId].currentViewName);
        }
      });
      initClinicalStudyData(studies[study.studyId], studyFiles);
    }
  });
}


export function openStudy(treeNode: DG.TreeViewGroup, standard: string,
  currentStudyName: string, currentViewName: string) {
  if (currentStudyName && !Object.keys(studies).includes(currentStudyName))
    grok.shell.error(`Study ${currentStudyName} doesn't exist`);
  else if (currentStudyName) {
    const cdiscStandardNode = treeNode.getOrCreateGroup(standard);
    const studyNode = cdiscStandardNode.getOrCreateGroup(currentStudyName);
    if (currentViewName && !Object.keys(VIEW_CREATE_FUNC).includes(currentViewName)) {
      grok.shell.warning(`${currentViewName} view doesn't exist, opening summary view`);
      currentViewName = SUMMARY_VIEW_NAME;
    } else if (!currentViewName)
      currentViewName = SUMMARY_VIEW_NAME;
    studies[currentStudyName].currentViewName = currentViewName;
    //this will trigger onNodeExpanding
    if (!studyNode.expanded)
      studyNode.expanded = true;
    else
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
      if (!dataLoaded) {
        studies[study.studyId].loadingStudyData = false;
        studyLoadedSubject.next({name: study.studyId, loaded: false});
      }
      studies[study.studyId].init();
      progressBar.close();
      grok.shell.info(`Data for study ${study.studyId} is ready`);
      studies[study.studyId].loadingStudyData = false;
      studyLoadedSubject.next({name: study.studyId, loaded: true});
    } catch (e: any) {
      studies[study.studyId].loadingStudyData = false;
      studyLoadedSubject.next({name: study.studyId, loaded: false});
      throw e;
    }
  }
}


export async function readClinicalData(study: ClinicalStudy, importedFiles?: DG.FileInfo[]): Promise<boolean> {
  const pb = DG.TaskBarProgressIndicator
    .create(`Reading data for ${study.config.name}...`);
  try {
    const studyFiles = importedFiles ?? await _package.files.list(`studies/${study.studyId}`);
    const domainsList = domains(study.studyId);

    //look for d42 file and read it in case it exists
    const d42DataFrames = studyFiles
      .filter((it) => removeExtension(it.fileName) === study.studyId && it.extension === 'd42');
    if (d42DataFrames.length) {
      const dfs = await grok.dapi.files.readBinaryDataFrames(d42DataFrames[0]);
      for (const df of dfs) {
        if (!studies[study.studyId].domains[df.name])
          studies[study.studyId].domains[df.name] = df;
      }
    } else { //if there is no .d42 file with dfs list, reading file by file, looking for xpt or csv
      for (let i = 0; i < studyFiles.length; i++) {
        const domainNameWithExt = studyFiles[i].fileName.toLowerCase();
        const domainNameWithoutExt = removeExtension(domainNameWithExt);
        pb.update(i / studyFiles.length * 100, `Reading ${domainNameWithExt}...`);
        if (!studies[study.studyId].domains[domainNameWithoutExt] && domainsList.includes(domainNameWithoutExt)) {
          const df = await readClinicalFile(studyFiles[i]);
          if (df) {
            df.name = domainNameWithoutExt;
            if (!studies[study.studyId].domains[domainNameWithoutExt])
              studies[study.studyId].domains[domainNameWithoutExt] = df;
          }
        }
      }

      //saving .d42 format for further fast reading
      const d42FileDirectory = `System:AppData/ClinicalCase/studies/${study.studyId}/${study.studyId}.d42`;
      const dfsList = studies[study.studyId].domains.all();
      grok.dapi.files.writeBinaryDataFrames(d42FileDirectory, dfsList)
        .then(() => grok.shell.info(`.d42 file has been saved for study ${study.studyId}`));

      //in case we import study - also save config file
      if (importedFiles) {
        //save define.xml, if we have it. Otherwise, save study.json
        const defineXml = importedFiles.filter((it) => it.name === defineXmlFileName);
        if (defineXml.length) {
          grok.dapi.files.writeAsText(`System:AppData/ClinicalCase/studies/${study.studyId}/define.xml`,
            await defineXml[0].readAsString());
        } else {
          grok.dapi.files.writeAsText(`System:AppData/ClinicalCase/studies/${study.studyId}/study.json`,
            JSON.stringify(study.config));
        }
      }
    }
    return true;
  } finally {
    pb.close();
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
        grok.dapi.files.writeAsText(`System:AppData/ClinicalCase/studies/${config.name}/${studyConfigJsonFileName}`,
          JSON.stringify(config));
        addStudyToBrowseTree(studies[config.name], treeNode, files);
      }
    }
    return config;
  } catch (e: any) {
    throw e;
  }
}

export async function createStudiesFromAppData(treeNode: DG.TreeViewGroup) {
  const folders = await _package.files.list(`studies`);
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

export class PackageFunctions {
  @grok.decorators.app({
    'browsePath': 'Clinical',
    'name': 'Clinical Case',
  })
  static async clinicalCaseApp(): Promise<DG.ViewBase | void> {
    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/img/clin_case_icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/ClinicalCase/README.md',
      description:
      '-  Visualize and explore your SDTM data\n' +
      '-  Find patterns and trends in your data\n' +
      '-  Explore data on patient specific and trial specific levels\n' +
      '-  Browse through AEs and related data\n' +
      '-  Validate your SDTM data',
    });

    const studiesHeader = ui.h1('Studies');
    const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
      .getOrCreateGroup('Clinical').getOrCreateGroup('Clinical Case');
    //this creates studies objects and tree view nodes
    await createStudiesFromAppData(clinicalCaseNode);
    const studiesDiv = ui.div();
    studiesDiv.append(createInitialSatistics(clinicalCaseNode, Object.values(studies).map((it) => it.config)));

    const importStudyDiv = ui.div('', {style: {position: 'relative'}});
    const importStudyButton = ui.button('Import study...', () => {
      const dialog = ui.dialog('Import study');
      let studyConfig: ClinStudyConfig | null = null;
      const filesInput = ui.input.files('Files', {
        onValueChanged: async () => {
          //check if configuration file is present in the folder
          try {
            dialog.getButton('OK').disabled = true;
            studyConfig = await createStudyWithConfig(filesInput.value, clinicalCaseNode, true);
            dialog.getButton('OK').disabled = false;
          } catch (e) {
            grok.shell.error(e);
          }
          dialog.getButton('OK').disabled = !filesInput.validate();
        },
      });
      dialog.add(filesInput.root).show({resizable: true})
        .onOK(() => {
          setTimeout(() => {
            ui.setUpdateIndicator(importStudyDiv, true, `Loading data for study ${studyConfig.name}`);
            const sub = studyLoadedSubject.subscribe((data) => {
              if (data.name === studyConfig.name) {
                sub.unsubscribe();
                if (data.loaded) {
                  ui.empty(studiesDiv);
                  studiesDiv.append(createInitialSatistics(clinicalCaseNode,
                    Object.values(studies).map((it) => it.config)));
                }
                ui.setUpdateIndicator(importStudyDiv, false);
              }
            });
            // eslint-disable-next-line max-len
            grok.dapi.files.writeAsText(`System:AppData/ClinicalCase/studies/${studyConfig.name}/${studyConfigJsonFileName}`,
              JSON.stringify(studyConfig));
            addStudyToBrowseTree(studies[studyConfig.name], clinicalCaseNode, filesInput.value);
            openStudy(clinicalCaseNode, studyConfig.standard, studyConfig.name, SUMMARY_VIEW_NAME);
          }, 100);
        });
    });
    importStudyButton.classList.add('clinical-case-import-study-button');
    importStudyDiv.append(importStudyButton);
    const view = DG.View.create();
    view.name = 'Clinical Case';
    view.path = CLINICAL_CASE_APP_PATH;
    view.root.append(ui.divV([
      appHeader,
      studiesHeader,
      studiesDiv,
      importStudyDiv,
    ]));
    return view;
  }

  @grok.decorators.func()
  static async clinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    if (!validationRulesList)
      validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);
    const url = new URL(window.location.href);
    const currentStudyAndViewPath = url.pathname.includes(`/browse${CLINICAL_CASE_APP_PATH}`) ?
      url.pathname.replace(`/browse${CLINICAL_CASE_APP_PATH}`, ``) : '';
    const studyAndView = getCurrentStudyAndView(currentStudyAndViewPath);
    await clinicalCaseAppTB(treeNode, studyAndView.standard, studyAndView.study, studyAndView.viewName);
  }

  @grok.decorators.func({
    'name': 'Get list of studies',
    'description': 'Return list of clinical and preclinical studies loaded into Clinical Case application',
  })
  static async getListOfStudies(
    @grok.decorators.param({type: 'string', options: {optional: true}}) name?: string,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'More detailed study information including species, drug, dosing'}}) description?: string,
    @grok.decorators.param({type: 'int', options: {optional: true}}) numSubjects?: number,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'string', options: {optional: true, description: '>, <, ='}}) numSubjectsOperator?: string,
    @grok.decorators.param({type: 'datetime', options: {optional: true}}) startDate?: dayjs.Dayjs,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'string', options: {optional: true, description: '>, <, ='}}) startDateOperator?: string,
    @grok.decorators.param({type: 'datetime', options: {optional: true}}) endDate?: dayjs.Dayjs,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'string', options: {optional: true, description: '>, <, ='}}) endDateOperator?: string,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'bool', options: {optional: true}}) ongoing?: boolean,
    // eslint-disable-next-line max-len
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'CDISC data format, either SDTM or SEND'}})
      standard?: CDISC_STANDARD): Promise<DG.Widget> {
    const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
      .getOrCreateGroup('Clinical').getOrCreateGroup('Clinical Case');
    await createStudiesFromAppData(clinicalCaseNode);
    const filterConfig: ClinStudyConfig = {
      name: name ? name : undefined,
      description: description ? description : undefined,
      totalSubjects: numSubjects,
      startDate,
      endDate,
      standard: standard ? standard : undefined,
    };
    let filteredStudies: ClinicalStudy[] = [];
    for (const study of Object.values(studies)) {
      let dateFieldProcessed = false;
      let matched = true;
      for (const key of Object.keys(filterConfig)) {
        if (filterConfig[key]) {
          if ((key === 'startDate' || key === 'endDate') && !dateFieldProcessed) {
            dateFieldProcessed = true;
            const compareDates = (studyValue?: dayjs.Dayjs, filterValue?: dayjs.Dayjs | null, operator?: string) => {
              if (!studyValue || !filterValue)
                return false;
              const op = (operator ?? '>').trim();
              if (op === '<')
                return studyValue.isBefore(filterValue, 'day');
              if (op === '=')
                return studyValue.isSame(filterValue, 'day');
              return studyValue.isAfter(filterValue, 'day');
            };
            let studyMatch = false;
            if (startDate) {
              const startMatches = compareDates(study.config.startDate ?
                dayjs(study.config.startDate) : undefined, startDate, startDateOperator);
              if (endDate) {
                const endMathes = compareDates(study.config.endDate ?
                  dayjs(study.config.endDate) : undefined, endDate, endDateOperator);
                studyMatch = startMatches && endMathes;
              } else
                studyMatch = startMatches;
            } else if (endDate) {
              studyMatch = compareDates(study.config.endDate ?
                dayjs(study.config.endDate) : undefined, endDate, endDateOperator);
            }
            if (!studyMatch) {
              matched = false;
              break;
            }
          } else if (typeof(filterConfig[key]) === 'string') {
            if (!study.config[key]?.toLowerCase().includes(filterConfig[key]?.toLowerCase())) {
              matched = false;
              break;
            }
          } else if (typeof(filterConfig[key]) === 'number') {
            const op = (numSubjectsOperator ?? '=').trim();
            switch (op) {
            case '>':
              matched = study.config[key] && study.config[key] > filterConfig[key];
              break;
            case '<':
              matched = study.config[key] && study.config[key] < filterConfig[key];
              break;
            case '=':
              matched = study.config[key] && study.config[key] === filterConfig[key];
              break;
            }
            if (!matched)
              break;
          } else {
            if (!study.config[key] === filterConfig[key]) {
              matched = false;
              break;
            }
          }
        }
      }
      if (matched)
        filteredStudies.push(study);
    }
    if (ongoing)
      filteredStudies = filteredStudies.filter((it) => !it.config.endDate);
    const studiesDiv = createInitialSatistics(clinicalCaseNode, filteredStudies.map((it) => it.config));
    return new DG.Widget(studiesDiv);
  }

  @grok.decorators.folderViewer({
    'name': 'clinicalCaseFolderLauncher',
  })
  static async clinicalCaseFolderLauncher(
    folder: DG.FileInfo,
    files: DG.FileInfo[]):
    Promise<DG.Widget | undefined> {
    if (files.some((f) => f.fileName.toLowerCase() === 'dm.csv')) {
      const res = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
      const table = DG.DataFrame.fromCsv(res);
      const studyId = table.columns.names().includes(STUDY_ID) ? table.get(STUDY_ID, 0) : 'undefined';
      return DG.Widget.fromRoot(ui.div([
        ui.panel([
          ui.divText('Folder contains SDTM data'),
          ui.divText(`Study ID: ${studyId}`)]),
        ui.button('Run ClinicalCase', async () => {
          studies[studyId] = new ClinicalStudy(studyId);
          await Promise.all(files.map(async (file) => {
            if (domains(studyId).includes(file.fileName.toLowerCase())) {
              const df = await grok.data.files.openTable(`${folder.fullPath}/${file.fileName.toLowerCase()}`);
              grok.shell.addTableView(df);
            }
          }));
          grok.functions.call('Clinicalcase:clinicalCaseApp');
        }),
      ]));
    }
  }

  @grok.decorators.fileHandler({
    'ext': 'xpt',
  })
  static async xptFileHandler(
    @grok.decorators.param({'type': 'list'}) file: DG.FileInfo): Promise<DG.DataFrame[]> {
    const res: DG.DataFrame = await scripts.readSas(file);
    res.name = file.name;
    return [res];
  }
}

export const SUPPORTED_VIEWS: {[key: string]: string[]} = {
  [CDISC_STANDARD.SDTM]: [SUMMARY_VIEW_NAME, TIMELINES_VIEW_NAME, LABORATORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME,
    ADVERSE_EVENTS_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME,
    CORRELATIONS_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME,
    VISITS_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, VALIDATION_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME,
    AE_BROWSER_VIEW_NAME],

  [CDISC_STANDARD.SEND]: [SUMMARY_VIEW_NAME, TIMELINES_VIEW_NAME, LABORATORY_VIEW_NAME, ANIMAL_PROFILE_VIEW_NAME,
    DISTRIBUTIONS_VIEW_NAME, CORRELATIONS_VIEW_NAME, TIME_PROFILE_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME,
    VALIDATION_VIEW_NAME],
};

export const VIEW_CREATE_FUNC: {[key: string]: (studyId: string, args?: any) => DG.ViewBase | ClinCaseTableView} = {

  [SUMMARY_VIEW_NAME]:
    (studyId, errorLinkHandler?: () => void) => new StudySummaryView(SUMMARY_VIEW_NAME, studyId, errorLinkHandler),
  [TIMELINES_VIEW_NAME]: (studyId) => new TimelinesView(TIMELINES_VIEW_NAME, studyId),
  [LABORATORY_VIEW_NAME]: (studyId) => new LaboratoryView(LABORATORY_VIEW_NAME, studyId),
  [PATIENT_PROFILE_VIEW_NAME]: (studyId) => new PatientProfileView(PATIENT_PROFILE_VIEW_NAME, studyId),
  [ANIMAL_PROFILE_VIEW_NAME]: (studyId) => new PatientProfileView(ANIMAL_PROFILE_VIEW_NAME, studyId),
  [ADVERSE_EVENTS_VIEW_NAME]: (studyId) => new AdverseEventsView(ADVERSE_EVENTS_VIEW_NAME, studyId),
  [AE_RISK_ASSESSMENT_VIEW_NAME]: (studyId) => new AERiskAssessmentView(AE_RISK_ASSESSMENT_VIEW_NAME, studyId),
  [SURVIVAL_ANALYSIS_VIEW_NAME]: (studyId) => new SurvivalAnalysisView(SURVIVAL_ANALYSIS_VIEW_NAME, studyId),
  [DISTRIBUTIONS_VIEW_NAME]: (studyId) => new BoxPlotsView(DISTRIBUTIONS_VIEW_NAME, studyId),
  [CORRELATIONS_VIEW_NAME]: (studyId) => new MatrixesView(CORRELATIONS_VIEW_NAME, studyId),
  [TIME_PROFILE_VIEW_NAME]: (studyId) => new TimeProfileView(TIME_PROFILE_VIEW_NAME, studyId),
  [TREE_MAP_VIEW_NAME]: (studyId) => new TreeMapView(TREE_MAP_VIEW_NAME, studyId),
  [MEDICAL_HISTORY_VIEW_NAME]: (studyId) => new MedicalHistoryView(MEDICAL_HISTORY_VIEW_NAME, studyId),
  [VISITS_VIEW_NAME]: (studyId) => new VisitsView(VISITS_VIEW_NAME, studyId),
  [STUDY_CONFIGURATIN_VIEW_NAME]: (studyId, addView?: boolean) =>
    new StudyConfigurationView(STUDY_CONFIGURATIN_VIEW_NAME, studyId, addView),
  [VALIDATION_VIEW_NAME]: (studyId) => new ValidationView(VALIDATION_VIEW_NAME, studyId),
  [QUESTIONNAIRES_VIEW_NAME]: (studyId) => new QuestionnaiesView(QUESTIONNAIRES_VIEW_NAME, studyId),
  [AE_BROWSER_VIEW_NAME]: (studyId, viewName) => createClinCaseTableView(studyId, viewName),
};


