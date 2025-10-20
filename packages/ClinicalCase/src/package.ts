/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StudySummaryView} from './views/study-summary-view';
import {TimelinesView} from './views/timelines-view';
import {PatientProfileView} from './views/patient-profile-view';
import {AdverseEventsView} from './views/adverse-events-view';
import {ValidationView} from './views/validation-view';
import {AdverseEventHandler} from './panels/adverse-event-handler';
import {LaboratoryView} from './views/laboratory-view';
import {AERiskAssessmentView} from './views/ae-risk-assessment-view';
import {SurvivalAnalysisView} from './views/survival-analysis-view';
import {BoxPlotsView} from './views/boxplots-view';
import {MatrixesView} from './views/matrixes-view';
import {TimeProfileView} from './views/time-profile-view';
import {STUDY_ID} from './constants/columns-constants';
import {TreeMapView} from './views/tree-map-view';
import {MedicalHistoryView} from './views/medical-history-view';
import {VisitsView} from './views/visits-view';
import {StudyConfigurationView} from './views/study-config-view';
import {ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, COHORT_VIEW_NAME,
  CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME,
  PATIENT_PROFILE_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, SUMMARY_VIEW_NAME,
  SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME,
  VALIDATION_VIEW_NAME, VISITS_VIEW_NAME} from './constants/view-names-constants';
import {createClinCaseTableView, TABLE_VIEWS} from './utils/views-creation-utils';
import {CohortView} from './views/cohort-view';
import {QuestionnaiesView} from './views/questionnaires-view';
import {ClinCaseTableView, ClinStudyConfig} from './utils/types';
import {domainsToValidate, StudyJsonName} from './constants/constants';
import {ClinicalStudy, studies} from './clinical-study';
import {ClinicalCaseViewBase} from './model/ClinicalCaseViewBase';
import '../css/clinical-case.css';
import {u2} from '@datagrok-libraries/utils/src/u2';

export * from './package.g';

export const _package = new DG.Package();

export let validationRulesList = null;

export const VIEWS: {[key: string]: {[key: string]: DG.ViewBase}} = {};

const loadingStudyData: {[key: string]: boolean} = {};

const CLINICAL_CASE_APP_PATH: string = '/apps/ClinicalCase';

let currentOpenedView: DG.ViewBase | null = null;

type CurrentStudyAndView = {
  study: string;
  viewName: string;
}

export let existingStudies: {[key: string]: ClinStudyConfig} | null = null;
const validationNodes: {[key: string]: DG.TreeViewNode} = {};

const domains = (studyId: string, exactDomains?: string[]) =>
  (exactDomains ?? Object.keys(studies[studyId].domains)).map((it) => `${it.toLocaleLowerCase()}.csv`);
export let c: DG.FuncCall;

let cliniclaCaseLaunched = false;

function getCurrentStudyAndView(path: string): CurrentStudyAndView {
  let currentStudy = '';
  let currentViewName = '';
  if (path && !path.startsWith('/browse')) {
    const pathSegments = path.split('/');
    currentStudy = pathSegments[1];
    currentViewName = pathSegments.length > 2 ? pathSegments[2] : SUMMARY_VIEW_NAME;
  }
  return {study: currentStudy, viewName: currentViewName};
}

async function clinicalCaseAppTB(treeNode: DG.TreeViewGroup,
  currentStudy: string, currentViewName: string) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  await loadStudies([]);
  openStudy(treeNode, currentStudy, currentViewName);
  loaderItem.remove();
}

export function openStudy(treeNode: DG.TreeViewGroup,
  currentStudy: string, currentViewName: string) {
  let initialViewSelected = false;

  const openStudyNode = async (study: ClinStudyConfig, node: DG.TreeViewGroup) => {
    const summaryView = !VIEWS[study.name] ? null : VIEWS[study.name][SUMMARY_VIEW_NAME];
    if (!summaryView) {
      if (loadingStudyData[study.name]) {
        grok.shell.warning(`Loading data for study ${study.name}`);
        return;
      }
      if (!!studies[study.name]) {
        const dataRead = await readClinicalData(study, domainsToValidate);
        if (!dataRead)
          return;
        studies[study.name].validate();
        studies[study.name].processSitesAndSubjectCount();
      }
    }
    loadView(study, SUMMARY_VIEW_NAME, node);
  };


  const loadView = async (study: ClinStudyConfig, viewName: string, parentNode: DG.TreeViewGroup) => {
    let view = VIEWS[study.name][viewName];
    let helper: any;
    if (!view) {
      if (Object.keys(TABLE_VIEWS).includes(viewName)) { //load table view
        const clinCaseTableView = VIEW_CREATE_FUNC[viewName](study.name, viewName) as ClinCaseTableView;
        view = clinCaseTableView.view;
        helper = clinCaseTableView.helper;
      } else { // load view
        if (viewName === SUMMARY_VIEW_NAME) {
          const errorLinkHandler = () => {
            parentNode.expanded = true;
            parentNode.currentItem = validationNodes[study.name];
          };
          view = VIEW_CREATE_FUNC[viewName](study.name, errorLinkHandler) as DG.ViewBase;
        } else
          view = VIEW_CREATE_FUNC[viewName](study.name) as DG.ViewBase;
      }
    }
    currentOpenedView?.close();
    currentOpenedView = grok.shell.addPreview(view);
    if (view.hasOwnProperty('loaded') && !(view as ClinicalCaseViewBase).loaded)
      (view as ClinicalCaseViewBase).load();
    else
      helper?.propertyPanel();
    view.path =
        `browse${CLINICAL_CASE_APP_PATH}/${study.name}/${viewName.replaceAll(' ', '')}`;
  };

  for (const [_, study] of Object.entries(existingStudies)) {
    const node = treeNode.getOrCreateGroup(study.friendlyName ?? study.name, null, false);
    if (!studies[study.name]) {
      studies[study.name] = new ClinicalStudy(study.name);
      VIEWS[study.name] = {};
    }

    node.onSelected.subscribe(async (_) => {
      await openStudyNode(study, node);
    });

    node.onNodeExpanding.subscribe(async (_) => {
      if (loadingStudyData[study.name])
        return;
      for (const viewName of Object.keys(VIEW_CREATE_FUNC)) {
        const viewNode = node.item(viewName);
        if (viewName === VALIDATION_VIEW_NAME)
          validationNodes[study.name] = viewNode;
        viewNode.onSelected.subscribe(() => {
          if (loadingStudyData[study.name]) {
            grok.shell.warning(`Loading data for study ${study.name}`);
            treeNode.currentItem = node;
            return;
          }
          loadView(study, viewName, node);
        });
      }
      await initClinicalStudy(study);
      const viewItem = node.items
        .find((node) => node.text === (initialViewSelected ? currentViewName : SUMMARY_VIEW_NAME))?.root;
      viewItem?.click();
      initialViewSelected = false;
    });
  }

  if (currentStudy && !Object.keys(existingStudies).includes(currentStudy))
    grok.shell.error(`Study ${currentStudy} doesn't exist`);
  else if (currentStudy) {
    const studyNode = treeNode.getOrCreateGroup(currentStudy);
    if (currentViewName && !Object.keys(VIEW_CREATE_FUNC).includes(currentViewName)) {
      grok.shell.warning(`${currentViewName} view doesn't exist, opening summary view`);
      currentViewName = SUMMARY_VIEW_NAME;
    } else if (!currentViewName)
      currentViewName = SUMMARY_VIEW_NAME;
    initialViewSelected = true;
    if (studyNode.expanded)
      openStudyNode(existingStudies[currentStudy], studyNode);
    else
      studyNode.expanded = true;
  }
}


export async function initClinicalStudy(study: ClinStudyConfig) {
  if (!studies[study.name].initCompleted) {
    try {
      loadingStudyData[study.name] = true;
      const progressBar = DG.TaskBarProgressIndicator.create(`Reading data for study ${study.name}`);
      await readClinicalData(study);
      studies[study.name].init();
      progressBar.close();
      grok.shell.info(`Data for study ${study.name} is ready`);
    } finally {
      delete loadingStudyData[study.name];
    }
  }
}


export async function readClinicalData(study: ClinStudyConfig, domainsToDownLoad?: string[]): Promise<boolean> {
  const studyFiles = await _package.files.list(`studies/${study.name}`);
  const domainsList = domains(study.name, domainsToDownLoad);
  await Promise.all(studyFiles.map(async (file) => {
    const domainName = file.fileName.toLowerCase();
    if (!studies[study.name].domains[domainName] && domainsList.includes(domainName)) {
      const df = await _package.files.readCsv(`studies/${study.name}/${domainName}`);
      studies[study.name].domains[domainName.replace('.csv', '')] = df;
    }
  }));
  if (!studies[study.name].domains.dm) {
    grok.shell.error(`No demographic data found for study ${study.name}`);
    return false;
  }
  return true;
}


export async function loadStudies(deletedCampaigns: string[]): Promise<{[name: string]: ClinStudyConfig}> {
  if (!existingStudies) {
    const studiesFolders = (await _package.files.list(`studies`))
      .filter((f) => deletedCampaigns.indexOf(f.name) === -1);
    const studiesNamesMap: {[name: string]: ClinStudyConfig} = {};
    for (const folder of studiesFolders) {
      try {
        const studyJson: ClinStudyConfig = JSON.parse(await _package.files
          .readAsText(`studies/${folder.name}/${StudyJsonName}`));
        studiesNamesMap[studyJson.name] = studyJson;
      } catch (e) {
        continue;
      }
    }
    existingStudies = studiesNamesMap;
  }
  return existingStudies;
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
    const existingStudies = await loadStudies([]);
    const existingStudiesNames = Object.keys(existingStudies);
    const studiesDiv = ui.divV([]);
    for (const studyName of existingStudiesNames) {
      const studyLink = ui.link(studyName, async () => {
        const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
          .getOrCreateGroup('Clinical').getOrCreateGroup('Clinical Case');
        openStudy(clinicalCaseNode, studyName, SUMMARY_VIEW_NAME);
      }, 'Click to run the study');
      studyLink.style.paddingBottom = '10px';
      studiesDiv.append(studyLink);
    }
    const view = DG.View.create();
    view.name = 'Clinical Case';
    view.path = CLINICAL_CASE_APP_PATH;
    view.root.append(ui.divV([
      appHeader,
      studiesHeader,
      studiesDiv,
    ]));
    cliniclaCaseLaunched = true;
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
    await clinicalCaseAppTB(treeNode, studyAndView.study, studyAndView.viewName);
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
}

export const VIEW_CREATE_FUNC: {[key: string]: (studyId: string, args?: any) => DG.ViewBase | ClinCaseTableView} = {

  [SUMMARY_VIEW_NAME]:
    (studyId, errorLinkHandler?: () => void) => new StudySummaryView(SUMMARY_VIEW_NAME, studyId, errorLinkHandler),
  [TIMELINES_VIEW_NAME]: (studyId) => new TimelinesView(TIMELINES_VIEW_NAME, studyId),
  [LABORATORY_VIEW_NAME]: (studyId) => new LaboratoryView(LABORATORY_VIEW_NAME, studyId),
  [PATIENT_PROFILE_VIEW_NAME]: (studyId) => new PatientProfileView(PATIENT_PROFILE_VIEW_NAME, studyId),
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


