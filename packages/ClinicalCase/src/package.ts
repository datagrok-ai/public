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
import {addView, createTableView, TABLE_VIEWS} from './utils/views-creation-utils';
import {CohortView} from './views/cohort-view';
import {QuestionnaiesView} from './views/questionnaires-view';
import {ClinCaseTableView, ClinStudyConfig} from './utils/types';
import {demoStudyId, domainsToValidate, StudyJsonName} from './constants/constants';
import {ClinicalStudy, studies} from './clinical-study';
import {ClinicalCaseViewBase} from './model/ClinicalCaseViewBase';

export const _package = new DG.Package();

export let validationRulesList = null;

export const VIEWS: {[key: string]: {[key: string]: DG.ViewBase}} = {};


const domains = (studyId: string, exactDomains?: string[]) =>
  (exactDomains ?? Object.keys(studies[studyId].domains)).map((it) => `${it.toLocaleLowerCase()}.csv`);
export let c: DG.FuncCall;

//tags: app
//name: Clinical Case
//output: view v
//meta.browsePath: Clinical
export async function clinicalCaseApp(): Promise<DG.ViewBase> {
  return DG.View.create();
}

//input: dynamic treeNode
//input: view browseView
export async function clinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup,
  browseView: DG.BrowseView) {// TODO: DG.BrowseView
  if (!validationRulesList)
    validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);
  await clinicalCaseAppTB(treeNode, browseView);
}

async function clinicalCaseAppTB(treeNode: DG.TreeViewGroup, browseView: DG.BrowseView) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  const existingStudies = await loadStudies([]);

  for (const [_, study] of Object.entries(existingStudies)) {
    const node = treeNode.group(study.friendlyName ?? study.name, null, false);
    if (!studies[study.name]) {
      studies[study.name] = new ClinicalStudy(study.name);
      VIEWS[study.name] = {};
    }

    node.onSelected.subscribe(async (_) => {
      const summaryView = !VIEWS[study.name] ? null : VIEWS[study.name][SUMMARY_VIEW_NAME];
      if (!summaryView) {
        const dataRead = await readClinicalData(study, domainsToValidate);
        if (!dataRead)
          return;
        //validation results are required for summary view
        studies[study.name].validate();
        studies[study.name].processSitesAndSubjectCount();
      }
      loadView(SUMMARY_VIEW_NAME, node);
    });

    const loadView = async (viewName: string, parentNode: DG.TreeViewGroup) => {
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
              parentNode.currentItem = validationNode;
            };
            view = VIEW_CREATE_FUNC[viewName](study.name, errorLinkHandler) as DG.ViewBase;
          } else
            view = VIEW_CREATE_FUNC[viewName](study.name) as DG.ViewBase;
        }
      }
      browseView.preview = view as DG.View;
      if (view.hasOwnProperty('loaded') && !(view as ClinicalCaseViewBase).loaded)
        (view as ClinicalCaseViewBase).load();
      else
        helper?.propertyPanel();
    };

    let validationNode: DG.TreeViewNode = null;
    node.onNodeExpanding.subscribe(async (_) => {
      for (const viewName of Object.keys(VIEW_CREATE_FUNC)) {
        const viewNode = node.item(viewName);
        if (viewName === VALIDATION_VIEW_NAME)
          validationNode = viewNode;
        viewNode.onSelected.subscribe(() => {
          loadView(viewName, node);
        });
      }
      await initClinicalStudy(study);
    });
  }
  loaderItem.remove();
}

export function createClinCaseTableView(studyId: string, viewName: string): ClinCaseTableView {
  const tableView = createTableView(
    studyId,
    TABLE_VIEWS[viewName].domainsAndColsToCheck,
    viewName,
    TABLE_VIEWS[viewName].helpUrl,
    TABLE_VIEWS[viewName].createViewHelper,
    TABLE_VIEWS[viewName].paramsForHelper,
  );
  return {view: tableView.view, helper: tableView.helper};
}

export async function initClinicalStudy(study: ClinStudyConfig) {
  if (!studies[study.name].initCompleted) {
    const progressBar = DG.TaskBarProgressIndicator.create(`Reading data for study ${study.name}`);
    await readClinicalData(study);
    studies[study.name].init();
    progressBar.close();
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
  return studiesNamesMap;
}

//input: string studyId {optional: true}
export async function createClinicalCaseViews(studyId?: string): Promise<any> {
  c = grok.functions.getCurrentCall();
  if (!validationRulesList)
    validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);

  if (!studyId) {
    const demoFiles = await grok.dapi.projects.filter('clin-demo-files-2').list();
    if (demoFiles.length) {
      studyId = demoStudyId;
      studies[studyId] = new ClinicalStudy(studyId);
      await (await grok.dapi.projects.find(demoFiles[0].id)).open();
      studies[studyId].initFromWorkspace();
    } else
      grok.shell.warning('Please load SDTM data or demo files');
  } else
    studies[studyId].init();

  if (!VIEWS[studyId])
    VIEWS[studyId] = {};

  VIEWS[studyId][SUMMARY_VIEW_NAME] = <StudySummaryView>addView(new StudySummaryView(SUMMARY_VIEW_NAME, studyId));
  VIEWS[studyId][VISITS_VIEW_NAME] = <VisitsView>addView(new VisitsView(VISITS_VIEW_NAME, studyId));
  //VIEWS[studyId][COHORT_VIEW_NAME] = <CohortView>addView(new CohortView(COHORT_VIEW_NAME, studyId));
  VIEWS[studyId][TIMELINES_VIEW_NAME] = <TimelinesView>addView(new TimelinesView(TIMELINES_VIEW_NAME, studyId));
  VIEWS[studyId][PATIENT_PROFILE_VIEW_NAME] =
  <PatientProfileView>addView(new PatientProfileView(PATIENT_PROFILE_VIEW_NAME, studyId));
  VIEWS[studyId][ADVERSE_EVENTS_VIEW_NAME] =
  <AdverseEventsView>addView(new AdverseEventsView(ADVERSE_EVENTS_VIEW_NAME, studyId));
  VIEWS[studyId][LABORATORY_VIEW_NAME] = <LaboratoryView>addView(new LaboratoryView(LABORATORY_VIEW_NAME, studyId));
  VIEWS[studyId][AE_RISK_ASSESSMENT_VIEW_NAME] =
  <AERiskAssessmentView>addView(new AERiskAssessmentView(AE_RISK_ASSESSMENT_VIEW_NAME, studyId));
  VIEWS[studyId][SURVIVAL_ANALYSIS_VIEW_NAME] =
  <SurvivalAnalysisView>addView(new SurvivalAnalysisView(SURVIVAL_ANALYSIS_VIEW_NAME, studyId));
  VIEWS[studyId][DISTRIBUTIONS_VIEW_NAME] = <BoxPlotsView>addView(new BoxPlotsView(DISTRIBUTIONS_VIEW_NAME, studyId));
  VIEWS[studyId][CORRELATIONS_VIEW_NAME] = <MatrixesView>addView(new MatrixesView(CORRELATIONS_VIEW_NAME, studyId));
  VIEWS[studyId][TIME_PROFILE_VIEW_NAME] =
  <TimeProfileView>addView(new TimeProfileView(TIME_PROFILE_VIEW_NAME, studyId));
  VIEWS[studyId][TREE_MAP_VIEW_NAME] = <TreeMapView>addView(new TreeMapView(TREE_MAP_VIEW_NAME, studyId));
  VIEWS[studyId][MEDICAL_HISTORY_VIEW_NAME] =
  <MedicalHistoryView>addView(new MedicalHistoryView(MEDICAL_HISTORY_VIEW_NAME, studyId));
  VIEWS[studyId][QUESTIONNAIRES_VIEW_NAME] =
  <QuestionnaiesView>addView(new QuestionnaiesView(QUESTIONNAIRES_VIEW_NAME, studyId));

  const tableViewHelpers = {};

  Object.keys(TABLE_VIEWS).forEach((it) => {
    const tableView = createTableView(
      TABLE_VIEWS[it].domainsAndColsToCheck,
      it,
      TABLE_VIEWS[it].helpUrl,
      TABLE_VIEWS[it].createViewHelper,
      TABLE_VIEWS[it].paramsForHelper,
    );
    VIEWS[studyId][it] = addView(tableView.view);
    tableViewHelpers[it] = tableView.helper;
  });

  DG.ObjectHandler.register(new AdverseEventHandler());

  const summary = VIEWS[studyId][SUMMARY_VIEW_NAME] as StudySummaryView;
  summary.load();
  const valView = addView(new ValidationView(VALIDATION_VIEW_NAME, studyId));
  summary.validationView = valView;
  VIEWS[studyId][VALIDATION_VIEW_NAME] = valView;

  VIEWS[studyId][STUDY_CONFIGURATIN_VIEW_NAME] =
    <StudyConfigurationView>addView(new StudyConfigurationView(STUDY_CONFIGURATIN_VIEW_NAME, studyId));

  setTimeout(() => {
    grok.shell.v = summary;
  }, 2000);

  const setObj = async (obj) => {
    grok.shell.o = await obj.propertyPanel();
  };

  const viewChangeSub = grok.events.onCurrentViewChanged.subscribe((v) => {
    setTimeout(() => {
      const studyViews = Object.keys(VIEWS[studyId]);
      if (!studyViews.length || studyViews.length === 1 && studyViews[AE_BROWSER_VIEW_NAME]) {
        viewChangeSub.unsubscribe();
        return;
      }
      const obj = VIEWS[studyId][grok.shell.v.name];
      if (obj) {
        if (obj.hasOwnProperty('loaded')) {
          if (!(obj as ClinicalCaseViewBase).loaded)
            (obj as ClinicalCaseViewBase).load();
        } else if (obj.type === DG.TYPE.TABLE_VIEW &&
          (obj as DG.TableView).dataFrame.name === studies[studyId].domains.ae.name)
          tableViewHelpers[AE_BROWSER_VIEW_NAME].propertyPanel();


        if ((obj as ClinicalCaseViewBase).loaded) {
          setObj(obj);
          if ((obj as ClinicalCaseViewBase).filterChanged) {
            (obj as ClinicalCaseViewBase).updateGlobalFilter();
            (obj as ClinicalCaseViewBase).filterChanged = false;
          }
        }
      }
    }, 100);
  });

  const viewClosedSub = grok.events.onViewRemoved.subscribe((v) => {
    if (VIEWS[studyId][v.name])
      delete VIEWS[studyId][v.name];
    if (!Object.keys(VIEWS[studyId]).length) {
      viewClosedSub.unsubscribe();
      return;
    }
  });

  if (studies[studyId].domains.dm) {
    const updateFilterProperty = (obj) => {
      if (obj.name === grok.shell.v.name) {
        setTimeout(() => {
          obj.updateGlobalFilter();
        });
      }
      obj.filterChanged = true;
    };
    studies[studyId].domains.dm.onFilterChanged.subscribe(()=> {
      Object.values(VIEWS[studyId]).forEach((it) => {
        if (it.hasOwnProperty('filterChanged'))
          updateFilterProperty(it);
        else {
          if (tableViewHelpers[it.name].updateGlobalFilter !== undefined)
            updateFilterProperty(tableViewHelpers[it.name]);
        }
      });
    });
  }
}

//name: clinicalCaseFolderLauncher
//tags: folderViewer
//input: file folder
//input: list<file> files
//output: widget res
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]):
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
  [STUDY_CONFIGURATIN_VIEW_NAME]: (studyId) => new StudyConfigurationView(STUDY_CONFIGURATIN_VIEW_NAME, studyId),
  [VALIDATION_VIEW_NAME]: (studyId) => new ValidationView(VALIDATION_VIEW_NAME, studyId),
  //[COHORT_VIEW_NAME]: (studyId) => new StudySummaryView(COHORT_VIEW_NAME, studyId),
  [QUESTIONNAIRES_VIEW_NAME]: (studyId) => new QuestionnaiesView(QUESTIONNAIRES_VIEW_NAME, studyId),
  [AE_BROWSER_VIEW_NAME]: (studyId, viewName) => createClinCaseTableView(studyId, viewName),
};


