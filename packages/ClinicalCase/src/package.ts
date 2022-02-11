/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as meta from './sdtm-meta';
import {study} from "./clinical-study";
import {StudySummaryView} from "./views/study-summary-view";
import {TimelinesView} from "./views/timelines-view";
import {PatientProfileView} from "./views/patient-profile-view";
import {AdverseEventsView} from "./views/adverse-events-view";
import {ValidationView} from './views/validation-view';
import {AdverseEventHandler} from './panels/adverse-event-handler';
import {LaboratoryView} from './views/laboratory-view';
import {AERiskAssessmentView} from './views/ae-risk-assessment-view';
import {SurvivalAnalysisView} from './views/survival-analysis-view';
import { BoxPlotsView } from './views/boxplots-view';
import { MatrixesView } from './views/matrixes-view';
import { TimeProfileView } from './views/time-profile-view';
import { AEBrowserHelper } from './helpers/ae-browser-helper';
import { AE_END_DAY, AE_START_DAY, SUBJECT_ID } from './columns-constants';
import { STUDY_ID } from './columns-constants';
import { createValidationErrorsDiv, updateDivInnerHTML } from './views/utils';
import { TreeMapView } from './views/tree-map-view';
import { MedicalHistoryView } from './views/medical-history-view';
import { VisitsView } from './views/visits-view';
import { StudyConfigurationView } from './views/study-config-view';
import { ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, SUMMARY_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, VALIDATION_VIEW_NAME, VISITS_VIEW_NAME } from './view-names-constants';
import { AE_TERM_FIELD, VIEWS_CONFIG } from './views-config';
import { ValidationHelper } from './helpers/validation-helper';

export let _package = new DG.Package();

export let validationRulesList = null;

let domains = Object.keys(study.domains).map(it => `${it.toLocaleLowerCase()}.csv`);


//name: Clinical Case
//tags: app
export async function clinicalCaseApp(): Promise<any> {
  let c: DG.FuncCall = grok.functions.getCurrentCall();
  validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);

  if (Object.keys(study.domains).every((name) => grok.shell.table(name) == null)) {
    let demoFiles = await grok.dapi.projects.filter('clin-demo-files-2').list();
    if (demoFiles.length) {
      await grok.dapi.projects.open('clin-demo-files-2');
    } else {
      grok.shell.warning('Please load SDTM data or demo files');
    }
  }

  study.initFromWorkspace();

  function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/' + view.name;
    grok.shell.addView(view);
    return view;
  }

  function createTableView(
    domainsAndColsToCheck: any,
    viewName: string,
    helpUrl: string,
    createViewHelper: (params: any) => any,
    paramsForHelper = null) {
    let tableView;
    let viewHelper;
    let validator = new ValidationHelper(domainsAndColsToCheck);
    if(validator.validate()){
      let { helper, df } = createViewHelper(paramsForHelper);
      tableView = DG.TableView.create(df);
      viewHelper = helper
    } else {
      tableView = DG.View.create();
      updateDivInnerHTML(tableView.root, createValidationErrorsDiv(validator.missingDomains, validator.missingColumnsInReqDomains, validator.missingColumnsInOptDomains)); 
    }
    tableView.name = viewName;
    if (helpUrl) {
      tableView.helpUrl = helpUrl;
    }
    return { helper: viewHelper, view: tableView };
  }

  function createAEBrowserHelper(timelinesView: TimelinesView): any {
    const aeBrowserDf = study.domains.ae.clone();
    const aeBrowserHelper = new AEBrowserHelper(aeBrowserDf);
    timelinesView.aeBrowserHelper = aeBrowserHelper;
    aeBrowserDf.onCurrentRowChanged.subscribe(() => {
      aeBrowserHelper.currentSubjId = aeBrowserDf.get(SUBJECT_ID, aeBrowserDf.currentRowIdx);
      aeBrowserHelper.currentAeDay = aeBrowserDf.get(AE_START_DAY, aeBrowserDf.currentRowIdx);
      aeBrowserHelper.propertyPanel();
    })
    return { helper: aeBrowserHelper, df: aeBrowserDf };
  }

  const views = [];

  views.push(<StudySummaryView>addView(new StudySummaryView(SUMMARY_VIEW_NAME)));
  views.push(<VisitsView>addView(new VisitsView(VISITS_VIEW_NAME)));
  const timelinesView = new TimelinesView(TIMELINES_VIEW_NAME);
  views.push(<TimelinesView>addView(timelinesView));
  views.push(<PatientProfileView>addView(new PatientProfileView(PATIENT_PROFILE_VIEW_NAME)));
  views.push(<AdverseEventsView>addView(new AdverseEventsView(ADVERSE_EVENTS_VIEW_NAME)));
  views.push(<LaboratoryView>addView(new LaboratoryView(LABORATORY_VIEW_NAME)));
  views.push(<AERiskAssessmentView>addView(new AERiskAssessmentView(AE_RISK_ASSESSMENT_VIEW_NAME)));
  views.push(<SurvivalAnalysisView>addView(new SurvivalAnalysisView(SURVIVAL_ANALYSIS_VIEW_NAME)));
  views.push(<BoxPlotsView>addView(new BoxPlotsView(DISTRIBUTIONS_VIEW_NAME)));
  views.push(<MatrixesView>addView(new MatrixesView(CORRELATIONS_VIEW_NAME)));
  views.push(<TimeProfileView>addView(new TimeProfileView(TIME_PROFILE_VIEW_NAME)));
  views.push(<TreeMapView>addView(new TreeMapView(TREE_MAP_VIEW_NAME)));
  views.push(<MedicalHistoryView>addView(new MedicalHistoryView(MEDICAL_HISTORY_VIEW_NAME)));


   const aeBrowserView = createTableView(
    {
      'req_domains': {
        'ae': {
          'req': [VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_TERM_FIELD], AE_START_DAY, AE_END_DAY]
        }
      }
    },
    AE_BROWSER_VIEW_NAME,
    'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/ae_browser.md',
    createAEBrowserHelper,
    timelinesView
  );
  views.push(addView(aeBrowserView.view)); 

  DG.ObjectHandler.register(new AdverseEventHandler());

  let summary = views.find(it => it.name === SUMMARY_VIEW_NAME);
  summary.load();
  let valView = addView(new ValidationView(summary.errorsByDomain, VALIDATION_VIEW_NAME));
  summary.validationView = valView;
  views.push(valView);

  views.push(<StudyConfigurationView>addView(new StudyConfigurationView(STUDY_CONFIGURATIN_VIEW_NAME)));

  setTimeout(() => {
    grok.shell.v = summary;
  }, 1000);

  let setObj = async (obj) => {
    grok.shell.o = await obj.propertyPanel();
  }

  grok.events.onCurrentViewChanged.subscribe((v) => {
    setTimeout(() => {
      const obj = views.find(it => it.name === grok.shell.v.name);
      if (obj) {
        if (obj.hasOwnProperty('loaded') && !obj.loaded) {
          obj.load();
        }
        if (obj.loaded) {
          setObj(obj);
        }
      }
    }, 100)
  });
}


//tags: folderViewer
//input: file folder
//input: list<file> files
//output: widget
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]): Promise<DG.Widget | undefined> {
  if (files.some((f) => f.fileName.toLowerCase() === 'dm.csv')) {
    let res = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
    let table = DG.DataFrame.fromCsv(res);
    let studyId = table.columns.names().includes(STUDY_ID) ? table.get(STUDY_ID, 0) : 'undefined';
    return DG.Widget.fromRoot(ui.div([
      ui.panel([
        ui.divText('Folder contains SDTM data'),
        ui.divText(`Study ID: ${studyId}`)]),
      ui.button('Run ClinicalCase', async () => {
        await Promise.all(files.map(async (file) => {
          if(domains.includes(file.fileName.toLowerCase())){
            const df = await grok.data.files.openTable(`${folder.fullPath}/${file.fileName.toLowerCase()}`);
            grok.shell.addTableView(df);
          }
        }));
        grok.functions.call("Clinicalcase:clinicalCaseApp");
      })
    ]));
  }
}


