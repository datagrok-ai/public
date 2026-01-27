/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StudySummaryView} from './views/study-summary-view';
import {TimelinesView} from './views/timelines-view';
import {PatientProfileView} from './views/patient-profile-view';
import {AdverseEventsView} from './views/adverse-events-view';
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
import {ConfigurationView} from './views/configuration-view';
import {ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME,
  ANIMAL_PROFILE_VIEW_NAME,
  CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MATRIX_TABLE_VIEW_NAME,
  MEDICAL_HISTORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME,
  CONFIGURATION_VIEW_NAME, SUMMARY_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME,
  MEASUREMENT_PROFILE_TABLE_VIEW_NAME,
  TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME,
  VALIDATION_VIEW_NAME, VISITS_VIEW_NAME} from './constants/view-names-constants';
import {createClinCaseTableView} from './utils/views-creation-utils';
//import {CohortView} from './views/cohort-view';
import {QuestionnaiesView} from './views/questionnaires-view';
import {CDISC_STANDARD, ClinCaseTableView, ClinStudyConfig} from './utils/types';
import {ClinicalStudy} from './clinical-study';
import '../css/clinical-case.css';
import {scripts} from './package-api';
import dayjs from 'dayjs';
import {createInitialSatistics} from './utils/initial-statistics-widget';
import {cdiscAppTB, CLINICAL_CASE_APP_PATH, createStudiesFromAppData,
  openApp, PRECLINICAL_CASE_APP_PATH,
  studies} from './utils/app-utils';
// Import renderer to ensure it's registered
import './utils/rule-violation-cell-renderer';

export * from './package.g';

export const _package = new DG.Package();

//export let study: ClinicalStudy = new ClinicalStudy();

type CurrentStudyAndView = {
  study: string;
  standard: string;
  viewName: string;
}

export let c: DG.FuncCall;

function getCurrentStudyAndView(path: string): CurrentStudyAndView {
  let currentStandard = '';
  let currentStudy = '';
  let currentViewName = '';
  if (path) {
    const pathSegments = path.split('/');
    currentStandard = pathSegments[1];
    currentStudy = pathSegments.length > 2 ? decodeURI(pathSegments[2]) : '';
    currentViewName = pathSegments.length > 3 ? decodeURI(pathSegments[3]) : SUMMARY_VIEW_NAME;
  }
  return {study: currentStudy, viewName: currentViewName, standard: currentStandard};
}


export class PackageFunctions {
  @grok.decorators.app({
    'name': 'Clinical Case',
    'icon': '/img/clin_case_icon.png',
  })
  static async clinicalCaseApp(): Promise<DG.ViewBase | void> {
    return await openApp(CDISC_STANDARD.SDTM);
  }

  @grok.decorators.app({
    'name': 'Preclinical Case',
    'icon': '/img/preclinical_case_icon.png',
  })
  static async PreclinicalCaseApp(): Promise<DG.ViewBase | void> {
    return await openApp(CDISC_STANDARD.SEND);
  }

  @grok.decorators.appTreeBrowser({app: 'Clinical Case'})
  static async clinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    const url = new URL(window.location.href);
    const currentStudyAndViewPath = url.pathname.includes(`${CLINICAL_CASE_APP_PATH}`) ?
      url.pathname.replace(`${CLINICAL_CASE_APP_PATH}`, `/${CDISC_STANDARD.SDTM}`) : '';
    const studyAndView = getCurrentStudyAndView(currentStudyAndViewPath);
    await cdiscAppTB(treeNode, CDISC_STANDARD.SDTM, studyAndView.study, studyAndView.viewName);
  }

  @grok.decorators.appTreeBrowser({app: 'Preclinical Case'})
  static async preclinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    const url = new URL(window.location.href);
    const currentStudyAndViewPath = url.pathname.includes(`${PRECLINICAL_CASE_APP_PATH}`) ?
      url.pathname.replace(`${PRECLINICAL_CASE_APP_PATH}`, `/${CDISC_STANDARD.SEND}`) : '';
    const studyAndView = getCurrentStudyAndView(currentStudyAndViewPath);
    await cdiscAppTB(treeNode, CDISC_STANDARD.SEND, studyAndView.study, studyAndView.viewName);
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
      .getOrCreateGroup(standard === CDISC_STANDARD.SEND ? 'Preclinical Case' : 'Clinical Case');
    await createStudiesFromAppData(clinicalCaseNode, standard ?? CDISC_STANDARD.SDTM);
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
          const domainsNames = Object.keys(studies[studyId].domains).filter((it) => it !== 'supp');
          await Promise.all(files.map(async (file) => {
            if (domainsNames.includes(file.fileName.toLowerCase())) {
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

  @grok.decorators.func({
    'name': 'Run CDISC CORE Validation',
    'description': 'Run CDISC CORE validation on datasets',
  })
  static async runCoreValidate(
    @grok.decorators.param({
      'name': 'standard',
      'description': 'CDISC standard (e.g., sendig, sdtmig)',
    }) standard: string,
    @grok.decorators.param({
      'name': 'dataPath',
      'description': 'Path to datasets directory (local path in container or remote System:AppData path)',
    }) dataPath: string,
    @grok.decorators.param({
      'name': 'version',
      'description': 'Standard version (e.g., 3.1)',
    }) version?: string,
    @grok.decorators.param({
      'name': 'outputFormat',
      'description': 'Output format (json, csv, xlsx)',
      'optional': true,
    }) outputFormat?: string,
    @grok.decorators.param({
      'name': 'options',
      'description': 'Additional validation options as JSON string',
      'optional': true,
    }) options?: any,
  ): Promise<string> {
    const result = await grok.functions.call('ClinicalCase:run_core_validate', {
      standard,
      version,
      data_path: dataPath,
      output_format: outputFormat || 'json',
      options: options ? JSON.stringify(options) : undefined,
    });
    return result;
  }
}

export const SUPPORTED_VIEWS: {[key: string]: string[]} = {
  [CDISC_STANDARD.SDTM]: [SUMMARY_VIEW_NAME, TIMELINES_VIEW_NAME, LABORATORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME,
    ADVERSE_EVENTS_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME,
    CORRELATIONS_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME,
    VISITS_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, VALIDATION_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME,
    AE_BROWSER_VIEW_NAME],

  [CDISC_STANDARD.SEND]: [SUMMARY_VIEW_NAME, TIMELINES_VIEW_NAME, LABORATORY_VIEW_NAME, ANIMAL_PROFILE_VIEW_NAME,
    DISTRIBUTIONS_VIEW_NAME, VALIDATION_VIEW_NAME, MATRIX_TABLE_VIEW_NAME, MEASUREMENT_PROFILE_TABLE_VIEW_NAME,
    CONFIGURATION_VIEW_NAME],
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
  [CONFIGURATION_VIEW_NAME]: (studyId) => new ConfigurationView(CONFIGURATION_VIEW_NAME, studyId),
  [VALIDATION_VIEW_NAME]: (studyId) => createClinCaseTableView(studyId, VALIDATION_VIEW_NAME),
  [QUESTIONNAIRES_VIEW_NAME]: (studyId) => new QuestionnaiesView(QUESTIONNAIRES_VIEW_NAME, studyId),
  [AE_BROWSER_VIEW_NAME]: (studyId, viewName) => createClinCaseTableView(studyId, viewName),
  [MATRIX_TABLE_VIEW_NAME]: (studyId) => createClinCaseTableView(studyId, MATRIX_TABLE_VIEW_NAME),
  [MEASUREMENT_PROFILE_TABLE_VIEW_NAME]: (studyId) =>
    createClinCaseTableView(studyId, MEASUREMENT_PROFILE_TABLE_VIEW_NAME),
};


