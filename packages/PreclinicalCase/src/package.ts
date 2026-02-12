/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StudySummaryView} from './views/study-summary-view';
import {ConfigurationView} from './views/configuration-view';
import {STUDY_ID} from './constants/columns-constants';
import {CONFIGURATION_VIEW_NAME, EVENTS_VIEW_NAME, MATRIX_TABLE_VIEW_NAME,
  MEASUREMENT_PROFILE_TABLE_VIEW_NAME, MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME,
  SUMMARY_VIEW_NAME, VALIDATION_VIEW_NAME} from './constants/view-names-constants';
import {createTableView} from './utils/views-creation-utils';
import {PreclinicalStudy} from './preclinical-study';
import {TableView} from './types/types';
import '../css/preclinical-case.css';
import {scripts} from './package-api';
import {cdiscAppTB, createStudiesFromAppData,
  openApp, PRECLINICAL_CASE_APP_PATH,
  studies} from './utils/app-utils';
// Import renderers to ensure they're registered
import './utils/rule-violation-cell-renderer';
import './utils/combined-measurements-cell-renderer';
import { createValidationView } from './views/validation-table-view';
import { createMatrixTableView } from './views/matrix-table-view';
import { createMeasurementProfileTableView } from './views/measurement-profile-table-view';
import { createMICrossDomainView } from './views/mi-cross-domain-analysis';
import { createEventsView } from './views/events-view';
// Import renderers to ensure they're registered                                                                                                                    
import './utils/rule-violation-cell-renderer';                                                                                                                      
import './utils/combined-measurements-cell-renderer';   
import './utils/y-axis-cell-renderer';

export * from './package.g';

export const _package = new DG.Package();

type CurrentStudyAndView = {
  study: string;
  viewName: string;
}

export let c: DG.FuncCall;

function getCurrentStudyAndView(path: string): CurrentStudyAndView {
  let currentStudy = '';
  let currentViewName = '';
  if (path) {
    const pathSegments = path.split('/');
    currentStudy = pathSegments.length > 1 ? decodeURI(pathSegments[1]) : '';
    currentViewName = pathSegments.length > 2 ? decodeURI(pathSegments[2]) : SUMMARY_VIEW_NAME;
  }
  return {study: currentStudy, viewName: currentViewName};
}


export class PackageFunctions {
  @grok.decorators.app({
    'name': 'Preclinical Case',
    'icon': '/img/preclinical_case_icon.png',
  })
  static async PreclinicalCaseApp(): Promise<DG.ViewBase | void> {
    return await openApp();
  }

  @grok.decorators.appTreeBrowser({app: 'Preclinical Case'})
  static async preclinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    const url = new URL(window.location.href);
    const currentStudyAndViewPath = url.pathname.includes(`${PRECLINICAL_CASE_APP_PATH}`) ?
      url.pathname.replace(`${PRECLINICAL_CASE_APP_PATH}`, '') : '';
    const studyAndView = getCurrentStudyAndView(currentStudyAndViewPath);
    await cdiscAppTB(treeNode, studyAndView.study, studyAndView.viewName);
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
    'description': 'Run CDISC CORE validation on SEND datasets',
  })
  static async runCoreValidate(
    @grok.decorators.param({
      'name': 'standard',
      'description': 'CDISC standard (e.g., sendig)',
    }) standard: string,
    @grok.decorators.param({
      'name': 'dataPath',
      'description': 'Path to datasets directory',
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
    const result = await grok.functions.call('PreclinicalCase:run_core_validate', {
      standard,
      version,
      data_path: dataPath,
      output_format: outputFormat || 'json',
      options: options ? JSON.stringify(options) : undefined,
    });
    return result;
  }
}

export const SUPPORTED_VIEWS: string[] = [
  SUMMARY_VIEW_NAME, VALIDATION_VIEW_NAME, MATRIX_TABLE_VIEW_NAME,
  MEASUREMENT_PROFILE_TABLE_VIEW_NAME, MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME,
  EVENTS_VIEW_NAME, CONFIGURATION_VIEW_NAME,
];

export const VIEW_CREATE_FUNC: {[key: string]: (studyId: string, args?: any) => DG.ViewBase | TableView} = {
  [SUMMARY_VIEW_NAME]:
    (studyId, errorLinkHandler?: () => void) => new StudySummaryView(SUMMARY_VIEW_NAME, studyId, errorLinkHandler),
  [CONFIGURATION_VIEW_NAME]: (studyId) => new ConfigurationView(CONFIGURATION_VIEW_NAME, studyId),
  [VALIDATION_VIEW_NAME]: (studyId) => createTableView(studyId, VALIDATION_VIEW_NAME, createValidationView),
  [MATRIX_TABLE_VIEW_NAME]: (studyId) => createTableView(studyId, MATRIX_TABLE_VIEW_NAME, createMatrixTableView),
  [MEASUREMENT_PROFILE_TABLE_VIEW_NAME]: (studyId) =>
    createTableView(studyId, MEASUREMENT_PROFILE_TABLE_VIEW_NAME, createMeasurementProfileTableView),
  [MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME]: (studyId) =>
    createTableView(studyId, MICROSCOPIC_FINDINGS_TABLE_VIEW_NAME, createMICrossDomainView),
  [EVENTS_VIEW_NAME]: (studyId) =>
    createTableView(studyId, EVENTS_VIEW_NAME, createEventsView),
};
