import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AE_BROWSER_VIEW_NAME, MATRIX_TABLE_VIEW_NAME, MEASUREMENT_PROFILE_TABLE_VIEW_NAME, TIMELINES_VIEW_NAME,
  VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import * as sdtmCols from '../constants/columns-constants';
import {AE_START_DAY_FIELD} from '../views-config';
import {AEBrowserHelper} from '../helpers/ae-browser-helper';
import {ValidationHelper} from '../helpers/validation-helper';
import {c} from '../package';
import {createValidationErrorsDiv, getRequiredColumnsByView} from './views-validation-utils';
import {hideValidationColumns, updateDivInnerHTML} from './utils';
import {ClinCaseTableView} from './types';
import {studies} from './app-utils';
import {createValidationView} from '../views/validation-table-view';
import {createMatrixTableView} from '../views/matrix-table-view';
import {createMeasurementProfileTableView} from '../views/measurement-profile-table-view';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';


export function createAEBrowserHelper(studyId: string): any {
  const aeBrowserDf = studies[studyId].domains.ae.clone();
  const aeBrowserHelper = new AEBrowserHelper(aeBrowserDf, studyId);
  const timelinesView = grok.shell.view(TIMELINES_VIEW_NAME) as any;
  if (timelinesView)
    timelinesView.aeBrowserHelper = aeBrowserHelper;

  aeBrowserDf.onCurrentRowChanged.subscribe(() => {
    aeBrowserHelper.currentSubjId = aeBrowserDf.get(sdtmCols.SUBJECT_ID, aeBrowserDf.currentRowIdx);
    aeBrowserHelper.currentAeDay = aeBrowserDf
      .get(studies[studyId].viewsConfig.config[AE_BROWSER_VIEW_NAME][AE_START_DAY_FIELD],
        aeBrowserDf.currentRowIdx);
    aeBrowserHelper.propertyPanel();
  });
  return {helper: aeBrowserHelper, df: aeBrowserDf};
}

export function addView(view: DG.ViewBase): DG.ViewBase {
  view.box = true;
  view.parentCall = c;
  view.path = '/' + view.name;
  grok.shell.addView(view);
  return view;
}

export function createClinCaseTableView(studyId: string, viewName: string): ClinCaseTableView {
  const tableView = createTableView(
    studyId,
    getRequiredColumnsByView(studyId),
    viewName,
    TABLE_VIEWS_META[viewName].helpUrl,
    TABLE_VIEWS_META[viewName].createViewHelper,
  );
  return {view: tableView.view, helper: tableView.helper};
}

export function createTableView(
  studyId: string,
  domainsAndColsToCheck: any,
  viewName: string,
  helpUrl: string,
  createViewHelper: (studyId: string, params: any) => {df: DG.DataFrame, helper?: any,
    onTableViewAddedFunc?: (tv: DG.TableView) => Promise<any>},
  paramsForHelper?: any) {
  let tableView: DG.TableView | DG.View;
  let viewHelper;
  const validator = new ValidationHelper(domainsAndColsToCheck, studyId);
  if (validator.validate()) {
    const {helper, df, onTableViewAddedFunc} = createViewHelper(studyId, paramsForHelper);
    tableView = DG.TableView.create(df, false);
    viewHelper = helper;
    if (onTableViewAddedFunc)
      onTableViewAddedFunc(tableView as DG.TableView);
    //wait for grid to become available to set validation columns invisible
    awaitCheck(() => (tableView as DG.TableView).grid !== null, `${viewName} hasn't been added`, 10000)
      .then(() => hideValidationColumns(tableView as DG.TableView))
      .catch(() => {});
  } else {
    tableView = DG.View.create();
    updateDivInnerHTML(tableView.root, createValidationErrorsDiv(validator.missingDomains,
      validator.missingColumnsInReqDomains, validator.missingColumnsInOptDomains));
  }
  tableView.name = viewName;
  if (helpUrl)
    tableView.helpUrl = helpUrl;
  grok.shell.windows.showHelp = true;

  return {helper: viewHelper, view: tableView};
}

export const TABLE_VIEWS_META = {
  [AE_BROWSER_VIEW_NAME]: {
    // eslint-disable-next-line max-len
    helpUrl: 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/ae_browser.md',
    createViewHelper: createAEBrowserHelper,
  },
  [VALIDATION_VIEW_NAME]: {
    createViewHelper: createValidationView,
  },
  [MATRIX_TABLE_VIEW_NAME]: {
    createViewHelper: createMatrixTableView,
  },
  [MEASUREMENT_PROFILE_TABLE_VIEW_NAME]: {
    createViewHelper: createMeasurementProfileTableView,
  },
};
