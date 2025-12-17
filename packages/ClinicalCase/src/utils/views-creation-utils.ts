import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {AE_BROWSER_VIEW_NAME, TIMELINES_VIEW_NAME, VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import * as sdtmCols from '../constants/columns-constants';
import {AE_START_DAY_FIELD} from '../views-config';
import {AEBrowserHelper} from '../helpers/ae-browser-helper';
import {ValidationHelper} from '../helpers/validation-helper';
import {c} from '../package';
import {createValidationErrorsDiv, getRequiredColumnsByView, setupValidationErrorColumns, setupValidationErrorIndicators} from './views-validation-utils';
import {updateDivInnerHTML} from './utils';
import {ClinCaseTableView} from './types';
import {studies} from './app-utils';
import {ValidationResult} from '../types/validation-result';

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

export function createValidationView(studyId: string): any {
  if (!studies[studyId].validated) {
    grok.shell.warning('Validation in progress...');
    studies[studyId].validationCompleted.subscribe(() => {
      grok.shell.info(`Validation for study ${studyId} completed`);
    });
    return;
  }

  const validationResults: ValidationResult = studies[studyId].validationResults;

  if (!validationResults) {
    grok.shell.error('No validation results available');
    return;
  }
  const issueDetailsDf = DG.DataFrame.fromObjects(validationResults.Issue_Details);
  const issueSummaryDf = DG.DataFrame.fromObjects(validationResults.Issue_Summary);

  grok.data.linkTables(issueSummaryDf, issueDetailsDf,
    ['dataset', 'core_id'], ['dataset', 'core_id'],
    [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]);

  // tv.dockManager.dock(issueDetailsDf.plot.grid(), 'right');
  const errorsGridDiv = ui.div();

  let dockNode: DG.DockNode | null = null;

  issueSummaryDf.onCurrentRowChanged.subscribe(() => {
    ui.empty(errorsGridDiv);
    if (issueSummaryDf.currentRowIdx === -1)
      return;
    if (!dockNode)
      dockNode = grok.shell.tv.dockManager.dock(errorsGridDiv, 'down');
    let grid: DG.Grid | null = null;
    //soem general rule violated - show row from issueDetailsDf
    const firstIssueDetailsIdx = issueDetailsDf.filter.findNext(-1, true);
    if (issueDetailsDf.filter.trueCount === 1 &&
      validationResults.Issue_Details[firstIssueDetailsIdx].row === '')
      grid = issueDetailsDf.plot.grid();
    else { //collect data from corresponding domain
      const domain: string = issueSummaryDf.get('dataset', issueSummaryDf.currentRowIdx);
      const domainWithoutExtension = domain.replace('.xpt', '').replace('.csv', '');
      const domainDf: DG.DataFrame = studies[studyId].domains[domainWithoutExtension];
      if (domainDf) {
        let columnName = '';
        //find the first variable name to scroll grid to corresponding column
        for (let i = 0; i < validationResults.Issue_Details[firstIssueDetailsIdx].variables.length; i++) {
          if (validationResults.Issue_Details[firstIssueDetailsIdx].values[i].toLowerCase() !== 'not in dataset') {
            columnName = validationResults.Issue_Details[firstIssueDetailsIdx].variables[i];
            break;
          }
        }
        validationResults.Issue_Details[issueDetailsDf.filter.findNext(-1, true)].variables
        if (!domainDf.col('~row'))
          domainDf.columns.addNewInt('~row').init((i) => i + 1);

        domainDf.filter.setAll(false);
        for (let i = 0; i < issueDetailsDf.rowCount; i++) {
          if (issueDetailsDf.filter.get(i)) {
            const idx = issueDetailsDf.get('row', i) - 1;
            domainDf.filter.set(idx, true);
          }
        }
        grid = domainDf.plot.grid();
        setupValidationErrorColumns(domainDf);
        setupValidationErrorIndicators(grid, domainDf, issueSummaryDf.get('core_id', issueSummaryDf.currentRowIdx));
        if (domainDf.col(columnName))
          grid.scrollToCell(columnName, 0);
      }
    }
    grid.root.style.width = '100%';
    grid.root.style.height = '95%';
    errorsGridDiv.append(grid.root);
  });
  return {df: issueSummaryDf};
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
  createViewHelper: (studyId: string, params: any) => any,
  paramsForHelper?: any) {
  let tableView: DG.TableView | DG.View;
  let viewHelper;
  const validator = new ValidationHelper(domainsAndColsToCheck, studyId);
  if (validator.validate()) {
    const {helper, df} = createViewHelper(studyId, paramsForHelper);
    tableView = DG.TableView.create(df, false);
    viewHelper = helper;
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
};
