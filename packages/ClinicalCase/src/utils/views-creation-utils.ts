import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {AE_BROWSER_VIEW_NAME, TIMELINES_VIEW_NAME, TRELLIS_PLOT_VIEW_NAME,
  VALIDATION_VIEW_NAME} from '../constants/view-names-constants';
import * as sdtmCols from '../constants/columns-constants';
import {AE_START_DAY_FIELD} from '../views-config';
import {AEBrowserHelper} from '../helpers/ae-browser-helper';
import {ValidationHelper} from '../helpers/validation-helper';
import {c} from '../package';
import {createValidationErrorsDiv, getRequiredColumnsByView, setupValidationErrorColumns,
  setupValidationErrorIndicators} from './views-validation-utils';
import {updateDivInnerHTML} from './utils';
import {CDISC_STANDARD, ClinCaseTableView} from './types';
import {studies} from './app-utils';
import {ValidationResult, IssueDetail} from '../types/validation-result';
import {SUBJECT_ID} from '../constants/columns-constants';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {validationFixFunctions} from './validation-fix-utils';
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

  // Add 'action' column to issue summary
  // Store core_id in the cell if a fix function exists, otherwise leave empty
  issueSummaryDf.columns.addNewString('action').init((i) => {
    const coreId = issueSummaryDf.get('core_id', i);
    return (coreId && validationFixFunctions[coreId]) ? coreId : '';
  });

  grok.data.linkTables(issueSummaryDf, issueDetailsDf,
    ['dataset', 'core_id'], ['dataset', 'core_id'],
    [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]);

  const errorsGridDiv = ui.div();
  let dockNode: DG.DockNode | null = null;
  let currentDomainDf: DG.DataFrame | null = null;
  let applyFixesButton: HTMLElement | null = null;
  let currentTableView: DG.TableView | null = null;
  let previewDockNode: DG.DockNode | null = null;


  issueSummaryDf.onCurrentRowChanged.subscribe(() => {
    ui.empty(errorsGridDiv);
    let header = '';
    if (issueSummaryDf.currentRowIdx === -1)
      return;
    if (!dockNode)
      dockNode = grok.shell.tv.dockManager.dock(errorsGridDiv, 'down');
    let grid: DG.Grid | null = null;
    //a general rule violated - show row from issueDetailsDf
    const firstIssueDetailsIdx = issueDetailsDf.filter.findNext(-1, true);
    if (issueDetailsDf.filter.trueCount === 1 &&
      validationResults.Issue_Details[firstIssueDetailsIdx].row === '') {
      grid = issueDetailsDf.plot.grid();
      header = 'Issue details';
      currentDomainDf = null;
    } else { //collect data from corresponding domain
      const domain: string = issueSummaryDf.get('dataset', issueSummaryDf.currentRowIdx);
      const domainWithoutExtension = domain.replace('.xpt', '').replace('.csv', '');
      let domainDf: DG.DataFrame | null = null;
      if (domainWithoutExtension.startsWith('supp')) {
        const domainIdx = studies[studyId].domains.supp.findIndex((it) => it.name === domainWithoutExtension);
        if (domainIdx !== -1)
          domainDf = studies[studyId].domains.supp[domainIdx];
      } else
        domainDf = studies[studyId].domains[domainWithoutExtension];
      if (domainDf) {
        //check if the violated rule related to metadata (variables in rules do not contain domain variables)
        const colNames = domainDf.columns.names().map((it) => it.toLowerCase());
        const ruleContainsDomainCols = validationResults.Issue_Details[firstIssueDetailsIdx].variables
          .filter((it) => colNames.includes(it.toLowerCase())).length;
        if (!ruleContainsDomainCols) {
          grid = issueDetailsDf.plot.grid();
          header = 'Issue details';
          currentDomainDf = null;
        } else {
          let columnName = '';
          //find the first variable name to scroll grid to corresponding column
          for (let i = 0; i < validationResults.Issue_Details[firstIssueDetailsIdx].variables.length; i++) {
            if (validationResults.Issue_Details[firstIssueDetailsIdx].values[i].toLowerCase() !== 'not in dataset') {
              columnName = validationResults.Issue_Details[firstIssueDetailsIdx].variables[i];
              break;
            }
          }
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
          header = domain;
          setupValidationErrorColumns(domainDf);
          setupValidationErrorIndicators(grid, domainDf, issueSummaryDf.get('core_id', issueSummaryDf.currentRowIdx));
          if (domainDf.col(columnName))
            grid.scrollToCell(columnName, 0);
          currentDomainDf = domainDf;
        }
      }
    }

    grid.root.prepend(ui.h1(header, {style: {margin: '0px 0px 10px 10px'}}));
    grid.root.style.width = '100%';
    grid.root.style.height = '95%';
    errorsGridDiv.append(grid.root);
  });

  // Handler for fix button click
  const handleFixClick = (
    fixFunction: (df: DG.DataFrame, issueDetails: IssueDetail[]) => {df: DG.DataFrame, colsToFix: string[],
      colsOrder: string[]}, rowIdx: number,
  ) => {
    if (rowIdx !== currentTableView.dataFrame.currentRowIdx) {
      grok.shell.warning('Select row to run fixes');
      return;
    }

    // Get affected issue details
    const affectedIssueDetails: IssueDetail[] = [];
    for (let i = 0; i < issueDetailsDf.rowCount; i++) {
      if (issueDetailsDf.filter.get(i))
        affectedIssueDetails.push(validationResults.Issue_Details[i]);
    }

    if (affectedIssueDetails.length === 0)
      return;


    // Create preview dataframe: clone filtered rows from original dataframe
    // Use the filter BitSet directly
    if (currentDomainDf.filter.trueCount === 0) {
      grok.shell.warning('No filtered rows available for preview');
      return;
    }


    const fixRes = fixFunction(currentDomainDf, affectedIssueDetails);
    const previewDf = fixRes.df;
    const colsOrder = fixRes.colsOrder;
    const colsToFix = fixRes.colsToFix;

    // Create preview grid and dock it on the right
    const previewGrid = previewDf.plot.grid();
    previewGrid.root.prepend(ui.h1('Fixes Preview', {style: {margin: '0px 0px 10px 10px'}}));
    previewGrid.root.style.width = '100%';
    previewGrid.root.style.height = '95%';
    previewGrid.columns.setOrder(colsOrder);

    // Close existing preview if any
    if (previewDockNode) {
      grok.shell.tv.dockManager.close(previewDockNode);
      previewDockNode = null;
    }

    // Dock preview grid on the right
    previewDockNode = grok.shell.tv.dockManager.dock(previewGrid.root, 'right');

    // Show 'Apply Fixes' button in ribbon
    if (!applyFixesButton && currentTableView) {
      applyFixesButton = ui.button('Apply Fixes', () => {
        //apply fixes
        for (const colTofix of colsToFix) {
          let counter = 0;
          for (let i = 0; i < currentDomainDf.filter.length; i++) {
            if (currentDomainDf.filter.get(i)) {
              currentDomainDf.col(colTofix).set(i, previewDf.get(`${colTofix}_fix`, counter));
              counter++;
            }
          }
        }
        //close fixesPreview
        grok.shell.tv.dockManager.close(previewDockNode);
        previewDockNode = null;
      });
      applyFixesButton.style.margin = '0 5px';

      // Add button to ribbon
      const panels = currentTableView.getRibbonPanels();
      if (panels.length === 0)
        panels.push([]);
      panels[0].push(applyFixesButton);
      currentTableView.setRibbonPanels(panels);
    }
  };

  // onTableViewAdded function to set up fix button rendering
  const onTableViewAdded = (tableView: DG.TableView) => {
    //need to wait for grid to be created
    awaitCheck(() => tableView.grid !== null, `Validation view hasn't been added`, 5000).then(() => {
      currentTableView = tableView;
      const fixCol = tableView.grid.columns.byName('action');
      fixCol.cellType = 'html';

      tableView.grid.onCellPrepare((gc) => {
        if (gc.isTableCell && gc.gridColumn.name === 'action') {
          // Get core_id from the fix column cell value
          const coreId = gc.cell.value;
          const fixFunction = coreId ? validationFixFunctions[coreId] : null;

          if (fixFunction) {
            const button = ui.button('Fix', () => {
              handleFixClick(fixFunction, gc.cell.rowIndex);
            });
            button.style.marginTop = '0px';
            gc.style.element = button;
          } else
            gc.style.element = ui.divText('');
        }
      });
    });
  };

  return {df: issueSummaryDf, onTableViewAddedFunc: onTableViewAdded};
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
    onTableViewAddedFunc?: (tv: DG.TableView) => any},
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
