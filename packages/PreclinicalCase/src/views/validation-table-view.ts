import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ValidationResult, IssueDetail} from '../types/validation-result';
import {awaitCheck} from '@datagrok-libraries/test/src/test';
import {studies} from '../utils/app-utils';
import {validationFixFunctions} from '../utils/validation-fix-utils';
import {setupValidationErrorColumns, setupValidationErrorIndicators} from '../utils/views-validation-utils';
import { StudyTableViewParams } from '../utils/views-creation-utils';

export function createValidationView(studyId: string): StudyTableViewParams {
  if (!studies[studyId].validated) {
    grok.shell.warning('Validation in progress...');
    studies[studyId].validationCompleted.subscribe(() => {
      grok.shell.info(`Validation for study ${studyId} completed`);
    });
    return {df: DG.DataFrame.create()};
  }

  const validationResults: ValidationResult | null = studies[studyId].validationResults;

  if (!validationResults) {
    grok.shell.error('No validation results available');
    return {df: DG.DataFrame.create()};
  }
  const issueDetailsDf = DG.DataFrame.fromObjects(validationResults.Issue_Details);
  const issueSummaryDf = DG.DataFrame.fromObjects(validationResults.Issue_Summary);

  issueSummaryDf!.columns.addNewString('action').init((i) => {
    const coreId = issueSummaryDf!.get('core_id', i);
    return (coreId && validationFixFunctions[coreId]) ? coreId : '';
  });

  grok.data.linkTables(issueSummaryDf!, issueDetailsDf!,
    ['dataset', 'core_id'], ['dataset', 'core_id'],
    [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]);

  const errorsGridDiv = ui.div();
  let dockNode: DG.DockNode | null = null;
  let currentDomainDf: DG.DataFrame | null = null;
  let applyFixesButton: HTMLElement | null = null;
  let currentTableView: DG.TableView | null = null;

  issueSummaryDf!.onCurrentRowChanged.subscribe(() => {
    ui.empty(errorsGridDiv);
    let header = '';
    if (issueSummaryDf!.currentRowIdx === -1)
      return;
    if (!dockNode)
      dockNode = grok.shell.tv.dockManager.dock(errorsGridDiv, 'down');
    let grid: DG.Grid | null = null;
    const firstIssueDetailsIdx = issueDetailsDf!.filter.findNext(-1, true);
    if (issueDetailsDf!.filter.trueCount === 1 &&
      validationResults.Issue_Details[firstIssueDetailsIdx].row === '') {
      grid = issueDetailsDf!.plot.grid();
      header = 'Issue details';
      currentDomainDf = null;
    } else {
      const domain: string = issueSummaryDf!.get('dataset', issueSummaryDf!.currentRowIdx);
      const domainWithoutExtension = domain.toLowerCase().replace('.xpt', '').replace('.csv', '');
      let domainDf: DG.DataFrame | null = null;
      if (domainWithoutExtension.startsWith('supp')) {
        const domainIdx = studies[studyId].domains.supp.findIndex((it) => it.name === domainWithoutExtension);
        if (domainIdx !== -1)
          domainDf = studies[studyId].domains.supp[domainIdx];
      } else 
        domainDf = (studies[studyId].domains as any)[domainWithoutExtension];
      if (domainDf) {
        const colNames = domainDf.columns.names().map((it) => it.toLowerCase());
        const ruleContainsDomainCols = validationResults.Issue_Details[firstIssueDetailsIdx].variables
          .filter((it) => colNames.includes(it.toLowerCase())).length;
        if (!ruleContainsDomainCols) {
          grid = issueDetailsDf!.plot.grid();
          header = 'Issue details';
          currentDomainDf = null;
        } else {
          let columnName = '';
          for (let i = 0; i < validationResults.Issue_Details[firstIssueDetailsIdx].variables.length; i++) {
            if (validationResults.Issue_Details[firstIssueDetailsIdx].values[i].toLowerCase() !== 'not in dataset') {
              columnName = validationResults.Issue_Details[firstIssueDetailsIdx].variables[i];
              break;
            }
          }
          if (!domainDf.col('~row'))
            domainDf.columns.addNewInt('~row').init((i) => i + 1);

          domainDf.filter.setAll(false);
          for (let i = 0; i < issueDetailsDf!.rowCount; i++) {
            if (issueDetailsDf!.filter.get(i)) {
              const idx = issueDetailsDf!.get('row', i) - 1;
              domainDf.filter.set(idx, true);
            }
          }
          grid = domainDf.plot.grid();
          header = domain;
          setupValidationErrorColumns(domainDf);
          setupValidationErrorIndicators(grid, domainDf, issueSummaryDf!.get('core_id', issueSummaryDf!.currentRowIdx));
          if (domainDf.col(columnName))
            grid.scrollToCell(columnName, 0);
          currentDomainDf = domainDf;
        }
      }
    }

    grid!.root.prepend(ui.h1(header, {style: {margin: '0px 0px 10px 10px'}}));
    grid!.root.style.width = '100%';
    grid!.root.style.height = '95%';
    errorsGridDiv.append(grid!.root);
  });

  const handleFixClick = (
    fixFunction: (df: DG.DataFrame, issueDetails: IssueDetail[]) => {df: DG.DataFrame, colsToFix: string[],
      colsOrder: string[]}, rowIdx: number,
  ) => {
    if (rowIdx !== currentTableView!.dataFrame.currentRowIdx) {
      grok.shell.warning('Select row to run fixes');
      return;
    }

    const affectedIssueDetails: IssueDetail[] = [];
    for (let i = 0; i < issueDetailsDf!.rowCount; i++) {
      if (issueDetailsDf!.filter.get(i))
        affectedIssueDetails.push(validationResults.Issue_Details[i]);
    }

    if (affectedIssueDetails.length === 0)
      return;

    if (currentDomainDf!.filter.trueCount === 0) {
      grok.shell.warning('No filtered rows available for preview');
      return;
    }

    const fixRes = fixFunction(currentDomainDf!, affectedIssueDetails);
    const previewDf = fixRes.df;
    const colsOrder = fixRes.colsOrder;
    const colsToFix = fixRes.colsToFix;

    const previewGrid = previewDf.plot.grid();
    previewGrid.root.prepend(ui.h1('Fixes Preview', {style: {margin: '0px 0px 10px 10px'}}));
    previewGrid.root.style.width = '100%';
    previewGrid.root.style.height = '95%';
    previewGrid.columns.setOrder(colsOrder);

    grok.shell.o = previewGrid.root;

    if (!applyFixesButton && currentTableView) {
      applyFixesButton = ui.button('Apply Fixes', () => {
        for (const colTofix of colsToFix) {
          let counter = 0;
          for (let i = 0; i < currentDomainDf!.filter.length; i++) {
            if (currentDomainDf!.filter.get(i)) {
              currentDomainDf!.col(colTofix)!.set(i, previewDf.get(`${colTofix}_fix`, counter));
              counter++;
            }
          }
        }
        if (grok.shell.o === previewGrid.root)
          grok.shell.o = null;
      });
      applyFixesButton.style.margin = '0 5px';
      currentTableView.setRibbonPanels([[applyFixesButton]]);
    }
  };

  const onTableViewAdded = (tableView: DG.TableView) => {
    awaitCheck(() => tableView.grid !== null, 'Validation view hasn\'t been added', 5000).then(() => {
      currentTableView = tableView;
      const fixCol = tableView.grid.columns.byName('action');
      fixCol!.cellType = 'html';

      tableView.grid.onCellPrepare((gc) => {
        if (gc.isTableCell && gc.gridColumn.name === 'action') {
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

  return {df: issueSummaryDf!, onTableViewAddedFunc: onTableViewAdded};
}
