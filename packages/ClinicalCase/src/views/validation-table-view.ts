import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ValidationResult, IssueDetail} from '../types/validation-result';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {reloadValidationView, studies} from '../utils/app-utils';
import {validationFixFunctions} from '../utils/validation-fix-utils';
import {setupValidationErrorColumns, setupValidationErrorIndicators} from '../utils/views-validation-utils';

/** True for cells that represent "no specific row" — empty string, null, NaN. */
function isEmptyCell(v: any): boolean {
  return v == null || v === '' || (typeof v === 'number' && Number.isNaN(v));
}

/** Reads an array-valued cell. CSVs stringify arrays as JSON; the d42 cache keeps the
 * JSON string verbatim. Fallback path (old JSON report) stored arrays directly.
 * @param {DG.DataFrame} df Dataframe holding the cell.
 * @param {string} col Column name.
 * @param {number} row Row index.
 * @return {string[]} Decoded array; empty when the cell is empty or unparseable. */
function getArrayCell(df: DG.DataFrame, col: string, row: number): string[] {
  const v = df.col(col)?.get(row);
  if (v == null) return [];
  if (Array.isArray(v)) return v;
  if (typeof v === 'string' && v.length) {
    try { return JSON.parse(v); } catch { return []; }
  }
  return [];
}

export function createValidationView(studyId: string): any {
  if (!studies[studyId].validated) {
    // validate() is async — either we're streaming cached files in (seconds) or
    // the container is still running CORE (minutes). Return a placeholder view
    // with a persistent message so the user can navigate away and come back to
    // see what state it's in. The auto-rebuild fires when validate() settles.
    const sub = studies[studyId].validationCompleted.subscribe(() => {
      sub.unsubscribe();
      reloadValidationView(studyId);
    });
    const placeholderDf = DG.DataFrame.create();
    const onTableViewAdded = async (tableView: DG.TableView) => {
      if (tableView.grid?.root)
        tableView.grid.root.style.display = 'none';

      // Distinguish "reading cached files" (seconds) from "container running CORE" (minutes).
      let cacheExists = true;
      try {
        const base = `System:AppData/ClinicalCase/${studies[studyId].config.standard!}/${studyId}`;
        cacheExists = (await Promise.all([
          grok.dapi.files.exists(`${base}/validation_results.d42`),
          grok.dapi.files.exists(`${base}/issue_summary.csv`),
        ])).some(Boolean);
      } catch (_) { /* fall back to the optimistic message */ }

      const text = cacheExists ?
        'Loading validation results…' :
        'Running CORE validation in the container. This can take several minutes — the view ' +
          'will update automatically when results are ready.';
      tableView.root.appendChild(ui.divText(text, {style: {padding: '24px', maxWidth: '640px'}}));
    };
    return {df: placeholderDf, onTableViewAddedFunc: onTableViewAdded};
  }

  const validationResults: ValidationResult = studies[studyId].validationResults;
  // Dataframes were built once in validate() — from the d42 cache, the per-table
  // CSVs, or as a last-resort fallback from the old fat JSON. Reuse them here
  // instead of paying DG.DataFrame.fromObjects on the main thread.
  const issueSummaryDf = studies[studyId].issueSummaryDf;
  const issueDetailsDf = studies[studyId].issueDetailsDf;

  if (!validationResults && !issueSummaryDf) {
    grok.shell.error('No validation results available');
    return;
  }

  if (!issueSummaryDf || !issueSummaryDf.rowCount) {
    grok.shell.info(`No validation issues found for study ${studyId}`);
    return;
  }

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
    if (issueDetailsDf.filter.trueCount === 1 && isEmptyCell(issueDetailsDf.get('row', firstIssueDetailsIdx))) {
      grid = issueDetailsDf.plot.grid();
      header = 'Issue details';
      currentDomainDf = null;
    } else { //collect data from corresponding domain
      const domain: string = issueSummaryDf.get('dataset', issueSummaryDf.currentRowIdx);
      // CORE reports dataset names in uppercase (e.g. "AE"), but the in-memory
      // domain map uses lowercase keys (ae, dm, lb, …) — normalise before lookup.
      const domainWithoutExtension = domain.replace('.xpt', '').replace('.csv', '').toLowerCase();
      let domainDf: DG.DataFrame | null = null;
      if (domainWithoutExtension.startsWith('supp')) {
        const domainIdx = studies[studyId].domains.supp.findIndex((it) => it.name === domainWithoutExtension);
        if (domainIdx !== -1)
          domainDf = studies[studyId].domains.supp[domainIdx];
      } else
        domainDf = studies[studyId].domains[domainWithoutExtension];
      if (domainDf) {
        //check if the violated rule related to metadata (variables in rules do not contain domain variables)
        const variables = getArrayCell(issueDetailsDf, 'variables', firstIssueDetailsIdx);
        const values = getArrayCell(issueDetailsDf, 'values', firstIssueDetailsIdx);
        const colNames = domainDf.columns.names().map((it) => it.toLowerCase());
        const ruleContainsDomainCols = variables.filter((it) => colNames.includes(it.toLowerCase())).length;
        if (!ruleContainsDomainCols) {
          grid = issueDetailsDf.plot.grid();
          header = 'Issue details';
          currentDomainDf = null;
        } else {
          let columnName = '';
          //find the first variable name to scroll grid to corresponding column
          for (let i = 0; i < variables.length; i++) {
            if ((values[i] ?? '').toLowerCase() !== 'not in dataset') {
              columnName = variables[i];
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

    // Fallback: if no domain was matched, still show the issue details grid
    // rather than crashing on `grid.root` below.
    if (!grid) {
      grid = issueDetailsDf.plot.grid();
      header = 'Issue details';
      currentDomainDf = null;
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

    // Get affected issue details. Reconstruct IssueDetail records from the dataframe
    // — the validationResults metadata JSON no longer carries the Issue_Details array.
    const affectedIssueDetails: IssueDetail[] = [];
    for (let i = 0; i < issueDetailsDf.rowCount; i++) {
      if (issueDetailsDf.filter.get(i))
        affectedIssueDetails.push(studies[studyId].issueDetailFromRow(i));
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

    grok.shell.o = previewGrid.root;

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
        if (grok.shell.o === previewGrid.root)
          grok.shell.o = null;
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
