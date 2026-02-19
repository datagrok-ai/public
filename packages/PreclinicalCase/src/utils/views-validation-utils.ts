import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as sdtmCols from '../constants/columns-constants';
import {updateDivInnerHTML} from './utils';
import {Subscription} from 'rxjs';

const ERROR_ICON_SIZE = 9;
const ERROR_ICON_MARGIN = 2;
const COLUMNS_WITH_VALIDATION_ERRORS_TAG = 'columnsWithErrors';

export function setupValidationErrorColumns(df: DG.DataFrame) {
  if (df.getTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG))
    return;

  const columnsWithErrors: {[colName: string]: {hasErrorsCol: string, errorsCol: string}} = {};

  for (const colName of df.columns.names()) {
    if (colName.endsWith(sdtmCols.COL_HAS_ERRORS_POSTFIX)) {
      const baseColName = colName.replace(sdtmCols.COL_HAS_ERRORS_POSTFIX, '');
      const baseCol = df.columns.byName(baseColName);
      const errorsCol = df.columns.byName(`${baseColName}${sdtmCols.ERRORS_POSTFIX}`);

      if (baseCol && errorsCol) {
        errorsCol.semType = 'sdisc-rule-violation';
        columnsWithErrors[baseColName] = {
          hasErrorsCol: colName,
          errorsCol: `${baseColName}${sdtmCols.ERRORS_POSTFIX}`,
        };
      }
    }
  }

  if (Object.keys(columnsWithErrors).length === 0)
    df.setTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG, '');

  df.setTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG, JSON.stringify(columnsWithErrors));
}

let currentErrorCell: {
    tableColName: string,
    tableRowIndex: number,
    errors: Array<{ruleID: string, message: string, value: string}>,
    iconX: number,
    iconY: number,
    iconRight: number,
    iconBottom: number
  } | null = null;

export function setupValidationErrorIndicators(grid: DG.Grid, df: DG.DataFrame, ruleId?: string): Subscription[] {
  const columnsWithErrorsString = df.getTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG);
  if (!columnsWithErrorsString)
    return [];

  const columnsWithErrors = JSON.parse(columnsWithErrorsString);

  const onCellRenderedSub = grid.onCellRendered.subscribe((args: DG.GridCellRenderArgs) => {
    const cell = args.cell;
    const tableColName = cell.tableColumn?.name;
    const tableRowIndex = cell.tableRowIndex!;

    if (!cell.isTableCell || !tableColName || tableRowIndex === undefined || tableRowIndex === null)
      return;

    const errorColsNames: {hasErrorsCol: string, errorsCol: string} = columnsWithErrors[tableColName];
    if (!errorColsNames)
      return;
    const hasErrorsCol = df.col(errorColsNames.hasErrorsCol);
    const errorsCol = df.col(errorColsNames.errorsCol);
    if (!hasErrorsCol || !errorsCol)
      return;

    if (tableRowIndex < 0 || tableRowIndex >= hasErrorsCol.length)
      return;

    const hasErrors = hasErrorsCol.get(tableRowIndex);

    if (hasErrors) {
      const errorsStr = errorsCol.get(tableRowIndex);
      let errors: Array<{ruleID: string, message: string, value: string}> = [];
      if (errorsStr) {
        try {
          errors = JSON.parse(errorsStr);
        }
        catch (e) {
          console.error('Failed to parse column errors:', e);
        }
      }

      if (ruleId) {
        const filteredRule = errors.filter((it) => it.ruleID === ruleId);
        if (!filteredRule.length)
          return;
        errors = filteredRule;
      }

      const bounds = args.bounds;
      const g = args.g;

      if (!g || !bounds)
        return;

      g.save();
      g.beginPath();
      g.rect(bounds.x, bounds.y, bounds.width, bounds.height);
      g.clip();

      const iconSize = ERROR_ICON_SIZE;
      const margin = ERROR_ICON_MARGIN;
      const iconX = bounds.x + bounds.width - iconSize - margin;
      const iconY = bounds.y + margin;

      g.strokeStyle = '#dc3545';
      g.lineWidth = 1;
      g.beginPath();
      g.arc(iconX + iconSize / 2, iconY + iconSize / 2, iconSize / 2 - 0.5, 0, 2 * Math.PI);
      g.stroke();

      g.fillStyle = '#dc3545';
      g.font = 'bold 6px Arial';
      g.textAlign = 'center';
      g.textBaseline = 'middle';
      g.fillText('!', iconX + iconSize / 2, iconY + iconSize / 2 + 0.5);

      g.restore();
    }
  });

  const onCellMouseEnterSub = grid.onCellMouseEnter.subscribe((cell: DG.GridCell) => {
    if (!cell || !cell.isTableCell || !cell.tableColumn?.name) {
      currentErrorCell = null;
      return;
    }

    const errorColsNames: {hasErrorsCol: string, errorsCol: string} = columnsWithErrors[cell.tableColumn.name];
    if (!errorColsNames)
      return;
    const hasErrorsCol = df.col(errorColsNames.hasErrorsCol);
    const errorsCol = df.col(errorColsNames.errorsCol);
    if (!hasErrorsCol || !errorsCol) {
      currentErrorCell = null;
      return;
    }

    const errorsStr = errorsCol.get(cell.tableRowIndex!);
    let errors: Array<{ruleID: string, message: string, value: string}> = [];
    if (errorsStr) {
      try {
        errors = JSON.parse(errorsStr);
      }
      catch (e) {
        console.error('Failed to parse column errors:', e);
        currentErrorCell = null;
        return;
      }
    }

    if (ruleId) {
      const filteredRule = errors.filter((it) => it.ruleID === ruleId);
      if (!filteredRule.length)
        return;
      errors = filteredRule;
    }

    const cellBounds = cell.bounds;
    if (!cellBounds) {
      currentErrorCell = null;
      return;
    }

    const iconSize = ERROR_ICON_SIZE;
    const margin = ERROR_ICON_MARGIN;
    const iconX = cellBounds.x + cellBounds.width - iconSize - margin;
    const iconY = cellBounds.y + margin;
    const iconRight = iconX + iconSize;
    const iconBottom = iconY + iconSize;

    currentErrorCell = {tableColName: cell.tableColumn!.name,
      tableRowIndex: cell.tableRowIndex!, errors, iconX, iconY, iconRight, iconBottom};
  });

  grid.overlay.addEventListener('mousemove', handleMouseMoveOverErrorCell);

  const onCellMouseLeave = grid.onCellMouseLeave.subscribe(() => {
    currentErrorCell = null;
    ui.tooltip.hide();
  });

  return [onCellRenderedSub, onCellMouseEnterSub, onCellMouseLeave];
}

export function handleMouseMoveOverErrorCell(e: MouseEvent) {
  if (!currentErrorCell) {
    ui.tooltip.hide();
    return;
  }

  const tolerance = 3;
  if (e.offsetX >= currentErrorCell.iconX - tolerance && e.offsetX <= currentErrorCell.iconRight + tolerance &&
        e.offsetY >= currentErrorCell.iconY - tolerance && e.offsetY <= currentErrorCell.iconBottom + tolerance) {
    const tooltipContent = createColumnValidationTooltip(currentErrorCell.errors);
    ui.tooltip.show(tooltipContent, e.clientX, e.clientY);
  } else
    ui.tooltip.hide();
}

export function createColumnValidationTooltip(
  errors: Array<{ruleID: string, message: string, value: string}>,
): HTMLElement | string {
  if (!errors || errors.length === 0)
    return 'No validation errors';

  const tooltipDiv = ui.div([], {style: {maxWidth: '400px'}});
  errors.forEach((error, index) => {
    const errorDiv = ui.div([],
      {style: {
        marginBottom: '8px',
        paddingBottom: '8px',
        borderBottom: index < errors.length - 1 ? '1px solid var(--grey-2)' : 'none',
      }});

    const ruleIdDiv = ui.div([
      ui.span(['Rule ID: '], {style: {fontWeight: 'bold'}}),
      ui.span([error.ruleID]),
    ]);
    errorDiv.append(ruleIdDiv);

    const messageDiv = ui.div([
      ui.span(['Message: '], {style: {fontWeight: 'bold'}}),
      ui.span([error.message]),
    ]);
    errorDiv.append(messageDiv);

    if (error.value) {
      const valueDiv = ui.div([
        ui.span(['Value: '], {style: {fontWeight: 'bold'}}),
        ui.span([error.value]),
      ]);
      errorDiv.append(valueDiv);
    }

    tooltipDiv.append(errorDiv);
  });

  return tooltipDiv;
}
