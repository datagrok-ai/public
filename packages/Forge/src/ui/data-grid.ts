import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const MAX_VISIBLE_ROWS = 15;
// By the grid's root: a grid found again from its element is a new wrapper of the same root.
const TOOLTIPS = new WeakMap<HTMLElement, GridTooltips>();

export interface GridTooltips {
  /** Header tooltips by column name. */
  headers?: {[column: string]: string};
  /** Cell tooltips by column name, one line per `\n`; without one, or when it returns '', the cell's text. */
  cells?: {[column: string]: (row: number) => string};
}

/** The platform grid over [df] for a pane or a section: read-only, no row header, as tall as its rows (then it
 * scrolls), as wide as its host, with [tooltips] on the headers and the cells. */
export function readOnlyGrid(df: DG.DataFrame, tooltips: GridTooltips = {}): DG.Grid {
  const grid = DG.Viewer.grid(df);
  grid.props.allowEdit = false;
  grid.props.showRowHeader = false;
  grid.props.showCurrentRowIndicator = false;
  grid.props.showAddNewRowIcon = false;
  // Inline, computed: one header row and the rows, so a short table shows no scroll bar.
  grid.root.style.width = '100%';
  grid.root.style.height = `${(Math.min(df.rowCount, MAX_VISIBLE_ROWS) + 1) * grid.props.rowHeight + 2}px`;
  TOOLTIPS.set(grid.root, tooltips);
  grid.onCellTooltip((cell, x, y) => {
    const text = gridTooltip(grid, cell);
    if (text === '')
      return false;
    ui.tooltip.show(ui.divV(text.split('\n').map((line) => ui.divText(line))), x, y);
    return true;
  });
  return grid;
}

/** The tooltip of [cell] in a {@link readOnlyGrid}; '' for none. */
export function gridTooltip(grid: DG.Grid, cell: DG.GridCell): string {
  const tooltips = TOOLTIPS.get(grid.root);
  const column = cell.tableColumn;
  if (tooltips === undefined || column === null)
    return '';
  if (cell.isColHeader)
    return tooltips.headers?.[column.name] ?? '';
  const row = cell.tableRowIndex;
  if (!cell.isTableCell || row === null)
    return '';
  return tooltips.cells?.[column.name]?.(row) || (column.isNone(row) ? '' : column.getString(row));
}

/** A string column of [values], a missing value empty. */
export function textColumn(name: string, values: (string | undefined)[]): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values.map((value) => value ?? ''));
}
