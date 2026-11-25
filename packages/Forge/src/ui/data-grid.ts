import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const MAX_VISIBLE_ROWS = 15;
// By the grid's root: a grid found again from its element is a new wrapper of the same root.
const TOOLTIPS = new WeakMap<HTMLElement, GridTooltips>();

export interface GridTooltips {
  /** Header tooltips by column name. */
  headers?: {[column: string]: string};
  /** Cell tooltips by column name, one line per `\n`. `null`: the cell has no tooltip at all; '' or no function: the
   * cell's text. */
  cells?: {[column: string]: (row: number) => string | null};
}

/** The platform grid over [df] for a pane or a section: read-only, no row header, as tall as its rows up to
 * [maxRows] (then it scrolls), as wide as its host, with [tooltips] on the headers and the cells. */
export function readOnlyGrid(df: DG.DataFrame, tooltips: GridTooltips = {},
  maxRows = MAX_VISIBLE_ROWS): DG.Grid {
  const grid = DG.Viewer.grid(df);
  grid.props.allowEdit = false;
  grid.props.showRowHeader = false;
  grid.props.showCurrentRowIndicator = false;
  grid.props.showAddNewRowIcon = false;
  // Inline, computed: the header (one row until it is drawn, then as tall as its wrapped names) and the rows, so a
  // short table shows no scroll bar.
  grid.root.style.width = '100%';
  const rows = Math.min(df.rowCount, maxRows) * grid.props.rowHeight + 2;
  const fit = (header: number) => {
    const height = `${Math.max(header, grid.props.rowHeight) + rows}px`;
    if (grid.root.style.height !== height)
      grid.root.style.height = height;
  };
  fit(grid.props.rowHeight);
  grid.sub(grid.onAfterDrawContent.subscribe(() => fit(grid.colHeaderHeight)));
  bindTooltips(grid, tooltips);
  return grid;
}

/** Shows [tooltips] on the headers and the cells of [grid]; a cell whose function returns `null` shows none. */
function bindTooltips(grid: DG.Grid, tooltips: GridTooltips): void {
  TOOLTIPS.set(grid.root, tooltips);
  grid.onCellTooltip((cell, x, y) => {
    const text = cellTooltip(grid, cell);
    if (text === null) {
      ui.tooltip.hide();
      return true;
    }
    if (text === '')
      return false;
    ui.tooltip.show(ui.divV(text.split('\n').map((line) => ui.divText(line))), x, y);
    return true;
  });
}

/** The tooltip of [cell] in a {@link readOnlyGrid}; '' for none. */
export function gridTooltip(grid: DG.Grid, cell: DG.GridCell): string {
  return cellTooltip(grid, cell) ?? '';
}

/** The tooltip of [cell]; '' for the platform's, null for none at all. */
function cellTooltip(grid: DG.Grid, cell: DG.GridCell): string | null {
  const tooltips = TOOLTIPS.get(grid.root);
  const column = cell.tableColumn;
  if (tooltips === undefined || column === null)
    return '';
  if (cell.isColHeader)
    return tooltips.headers?.[column.name] ?? '';
  const row = cell.tableRowIndex;
  if (!cell.isTableCell || row === null)
    return '';
  const text = tooltips.cells?.[column.name]?.(row);
  if (text === null)
    return null;
  return text || (column.isNone(row) ? '' : column.getString(row));
}

/** A float column of [values] shown with three decimals. */
export function numberColumn(name: string, values: (number | null)[]): DG.Column {
  const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
  col.meta.format = '0.000';
  return col;
}
