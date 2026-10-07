import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {take} from 'rxjs/operators';
import dayjs from 'dayjs';

export const colors = {'passed': '#3CB173', 'failed': '#EB6767', 'skipped': '#FFA24A'};

export function getTime(date: Date, format: string = 'en-GB'): string {
  return date.toLocaleString(format, {hour12: false, timeZone: 'GMT'}).replace(',', '');
}

export function getDate(date: Date): string {
  return date.toLocaleDateString('en-US', {month: '2-digit', day: '2-digit', year: 'numeric'});
}

let usersCache: { [name: string]: DG.User } | null = null;

export async function loadUsers(): Promise<{ [name: string]: DG.User }> {
  if (!usersCache) {
    usersCache = {};
    for (const user of await grok.dapi.users.list())
      usersCache[user.friendlyName] = user;
  }
  return usersCache;
}

export function showEventDetails(table: DG.DataFrame): void {
  const rowIdx = table.currentRowIdx;
  if (rowIdx < 0)
    return;
  const eventId = table.getCol('id').get(rowIdx);
  if (!eventId)
    return;
  const accordion = DG.Accordion.create();
  accordion.addPane('Details', () => ui.wait(async () => {
    const t: DG.DataFrame = await grok.functions.call('UsageAnalysis:LogEventParameters', {eventId});
    if (t.rowCount === 0)
      return ui.divText('No details available');
    const names = t.getCol('param_name').toList();
    const values = t.getCol('value').toList();
    const map: {[key: string]: string} = {};
    for (let i = 0; i < names.length; i++)
      map[names[i]] = values[i];
    return ui.tableFromMap(map);
  }), true);

  grok.shell.o = accordion.root;
}

export function setupUserIconRenderer(grid: DG.Grid, users: { [name: string]: DG.User }, columnNames: string[]): void {
  for (const name of columnNames) {
    const col = grid.col(name);
    if (col) {
      col.width = 25;
      col.cellType = 'html';
    }
  }
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell || !columnNames.includes(gc.gridColumn.name) || !gc.cell.value) return;
    const user = users[gc.cell.value];
    if (!user) return;
    const icon = DG.ObjectHandler.forEntity(user)?.renderIcon(user.dart);
    if (icon) {
      icon.style.top = 'calc(50% - 8px)';
      icon.style.left = 'calc(50% - 8px)';
      gc.style.element = ui.tooltip.bind(icon, () => {
        return DG.ObjectHandler.forEntity(user)?.renderTooltip(user.dart)!;
      });
    }
  });
}

/** A grid wider than its view opens scrolled right past its first column: scroll back on the first draw. */
export function scrollToStartOnFirstDraw(grid: DG.Grid): void {
  grid.onAfterDrawContent.pipe(take(1)).subscribe(() => grid.horzScroll.scrollTo(0));
}

/** Calls [handler] with the table row a context menu was opened on (a right-click does not move the current row). */
export function onRowContextMenu(grid: DG.Grid, handler: (menu: DG.Menu, row: number) => void): void {
  let row = -1;
  grid.root.addEventListener('mousedown', (e) => {
    if (e.button !== 2)
      return;
    const r = grid.root.getBoundingClientRect();
    row = grid.hitTest(e.clientX - r.left, e.clientY - r.top)?.tableRowIndex ?? -1;
  }, true);
  grid.onContextMenu.subscribe((menu) => {
    if (row >= 0)
      handler(menu, row);
  });
}

const GRID_TIME_FORMAT = 'yyyy-MM-dd HH:mm:ss UTC';

/** Shows every date column of [grid] without milliseconds, marked UTC: the grid keeps the platform's UTC. */
export function formatGridTimes(grid: DG.Grid): void {
  for (const col of grid.dataFrame.columns.toList()) {
    if (col.type === DG.TYPE.DATE_TIME && grid.col(col.name))
      grid.col(col.name)!.format = GRID_TIME_FORMAT;
  }
}

/** A date column's value as [formatGridTimes] shows it; empty for none. */
export function formatTime(value: any): string {
  return value == null ? '' : `${dayjs(value).toISOString().replace('T', ' ').substring(0, 19)} UTC`;
}

/** A line that says why an action is disabled: a disabled button shows no tooltip. */
export function problemLine(): HTMLDivElement {
  const line = ui.divText('', 'ua-problem');
  ui.setDisplay(line, false);
  return line;
}

/** Disables [button] while there is a [problem] and shows it in [line]. */
export function showProblem(button: HTMLButtonElement, line: HTMLElement, problem: string | null): void {
  button.disabled = problem != null;
  line.textContent = problem ?? '';
  ui.setDisplay(line, problem != null);
}

/** A centred message with a hint, for a list with nothing to show. */
export function emptyState(message: string, hint: string): HTMLDivElement {
  return ui.divV([ui.divText(message), ui.divText(hint, 'ua-empty-hint')], 'ua-empty');
}

export function rowsTable(t: DG.DataFrame, empty: string, row: (i: number) => any[], headers: string[]): HTMLElement {
  return t.rowCount ? ui.table([...Array(t.rowCount).keys()], row, headers) : ui.divText(empty);
}
