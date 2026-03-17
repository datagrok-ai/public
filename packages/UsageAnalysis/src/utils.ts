import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

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
