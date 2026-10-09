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

export function showEventDetails(table: DG.DataFrame): void {
  const rowIdx = table.currentRowIdx;
  if (rowIdx < 0)
    return;
  const eventId = table.getCol('id').get(rowIdx);
  if (!eventId)
    return;
  const requestId: string | null = table.col('request_id')?.get(rowIdx) ?? null;
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

  grok.shell.o = requestId ? ui.divV([ui.button('Timeline', () => openTimeline('action', requestId.split('.')[0]),
    'What happened in the user action of this event'), accordion.root]) : accordion.root;
}

/** Opens what happened in one action or request, oldest first, as a table view. */
export async function openTimeline(key: 'action' | 'request', id: string): Promise<void> {
  const progress = DG.TaskBarProgressIndicator.create('Loading the timeline...');
  try {
    const rows = await grok.dapi.log.getTimeline({[key]: id});
    if (!rows.length) {
      grok.shell.info(`Nothing recorded for ${key} ${id}`);
      return;
    }
    const t = DG.DataFrame.fromObjects(rows)!;
    t.name = `Timeline ${id}`;
    grok.shell.addTableView(t);
  }
  catch (e: any) {
    grok.shell.error(`Timeline: ${e?.message ?? e}`);
  }
  finally {
    progress.close();
  }
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
