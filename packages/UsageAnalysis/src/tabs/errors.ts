import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import {emptyState, formatGridTimes, formatTime, onRowContextMenu, problemLine, rowsTable, scrollToStartOnFirstDraw,
  showProblem} from '../utils';
import {funcs, queries} from '../package-api';

import '../../css/usage_analysis.css';

export const DIMENSIONS = ['signature', 'package', 'version', 'user', 'group', 'service', 'route', 'server',
  'connection', 'function'];
const TEXT_FILTERS = ['signature', 'package', 'version', 'user', 'route', 'server', 'connection', 'function'];
const SINCE: {[label: string]: string} = {'1 hour': '1h', '24 hours': '24h', '7 days': '7d', '30 days': '30d',
  '90 days': '90d'};
const RANGE = 'From - To';
const MAX_TREND = 60;
const DRILL_LIMIT = 500;
const SHOWN = 5;
const ROUTE = /^([A-Za-z]+\s+)?\//;
const FORMULA = /^[=+\-@\t\r]/;
const HEADERS: {[column: string]: string} = {count: 'occurrences', firstVersion: 'first seen in',
  firstSeen: 'first seen', lastSeen: 'last seen', newInRange: 'new in range', mttrMinutes: 'MTTR, min',
  requestId: 'request id', state: 'alert'};
const ERROR_COLUMNS = ['error', 'topError'];
const TIME_COLUMNS = ['time', 'firstSeen', 'lastSeen', 'stateUntil'];
const INT_COLUMNS = ['count', 'users', 'sessions', 'signatures', 'newInRange', 'incidents', 'ms'];
const USERS_SHOWN = 10;
const NOT_HINTED = ['since', 'from', 'to', 'by', 'trend'];

type Spec = {[name: string]: string | number | boolean | undefined};

/** The platform's errors as data (`grok.dapi.log.getErrors`, the `GET /errors` query): occurrences, or figures
 * by up to three dimensions, a drill-down per row and exports. Needs ViewTelemetry. */
export class ErrorsView extends UaView {
  since = ui.input.choice('Since', {value: '7 days', items: [...Object.keys(SINCE), RANGE], nullable: false});
  from = ui.input.date('From');
  to = ui.input.date('To');
  group!: DG.InputBase<string>;
  by = ui.input.multiChoice<string>('Group by', {value: ['signature'], items: DIMENSIONS});
  service = ui.input.choice('Service', {value: '', items: ['', 'server', 'client']});
  text: {[name: string]: DG.InputBase<string>} = Object.fromEntries(TEXT_FILTERS.map((f) =>
    [f, ui.input.string(f[0].toUpperCase() + f.substring(1))]));
  minUsers = ui.input.int('Min users');
  minCount = ui.input.int('Min count');
  regressed = ui.input.bool('Regressed', {tooltipText: 'Only signatures that came back in the window after a fix, or on the version their mute waited for'});
  trend = ui.input.choice('Trend', {value: 'day', items: ['day', 'hour'], nullable: false});
  applyButton = ui.bigButton('Apply', () => this.load());
  applyProblem = problemLine();
  host: HTMLDivElement = ui.box();
  table?: DG.DataFrame;
  shownGrid?: DG.Grid;
  private runs = 0;
  /** The app's `?error=<stack hash>` parameter, which alert links carry (`/apps/usage/errors?error=...`): the
   * platform passes it to the app function, not in the URL. */
  static urlError?: string;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Errors';
  }

  async initViewers(path?: string): Promise<void> {
    const allUsers = DG.Group.defaultGroupsIds['All users'];
    const groups = (await grok.dapi.groups.list()).filter((g) => !g.personal && g.id !== allUsers)
      .map((g) => g.friendlyName).sort();
    this.group = ui.typeAhead('Group', {source: {local: groups}, minLength: 1, limit: 30, highlight: true});
    this.group.setTooltip('Only errors of this group\'s members; empty for everyone. A group filter leaves out ' +
      'the errors without a user');
    this.text.signature.setTooltip('A stack hash or its first 6+ characters');
    this.text.route.setTooltip('"<METHOD> /path" or "/path", without /api');
    this.text.connection.setTooltip('Namespace:Name or the connection id');
    this.text.function.setTooltip('The function\'s nqName');
    if (ErrorsView.urlError)
      this.text.signature.value = ErrorsView.urlError;
    const main: DG.InputBase[] = [this.since, this.from, this.to, this.group, this.by];
    const more: DG.InputBase[] = [this.service, ...Object.values(this.text), this.minUsers, this.minCount,
      this.regressed, this.trend];
    for (const input of [...main, ...more])
      input.onChanged.subscribe(() => this.refresh());
    const form = ui.narrowForm(main);
    ui.tooltip.bind(this.applyButton, 'Show the errors');
    const moreFilters = ui.accordion();
    moreFilters.addPane('More filters', () => ui.narrowForm(more), false);
    form.append(this.applyButton, this.applyProblem, moreFilters.root);
    this.uaToolbox.addTabPane(this.name, form);

    const exportButton = ui.button('Export', (e: MouseEvent) => this.exportMenu(e));
    this.root.append(ui.divV([ui.divH([exportButton], 'ua-toolbar'), this.host],
      'ui-box'));
    this.refresh();
    this.load();
  }

  refresh(): void {
    const range = this.since.value === RANGE;
    ui.setDisplay(this.from.root, range);
    ui.setDisplay(this.to.root, range);
    this.trend.enabled = (this.by.value ?? []).length > 0;
    showProblem(this.applyButton, this.applyProblem, this.problem());
  }

  problem(): string | null {
    if ((this.by.value ?? []).length > 3)
      return 'Group by takes up to three dimensions';
    if (this.since.value === RANGE) {
      if (!this.from.value || !this.to.value)
        return 'Choose From and To';
      if (!this.from.value.isBefore(this.to.value))
        return 'From must be before To';
    }
    return null;
  }

  spec(): Spec {
    const by = DIMENSIONS.filter((d) => (this.by.value ?? []).includes(d));
    return {
      ...(this.since.value === RANGE ?
        {from: this.from.value!.startOf('day').toISOString(), to: this.to.value!.endOf('day').toISOString()} :
        {since: SINCE[this.since.value!]}),
      ...Object.fromEntries(TEXT_FILTERS.map((f) => [f, this.text[f].value?.trim() || undefined])),
      group: this.group.value?.trim() || undefined,
      service: this.service.value || undefined,
      minUsers: this.minUsers.value || undefined,
      minCount: this.minCount.value || undefined,
      regressed: this.regressed.value || undefined,
      by: by.join(',') || undefined,
      trend: by.length ? this.trend.value! : undefined,
    };
  }

  load(): void {
    if (this.problem() != null)
      return;
    const spec = this.spec();
    const run = ++this.runs;
    this.table = undefined;
    this.shownGrid = undefined;
    this.refresh();
    grok.shell.o = null;
    ui.empty(this.host);
    this.host.append(ui.waitBox(async () => {
      try {
        const t = ErrorsView.frame(await grok.dapi.log.getErrors(spec));
        if (run !== this.runs)
          return ui.div();
        t.name = 'Errors';
        this.table = t;
        if (t.rowCount === 0) {
          this.refresh();
          return emptyState('No errors match', ErrorsView.emptyHint(spec));
        }
        this.shownGrid = this.grid(t, spec);
        this.refresh();
        if (ErrorsView.urlError) {
          ErrorsView.urlError = undefined;
          t.currentRowIdx = 0;
        }
        return this.shownGrid.root;
      }
      catch (e: any) {
        return ui.divText(`Errors: ${e?.message ?? e}`, 'd4-viewer-error');
      }
    }));
  }

  grid(t: DG.DataFrame, spec: Spec): DG.Grid {
    const by = spec.by ? (spec.by as string).split(',') : [];
    const bySignature = by.includes('signature');
    const grid = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
    if (by.length) {
      const order = [...by, ...(bySignature ? ['topError'] : []), 'count', 'users', 'firstVersion', 'firstSeen',
        'lastSeen', 'trend', ...(bySignature ? ['state'] : ['signatures', 'topError']), 'newInRange', 'sessions',
        'incidents', 'mttrMinutes', 'regressed'];
      grid.columns.setVisible(order);
      grid.columns.setOrder(order);
      grid.col('trend')!.width = 120;
      grid.col('count')!.width = 80;
      grid.col('topError')!.width = 400;
      grid.col('topError')!.name = bySignature ? 'error' : 'top error';
    }
    else {
      grid.col('error')!.width = 400;
      onRowContextMenu(grid, (menu, i) => {
        const request = t.get('requestId', i);
        if (request)
          menu.item('Timeline', () => this.openTimeline('request', request));
      });
    }
    if (grid.col('signature'))
      grid.col('signature')!.width = 70;
    formatGridTimes(grid);
    for (const [name, header] of Object.entries(HEADERS))
      if (grid.col(name))
        grid.col(name)!.name = header;
    grid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;
      const name = gc.gridColumn.column?.name ?? '';
      if (ERROR_COLUMNS.includes(name))
        gc.style.textWrap = 'none';
      else if (name === 'signature')
        gc.customText = ErrorsView.shortSignature(gc.cell.value);
      else if (name === 'state')
        gc.customText = ErrorsView.stateText(t, gc.cell.rowIndex);
    });
    grid.onCellTooltip((gc, x, y) => {
      if (!gc.isTableCell || !ERROR_COLUMNS.includes(gc.gridColumn.column?.name ?? '') || !gc.cell.value)
        return false;
      ui.tooltip.show(ui.divText(gc.cell.value, 'ua-error-text'), x, y);
      return true;
    });
    grid.onCellRender.subscribe((args) => {
      if (!args.cell.isTableCell || args.cell.gridColumn.column?.name !== 'trend')
        return;
      const counts = ErrorsView.buckets(args.cell.cell.value);
      const max = Math.max(0, ...counts);
      const b = args.bounds;
      const w = (b.width - 8) / Math.max(1, counts.length);
      args.g.fillStyle = DG.Color.toHtml(DG.Color.getCategoricalColor(0));
      for (let i = 0; i < counts.length; i++) {
        const h = max ? Math.max(1, counts[i] / max * (b.height - 8)) : 1;
        args.g.fillRect(b.x + 4 + i * w, b.y + b.height - 4 - h, Math.max(1, w - 1), h);
      }
      args.preventDefault();
    });
    t.onCurrentRowChanged.subscribe(() => {
      if (t.currentRowIdx >= 0)
        this.showRow(t, spec, t.currentRowIdx);
    });
    scrollToStartOnFirstDraw(grid);
    return grid;
  }

  /** What to change when [spec] finds no errors: the filters it sets beyond Since and Group by. */
  static emptyHint(spec: Spec): string {
    const set = Object.entries(spec).filter(([k, v]) => v != null && !NOT_HINTED.includes(k)).map(([k, v]) =>
      `${k[0].toUpperCase()}${k.substring(1).replace(/[A-Z]/g, (c) => ` ${c.toLowerCase()}`)} is ${v}`);
    return set.length ? `${set.join('; ')}. Clear ${set.length > 1 ? 'them' : 'it'}, or choose a longer Since` :
      'Choose a longer Since';
  }

  static shortSignature(value: string | null): string {
    return value ? value.replace(/-/g, '').substring(0, 6) : '';
  }

  /** The space-separated bucket counts of a trend, summed down to at most [MAX_TREND] buckets. */
  static buckets(value: string | null): number[] {
    const counts = (value ?? '').split(' ').filter((s) => s).map(Number);
    const step = Math.ceil(counts.length / MAX_TREND);
    const buckets: number[] = [];
    for (let i = 0; i < counts.length; i += step)
      buckets.push(counts.slice(i, i + step).reduce((a, b) => a + b, 0));
    return buckets;
  }

  /** The problem's status when it is not active, else its latest alert's status. */
  static stateText(t: DG.DataFrame, i: number): string {
    const state = t.get('state', i);
    if (state === 'not-a-problem')
      return 'not a problem';
    if (state !== 'muted')
      return state ?? '';
    const version = t.get('stateVersion', i);
    const until = t.get('stateUntil', i);
    return version ? `muted → ${version}` : until ? `muted until ${formatTime(until)}` : 'muted';
  }

  /** The filters that narrow [spec] to row [i]: its dimension values, or an occurrence's signature. */
  static rowFilters(t: DG.DataFrame, spec: Spec, i: number): Spec {
    const by = spec.by ? (spec.by as string).split(',') : ['signature'];
    const filters: Spec = {};
    for (const d of by) {
      const value = t.get(d, i);
      if (value != null && value !== '' && (d !== 'route' || ROUTE.test(value)))
        filters[d] = String(value);
    }
    return filters;
  }

  /** `a41f9c · package Chem · 142 occurrences, 4 users`: what row [i] of [t] stands for. */
  static rowTitle(t: DG.DataFrame, spec: Spec, i: number): string {
    const filters = ErrorsView.rowFilters(t, spec, i);
    const parts = Object.entries(filters).map(([d, v]) => d === 'signature' ? ErrorsView.shortSignature(v as string) :
      `${d} ${v}`);
    if (spec.by)
      parts.push(`${t.get('count', i)} occurrences, ${t.get('users', i)} users`);
    else
      parts.push(formatTime(t.get('time', i)));
    return parts.join(' · ');
  }

  showRow(t: DG.DataFrame, spec: Spec, i: number): void {
    const filters = ErrorsView.rowFilters(t, spec, i);
    const drill: Spec = {...spec, ...filters, by: undefined, trend: undefined, limit: DRILL_LIMIT};
    const occurrences = grok.dapi.log.getErrors(drill).then((rows) => ErrorsView.frame(rows));
    const error: string = t.col('topError') ? t.get('topError', i) : t.get('error', i);
    const usersSpec: Spec = {...spec, ...filters, by: 'user', trend: undefined, minUsers: undefined,
      minCount: undefined};
    const acc = DG.Accordion.create();
    acc.addPane(ErrorsView.rowTitle(t, spec, i), () => ui.divV([
      spec.by ? ui.wait(async () => {
        const u = ErrorsView.frame(await grok.dapi.log.getErrors(usersSpec));
        return ui.divText(`Users: ${ErrorsView.usersText(u.col('user')?.toList() ?? [])}`);
      }) : null,
      ui.divText(error ?? '', 'ua-error-text'),
      ui.wait(async () => {
        const o = await occurrences;
        const signature = o.rowCount ? o.get('signature', 0) : null;
        const sample = signature ? await queries.errorSample(signature) : null;
        if (!sample?.rowCount || !sample.get('stack', 0))
          return ui.divText('No stack trace');
        const pre = ui.element('pre', 'ua-stack');
        pre.textContent = sample.get('stack', 0);
        return ui.divV([ui.divText(`Stack trace of the latest occurrence, ${formatTime(sample.get('time', 0))}`), pre]);
      }),
    ]), true);
    if (!spec.by) {
      const details: {[key: string]: any} = {};
      for (const c of t.columns.toList())
        details[HEADERS[c.name] ?? c.name] = c.type === DG.TYPE.DATE_TIME ? formatTime(c.get(i)) : c.getString(i);
      const request = t.get('requestId', i);
      acc.addPane('Occurrence', () => ui.divV([ui.tableFromMap(details),
        request ? ui.link('Timeline', () => this.openTimeline('request', request)) : null,
      ]));
    }
    acc.addPane('Occurrences', () => ui.wait(async () => {
      const o = await occurrences;
      const rows = [...Array(Math.min(o.rowCount, SHOWN)).keys()];
      const table = ui.table(rows, (r) => [formatTime(o.get('time', r)), o.get('user', r) ?? '',
        ErrorsView.shortSignature(o.get('signature', r)), o.get('error', r),
        o.get('requestId', r) ? ui.link('timeline', () => this.openTimeline('request', o.get('requestId', r))) : ''],
      ['time', 'user', 'signature', 'error', '']);
      o.name = `Errors of ${Object.values(filters).join(' · ')}`;
      return ui.divV([table, o.rowCount > SHOWN || o.rowCount === 0 ?
        ui.divText(o.rowCount === 0 ? 'No occurrences' : `${o.rowCount}${o.rowCount === DRILL_LIMIT ? '+' : ''} in all`) :
        null, o.rowCount ? ui.link('Open as a table', () => grok.shell.addTableView(o)) : null]);
    }));
    acc.addPane('Sessions', () => ui.wait(async () => {
      const o = await occurrences;
      const signatures = ErrorsView.values(o, 'signature');
      const users = ErrorsView.values(o, 'user');
      if (o.rowCount === 0 || users.length === 0)
        return ui.divText('No sessions');
      const times = o.col('time')!.toList().filter((v) => v != null).map((v) => dayjs(v).valueOf());
      const s = await queries.errorSessions(signatures, users, new Date(Math.min(...times)).toISOString(),
        new Date(Math.max(...times) + 1000).toISOString());
      return rowsTable(s, 'No sessions', (r) => [s.get('user', r), formatTime(s.get('first', r)), s.get('count', r),
        ui.link('timeline', () => this.openTimeline('session', s.get('session', r)))], ['user', 'first', 'errors', '']);
    }));
    acc.addPane('Reports', () => ui.wait(async () => {
      const signatures = ErrorsView.values(await occurrences, 'signature');
      const r = signatures.length ? await queries.errorReports(signatures) : DG.DataFrame.create();
      return rowsTable(r, 'No reports', (k) => [
        ui.link(`#${r.get('number', k)}`, async () => grok.shell.addView(await funcs.reportsApp(`/${r.get('number', k)}`))),
        formatTime(r.get('created_on', k)), r.get('is_auto', k) ? 'auto' : r.get('reporter', k) ?? '',
        r.get('is_resolved', k) ? 'resolved' : 'open', r.get('description', k) ?? ''],
      ['report', 'created', 'by', 'status', 'description']);
    }));
    acc.addPane('Alerts', () => ui.wait(async () => {
      const signatures = ErrorsView.values(await occurrences, 'signature');
      const a = signatures.length ? await queries.errorAlerts(signatures) : DG.DataFrame.create();
      return rowsTable(a, 'No alerts', (k) => [a.get('kind', k), a.get('status', k), a.get('problem', k) ?? '',
        formatTime(a.get('opened_at', k)), a.get('cleared_at', k) ? formatTime(a.get('cleared_at', k)) : '', a.get('summary', k) ?? ''],
      ['kind', 'status', 'problem', 'opened', 'cleared', 'summary']);
    }));
    grok.shell.o = acc.root;
  }

  /** Up to [USERS_SHOWN] logins, then `+N more`; `none` for errors without a user. */
  static usersText(users: (string | null)[]): string {
    const logins = users.filter((u) => u);
    if (logins.length === 0)
      return 'none';
    const more = logins.length - USERS_SHOWN;
    return logins.slice(0, USERS_SHOWN).join(', ') + (more > 0 ? ` +${more} more` : '');
  }

  static values(t: DG.DataFrame, name: string): string[] {
    const col = t.col(name);
    return col ? col.categories.filter((v) => v !== '') : [];
  }

  exportMenu(e: MouseEvent): void {
    const noTable = () => this.shownGrid ? null : 'Nothing to export yet';
    const arrow = DG.Func.find({package: 'Arrow', name: 'toParquet'}).length > 0;
    const file = (ext: string) => `errors-${new Date().toISOString().substring(0, 10)}.${ext}`;
    DG.Menu.popup()
      .item('CSV', () => DG.Utils.download(file('csv'), ErrorsView.toCsv(ErrorsView.exportTable(this.shownGrid!)),
        'text/csv'), null, {isEnabled: noTable})
      .item('JSON', () => DG.Utils.download(file('json'),
        JSON.stringify(ErrorsView.exportTable(this.shownGrid!).toJson()), 'application/json'), null, {isEnabled: noTable})
      .item('Parquet', async () => {
        try {
          const bytes: Uint8Array<ArrayBuffer> = await grok.functions.call('Arrow:toParquet',
            {table: ErrorsView.exportTable(this.shownGrid!)});
          DG.Utils.download(file('parquet'), bytes);
        }
        catch (err: any) {
          grok.shell.error(`Parquet: ${err?.message ?? err}`);
        }
      }, null, {isEnabled: () => noTable() ?? (arrow ? null : 'Install the Arrow package to export Parquet')})
      .show({causedBy: e});
  }

  /** The columns [grid] shows, in its order, under its headers, with the texts it shows for the signature and the
   * alert. */
  static exportTable(grid: DG.Grid): DG.DataFrame {
    const t = grid.dataFrame;
    const rows = [...Array(t.rowCount).keys()];
    const columns: DG.Column[] = [];
    for (let i = 0; i < grid.columns.length; i++) {
      const gc = grid.columns.byIndex(i)!;
      const col = gc.column;
      if (!gc.visible || !col)
        continue;
      const copy = col.name === 'signature' ?
        DG.Column.fromStrings(gc.name, rows.map((r) => ErrorsView.shortSignature(col.get(r)))) :
        col.name === 'state' ? DG.Column.fromStrings(gc.name, rows.map((r) => ErrorsView.stateText(t, r))) :
          col.clone();
      copy.name = gc.name;
      columns.push(copy);
    }
    return DG.DataFrame.fromColumns(columns);
  }

  /** [t] as CSV, like the server's: a cell a spreadsheet would run as a formula starts with `'`. */
  static toCsv(t: DG.DataFrame): string {
    const safe = t.clone();
    for (const col of safe.columns.toList()) {
      if (col.type !== DG.TYPE.STRING)
        continue;
      const values: (string | null)[] = col.toList();
      col.init((i) => values[i] && FORMULA.test(values[i]!) ? `'${values[i]}` : values[i]);
    }
    return safe.toCsv();
  }

  /** Rows of `grok.dapi.log.getErrors` or `getTimeline` as a table: times as dates, counts as integers,
   * lists (a trend's bucket counts) as their values separated by spaces. */
  static frame(rows: {[key: string]: any}[]): DG.DataFrame {
    if (rows.length === 0)
      return DG.DataFrame.create();
    return DG.DataFrame.fromColumns([...new Set(rows.flatMap((r) => Object.keys(r)))].map((name) => {
      const values = rows.map((r) => r[name] ?? null);
      if (TIME_COLUMNS.includes(name))
        return DG.Column.dateTime(name, rows.length).init((i) => values[i] == null ? null : dayjs(values[i]));
      if (INT_COLUMNS.includes(name))
        return DG.Column.fromList(DG.TYPE.INT, name, values);
      if (name === 'mttrMinutes')
        return DG.Column.fromList(DG.TYPE.FLOAT, name, values);
      if (name === 'regressed')
        return DG.Column.fromList(DG.TYPE.BOOL, name, values);
      return DG.Column.fromList(DG.TYPE.STRING, name,
        values.map((v) => v == null ? null : Array.isArray(v) ? v.join(' ') : `${v}`));
    }));
  }
}
