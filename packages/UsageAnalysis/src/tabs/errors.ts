import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {UaView} from './ua';
import {TimelineView} from './timeline';
import {UaToolbox} from '../ua-toolbox';
import {emptyState, formatGridTimes, formatTime, onRowContextMenu, problemLine, scrollToStartOnFirstDraw,
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
const DEFAULT_FOLDER = 'System:AppData/Ops/errors/';
const DAYS: {[day: string]: string} = {SUN: '0', MON: '1', TUE: '2', WED: '3', THU: '4', FRI: '5', SAT: '6', DAILY: '*',
  WEEKDAYS: '1-5'};
const DAY_NAMES: {[day: string]: string} = {'0': 'Sundays', '1': 'Mondays', '2': 'Tuesdays', '3': 'Wednesdays',
  '4': 'Thursdays', '5': 'Fridays', '6': 'Saturdays', '*': 'Daily', '1-5': 'Weekdays'};
const USERS_SHOWN = 10;

type Spec = {[name: string]: string | number | undefined};

/** The platform's errors as data (`ErrorStats`, the `GET /errors` query): occurrences, or figures by up to three
 * dimensions, a drill-down per row, exports and scheduled exports (`ErrorsSaveJob`). Needs ViewTelemetry. */
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
  saveButton = ui.button('Save as job...', () => this.saveJobDialog());
  saveProblemLine = problemLine();
  host: HTMLDivElement = ui.box();
  table?: DG.DataFrame;
  shownGrid?: DG.Grid;
  shownSpec?: Spec;
  private runs = 0;

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
    ui.tooltip.bind(this.saveButton, 'Save this view as an export job, optionally scheduled');
    this.root.append(ui.divV([ui.divH([exportButton, this.saveButton, this.saveProblemLine], 'ua-toolbar'), this.host],
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
    showProblem(this.saveButton, this.saveProblemLine, this.saveProblem());
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
    const spec: Spec = this.since.value === RANGE ?
      {from: this.from.value!.startOf('day').toISOString(), to: this.to.value!.endOf('day').toISOString()} :
      {since: SINCE[this.since.value!]};
    for (const f of TEXT_FILTERS) {
      const value = this.text[f].value?.trim();
      if (value)
        spec[f] = value;
    }
    if (this.group.value?.trim())
      spec.group = this.group.value.trim();
    if (this.service.value)
      spec.service = this.service.value;
    if (this.minUsers.value)
      spec.minUsers = this.minUsers.value;
    if (this.minCount.value)
      spec.minCount = this.minCount.value;
    if (this.regressed.value)
      spec.regressed = 'true';
    if (by.length) {
      spec.by = by.join(',');
      spec.trend = this.trend.value!;
    }
    return spec;
  }

  load(): void {
    if (this.problem() != null)
      return;
    const spec = this.spec();
    const run = ++this.runs;
    this.table = undefined;
    this.shownGrid = undefined;
    this.shownSpec = undefined;
    this.refresh();
    grok.shell.o = null;
    ui.empty(this.host);
    this.host.append(ui.waitBox(async () => {
      try {
        const t: DG.DataFrame = await grok.functions.call('ErrorStats', {spec: JSON.stringify(spec)});
        if (run !== this.runs)
          return ui.div();
        t.name = 'Errors';
        this.table = t;
        this.shownSpec = spec;
        if (t.rowCount === 0) {
          this.refresh();
          return emptyState('No errors match', ErrorsView.emptyHint(spec));
        }
        this.shownGrid = this.grid(t, spec);
        this.refresh();
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
          menu.item('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'request', request));
      });
    }
    if (grid.col('signature'))
      grid.col('signature')!.width = 70;
    formatGridTimes(grid);
    for (const [name, header] of Object.entries(HEADERS)) {
      if (grid.col(name))
        grid.col(name)!.name = header;
    }
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
    const set: string[] = [];
    const by = spec.by ? (spec.by as string).split(',') : [];
    if (spec.group)
      set.push(`Group is ${spec.group}`);
    if (spec.minUsers) {
      set.push(`Min users is ${spec.minUsers}` + (by.includes('user') && Number(spec.minUsers) > 1 ?
        ' — grouping by user leaves one user per row' : ''));
    }
    if (spec.minCount)
      set.push(`Min count is ${spec.minCount}`);
    if (spec.service)
      set.push(`Service is ${spec.service}`);
    for (const f of TEXT_FILTERS) {
      if (spec[f])
        set.push(`${f[0].toUpperCase()}${f.substring(1)} is "${spec[f]}"`);
    }
    if (spec.regressed)
      set.push('Regressed is on');
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
    const occurrences: Promise<DG.DataFrame> = grok.functions.call('ErrorStats', {spec: JSON.stringify(drill)});
    const error: string = t.col('topError') ? t.get('topError', i) : t.get('error', i);
    const usersSpec: Spec = {...spec, ...filters, by: 'user', trend: undefined, minUsers: undefined,
      minCount: undefined};
    const acc = DG.Accordion.create();
    acc.addPane(ErrorsView.rowTitle(t, spec, i), () => ui.divV([
      spec.by ? ui.wait(async () => {
        const u: DG.DataFrame = await grok.functions.call('ErrorStats', {spec: JSON.stringify(usersSpec)});
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
        request ? ui.link('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'request', request)) : null,
      ]));
    }
    acc.addPane('Occurrences', () => ui.wait(async () => {
      const o = await occurrences;
      const rows = [...Array(Math.min(o.rowCount, SHOWN)).keys()];
      const table = ui.table(rows, (r) => [formatTime(o.get('time', r)), o.get('user', r) ?? '',
        ErrorsView.shortSignature(o.get('signature', r)), o.get('error', r),
        o.get('requestId', r) ? ui.link('timeline', () =>
          TimelineView.open(this.uaToolbox.viewHandler, 'request', o.get('requestId', r))) : ''],
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
      if (s.rowCount === 0)
        return ui.divText('No sessions');
      return ui.table([...Array(s.rowCount).keys()], (r) => [s.get('user', r), formatTime(s.get('first', r)),
        s.get('count', r), ui.link('timeline', () =>
          TimelineView.open(this.uaToolbox.viewHandler, 'session', s.get('session', r)))],
      ['user', 'first', 'errors', '']);
    }));
    acc.addPane('Reports', () => ui.wait(async () => {
      const signatures = ErrorsView.values(await occurrences, 'signature');
      const r = signatures.length ? await queries.errorReports(signatures) : DG.DataFrame.create();
      if (r.rowCount === 0)
        return ui.divText('No reports');
      return ui.table([...Array(r.rowCount).keys()], (k) => [
        ui.link(`#${r.get('number', k)}`, async () => grok.shell.addView(await funcs.reportsApp(`/${r.get('number', k)}`))),
        formatTime(r.get('created_on', k)), r.get('is_auto', k) ? 'auto' : r.get('reporter', k) ?? '',
        r.get('is_resolved', k) ? 'resolved' : 'open', r.get('description', k) ?? ''],
      ['report', 'created', 'by', 'status', 'description']);
    }));
    acc.addPane('Alerts', () => ui.wait(async () => {
      const signatures = ErrorsView.values(await occurrences, 'signature');
      const a = signatures.length ? await queries.errorAlerts(signatures) : DG.DataFrame.create();
      if (a.rowCount === 0)
        return ui.divText('No alerts');
      return ui.table([...Array(a.rowCount).keys()], (k) => [a.get('kind', k), a.get('status', k), a.get('problem', k) ?? '',
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
          const bytes: Uint8Array = await grok.functions.call('Arrow:toParquet',
            {table: ErrorsView.exportTable(this.shownGrid!)});
          DG.Utils.download(file('parquet'), bytes as Uint8Array<ArrayBuffer>);
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

  saveProblem(): string | null {
    if (!this.shownSpec)
      return 'Nothing to save yet';
    return this.shownSpec.since ? null : 'A saved job runs over Since, not From - To';
  }

  /** `since 7d · by signature · trend day · group Chemists`: what a job of [spec] exports. */
  static describe(spec: Spec): string {
    return Object.entries(spec).filter(([_, v]) => v != null && v !== '').map(([k, v]) => `${k} ${v}`).join(' · ');
  }

  saveJobDialog(): void {
    if (this.saveProblem())
      return;
    const spec = {...this.shownSpec};
    const saves = ui.input.string('Saves', {value: ErrorsView.describe(spec)});
    saves.readOnly = true;
    const name = ui.input.string('Name');
    const format = ui.input.choice('Format', {value: 'csv', items: ['csv', 'json'], nullable: false,
      tooltipText: 'The server writes CSV or JSON. For Parquet, use Export > Parquet in the browser'});
    const path = ui.input.string('Path', {value: DEFAULT_FOLDER,
      tooltipText: '<connection>/<path>; a path ending in / gets <name>-{date}.<format>; {date} is the run\'s UTC date'});
    const schedule = ui.input.string('Schedule, UTC', {placeholder: 'MON 07:00',
      tooltipText: 'Empty for no schedule, "MON 07:00", "DAILY 07:00", "WEEKDAYS 07:00" or a five-field cron, in UTC'});
    const dialogProblem = (): string | null => {
      if (!name.value?.trim())
        return 'Enter the name';
      if (!path.value?.trim() || path.value.trim().indexOf('/') <= 0)
        return 'Enter the path as <connection>/<path>';
      return schedule.value?.trim() && ErrorsView.cron(schedule.value) == null ?
        'The schedule is "MON 07:00", "DAILY 07:00", "WEEKDAYS 07:00" or a five-field cron' : null;
    };
    const line = problemLine();
    const dialog = ui.dialog('Save as job');
    dialog.add(saves);
    for (const input of [name, format, path, schedule]) {
      dialog.add(input);
      input.onChanged.subscribe(() => showProblem(dialog.getButton('OK'), line, dialogProblem()));
    }
    dialog.add(line);
    dialog.onOK(async () => {
      const folder = path.value.trim();
      const target = folder.endsWith('/') ?
        `${folder}${name.value.trim().toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '')}-{date}.${format.value}` :
        folder;
      const cron = schedule.value?.trim() ? ErrorsView.cron(schedule.value)! : '';
      try {
        const job = JSON.parse(await grok.functions.call('ErrorsSaveJob',
          {name: name.value.trim(), spec: JSON.stringify(spec), format: format.value, path: target, cron}));
        grok.shell.info(`Saved job "${job.name}": ${ErrorsView.scheduleText(job.cron)} → ${job.path}`);
      }
      catch (e: any) {
        grok.shell.error(`Save as job: ${e?.message ?? e}`);
      }
    });
    dialog.show();
    showProblem(dialog.getButton('OK'), line, dialogProblem());
    name.input.focus();
  }

  /** [cron] in words, `Mondays 08:00 UTC`; `no schedule` for none. */
  static scheduleText(cron: string | null | undefined): string {
    if (!cron?.trim())
      return 'no schedule';
    const f = cron.trim().split(/\s+/);
    const pad = (v: string) => v.padStart(2, '0');
    return f.length === 5 && /^\d+$/.test(f[0]) && /^\d+$/.test(f[1]) && f[2] === '*' && f[3] === '*' &&
      DAY_NAMES[f[4]] ? `${DAY_NAMES[f[4]]} ${pad(f[1])}:${pad(f[0])} UTC` : `cron "${cron.trim()}" UTC`;
  }

  /** `MON 07:00`, `DAILY 07:00`, `WEEKDAYS 07:00` or a five-field cron → cron; null when it is neither. */
  static cron(schedule: string): string | null {
    const s = schedule.trim();
    const m = /^([A-Za-z]+)\s+(\d{1,2}):(\d{2})$/.exec(s);
    if (m) {
      const day = DAYS[m[1].toUpperCase()];
      return day === undefined || +m[2] > 23 || +m[3] > 59 ? null : `${+m[3]} ${+m[2]} * * ${day}`;
    }
    return s.split(/\s+/).length === 5 ? s : null;
  }
}
