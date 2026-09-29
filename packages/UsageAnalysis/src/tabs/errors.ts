import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {UaView} from './ua';
import {TimelineView} from './timeline';
import {UaToolbox} from '../ua-toolbox';
import {onRowContextMenu, scrollToStartOnFirstDraw} from '../utils';
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
const SHOWN = 20;
const ROUTE = /^([A-Za-z]+\s+)?\//;
const HEADERS: {[column: string]: string} = {count: 'occurrences', firstVersion: 'first seen in', firstSeen: 'first seen',
  lastSeen: 'last seen', newInRange: 'new', mttrMinutes: 'MTTR, min', requestId: 'request'};
const DEFAULT_FOLDER = 'System:AppData/Ops/errors/';
const DAYS: {[day: string]: string} = {SUN: '0', MON: '1', TUE: '2', WED: '3', THU: '4', FRI: '5', SAT: '6', DAILY: '*',
  WEEKDAYS: '1-5'};

type Spec = {[name: string]: string | number | undefined};

/** The platform's errors as data (`ErrorStats`, the `GET /errors` query): occurrences, or figures by up to three
 * dimensions, a drill-down per row, exports and scheduled exports (`ErrorsSaveJob`). Needs ViewTelemetry. */
export class ErrorsView extends UaView {
  since = ui.input.choice('Since', {value: '7 days', items: [...Object.keys(SINCE), RANGE], nullable: false});
  from = ui.input.date('From');
  to = ui.input.date('To');
  group = ui.input.choice<string>('Group', {value: '', items: ['']});
  by = ui.input.multiChoice<string>('Group by', {value: ['signature'], items: DIMENSIONS});
  service = ui.input.choice('Service', {value: '', items: ['', 'server', 'client']});
  text: {[name: string]: DG.InputBase<string>} = Object.fromEntries(TEXT_FILTERS.map((f) =>
    [f, ui.input.string(f[0].toUpperCase() + f.substring(1))]));
  minUsers = ui.input.int('Min users');
  minCount = ui.input.int('Min count');
  regressed = ui.input.bool('Regressed', {tooltipText: 'Only signatures unmuted by a newer version in the window'});
  trend = ui.input.choice('Trend', {value: 'day', items: ['day', 'hour'], nullable: false});
  applyButton = ui.bigButton('Apply', () => this.load());
  host: HTMLDivElement = ui.box();
  table?: DG.DataFrame;
  shownSpec?: Spec;
  private runs = 0;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Errors';
  }

  async initViewers(path?: string): Promise<void> {
    const groups = (await grok.dapi.groups.list()).filter((g) => !g.personal).map((g) => g.friendlyName).sort();
    this.group.items = ['', ...groups];
    this.text.signature.setTooltip('A stack hash or its first 6+ characters');
    this.text.route.setTooltip('"<METHOD> /path" or "/path", without /api');
    this.text.connection.setTooltip('Namespace:Name or the connection id');
    this.text.function.setTooltip('The function\'s nqName');
    const inputs: DG.InputBase[] = [this.since, this.from, this.to, this.group, this.by, this.service,
      ...Object.values(this.text), this.minUsers, this.minCount, this.regressed, this.trend];
    for (const input of inputs)
      input.onChanged.subscribe(() => this.refresh());
    const form = ui.narrowForm(inputs);
    form.append(this.applyButton);
    ui.tooltip.bind(this.applyButton, () => this.problem() ?? 'Show the errors');
    this.uaToolbox.addTabPane(this.name, form);

    const exportButton = ui.button('Export', (e: MouseEvent) => this.exportMenu(e));
    const saveButton = ui.button('Save as job...', () => this.saveJobDialog());
    ui.tooltip.bind(saveButton, () => this.saveProblem() ?? 'Save this view as an export job, optionally scheduled');
    this.root.append(ui.divV([ui.divH([exportButton, saveButton], 'ua-toolbar'), this.host], 'ui-box'));
    this.refresh();
    this.load();
  }

  refresh(): void {
    const range = this.since.value === RANGE;
    ui.setDisplay(this.from.root, range);
    ui.setDisplay(this.to.root, range);
    this.trend.enabled = (this.by.value ?? []).length > 0;
    this.applyButton.disabled = this.problem() != null;
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
    if (this.group.value)
      spec.group = this.group.value;
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
    ui.empty(this.host);
    this.host.append(ui.waitBox(async () => {
      try {
        const t: DG.DataFrame = await grok.functions.call('ErrorStats', {spec: JSON.stringify(spec)});
        if (run !== this.runs)
          return ui.div();
        t.name = 'Errors';
        this.table = t;
        this.shownSpec = spec;
        return t.rowCount === 0 ? ui.divText('No errors match') : this.grid(t, spec).root;
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
      grid.col('topError')!.width = 300;
      grid.col('topError')!.name = bySignature ? 'error' : 'top error';
    }
    else {
      grid.col('error')!.width = 300;
      onRowContextMenu(grid, (menu, i) => {
        const request = t.get('requestId', i);
        if (request)
          menu.item('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'request', request));
      });
    }
    if (grid.col('signature'))
      grid.col('signature')!.width = 70;
    for (const [name, header] of Object.entries(HEADERS)) {
      if (grid.col(name))
        grid.col(name)!.name = header;
    }
    grid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;
      const name = gc.gridColumn.column?.name;
      if (name === 'signature')
        gc.customText = ErrorsView.shortSignature(gc.cell.value);
      else if (name === 'state')
        gc.customText = ErrorsView.stateText(t, gc.cell.rowIndex);
    });
    grid.onCellRender.subscribe((args) => {
      if (!args.cell.isTableCell || args.cell.gridColumn.column?.name !== 'trend')
        return;
      const counts = ErrorsView.buckets(args.cell.cell.value);
      const max = Math.max(0, ...counts);
      const b = args.bounds;
      const w = (b.width - 8) / Math.max(1, counts.length);
      args.g.fillStyle = DG.Color.toHtml(DG.Color.blue);
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

  static stateText(t: DG.DataFrame, i: number): string {
    const state = t.get('state', i);
    if (state !== 'muted')
      return state ?? '';
    const version = t.get('stateVersion', i);
    const until = t.col('stateUntil')!.getString(i);
    return version ? `muted → ${version}` : until ? `muted until ${until}` : 'muted';
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

  showRow(t: DG.DataFrame, spec: Spec, i: number): void {
    const filters = ErrorsView.rowFilters(t, spec, i);
    const drill: Spec = {...spec, ...filters, by: undefined, trend: undefined, limit: DRILL_LIMIT};
    const occurrences: Promise<DG.DataFrame> = grok.functions.call('ErrorStats', {spec: JSON.stringify(drill)});
    const acc = DG.Accordion.create();
    if (!spec.by) {
      const details: {[key: string]: any} = {};
      for (const c of t.columns.toList())
        details[HEADERS[c.name] ?? c.name] = c.getString(i);
      const request = t.get('requestId', i);
      acc.addPane('Occurrence', () => ui.divV([ui.tableFromMap(details),
        request ? ui.link('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'request', request)) : null,
      ]), true);
    }
    acc.addPane('Occurrences', () => ui.wait(async () => {
      const o = await occurrences;
      const rows = [...Array(Math.min(o.rowCount, SHOWN)).keys()];
      const table = ui.table(rows, (r) => [o.col('time')!.getString(r), o.get('user', r) ?? '',
        ErrorsView.shortSignature(o.get('signature', r)), o.get('error', r),
        o.get('requestId', r) ? ui.link('timeline', () =>
          TimelineView.open(this.uaToolbox.viewHandler, 'request', o.get('requestId', r))) : ''],
      ['time', 'user', 'signature', 'error', '']);
      o.name = `Errors of ${Object.values(filters).join(' · ')}`;
      return ui.divV([table, o.rowCount > SHOWN || o.rowCount === 0 ?
        ui.divText(o.rowCount === 0 ? 'No occurrences' : `${o.rowCount}${o.rowCount === DRILL_LIMIT ? '+' : ''} in all`) :
        null, o.rowCount ? ui.link('Open as a table', () => grok.shell.addTableView(o)) : null]);
    }), true);
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
      return ui.table([...Array(s.rowCount).keys()], (r) => [s.get('user', r), s.col('first')!.getString(r),
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
        r.col('created_on')!.getString(k), r.get('is_auto', k) ? 'auto' : r.get('reporter', k) ?? '',
        r.get('is_resolved', k) ? 'resolved' : 'open', r.get('description', k) ?? ''],
      ['report', 'created', 'by', 'status', 'description']);
    }));
    acc.addPane('Alerts', () => ui.wait(async () => {
      const signatures = ErrorsView.values(await occurrences, 'signature');
      const a = signatures.length ? await queries.errorAlerts(signatures) : DG.DataFrame.create();
      if (a.rowCount === 0)
        return ui.divText('No alerts');
      return ui.table([...Array(a.rowCount).keys()], (k) => [a.get('kind', k), a.get('status', k),
        a.col('opened_at')!.getString(k), a.get('summary', k) ?? ''], ['kind', 'status', 'opened', 'summary']);
    }));
    grok.shell.o = acc.root;
  }

  static values(t: DG.DataFrame, name: string): string[] {
    const col = t.col(name);
    return col ? col.categories.filter((v) => v !== '') : [];
  }

  exportMenu(e: MouseEvent): void {
    const noTable = () => this.table ? null : 'Nothing to export yet';
    const arrow = DG.Func.find({package: 'Arrow', name: 'toParquet'}).length > 0;
    DG.Menu.popup()
      .item('CSV', () => DG.Utils.download('errors.csv', this.table!.toCsv(), 'text/csv'), null, {isEnabled: noTable})
      .item('JSON', () => DG.Utils.download('errors.json', JSON.stringify(this.table!.toJson()), 'application/json'),
        null, {isEnabled: noTable})
      .item('Parquet', async () => {
        try {
          const bytes: Uint8Array = await grok.functions.call('Arrow:toParquet', {table: this.table});
          DG.Utils.download('errors.parquet', bytes as Uint8Array<ArrayBuffer>);
        }
        catch (err: any) {
          grok.shell.error(`Parquet: ${err?.message ?? err}`);
        }
      }, null, {isEnabled: () => noTable() ?? (arrow ? null : 'Install the Arrow package to export Parquet')})
      .show({causedBy: e});
  }

  saveProblem(): string | null {
    if (!this.shownSpec)
      return 'Nothing to save yet';
    return this.shownSpec.since ? null : 'A saved job runs over Since, not From - To';
  }

  saveJobDialog(): void {
    const problem = this.saveProblem();
    if (problem) {
      grok.shell.warning(problem);
      return;
    }
    const spec = {...this.shownSpec};
    const name = ui.input.string('Name');
    const format = ui.input.choice('Format', {value: 'csv', items: ['csv', 'json'], nullable: false});
    const path = ui.input.string('Path', {value: DEFAULT_FOLDER,
      tooltipText: '<connection>/<path>; a path ending in / gets <name>-{date}.<format>; {date} is the run\'s UTC date'});
    const schedule = ui.input.string('Schedule', {
      tooltipText: 'Empty for no schedule, "MON 07:00", "DAILY 07:00", "WEEKDAYS 07:00" or a five-field cron (UTC)'});
    const dialogProblem = (): string | null => {
      if (!name.value?.trim())
        return 'Enter the name';
      if (!path.value?.trim() || path.value.trim().indexOf('/') <= 0)
        return 'Enter the path as <connection>/<path>';
      return schedule.value?.trim() && ErrorsView.cron(schedule.value) == null ?
        'The schedule is "MON 07:00", "DAILY 07:00", "WEEKDAYS 07:00" or a five-field cron' : null;
    };
    const dialog = ui.dialog('Save as job');
    for (const input of [name, format, path, schedule]) {
      dialog.add(input);
      input.onChanged.subscribe(() => dialog.getButton('OK').disabled = dialogProblem() != null);
    }
    dialog.onOK(async () => {
      const folder = path.value.trim();
      const target = folder.endsWith('/') ?
        `${folder}${name.value.trim().toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '')}-{date}.${format.value}` :
        folder;
      const cron = schedule.value?.trim() ? ErrorsView.cron(schedule.value)! : '';
      try {
        const job = JSON.parse(await grok.functions.call('ErrorsSaveJob',
          {name: name.value.trim(), spec: JSON.stringify(spec), format: format.value, path: target, cron}));
        grok.shell.info(`Saved job "${job.name}" ${job.cron ?? '(no schedule)'} → ${job.path}`);
      }
      catch (e: any) {
        grok.shell.error(`Save as job: ${e?.message ?? e}`);
      }
    });
    dialog.show();
    ui.tooltip.bind(dialog.getButton('OK'), () => dialogProblem() ?? 'Save the job');
    dialog.getButton('OK').disabled = true;
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
