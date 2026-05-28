import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';

import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';

type CardColor = 'green' | 'orange' | 'red' | 'info';
type QueriesMode = 'slowest' | 'most-called' | 'worst-hit';
type PssMode = 'modern' | 'legacy' | 'unavailable';

const DASHBOARD_LIMIT = 10;
const FULL_VIEW_LIMIT = 100000;
const EMAIL_LIMIT = 100;

const THRESH = {
  cacheHit: {green: 99, orange: 97},
  tableHealthCount: {green: 0, orange: 2},
  idleXactSec: {green: 60, orange: 300},
  queryMeanMs: {green: 50, orange: 1000},
  queryHitPct: {green: 90, orange: 70},
  cellDeadPct: {green: 20, orange: 40},
  latencyMs: {green: 500, orange: 2000},
  routeP95Ms: {green: 500, orange: 2000},
  routeP99Ms: {green: 2000, orange: 5000},
  routeErrPct: {green: 10, orange: 50},
  diskUsedPct: {orange: 75, red: 90},
  diskFreeBytes: {orange: 5 * 1024 * 1024 * 1024, red: 2 * 1024 * 1024 * 1024},
};

interface DiskStat {
  path: string;
  mount: string;
  totalBytes: number;
  usedBytes: number;
  freeBytes: number;
  usedPct: number;
}

const COLOR_GREEN = 0xFF3CB173;
const COLOR_ORANGE = 0xFFFFA24A;
const COLOR_RED = 0xFFEB6767;

const PSS_QUERIES: Record<QueriesMode, {modern: (n: number) => Promise<DG.DataFrame>, legacy: (n: number) => Promise<DG.DataFrame>, fullViewQuery: {modern: string, legacy: string}}> = {
  'slowest': {
    modern: queries.metricsTopSlowestQueries,
    legacy: queries.metricsTopSlowestQueriesPg12,
    fullViewQuery: {modern: 'MetricsTopSlowestQueries', legacy: 'MetricsTopSlowestQueriesPg12'},
  },
  'most-called': {
    modern: queries.metricsTopMostCalledQueries,
    legacy: queries.metricsTopMostCalledQueriesPg12,
    fullViewQuery: {modern: 'MetricsTopMostCalledQueries', legacy: 'MetricsTopMostCalledQueriesPg12'},
  },
  'worst-hit': {
    modern: queries.metricsWorstCacheHitQueries,
    legacy: queries.metricsWorstCacheHitQueriesPg12,
    fullViewQuery: {modern: 'MetricsWorstCacheHitQueries', legacy: 'MetricsWorstCacheHitQueriesPg12'},
  },
};

function formatTime(d: Date): string {
  return d.toTimeString().substring(0, 5);
}

function thresholdColor(v: number, t: {green: number, orange: number}, higherIsBetter: boolean): number {
  const isGreen = higherIsBetter ? v >= t.green : v < t.green;
  const isOrange = higherIsBetter ? v >= t.orange : v < t.orange;
  if (isGreen)
    return COLOR_GREEN;
  if (isOrange)
    return COLOR_ORANGE;
  return COLOR_RED;
}

function formatRange(start: dayjs.Dayjs, end: dayjs.Dayjs): string {
  const fmt = 'YYYY-MM-DD HH:mm';
  return `${start.format(fmt)} → ${end.format(fmt)}`;
}

function formatAgo(d: dayjs.Dayjs | null): string {
  if (!d)
    return '—';
  const ms = Date.now() - d.valueOf();
  const s = Math.floor(ms / 1000);
  if (s < 60)
    return `${s} s ago`;
  const m = Math.floor(s / 60);
  if (m < 60)
    return `${m} m ago`;
  const h = Math.floor(m / 60);
  if (h < 24)
    return `${h} h ago`;
  return `${Math.floor(h / 24)} d ago`;
}

function formatBytes(n: number | null | undefined): string {
  if (n == null) return '—';
  if (n < 1024) return `${n} B`;
  const units = ['KB', 'MB', 'GB', 'TB', 'PB'];
  let v = n / 1024;
  let i = 0;
  while (v >= 1024 && i < units.length - 1) { v /= 1024; i++; }
  return `${v >= 100 ? v.toFixed(0) : v.toFixed(1)} ${units[i]}`;
}

function formatCount(n: number | null | undefined): string {
  if (n == null) return '—';
  if (n < 1000) return `${n}`;
  const units = ['K', 'M', 'B', 'T'];
  let v = n / 1000;
  let i = 0;
  while (v >= 1000 && i < units.length - 1) { v /= 1000; i++; }
  return `${v >= 100 ? v.toFixed(0) : v.toFixed(1)}${units[i]}`;
}

export class MetricsView extends UaView {
  private cards = new Map<string, {root: HTMLElement, value: HTMLElement, sub: HTMLElement}>();
  private asOfLabel!: HTMLElement;
  private queriesMode: QueriesMode = 'slowest';
  private queriesGridHost!: HTMLElement;
  private largestTablesHost!: HTMLElement;
  private tableHealthHost!: HTMLElement;
  private httpRoutesHost!: HTMLElement;
  private refreshing = false;
  private pssMode: PssMode | null = null;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Metrics';
  }

  async initViewers(_path?: string): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.asOfLabel = ui.divText('as of —', 'ua-metrics-asof');
    this.queriesGridHost = ui.div([], 'ua-metrics-grid-host ua-metrics-queries-host');
    this.largestTablesHost = ui.div([], 'ua-metrics-grid-host');
    this.tableHealthHost = ui.div([], 'ua-metrics-grid-host');
    this.httpRoutesHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn');
    const sendBtn = ui.bigButton('Send to Datagrok', () => this.sendToDatagrok());
    sendBtn.prepend(ui.iconFA('envelope'), ' ');
    sendBtn.classList.add('ua-metrics-btn-secondary');

    const header = ui.divH([
      this.asOfLabel,
      ui.div([], {style: {flex: '1'}}),
      refreshBtn,
      sendBtn,
    ], 'ua-metrics-header');

    const cardsRow1 = ui.divH([
      this.makeCard('db', 'Database'),
      this.makeCard('tblHealth', 'Table health'),
      this.makeCard('storage', 'Storage'),
      this.makeCard('connections', 'Connections'),
      this.makeCard('disk', 'Disk free'),
    ], 'ua-metrics-cards-row');

    const cardsRow2 = ui.divH([
      this.makeCard('errors', 'Errors'),
      this.makeCard('latency', 'Latency'),
      this.makeCard('sessions', 'Sessions'),
    ], 'ua-metrics-cards-row');

    const queriesPanel = this.buildQueriesPanel();
    const largestTablesPanel = MetricsView.buildPanel('Largest tables', this.largestTablesHost,
      () => this.openFullView(queries.metricsLargestTables, 'MetricsLargestTables', 'Largest tables'));
    const tableHealthPanel = MetricsView.buildPanel('Table health', this.tableHealthHost,
      () => this.openFullView(queries.metricsTableHealth, 'MetricsTableHealth', 'Table health'));
    const httpRoutesPanel = MetricsView.buildPanel('HTTP routes', this.httpRoutesHost,
      () => this.openFullHttpRoutes());

    const tablesRow = ui.divH([largestTablesPanel, tableHealthPanel], 'ua-metrics-tables-row');

    this.root.append(ui.div([
      header,
      cardsRow1,
      cardsRow2,
      queriesPanel,
      tablesRow,
      httpRoutesPanel,
    ], 'ua-metrics-root'));

    this.uaToolbox.filterStream.subscribe(() => this.refreshWindowCards());

    await this.refresh();
  }

  private makeCard(id: string, title: string): HTMLElement {
    const value = ui.divText('—', 'ua-metrics-card-value');
    const sub = ui.divText('', 'ua-metrics-card-sub');
    const root = ui.div([
      ui.divText(title, 'ua-metrics-card-title'),
      value,
      sub,
    ], 'ua-metrics-card');
    this.cards.set(id, {root, value, sub});
    return root;
  }

  private card(id: string): {root: HTMLElement, value: HTMLElement, sub: HTMLElement} {
    return this.cards.get(id)!;
  }

  private setCard(id: string, value: string, sub: string, color: CardColor): void {
    const card = this.card(id);
    card.value.textContent = value;
    card.sub.textContent = sub;
    card.sub.className = `ua-metrics-card-sub ua-metrics-${color}`;
  }

  private static buildPanel(title: string, host: HTMLElement, onOpenFull?: () => void, extraHeader?: HTMLElement): HTMLElement {
    const headerChildren: HTMLElement[] = [ui.divText(title, 'ua-metrics-panel-title')];
    if (extraHeader)
      headerChildren.push(extraHeader);
    headerChildren.push(ui.div([], {style: {flex: '1'}}));
    if (onOpenFull) {
      const addIcon = ui.iconFA('plus', onOpenFull, 'Add to workspace');
      addIcon.classList.add('ua-metrics-add-icon');
      headerChildren.push(addIcon);
    }
    const header = ui.divH(headerChildren, 'ua-metrics-panel-header');
    return ui.div([header, host], 'ua-metrics-panel');
  }

  private buildQueriesPanel(): HTMLElement {
    const modes: Array<[string, QueriesMode]> = [
      ['slowest', 'slowest'],
      ['most called', 'most-called'],
      ['worst cache hit', 'worst-hit'],
    ];
    const buttons = modes.map(([label, mode]) => {
      const btn = ui.toggleButton(label, () => {
        this.queriesMode = mode;
        this.loadQueries();
      });
      if (mode === this.queriesMode)
        btn.classList.add('d4-current');
      return btn;
    });
    const toggle = ui.toggleButtonGroup(buttons);
    toggle.classList.add('ua-metrics-mode-toggle');
    return MetricsView.buildPanel('Queries', this.queriesGridHost,
      () => {
        const q = PSS_QUERIES[this.queriesMode];
        const name = this.pssMode === 'legacy' ? q.fullViewQuery.legacy : q.fullViewQuery.modern;
        const fn = this.pssMode === 'legacy' ? q.legacy : q.modern;
        this.openFullView(fn, name, 'pg_stat_statements');
      }, toggle);
  }

  private async refresh(): Promise<void> {
    if (this.refreshing)
      return;
    this.refreshing = true;
    try {
      this.asOfLabel.textContent = `as of ${formatTime(new Date())}`;
      await Promise.all([
        this.loadDbStats(),
        this.loadTableHealthSummary(),
        this.loadConnections(),
        this.loadQueries(),
        this.loadLargestTables(),
        this.loadTableHealth(),
        this.refreshWindowCards(),
        this.loadStorage(),
        this.loadDisk(),
      ]);
    } finally {
      this.refreshing = false;
    }
  }

  private async refreshWindowCards(): Promise<void> {
    await Promise.all([this.loadErrors(), this.loadSessions(), this.loadLatency(), this.loadHttpRoutes()]);
  }

  private static async safeCall<T>(fn: () => Promise<T>, label: string): Promise<T | null> {
    try {
      return await fn();
    } catch (e) {
      console.warn(`[Metrics] ${label} failed:`, e);
      return null;
    }
  }

  private async loadDbStats(): Promise<void> {
    const card = this.card('db');
    ui.tooltip.bind(card.root, null);
    const df = await MetricsView.safeCall(() => queries.metricsDbStats(), 'MetricsDbStats');
    if (!df || df.rowCount === 0) {
      this.setCard('db', '—', 'unavailable', 'info');
      return;
    }
    const sizePretty = df.get('db_size_pretty', 0) as string;
    const hit = Number(df.get('cache_hit_pct', 0) ?? 0);
    this.setCard('db', sizePretty ?? '—', `hit ${hit.toFixed(0)}%`,
      hit >= THRESH.cacheHit.green ? 'green' : hit >= THRESH.cacheHit.orange ? 'orange' : 'red');
    ui.tooltip.bind(card.root, () => MetricsView.buildDbStatsTooltip(df, hit, sizePretty));
  }

  private static buildDbStatsTooltip(df: DG.DataFrame, hit: number, sizePretty: string): HTMLElement {
    const statsReset = df.get('stats_reset', 0) as dayjs.Dayjs | null;
    const lines: HTMLElement[] = [];
    lines.push(ui.divText(`${sizePretty ?? '—'} · cache hit ${hit.toFixed(2)}%`));
    lines.push(ui.divText(
      'Cache hit = share of heap-block reads served from shared_buffers instead of disk. ' +
      '>99% is healthy for OLTP; lower means more disk I/O.',
      'ua-metrics-tt-dim'));
    if (statsReset)
      lines.push(ui.divText(`Stats since ${statsReset.format('YYYY-MM-DD HH:mm')} · ${formatAgo(statsReset)}`,
        'ua-metrics-tt-dim'));

    const hasOffenders = df.rowCount > 0 && df.get('offender_table', 0) != null;
    if (hasOffenders) {
      lines.push(ui.div([], {style: {height: '6px'}}));
      lines.push(ui.divText('Top tables by cache misses (hit ratio < 95%):'));
      const rows: HTMLElement[] = [
        ui.divH([
          ui.divText('table', 'ua-metrics-tt-col-prefix ua-metrics-tt-head'),
          ui.divText('hit %', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
          ui.divText('disk reads', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
        ]),
      ];
      for (let i = 0; i < df.rowCount; i++) {
        const name = df.get('offender_table', i);
        if (name == null)
          continue;
        rows.push(ui.divH([
          ui.divText(name as string, 'ua-metrics-tt-col-prefix'),
          ui.divText(`${Number(df.get('offender_hit_pct', i) ?? 0).toFixed(1)}%`, 'ua-metrics-tt-col-num'),
          ui.divText(formatCount(Number(df.get('offender_disk_reads', i) ?? 0)), 'ua-metrics-tt-col-num'),
        ]));
      }
      lines.push(ui.div(rows, 'ua-metrics-tt-table'));
      lines.push(ui.divText('Hit % includes heap + index blocks. Filtered to tables with ≥1K disk reads.',
        'ua-metrics-tt-dim'));
    }
    else {
      lines.push(ui.div([], {style: {height: '6px'}}));
      lines.push(ui.divText('No user tables with notable traffic below 95% hit ratio.',
        'ua-metrics-tt-dim'));
    }
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async loadTableHealthSummary(): Promise<void> {
    const card = this.card('tblHealth');
    ui.tooltip.bind(card.root, null);
    const df = await MetricsView.safeCall(() => queries.metricsTableHealthSummary(DASHBOARD_LIMIT),
      'MetricsTableHealthSummary');
    if (!df) {
      this.setCard('tblHealth', '—', 'unavailable', 'info');
      return;
    }
    const count = df.rowCount === 0 ? 0 : ((df.get('unhealthy_count', 0) as number) ?? 0);
    const maxDead = df.rowCount === 0 ? 0 : ((df.get('max_dead_pct', 0) as number) ?? 0);
    this.setCard('tblHealth', `${count}`, count > 0 ? `dead ${maxDead}%` : 'healthy',
      count <= THRESH.tableHealthCount.green ? 'green' : count <= THRESH.tableHealthCount.orange ? 'orange' : 'red');
    ui.tooltip.bind(card.root, () => MetricsView.buildTableHealthTooltip(count, maxDead, df));
  }

  private static buildTableHealthTooltip(count: number, maxDead: number, df: DG.DataFrame): HTMLElement {
    const lines: HTMLElement[] = [];
    if (count === 0) {
      lines.push(ui.divText('All user tables under the dead-tuple threshold.'));
      lines.push(ui.divText('"Unhealthy" = ≥10K live rows and >40% dead tuples.', 'ua-metrics-tt-dim'));
      return ui.div(lines, 'ua-metrics-tt');
    }
    lines.push(ui.divText(
      `${count} table${count === 1 ? '' : 's'} above 40% dead tuples (max ${maxDead}%).`));
    lines.push(ui.divText('Bloat slows scans and stales planner statistics.', 'ua-metrics-tt-dim'));
    if (df.rowCount > 0) {
      const rows: HTMLElement[] = [
        ui.divH([
          ui.divText('table', 'ua-metrics-tt-col-prefix ua-metrics-tt-head'),
          ui.divText('dead %', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
          ui.divText('last vacuum', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
        ]),
      ];
      for (let i = 0; i < df.rowCount; i++) {
        const lastVac = df.get('last_vacuum', i) as dayjs.Dayjs | null;
        rows.push(ui.divH([
          ui.divText(df.get('table_name', i) as string, 'ua-metrics-tt-col-prefix'),
          ui.divText(`${df.get('dead_pct', i)}%`, 'ua-metrics-tt-col-num'),
          ui.divText(formatAgo(lastVac), 'ua-metrics-tt-col-num'),
        ]));
      }
      lines.push(ui.div(rows, 'ua-metrics-tt-table'));
      if (count > df.rowCount)
        lines.push(ui.divText(`showing ${df.rowCount} of ${count}`, 'ua-metrics-tt-dim'));
    }
    lines.push(ui.div([], {style: {height: '6px'}}));
    lines.push(ui.divText(
      'Run VACUUM ANALYZE; if last vacuum is recent, lower autovacuum_vacuum_scale_factor for these tables.',
      'ua-metrics-tt-dim'));
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async loadConnections(): Promise<void> {
    const card = this.card('connections');
    ui.tooltip.bind(card.root, null);
    const summary = await MetricsView.safeCall(() => queries.metricsConnections(), 'MetricsConnections');
    if (!summary || summary.rowCount === 0) {
      this.setCard('connections', '—', 'unavailable', 'info');
      return;
    }
    const total = (summary.get('total', 0) as number) ?? 0;
    const idleX = (summary.get('idle_in_xact', 0) as number) ?? 0;
    const oldestSec = (summary.get('oldest_idle_xact_seconds', 0) as number) ?? 0;
    this.setCard('connections', `${total}`, idleX > 0 ? `${idleX} idle xact` : 'no idle xact',
      oldestSec >= THRESH.idleXactSec.orange ? 'red' : oldestSec >= THRESH.idleXactSec.green ? 'orange' : 'green');

    let offenders: DG.DataFrame | null = null;
    let pending = true;
    ui.tooltip.bind(card.root, () => MetricsView.buildConnectionsTooltip(summary, offenders, pending));
    setTimeout(async () => {
      offenders = await MetricsView.safeCall(
        () => queries.metricsConnectionsOffenders(DASHBOARD_LIMIT, THRESH.idleXactSec.green, 30),
        'MetricsConnectionsOffenders');
      pending = false;
    }, 500);
  }

  private static buildConnectionsTooltip(
    summary: DG.DataFrame, offenders: DG.DataFrame | null, pending: boolean): HTMLElement {
    const total = (summary.get('total', 0) as number) ?? 0;
    const active = (summary.get('active', 0) as number) ?? 0;
    const idle = (summary.get('idle', 0) as number) ?? 0;
    const idleX = (summary.get('idle_in_xact', 0) as number) ?? 0;
    const waiting = (summary.get('waiting_on_lock', 0) as number) ?? 0;
    const oldestSec = (summary.get('oldest_idle_xact_seconds', 0) as number) ?? 0;

    const lines: HTMLElement[] = [];
    const headline = idleX > 0
      ? `${total} connections · ${idleX} idle in transaction · oldest ${DG.Utils.formatDuration(oldestSec)} ago`
      : `${total} connections`;
    lines.push(ui.divText(headline));
    lines.push(ui.divText(
      `${active} active · ${idle} idle · ${idleX} idle in transaction · ${waiting} waiting on lock`,
      'ua-metrics-tt-dim'));

    if (pending) {
      lines.push(ui.div([], {style: {height: '4px'}}));
      lines.push(ui.divText('Loading session details…', 'ua-metrics-tt-dim'));
      return ui.div(lines, 'ua-metrics-tt');
    }
    if (!offenders) {
      lines.push(ui.div([], {style: {height: '4px'}}));
      lines.push(ui.divText('Session details unavailable.', 'ua-metrics-tt-dim'));
      return ui.div(lines, 'ua-metrics-tt');
    }
    const rowCount = offenders.rowCount;
    if (rowCount === 0) {
      lines.push(ui.div([], {style: {height: '4px'}}));
      lines.push(ui.divText('All sessions look healthy.'));
      return ui.div(lines, 'ua-metrics-tt');
    }

    const states: string[] = [];
    const waitEvtTypes: string[] = [];
    let maxBlocks = 0;
    for (let i = 0; i < rowCount; i++) {
      states.push(offenders.get('state', i) as string ?? '');
      waitEvtTypes.push(offenders.get('wait_event_type', i) as string ?? '');
      maxBlocks = Math.max(maxBlocks, (offenders.get('blocks', i) as number) ?? 0);
    }
    const hasIdleXact = states.some((s) => s === 'idle in transaction');
    const hasActiveLock = states.some((s, i) => s === 'active' && waitEvtTypes[i] === 'Lock');
    const showBlocks = maxBlocks > 0;

    if (hasIdleXact) {
      lines.push(ui.div([], {style: {height: '4px'}}));
      lines.push(ui.divText(
        'Idle in transaction = client opened BEGIN but hasn\'t sent COMMIT/ROLLBACK. ' +
        'These block VACUUM and hold row locks until they close.',
        'ua-metrics-tt-dim'));
    }

    const header: HTMLElement[] = [
      ui.divText('state', 'ua-metrics-tt-col-state ua-metrics-tt-head'),
      ui.divText('age', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
      ui.divText('app', 'ua-metrics-tt-col-app ua-metrics-tt-head'),
      ui.divText('user', 'ua-metrics-tt-col-user ua-metrics-tt-head'),
      ui.divText('query', 'ua-metrics-tt-col-query ua-metrics-tt-head'),
    ];
    if (showBlocks)
      header.push(ui.divText('blocks', 'ua-metrics-tt-col-num ua-metrics-tt-head'));
    const rows: HTMLElement[] = [ui.divH(header)];
    for (let i = 0; i < rowCount; i++) {
      const state = offenders.get('state', i) as string ?? '';
      const wet = offenders.get('wait_event_type', i) as string ?? '';
      const stateLabel = state === 'idle in transaction' ? 'idle in xact'
        : state === 'active' && wet === 'Lock' ? 'active (Lock)'
        : state;
      const cells: HTMLElement[] = [
        ui.divText(stateLabel, 'ua-metrics-tt-col-state'),
        ui.divText(DG.Utils.formatDuration((offenders.get('age_sec', i) as number) ?? 0),
          'ua-metrics-tt-col-num'),
        ui.divText((offenders.get('application_name', i) as string) || '—', 'ua-metrics-tt-col-app'),
        ui.divText((offenders.get('usename', i) as string) || '—', 'ua-metrics-tt-col-user'),
        ui.divText((offenders.get('query', i) as string) || '—', 'ua-metrics-tt-col-query'),
      ];
      if (showBlocks)
        cells.push(ui.divText(`${(offenders.get('blocks', i) as number) ?? 0}`, 'ua-metrics-tt-col-num'));
      rows.push(ui.divH(cells));
    }
    lines.push(ui.div(rows, 'ua-metrics-tt-table'));

    if (hasActiveLock) {
      lines.push(ui.divText(
        'active (Lock) = running a query but waiting for another session\'s lock.',
        'ua-metrics-tt-dim'));
    }
    if (showBlocks) {
      lines.push(ui.divText(
        '"blocks N" = N other sessions are queued behind this one\'s locks.',
        'ua-metrics-tt-dim'));
    }

    lines.push(ui.div([], {style: {height: '6px'}}));
    lines.push(ui.divText(
      'Run SELECT pg_terminate_backend(pid) to drop a session; ' +
      'fix missing COMMIT/ROLLBACK in the app for repeat offenders.',
      'ua-metrics-tt-dim'));
    return ui.div(lines, 'ua-metrics-tt ua-metrics-tt-wide');
  }

  private async loadErrors(): Promise<void> {
    const filter = this.uaToolbox.getFilter();
    const card = this.card('errors');
    ui.tooltip.bind(card.root, null);
    const df = await MetricsView.safeCall(() => queries.metricsErrorsCount(filter.date!), 'MetricsErrorsCount');
    if (!df || df.rowCount === 0) {
      this.setCard('errors', '0', 'no errors', 'info');
      return;
    }
    const now = Number(df.get('errors_now', 0) ?? 0);
    const prev = Number(df.get('errors_prev', 0) ?? 0);
    const delta = now - prev;
    const sub = delta === 0 ? 'no change' : `${delta > 0 ? '+' : ''}${delta} vs prev`;
    const color: CardColor = now === 0 ? 'info'
      : delta > 0 ? 'red'
      : delta < 0 ? 'green'
      : 'orange';
    this.setCard('errors', `${now}`, sub, color);
    const winStart = df.get('window_start', 0) as dayjs.Dayjs | null;
    const winEnd = df.get('window_end', 0) as dayjs.Dayjs | null;
    const prevStart = df.get('prev_window_start', 0) as dayjs.Dayjs | null;
    ui.tooltip.bind(card.root, () => MetricsView.buildErrorsTooltip({now, prev, delta, color, winStart, winEnd, prevStart}));
  }

  private static buildErrorsTooltip(d: {
    now: number, prev: number, delta: number, color: CardColor,
    winStart: dayjs.Dayjs | null, winEnd: dayjs.Dayjs | null, prevStart: dayjs.Dayjs | null,
  }): HTMLElement {
    const lines: HTMLElement[] = [];
    lines.push(ui.divText(
      `${d.now.toLocaleString()} error event${d.now === 1 ? '' : 's'} logged in this window.`));
    lines.push(ui.divText('An error event = one row in events with event_types.source = \'error\'.',
      'ua-metrics-tt-dim'));
    if (d.winStart && d.winEnd)
      lines.push(ui.divText(formatRange(d.winStart, d.winEnd), 'ua-metrics-tt-dim'));
    if (d.prevStart && d.winStart) {
      lines.push(ui.div([], {style: {height: '6px'}}));
      const pct = d.prev > 0 ? ` (${d.delta > 0 ? '+' : ''}${Math.round((d.delta / d.prev) * 100)}%)` : '';
      const change = d.delta === 0 ? 'no change'
        : `${d.delta > 0 ? '+' : ''}${d.delta}${pct} vs previous window`;
      lines.push(ui.divText(`Previous window: ${d.prev.toLocaleString()} error${d.prev === 1 ? '' : 's'} · ${change}`));
      lines.push(ui.divText(formatRange(d.prevStart, d.winStart), 'ua-metrics-tt-dim'));
      lines.push(ui.div([], {style: {height: '6px'}}));
      const trend = d.color === 'red' ? 'Red = more errors than the previous window.'
        : d.color === 'green' ? 'Green = fewer errors than the previous window.'
        : d.color === 'orange' ? 'Orange = same error count as the previous window.'
        : '"vs prev" compares to the immediately preceding window of equal length.';
      lines.push(ui.divText(trend, 'ua-metrics-tt-dim'));
    }
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async loadSessions(): Promise<void> {
    const filter = this.uaToolbox.getFilter();
    const card = this.card('sessions');
    ui.tooltip.bind(card.root, null);
    const df = await MetricsView.safeCall(() => queries.metricsSessionsCount(filter.date!), 'MetricsSessionsCount');
    if (!df || df.rowCount === 0) {
      this.setCard('sessions', '—', 'unavailable', 'info');
      return;
    }
    const now = Number(df.get('sessions_now', 0) ?? 0);
    const prev = Number(df.get('sessions_prev', 0) ?? 0);
    const delta = now - prev;
    const sub = delta === 0 ? 'no change' : `${delta > 0 ? '+' : ''}${delta} vs prev`;
    this.setCard('sessions', `${now}`, sub, 'info');
    const winStart = df.get('window_start', 0) as dayjs.Dayjs | null;
    const winEnd = df.get('window_end', 0) as dayjs.Dayjs | null;
    const prevStart = df.get('prev_window_start', 0) as dayjs.Dayjs | null;
    ui.tooltip.bind(card.root, () => MetricsView.buildSessionsTooltip({now, prev, delta, winStart, winEnd, prevStart}));
  }

  private static buildSessionsTooltip(d: {
    now: number, prev: number, delta: number,
    winStart: dayjs.Dayjs | null, winEnd: dayjs.Dayjs | null, prevStart: dayjs.Dayjs | null,
  }): HTMLElement {
    const lines: HTMLElement[] = [];
    lines.push(ui.divText(
      `${d.now.toLocaleString()} distinct user session${d.now === 1 ? '' : 's'} started in this window.`));
    lines.push(ui.divText('A session = one users_sessions row (one browser tab login).',
      'ua-metrics-tt-dim'));
    if (d.winStart && d.winEnd)
      lines.push(ui.divText(formatRange(d.winStart, d.winEnd), 'ua-metrics-tt-dim'));
    if (d.prevStart && d.winStart) {
      lines.push(ui.div([], {style: {height: '6px'}}));
      const pct = d.prev > 0 ? ` (${d.delta > 0 ? '+' : ''}${Math.round((d.delta / d.prev) * 100)}%)` : '';
      const change = d.delta === 0 ? 'no change'
        : `${d.delta > 0 ? '+' : ''}${d.delta}${pct} vs previous window`;
      lines.push(ui.divText(`Previous window: ${d.prev.toLocaleString()} session${d.prev === 1 ? '' : 's'} · ${change}`));
      lines.push(ui.divText(formatRange(d.prevStart, d.winStart), 'ua-metrics-tt-dim'));
      lines.push(ui.div([], {style: {height: '6px'}}));
      lines.push(ui.divText(
        '"vs prev" compares to the immediately preceding window of equal length.',
        'ua-metrics-tt-dim'));
    }
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async loadStorage(): Promise<void> {
    const card = this.card('storage');
    ui.tooltip.bind(card.root, null);
    const stats = await MetricsView.safeCall(() => grok.dapi.info.getStorageStats(), 'getStorageStats');
    if (!stats || Object.keys(stats).length === 0) {
      this.setCard('storage', '—', 'unavailable', 'info');
      return;
    }
    const more = stats.truncated ? '>' : '';
    const sub = `${stats.type} · ${more}${formatCount(stats.objectCount)} obj`;
    this.setCard('storage', `${more}${formatBytes(stats.totalBytes)}`, sub, 'info');
    const coords = await this.loadStorageCoords();
    ui.tooltip.bind(card.root, () => MetricsView.buildStorageTooltip(stats, coords));
  }

  private async loadStorageCoords(): Promise<string[]> {
    const conn = await MetricsView.safeCall(
      () => grok.dapi.connections.filter('namespace = "System:" and shortName = "AppData"').first(),
      'AppData connection');
    if (!conn) return [];
    const p = conn.parameters ?? {};
    return [p['bucket'], p['region']].filter((x) => x != null && x !== '');
  }

  private static buildStorageTooltip(s: any, coords: string[]): HTMLElement {
    const lines: HTMLElement[] = [];
    const headerCoords = coords.length > 0 ? coords : (s.root ? [s.root] : []);
    const more = s.truncated ? '>' : '';
    lines.push(ui.divText(`${s.type}${headerCoords.length > 0 ? ` · ${headerCoords.join(' · ')}` : ''}`));
    lines.push(ui.divText(`${more}${formatBytes(s.totalBytes)} · ${more}${formatCount(s.objectCount)} objects`));
    if (s.collectedAt) {
      const collected = dayjs(s.collectedAt);
      lines.push(ui.divText(`as of ${collected.format('YYYY-MM-DD HH:mm')} · ${formatAgo(collected)}`,
        'ua-metrics-tt-dim'));
    }
    const top = (s.topPrefixes as any[]) ?? [];
    if (top.length > 0) {
      const rows: HTMLElement[] = [
        ui.divH([
          ui.divText('prefix', 'ua-metrics-tt-col-prefix ua-metrics-tt-head'),
          ui.divText('size', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
          ui.divText('objects', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
          ui.divText('%', 'ua-metrics-tt-col-pct ua-metrics-tt-head'),
        ]),
      ];
      for (const p of top) {
        const pct = s.totalBytes > 0 ? `${((p.bytes / s.totalBytes) * 100).toFixed(1)}%` : '—';
        rows.push(ui.divH([
          ui.divText(p.name, 'ua-metrics-tt-col-prefix'),
          ui.divText(p.truncated ? '>?' : formatBytes(p.bytes), 'ua-metrics-tt-col-num'),
          ui.divText(p.truncated ? '>?' : formatCount(p.objectCount), 'ua-metrics-tt-col-num'),
          ui.divText(p.truncated ? '—' : pct, 'ua-metrics-tt-col-pct'),
        ], p.truncated ? 'ua-metrics-tt-dim' : ''));
      }
      lines.push(ui.div(rows, 'ua-metrics-tt-table'));
    }
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async loadDisk(): Promise<void> {
    const card = this.card('disk');
    ui.tooltip.bind(card.root, null);
    const s = await MetricsView.safeCall(async () => {
      const raw = await grok.functions.call('DiskStats') as string;
      return JSON.parse(raw) as DiskStat | null;
    }, 'DiskStats');
    if (!s) {
      this.setCard('disk', '—', 'unavailable', 'info');
      return;
    }
    this.setCard('disk', formatBytes(s.freeBytes), `${s.usedPct}% used`, MetricsView.diskColor(s));
    ui.tooltip.bind(card.root, () => MetricsView.buildDiskTooltip(s));
  }

  private static diskColor(s: DiskStat): CardColor {
    if (s.usedPct >= THRESH.diskUsedPct.red || s.freeBytes < THRESH.diskFreeBytes.red)
      return 'red';
    if (s.usedPct >= THRESH.diskUsedPct.orange || s.freeBytes < THRESH.diskFreeBytes.orange)
      return 'orange';
    return 'green';
  }

  private static buildDiskTooltip(s: DiskStat): HTMLElement {
    const mount = s.mount && s.mount !== s.path ? ` · mount ${s.mount}` : '';
    return ui.div([
      ui.divText(s.path),
      ui.divText(
        `${formatBytes(s.freeBytes)} free of ${formatBytes(s.totalBytes)} · ${s.usedPct}% used${mount}`,
        'ua-metrics-tt-dim'),
    ], 'ua-metrics-tt');
  }

  private async loadLatency(): Promise<void> {
    const filter = this.uaToolbox.getFilter();
    const df = await MetricsView.safeCall(() => grok.data.query('UsageAnalysis:MetricsLatency', {date: filter.date!}),
      'MetricsLatency');
    const card = this.card('latency');
    ui.tooltip.bind(card.root, null);
    if (!df || df.rowCount === 0 || (df.get('count_now', 0) as number) === 0) {
      this.setCard('latency', '—', 'no traffic', 'info');
      return;
    }
    const p50 = (df.get('p50_now', 0) as number) ?? 0;
    const p95 = (df.get('p95_now', 0) as number) ?? 0;
    const p99 = (df.get('p99_now', 0) as number) ?? 0;
    const p95Prev = (df.get('p95_prev', 0) as number) ?? 0;
    const countNow = (df.get('count_now', 0) as number) ?? 0;
    const countPrev = (df.get('count_prev', 0) as number) ?? 0;
    const winStart = df.get('window_start', 0) as dayjs.Dayjs | null;
    const winEnd = df.get('window_end', 0) as dayjs.Dayjs | null;
    const prevStart = df.get('prev_window_start', 0) as dayjs.Dayjs | null;

    const delta = p95 - p95Prev;
    const sub = `p95 ms · ${delta === 0 ? 'no change' : `${delta > 0 ? '+' : ''}${delta} vs prev`}`;
    const color: CardColor = p95 === 0 ? 'info'
      : delta > 0 ? 'red'
      : delta < 0 ? 'green'
      : 'orange';
    this.setCard('latency', `${p95}`, sub, color);

    ui.tooltip.bind(card.root, () => MetricsView.buildLatencyTooltip({
        p50, p95, p99, p95Prev, countNow, countPrev, winStart, winEnd, prevStart, color,
      }));
  }

  private static buildLatencyTooltip(d: {
    p50: number, p95: number, p99: number, p95Prev: number,
    countNow: number, countPrev: number,
    winStart: dayjs.Dayjs | null, winEnd: dayjs.Dayjs | null, prevStart: dayjs.Dayjs | null,
    color: CardColor,
  }): HTMLElement {
    const lines: HTMLElement[] = [];
    if (d.color === 'red')
      lines.push(ui.divText(`p95 above ${THRESH.latencyMs.orange} ms — most users see noticeable delay.`));
    else
      lines.push(ui.divText(`95% of requests finish under ${d.p95} ms.`));
    lines.push(ui.divText(`Based on ${d.countNow.toLocaleString()} request${d.countNow === 1 ? '' : 's'}`));
    if (d.winStart && d.winEnd)
      lines.push(ui.divText(formatRange(d.winStart, d.winEnd), 'ua-metrics-tt-dim'));
    if (d.countPrev > 0 && d.prevStart && d.winStart) {
      lines.push(ui.div([], {style: {height: '6px'}}));
      lines.push(ui.divText(`Previous window: p95 = ${d.p95Prev} ms (${d.countPrev.toLocaleString()} request${d.countPrev === 1 ? '' : 's'})`));
      lines.push(ui.divText(formatRange(d.prevStart, d.winStart), 'ua-metrics-tt-dim'));
    }
    lines.push(ui.div([], {style: {height: '6px'}}));
    lines.push(ui.divText(`p50 ${d.p50} ms · p99 ${d.p99} ms`, 'ua-metrics-tt-dim'));
    return ui.div(lines, 'ua-metrics-tt');
  }

  private async detectPssMode(): Promise<PssMode> {
    if (this.pssMode != null)
      return this.pssMode;
    const df = await MetricsView.safeCall(() => queries.metricsPgStatStatementsVersion(), 'MetricsPgStatStatementsVersion');
    if (!df || df.rowCount === 0) {
      this.pssMode = 'unavailable';
      return this.pssMode;
    }
    const v = String(df.get('extversion', 0) ?? '');
    const parts = v.split('.');
    const major = parseInt(parts[0], 10) || 0;
    const minor = parseInt(parts[1], 10) || 0;
    this.pssMode = (major > 1 || (major === 1 && minor >= 8)) ? 'modern' : 'legacy';
    return this.pssMode;
  }

  private static resetHost(host: HTMLElement): void {
    host.innerHTML = '';
    host.append(ui.loader());
  }

  private async loadQueries(): Promise<void> {
    MetricsView.resetHost(this.queriesGridHost);
    const mode = await this.detectPssMode();
    if (mode === 'unavailable') {
      this.queriesGridHost.innerHTML = '';
      this.queriesGridHost.append(MetricsView.degradedMessage('pg_stat_statements unavailable. ' +
        'Ensure the extension is installed and the System:Datagrok role has pg_read_all_stats.'));
      return;
    }
    const q = PSS_QUERIES[this.queriesMode];
    const fn = mode === 'legacy' ? q.legacy : q.modern;
    const df = await MetricsView.safeCall(() => fn(DASHBOARD_LIMIT), 'pg_stat_statements query');
    this.queriesGridHost.innerHTML = '';
    if (!df) {
      this.queriesGridHost.append(MetricsView.degradedMessage('pg_stat_statements query failed.'));
      return;
    }
    const grid = DG.Viewer.grid(df, MetricsView.gridOptions('13px monospace'));
    const queryCol = grid.col('query');
    if (queryCol)
      queryCol.width = 700;
    MetricsView.formatMsColumns(grid, ['total_ms', 'mean_ms']);
    grid.onCellPrepare((gc) => MetricsView.colorQueriesCell(gc));
    MetricsView.fitGrid(grid, this.queriesGridHost);
  }

  private async loadLargestTables(): Promise<void> {
    MetricsView.resetHost(this.largestTablesHost);
    const df = await MetricsView.safeCall(() => queries.metricsLargestTables(DASHBOARD_LIMIT), 'MetricsLargestTables');
    this.largestTablesHost.innerHTML = '';
    if (!df) {
      this.largestTablesHost.append(MetricsView.degradedMessage('Largest tables unavailable.'));
      return;
    }
    const grid = DG.Viewer.grid(df, MetricsView.gridOptions());
    const totalBytesCol = grid.col('total_bytes');
    if (totalBytesCol)
      totalBytesCol.visible = false;
    MetricsView.fitGrid(grid, this.largestTablesHost);
  }

  private async loadTableHealth(): Promise<void> {
    MetricsView.resetHost(this.tableHealthHost);
    const df = await MetricsView.safeCall(() => queries.metricsTableHealth(DASHBOARD_LIMIT), 'MetricsTableHealth');
    this.tableHealthHost.innerHTML = '';
    if (!df) {
      this.tableHealthHost.append(MetricsView.degradedMessage('Table health unavailable.'));
      return;
    }
    const grid = DG.Viewer.grid(df, MetricsView.gridOptions());
    const lastVacCol = grid.col('last_vacuum');
    if (lastVacCol)
      lastVacCol.name = 'last vacuum';
    grid.onCellPrepare((gc) => MetricsView.renderTableHealthCell(gc));
    MetricsView.fitGrid(grid, this.tableHealthHost);
  }

  private async loadHttpRoutes(): Promise<void> {
    MetricsView.resetHost(this.httpRoutesHost);
    const filter = this.uaToolbox.getFilter();
    const df = await MetricsView.safeCall(
      () => grok.data.query('UsageAnalysis:MetricsHttpRoutes', {date: filter.date!, limit: DASHBOARD_LIMIT}),
      'MetricsHttpRoutes');
    this.httpRoutesHost.innerHTML = '';
    if (!df) {
      this.httpRoutesHost.append(MetricsView.degradedMessage('HTTP routes unavailable.'));
      return;
    }
    if (df.rowCount === 0) {
      this.httpRoutesHost.append(MetricsView.degradedMessage('No http-request events in this window.'));
      return;
    }
    const grid = DG.Viewer.grid(df, MetricsView.gridOptions());
    const routeCol = grid.col('route');
    if (routeCol)
      routeCol.width = 360;
    MetricsView.formatMsColumns(grid, ['p50', 'p95', 'p99']);
    for (const p of ['p50', 'p95', 'p99']) {
      const c = grid.col(p);
      if (c)
        c.name = `${p}, ms`;
    }
    const errPctCol = grid.col('err_pct');
    if (errPctCol)
      errPctCol.format = '#0.##';
    grid.onCellPrepare((gc) => MetricsView.colorHttpRoutesCell(gc));
    MetricsView.fitGrid(grid, this.httpRoutesHost);
  }

  private async openFullHttpRoutes(): Promise<void> {
    const filter = this.uaToolbox.getFilter();
    const progress = DG.TaskBarProgressIndicator.create('Loading HTTP routes...');
    try {
      const df = await MetricsView.safeCall(
        () => grok.data.query('UsageAnalysis:MetricsHttpRoutes', {date: filter.date!, limit: FULL_VIEW_LIMIT}),
        'MetricsHttpRoutes');
      if (!df) {
        grok.shell.warning('HTTP routes unavailable.');
        return;
      }
      df.name = 'HTTP routes';
      grok.shell.addTableView(df);
    } finally {
      progress.close();
    }
  }

  private static degradedMessage(text: string): HTMLElement {
    return ui.div([ui.iconFA('info-circle'), ' ', text], 'ua-metrics-degraded');
  }

  private static formatMsColumns(grid: DG.Viewer, names: string[]): void {
    for (const name of names) {
      const col = (grid as any).col(name);
      if (col)
        col.format = '#,##0';
    }
  }

  private static fitGrid(grid: DG.Viewer, host: HTMLElement): void {
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    host.append(grid.root);
  }

  private static gridOptions(font?: string): any {
    const opts: any = {
      'showColumnLabels': true,
      'showRowHeader': true,
      'showColumnGridlines': false,
      'allowRowSelection': false,
      'allowBlockSelection': false,
      'allowRowResizing': true,
      'showCurrentCellOutline': false,
    };
    if (font)
      opts['defaultCellFont'] = font;
    return opts;
  }

  private static colorQueriesCell(gc: DG.GridCell): void {
    const name = gc.gridColumn.name;
    const v = gc.cell.value as number;
    if (v == null)
      return;
    if (name === 'mean_ms')
      gc.style.textColor = thresholdColor(v, THRESH.queryMeanMs, false);
    else if (name === 'hit_pct')
      gc.style.textColor = thresholdColor(v, THRESH.queryHitPct, true);
  }

  private static colorHttpRoutesCell(gc: DG.GridCell): void {
    const name = gc.gridColumn.name;
    const v = gc.cell.value as number;
    if (v == null)
      return;
    if (name === 'p95, ms')
      gc.style.textColor = thresholdColor(v, THRESH.routeP95Ms, false);
    else if (name === 'p99, ms')
      gc.style.textColor = thresholdColor(v, THRESH.routeP99Ms, false);
    else if (name === 'err_pct')
      gc.style.textColor = thresholdColor(v, THRESH.routeErrPct, false);
  }

  private static renderTableHealthCell(gc: DG.GridCell): void {
    const name = gc.gridColumn.name;
    if (name === 'dead_pct') {
      const v = gc.cell.value as number;
      if (v != null)
        gc.style.textColor = thresholdColor(v, THRESH.cellDeadPct, false);
    }
    else if (name === 'last vacuum') {
      const v = gc.cell.value as dayjs.Dayjs | null;
      gc.style.element = ui.divText(formatAgo(v));
    }
  }

  private async openFullView(fn: (limit: number) => Promise<DG.DataFrame>, queryName: string, title: string): Promise<void> {
    const progress = DG.TaskBarProgressIndicator.create(`Loading ${title}...`);
    try {
      const df = await MetricsView.safeCall(() => fn(FULL_VIEW_LIMIT), queryName);
      if (!df) {
        grok.shell.warning(`${title} unavailable.`);
        return;
      }
      df.name = title;
      grok.shell.addTableView(df);
    } finally {
      progress.close();
    }
  }

  private async sendToDatagrok(): Promise<void> {
    const progress = DG.TaskBarProgressIndicator.create('Preparing snapshot...');
    let files: DG.FileInfo[];
    try {
      files = await this.collectSnapshotAttachments();
    } catch (e) {
      grok.shell.error(`Failed to prepare snapshot: ${e}`);
      return;
    } finally {
      progress.close();
    }

    const reportEmail = await MetricsView.safeCall(() => grok.dapi.admin.getReportEmail(), 'getReportEmail');
    const subject = `Datagrok metrics — ${window.location.host} — ${dayjs().format('YYYY-MM-DD HH:mm')}`;
    ui.composeEmail({
      title: 'Send to Datagrok', subject, to: reportEmail ? reportEmail.split(/[,;\s]+/).filter((s) => s) : [],
      bcc: [], attachments: files,
      onSend: async (r) => {
        await grok.dapi.admin.sendEmail({
          subject: r.subject,
          to: r.to,
          bcc: r.bcc,
          html: r.html,
          text: r.text,
          attachments: r.attachments.map((fi) => ({name: fi.name, data: fi.data})),
        });
        grok.shell.info('Metrics snapshot sent.');
      },
    }).show();
  }

  private async collectSnapshotAttachments(): Promise<DG.FileInfo[]> {
    const filter = this.uaToolbox.getFilter();
    const date = filter.date!;
    const pssMode = await this.detectPssMode();
    const useLegacy = pssMode === 'legacy';

    type Src = {name: string, fn: () => Promise<DG.DataFrame | string | null>};
    const sources: Src[] = [
      {name: 'db_summary.csv',           fn: () => queries.metricsDbStats()},
      {name: 'table_health.csv',         fn: () => queries.metricsTableHealth(EMAIL_LIMIT)},
      {name: 'largest_tables.csv',       fn: () => queries.metricsLargestTables(20)},
      {name: 'connections.csv',          fn: () => queries.metricsConnections()},
      {name: 'connections_offenders.csv', fn: () => queries.metricsConnectionsOffenders(20, 60, 30)},
      {name: 'http_routes.csv',          fn: () => queries.metricsHttpRoutes(date, EMAIL_LIMIT)},
      {name: 'errors.csv',               fn: () => queries.metricsErrorsCount(date)},
      {name: 'latency.csv',              fn: () => queries.metricsLatency(date)},
      {name: 'sessions.csv',             fn: () => queries.metricsSessionsCount(date)},
    ];
    if (pssMode !== 'unavailable') {
      sources.push(
        {name: 'slowest_queries.csv',    fn: () => useLegacy ? queries.metricsTopSlowestQueriesPg12(EMAIL_LIMIT) : queries.metricsTopSlowestQueries(EMAIL_LIMIT)},
        {name: 'most_called_queries.csv', fn: () => useLegacy ? queries.metricsTopMostCalledQueriesPg12(EMAIL_LIMIT) : queries.metricsTopMostCalledQueries(EMAIL_LIMIT)},
        {name: 'worst_cache_hit.csv',    fn: () => useLegacy ? queries.metricsWorstCacheHitQueriesPg12(EMAIL_LIMIT) : queries.metricsWorstCacheHitQueries(EMAIL_LIMIT)},
      );
    }
    sources.push(
      {name: 'storage.json',             fn: () => grok.functions.call('StorageStats') as Promise<string>},
      {name: 'disk.json',                fn: () => grok.functions.call('DiskStats') as Promise<string>},
    );

    const results = await Promise.all(sources.map((s) => MetricsView.safeCall(s.fn, s.name)));

    const files: DG.FileInfo[] = [];
    for (let i = 0; i < sources.length; i++) {
      const r = results[i];
      if (r == null)
        continue;
      const csv = typeof r === 'string' ? r : (r as DG.DataFrame).toCsv();
      files.push(DG.FileInfo.fromString(sources[i].name, csv));
    }
    return files;
  }
}
