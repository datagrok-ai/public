/** See {@link ./METRICS.md} for developer documentation. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';

import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';
import {waitForWidth} from '../release/data';

type CardColor = 'green' | 'orange' | 'red' | 'info';
type QueriesMode = 'slowest' | 'mostCalled' | 'worstCacheHit';
type GrowthArea = 'database' | 'storage' | 'gc';

const DASHBOARD_LIMIT = 10;
const FULL_VIEW_LIMIT = 100000;
const EMAIL_LIMIT = 100;
const GROWTH_DAYS = [7, 30, 90];
const TOOLTIP_TOP = 5;

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
  growthPct: {green: 10, orange: 25},
};

interface DiskStat {
  path: string;
  mount: string;
  totalBytes: number;
  usedBytes: number;
  freeBytes: number;
  usedPct: number;
}

const COLOR_GREEN = DG.Color.success;
const COLOR_ORANGE = 0xFFFFA24A;
const COLOR_RED = DG.Color.failure;

function num(df: DG.DataFrame, col: string, row = 0): number {
  return Number(df.get(col, row) ?? 0);
}

function dts(df: DG.DataFrame, col: string, row = 0): dayjs.Dayjs | null {
  return df.get(col, row) as dayjs.Dayjs | null;
}

function formatTime(d: Date): string {
  return d.toTimeString().substring(0, 5);
}

function thresholdBand(v: number, t: {green: number, orange: number}, higherIsBetter: boolean): CardColor {
  const isGreen = higherIsBetter ? v >= t.green : v < t.green;
  const isOrange = higherIsBetter ? v >= t.orange : v < t.orange;
  return isGreen ? 'green' : isOrange ? 'orange' : 'red';
}

function thresholdColor(v: number, t: {green: number, orange: number}, higherIsBetter: boolean): number {
  const band = thresholdBand(v, t, higherIsBetter);
  return band === 'green' ? COLOR_GREEN : band === 'orange' ? COLOR_ORANGE : COLOR_RED;
}

function deltaSub(delta: number): string {
  return delta === 0 ? 'no change' : `${delta > 0 ? '+' : ''}${delta} vs prev`;
}

function deltaColor(now: number, delta: number): CardColor {
  return now === 0 ? 'info' : delta > 0 ? 'red' : delta < 0 ? 'green' : 'orange';
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

function formatSigned(n: number | null | undefined, format: (n: number) => string): string {
  return n == null ? '—' : `${n > 0 ? '+' : n < 0 ? '-' : ''}${format(Math.abs(n))}`;
}

/** Growth of [s] over 7 days, percent of the value a week ago; null when there is no sample that old. */
function growthPct(s: DG.GrowthMetricsSeries | null): number | null {
  if (s == null || s.d7 == null)
    return null;
  const before = s.now - s.d7;
  return before > 0 ? s.d7 * 100 / before : null;
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
  private metrics: DG.ServerMetrics | null = null;
  private growthHost!: HTMLElement;
  private growth: DG.GrowthMetrics | null = null;
  private growthArea: GrowthArea = 'database';
  private growthDays = 30;
  private onlyGrowing = false;

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
    this.growthHost = ui.div([], 'ua-metrics-grid-host ua-metrics-growth-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');

    const header = ui.divH([
      this.asOfLabel,
      ui.div([], {style: {flex: '1'}}),
      refreshBtn,
    ], 'ua-metrics-header');

    const cardsRow1 = ui.divH([
      this.makeCard('db', 'Database'),
      this.makeCard('tblHealth', 'Table health'),
      this.makeCard('storage', 'Storage'),
      this.makeCard('connections', 'Connections'),
      this.makeCard('queue', 'Queue'),
      this.makeCard('disk', 'Disk free'),
    ], 'ua-metrics-cards-row');

    const cardsRow2 = ui.divH([
      this.makeCard('errors', 'Errors'),
      this.makeCard('latency', 'Latency'),
      this.makeCard('sessions', 'Sessions'),
    ], 'ua-metrics-cards-row');

    const cardsRow3 = ui.divH([
      this.makeCard('dbGrowth', 'DB growth 7d'),
      this.makeCard('storageGrowth', 'Storage growth 7d'),
      this.makeCard('gcDeleted', 'GC deleted'),
      this.makeCard('gcPending', 'GC pending'),
      this.makeCard('growthAlerts', 'Growth alerts'),
    ], 'ua-metrics-cards-row');

    const queriesPanel = this.buildQueriesPanel();
    const largestTablesPanel = MetricsView.buildPanel('Largest tables', this.largestTablesHost,
      () => this.openFullView(queries.metricsLargestTables, 'MetricsLargestTables', 'Largest tables'));
    const tableHealthPanel = MetricsView.buildPanel('Table health', this.tableHealthHost,
      () => this.openFullView(queries.metricsTableHealth, 'MetricsTableHealth', 'Table health'));
    const httpRoutesPanel = MetricsView.buildPanel('HTTP routes', this.httpRoutesHost,
      () => this.openFullView(async (limit) => MetricsView.routesFrame((await this.fetchMetrics(limit)).http.routes),
        'getMetrics', 'HTTP routes'));

    const tablesRow = ui.divH([largestTablesPanel, tableHealthPanel], 'ua-metrics-tables-row');
    const growthPanel = this.buildGrowthPanel();

    this.root.append(ui.div([
      header,
      cardsRow1,
      cardsRow2,
      cardsRow3,
      growthPanel,
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

  private static buildPanel(title: string, host: HTMLElement, onOpenFull?: () => void,
    extraHeader?: HTMLElement, actions?: HTMLElement): HTMLElement {
    const headerChildren: HTMLElement[] = [ui.divText(title, 'ua-metrics-panel-title')];
    if (extraHeader)
      headerChildren.push(extraHeader);
    headerChildren.push(ui.div([], {style: {flex: '1'}}));
    if (actions)
      headerChildren.push(actions);
    else if (onOpenFull) {
      const addIcon = ui.iconFA('plus', onOpenFull, 'Add to workspace');
      addIcon.classList.add('ua-metrics-add-icon');
      headerChildren.push(addIcon);
    }
    const header = ui.divH(headerChildren, 'ua-metrics-panel-header');
    return ui.div([header, host], 'ua-metrics-panel');
  }

  private buildQueriesPanel(): HTMLElement {
    const modes: Array<[string, QueriesMode, string]> = [
      ['slowest', 'slowest', 'Highest total execution time — biggest aggregate load.'],
      ['most called', 'mostCalled', 'Highest call count — small per-call costs add up.'],
      ['worst cache hit', 'worstCacheHit',
        'Lowest cache hit ratio (queries with ≥ 100 calls) — disk reads instead of memory; index or RAM candidates.'],
    ];
    const buttons = modes.map(([label, mode, tooltip]) => {
      const btn = ui.toggleButton(label, () => {
        this.queriesMode = mode;
        this.loadQueries();
      }, tooltip);
      if (mode === this.queriesMode)
        btn.classList.add('d4-current');
      return btn;
    });
    const toggle = ui.toggleButtonGroup(buttons);
    toggle.classList.add('ua-metrics-mode-toggle');

    const more = ui.iconFA('ellipsis-h', (e: MouseEvent) => {
      e.stopImmediatePropagation();
      const menu = DG.Menu.popup();
      menu.item('Add to workspace', () => this.openFullView(
        async (limit) => MetricsView.statementsFrame((await this.fetchMetrics(limit)).database.statements[this.queriesMode]),
        'getMetrics', 'pg_stat_statements'));
      menu.item('Reset stats', () => this.confirmResetPgStats());
      menu.show({causedBy: e});
    }, 'More actions');
    more.classList.add('ua-metrics-add-icon');

    return MetricsView.buildPanel('Queries', this.queriesGridHost, undefined, toggle, more);
  }

  private confirmResetPgStats(): void {
    const body = ui.div([
      ui.divText(
        'This calls pg_stat_statements_reset() — it clears accumulated query stats for ' +
        'all administrators, not just your view. It can\'t be undone.'),
      ui.divText('The slowest-query rankings will be empty until traffic rebuilds them.'),
    ], 'ua-metrics-reset-body');

    const dialog = ui.dialog({title: 'Reset query statistics?'}).add(body);
    dialog.onOK(async () => {
      try {
        await grok.data.query('UsageAnalysis:MetricsResetPgStatStatements');
        grok.shell.info('pg_stat_statements stats reset.');
        await this.loadMetrics();
      } catch (e) {
        grok.shell.error(`Failed to reset stats: ${e}`);
      }
    });
    dialog.show();

    const okBtn = dialog.getButton('OK');
    if (okBtn) {
      okBtn.textContent = 'Reset stats';
      okBtn.classList.add('ua-metrics-btn-danger');
    }
  }

  async exportFiles(): Promise<DG.FileInfo[]> {
    return this.collectSnapshotAttachments();
  }

  private async refresh(): Promise<void> {
    if (this.refreshing)
      return;
    this.refreshing = true;
    try {
      this.asOfLabel.textContent = `as of ${formatTime(new Date())}`;
      await Promise.all([
        this.loadTableHealthSummary(),
        this.loadLargestTables(),
        this.loadTableHealth(),
        this.refreshWindowCards(),
        this.loadStorage(),
        this.loadDisk(),
        this.loadGrowth(),
      ]);
    } finally {
      this.refreshing = false;
    }
  }

  private async refreshWindowCards(): Promise<void> {
    await Promise.all([this.loadErrors(), this.loadSessions(), this.loadMetrics()]);
  }

  private fetchMetrics(limit: number): Promise<DG.ServerMetrics> {
    return grok.dapi.admin.getMetrics({date: this.uaToolbox.getFilter().date!, limit});
  }

  private async loadMetrics(): Promise<void> {
    MetricsView.resetHost(this.queriesGridHost);
    MetricsView.resetHost(this.httpRoutesHost);
    this.metrics = await MetricsView.safeCall(() => this.fetchMetrics(DASHBOARD_LIMIT), 'getMetrics');
    this.loadLatency();
    this.loadHttpRoutes();
    this.loadQueue();
    this.loadConnections();
    this.loadQueries();
    await this.loadDbStats();
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
    const db = this.metrics?.database;
    if (!db) {
      this.setCard('db', '—', 'unavailable', 'info');
      return;
    }
    const hit = db.cacheHitPct;
    this.setCard('db', formatBytes(db.sizeBytes), `hit ${hit.toFixed(0)}%`,
      thresholdBand(hit, THRESH.cacheHit, true));
    const offenders = await MetricsView.safeCall(() => queries.metricsCacheMissTables(), 'MetricsCacheMissTables');
    ui.tooltip.bind(card.root, () => MetricsView.buildDbStatsTooltip(db, offenders));
  }

  private static buildDbStatsTooltip(db: DG.ServerMetrics['database'], offenders: DG.DataFrame | null): HTMLElement {
    const lines: HTMLElement[] = [];
    lines.push(ui.divText(`${formatBytes(db.sizeBytes)} · cache hit ${db.cacheHitPct.toFixed(2)}%`));
    lines.push(ui.divText(
      'Cache hit = share of heap-block reads served from shared_buffers instead of disk. ' +
      '>99% is healthy for OLTP; lower means more disk I/O.',
      'ua-metrics-tt-dim'));
    if (db.statsReset) {
      const statsReset = dayjs(db.statsReset);
      lines.push(ui.divText(`Stats since ${statsReset.format('YYYY-MM-DD HH:mm')} · ${formatAgo(statsReset)}`,
        'ua-metrics-tt-dim'));
    }

    lines.push(ui.div([], {style: {height: '6px'}}));
    if (!offenders)
      lines.push(ui.divText('Cache-miss tables unavailable.', 'ua-metrics-tt-dim'));
    else if (offenders.rowCount === 0)
      lines.push(ui.divText('No user tables with notable traffic below 95% hit ratio.', 'ua-metrics-tt-dim'));
    else {
      lines.push(ui.divText('Top tables by cache misses (hit ratio < 95%):'));
      const rows: HTMLElement[] = [
        ui.divH([
          ui.divText('table', 'ua-metrics-tt-col-prefix ua-metrics-tt-head'),
          ui.divText('hit %', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
          ui.divText('disk reads', 'ua-metrics-tt-col-num ua-metrics-tt-head'),
        ]),
      ];
      for (let i = 0; i < offenders.rowCount; i++) {
        rows.push(ui.divH([
          ui.divText(offenders.get('offender_table', i) as string, 'ua-metrics-tt-col-prefix'),
          ui.divText(`${Number(offenders.get('offender_hit_pct', i) ?? 0).toFixed(1)}%`, 'ua-metrics-tt-col-num'),
          ui.divText(formatCount(Number(offenders.get('offender_disk_reads', i) ?? 0)), 'ua-metrics-tt-col-num'),
        ]));
      }
      lines.push(ui.div(rows, 'ua-metrics-tt-table'));
      lines.push(ui.divText('Hit % includes heap + index blocks. Filtered to tables with ≥1K disk reads.',
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

  private loadConnections(): void {
    const card = this.card('connections');
    ui.tooltip.bind(card.root, null);
    const c = this.metrics?.database.connections;
    if (!c) {
      this.setCard('connections', '—', 'unavailable', 'info');
      return;
    }
    const idleX = c.idleInTransaction;
    this.setCard('connections', `${c.total}`, idleX > 0 ? `${idleX} idle xact` : 'no idle xact',
      thresholdBand(c.oldestIdleTransactionSeconds, THRESH.idleXactSec, false));

    let offenders: DG.DataFrame | null = null;
    let pending = true;
    ui.tooltip.bind(card.root, () => MetricsView.buildConnectionsTooltip(c, offenders, pending));
    setTimeout(async () => {
      offenders = await MetricsView.safeCall(
        () => queries.metricsConnectionsOffenders(DASHBOARD_LIMIT, THRESH.idleXactSec.green, 30),
        'MetricsConnectionsOffenders');
      pending = false;
    }, 500);
  }

  private static buildConnectionsTooltip(
    c: DG.ServerMetricsConnections, offenders: DG.DataFrame | null, pending: boolean): HTMLElement {
    const idleX = c.idleInTransaction;
    const lines: HTMLElement[] = [];
    const headline = idleX > 0
      ? `${c.total} connections · ${idleX} idle in transaction · oldest ${DG.Utils.formatDuration(c.oldestIdleTransactionSeconds)} ago`
      : `${c.total} connections`;
    lines.push(ui.divText(headline));
    lines.push(ui.divText(
      `${c.active} active · ${c.idle} idle · ${idleX} idle in transaction · ${c.waitingOnLock} waiting on lock`,
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

  private loadQueue(): void {
    const q = this.metrics?.queue;
    if (!q) {
      this.setCard('queue', '—', 'unavailable', 'info');
      return;
    }
    this.setCard('queue', `${q.queued + q.running}`, `${q.queued} queued · ${q.running} running`,
      q.queued === 0 ? 'green' : 'orange');
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
    const now = num(df, 'errors_now');
    const prev = num(df, 'errors_prev');
    const delta = now - prev;
    const color = deltaColor(now, delta);
    this.setCard('errors', `${now}`, deltaSub(delta), color);
    const winStart = dts(df, 'window_start');
    const winEnd = dts(df, 'window_end');
    const prevStart = dts(df, 'prev_window_start');
    ui.tooltip.bind(card.root, () => MetricsView.buildErrorsTooltip({now, prev, delta, color, winStart, winEnd, prevStart}));
  }

  // The "Previous window: N nouns · ±delta (pct) vs previous window" block shared by the
  // errors and sessions cards. `trend` is the closing dim explainer (color-specific for errors).
  private static prevCountBlock(d: {
    noun: string, prev: number, delta: number, trend: string,
    prevStart: dayjs.Dayjs | null, winStart: dayjs.Dayjs | null,
  }): HTMLElement[] {
    if (!(d.prevStart && d.winStart))
      return [];
    const pct = d.prev > 0 ? ` (${d.delta > 0 ? '+' : ''}${Math.round((d.delta / d.prev) * 100)}%)` : '';
    const change = d.delta === 0 ? 'no change'
      : `${d.delta > 0 ? '+' : ''}${d.delta}${pct} vs previous window`;
    return [
      ui.div([], {style: {height: '6px'}}),
      ui.divText(`Previous window: ${d.prev.toLocaleString()} ${d.noun}${d.prev === 1 ? '' : 's'} · ${change}`),
      ui.divText(formatRange(d.prevStart, d.winStart), 'ua-metrics-tt-dim'),
      ui.div([], {style: {height: '6px'}}),
      ui.divText(d.trend, 'ua-metrics-tt-dim'),
    ];
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
    const trend = d.color === 'red' ? 'Red = more errors than the previous window.'
      : d.color === 'green' ? 'Green = fewer errors than the previous window.'
      : d.color === 'orange' ? 'Orange = same error count as the previous window.'
      : '"vs prev" compares to the immediately preceding window of equal length.';
    lines.push(...MetricsView.prevCountBlock(
      {noun: 'error', prev: d.prev, delta: d.delta, prevStart: d.prevStart, winStart: d.winStart, trend}));
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
    const now = num(df, 'sessions_now');
    const prev = num(df, 'sessions_prev');
    const delta = now - prev;
    this.setCard('sessions', `${now}`, deltaSub(delta), 'info');
    const winStart = dts(df, 'window_start');
    const winEnd = dts(df, 'window_end');
    const prevStart = dts(df, 'prev_window_start');
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
    lines.push(...MetricsView.prevCountBlock({
      noun: 'session', prev: d.prev, delta: d.delta, prevStart: d.prevStart, winStart: d.winStart,
      trend: '"vs prev" compares to the immediately preceding window of equal length.'}));
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

  private loadLatency(): void {
    const card = this.card('latency');
    ui.tooltip.bind(card.root, null);
    const m = this.metrics;
    if (!m) {
      this.setCard('latency', '—', 'unavailable', 'info');
      return;
    }
    const {p50, p95, p99, count: countNow} = m.http.now;
    if (countNow === 0) {
      this.setCard('latency', '—', 'no traffic', 'info');
      return;
    }
    const p95Prev = m.http.previous.p95;
    const countPrev = m.http.previous.count;
    const winStart = dayjs(m.window.start);
    const winEnd = dayjs(m.window.end);
    const prevStart = dayjs(m.window.previousStart);

    const delta = p95 - p95Prev;
    const color = deltaColor(p95, delta);
    this.setCard('latency', `${p95}`, `p95 ms · ${deltaSub(delta)}`, color);

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

  private static resetHost(host: HTMLElement): void {
    host.innerHTML = '';
    host.append(ui.loader());
  }

  private loadQueries(): void {
    this.queriesGridHost.innerHTML = '';
    const statements = this.metrics?.database.statements;
    if (!statements) {
      this.queriesGridHost.append(MetricsView.degradedMessage('Query statistics unavailable.'));
      return;
    }
    if (!statements.available) {
      this.queriesGridHost.append(MetricsView.degradedMessage('pg_stat_statements unavailable. ' +
        'Ensure the extension is installed and the System:Datagrok role has pg_read_all_stats.'));
      return;
    }
    const grid = DG.Viewer.grid(MetricsView.statementsFrame(statements[this.queriesMode]),
      MetricsView.gridOptions('13px monospace'));
    const queryCol = grid.col('query');
    if (queryCol)
      queryCol.width = 700;
    MetricsView.formatMsColumns(grid, ['calls', 'total_ms', 'mean_ms']);
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
    const rowsCol = grid.col('#rows');
    if (rowsCol)
      rowsCol.format = '#,##0';
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

  private loadHttpRoutes(): void {
    this.httpRoutesHost.innerHTML = '';
    const routes = this.metrics?.http.routes;
    if (!routes) {
      this.httpRoutesHost.append(MetricsView.degradedMessage('HTTP routes unavailable.'));
      return;
    }
    if (routes.length === 0) {
      this.httpRoutesHost.append(MetricsView.degradedMessage('No requests recorded in this window.'));
      return;
    }
    const grid = DG.Viewer.grid(MetricsView.routesFrame(routes), MetricsView.gridOptions());
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
    const g = grid as DG.Grid;
    g.onCellTooltip((gc, x, y) => {
      if (!gc.isTableCell || gc.gridColumn.name !== 'err_pct' || !(gc.cell.value > 0))
        return false;
      ui.tooltip.show('Click to see the failed requests, their errors and log', x, y);
      return true;
    });
    g.onCellClick.subscribe((gc) => {
      if (gc.isTableCell && gc.gridColumn.name === 'err_pct' && gc.cell.value > 0)
        this.openRouteFailures(routes[gc.tableRowIndex!]);
    });
    MetricsView.fitGrid(grid, this.httpRoutesHost);
  }

  private buildGrowthPanel(): HTMLElement {
    const areas: Array<[string, GrowthArea, string]> = [
      ['database', 'database',
        'Platform database tables, and their categories: events by source and type, entities by type, ...'],
      ['storage', 'storage', 'File storage by location and kind of data, and each server\'s disk.'],
      ['gc', 'gc', 'What the garbage collector deleted and has left to delete.'],
    ];
    const toggle = ui.toggleButtonGroup(areas.map(([label, area, tooltip]) => {
      const btn = ui.toggleButton(label, () => {
        this.growthArea = area;
        this.renderGrowth();
      }, tooltip);
      if (area === this.growthArea)
        btn.classList.add('d4-current');
      return btn;
    }));
    toggle.classList.add('ua-metrics-mode-toggle');
    const growing = ui.input.bool('Only growing', {value: this.onlyGrowing,
      tooltipText: 'Only what grew over the last 7 days', onValueChanged: (v) => {
        this.onlyGrowing = v;
        this.renderGrowth();
      }});
    const days = ui.input.choice<number>('Days', {value: this.growthDays, items: GROWTH_DAYS, nullable: false,
      tooltipText: 'Days of the trend', onValueChanged: (v) => {
        this.growthDays = v ?? 30;
        this.loadGrowth();
      }});
    const controls = ui.divH([toggle, growing.root, days.root], 'ua-metrics-growth-controls');
    return MetricsView.buildPanel('Data growth', this.growthHost,
      () => this.openFullView(async () => this.growthFrame(this.growthRows()), 'getGrowthMetrics', 'Data growth'),
      controls);
  }

  private async loadGrowth(): Promise<void> {
    MetricsView.resetHost(this.growthHost);
    this.growth = await MetricsView.safeCall(() => grok.dapi.admin.getGrowthMetrics({days: this.growthDays}),
      'getGrowthMetrics');
    this.loadGrowthCards();
    this.renderGrowth();
  }

  /** The rows of the current area, the growing ones only when asked, fastest-growing first. */
  private growthRows(): DG.GrowthMetricsRow[] {
    const rows = (this.growth?.rows ?? []).filter((r) => r.area === this.growthArea &&
      (!this.onlyGrowing || (r.bytes?.d7 ?? 0) > 0 || (r.count?.d7 ?? 0) > 0));
    return rows.sort((a, b) => ((b.bytes?.d7 ?? -Infinity) - (a.bytes?.d7 ?? -Infinity)) ||
      ((b.bytes?.now ?? b.count?.now ?? 0) - (a.bytes?.now ?? a.count?.now ?? 0)));
  }

  private static growthName(r: DG.GrowthMetricsRow): string {
    if (r.area !== 'storage')
      return r.dimension == null ? r.key : `${r.table} · ${r.dimension}=${r.value ?? ''}`;
    return `${r.location} · ${r.nature ?? 'all'}${r.value != null ? `/${r.value}` : ''}`;
  }

  /** One row per [rows], typed: sizes and counts are float columns, since a size passes 2^31. */
  private growthFrame(rows: DG.GrowthMetricsRow[]): DG.DataFrame {
    const text = (name: string, f: (r: DG.GrowthMetricsRow) => string | null) =>
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, rows.map((r) => f(r) ?? ''));
    const number = (name: string, f: (r: DG.GrowthMetricsRow) => number | null | undefined) =>
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, rows.map((r) => f(r) ?? null));
    const columns = [
      text('name', (r) => MetricsView.growthName(r)), text('area', (r) => r.area), text('location', (r) => r.location),
      text('nature', (r) => r.nature), text('table', (r) => r.table), text('dimension', (r) => r.dimension),
      text('value', (r) => r.value), text('key', (r) => r.key),
      number('size', (r) => r.bytes?.now), number('count', (r) => r.count?.now),
      number('d1_bytes', (r) => r.bytes?.d1), number('d7_bytes', (r) => r.bytes?.d7),
      number('d30_bytes', (r) => r.bytes?.d30),
      number('d7_pct', (r) => growthPct(r.bytes ?? r.count)), number('d7_count', (r) => r.count?.d7),
      number('written_7d', (r) => r.count?.d7 == null ? null : r.count.d7 + (r.gc?.deleted7d ?? 0)),
      number('gc_last', (r) => r.gc?.deletedLast), number('gc_7d', (r) => r.gc?.deleted7d),
      number('gc_pending', (r) => r.gc?.pending),
      number('free_bytes', (r) => r.freeBytes), number('used_pct', (r) => r.usedPct),
      text('alert', (r) => r.alerts.filter((a) => a.holding).map((a) => a.name).join(', ')),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'estimated',
        rows.map((r) => (r.bytes?.estimated || r.count?.estimated) ?? false)),
    ];
    for (let i = 0; i < this.growthDays; i++)
      columns.push(number(`t${i}`, (r) => (r.bytes ?? r.count)?.trend[i]));
    const df = DG.DataFrame.fromColumns(columns);
    df.name = 'Data growth';
    return df;
  }


  private renderGrowth(): void {
    this.growthHost.innerHTML = '';
    if (!this.growth) {
      this.growthHost.append(MetricsView.degradedMessage('Data growth unavailable.'));
      return;
    }
    if (this.growth.collectedAt == null) {
      this.growthHost.append(MetricsView.degradedMessage('Collecting the first samples…'));
      return;
    }
    const rows = this.growthRows();
    if (rows.length === 0) {
      this.growthHost.append(MetricsView.degradedMessage(this.onlyGrowing ?
        'Nothing grew over the last 7 days.' : 'No samples in this area yet.'));
      return;
    }
    const df = this.growthFrame(rows);
    const grid = DG.Viewer.grid(df, MetricsView.gridOptions());
    for (const name of ['count', 'd7_count', 'written_7d', 'gc_last', 'gc_7d', 'gc_pending']) {
      const col = grid.col(name);
      if (col)
        col.format = '#,##0';
    }
    const pct = grid.col('d7_pct');
    if (pct)
      pct.format = '#0.#';
    const nameCol = grid.col('name');
    if (nameCol)
      nameCol.width = 260;
    grid.onCellPrepare((gc) => MetricsView.renderGrowthCell(gc, rows));
    MetricsView.bindGrowthAlerts(grid, rows);
    const trendCols = Array.from({length: this.growthDays}, (_, i) => `t${i}`);
    const area = this.growthArea;
    const visible = ['name', 'size', 'count', 'd7_bytes', 'd7_pct', 'trend',
      ...(area === 'storage' ? ['free_bytes', 'used_pct'] : ['written_7d', 'gc_7d', 'gc_pending']), 'alert'];
    waitForWidth(grid.root).then(() => {
      const trend = grid.columns.add({gridColumnName: 'trend', cellType: 'sparkline'});
      trend.settings = {columnNames: trendCols};
      grid.columns.setVisible(visible.filter((n) => grid.col(n) != null));
    });
    MetricsView.fitGrid(grid, this.growthHost);
  }

  private static renderGrowthCell(gc: DG.GridCell, rows: DG.GrowthMetricsRow[]): void {
    if (!gc.isTableCell)
      return;
    const name = gc.gridColumn.name;
    const v = gc.cell.value;
    const row = rows[gc.tableRowIndex!];
    if (row == null)
      return;
    const estimated = row.count?.estimated ? '~' : '';
    if (name === 'size' || name === 'd1_bytes' || name === 'd7_bytes' || name === 'd30_bytes' || name === 'free_bytes')
      gc.customText = v == null ? '' :
        name === 'size' || name === 'free_bytes' ? formatBytes(v) : formatSigned(v, formatBytes);
    else if (name === 'count' && v != null)
      gc.customText = `${(v as number).toLocaleString()}${estimated}`;
    else if (name === 'd7_pct' && v != null) {
      gc.customText = `${v > 0 ? '+' : ''}${(v as number).toFixed(v >= 10 || v <= -10 ? 0 : 1)}%`;
      gc.style.textColor = thresholdColor(v as number, THRESH.growthPct, false);
    }
    else if (name === 'alert') {
      gc.customText = row.alerts.some((a) => a.holding) ? '●' : '';
      gc.style.textColor = COLOR_RED;
    }
  }

  /** The alerts holding for the row of [gc] as a tooltip; a click opens the first one's problem in Alerts (its id is
   *  only there for `ManageAlerts` holders). */
  private static bindGrowthAlerts(grid: DG.Grid, rows: DG.GrowthMetricsRow[]): void {
    const holding = (gc: DG.GridCell) => !gc.isTableCell || gc.gridColumn.name !== 'alert' ? [] :
      (rows[gc.tableRowIndex!]?.alerts ?? []).filter((a) => a.holding);
    grid.onCellTooltip((gc, x, y) => {
      const alerts = holding(gc);
      if (alerts.length === 0)
        return false;
      const linked = alerts.some((a) => a.problemId != null);
      ui.tooltip.show(ui.div([...alerts.map((a) => ui.divText(a.name)),
        ...(linked ? [ui.divText('Click to open the problem in Alerts.', 'ua-metrics-tt-dim')] : [])],
      'ua-metrics-tt'), x, y);
      return true;
    });
    grid.onCellClick.subscribe((gc) => {
      const id = holding(gc).find((a) => a.problemId != null)?.problemId;
      if (id)
        grok.shell.route(`/alerts/${id}`);
    });
  }

  private loadGrowthCards(): void {
    for (const id of ['dbGrowth', 'storageGrowth', 'gcDeleted', 'gcPending', 'growthAlerts'])
      ui.tooltip.bind(this.card(id).root, null);
    const g = this.growth;
    if (!g) {
      for (const id of ['dbGrowth', 'storageGrowth', 'gcDeleted', 'gcPending', 'growthAlerts'])
        this.setCard(id, '—', 'unavailable', 'info');
      return;
    }
    if (g.collectedAt == null) {
      for (const id of ['dbGrowth', 'storageGrowth', 'gcDeleted', 'gcPending', 'growthAlerts'])
        this.setCard(id, '—', 'collecting…', 'info');
      return;
    }
    const total = g.rows.find((r) => r.area === 'database' && r.key === '(total)');
    this.setGrowthCard('dbGrowth', total?.bytes ?? null, `now ${formatBytes(total?.bytes?.now)}`,
      g.rows.filter((r) => r.area === 'database' && r.dimension == null && r.key !== '(total)'));
    const shared = g.rows.find((r) => r.area === 'storage' && r.key === 'shared');
    this.setGrowthCard('storageGrowth', shared?.bytes ?? null,
      `${g.storage?.type ?? 'storage'} · ${formatBytes(shared?.bytes?.now)}`,
      g.rows.filter((r) => r.area === 'storage' && r.location === 'shared' && r.nature != null && r.value == null));

    const gcRows = g.rows.filter((r) => r.area === 'gc' && r.gc != null);
    const sum = (f: (r: DG.GrowthMetricsRow) => number | null) => gcRows.reduce((s, r) => s + (f(r) ?? 0), 0);
    const deleted = sum((r) => r.gc!.deletedLast);
    const lastRun = g.gc.lastRun ? dayjs(g.gc.lastRun) : null;
    this.setCard('gcDeleted', g.gc.lastRun ? formatCount(deleted) : '—',
      g.gc.enabled === false ? 'GC disabled' : lastRun ? `last run ${lastRun.format('MM-DD HH:mm')}` : 'no run yet',
      g.gc.enabled === false ? 'orange' : 'info');
    this.bindTopTooltip('gcDeleted', gcRows, (r) => r.gc!.deleted7d, (n) => formatCount(n),
      'Deleted in 7 days, by key');
    const pending = sum((r) => r.gc!.pending);
    const biggest = [...gcRows].sort((a, b) => (b.gc!.pending ?? 0) - (a.gc!.pending ?? 0))[0];
    this.setCard('gcPending', g.gc.pendingAt ? formatCount(pending) : '—',
      biggest?.gc?.pending ? biggest.key : g.gc.pendingAt ? 'nothing pending' : 'not counted yet', 'info');
    this.bindTopTooltip('gcPending', gcRows, (r) => r.gc!.pending, (n) => formatCount(n), 'Left to delete, by key');

    const holding: Array<[string, string, DG.GrowthMetricsAlert]> = [];
    for (const r of g.rows)
      for (const a of r.alerts.filter((a) => a.holding))
        holding.push([a.name, MetricsView.growthName(r), a]);
    const open = holding.some(([, , a]) => a.alertStatus === 'open' || a.alertStatus === 'acknowledged');
    const names = [...new Set(holding.map(([name]) => name))].join(', ');
    this.setCard('growthAlerts', `${holding.length}`, holding.length === 0 ? 'none holding' : names,
      holding.length === 0 ? 'green' : open ? 'red' : 'orange');
    const lines = holding.length === 0 ?
      [ui.divText('No metric alert holds.'),
        ui.divText('Settings > Alerts > Data growth sets them.', 'ua-metrics-tt-dim')] :
      [...holding.slice(0, TOOLTIP_TOP * 2).map(([name, row]) => ui.divText(`${name}: ${row}`)),
        ui.divText('Red: an alert is open. Browse > Platform > Observability > Alerts acts on them.', 'ua-metrics-tt-dim')];
    ui.tooltip.bind(this.card('growthAlerts').root, () => ui.div(lines, 'ua-metrics-tt'));
  }

  private setGrowthCard(id: string, s: DG.GrowthMetricsSeries | null, sub: string, parts: DG.GrowthMetricsRow[]): void {
    if (s == null || s.d7 == null) {
      this.setCard(id, '—', s == null ? 'collecting…' : `${sub} · collecting…`, 'info');
      return;
    }
    const pct = growthPct(s);
    const pctText = pct == null ? '' : `${pct > 0 ? '+' : ''}${pct.toFixed(0)}% · `;
    this.setCard(id, formatSigned(s.d7, formatBytes), `${pctText}${sub}`,
      pct == null ? 'info' : thresholdBand(pct, THRESH.growthPct, false));
    this.bindTopTooltip(id, parts, (r) => r.bytes?.d7 ?? null, (n) => formatSigned(n, formatBytes),
      'Grew most in 7 days');
  }

  /** Binds to card [id] the [TOOLTIP_TOP] rows of [rows] with the largest [value]. */
  private bindTopTooltip(id: string, rows: DG.GrowthMetricsRow[], value: (r: DG.GrowthMetricsRow) => number | null,
    format: (n: number) => string, title: string): void {
    const top = rows.filter((r) => (value(r) ?? 0) !== 0)
      .sort((a, b) => (value(b) ?? 0) - (value(a) ?? 0)).slice(0, TOOLTIP_TOP);
    ui.tooltip.bind(this.card(id).root, () => ui.div([
      ui.divText(title),
      top.length === 0 ? ui.divText('Nothing yet.', 'ua-metrics-tt-dim') : ui.div(top.map((r) => ui.divH([
        ui.divText(MetricsView.growthName(r), 'ua-metrics-tt-col-prefix'),
        ui.divText(format(value(r)!), 'ua-metrics-tt-col-num'),
      ])), 'ua-metrics-tt-table'),
    ], 'ua-metrics-tt'));
  }

  private static frame(rows: object[]): DG.DataFrame {
    return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.create();
  }

  private static routesFrame(routes: DG.ServerMetricsRoute[]): DG.DataFrame {
    return MetricsView.frame(routes.map((r) =>
      ({route: `${r.method} ${r.route}`, count: r.count, p50: r.p50, p95: r.p95, p99: r.p99, err_pct: r.errPct})));
  }

  private static statementsFrame(rows: DG.ServerMetricsStatement[]): DG.DataFrame {
    return MetricsView.frame(rows.map((s) =>
      ({query: s.query, calls: s.calls, total_ms: s.totalMs, mean_ms: s.meanMs, hit_pct: s.hitPct})));
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

  /** The failed requests of [r] in the window, each with its nearest error; the selected request's log docked below. */
  private async openRouteFailures(r: DG.ServerMetricsRoute): Promise<void> {
    const title = `${r.method} ${r.route} failures`;
    const df = await MetricsView.safeCall(
      () => queries.metricsRouteFailures(r.method, r.route, this.uaToolbox.getFilter().date!), 'MetricsRouteFailures');
    if (!df) {
      grok.shell.warning(`${title} unavailable.`);
      return;
    }
    df.name = title;
    const tv = grok.shell.addTableView(df);
    const log = DG.Viewer.grid(DG.DataFrame.create(0), MetricsView.gridOptions());
    tv.dockManager.dock(log, DG.DOCK_TYPE.DOWN, null, 'Request log', 0.4);
    const showLog = async (): Promise<void> => {
      const i = df.currentRowIdx;
      const session = i < 0 ? null : df.get('session_id', i) as string | null;
      if (!session)
        return;
      const logDf = await MetricsView.safeCall(
        () => queries.metricsRequestLog(session, (df.get('time', i) as dayjs.Dayjs).toISOString()), 'MetricsRequestLog');
      if (logDf)
        log.dataFrame = logDf;
    };
    tv.subs.push(df.onCurrentRowChanged.subscribe(() => showLog()));
    if (df.rowCount > 0)
      df.currentRowIdx = 0;
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

  private async collectSnapshotAttachments(): Promise<DG.FileInfo[]> {
    const date = this.uaToolbox.getFilter().date!;
    const metrics = MetricsView.safeCall(() => this.fetchMetrics(EMAIL_LIMIT), 'getMetrics');
    const storage = grok.dapi.info.getStorageStats();

    type Src = {name: string, fn: () => Promise<DG.DataFrame>};
    const sources: Src[] = [
      {name: 'table_health.csv',          fn: () => queries.metricsTableHealth(EMAIL_LIMIT)},
      {name: 'table_health_summary.csv',  fn: () => queries.metricsTableHealthSummary(EMAIL_LIMIT)},
      {name: 'largest_tables.csv',        fn: () => queries.metricsLargestTables(20)},
      {name: 'cache_miss_tables.csv',     fn: () => queries.metricsCacheMissTables()},
      {name: 'connections_offenders.csv', fn: () => queries.metricsConnectionsOffenders(20, 60, 30)},
      {name: 'errors.csv',                fn: () => queries.metricsErrorsCount(date)},
      {name: 'sessions.csv',              fn: () => queries.metricsSessionsCount(date)},
      {name: 'storage.csv',               fn: () => storage.then(({topPrefixes, ...totals}) => DG.DataFrame.fromObjects([totals])!)},
      {name: 'storage_prefixes.csv',      fn: () => storage.then((s) => DG.DataFrame.fromObjects(s.topPrefixes)!)},
      {name: 'disk.csv',                  fn: async () => DG.DataFrame.fromObjects([JSON.parse(await grok.functions.call('DiskStats'))])!},
    ];

    const [m, results] = await Promise.all([metrics, Promise.all(sources.map((s) => MetricsView.safeCall(s.fn, s.name)))]);

    const files: DG.FileInfo[] = [];
    for (let i = 0; i < sources.length; i++) {
      const r = results[i];
      if (r == null)
        continue;
      files.push(DG.FileInfo.fromString(sources[i].name, r.toCsv()));
    }
    if (m) {
      const {window: w, http, queue, database: db} = m;
      const frames: {[name: string]: DG.DataFrame} = {
        'db_summary.csv': MetricsView.frame([{size_bytes: db.sizeBytes, cache_hit_pct: db.cacheHitPct, stats_reset: db.statsReset}]),
        'connections.csv': MetricsView.frame([db.connections]),
        'queue.csv': MetricsView.frame([queue]),
        'http_routes.csv': MetricsView.routesFrame(http.routes),
        'latency.csv': MetricsView.frame([{...http.now, p95_prev: http.previous.p95, count_prev: http.previous.count,
          window_start: w.start, window_end: w.end, prev_window_start: w.previousStart}]),
      };
      if (db.statements.available) {
        frames['slowest_queries.csv'] = MetricsView.statementsFrame(db.statements.slowest);
        frames['most_called_queries.csv'] = MetricsView.statementsFrame(db.statements.mostCalled);
        frames['worst_cache_hit.csv'] = MetricsView.statementsFrame(db.statements.worstCacheHit);
      }
      for (const name of Object.keys(frames))
        files.push(DG.FileInfo.fromString(name, frames[name].toCsv()));
    }
    return files;
  }

}
