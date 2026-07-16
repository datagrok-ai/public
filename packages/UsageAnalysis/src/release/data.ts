import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {BehaviorSubject, Subject} from 'rxjs';

import {queries} from '../package-api';

export const DEV_HOST = 'dev.datagrok.ai';
export const JIRA_BROWSE = 'https://reddata.atlassian.net/browse/';
export const JENKINS_TEST_JOB = 'https://ops.datagrok.ai/job/test-build/';
const SLOW_FACTOR = 1.25;
const MIN_SLOW_MS = 1000; // ignore sub-second tests — a +% on a few-ms test is noise, not a regression

/** Instances the dashboard can inspect; each maps to a `test_runs.instance` host (see release_tests.sql). */
export const ENV_CHOICES = ['dev', 'release', 'public', 'release-ec2'];
export const DEFAULT_ENV = 'dev';

/** Dashboard-wide state shared by every tab. The ribbon env picker (ReleaseHandler) writes `env`;
 * env-scoped tabs (Overview, Tests) read `env.value` and re-fetch when it changes. The ribbon Refresh
 * button fires `refresh`; every tab re-fetches on it. */
export class ReleaseContext {
  readonly env: BehaviorSubject<string>;
  readonly refresh = new Subject<void>();
  constructor(env: string = DEFAULT_ENV) {
    this.env = new BehaviorSubject(env);
  }
}

// Version-bound test muting. A test is muted for a specific release version; the mute list lives in its
// own sticky-meta schema (semtype=autotest) so it never touches the shared 'Autotests' schema.
export const RELEASE_MUTES_SCHEMA = 'ReleaseMutes';
export const MUTED_VERSIONS = 'mutedVersions';           // failing-alert mute (version-bound)
export const STALE_MUTED_VERSIONS = 'staleMutedVersions'; // stale-alert mute (version-bound, independent)
export const MUTE_ON = '🔕';  // muted for this release
export const MUTE_OFF = '🔔'; // active (click to mute)
export const STALE_DAYS = 7;  // "not run for a week or more"
let _mutesSchemaPromise: Promise<any> | null = null;

/** The release-mute sticky-meta schema (fields: mutedVersions, staleMutedVersions), created/migrated on
 * first use. Promise-cached so the Overview and Tests tabs (both fetch on load) can't race into two
 * createSchema calls; recovers if one loses the race, and adds staleMutedVersions to a pre-existing schema. */
export function getReleaseMutesSchema(): Promise<any> {
  _mutesSchemaPromise ??= (async () => {
    let schemas = await grok.dapi.stickyMeta.getSchemas();
    let schema = schemas.find((s) => s.name === RELEASE_MUTES_SCHEMA);
    if (schema)
      return await ensureStaleField(schema);
    try {
      // Entity-type name must be globally unique across schemas (the shared 'Autotests' schema owns
      // 'autotest'), so use a distinct name while still matching the same autotest-semtype column.
      return await grok.dapi.stickyMeta.createSchema(RELEASE_MUTES_SCHEMA,
        [{name: 'releaseMuteAutotest', matchBy: 'semtype=autotest'}],
        [{name: MUTED_VERSIONS, type: DG.TYPE.STRING}, {name: STALE_MUTED_VERSIONS, type: DG.TYPE.STRING}]);
    } catch (e) {
      // A concurrent caller created it first — re-fetch before giving up.
      schemas = await grok.dapi.stickyMeta.getSchemas();
      schema = schemas.find((s) => s.name === RELEASE_MUTES_SCHEMA);
      if (schema)
        return await ensureStaleField(schema);
      _mutesSchemaPromise = null; // allow a later retry on a genuine failure
      throw e;
    }
  })();
  return _mutesSchemaPromise;
}

/** Adds the staleMutedVersions property to a pre-existing ReleaseMutes schema (migration). */
async function ensureStaleField(schema: any): Promise<any> {
  if (schema.properties.some((p: any) => p.name === STALE_MUTED_VERSIONS))
    return schema;
  try {
    schema.properties = [...schema.properties, DG.EntityProperty.create(STALE_MUTED_VERSIONS, DG.TYPE.STRING)];
    await grok.dapi.stickyMeta.saveSchema(schema);
  } catch (_) { /* field may already exist from a concurrent migration — safe to ignore */ }
  return schema;
}

/** Extracts the release version (e.g. '1.28.0') from a build name like '2026-07-15-1.28.0.BE-3'. */
export function versionFromBuild(build: string): string {
  const m = /(\d+\.\d+\.\d+)/.exec(build ?? '');
  return m ? m[1] : '';
}

/** Parses a stored 'a;b;c' mute list into a trimmed, non-empty version array. */
export function parseMutedVersions(raw: string | null | undefined): string[] {
  return (raw ?? '').split(';').map((s) => s.trim()).filter((s) => s);
}

/** True when `version` is present in the stored mute list. */
export function isMutedForVersion(raw: string | null | undefined, version: string): boolean {
  return !!version && parseMutedVersions(raw).includes(version);
}

// Grid status-cell colors (ARGB — grid cells take numeric colors, not design tokens; see metrics.ts/errors.ts).
export const STATUS_BACK: Record<string, number> = {passed: 0xFFE6F4EA, failed: 0xFFFDE7E7, skipped: 0xFFF3F3F3};
export const NOT_RUN_BACK = 0xFFFAFAFA;
export const COLOR_FAIL_TEXT = 0xFF7B2D24;
// Dimmed palette for a carried-forward result (test didn't run in the latest build → show the last known one, muted).
export const DIM_STATUS_BACK: Record<string, number> = {passed: 0xFFF0F7F1, failed: 0xFFF7ECEC, skipped: 0xFFF6F6F6};
export const COLOR_DIM_TEXT = 0xFF9AA0A6;

/** Resolves once the element has a non-zero width (or after ~1s). Sparkline grid columns throw a
 * layout "Wrong range" error if added before the grid is laid out — wait for a real width first. */
export async function waitForWidth(el: HTMLElement, tries = 60): Promise<void> {
  for (let i = 0; i < tries && el.clientWidth === 0; i++)
    await new Promise((r) => requestAnimationFrame(r));
}

/** Tints per-build/per-batch status cells (passed/failed/skipped/did-not-run) in a results grid. */
export function colorStatusCell(gc: DG.GridCell, statusCols: string[]): void {
  if (!gc.isTableCell || !statusCols.includes(gc.gridColumn.name))
    return;
  const v = gc.cell.value as string;
  gc.style.backColor = STATUS_BACK[v] ?? NOT_RUN_BACK;
  if (v === 'failed')
    gc.style.textColor = COLOR_FAIL_TEXT;
}

/** A '+' icon that opens the current dataframe as a table view in the workspace. */
export function openInWorkspaceIcon(getDf: () => DG.DataFrame | null): HTMLElement {
  const icon = ui.iconFA('plus', () => {
    const df = getDf();
    if (df && df.rowCount)
      grok.shell.addTableView(df);
  }, 'Open as a table in the workspace');
  icon.classList.add('ua-release-open');
  return icon;
}

/** The Release dashboard is dev-only. True when the current server is dev.datagrok.ai. */
export function isDevHost(): boolean {
  const root = (grok.dapi.root ?? '').toLowerCase();
  if (root.includes(DEV_HOST))
    return true;
  return typeof window !== 'undefined' && window.location.hostname === DEV_HOST;
}

/** Default next-release version to prefill the tickets filter (current client build). */
export function defaultNextVersion(): string {
  return grok.shell.build?.client?.version ?? '';
}

export function median(xs: number[]): number {
  if (xs.length === 0)
    return 0;
  const s = [...xs].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

/** One-row-per-test pivot of ReleaseTests plus the derived per-build column names. */
export interface ReleasePivot {
  df: DG.DataFrame;
  buildIndexes: number[]; // ascending: first = oldest, last = latest
  latest: number;
  statusCols: string[];
  durationCols: string[];
  successCols: string[];
  resultCols: string[];
  version: string; // release version the latest build belongs to (for version-bound muting)
}

export interface TestAlerts {
  total: number;
  passed: number;
  failed: number;
  skipped: number;
  notRun: number;
  flaky: number;
  slower: number;
  passRate: number; // 0..100 over the latest build: passed / (passed + failed), excluding skipped & did-not-run
  stale: number; // tests not run for >= STALE_DAYS (unmuted)
  failingByPkg: Record<string, string[]>;
  slowerByPkg: Record<string, string[]>;
  flakyByPkg: Record<string, string[]>;
  staleByPkg: Record<string, string[]>;
}

function pushByPkg(map: Record<string, string[]>, pkg: string, test: string): void {
  (map[pkg] ??= []).push(test);
}

function countGroups(map: Record<string, string[]>): number {
  return Object.values(map).reduce((n, xs) => n + xs.length, 0);
}

export interface StressRegression {
  regressed: boolean;
  detail: string;
}

function renameCol(df: DG.DataFrame, from: string, to: string): void {
  const c = df.columns.byName(from);
  if (c) {
    c.name = to;
    c.setTag('friendlyName', to);
  }
}

/** Runs ReleaseTests and pivots it client-side into one row per test with per-build columns. */
export async function fetchReleaseTests(instanceFilter: string, lastBuildsNum: number): Promise<ReleasePivot | null> {
  const raw: DG.DataFrame = await queries.releaseTests(instanceFilter, lastBuildsNum);
  if (raw.rowCount === 0)
    return null;

  const biCol = raw.getCol('build_index');
  const buildIndexes = Array.from(new Set(Array.from({length: biCol.length}, (_, i) => biCol.get(i))))
    .filter((v) => v != null).map(Number).sort((a, b) => a - b);
  const latest = buildIndexes[buildIndexes.length - 1];

  const pivot = raw.groupBy(['package', 'test', 'owner'])
    .pivot('build_index')
    .first('status')
    .avg('duration')
    .first('flaking')
    .first('result')
    .first('build_commit')
    .aggregate();

  // Per-build success (1 pass / 0 fail / gap for skipped & did-not-run) for the success sparkline.
  for (const bi of buildIndexes) {
    const sc = pivot.getCol(`${bi} first(status)`);
    pivot.columns.addNewFloat(`${bi} ok`).init((i) => {
      const s = sc.get(i);
      return s === 'passed' ? 1 : s === 'failed' ? 0 : null;
    });
  }

  const latestCommit = pivot.columns.addNewString('latest_commit');
  latestCommit.init((i) => pivot.get(`${latest} first(build_commit)`, i) ?? '');
  const latestResult = pivot.columns.addNewString('latest_result');

  for (const bi of buildIndexes) {
    renameCol(pivot, `${bi} first(status)`, `${bi}`);
    renameCol(pivot, `${bi} avg(duration)`, `${bi} ms`);
    renameCol(pivot, `${bi} first(result)`, `${bi} result`);
    renameCol(pivot, `${bi} first(flaking)`, `${bi} flaking`);
    pivot.columns.remove(`${bi} first(build_commit)`);
  }
  latestResult.init((i) => pivot.get(`${latest} result`, i) ?? '');

  const statusCols = buildIndexes.map((bi) => `${bi}`);
  const durationCols = buildIndexes.map((bi) => `${bi} ms`);
  const successCols = buildIndexes.map((bi) => `${bi} ok`);
  const resultCols = buildIndexes.map((bi) => `${bi} result`);

  // Release version + per-test baselines (last-run date, previous-release avg duration), read from the raw
  // result up front so the derived flags below (slower, stale) can use them.
  const buildCol = raw.getCol('build');
  let latestBuildName = '';
  for (let i = 0; i < raw.rowCount; i++)
    if (Number(biCol.get(i)) === latest) { latestBuildName = buildCol.get(i) ?? ''; break; }
  const version = versionFromBuild(latestBuildName);
  const nowMs = Date.now();
  const rawTestCol = raw.getCol('test');
  const lastRunByTest = new Map<string, number>();
  const rawLastRun = raw.columns.byName('last_run');
  if (rawLastRun) {
    for (let i = 0; i < raw.rowCount; i++) {
      const t = rawTestCol.get(i);
      if (!lastRunByTest.has(t)) {
        const v: any = rawLastRun.get(i);
        lastRunByTest.set(t, v == null ? NaN : (v.valueOf ? Number(v.valueOf()) : Number(v)));
      }
    }
  }
  const prevReleaseByTest = new Map<string, number>();
  const rawPrev = raw.columns.byName('prev_release_ms');
  if (rawPrev) {
    for (let i = 0; i < raw.rowCount; i++) {
      const t = rawTestCol.get(i);
      if (!prevReleaseByTest.has(t)) { const v: any = rawPrev.get(i); prevReleaseByTest.set(t, v == null ? NaN : Number(v)); }
    }
  }

  // Carry the last known result forward when the test didn't run in the latest build.
  const descBuilds = [...buildIndexes].reverse(); // latest → oldest
  pivot.columns.addNewString('effective_status').init((i) => {
    for (const bi of descBuilds) {
      const s = pivot.get(`${bi}`, i);
      if (s && s !== 'did not run')
        return s;
    }
    return '';
  });
  pivot.columns.addNewBool('carried_forward').init((i) => {
    const latestS = pivot.get(`${latest}`, i);
    return (latestS == null || latestS === 'did not run') && pivot.get('effective_status', i) !== '';
  });

  pivot.columns.addNewBool('failing').init((i) => pivot.get('effective_status', i) === 'failed');
  pivot.columns.addNewBool('stable').init((i) => statusCols.every((c) => pivot.get(c, i) === 'passed'));
  pivot.columns.addNewBool('flaky').init((i) => {
    const statuses = statusCols.map((c) => pivot.get(c, i));
    const flakingFlag = buildIndexes.some((bi) => pivot.get(`${bi} flaking`, i) === true);
    return flakingFlag || (statuses.includes('passed') && statuses.includes('failed'));
  });
  // Slow-test baselines: bleeding-edge recent average (prior window passed runs) and the previous-release
  // average (from the query). A test is "slower" if the latest run exceeds SLOW_FACTOR x either baseline.
  const baselineTestCol = pivot.getCol('test');
  pivot.columns.addNewInt('prev_release_ms').init((i) => {
    const v = prevReleaseByTest.get(baselineTestCol.get(i));
    return (v == null || isNaN(v)) ? 0 : v;
  });
  pivot.columns.addNewInt('be_avg_ms').init((i) => {
    const prior = buildIndexes.slice(0, -1)
      .filter((bi) => pivot.get(`${bi}`, i) === 'passed')
      .map((bi) => Number(pivot.get(`${bi} ms`, i))).filter((v) => v && !isNaN(v));
    return prior.length ? Math.round(prior.reduce((a, b) => a + b, 0) / prior.length) : 0;
  });
  pivot.columns.addNewBool('slower').init((i) => {
    if (pivot.get(`${latest}`, i) !== 'passed')
      return false;
    const latestMs = Number(pivot.get(`${latest} ms`, i));
    if (!latestMs || isNaN(latestMs) || latestMs < MIN_SLOW_MS)
      return false;
    const be = pivot.get('be_avg_ms', i);
    const prev = pivot.get('prev_release_ms', i);
    return (be > 0 && latestMs > be * SLOW_FACTOR) || (prev > 0 && latestMs > prev * SLOW_FACTOR);
  });
  pivot.columns.addNewString('slow_detail').init((i) => {
    if (pivot.get('slower', i) !== true)
      return '';
    const latestMs = Number(pivot.get(`${latest} ms`, i));
    const be = pivot.get('be_avg_ms', i);
    const prev = pivot.get('prev_release_ms', i);
    const pct = (base: number) => { const p = Math.round((latestMs / base - 1) * 100); return `${p >= 0 ? '+' : ''}${p}%`; };
    const parts: string[] = [];
    if (be > 0) parts.push(`${pct(be)} vs recent avg ${be}ms`);
    if (prev > 0) parts.push(`${pct(prev)} vs prev release ${prev}ms`);
    return `${latestMs}ms — ${parts.join(', ')}`;
  });

  const schemas = await grok.dapi.stickyMeta.getSchemas();
  const schema = schemas.find((s) => s.name === 'Autotests');
  pivot.getCol('test').semType = 'autotest';
  pivot.getCol('owner').semType = 'User';
  pivot.getCol('owner').setTag('cell.renderer', 'User');
  if (schema) {
    const meta = await grok.dapi.stickyMeta.getAllValues(schema, pivot.getCol('test'));
    for (const name of ['ignore?', 'ignoreReason', 'lastResolved']) {
      const c = meta.col(name);
      if (c)
        pivot.columns.add(c);
    }
  }
  // Version-bound mutes (own schema). 'mutedForVersion' drives the grid icon; 'muted' = global ignore OR
  // muted-for-this-version and gates needs_attention and every failing/flaky count.
  const mutesSchema = await getReleaseMutesSchema();
  const mutes = await grok.dapi.stickyMeta.getAllValues(mutesSchema, pivot.getCol('test'));
  const mvCol = mutes.col(MUTED_VERSIONS);
  if (mvCol)
    pivot.columns.add(mvCol);
  const mutedVersionsCol = pivot.columns.byName(MUTED_VERSIONS);
  const ignoreCol = pivot.columns.byName('ignore?');
  pivot.columns.addNewBool('mutedForVersion')
    .init((i) => isMutedForVersion(mutedVersionsCol?.get(i), version));
  pivot.columns.addNewBool('muted')
    .init((i) => (ignoreCol != null && ignoreCol.get(i) === true) || pivot.get('mutedForVersion', i) === true);
  // Mute toggle column: the glyph is the cell's own value so it always stays aligned with its row
  // (reading another column by grid-row index misaligns once the grid sorts). tests.ts toggles it on click.
  pivot.columns.addNewString('mute').init((i) => pivot.get('mutedForVersion', i) === true ? MUTE_ON : MUTE_OFF);
  pivot.columns.addNewBool('needs_attention')
    .init((i) => pivot.get('failing', i) && !pivot.get('muted', i));

  // Staleness: tests not run for >= STALE_DAYS. Its mute is separate (staleMutedVersions) from the failing mute.
  const testColP = pivot.getCol('test');
  const smvCol = mutes.col(STALE_MUTED_VERSIONS);
  if (smvCol)
    pivot.columns.add(smvCol);
  const staleMutedVersionsCol = pivot.columns.byName(STALE_MUTED_VERSIONS);
  pivot.columns.addNewInt('stale_days').init((i) => {
    const lr = lastRunByTest.get(testColP.get(i));
    return (lr == null || isNaN(lr)) ? 0 : Math.floor((nowMs - lr) / 86400000);
  });
  pivot.columns.addNewBool('stale').init((i) => pivot.get('stale_days', i) >= STALE_DAYS);
  pivot.columns.addNewBool('staleMutedForVersion').init((i) => isMutedForVersion(staleMutedVersionsCol?.get(i), version));
  pivot.columns.addNewBool('stale_alert')
    .init((i) => pivot.get('stale', i) === true && pivot.get('staleMutedForVersion', i) !== true);
  // Toggle glyph only for stale tests (blank otherwise).
  pivot.columns.addNewString('staleMute')
    .init((i) => pivot.get('stale', i) !== true ? '' : (pivot.get('staleMutedForVersion', i) === true ? MUTE_ON : MUTE_OFF));

  pivot.name = 'Release tests';
  return {df: pivot, buildIndexes, latest, statusCols, durationCols, successCols, resultCols, version};
}

/** Rolls the pivot up into the counts and package-grouped alert lists the Overview needs. */
export function computeTestAlerts(p: ReleasePivot): TestAlerts {
  const df = p.df;
  const a: TestAlerts = {total: 0, passed: 0, failed: 0, skipped: 0, notRun: 0, flaky: 0, slower: 0,
    passRate: 0, stale: 0, failingByPkg: {}, slowerByPkg: {}, flakyByPkg: {}, staleByPkg: {}};
  const testCol = df.getCol('test');
  const pkgCol = df.getCol('package');
  const mutedCol = df.columns.byName('muted');
  let rawFailed = 0;
  for (let i = 0; i < df.rowCount; i++) {
    const status = df.get('effective_status', i) || 'did not run';
    const muted = !!(mutedCol && mutedCol.get(i) === true);
    const pkg = pkgCol.get(i) ?? '';
    const test = testCol.get(i);
    a.total++;
    if (status === 'passed') a.passed++;
    else if (status === 'failed') rawFailed++;
    else if (status === 'skipped') a.skipped++;
    else a.notRun++;
    if (status === 'failed' && !muted)
      pushByPkg(a.failingByPkg, pkg, test);
    if (df.get('flaky', i) === true && !muted)
      pushByPkg(a.flakyByPkg, pkg, test);
    if (df.get('slower', i) === true && !muted)
      pushByPkg(a.slowerByPkg, pkg, `${test}: ${df.get('slow_detail', i) ?? ''}`);
    if (df.get('stale_alert', i) === true)
      pushByPkg(a.staleByPkg, pkg, test);
  }
  a.failed = countGroups(a.failingByPkg);
  a.flaky = countGroups(a.flakyByPkg);
  a.slower = countGroups(a.slowerByPkg);
  a.stale = countGroups(a.staleByPkg);
  // Success rate excludes skipped and did-not-run: passed / (passed + failed).
  // floor so it never rounds up to 100% while any test is failing.
  const ran = a.passed + rawFailed;
  a.passRate = ran > 0 ? Math.floor((a.passed / ran) * 100) : 0;
  return a;
}

/** One-row-per-test pivot of the Manual Tests query (Test Track batches). batch_index 1 = newest. */
export interface ManualPivot {
  df: DG.DataFrame;
  chronological: number[]; // batch indexes ordered oldest -> newest (for left-to-right sparklines)
  latest: number; // batch_index of the newest batch (= 1)
  statusCols: string[];
  successCols: string[];
  latestBatchName: string;
  passed: number;
  failed: number;
  total: number;
}

export async function fetchReleaseManual(lastBatches: number): Promise<ManualPivot | null> {
  const raw: DG.DataFrame = await queries.releaseManualTests(lastBatches);
  if (raw.rowCount === 0)
    return null;
  const biCol = raw.getCol('batch_index');
  const idxs = Array.from(new Set(Array.from({length: biCol.length}, (_, i) => biCol.get(i))))
    .filter((v) => v != null).map(Number).sort((a, b) => a - b);
  const latest = idxs[0]; // ROW_NUMBER ... ORDER BY latest_batch_run DESC => 1 is newest
  const chronological = [...idxs].reverse();
  let latestBatchName = '';
  for (let i = 0; i < raw.rowCount; i++)
    if (Number(biCol.get(i)) === latest) { latestBatchName = raw.get('batch_name', i) ?? ''; break; }

  const pivot = raw.groupBy(['test']).pivot('batch_index').first('status').first('result').aggregate();
  for (const bi of idxs) {
    const sc = pivot.getCol(`${bi} first(status)`);
    pivot.columns.addNewFloat(`${bi} ok`).init((i) => { const s = sc.get(i); return s === 'passed' ? 1 : s === 'failed' ? 0 : null; });
  }
  const latestResult = pivot.columns.addNewString('latest_result');
  for (const bi of idxs) {
    renameCol(pivot, `${bi} first(status)`, `${bi}`);
    renameCol(pivot, `${bi} first(result)`, `${bi} result`);
  }
  latestResult.init((i) => pivot.get(`${latest} result`, i) ?? '');
  pivot.columns.addNewBool('failing').init((i) => pivot.get(`${latest}`, i) === 'failed');
  pivot.getCol('test').semType = 'autotest';
  pivot.name = 'Manual tests';

  let passed = 0; let failed = 0;
  for (let i = 0; i < pivot.rowCount; i++) {
    const s = pivot.get(`${latest}`, i);
    if (s === 'passed') passed++;
    else if (s === 'failed') failed++;
  }
  return {df: pivot, chronological, latest, statusCols: chronological.map((bi) => `${bi}`),
    successCols: chronological.map((bi) => `${bi} ok`), latestBatchName, passed, failed, total: pivot.rowCount};
}

/** Compares the latest build's worst-case p95 to the median of prior builds (from StressTestsSummary). */
export function stressRegression(df: DG.DataFrame, factor = 1.2): StressRegression {
  const n = df.rowCount;
  if (n === 0)
    return {regressed: false, detail: 'No stress runs reported.'};
  const buildCol = df.getCol('build');
  const startedCol = df.columns.byName('started') ?? df.columns.byName('build_date');
  const p95Col = df.columns.byName('p95_ms') ?? df.columns.byName('median_ms');
  if (!p95Col)
    return {regressed: false, detail: 'No duration metric available.'};

  const byBuild = new Map<string, {started: number, p95: number}>();
  for (let i = 0; i < n; i++) {
    const b = buildCol.get(i);
    const p = Number(p95Col.get(i) ?? 0);
    const raw = startedCol ? startedCol.get(i) : i;
    const started = raw && (raw as any).valueOf ? (raw as any).valueOf() : Number(raw) || 0;
    const cur = byBuild.get(b);
    byBuild.set(b, {started: Math.max(started, cur?.started ?? 0), p95: Math.max(p, cur?.p95 ?? 0)});
  }
  const builds = [...byBuild.entries()].map(([build, v]) => ({build, ...v})).sort((x, y) => y.started - x.started);
  if (builds.length < 3)
    return {regressed: false, detail: 'Not enough stress history to compare.'};
  const latest = builds[0];
  const priorMed = median(builds.slice(1).map((b) => b.p95));
  const pct = priorMed > 0 ? Math.round((latest.p95 / priorMed - 1) * 100) : 0;
  const regressed = priorMed > 0 && latest.p95 > priorMed * factor;
  return {regressed, detail: regressed ?
    `Latest p95 ${Math.round(latest.p95)}ms is ${pct}% above prior median ${Math.round(priorMed)}ms` :
    `Latest p95 within range (${pct >= 0 ? '+' : ''}${pct}% vs prior median)`};
}
