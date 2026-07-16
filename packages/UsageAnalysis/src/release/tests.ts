import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {fetchReleaseTests, computeTestAlerts, ReleasePivot, ReleaseContext, JENKINS_TEST_JOB, colorStatusCell,
  COLOR_FAIL_TEXT, DIM_STATUS_BACK, COLOR_DIM_TEXT, NOT_RUN_BACK, waitForWidth, openInWorkspaceIcon,
  getReleaseMutesSchema, MUTED_VERSIONS, STALE_MUTED_VERSIONS, MUTE_ON, MUTE_OFF,
  parseMutedVersions, isMutedForVersion} from './data';

export class TestsView extends UaView {
  private lastBuildsInput!: DG.InputBase;
  private summaryHost!: HTMLElement;
  private gridHost!: HTMLElement;
  private schema: any = null;
  private pivot: ReleasePivot | null = null;
  private summaryDf: DG.DataFrame | null = null;
  private detailDf: DG.DataFrame | null = null;
  private grid: DG.Grid | null = null;
  private persisting = false;
  private pendingPersist = false;
  private persistingMutes = false;
  private pendingMutes = false;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tests';
    this.ctx.env.subscribe(() => { if (this.initialized) this.refresh(); });
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics ua-metrics-fill';
    this.lastBuildsInput = ui.input.int('Builds', {value: 5, min: 2, max: 20, onValueChanged: () => this.refresh()});
    this.summaryHost = ui.div([], 'ua-metrics-grid-host');
    this.gridHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([
      this.lastBuildsInput.root,
      ui.div([], {style: {flex: '1'}}), refreshBtn,
    ], 'ua-metrics-header');

    const summaryPanel = ui.div([
      ui.divH([ui.divText('Results by package (latest build)', 'ua-metrics-panel-title'),
        openInWorkspaceIcon(() => this.summaryDf)], 'ua-metrics-panel-header'),
      this.summaryHost,
    ], 'ua-metrics-panel');
    const gridPanel = ui.div([
      ui.divH([ui.divText('Tests (last builds — status, success & duration sparklines)', 'ua-metrics-panel-title'),
        openInWorkspaceIcon(() => this.detailDf)], 'ua-metrics-panel-header'),
      this.gridHost,
    ], 'ua-metrics-panel ua-metrics-panel-grow');

    this.root.append(ui.div([header, summaryPanel, gridPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    this.summaryHost.innerHTML = '';
    this.summaryHost.append(ui.loader());
    this.gridHost.innerHTML = '';
    try {
      this.pivot = await fetchReleaseTests(this.ctx.env.value, this.lastBuildsInput.value ?? 5);
      if (!this.pivot) {
        this.summaryHost.innerHTML = '';
        this.summaryHost.append(ui.divText('No package test runs reported for this instance.', 'ua-metrics-degraded'));
        return;
      }
      this.buildSummary(this.pivot);
      this.buildGrid(this.pivot);
    } catch (e) {
      this.summaryHost.innerHTML = '';
      this.summaryHost.append(ui.divText(`ReleaseTests failed: ${e}`, 'ua-metrics-degraded'));
    }
  }

  private buildSummary(p: ReleasePivot): void {
    const a = computeTestAlerts(p);
    const df = p.df;
    const packageCol = df.getCol('package');
    const mutedCol = df.columns.byName('muted');
    const perPkg = new Map<string, {passed: number, failed: number, skipped: number, notRun: number}>();
    for (let i = 0; i < df.rowCount; i++) {
      const pkg = packageCol.get(i) ?? '';
      const s = df.get('effective_status', i) || 'did not run';
      const muted = !!(mutedCol && mutedCol.get(i) === true);
      const row = perPkg.get(pkg) ?? {passed: 0, failed: 0, skipped: 0, notRun: 0};
      if (s === 'passed') row.passed++;
      else if (s === 'failed' && !muted) row.failed++;
      else if (s === 'skipped') row.skipped++;
      else if (s !== 'failed') row.notRun++;
      perPkg.set(pkg, row);
    }
    const pkgs = [...perPkg.entries()].sort((x, y) => y[1].failed - x[1].failed || x[0].localeCompare(y[0]));
    const summary = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('package', pkgs.map((e) => e[0])),
      DG.Column.fromInt32Array('passed', Int32Array.from(pkgs.map((e) => e[1].passed))),
      DG.Column.fromInt32Array('failed', Int32Array.from(pkgs.map((e) => e[1].failed))),
      DG.Column.fromInt32Array('skipped', Int32Array.from(pkgs.map((e) => e[1].skipped))),
      DG.Column.fromInt32Array('not run', Int32Array.from(pkgs.map((e) => e[1].notRun))),
      DG.Column.fromInt32Array('pass rate', Int32Array.from(pkgs.map((e) => {
        const ran = e[1].passed + e[1].failed; // success rate excludes skipped & not-run
        return ran > 0 ? Math.floor((e[1].passed / ran) * 100) : 0;
      }))),
    ]);
    summary.name = `By package — ${a.passed}/${a.total} passed, ${a.failed} failing`;
    this.summaryDf = summary;
    const grid = DG.Viewer.grid(summary, {showColumnGridlines: false, allowBlockSelection: false});
    grid.onCellPrepare((gc) => {
      if (!gc.isTableCell) return;
      if (gc.gridColumn.name === 'failed' && (gc.cell.value ?? 0) > 0)
        gc.style.textColor = COLOR_FAIL_TEXT;
      if (gc.gridColumn.name === 'pass rate') {
        const v = gc.cell.value ?? 0;
        gc.style.textColor = v >= 95 ? 0xFF2CA02C : v >= 80 ? 0xFFFFA500 : 0xFFFF5A00;
      }
    });
    this.summaryHost.innerHTML = '';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.summaryHost.append(grid.root);
  }

  private buildGrid(p: ReleasePivot): void {
    const df = p.df;
    this.detailDf = df;
    const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false, showCurrentCellOutline: false});
    this.grid = grid;

    this.gridHost.innerHTML = '';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.gridHost.append(grid.root);

    grid.onCellPrepare((gc) => {
      colorStatusCell(gc, p.statusCols);
      // Compact per-build status cells: '+' passed / '−' failed / '·' skipped / blank did-not-run.
      if (gc.isTableCell && p.statusCols.includes(gc.gridColumn.name)) {
        const v = gc.cell.value as string;
        gc.customText = v === 'passed' ? '+' : v === 'failed' ? '−' : v === 'skipped' ? '·' : '';
      }
      // Latest build didn't run this test → show the last known result in that cell, dimmed.
      if (gc.isTableCell && gc.gridColumn.name === `${p.latest}` && df.get('carried_forward', gc.cell.rowIndex) === true) {
        const eff = df.get('effective_status', gc.cell.rowIndex) as string;
        gc.customText = eff === 'passed' ? '+' : eff === 'failed' ? '−' : eff === 'skipped' ? '·' : '';
        gc.style.backColor = DIM_STATUS_BACK[eff] ?? NOT_RUN_BACK;
        gc.style.textColor = COLOR_DIM_TEXT;
      }
      // The 'mute' column value already holds the 🔔/🔕 glyph (kept aligned per row); click toggles it.
    });
    grid.onCellClick.subscribe((gc) => {
      if (!gc.isTableCell)
        return;
      if (gc.gridColumn.name === 'mute')
        this.toggleMute(df, gc.cell.rowIndex, p);
      else if (gc.gridColumn.name === 'staleMute' && df.get('stale', gc.cell.rowIndex) === true)
        this.toggleStaleMute(df, gc.cell.rowIndex, p);
    });
    df.onCurrentRowChanged.subscribe(() => {
      if (df.currentRowIdx >= 0)
        grok.shell.o = this.buildDrilldown(df, df.currentRowIdx, p);
    });
    // Persist mute edits (ignore? / ignoreReason) back to the Autotests sticky-meta schema.
    df.onValuesChanged.subscribe(() => this.persistMute(df));

    // Sparklines + ordering must wait until the grid has a real width, else the grid's layout throws
    // a "Wrong range" error. Degrade gracefully (raw columns) if it still can't add them.
    waitForWidth(grid.root).then(() => {
      try {
        const passing = grid.columns.add({gridColumnName: 'passing', cellType: 'sparkline'});
        passing.settings = {columnNames: p.successCols};
        const duration = grid.columns.add({gridColumnName: 'duration', cellType: 'sparkline'});
        duration.settings = {columnNames: p.durationCols};
      } catch (_) { /* grid not ready for sparklines — leave the raw per-build columns visible */ }

      // Show only these, in this order: test → package → passing history (sparkline + per-build status
      // cells) → duration sparkline → flags. `setVisible` whitelists + orders in one deterministic call.
      const visible = ['mute', 'test', 'package', 'passing', ...p.statusCols, 'duration',
        'needs_attention', 'stale_days', 'staleMute', 'flaky', 'slower', 'owner', 'ignore?', 'ignoreReason']
        .filter((n) => grid.col(n) != null);
      grid.columns.setVisible(visible);
      grid.columns.setOrder(visible);

      const widths: Record<string, number> = {mute: 40, test: 340, package: 130, passing: 100, duration: 100,
        needs_attention: 92, stale_days: 64, staleMute: 46, flaky: 54, slower: 54, owner: 110, 'ignore?': 54,
        ignoreReason: 150};
      for (const name of Object.keys(widths)) {
        const c = grid.col(name);
        if (c)
          c.width = widths[name];
      }
      for (const s of p.statusCols) {
        const c = grid.col(s);
        if (c)
          c.width = 36;
      }
      // Surface the problem tests first so the grid mirrors the Overview alerts.
      grid.sort(['needs_attention', 'flaky', 'slower', 'package', 'test'], [false, false, false, true, true]);
    });
  }

  private buildDrilldown(df: DG.DataFrame, i: number, p: ReleasePivot): HTMLElement {
    const test = df.get('test', i);
    const pkg = df.get('package', i);
    const status = df.get(`${p.latest}`, i);
    const result = df.get('latest_result', i) ?? '';
    const commit = df.get('latest_commit', i) ?? '';
    const resultBox = ui.div([ui.divText(result || 'No output for the latest build.')], 'ua-metrics-degraded');
    resultBox.style.whiteSpace = 'pre-wrap';
    resultBox.style.maxHeight = '320px';
    resultBox.style.overflow = 'auto';
    return ui.divV([
      ui.h2(test),
      ui.tableFromMap({Package: pkg, 'Latest status': status, Commit: commit || '—'}),
      ui.link('Open test-build in Jenkins', JENKINS_TEST_JOB),
      ui.divText('Latest build output', 'ua-metrics-panel-title'),
      resultBox,
    ]);
  }

  private async persistMute(df: DG.DataFrame): Promise<void> {
    const ignore = df.columns.byName('ignore?');
    if (!ignore)
      return;
    if (this.persisting) {
      this.pendingPersist = true; // an edit landed mid-write — coalesce into one more pass
      return;
    }
    this.persisting = true;
    try {
      this.schema ??= (await grok.dapi.stickyMeta.getSchemas()).find((s) => s.name === 'Autotests');
      if (!this.schema)
        return;
      const reason = df.columns.byName('ignoreReason');
      const n = df.rowCount;
      do {
        this.pendingPersist = false;
        const values = DG.DataFrame.fromColumns([
          DG.Column.fromList('bool', 'ignore?', Array.from({length: n}, (_, i) => ignore.get(i) === true)),
          DG.Column.fromStrings('ignoreReason', Array.from({length: n}, (_, i) => reason?.get(i) ?? '')),
        ]);
        await grok.dapi.stickyMeta.setAllValues(this.schema, df.getCol('test'), values);
      } while (this.pendingPersist);
    } finally {
      this.persisting = false;
    }
  }

  /** Toggles the current release version in the clicked test's mute list, updates derived flags, and persists. */
  private async toggleMute(df: DG.DataFrame, row: number, p: ReleasePivot): Promise<void> {
    if (!p.version) {
      grok.shell.warning('No release version detected for the latest build — cannot mute.');
      return;
    }
    const mvCol = df.columns.byName(MUTED_VERSIONS) ?? df.columns.addNewString(MUTED_VERSIONS);
    const versions = parseMutedVersions(mvCol.get(row));
    const at = versions.indexOf(p.version);
    if (at >= 0)
      versions.splice(at, 1);
    else
      versions.push(p.version);
    mvCol.set(row, versions.join(';'));

    const mutedForVersion = isMutedForVersion(mvCol.get(row), p.version);
    const globallyIgnored = df.columns.byName('ignore?')?.get(row) === true;
    const muted = mutedForVersion || globallyIgnored;
    df.getCol('mutedForVersion').set(row, mutedForVersion);
    df.getCol('muted').set(row, muted);
    df.getCol('mute').set(row, mutedForVersion ? MUTE_ON : MUTE_OFF);
    df.getCol('needs_attention').set(row, df.get('failing', row) === true && !muted);
    this.grid?.invalidate(); // repaint the toggled glyph immediately
    this.buildSummary(p); // muting changes the per-package pass/fail counts
    await this.persistReleaseMutes(df);
  }

  /** Toggles the current release version in the clicked test's STALE-mute list (independent of the failing mute). */
  private async toggleStaleMute(df: DG.DataFrame, row: number, p: ReleasePivot): Promise<void> {
    if (!p.version)
      return;
    const smvCol = df.columns.byName(STALE_MUTED_VERSIONS) ?? df.columns.addNewString(STALE_MUTED_VERSIONS);
    const versions = parseMutedVersions(smvCol.get(row));
    const at = versions.indexOf(p.version);
    if (at >= 0)
      versions.splice(at, 1);
    else
      versions.push(p.version);
    smvCol.set(row, versions.join(';'));
    const staleMuted = isMutedForVersion(smvCol.get(row), p.version);
    df.getCol('staleMutedForVersion').set(row, staleMuted);
    df.getCol('stale_alert').set(row, df.get('stale', row) === true && !staleMuted);
    df.getCol('staleMute').set(row, df.get('stale', row) !== true ? '' : (staleMuted ? MUTE_ON : MUTE_OFF));
    this.grid?.invalidate();
    await this.persistReleaseMutes(df);
  }

  /** Persists mutedVersions + staleMutedVersions back to the ReleaseMutes sticky-meta schema (coalesced). */
  private async persistReleaseMutes(df: DG.DataFrame): Promise<void> {
    const mvCol = df.columns.byName(MUTED_VERSIONS);
    const smvCol = df.columns.byName(STALE_MUTED_VERSIONS);
    if (!mvCol && !smvCol)
      return;
    if (this.persistingMutes) {
      this.pendingMutes = true;
      return;
    }
    this.persistingMutes = true;
    try {
      const schema = await getReleaseMutesSchema();
      const n = df.rowCount;
      do {
        this.pendingMutes = false;
        const cols: DG.Column[] = [];
        if (mvCol)
          cols.push(DG.Column.fromStrings(MUTED_VERSIONS, Array.from({length: n}, (_, i) => mvCol.get(i) ?? '')));
        if (smvCol)
          cols.push(DG.Column.fromStrings(STALE_MUTED_VERSIONS, Array.from({length: n}, (_, i) => smvCol.get(i) ?? '')));
        await grok.dapi.stickyMeta.setAllValues(schema, df.getCol('test'), DG.DataFrame.fromColumns(cols));
      } while (this.pendingMutes);
    } finally {
      this.persistingMutes = false;
    }
  }
}
