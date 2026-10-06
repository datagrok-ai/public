import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {fetchReleaseBenchmarks, computeBenchmarkAlerts, BenchmarkPivot, ReleaseContext, BENCH_SLOW_PCT,
  JENKINS_TEST_JOB, colorStatusCell, COLOR_FAIL_TEXT, waitForWidth, openInWorkspaceIcon} from './data';
import {combineLatest} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

const SLOW_BACK = 0xFFFFF3E0;
const SLOW_TEXT = 0xFFB35A00;

export class BenchmarksView extends UaView {
  private lastBuildsInput!: DG.InputBase;
  private slowPctInput!: DG.InputBase;
  private summaryHost!: HTMLElement;
  private gridHost!: HTMLElement;
  private pivot: BenchmarkPivot | null = null;
  private summaryDf: DG.DataFrame | null = null;
  private detailDf: DG.DataFrame | null = null;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Benchmarks';
    combineLatest([this.ctx.env, this.ctx.branch]).pipe(debounceTime(0))
      .subscribe(() => { if (this.initialized) this.refresh(); });
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics ua-metrics-fill';
    this.lastBuildsInput = ui.input.int('Builds', {value: 5, min: 2, max: 20, onValueChanged: () => this.refresh()});
    this.slowPctInput = ui.input.int('Slower by %', {value: BENCH_SLOW_PCT, min: 1, max: 500,
      onValueChanged: () => this.refresh()});
    this.slowPctInput.setTooltip('Flag a benchmark when the latest run exceeds its baseline by this much');
    this.summaryHost = ui.div([], 'ua-metrics-grid-host');
    this.gridHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([
      this.lastBuildsInput.root, this.slowPctInput.root,
      ui.div([], {style: {flex: '1'}}), refreshBtn,
    ], 'ua-metrics-header');

    const summaryPanel = ui.div([
      ui.divH([ui.divText('Benchmarks by package (latest build)', 'ua-metrics-panel-title'),
        openInWorkspaceIcon(() => this.summaryDf)], 'ua-metrics-panel-header'),
      this.summaryHost,
    ], 'ua-metrics-panel');
    const gridPanel = ui.div([
      ui.divH([ui.divText('Benchmarks (duration per build and per release line)', 'ua-metrics-panel-title'),
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
      this.pivot = await fetchReleaseBenchmarks(this.ctx.env.value, this.lastBuildsInput.value ?? 5,
        this.ctx.branch.value, this.slowPctInput.value ?? BENCH_SLOW_PCT);
      if (!this.pivot) {
        this.summaryHost.innerHTML = '';
        this.summaryHost.append(ui.divText('No benchmark runs reported for this instance and branch. ' +
          'The Benchmarks stage of Build-Deploy-Datagrok publishes them.', 'ua-metrics-degraded'));
        return;
      }
      this.buildSummary(this.pivot);
      this.buildGrid(this.pivot);
    } catch (e) {
      this.summaryHost.innerHTML = '';
      this.summaryHost.append(ui.divText(`ReleaseBenchmarks failed: ${e}`, 'ua-metrics-degraded'));
    }
  }

  private buildSummary(p: BenchmarkPivot): void {
    const a = computeBenchmarkAlerts(p);
    const df = p.df;
    const pkgCol = df.getCol('package');
    const perPkg = new Map<string, {total: number, slower: number, failed: number, worst: number, ms: number}>();
    for (let i = 0; i < df.rowCount; i++) {
      const pkg = pkgCol.get(i) ?? '';
      const row = perPkg.get(pkg) ?? {total: 0, slower: 0, failed: 0, worst: 0, ms: 0};
      row.total++;
      row.ms += df.get('latest_ms', i) ?? 0;
      if (df.get('slower', i) === true) {
        row.slower++;
        row.worst = Math.max(row.worst, df.get('slow_pct', i) ?? 0);
      }
      if (df.get('failing', i) === true)
        row.failed++;
      perPkg.set(pkg, row);
    }
    const pkgs = [...perPkg.entries()].sort((x, y) => y[1].slower - x[1].slower || x[0].localeCompare(y[0]));
    const summary = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('package', pkgs.map((e) => e[0])),
      DG.Column.fromInt32Array('benchmarks', Int32Array.from(pkgs.map((e) => e[1].total))),
      DG.Column.fromInt32Array('slower', Int32Array.from(pkgs.map((e) => e[1].slower))),
      DG.Column.fromInt32Array('failed', Int32Array.from(pkgs.map((e) => e[1].failed))),
      DG.Column.fromInt32Array('worst %', Int32Array.from(pkgs.map((e) => e[1].worst))),
      DG.Column.fromInt32Array('total ms', Int32Array.from(pkgs.map((e) => e[1].ms))),
    ]);
    summary.name = `By package — ${a.slower} slower, ${a.failed} failing of ${a.total}`;
    this.summaryDf = summary;
    const grid = DG.Viewer.grid(summary, {showColumnGridlines: false, allowBlockSelection: false});
    grid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;
      if ((gc.gridColumn.name === 'slower' || gc.gridColumn.name === 'worst %') && (gc.cell.value ?? 0) > 0)
        gc.style.textColor = SLOW_TEXT;
      if (gc.gridColumn.name === 'failed' && (gc.cell.value ?? 0) > 0)
        gc.style.textColor = COLOR_FAIL_TEXT;
    });
    this.summaryHost.innerHTML = '';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.summaryHost.append(grid.root);
  }

  private buildGrid(p: BenchmarkPivot): void {
    const df = p.df;
    this.detailDf = df;
    const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false,
      showCurrentCellOutline: false});

    this.gridHost.innerHTML = '';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.gridHost.append(grid.root);

    grid.onCellPrepare((gc) => {
      colorStatusCell(gc, p.statusCols);
      if (gc.isTableCell && p.statusCols.includes(gc.gridColumn.name)) {
        const v = gc.cell.value as string;
        gc.customText = v === 'passed' ? '+' : v === 'failed' ? '−' : v === 'skipped' ? '·' : '';
      }
      if (gc.isTableCell && (gc.gridColumn.name === 'slow_pct' || gc.gridColumn.name === 'latest_ms') &&
          df.get('slower', gc.cell.rowIndex) === true) {
        gc.style.backColor = SLOW_BACK;
        gc.style.textColor = SLOW_TEXT;
      }
    });
    df.onCurrentRowChanged.subscribe(() => {
      if (df.currentRowIdx >= 0)
        grok.shell.o = this.buildDrilldown(df, df.currentRowIdx, p);
    });

    waitForWidth(grid.root).then(() => {
      try {
        const trend = grid.columns.add({gridColumnName: 'per build', cellType: 'sparkline'});
        trend.settings = {columnNames: p.durationCols};
        if (p.versionCols.length > 1) {
          const byVersion = grid.columns.add({gridColumnName: 'per version', cellType: 'sparkline'});
          byVersion.settings = {columnNames: p.versionCols};
        }
      } catch (_) { /* grid not laid out yet — leave the raw per-build/per-version columns visible */ }

      const visible = ['test', 'package', 'per build', ...p.statusCols, 'per version', 'latest_ms',
        'baseline_ms', 'prev_release_ms', 'slow_pct', 'slower', 'failing', 'owner']
        .filter((n) => grid.col(n) != null);
      grid.columns.setVisible(visible);
      grid.columns.setOrder(visible);

      const widths: Record<string, number> = {test: 340, package: 130, 'per build': 100, 'per version': 100,
        latest_ms: 80, baseline_ms: 90, prev_release_ms: 110, slow_pct: 70, slower: 54, failing: 54, owner: 110};
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
      grid.sort(['slower', 'slow_pct', 'package', 'test'], [false, false, true, true]);
    });
  }

  private buildDrilldown(df: DG.DataFrame, i: number, p: BenchmarkPivot): HTMLElement {
    const result = df.get('latest_result', i) ?? '';
    const resultBox = ui.div([ui.divText(result || 'No output for the latest build.')], 'ua-metrics-degraded');
    resultBox.style.whiteSpace = 'pre-wrap';
    resultBox.style.maxHeight = '320px';
    resultBox.style.overflow = 'auto';
    const perVersion: Record<string, any> = {};
    for (const v of p.versions)
      perVersion[v] = df.get(`${v} ms`, i) ?? '—';
    const pct = df.get('slow_pct', i) ?? 0;
    return ui.divV([
      ui.h2(df.get('test', i)),
      ui.tableFromMap({
        Package: df.get('package', i),
        'Latest status': df.get(`${p.latest}`, i) ?? 'did not run',
        'Latest ms': df.get('latest_ms', i) || '—',
        'Baseline ms (median of prior builds)': df.get('baseline_ms', i) || '—',
        'Previous release ms': df.get('prev_release_ms', i) || '—',
        'Change': `${pct >= 0 ? '+' : ''}${pct}%`,
        Commit: df.get('latest_commit', i) || '—',
      }),
      ui.divText('Average per release line', 'ua-metrics-panel-title'),
      ui.tableFromMap(perVersion),
      ui.link('Open Jenkins', JENKINS_TEST_JOB),
      ui.divText('Latest build output', 'ua-metrics-panel-title'),
      resultBox,
    ]);
  }
}
