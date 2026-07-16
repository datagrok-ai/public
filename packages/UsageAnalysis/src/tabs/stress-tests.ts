import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';
import {stressRegression, openInWorkspaceIcon, ReleaseContext} from '../release/data';

const LAST_BUILDS = 20;

export class StressView extends UaView {
  private byBuildHost!: HTMLElement;
  private summaryHost!: HTMLElement;
  private boxHost!: HTMLElement;
  private boxHeader!: HTMLElement;
  private regressionHost!: HTMLElement;
  private summaryDf: DG.DataFrame | null = null;

  constructor(uaToolbox?: UaToolbox, ctx?: ReleaseContext) {
    super(uaToolbox);
    this.name = 'Stress';
    ctx?.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.byBuildHost = ui.div([], 'ua-metrics-grid-host ua-metrics-chart-host');
    this.summaryHost = ui.div([], 'ua-metrics-grid-host');
    this.boxHost = ui.div([], 'ua-metrics-grid-host ua-metrics-chart-host');
    this.boxHeader = ui.divText('Last build', 'ua-metrics-panel-title');

    this.regressionHost = ui.div([], 'ua-metrics-asof');
    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.regressionHost, ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    const byBuildPanel = ui.div([
      ui.divH([ui.divText(`Median duration by build (last ${LAST_BUILDS}, split by threads)`, 'ua-metrics-panel-title')],
        'ua-metrics-panel-header'),
      this.byBuildHost,
    ], 'ua-metrics-panel');
    const summaryPanel = ui.div([
      ui.divH([ui.divText('Metrics by build × threads', 'ua-metrics-panel-title'),
        openInWorkspaceIcon(() => this.summaryDf)], 'ua-metrics-panel-header'),
      this.summaryHost,
    ], 'ua-metrics-panel');
    const boxPanel = ui.div([
      ui.divH([this.boxHeader], 'ua-metrics-panel-header'),
      this.boxHost,
    ], 'ua-metrics-panel');

    this.root.append(ui.div([header, byBuildPanel, boxPanel, summaryPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    await Promise.all([this.loadSummary(), this.loadBox()]);
  }

  private static mount(host: HTMLElement, el: HTMLElement, height: string): void {
    el.style.width = '100%';
    el.style.height = height;
    host.innerHTML = '';
    host.append(el);
  }

  private async loadSummary(): Promise<void> {
    this.byBuildHost.innerHTML = '';
    this.byBuildHost.append(ui.loader());
    this.summaryHost.innerHTML = '';
    try {
      const df = await queries.stressTestsSummary(LAST_BUILDS);
      this.summaryDf = df;
      if (df.rowCount === 0) {
        this.byBuildHost.innerHTML = '';
        this.byBuildHost.append(ui.divText('No stress runs reported yet.', 'ua-metrics-degraded'));
        return;
      }
      const reg = stressRegression(df);
      this.regressionHost.textContent = reg.detail;
      this.regressionHost.className = `ua-metrics-asof ${reg.regressed ? 'ua-metrics-red' : 'ua-metrics-green'}`;
      df.columns.addNewString('threads ').init((i) => `${df.get('threads', i)} threads`);
      const line = DG.Viewer.lineChart(df, {
        xColumnName: 'build',
        yColumnNames: ['median_ms'],
        splitColumnNames: ['threads '],
      });
      StressView.mount(this.byBuildHost, line.root, '320px');

      const grid = DG.Viewer.grid(df, {'showColumnGridlines': false, 'allowBlockSelection': false});
      const threadsSplitCol = grid.col('threads ');
      if (threadsSplitCol)
        threadsSplitCol.visible = false;
      StressView.mount(this.summaryHost, grid.root, '100%');
    } catch (e) {
      this.byBuildHost.innerHTML = '';
      this.byBuildHost.append(ui.divText(`StressTestsSummary failed: ${e}`, 'ua-metrics-degraded'));
    }
  }

  private async loadBox(): Promise<void> {
    this.boxHost.innerHTML = '';
    this.boxHost.append(ui.loader());
    try {
      const df = await queries.stressTestsRaw(null);
      this.boxHost.innerHTML = '';
      if (df.rowCount === 0) {
        this.boxHost.append(ui.divText('No stress runs in the latest build.', 'ua-metrics-degraded'));
        return;
      }
      // Zero-padded thread count so the categorical axis sorts numerically (008 < 016 < 032).
      const threads = df.getCol('threads');
      df.columns.addNewString('threads_string').init((i) => `${threads.get(i) ?? ''}`.padStart(3, '0'));
      this.boxHeader.textContent =
        `Last build: duration distribution by threads (${df.rowCount} runs, ${df.getCol('test').categories.length} tests)`;
      const box = DG.Viewer.boxPlot(df, {
        valueColumnName: 'ms',
        categoryColumnNames: ['threads_string'],
        markerColorColumnName: 'passed',
        axisType: 'logarithmic',
      });
      StressView.mount(this.boxHost, box.root, '360px');
    } catch (e) {
      this.boxHost.innerHTML = '';
      this.boxHost.append(ui.divText(`StressTestsRaw failed: ${e}`, 'ua-metrics-degraded'));
    }
  }
}
