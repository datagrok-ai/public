import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';

const LAST_BUILDS = 20;

export class StressView extends UaView {
  private byBuildHost!: HTMLElement;
  private summaryHost!: HTMLElement;
  private scatterHost!: HTMLElement;
  private scatterHeader!: HTMLElement;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Stress';
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.byBuildHost = ui.div([], 'ua-metrics-grid-host');
    this.summaryHost = ui.div([], 'ua-metrics-grid-host');
    this.scatterHost = ui.div([], 'ua-metrics-grid-host');
    this.scatterHeader = ui.divText('Last build', 'ua-metrics-panel-title');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    const byBuildPanel = ui.div([
      ui.divH([ui.divText(`Median duration by build (last ${LAST_BUILDS}, split by threads)`, 'ua-metrics-panel-title')],
        'ua-metrics-panel-header'),
      this.byBuildHost,
    ], 'ua-metrics-panel');
    const summaryPanel = ui.div([
      ui.divH([ui.divText('Metrics by build × threads', 'ua-metrics-panel-title')], 'ua-metrics-panel-header'),
      this.summaryHost,
    ], 'ua-metrics-panel');
    const scatterPanel = ui.div([
      ui.divH([this.scatterHeader], 'ua-metrics-panel-header'),
      this.scatterHost,
    ], 'ua-metrics-panel');

    this.root.append(ui.div([header, byBuildPanel, scatterPanel, summaryPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    await Promise.all([this.loadSummary(), this.loadScatter()]);
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
      if (df.rowCount === 0) {
        this.byBuildHost.innerHTML = '';
        this.byBuildHost.append(ui.divText('No stress runs reported yet.', 'ua-metrics-degraded'));
        return;
      }
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

  private async loadScatter(): Promise<void> {
    this.scatterHost.innerHTML = '';
    this.scatterHost.append(ui.loader());
    try {
      const df = await queries.stressTestsRaw(null);
      this.scatterHost.innerHTML = '';
      if (df.rowCount === 0) {
        this.scatterHost.append(ui.divText('No stress runs in the latest build.', 'ua-metrics-degraded'));
        return;
      }
      this.scatterHeader.textContent =
        `Last build: duration vs threads (${df.rowCount} runs, ${df.getCol('test').categories.length} tests)`;
      const scatter = DG.Viewer.scatterPlot(df, {
        xColumnName: 'threads',
        yColumnName: 'ms',
        colorColumnName: 'passed',
        yAxisType: 'logarithmic',
      });
      StressView.mount(this.scatterHost, scatter.root, '360px');
    } catch (e) {
      this.scatterHost.innerHTML = '';
      this.scatterHost.append(ui.divText(`StressTestsRaw failed: ${e}`, 'ua-metrics-degraded'));
    }
  }
}
