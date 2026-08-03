import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {fetchReleaseManual, ManualPivot, ReleaseContext, colorStatusCell, waitForWidth, openInWorkspaceIcon} from './data';

export class ManualTestsView extends UaView {
  private batchesInput!: DG.InputBase;
  private gridHost!: HTMLElement;
  private titleLabel!: HTMLElement;
  private manualDf: DG.DataFrame | null = null;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Manual';
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.batchesInput = ui.input.int('Batches', {value: 5, min: 2, max: 20, onValueChanged: () => this.refresh()});
    this.titleLabel = ui.divText('Manual test runs (Test Track batches)', 'ua-metrics-panel-title');
    this.gridHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.batchesInput.root, ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');
    const panel = ui.div([ui.divH([this.titleLabel, openInWorkspaceIcon(() => this.manualDf)], 'ua-metrics-panel-header'),
      this.gridHost], 'ua-metrics-panel');
    this.root.append(ui.div([header, panel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    this.gridHost.innerHTML = '';
    this.gridHost.append(ui.loader());
    try {
      const p = await fetchReleaseManual(this.batchesInput.value ?? 5);
      if (!p) {
        this.gridHost.innerHTML = '';
        this.gridHost.append(ui.divText('No manual test runs reported.', 'ua-metrics-degraded'));
        return;
      }
      this.titleLabel.textContent =
        `Latest batch: ${p.latestBatchName || '—'} — ${p.passed} passed, ${p.failed} failing of ${p.total}`;
      this.buildGrid(p);
    } catch (e) {
      this.gridHost.innerHTML = '';
      this.gridHost.append(ui.divText(`Manual Tests failed: ${e}`, 'ua-metrics-degraded'));
    }
  }

  private buildGrid(p: ManualPivot): void {
    const df = p.df;
    this.manualDf = df;
    const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false, showCurrentCellOutline: false});
    this.gridHost.innerHTML = '';
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.gridHost.append(grid.root);

    grid.onCellPrepare((gc) => colorStatusCell(gc, p.statusCols));
    df.onCurrentRowChanged.subscribe(() => {
      if (df.currentRowIdx >= 0)
        grok.shell.o = this.buildDrilldown(df, df.currentRowIdx);
    });

    // Sparkline needs the grid laid out first (see waitForWidth) — degrade gracefully otherwise.
    waitForWidth(grid.root).then(() => {
      try {
        const success = grid.columns.add({gridColumnName: 'success', cellType: 'sparkline'});
        success.settings = {columnNames: p.successCols};
      } catch (_) { /* grid not ready for the sparkline — leave raw per-batch columns */ }
      const order = ['test', 'success', ...p.statusCols, 'failing'];
      grid.columns.setOrder(order.filter((n) => grid.col(n) != null));
      for (const name of [...p.successCols, 'latest_result', ...p.chronological.map((bi) => `${bi} result`)]) {
        const c = grid.col(name);
        if (c)
          c.visible = false;
      }
    });
  }

  private buildDrilldown(df: DG.DataFrame, i: number): HTMLElement {
    const test = df.get('test', i);
    const result = df.get('latest_result', i) ?? '';
    const box = ui.div([ui.divText(result || 'No output for the latest batch.')], 'ua-metrics-degraded');
    box.style.whiteSpace = 'pre-wrap';
    box.style.maxHeight = '320px';
    box.style.overflow = 'auto';
    return ui.divV([ui.h2(test), ui.divText('Latest batch output', 'ua-metrics-panel-title'), box]);
  }
}
