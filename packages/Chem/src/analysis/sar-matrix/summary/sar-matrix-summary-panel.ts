import * as ui from 'datagrok-api/ui';

import {SummaryMaking} from './sar-matrix-summary-making';
import {SummaryOverview} from './sar-matrix-summary-overview';
import {SummaryEffects} from './sar-matrix-summary-effects';
import {SummaryCollector, SummaryData} from './sar-matrix-summary-data';
import {chipBadge, count} from '../sar-matrix-ui-common';
import {SummaryScale} from './sar-matrix-summary-scale';
import {SummaryKit} from './sar-matrix-summary-kit';

export class SummaryPanel extends SummaryKit {
  readonly root = ui.divV([], 'chem-sar-main chem-sar-sum-root');
  private data: SummaryData | null = null;
  private paintTimer = 0;
  /** Collecting holds the main thread, so it runs after the loader paints. */
  private collectTimer = 0;

  /** A fault message stands until invalidated, or collecting would replace it with an empty-state note. */
  private message = false;

  private readonly making = new SummaryMaking(this);
  private readonly overview = new SummaryOverview(this);
  readonly effects = new SummaryEffects(this);
  private readonly scale = new SummaryScale(this);

  activateSummaryTab(): void {
    if (this.host.computing) {
      this.reset();
      this.root.appendChild(ui.divV([ui.loader(), ui.divText('Building SAR matrices...')],
        'chem-sar-empty-note'));
      return;
    }
    // Occupancy alone cannot say the tab is current: the loader is a child too.
    if (this.message || this.collectTimer !== 0 ||
      (this.data !== null && this.root.childElementCount > 0))
      return;
    if (this.data !== null) {
      this.render();
      return;
    }
    ui.empty(this.root);
    this.root.appendChild(ui.divV([ui.loader(), ui.divText('Reading the analysis...')],
      'chem-sar-empty-note'));
    this.collectTimer = window.setTimeout(() => {
      this.collectTimer = 0;
      if (this.host.computing || this.message)
        return;
      this.data = new SummaryCollector(this.host, this.tierFilter).collect();
      this.render();
    }, 0);
  }

  invalidate(): void {
    this.reset();
    this.message = false;
  }

  private setTierFilter(tier: number | null): void {
    if (this.tierFilter === tier)
      return;
    this.tierFilter = tier;
    this.data = new SummaryCollector(this.host, this.tierFilter).collect();
    this.render(true);
  }

  showMessage(text: string): void {
    this.reset();
    this.message = true;
    this.root.appendChild(ui.divText(text, 'chem-sar-empty-note'));
  }

  release(): void {
    this.invalidate();
  }

  /** Releases every timer and node, which a dropped reference does not free. */
  private reset(): void {
    window.clearTimeout(this.paintTimer);
    window.clearTimeout(this.collectTimer);
    this.collectTimer = 0;
    this.data = null;
    // A new analysis has its own tiers; a stale one would filter every matrix out.
    this.tierFilter = null;
    this.pendingPaints = [];
    this.stopPainting();
    ui.empty(this.root);
  }

  /** @param keepPlace Keep the reader's scroll position, for a re-read at another tier. */
  private render(keepPlace = false): void {
    const host = this.host;
    const scrollers = (): HTMLElement[] => [...this.root.querySelectorAll<HTMLElement>(
      '.chem-sar-sum-pane, .chem-sar-sum-flowcol, .chem-sar-sum-making-block')];
    const places = keepPlace ? scrollers().map((el) => el.scrollTop) : [];
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.stopPainting();
    ui.empty(this.root);

    if (host.matrices.length === 0) {
      this.root.appendChild(ui.divText(host.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }
    const data = this.data!;
    // A short analog list would leave the right column empty while the left one scrolls, so the answers join it.
    const keySeries = this.overview.keySeries(data);
    const making = this.making.build(data);
    const answers = this.overview.answers(data);
    const short = this.making.isShort(data);
    const right = short ? ui.divV([making, answers]) : making;
    right.classList.add('chem-sar-sum-flowcol');
    const left = ui.divV(short ? [keySeries] : [keySeries, answers], 'chem-sar-sum-flowcol');
    this.root.appendChild(ui.divV([this.overview.head(data), ui.divH([left, right], 'chem-sar-sum-flow')],
      'chem-sar-sum-pane'));
    scrollers().forEach((el, i) => el.scrollTop = places[i] ?? 0);
    this.paintTimer = window.setTimeout(() => this.flushPaints(), 0);
  }

  scaleRow(data: SummaryData, coverage: HTMLElement): HTMLElement {
    return ui.divH([this.scale.build(data), coverage], 'chem-sar-sum-top');
  }

  levelInput(): HTMLElement {
    const tiers = this.countTiers();
    if (tiers.length <= 1) {
      const families = new Set(this.host.matrixRoots).size;
      return this.hint(`${count(this.host.matrices.length)} series · ${count(families)} families`,
        'Matrices over fold lineages. The same compounds re-cut, not this many findings.');
    }
    const chips = [this.levelChip(null, 'All', this.host.matrixTiers.length),
      ...tiers.map(({tier, n}) => this.levelChip(tier, `L${tier}`, n))];
    return ui.divH(chips, 'chem-sar-sum-levels');
  }

  private levelChip(tier: number | null, label: string, n: number): HTMLElement {
    const chip = chipBadge(`${label} · ${count(n)}`, (tier === null ? `Read every ranking over all ${count(n)} ` +
      'series.' : `Read every ranking over the ${count(n)} series cut at level ${tier} only.`) +
      ' L1 is the finest; each level up is one cut coarser.', 'chem-sar-sum-level');
    chip.classList.toggle('chem-sar-sum-level-on', this.tierFilter === tier);
    chip.onclick = () => this.setTierFilter(tier);
    return chip;
  }

  private countTiers(): {tier: number, n: number}[] {
    const byTier = new Map<number, number>();
    for (const tier of this.host.matrixTiers)
      byTier.set(tier, (byTier.get(tier) ?? 0) + 1);
    return [...byTier.entries()].sort((a, b) => a[0] - b[0]).map(([tier, n]) => ({tier, n}));
  }
}
