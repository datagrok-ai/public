import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {SummaryMaking} from './sar-matrix-summary-making';
import {SummaryOverview} from './sar-matrix-summary-overview';
import {SummaryEffects} from './sar-matrix-summary-effects';
import {SummaryCollector, SummaryData} from './sar-matrix-summary-data';
import {count, scrollWithin} from '../sar-matrix-ui-common';
import {PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING, PANES} from './sar-matrix-summary-common';
import {SummaryScale} from './sar-matrix-summary-scale';
import {SummaryKit} from './sar-matrix-summary-kit';

/**
 * The Summary tab's frame: lifecycle, the section tabs and which pane is on screen. The scale is pinned
 * above every section because it invalidates every ranking under it. Each section is its own class; the
 * building blocks they share are in `SummaryKit`.
 *
 * The analog grid is its one Dart-backed object, built only when its section opens and closed on every
 * teardown path.
 */
export class SummaryPanel extends SummaryKit {
  readonly root = ui.divV([], 'chem-sar-main');
  private data: SummaryData | null = null;
  /** Depictions are painted one tick after the DOM lands, so a tab switch is not held up by RDKit. */
  private paintTimer = 0;
  /** Collecting holds the main thread, so it runs after the loader paints. */
  private collectTimer = 0;

  /** Set by the answer tile so the R-group leaderboard opens its top row when the Effects pane builds. */
  expandTopRGroup = false;
  /** Which Effects view is shown: a role card index, or past the last the per-series evidence. Kept across
   *  rebuilds. */
  effectsTab = 0;

  /** Which halves of the landing band the reader has shut, kept across rebuilds. */
  readonly foldedBands = new Set<string>();
  /** A fault message stands until invalidated, or collecting would replace it with an empty-state note. */
  private message = false;
  private sections: DG.TabControl | null = null;
  private sectionSub: Subscription | null = null;
  private readonly paneHosts = new Map<string, HTMLElement>();
  private currentPane = PANE_OVERVIEW;

  private readonly making = new SummaryMaking(this);
  private readonly overview = new SummaryOverview(this);
  private readonly effects = new SummaryEffects(this);
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

  /** Re-read the analysis at another fold tier. The matrices are not rebuilt, only walked differently. */
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

  /** Releases every timer, node and the grid, which a dropped reference does not free. */
  private reset(): void {
    window.clearTimeout(this.paintTimer);
    window.clearTimeout(this.collectTimer);
    this.collectTimer = 0;
    this.making.close();
    this.data = null;
    this.dropSections();
    this.currentPane = PANE_OVERVIEW;
    // A new analysis has its own tiers; a stale one would filter every matrix out.
    this.tierFilter = null;
    this.pendingPaints = [];
    ui.empty(this.root);
  }

  /** @param keepPlace Keep the section and view the reader is on: set when re-reading at another tier. */
  private render(keepPlace = false): void {
    const data = this.data!;
    const host = this.host;
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.making.close();
    this.dropSections();
    const pane = keepPlace ? this.currentPane : PANE_OVERVIEW;
    if (!keepPlace) {
      this.expandTopRGroup = false;
      this.effectsTab = 0;
    }
    ui.empty(this.root);

    if (host.matrices.length === 0) {
      // No tabs: three panes over nothing is three ways to read one message.
      this.root.appendChild(ui.divText(host.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }

    this.root.appendChild(ui.divH([this.scale.build(data), this.levelInput()], 'chem-sar-sum-top'));
    const sections = ui.tabControl(null, {vertical: false});
    sections.root.classList.add('chem-sar-sum-sections');
    for (const name of PANES) {
      const host = ui.div([], 'chem-sar-sum-pane');
      this.paneHosts.set(name, host);
      sections.addPane(name, () => host);
    }
    this.sections = sections;
    this.root.appendChild(sections.root);
    this.showPane(pane);
    this.sectionSub = sections.onTabChanged.subscribe(() => {
      const name = sections.currentPane?.name;
      if (name !== undefined && name !== this.currentPane)
        this.showPane(name);
    });
  }

  private dropSections(): void {
    this.sectionSub?.unsubscribe();
    this.sectionSub = null;
    this.sections = null;
    this.paneHosts.clear();
  }

  /** The fold tier every section reads, or the series count where there is only one tier. */
  private levelInput(): HTMLElement {
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
    const chip = ui.divText(`${label} · ${count(n)}`, 'chem-sar-chip-badge chem-sar-sum-level');
    chip.classList.toggle('chem-sar-sum-level-on', this.tierFilter === tier);
    ui.tooltip.bind(chip, () => (tier === null ? `Read every ranking over all ${count(n)} series.` :
      `Read every ranking over the ${count(n)} series cut at level ${tier} only.`) +
      ' L1 is the finest; each level up is one cut coarser.');
    chip.onclick = () => this.setTierFilter(tier);
    return chip;
  }

  private countTiers(): {tier: number, n: number}[] {
    const byTier = new Map<number, number>();
    for (const tier of this.host.matrixTiers)
      byTier.set(tier, (byTier.get(tier) ?? 0) + 1);
    return [...byTier.entries()].sort((a, b) => a[0] - b[0]).map(([tier, n]) => ({tier, n}));
  }

  /** Build one section and discard the last, so an unvisited Worth making never creates its grid. */
  showPane(name: string, anchor?: string): void {
    const paneHost = this.paneHosts.get(name);
    const sections = this.sections;
    const data = this.data;
    if (paneHost === undefined || sections === null || data === null)
      return;
    // Unconditional: reopening Worth making builds a second grid.
    this.making.close();
    this.currentPane = name;
    if (sections.currentPane?.name !== name)
      sections.currentPane = sections.getPane(name);
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    for (const host of this.paneHosts.values())
      ui.empty(host);
    // An anchor on this segment names a card, and cards live on tabs, so the tab is chosen first.
    if (name === PANE_EFFECTS && anchor !== undefined)
      this.effectsTab = this.effectsTabOf(data, anchor);
    paneHost.appendChild(name === PANE_EFFECTS ? this.effects.build(data) :
      name === PANE_MAKING ? this.making.build(data) : this.overview.build(data));

    this.paintTimer = window.setTimeout(() => {
      this.flushPaints();
      if (anchor === undefined)
        return;
      const target = paneHost.querySelector(`[data-anchor="${CSS.escape(anchor)}"]`);
      if (target instanceof HTMLElement)
        scrollWithin(paneHost, target);
    }, 0);
  }

  /** The Effects view an anchor lives on: a role anchor carries its index; the rest are per-series. */
  private effectsTabOf(data: SummaryData, anchor: string): number {
    const role = /^role-(\d+)$/.exec(anchor);
    if (role !== null)
      return Number(role[1]);
    if (data.axisRole === null)
      return 0;
    return this.roleFitRefusal(data) === null ? data.roleFit!.roles.length : 1;
  }
}
