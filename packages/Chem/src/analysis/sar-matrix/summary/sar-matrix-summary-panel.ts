import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {SummaryMaking} from './sar-matrix-summary-making';
import {SummaryOverview} from './sar-matrix-summary-overview';
import {SummaryEffects} from './sar-matrix-summary-effects';
import {SummaryCollector, SummaryData} from './sar-matrix-summary-data';
import {count, scrollWithin, tipText} from '../sar-matrix-ui-common';
import {NARROW_PX, XNARROW_PX, SHORT_PX, PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING, PANE_METHOD,
  PANES} from './sar-matrix-summary-common';
import {SummaryMethod} from './sar-matrix-summary-method';
import {SummaryKit} from './sar-matrix-summary-kit';

/**
 * The Summary tab: the landing screen of the analysis. The scale and direction are pinned above
 * everything, because they invalidate every ranking under them; under that a segmented control opens
 * one of four panes, of which the first — the answers, which series to open, and what the analysis
 * covers — fits without scrolling at a normal dock size. Each row lands on the cell it describes.
 *
 * This class is the frame — lifecycle, the segment bar, which pane is on screen. Each segment is its
 * own class, and the building blocks they draw with are in `SummaryKit`.
 *
 * Holds one Dart-backed object, the analog grid, and it is built only when its own segment is opened;
 * every teardown path closes it. It holds no view and no `TableView`, so it can reach no dock node.
 */
export class SummaryPanel extends SummaryKit {
  readonly root = ui.divV([], 'chem-sar-main');
  private data: SummaryData | null = null;
  /** Depictions are painted one tick after the DOM lands, so a tab switch is not held up by RDKit. */
  private paintTimer = 0;
  /** Collecting over every cell of every matrix holds the main thread, so it runs off the activation
   *  stack — which also lets the loader paint. */
  private collectTimer = 0;

  /** Set by the answer tile so the R-group leaderboard opens its top row when the Effects pane builds. */
  expandTopRGroup = false;
  /** Which of the Effects segment's own tabs is open: an index into its role cards, or past the last
   *  of them the per-series evidence. Held across pane rebuilds so leaving the segment and coming back
   *  does not lose the component the reader was on. */
  effectsTab = 0;
  /** Whether the poorly-fitting list starts open: a reader who followed a "the fit does not hold"
   *  link came for that list, not for the explanation above it. */
  openTrustList = false;

  /** Which halves of the landing band the reader has shut, by key. Held across pane rebuilds, so a
   *  collapsed half stays collapsed when the segment is left and come back to. */
  readonly foldedBands = new Set<string>();
  /** The one line that reports the transfer scan, kept so its state can be refreshed without
   *  rebuilding the pane around it. */
  transferLine: HTMLElement | null = null;
  /** A fault line stands until something invalidates it: collecting over it would replace the reason
   *  the analysis failed with a generic empty-state note. */
  private message = false;
  private paneHost: HTMLElement | null = null;
  private currentPane = PANE_OVERVIEW;
  private readonly segButtons = new Map<string, HTMLElement>();
  /** The one scrolling region of the pane on screen, or null where the pane does not scroll. Moving
   *  this rather than calling `scrollIntoView` is what keeps a scroll inside the panel. */
  scroller: HTMLElement | null = null;
  private sizeSub: Subscription | null = null;
  /** The width breakpoints the pane on screen was built against; empty before the first build. */
  private widthKey = '';

  private readonly making = new SummaryMaking(this);
  private readonly overview = new SummaryOverview(this);
  private readonly effects = new SummaryEffects(this);
  private readonly method = new SummaryMethod(this);

  activateSummaryTab(): void {
    if (this.host.computing) {
      this.reset();
      this.root.appendChild(ui.divV([ui.loader(), ui.divText('Building SAR matrices...')],
        'chem-sar-empty-note'));
      return;
    }
    // Occupancy alone cannot say the tab is current: the loader above is a child too, and the
    // re-activation that follows a compute would then leave it standing forever.
    if (this.message || this.collectTimer !== 0 ||
      (this.data !== null && this.root.childElementCount > 0)) {
      // Detection runs on the Transfer tab and this pane is re-shown rather than rebuilt, so the one
      // line that reports it would otherwise keep inviting a scan that has already run.
      this.syncTransferLine();
      return;
    }
    if (this.data !== null) {
      this.render();
      return;
    }
    // Letting the loader paint first also lets the dock finish laying out its tab strip before this
    // pane's content lands.
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

  /** Re-read the analysis at another fold tier and rebuild the segment on screen. Nothing is recomputed
   *  in the viewer: the matrices are already built, and this only changes which of them are walked. */
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
    this.sizeSub?.unsubscribe();
    this.sizeSub = null;
  }

  /** Every timer, every Dart-backed object and every rendered node this panel owns. The grid is the
   *  one thing a dropped reference does not release, so it is closed on all three teardown paths. */
  private reset(): void {
    window.clearTimeout(this.paintTimer);
    window.clearTimeout(this.collectTimer);
    this.collectTimer = 0;
    this.making.close();
    this.data = null;
    this.paneHost = null;
    this.scroller = null;
    this.currentPane = PANE_OVERVIEW;
    // A new analysis has its own tiers. Kept, a tier the new run does not hold filters every matrix
    // out, which reads as an empty dataset — and the chip bar that would clear it is built from the
    // tiers that survived the filter, so there is none.
    this.tierFilter = null;
    this.widthKey = '';
    this.segButtons.clear();
    this.pendingPaints = [];
    ui.empty(this.root);
  }

  /**
   * @param keepPlace Keep the segment and sub-tab the reader is on. Set when the same analysis is
   * re-read — changing the fold tier — and left alone for a new one.
   */
  private render(keepPlace = false): void {
    const data = this.data!;
    const host = this.host;
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.transferLine = null;
    this.making.close();
    this.segButtons.clear();
    this.paneHost = null;
    this.scroller = null;
    // A new analysis is a new landing: a reader returning to a different run must not arrive mid-way
    // into someone else's orientation. Re-reading the same one at another tier is the opposite case —
    // the reader asked a question about the segment they are looking at and has to stay on it.
    const pane = keepPlace ? this.currentPane : PANE_OVERVIEW;
    if (!keepPlace) {
      this.expandTopRGroup = false;
      this.effectsTab = 0;
      this.openTrustList = false;
    }
    ui.empty(this.root);

    if (host.matrices.length === 0) {
      // No segment bar: four panes over nothing is four ways to read one message.
      this.root.appendChild(ui.divText(host.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }

    this.root.appendChild(this.method.scaleBand(data));
    this.root.appendChild(this.buildSegBar());
    this.paneHost = ui.div([], 'chem-sar-sum-pane');
    this.root.appendChild(this.paneHost);
    // Before the first build, not after: a pane decides here whether it draws a mark or its text
    // fallback, and it can only read a class that is already on the root.
    this.widthKey = '';
    this.applySize();
    this.showPane(pane);

    if (this.sizeSub === null) {
      // The pane is docked, so its width is not the viewport's: a media query would key on the wrong
      // number entirely.
      this.sizeSub = DG.debounce(ui.onSizeChanged(this.root), 150).subscribe(() => this.applySize());
    }
  }

  /**
   * The breakpoint classes, and a rebuild of the pane on screen whenever a width one crosses.
   *
   * Height is read as well as width: the Overview's no-scroll look is an `overflow: hidden`, and a
   * dock that is wide but short would otherwise clip its bottom row with no scrollbar to say so.
   * The rebuild is what keeps a mark and its text fallback from both vanishing — CSS can hide the
   * reason lane, only a rebuild can put the phrases back.
   */
  private applySize(): void {
    const w = this.root.clientWidth;
    if (w === 0)
      return;
    const narrow = this.root.classList.toggle('chem-sar-sum-narrow', w < NARROW_PX);
    const xnarrow = this.root.classList.toggle('chem-sar-sum-xnarrow', w < XNARROW_PX);
    this.root.classList.toggle('chem-sar-sum-short', this.root.clientHeight < SHORT_PX);
    const key = `${narrow}|${xnarrow}`;
    const crossed = this.widthKey !== '' && this.widthKey !== key;
    this.widthKey = key;
    if (crossed && this.paneHost !== null)
      this.showPane(this.currentPane);
  }

  private tierChip(tier: number | null, label: string, n: number): HTMLElement {
    // The tier name goes in its own element so a reader or a test can address it without the count
    // beside it, which changes with the data.
    const chip = ui.divH([ui.divText(label, 'chem-sar-sum-tier-name'),
      ui.divText(`· ${count(n)}`, 'chem-sar-sum-tier-n')], 'chem-sar-chip-badge chem-sar-sum-role');
    chip.classList.toggle('chem-sar-sum-role-on', this.tierFilter === tier);
    ui.tooltip.bind(chip, () => (tier === null ?
      `Read every ranking on this tab over all ${count(n)} series at once. A compound cut at several ` +
      'breadths is then counted at each of them.' :
      `Read every ranking on this tab over the ${count(n)} series cut at level ${tier} only — the ` +
      'components, the swaps, the cores, Start here and Worth making.') +
      ' L1 is a leaf series; each level above it holds the same compounds over cores that agree one ' +
      'further cut deeper. Nothing is rebuilt: the series already exist either way.');
    chip.onclick = () => this.setTierFilter(tier);
    return chip;
  }

  /**
   * The second level of tabs, built from plain divs rather than `ui.tabControl`.
   *
   * A nested tab control renders the platform's own tab chrome 26px under the viewer's, and two
   * identical strips is exactly the confusion a second level has to avoid. This differs in shape, in
   * active state, in size and position, and in carrying a readout on its right — which a tab strip
   * never does and a toolbar always does.
   */
  private buildSegBar(): HTMLElement {
    const segs = PANES.map((name) => {
      const el = ui.divText(name, 'chem-sar-sum-seg');
      el.onclick = () => this.showPane(name);
      this.segButtons.set(name, el);
      return el;
    });
    // On the segment bar rather than on one pane: the tier decides what every segment is reading, so it
    // has to be visible and changeable from all of them. A decomposition handed over as columns is all
    // one tier, and there the chips would be a control with one setting.
    const readout: HTMLElement[] = [];
    const tierCounts = this.countTiers();
    if (tierCounts.length > 1) {
      readout.push(this.tierChip(null, 'All', this.host.matrixTiers.length));
      for (const {tier, n} of tierCounts)
        readout.push(this.tierChip(tier, `L${tier}`, n));
    } else {
      const families = new Set(this.host.matrixRoots).size;
      readout.push(tipText(`${count(this.host.matrices.length)} series · ${count(families)} families`,
        'chem-sar-sum-seg-readout', 'Matrices over fold lineages. The same compounds re-cut, not this many ' +
        'findings.'));
    }
    return ui.divH([ui.divH(segs, 'chem-sar-sum-seg-group'),
      ui.divH(readout, 'chem-sar-sum-seg-tiers')], 'chem-sar-sum-seg-bar');
  }

  private countTiers(): {tier: number, n: number}[] {
    const byTier = new Map<number, number>();
    for (const tier of this.host.matrixTiers)
      byTier.set(tier, (byTier.get(tier) ?? 0) + 1);
    return [...byTier.entries()].sort((a, b) => a[0] - b[0]).map(([tier, n]) => ({tier, n}));
  }

  /** Build one segment and discard the last. Panes are never held: an unvisited Worth-making segment
   *  then never creates a DataFrame, and a revisited one is rebuilt from data that has not changed. */
  showPane(name: string, anchor?: string): void {
    const paneHost = this.paneHost;
    const data = this.data;
    if (paneHost === null || data === null)
      return;
    // Unconditional: re-opening the segment that holds it builds a second grid, and only closing the
    // first releases the Dart-backed object behind it.
    this.making.close();
    this.currentPane = name;
    for (const [key, el] of this.segButtons)
      el.classList.toggle('chem-sar-sum-seg-on', key === name);
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.scroller = null;
    ui.empty(paneHost);
    // An anchor on this segment names a card, and cards live on tabs, so the tab is chosen first.
    if (name === PANE_EFFECTS && anchor !== undefined)
      this.effectsTab = this.effectsTabOf(data, anchor);
    paneHost.classList.toggle('chem-sar-sum-pane-fixed', name === PANE_OVERVIEW);
    paneHost.appendChild(name === PANE_EFFECTS ? this.effects.build(data) :
      name === PANE_MAKING ? this.making.build(data) :
        name === PANE_METHOD ? this.method.build(data) : this.overview.build(data));

    this.paintTimer = window.setTimeout(() => {
      this.flushPaints();
      if (anchor === undefined)
        return;
      const target = paneHost.querySelector(`[data-anchor="${anchor}"]`);
      if (target instanceof HTMLElement && this.scroller !== null)
        scrollWithin(this.scroller, target);
    }, 0);
  }

  /** What the transfer scan has found, said in the one line that offers it. Detection is lazy and runs
   *  on its own tab, so this is read at activation rather than when the pane was built. */
  syncTransferLine(): void {
    const line = this.transferLine;
    if (line === null)
      return;
    const scan = this.host.transferSummary;
    line.textContent = !scan.scanned ?
      'SAR transfers — pairs of cores whose trends run in parallel · open the tab to detect →' :
      scan.count === 0 ?
        'No SAR transfers: no two cores share enough R-groups, on scaffolds alike enough · open the ' +
        'tab to change the similarity →' :
        `${count(scan.count)} SAR transfers — pairs of cores whose trends run in parallel · open the ` +
        'tab →';
  }

  /** Switch to the methodology segment with the trust list on screen — where every "the fit does not
   *  hold here" on the tab resolves to. */
  revealTrust(): void {
    this.openTrustList = true;
    this.showPane(PANE_METHOD, 'trust');
  }

  /** The Effects tab an anchor lives on: a role anchor carries its own index, and everything else is
   *  on the per-series tab, which sits after the role cards. */
  private effectsTabOf(data: SummaryData, anchor: string): number {
    const role = /^role-(\d+)$/.exec(anchor);
    if (role !== null)
      return Number(role[1]);
    if (data.axisRole === null)
      return 0;
    return this.roleFitRefusal(data) === null ? data.roleFit!.roles.length : 1;
  }
}
