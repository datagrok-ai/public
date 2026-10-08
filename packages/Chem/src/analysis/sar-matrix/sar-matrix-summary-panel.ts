import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {_package} from '../../package';
import {ROLE_FIT_MAX_SWEEPS} from './sar-matrix-role-fit';
import {SarMatrix} from './sar-matrix-types';
import {SummaryMaking} from './sar-matrix-summary-making';
import {SummaryOverview} from './sar-matrix-summary-overview';
import {SummaryEffects} from './sar-matrix-summary-effects';
import {MIN_SUPPORT, SeriesStat, SummaryCollector, SummaryData, SummaryHost, SUM_ROWS,
  TRUST_R2} from './sar-matrix-summary-data';

export {SummaryHost} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, CORE_BG_ARGB, MatrixCellRef, paintMoleculeOnColor, STRIP_MOL_H, STRIP_MOL_W,
  count} from './sar-matrix-ui-common';
import {TRUST_LIST_MAX, NARROW_PX, XNARROW_PX, SHORT_PX, PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING, PANE_METHOD,
  PANES} from './sar-matrix-summary-common';

/**
 * The Summary tab: the landing screen of the analysis. The scale and direction are pinned above
 * everything, because they invalidate every ranking under them; under that a segmented control opens
 * one of four panes, of which the first — the answers, which series to open, and what the analysis
 * covers — fits without scrolling at a normal dock size. Each row lands on the cell it describes.
 *
 * Holds one Dart-backed object, the analog grid, and it is built only when its own segment is opened;
 * every teardown path closes it. It holds no view and no `TableView`, so it can reach no dock node.
 */
export class SummaryPanel {
  readonly root = ui.divV([], 'chem-sar-main');
  private data: SummaryData | null = null;
  /** Depictions are painted one tick after the DOM lands, so a tab switch is not held up by RDKit. */
  private paintTimer = 0;
  /** Collecting over every cell of every matrix holds the main thread, so it runs off the activation
   *  stack — which also lets the loader paint. */
  private collectTimer = 0;
  pendingPaints: (() => void)[] = [];

  /** Set by the answer tile so the R-group leaderboard opens its top row when the Effects pane builds. */
  expandTopRGroup = false;
  /** Which of the Effects segment's own tabs is open: an index into its role cards, or past the last
   *  of them the per-series evidence. Held across pane rebuilds so leaving the segment and coming back
   *  does not lose the component the reader was on. */
  effectsTab = 0;
  /** Whether the poorly-fitting list starts open: a reader who followed a "the fit does not hold"
   *  link came for that list, not for the explanation above it. */
  private openTrustList = false;
  /**
   * The fold tier every ranking on this tab is read at, or null for all of them together.
   *
   * Not a cut of the data: the matrices are built once and all of them stay. This only decides which
   * of them the tab's rankings walk — a tier holds the same compounds as the tier below it over cores
   * cut one bond broader, so reading at one tier answers "what does the SAR look like at this breadth"
   * without the other breadths mixed in.
   */
  private tierFilter: number | null = null;

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

  constructor(readonly host: SummaryHost) {}

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
    this.making.close();
    this.data = new SummaryCollector(this.host, this.tierFilter).collect();
    this.render(true);
  }

  showMessage(text: string): void {
    this.reset();
    this.message = true;
    this.root.appendChild(ui.divText(text, 'chem-sar-empty-note'));
  }

  release(): void {
    this.reset();
    this.message = false;
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

  /** `scrollIntoView` scrolls every scrollable ancestor, the dock container included, so it can move
   *  the layout around the panel; only this pane's own scroller may move. */
  private scrollTo(el: HTMLElement): void {
    const s = this.scroller;
    if (s === null)
      return;
    s.scrollTop += el.getBoundingClientRect().top - s.getBoundingClientRect().top;
  }

  // ---- Collection -----------------------------------------------------------------------------

  // ---- Layout ---------------------------------------------------------------------------------

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
    this.currentPane = PANE_OVERVIEW;
    ui.empty(this.root);

    if (host.matrices.length === 0) {
      // No segment bar: four panes over nothing is four ways to read one message.
      this.root.appendChild(ui.divText(host.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }

    this.root.appendChild(this.buildScaleBand(data));
    this.root.appendChild(this.buildSegBar(data));
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

  /**
   * The second level of tabs, built from plain divs rather than `ui.tabControl`.
   *
   * A nested tab control renders the platform's own tab chrome 26px under the viewer's, and two
   * identical strips is exactly the confusion a second level has to avoid. This differs in shape, in
   * active state, in size and position, and in carrying a readout on its right — which a tab strip
   * never does and a toolbar always does.
   */
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

  private buildSegBar(data: SummaryData): HTMLElement {
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
    if (data.tierCounts.length > 1) {
      const total = data.tierCounts.reduce((sum, {n}) => sum + n, 0);
      readout.push(this.tierChip(null, 'All', total));
      for (const {tier, n} of data.tierCounts)
        readout.push(this.tierChip(tier, `L${tier}`, n));
    } else {
      const line = ui.divText(`${count(this.host.matrices.length)} series · ` +
        `${count(data.families)} families`, 'chem-sar-sum-seg-readout');
      ui.tooltip.bind(line, 'Matrices over fold lineages. The same compounds re-cut, not this ' +
        'many findings.');
      readout.push(line);
    }
    return ui.divH([ui.divH(segs, 'chem-sar-sum-seg-group'),
      ui.divH(readout, 'chem-sar-sum-seg-tiers')], 'chem-sar-sum-seg-bar');
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
    // An anchor on this segment names a card, and a card now lives on a tab rather than at a scroll
    // offset, so it selects the tab before the pane is built.
    if (name === PANE_EFFECTS && anchor !== undefined)
      this.effectsTab = this.effectsTabOf(data, anchor);
    paneHost.classList.toggle('chem-sar-sum-pane-fixed', name === PANE_OVERVIEW);
    paneHost.appendChild(name === PANE_EFFECTS ? this.effects.build(data) :
      name === PANE_MAKING ? this.making.build(data) :
        name === PANE_METHOD ? this.buildMethodPane(data) : this.overview.build(data));

    this.paintTimer = window.setTimeout(() => {
      this.flushPaints();
      if (anchor === undefined)
        return;
      const target = paneHost.querySelector(`[data-anchor="${anchor}"]`);
      if (target instanceof HTMLElement)
        this.scrollTo(target);
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

  anchored(el: HTMLElement, anchor: string): HTMLElement {
    el.dataset.anchor = anchor;
    return el;
  }

  // ---- The four panes ---------------------------------------------------------------------------

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

  private buildMethodPane(data: SummaryData): HTMLElement {
    const parts: HTMLElement[] = [this.buildChips(data),
      this.anchored(this.buildTrustSection(data), 'trust')];
    const scroll = ui.divV(parts, 'chem-sar-sum-scroll');
    this.scroller = scroll;
    return scroll;
  }

  // ---- Band A: scale, then method ---------------------------------------------------------------

  /**
   * The two things that invalidate every ranking below, and nothing else — this band does not scroll,
   * so anything pinned here is width the answers never get back.
   */
  private buildScaleBand(data: SummaryData): HTMLElement {
    const host = this.host;
    const lines: HTMLElement[] = [];
    const line = (text: string): HTMLElement => ui.divText(text, 'chem-sar-sum-orient-line');

    const transform = host.scalingLabel === 'raw' ? 'untransformed' : `${host.scalingLabel} applied`;
    const direction = host.higherIsBetter ? 'higher is better' : 'lower is better';
    // Never printed as 0 when no series has a fit: a floor of zero reads as "everything is resolved".
    const floor = data.modelError === null ? '' :
      ` · ± ${host.formatActivity(data.modelError)} model error`;
    // Log-ness of a raw column is inferred from the declared direction, not read off the data, and
    // every fold and log-unit claim on the tab rests on it. Percent inhibition, ΔTm and ΔG are raw and
    // higher-is-better too, so the inference has to be on screen rather than in a getter.
    // Stated on the hover rather than on the band: it qualifies the fold figures further down, and on
    // the band it read as one more property of the column.
    const assumedLog = host.scalingLabel === 'raw' && host.activityIsLog;
    const scaleLine = line(`${host.activityColumnName || 'Activity'} · ${transform} · ${direction}`);
    ui.tooltip.bind(scaleLine, () => 'The scale and direction every ranking here depends on.' +
      (!assumedLog ? '' : ' Fold changes on this tab assume the column is already a log, inferred from ' +
      'its being untransformed and higher-is-better — a percent inhibition, ΔTm or ΔG column is the ' +
      'same shape and is not a log.') +
      (data.modelError === null ? '' : ' Model error, not assay error: no replicates here, so it is an ' +
      'upper bound on the assay.'));
    const head: HTMLElement[] = [scaleLine];
    if (data.minObserved !== null)
      head.push(this.rangeRule(data));
    else
      head.push(line('nothing observed'));
    if (floor !== '')
      head.push(line(floor.replace(/^ · /, '')));
    lines.push(ui.divH(head, 'chem-sar-sum-orient-row'));
    // Only where this analysis did transform: a column it left alone is on the scale it was measured
    // on, and log permeability or ΔG is negative throughout without anything being wrong.
    if (data.maxObserved !== null && data.maxObserved < 0 && host.scalingLabel !== 'raw') {
      // Redundant with the rule's own zero mark, and deliberately so: a double transform inverts every
      // ranking on the tab, which is the one consequence worth saying twice.
      lines.push(ui.divText('Every observed value is negative — check Scaling.',
        'chem-sar-sum-orient-line chem-sar-sum-alarm'));
    }
    return ui.divV(lines, 'chem-sar-sum-orient');
  }

  // ---- Marks ------------------------------------------------------------------------------------

  /**
   * The observed range as a rule, spanning exactly lo..hi. It must not extend to zero: on a column that
   * is negative throughout — log solubility, ΔG — the right-hand label would sit at the track's end
   * where the value is zero and not the label's number.
   */
  private rangeRule(data: SummaryData): HTMLElement {
    const host = this.host;
    const lo = data.minObserved!;
    const hi = data.maxObserved!;
    const span = hi - lo || 1;
    const track = ui.div([], 'chem-sar-sum-rule');
    const bar = ui.div([], 'chem-sar-sum-rule-bar');
    bar.style.left = '0%';
    bar.style.width = '100%';
    track.appendChild(bar);
    // Only where the measurements actually straddle it, which is the only place the mark can be read
    // off the rule at all.
    const showZero = lo < 0 && hi > 0;
    if (showZero) {
      const zero = ui.div([], 'chem-sar-sum-rule-zero');
      zero.style.left = `${((0 - lo) / span) * 100}%`;
      track.appendChild(zero);
    }
    const box = ui.divH([
      ui.divText(host.formatActivity(lo), 'chem-sar-sum-rule-end'),
      track,
      ui.divText(host.formatActivity(hi), 'chem-sar-sum-rule-end'),
    ], 'chem-sar-sum-rule-box');
    ui.tooltip.bind(box, () => `Observed ${host.formatActivity(lo)} to ${host.formatActivity(hi)} over ` +
      `${count(data.measuredCells)} measured cells. The rule spans what was measured, end to end.` +
      (showZero ? ' The tick is zero.' : ''));
    return box;
  }

  /** MK-E: whether this series' fit was checked, and whether it held. No colour — `confidence` is null
   *  where too few cells could be cross-validated or every measured value was identical, and a traffic
   *  light would make "unchecked" read as a failing grade. */
  trustDot(matrix: SarMatrix): HTMLElement {
    const conf = matrix.confidence;
    const dot = ui.div([], 'chem-sar-sum-dot' + (!conf ? ' chem-sar-sum-dot-none' :
      conf.r2 < TRUST_R2 ? ' chem-sar-sum-dot-open' : ''));
    ui.tooltip.bind(dot, () => !conf ?
      'Unchecked: this series has fewer than four measured cells that could be held out, or every ' +
      'measured value in it is identical. Its predictions are unverified, not wrong.' :
      `R² ${conf.r2.toFixed(2)} ± ${this.host.formatActivity(conf.rmse)} predicting cells held out of ` +
      `the fit, over ${conf.n} of ${conf.total} measured — the fit ` +
      `${conf.r2 >= TRUST_R2 ? 'holds' : 'does not hold'} at R² ${TRUST_R2}. A property of the model, ` +
      'not of the data.');
    return dot;
  }

  private buildChips(data: SummaryData): HTMLElement {
    const chip = (text: string, tip: string, cls = ''): HTMLElement => {
      const el = ui.divText(text, `chem-sar-chip-badge ${cls}`.trim());
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`${count(data.trustedCells)} predictions pass the trust gate · ` +
        `${count(data.trustedStructures)} distinct structures`,
      // The second constant is concatenated rather than interpolated: the bundler folds two adjacent
      // template operands that each carry a compile-time constant and drops the first one's trailing
      // text, which loses a clause without failing the build.
      `Gate: ${MIN_SUPPORT} supporting compounds on both axes, in a series whose own R² ≥ ` +
      TRUST_R2 + '. The distinct-structure count is the actionable one — it deduplicates the same ' +
      'analog predicted in several tiers.'),
    ];
    if (data.trustedNoStructure > 0) {
      items.push(chip(`${count(data.trustedNoStructure)} predictions have a value and no structure`,
        'These pass the same trust gate as the analogs in Do next, but the core carries an attachment ' +
        'point none of the picked fragment columns fills, so no structure can be completed over it. ' +
        'Add that column to the R-group columns and they appear in Do next.', 'chem-sar-chip-partial'));
    }
    return ui.divH(items, 'chem-sar-sum-chips');
  }

  // ---- Band A′: the three answers -------------------------------------------------------------

  /** Whether the front-runner may be named alone: its margin over the runner-up has to clear the
   *  screen's own model error, or the two are one answer with two names. */
  leads(first: number | null, second: number | null, modelError: number | null): boolean {
    if (second === null)
      return true;
    if (first === null)
      return false;
    return first - second > (modelError ?? 0);
  }

  /** The substituent a tile's conclusion is about, drawn. Undefined where there is no one structure —
   *  a tile that names two groups or declines to rank has nothing to draw. */
  answerArt(smiles: string | null): HTMLElement | undefined {
    return smiles ? this.depiction(smiles, STRIP_MOL_W, STRIP_MOL_H) : undefined;
  }

  // ---- Rows -----------------------------------------------------------------------------------

  /** Canvas now, RDKit later: a synchronous pass over every structure on the tab would stall the tab
   *  switch that asked for them. */
  depiction(smiles: string | null, w: number, h: number): HTMLElement {
    const canvas = ui.canvas(w, h);
    canvas.classList.add('chem-sar-card-core');
    if (smiles)
      this.pendingPaints.push(() => paintMoleculeOnColor(canvas, smiles, w, h, CORE_BG_ARGB));
    return canvas;
  }

  flushPaints(): void {
    const paints = this.pendingPaints;
    this.pendingPaints = [];
    for (const paint of paints)
      paint();
  }

  badge(text: string, tip: string, partial = false): HTMLElement {
    const el = ui.divText(text, `chem-sar-chip-badge${partial ? ' chem-sar-chip-partial' : ''}`);
    ui.tooltip.bind(el, tip);
    return el;
  }

  /** Signed, and finer than `formatActivity`: fitted effects and pooled differences live in tenths,
   *  where one decimal renders three differently-ranked R-groups as three rows all reading "+0.1". */
  formatEffect(value: number): string {
    const abs = Math.abs(value);
    return `${value < 0 ? '−' : '+'}${abs < 10 ? abs.toFixed(2) : abs.toFixed(1)}`;
  }

  /** A measured difference in the reader's own units: a signed log difference where the scale is a
   *  log, and the fold it corresponds to where it is not. A loss is written as its reciprocal fold:
   *  10^value rounds every loss past one log to "0.0×", which destroys the number it is reporting. */
  formatDelta(value: number): string {
    if (this.host.activityIsLog)
      return this.formatEffect(value);
    const fold = Math.pow(10, Math.abs(value));
    const shown = fold < 10 ? fold.toFixed(1) : fold.toFixed(0);
    return value < 0 ? `1/${shown}×` : `${shown}×`;
  }

  /** Every clickable row on the tab. `.chem-sar-card` carries the hover, the left border and the row
   *  metric the two text lines and the value column are all sized against. */
  summaryRow(opts: {
    depiction: HTMLElement | null,
    name: string,
    badges: HTMLElement[],
    desc: string,
    /** A mark lane between the text column and the value, at a fixed width so it reads down a list. */
    mark?: HTMLElement,
    value: string,
    caption: string,
    valueTip?: string,
    faint?: boolean,
    cart?: MatrixCellRef,
    tip: string,
    onClick: () => void,
    /** Makes the chevron a toggle of its own rather than a repeat of the row's click. */
    onChevron?: () => void,
  }): HTMLElement {
    const value = ui.divText(opts.value,
      `chem-sar-sum-value${opts.faint ? ' chem-sar-sum-faint' : ''}`);
    const caption = ui.divText(opts.caption, 'chem-sar-sum-cap');
    const valueBox = ui.divV([value, caption], 'chem-sar-sum-valuebox');
    if (opts.valueTip)
      ui.tooltip.bind(valueBox, () => opts.valueTip!);
    const parts: HTMLElement[] = [];
    if (opts.depiction !== null)
      parts.push(opts.depiction);
    parts.push(ui.divV([
      ui.divH([ui.divText(opts.name, 'chem-sar-card-name'), ...opts.badges], 'chem-sar-card-title'),
      ui.divText(opts.desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body'));
    if (opts.mark !== undefined)
      parts.push(opts.mark);
    parts.push(valueBox);
    if (opts.cart !== undefined) {
      const ref = opts.cart;
      const cart = ui.iconFA('cart-plus', (e: MouseEvent) => {
        e.stopPropagation();
        this.host.addCellsToMakeList([ref], 'This cell has no structure to add.');
      }, 'Add this compound to the Make list');
      cart.classList.add('chem-sar-struct-icon', 'chem-sar-cart-icon');
      parts.push(cart);
    }
    const toggle = opts.onChevron;
    const go = toggle === undefined ? ui.iconFA('chevron-right') :
      ui.iconFA('chevron-down', (e: MouseEvent) => {
        e.stopPropagation();
        toggle();
      }, 'What is worth opening this series for');
    go.classList.add('chem-sar-sum-go');
    parts.push(go);
    const row = ui.divH(parts, 'chem-sar-card chem-sar-sum-row');
    ui.tooltip.bind(row, () => opts.tip);
    row.onclick = opts.onClick;
    return row;
  }

  /** A card with no qualifying rows still renders: a vanishing card reflows the grid and hides the
   *  fact that the analysis found nothing of that kind. */
  reason(text: string): HTMLElement[] {
    return [ui.divText(text, 'chem-sar-cp-hint')];
  }

  // ---- Band B: start here ---------------------------------------------------------------------

  /**
   * The best unmade compound the model names in this series, against the best one already measured.
   *
   * An additive FLOOR, not a ceiling on the chemistry — a cliff is exactly what this model cannot see,
   * so a gain near zero means "no further additive gain", never "no further potency".
   */

  // ---- Band C: the three cards ----------------------------------------------------------------

  // ---- Band B′: one fit over every component column ---------------------------------------------

  /** Why the whole fit cannot be read, or null when it can. */
  roleFitRefusal(data: SummaryData): string | null {
    const fit = data.roleFit;
    if (fit === null) {
      return `No component column carries two values with ${MIN_SUPPORT} or more compounds each, so ` +
        'there is nothing to compare.';
    }
    if (!fit.converged) {
      return `The fit did not settle within ${ROLE_FIT_MAX_SWEEPS} passes. Each component's offsets are ` +
        'still shrunk toward zero by a different amount and cannot be compared, so nothing is ranked.';
    }
    if (fit.cvR2 === null) {
      return 'Too few compounds to cross-validate a fit over ' +
        `${count(fit.roles.reduce((sum, role) => sum + role.fitted, 0))} component values. Unchecked, ` +
        'not wrong — nothing is ranked.';
    }
    if (fit.cvR2 < TRUST_R2) {
      return `Cross-validated R² ${fit.cvR2.toFixed(2)} — below ${TRUST_R2}, the additive reading does ` +
        'not hold on this table: knowing every component does not predict a compound the fit had not ' +
        'seen. The offsets are not ranked.';
    }
    return null;
  }

  // ---- Cores ----------------------------------------------------------------------------------

  /** Series a core leaderboard may rank at all: a matrix has to BE one core, since elsewhere a core is
   *  one row of one matrix and its fitted mean compares groupings, not scaffolds. */
  rankableCores(data: SummaryData): SeriesStat[] {
    return data.coresAreSeries ?
      data.series.filter((stat) => stat.realCells > 0 && stat.cpd >= MIN_SUPPORT) : [];
  }

  /**
   * One tier's cores, best first.
   *
   * Never two tiers at once: a folded tier holds the same compounds over a core cut one bond broader,
   * so a mixed ranking enters one compound twice under two cores and reads as two scaffolds agreeing.
   */
  topCores(data: SummaryData): SeriesStat[] {
    const tier = this.tierFilter ?? this.defaultCoreTier(data);
    if (tier === null)
      return [];
    const dir = this.host.higherIsBetter ? 1 : -1;
    return this.rankableCores(data).filter((stat) => stat.tier === tier)
      .sort((a, b) => dir * (b.typical - a.typical) || (a.matrix.id < b.matrix.id ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  // ---- Worth making ------------------------------------------------------------------------------

  // ---- Trust section --------------------------------------------------------------------------

  /**
   * What the two R² on this tab mean and what the gate does with them, then the series at each end of it.
   * The poorly-fitting list opens itself when something elsewhere on the tab links here, since that is
   * the one arrival where the reader came for the list rather than the explanation.
   */
  private buildTrustSection(data: SummaryData): HTMLElement {
    const host = this.host;
    const parts: HTMLElement[] = [
      ui.divText('Fit quality', 'chem-sar-sum-card-title'),
      // Two quantities with one name on one screen, and nothing else says they are not the same number.
      this.prose('Two different R² appear on this tab: one for the whole table, one for each series.',
        'Both are scored by predicting compounds that were left out of the fit, so they say how well ' +
        'the model predicts rather than how well it describes what it was shown. In a series, 1 means ' +
        'its substituent effects add perfectly, 0 means the model does no better than simply using that ' +
        'series\' average, and below 0 it does worse than that average.'),
      this.prose(`Trust gate: ${MIN_SUPPORT} measured compounds on both axes, in a series whose own R² ` +
        'is at least ' + TRUST_R2 + '.',
      'Predictions that fail it are still computed and still drawn in the matrix. They are left out of ' +
        'the counts above, and out of Worth making unless nothing at all clears the gate.'),
    ];
    const worst = [...data.lowR2].sort((a, b) =>
      (a.confidence!.r2 - b.confidence!.r2) || (a.id < b.id ? -1 : 1));
    if (worst.length === 0) {
      parts.push(this.prose('Every cross-validated fit holds up.',
        'Substituent effects add in each of them, so their predictions can be read as estimates.'));
    } else {
      parts.push(this.foldable(`Where the additive model does not hold · ${count(worst.length)}`,
        `${count(data.lowR2Virtual)} predicted cells inherit a fit that reproduces its own measured ` +
        'cells poorly.', this.trustList(worst), this.openTrustList));
    }
    const held = host.matrices.filter((matrix) => matrix.confidence != null &&
      matrix.confidence.r2 >= TRUST_R2)
      .sort((a, b) => (b.confidence!.r2 - a.confidence!.r2) || (a.id < b.id ? -1 : 1));
    if (held.length > 0) {
      parts.push(this.foldable(`Where it holds best · ${count(held.length)}`,
        'Best first. A prediction in one of these rests on a fit that reproduced its own measured cells.',
        this.trustList(held), false));
    }
    return ui.divV(parts, 'chem-sar-sum-trust');
  }

  /** A statement on the pane and its qualification on the hover. Method is read for the one fact it is
   *  opened for, and a paragraph per fact buries that fact in the others. */
  prose(text: string, tip: string): HTMLElement {
    const el = ui.divText(text, 'chem-sar-sum-prose');
    ui.tooltip.bind(el, () => tip);
    return el;
  }

  /** The same, in the faint style a row uses for its trailing clause. */
  hint(text: string, tip: string): HTMLElement {
    const el = ui.divText(text, 'chem-sar-cp-hint');
    ui.tooltip.bind(el, () => tip);
    return el;
  }

  /** A chevron heading that shows and hides `body`; the caller owns whatever else sits on the row. */
  foldHead(title: string, body: HTMLElement, open: boolean, cls = '',
    onToggle?: (open: boolean) => void): HTMLElement {
    const icon = ui.iconFA('chevron-right');
    icon.classList.add('chem-sar-sum-fold-icon');
    icon.classList.toggle('chem-sar-sum-fold-open', open);
    body.style.display = open ? '' : 'none';
    const head = ui.divH([icon, ui.divText(title, 'chem-sar-sum-fold-title')],
      `chem-sar-sum-fold ${cls}`.trim());
    head.onclick = () => {
      const show = body.style.display === 'none';
      body.style.display = show ? '' : 'none';
      icon.classList.toggle('chem-sar-sum-fold-open', show);
      onToggle?.(show);
    };
    return head;
  }

  /** The lead line above the fold stays, since a count is worth reading without opening anything. */
  private foldable(title: string, lead: string, body: HTMLElement, open: boolean): HTMLElement {
    return ui.divV([this.foldHead(title, body, open), ui.divText(lead, 'chem-sar-cp-hint'), body]);
  }

  private trustList(matrices: SarMatrix[]): HTMLElement {
    const host = this.host;
    const list = ui.div(matrices.slice(0, TRUST_LIST_MAX).map((matrix) => {
      const conf = matrix.confidence!;
      return this.summaryRow({
        // "Series 7" names a series without showing one. Every other list on the tab draws the
        // chemistry it is talking about, and a reader deciding whether to trust a series' predictions
        // wants to see which scaffold they are about.
        depiction: this.depiction(matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
        name: matrix.label,
        badges: [this.badge(`R² ${conf.r2.toFixed(2)} ± ${host.formatActivity(conf.rmse)}`,
          'How well this series\' model predicts a measured cell held out of its own fit, and the ' +
          'typical size of that error. Near 0 or below, substituent effects do not add up here.',
          conf.r2 < TRUST_R2)],
        desc: `${conf.n} of ${conf.total} measured cells checked`,
        value: `${count(matrix.virtualCount)}`,
        caption: 'predicted cells',
        tip: `Open ${matrix.label} in the SAR Matrix`,
        onClick: () => host.revealMatrix(matrix),
      });
    }));
    if (matrices.length > TRUST_LIST_MAX)
      list.appendChild(ui.divText(`+${matrices.length - TRUST_LIST_MAX} more`, 'chem-sar-cp-hint'));
    return list;
  }

  cellLocation(matrix: SarMatrix, ri: number, ci: number): string {
    return `${matrix.label} · ${matrix.rows[ri].label} × ${matrix.columns[ci].position}`;
  }

  openTip(matrix: SarMatrix, ri: number, ci: number): string {
    return `Open ${this.cellLocation(matrix, ri, ci)} in the SAR Matrix`;
  }

  /** The tier ranked when the reader has not picked one: the one holding the most cores, which is the
   *  cut depth this library actually recurs at. */
  defaultCoreTier(data: SummaryData): number | null {
    const tiers = this.coreTiers(data);
    if (tiers.length === 0)
      return null;
    const counts = new Map(tiers.map((tier) =>
      [tier, this.rankableCores(data).filter((stat) => stat.tier === tier).length]));
    return tiers.reduce((best, tier) => counts.get(tier)! > counts.get(best)! ? tier : best, tiers[0]);
  }

  /** The fold tiers that carry rankable cores, shallowest first. */
  coreTiers(data: SummaryData): number[] {
    return [...new Set(this.rankableCores(data).map((stat) => stat.tier))].sort((a, b) => a - b);
  }
}
