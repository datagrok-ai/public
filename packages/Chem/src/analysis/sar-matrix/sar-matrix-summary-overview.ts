/* The Summary's Overview segment, the one the analysis opens on: the answers, the component band
   and its measured swaps, and Start here. Reaches the panel only through its kit. */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {RoleSummary} from './sar-matrix-role-fit';
import {SarMatrix} from './sar-matrix-types';
import {ANALOG_LIST_MAX, MIN_SUPPORT, RGROUP_MIN_SERIES, SeriesStat, SummaryData, SwapPool, SWAP_MIN_PAIRS,
  SWAP_MIN_SERIES, TRUST_R2} from './sar-matrix-summary-data';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, count, STRIP_MOL_H, STRIP_MOL_W,
  TAB_TRANSFER} from './sar-matrix-ui-common';
import {FINDING_ROWS, GAIN_TICKS, PANE_EFFECTS, PANE_MAKING, REASON_GLYPHS, REASON_WORDS, orient, Answer,
  shortSmiles, swapAnchor, roleRanks} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

/** What this segment reads from the panel: the shared rendering helpers and the state segments share. */
export type SummaryOverviewKit = Pick<SummaryPanel, 'answerArt' | 'badge' | 'coreTiers' | 'defaultCoreTier' |
  'depiction' | 'expandTopRGroup' | 'flushPaints' | 'foldHead' | 'foldedBands' | 'formatDelta' |
  'formatEffect' | 'hint' | 'host' | 'leads' | 'rankableCores' | 'reason' | 'revealTrust' | 'roleFitRefusal' |
  'root' | 'scroller' | 'showPane' | 'summaryRow' | 'syncTransferLine' | 'topCores' | 'transferLine' |
  'trustDot'>;

export class SummaryOverview {
  constructor(private readonly kit: SummaryOverviewKit) {}

  build(data: SummaryData): HTMLElement {
    const cores = this.kit.topCores(data);
    const list = ui.divV(this.buildStartHere(data), 'chem-sar-sum-list');
    this.kit.scroller = list;
    return ui.divV([
      this.coverageBar(data),
      this.setupLine(data),
      this.buildAnswers(data, cores),
      this.buildListHeader(data),
      list,
      this.buildTotals(data),
    ], 'chem-sar-sum-overview');
  }

  /** The five reason-lane glyphs printed once, directly above the lane they label and always on
   *  screen — a mark whose legend is a scroll away is a legend hunt. */
  private buildListHeader(data: SummaryData): HTMLElement {
    const label = ui.divText('Start here', 'chem-sar-sum-card-title');
    if (!this.lanesFit(data))
      return ui.divH([label], 'chem-sar-sum-list-head');
    const keys = REASON_GLYPHS.map((glyph, i) =>
      ui.divH([ui.divText(glyph, 'chem-sar-sum-slot chem-sar-sum-slot-on'),
        ui.divText(REASON_WORDS[i], 'chem-sar-sum-legend-word')], 'chem-sar-sum-legend-key'));
    return ui.divH([label, ui.divH(keys, 'chem-sar-sum-legend')], 'chem-sar-sum-list-head');
  }

  /** Three totals that each open the segment answering them. Magnitudes that only describe the run
   *  stay on Method: a landing screen carries decisions, not census figures. */
  private buildTotals(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const ranked = data.analogs.all.length;
    const thin = data.analogsThin.all.length;
    const gated = data.analogs.all.concat(data.analogsThin.all);
    const shown = gated.length > 0 ? gated : data.analogsAny.all;
    const structures = new Set(shown.map((r) => r.key)).size;
    // A capped list printed as a total states the cap as the finding, so the cap is named wherever the
    // length is.
    const capped = gated.length > 0 && data.analogOverflow > 0 ? `top ${count(ranked + thin)} of ` +
      `${count(ranked + thin + data.analogOverflow)}` : count(shown.length);
    const making = ui.divText(shown.length === 0 ? 'Nothing is predicted above what is already made →' :
      gated.length === 0 ? `${capped} best candidates, none past the trust gate →` :
        `${capped} worth making · ${count(structures)} distinct structures →`, 'chem-sar-sum-total');
    making.onclick = () => this.kit.showPane(PANE_MAKING);
    ui.tooltip.bind(making, () => 'Predicted analogs this dataset has no row for, ranked on the gain ' +
      'they buy over the best compound their own series has already made.' + (data.analogOverflow === 0 ?
      '' : ` The ranked and the thin list hold ${count(ANALOG_LIST_MAX)} rows each: ` +
      `${count(data.analogOverflow)} further structures cleared the same gate and are in neither.`));

    const trust = ui.divText(`${count(data.fitHolds)} fits hold · ${count(data.unchecked)} unchecked · ` +
      `${count(data.lowR2.length)} do not →`, 'chem-sar-sum-total');
    trust.onclick = () => this.kit.revealTrust();
    ui.tooltip.bind(trust, 'A verdict on each series\' fit, not a partition of your library — one ' +
      'compound sits in several series and is routinely on both sides. Unchecked means unverified, ' +
      'not wrong.');

    const parts = [making, trust];
    if (host.unscalableCount > 0) {
      const bad = ui.divText(`${count(host.unscalableCount)} values cannot be scaled by ` +
        `${host.scalingLabel}`, 'chem-sar-sum-total chem-sar-chip-partial');
      ui.tooltip.bind(bad, () => `${host.unscalableCount} of the "${host.activityColumnName}" values ` +
        `were assayed but cannot be scaled by ${host.scalingLabel}, which needs positive numbers. ` +
        'They are excluded from every matrix — set Scaling to "none" to use them as they are.');
      parts.push(bad);
    }
    return ui.divH(parts, 'chem-sar-sum-totals');
  }

  /**
   * MK-cov: what fraction of the table the analysis sees, in four segments in a fixed order, over a
   * label line naming them in the same order — so the label is the legend.
   */
  private coverageBar(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const total = Math.max(1, host.hostRowCount);
    // The two ways a row misses out are split because they call for different things: one is a plate
    // to run, the other a library too diverse for any core to recur. Lumped, a sparse assay and a
    // diverse library read the same.
    const unplaced = Math.max(0, host.assayedCount - data.compounds - host.unscalableCount);
    const never = Math.max(0, total - host.assayedCount - data.untested);
    const segs = [
      {n: data.compounds, word: 'measured', cls: 'chem-sar-sum-cov-real'},
      {n: data.untested, word: 'held, never assayed', cls: 'chem-sar-sum-cov-held'},
      {n: host.unscalableCount, word: 'assayed, cannot be scaled', cls: 'chem-sar-sum-cov-bad'},
      {n: unplaced, word: 'assayed, no analog to pair with', cls: 'chem-sar-sum-cov-out'},
      {n: never, word: 'never assayed', cls: 'chem-sar-sum-cov-none'},
      // A zero-count segment is omitted rather than drawn at the minimum width, which would state a
      // quantity that is not there.
    ].filter((s) => s.n > 0);
    const bar = ui.divH(segs.map((s) => {
      const el = ui.div([], `chem-sar-sum-cov-seg ${s.cls}`);
      el.style.flexGrow = `${s.n}`;
      return el;
    }), 'chem-sar-sum-cov');
    const label = ui.divText(segs.map((s) => `${count(s.n)} ${s.word}`).join(' · '),
      'chem-sar-sum-cov-label');
    const box = ui.divV([bar, label], 'chem-sar-sum-cov-box');
    ui.tooltip.bind(box, () => `${count(total)} table rows. Rows, not compounds: the denominator is the ` +
      'table\'s raw row count, so this ignores any filter on the source grid. A segment narrower than ' +
      'a couple of pixels is drawn at that floor and stops being proportional, which is why every ' +
      'count is printed.');
    return box;
  }

  /** MK-F: how many of a series' own typical prediction errors a predicted gain is worth. Unfilled
   *  slots are drawn, not omitted, so "one of four" never looks like "one of one". */
  private tickRun(gain: number, rmse: number | null): HTMLElement {
    const filled = rmse !== null && rmse > 0 ?
      Math.max(0, Math.min(GAIN_TICKS, Math.floor(gain / rmse))) : 0;
    const ticks: HTMLElement[] = [];
    for (let i = 0; i < GAIN_TICKS; i++)
      ticks.push(ui.div([], `chem-sar-sum-tick${i < filled ? ' chem-sar-sum-tick-on' : ''}`));
    const run = ui.divH(ticks, 'chem-sar-sum-ticks');
    ui.tooltip.bind(run, () => rmse === null || rmse <= 0 ?
      'This series has no cross-validated error to measure the gain against.' :
      `The predicted gain is ${(gain / rmse).toFixed(1)} times this series' own leave-one-out error. ` +
      'Under one error the model cannot tell this analog from the best compound already made.');
    return run;
  }

  /** MK-B: the five fixed reasons a series is worth opening, won or not won, always in the same order
   *  — the lane's shape is what a reader compares down the column, which five variable-length phrases
   *  can never be. */
  private reasonLane(reasons: number[]): HTMLElement {
    const won = new Set(reasons);
    const slots = REASON_GLYPHS.map((glyph, i) => ui.divText(won.has(i) ? glyph : '',
      `chem-sar-sum-slot${won.has(i) ? ' chem-sar-sum-slot-on' : ''}`));
    const lane = ui.divH(slots, 'chem-sar-sum-lane');
    ui.tooltip.bind(lane, () => reasons.map((i) => REASON_WORDS[i]).join(' · '));
    return lane;
  }

  /**
   * What the analysis answers, in words: which R-group, which core, and which swap buys the most.
   *
   * Every number here is already on a card below, so a reader who stops at this band has the finding
   * and one who does not finds the evidence under it. Each tile carries its own negative, because a
   * leaderboard with no stated limit is read as a ranking of the chemistry rather than of what was made.
   */
  private buildAnswers(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const findings = this.buildFindings(data);
    if (findings !== null)
      return findings;
    return ui.divH([
      this.answerTile(`Best ${data.axisRole ?? 'R-group'}`, this.rgroupAnswer(data), 'rgroup',
        () => this.kit.expandTopRGroup = true),
      this.answerTile(`Best ${data.coreRole ?? 'core'}`, this.coreAnswer(data, cores), 'cores'),
      this.answerTile('What improves potency', this.potencyAnswer(data), 'swaps'),
    ], 'chem-sar-sum-answers');
  }

  /**
   * One row per component column and one for the best measured swap, in a single full-width band.
   *
   * Every component, not two of them: tiles for the axis and the core left the rest to be hunted for,
   * which sent readers to the axis switch on Method — a rebuild of every matrix — to see a ranking that
   * was already computed. One band rather than a strip beside a tile, because the tile's column was too
   * narrow for a structure and fell back to a truncated SMILES, which reads identically for every analog
   * of one scaffold.
   *
   * Null outside fragment-columns mode and wherever the fit is refused — a fragmented substituent label
   * is local to its series, so there is no component to name and the three tiles stay.
   */
  private buildFindings(data: SummaryData): HTMLElement | null {
    const fit = data.roleFit;
    if (fit === null || this.kit.roleFitRefusal(data) !== null)
      return null;
    const rows = fit.roles.map((role, index) => this.componentRow(role, index));
    const body: HTMLElement[] = [
      ui.divText(`spans = the range of ${this.kit.host.activityColumnName} across that component's ` +
        `values · offsets are against the library mean of ${this.kit.host.formatActivity(fit.mean)}`,
      'chem-sar-cp-hint'),
      ...rows.slice(0, FINDING_ROWS),
    ];
    // A six-component table would otherwise push the swap row, and the band under it, off the screen the
    // reader lands on. The ranking is by spread, so what folds away is what moves the endpoint least.
    if (rows.length > FINDING_ROWS)
      body.push(this.moreComponents(rows.slice(FINDING_ROWS)));
    // One row per component, not one for whichever component the matrix columns happened to enumerate:
    // a matched pair is a grouping on the other components, so every component's pairs are in hand.
    const swapRows = fit.roles.map((role) => this.swapRow(data, role.name));
    const withPairs = fit.roles.filter((role) => (data.swapsByRole.get(role.name) ?? []).length > 0);
    // Each half folds to its own headline, so collapsing one to read the other never costs the answer.
    const lead = fit.roles[0];
    return ui.divV([
      this.bandGroup('comp', `What to change · ${count(fit.roles.length)} components`,
        `${lead.name} leads, spans ${lead.spread.toFixed(2)}`, ui.divV(body)),
      this.bandGroup('swap', `Best measured swap · ${count(fit.roles.length)} components`,
        withPairs.length === 0 ? 'none pooled' :
          `${withPairs[0].name} ${this.swapHeadline(data, withPairs[0].name)}`, ui.divV(swapRows)),
    ], 'chem-sar-sum-comp');
  }

  /**
   * Which column is the core, which one the matrix columns enumerate, and which fold into the row.
   *
   * The band below names components without saying what each one is to the analysis, and the three are
   * not interchangeable: the core is the scaffold every row is drawn from, the one across the columns is
   * the only one whose values sit side by side in the grid, and the rest are part of a row's identity.
   */
  private setupLine(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    if (data.axisRole === null) {
      return this.kit.hint('Cores and substituents found by cutting the molecules — no component columns ' +
        'were given.', 'Every core here is a fragment the cutting produced, so its substituent labels ' +
        'are local to the series they came from and mean nothing across series.');
    }
    const folded = host.roleColumns.filter((name) => name !== data.axisRole);
    const parts = [`core: ${data.coreRole}`, `across the matrix columns: ${data.axisRole}`];
    if (folded.length > 0)
      parts.push(`folded into the row: ${folded.join(', ')}`);
    return this.kit.hint(parts.join(' · '),
      'The core is the scaffold every row is drawn from. The component across the columns is the one ' +
      'whose values sit side by side in the SAR Matrix grid; the rest are part of what identifies a ' +
      'row. Every component is ranked and has its swaps pooled either way — the choice only decides ' +
      'the grid\'s layout.');
  }

  /** One foldable half of the landing band: a heading that states its own answer when shut, so the band
   *  can be narrowed to the half the reader wants without hiding what the other half concluded. */
  private bandGroup(key: string, title: string, headline: string, body: HTMLElement): HTMLElement {
    const open = !this.kit.foldedBands.has(key);
    const lead = ui.divText(headline, 'chem-sar-cp-hint');
    lead.style.display = open ? 'none' : '';
    const head = this.kit.foldHead(title, body, open, 'chem-sar-sum-band-head', (show) => {
      lead.style.display = show ? 'none' : '';
      if (show)
        this.kit.foldedBands.delete(key);
      else
        this.kit.foldedBands.add(key);
    });
    head.appendChild(lead);
    return ui.divV([head, body]);
  }

  /** The best pooled swap of one component, or of whichever component the matrix columns enumerated
   *  when the components are not columns and there is only the one pool. */
  private topSwap(data: SummaryData, role?: string): SwapPool | undefined {
    return role === undefined ? data.swaps[0] : (data.swapsByRole.get(role) ?? [])[0];
  }

  /** The swap's conclusion in one phrase, for the group heading when it is shut. */
  private swapHeadline(data: SummaryData, role?: string): string {
    const top = this.topSwap(data, role);
    if (top === undefined)
      return 'none pooled';
    const {worst, widest} = orient(top);
    return worst > 0 ? `≥ ${this.kit.formatDelta(worst)} in ${top.n} pairs` :
      `${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)}`;
  }

  /** One component: how much it moves the endpoint, and the best value of it.
   *
   * The range is written out rather than drawn against the widest component. A bar whose scale is the
   * other rows reads as a fraction of something, and with three components within 0.03 of one another it
   * draws a ranking the numbers do not support. */
  private componentRow(role: RoleSummary, index: number): HTMLElement {
    const top = roleRanks(role) ? role.levels[0] : undefined;
    const parts: HTMLElement[] = [
      ui.divText(role.name, 'chem-sar-sum-comp-name'),
      ui.divText(`spans ${role.spread.toFixed(2)}`, 'chem-sar-sum-comp-spread'),
    ];
    if (top === undefined)
      parts.push(ui.divH([ui.divText('not ranked', 'chem-sar-cp-hint')], 'chem-sar-sum-comp-slot'));
    else {
      parts.push(ui.divH([top.value === '' ? this.kit.hint('nothing at this position',
        `The compounds that leave ${role.name} empty score best. The column is blank for them, so there ` +
        'is no group to draw — in a decomposition that is hydrogen, in a list of named components it ' +
        'means nobody recorded one.') :
        this.kit.depiction(top.value, STRIP_MOL_W, STRIP_MOL_H)], 'chem-sar-sum-comp-slot'));
      parts.push(ui.divText(this.kit.formatEffect(top.coef), 'chem-sar-sum-comp-best'));
      parts.push(ui.divText(`over ${count(top.n)} compounds`, 'chem-sar-cp-hint'));
    }
    const row = ui.divH(parts, 'chem-sar-sum-comp-row');
    ui.tooltip.bind(row, () => `${role.name} moves ${this.kit.host.activityColumnName} by ` +
      `${role.spread.toFixed(2)} across its values. Click for its full ranking.`);
    row.onclick = () => this.kit.showPane(PANE_EFFECTS, `role-${index}`);
    return row;
  }

  /** The components that move the endpoint least, folded away behind their own count. */
  private moreComponents(rows: HTMLElement[]): HTMLElement {
    const body = ui.divV(rows);
    const head = this.kit.foldHead(`${rows.length} more`, body, false, 'chem-sar-sum-comp-more');
    return ui.divV([head, body]);
  }

  /**
   * The best measured swap, on the same row grammar as the components above it. Structures and two
   * numbers rather than a sentence, which has no natural length and clips wherever the band is narrow.
   */
  private swapRow(data: SummaryData, role?: string): HTMLElement {
    // Named like a component row, because that is what it reports on: the pairs are grouped on the other
    // components, so each component has its own and the row has to say which.
    const parts: HTMLElement[] = [
      ui.divText(role ?? 'Best swap', 'chem-sar-sum-comp-name'),
      ui.divText(role === undefined ? '' : 'swapped', 'chem-sar-sum-comp-spread'),
    ];
    const top = this.topSwap(data, role);
    const only = role === undefined ? '' :
      ` Pairs alike in every component but ${role}, pooled across the rest of the table.`;
    let tip: string;
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      parts.push(this.pairArt(from, to, 'chem-sar-sum-comp-slot chem-sar-sum-comp-pair'));
      parts.push(ui.divText(worst > 0 ? `≥ ${this.kit.formatDelta(worst)}` :
        `mean ${this.kit.formatDelta(mean)}`, 'chem-sar-sum-comp-best'));
      parts.push(ui.divText(`${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)} · ` +
        `${count(top.n)} pairs · ${count(top.roots.size)} ${role === undefined ? 'series' : 'contexts'}`,
      'chem-sar-cp-hint'));
      tip = (worst > 0 ? 'This swap improved potency in every pair it was measured in.' :
        'No swap improved potency in every pair it was measured in; this one has the best floor.') +
        ' Both compounds were made and measured; nothing here is fitted or predicted.' + only +
        ' Click for the full pool.';
    } else {
      const best = data.swapBest;
      parts.push(ui.divH([ui.divText(data.swapCandidates === 0 && role === undefined ?
        'no row carries two measured values' :
        `nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} ` +
        `${role === undefined ? 'series' : 'contexts'}`, 'chem-sar-cp-hint')],
      'chem-sar-sum-comp-slot'));
      if (role === undefined && best?.best != null) {
        parts.push(ui.divText(`largest single move ${this.kit.formatDelta(Math.abs(best.best.delta))} ` +
          `in ${best.best.matrix.label}`, 'chem-sar-cp-hint'));
      }
      tip = 'A swap is pooled only where the same substitution was measured against several different ' +
        'backgrounds, so one pair cannot carry it.' + only + ' Click for what was rejected.';
    }
    const row = ui.divH(parts, 'chem-sar-sum-comp-row');
    ui.tooltip.bind(row, () => tip);
    row.onclick = () => this.kit.showPane(PANE_EFFECTS, role === undefined ? 'swaps' : swapAnchor(role));
    return row;
  }

  /** One conclusion: the name, the number, its bar against the model error, and one clause. The
   *  evidence is a segment away, so nothing here has to carry it. */
  private answerTile(title: string, text: Answer, anchor: string, onOpen?: () => void): HTMLElement {
    const go = ui.iconFA('chevron-right');
    go.classList.add('chem-sar-sum-go');
    const body: HTMLElement[] = [
      ui.divH([ui.divText(title, 'chem-sar-sum-card-title'), go], 'chem-sar-sum-answer-head'),
    ];
    if (text.art !== undefined)
      body.push(text.art);
    body.push(ui.divH([ui.divText(text.answer, 'chem-sar-sum-answer-value')], 'chem-sar-sum-answer-line'));
    // Not `chem-sar-card-desc`: this is a sentence, and that class clips to one line.
    body.push(ui.divText(text.negative, 'chem-sar-cp-hint'));
    const tile = ui.divV(body, 'chem-sar-sum-answer');
    // The answer itself, not only the routing: the value line clips, and a conclusion that fits nowhere
    // on the tile has to be readable somewhere.
    ui.tooltip.bind(tile, () => `${text.answer} — ${text.negative}. The evidence behind this answer is ` +
      'on the Effects segment; click to go to it.');
    tile.onclick = () => {
      onOpen?.();
      this.kit.showPane(PANE_EFFECTS, anchor);
    };
    return tile;
  }

  private rgroupAnswer(data: SummaryData): Answer {
    if (!this.kit.host.activityIsLog) {
      // There is no mark for "not answerable", so this branch stays a sentence.
      return {answer: 'Not answerable on a raw scale',
        negative: 'a fitted effect is a difference in assay units and says nothing about how large ' +
          'the change is — set Scaling to lg or −lg, or declare the column higher-is-better if it is ' +
          'already a pIC50'};
    }
    const loser = data.rgroupLosers[0];
    const negative = loser === undefined ?
      'over each series\' own most-common substituent, not over your whole library' :
      `stop making ${shortSmiles(loser.subst)} — last in ${loser.k} of ${loser.tried} tried`;
    const rows = data.rgroups;
    if (rows.length === 0) {
      return {answer: `None in ${RGROUP_MIN_SERIES} series`,
        negative: 'no substituent comes first at its position in that many series whose additive fit holds'};
    }
    const top = rows[0];
    const partner = data.axisRole === null ? 'series' : 'core';
    if (!this.kit.leads(top.magnitude, rows[1]?.magnitude ?? null, data.modelError)) {
      return {answer: `Two lead equally`,
        negative: `the ${partner} decides — ${top.k} of ${top.tried} against ${rows[1].k} of ` +
          `${rows[1].tried}`,
        art: ui.divH([this.kit.depiction(top.subst, STRIP_MOL_W, STRIP_MOL_H),
          ui.divText('or', 'chem-sar-sum-comp-arrow'),
          this.kit.depiction(rows[1].subst, STRIP_MOL_W, STRIP_MOL_H)], 'chem-sar-sum-comp-pair')};
    }
    return {answer: (top.magnitude === null ? '' : `${this.kit.formatEffect(top.magnitude)} · `) +
      `first on ${top.k} of ${top.tried}`, negative, art: this.kit.answerArt(top.subst)};
  }

  private coreAnswer(data: SummaryData, cores: SeriesStat[]): Answer {
    const host = this.kit.host;
    if (!data.coresAreSeries) {
      // A core recurs only inside its own lineage, so there is no ranking a mark could carry.
      return {answer: 'Not comparable across series',
        // Two lines at tile width; a third is clipped, and a clipped sentence is what the tile is for.
        negative: 'a core recurs only inside its own lineage — nothing to pool across series'};
    }
    if (cores.length === 0)
      return {answer: `None with ${MIN_SUPPORT} compounds`, negative: 'mean of what was made'};
    // Which cut depth was ranked, because only cores cut alike are comparable and the card lets the
    // reader change it.
    const depth = this.kit.coreTiers(data).length > 1 ? `among the L${cores[0].tier} cores; ` : '';
    const negative = `${depth}mean of what was made — a core only ever paired with good caps scores high`;
    const dir = host.higherIsBetter ? 1 : -1;
    const top = cores[0];
    if (!this.kit.leads(dir * top.typical, cores.length > 1 ? dir * cores[1].typical : null, data.modelError)) {
      return {answer: `No core runs clear`,
        negative: `${depth}the top ${cores.length} sit inside ` +
          `${host.formatActivity(data.modelError ?? 0)} model error of each other`};
    }
    return {answer: `${top.matrix.label} · typical ${host.formatActivity(top.typical)} over ` +
      `${count(top.cpd)} cpd`, negative,
    art: this.kit.answerArt(top.matrix.rows[0]?.keySmiles ?? null)};
  }

  private potencyAnswer(data: SummaryData): Answer {
    const negative = 'both compounds measured, nothing fitted — one R-group swapped, the rest identical';
    const top = data.swaps[0];
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      const art = this.pairArt(from, to, 'chem-sar-sum-comp-pair');
      if (worst > 0)
        return {answer: `≥ ${this.kit.formatDelta(worst)} in ${top.n} pairs`, negative, art};
      // The pool is ranked on its floor, and a negative floor under this heading would state the
      // reverse of the question: the swap lost in at least one pair it was measured in.
      return {answer: `${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)} over ${top.n} pairs`,
        negative: `no swap improves potency in every pair it was measured in — this one is the best ` +
          `floor; mean ${this.kit.formatDelta(mean)}`, art};
    }
    // Nothing pooled: what the gate rejected is the finding, and there is nothing to draw.
    const gate = data.swapCandidates === 0 ? 'No row carries two measured substituents' :
      `Nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} series`;
    const best = data.swapBest;
    if (best === null || best.best === null)
      return {answer: gate, negative};
    return {answer: gate, negative: `the largest single measured move is ` +
      `${this.kit.formatDelta(Math.abs(best.best.delta))} in ${best.best.matrix.label}`};
  }

  /** Two structures and the arrow between them: the shape every swap is drawn in. */
  private pairArt(from: string, to: string, cls: string): HTMLElement {
    return ui.divH([this.kit.depiction(from, STRIP_MOL_W, STRIP_MOL_H),
      ui.divText('→', 'chem-sar-sum-comp-arrow'),
      this.kit.depiction(to, STRIP_MOL_W, STRIP_MOL_H)], cls);
  }

  /**
   * Whether the reason lane is drawn on this render — the header's key, the per-row lane and the text
   * fallback all read this one answer, so the mark can never lose its legend or the phrases with it.
   *
   * A lane that wins every slot carries no information, so with one series the phrases come back as
   * text: a constant mark is decoration. Below the narrow breakpoint five 14px slots and a readable row
   * cannot both fit, and `applySize` rebuilds the pane when that breakpoint is crossed.
   */
  private lanesFit(data: SummaryData): boolean {
    return data.startHere.length > 1 && !this.kit.root.classList.contains('chem-sar-sum-xnarrow');
  }

  private buildStartHere(data: SummaryData): HTMLElement[] {
    const host = this.kit.host;
    const parts: HTMLElement[] = [];
    if (data.startHere.length === 0)
      parts.push(...this.kit.reason('No series holds a measured compound.'));
    const lanes = this.lanesFit(data);
    for (const row of data.startHere) {
      const stat = row.primary;
      const {matrix} = stat;
      const rowKey = stat.best !== null ? matrix.rows[stat.best.ri].keySmiles : matrix.rows[0]?.keySmiles;
      const neighbours = stat.best === null ? 0 :
        host.observedNeighbours(matrix, stat.best.ri, stat.best.ci);
      const expand = ui.divV([]);
      expand.style.display = 'none';
      let filled = false;
      parts.push(this.kit.summaryRow({
        depiction: this.kit.depiction(rowKey ?? null, CARD_CORE_W, CARD_CORE_H),
        name: matrix.label,
        badges: [
          this.kit.badge(`L${stat.tier}`, stat.tier === 1 ?
            'A leaf series: no finer series sit under this one.' :
            `A coarser series, holding the compounds of the L${stat.tier - 1} series below it, whose ` +
            'cores agree one further cut deeper.'),
          this.kit.trustDot(matrix),
          this.gainBadge(stat),
        ],
        mark: lanes ? this.reasonLane(row.reasons) : undefined,
        // Cells the decomposition cannot express are never predicted and never offered, so they are no
        // part of the denominator — but only a matrix whose assembly marked them can say so.
        desc: `${count(stat.cpd)} cpd · ${count(stat.realCells)} of ` +
          `${count(stat.totalCells - stat.impossibleCells)} ` +
          `${stat.impossibleCells > 0 ? 'makeable cells' : 'cells'} measured · ` +
          `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}` +
          (stat.realCells ? ` · typical ${host.formatActivity(stat.typical)}` : ''),
        value: stat.best === null ? '—' : host.formatActivity(stat.best.value),
        caption: 'best',
        valueTip: stat.best === null ? 'This series holds no measured compound.' :
          `${host.formatActivity(stat.best.value)} measured. ${neighbours} measured cells share this ` +
          'core or this substituent.',
        tip: `Open ${matrix.label} in the SAR Matrix. "Typical" is this series' fitted mean — the mean ` +
          'of what was made, count-weighted, so a series whose chemists only made their good ideas ' +
          'scores high and one heavily resampled core dominates it. The gap between typical and best ' +
          'is the lucky-singleton signal.' + (stat.impossibleCells > 0 ? '' :
          ' No cell of this series was marked unmakeable, which can mean the decomposition never ' +
          'assessed them rather than that the grid is complete.'),
        onClick: () => stat.best === null ? host.revealMatrix(matrix) :
          host.revealCell(matrix, stat.best.ri, stat.best.ci),
        onChevron: () => {
          if (!filled) {
            filled = true;
            for (const part of this.startExpandParts(stat))
              expand.appendChild(part);
            // Rows added after the paint tick has run, so these depictions need one of their own.
            this.kit.flushPaints();
          }
          expand.style.display = expand.style.display === 'none' ? '' : 'none';
        },
      }));
      if (!lanes)
        parts.push(ui.divText(row.reasons.map((i) => REASON_WORDS[i]).join(' · '), 'chem-sar-cp-hint'));
      parts.push(expand);
    }
    const transfers = ui.divText('', 'chem-sar-cp-hint');
    this.kit.transferLine = transfers;
    this.kit.syncTransferLine();
    transfers.classList.add('chem-sar-sum-click');
    transfers.onclick = () => host.showTab(TAB_TRANSFER);
    parts.push(transfers);
    return parts;
  }

  /**
   * Which gate left this series with nothing to aim at. The gate is one condition — a predicted cell
   * needs MIN_SUPPORT measured compounds on both axes AND a series whose leave-one-out R² holds — but
   * it fails in ways that call for different things: a setting to change, more compounds, or nothing
   * at all because the grid is already full.
   */
  private noGainReason(stat: SeriesStat): {label: string, tip: string} {
    const r2 = stat.matrix.confidence?.r2 ?? null;
    if (!this.kit.host.predictVirtual) {
      return {label: 'predictions off', tip: '"Predict virtual analogs" is off, so this series has no ' +
        'predicted cells to compare against. Turn it on in the viewer\'s properties.'};
    }
    if (stat.virtualCells === 0) {
      return {label: 'nothing unfilled', tip: 'Every cell this series can express is already measured ' +
        'or cannot be made, so there is no unfilled cell to predict.'};
    }
    if (r2 === null) {
      return {label: 'fit unchecked', tip: `This series carries ${count(stat.virtualCells)} predicted ` +
        'cells, but its additive fit was never cross-validated — too few measured cells to leave one ' +
        'out — so none of them can be trusted. Unchecked, not wrong.'};
    }
    if (r2 < TRUST_R2) {
      return {label: 'fit fails', tip: `This series carries ${count(stat.virtualCells)} predicted ` +
        `cells, but its own R² is ${r2.toFixed(2)}, below ${TRUST_R2}: substituent effects ` +
        'do not add up here, so its predictions are arithmetic rather than estimates.'};
    }
    return {label: 'thin support', tip: `All ${count(stat.virtualCells)} predicted cells here rest on ` +
      `fewer than ${MIN_SUPPORT} measured compounds on one of their two axes. The fit holds; the ` +
      'evidence behind each individual cell does not.'};
  }

  private gainBadge(stat: SeriesStat): HTMLElement {
    const host = this.kit.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const reach = stat.bestVirtual !== null && stat.best !== null ?
      dir * (stat.bestVirtual.value - stat.best.value) : null;
    let tip = 'The best unfilled cell the model names, against the best compound already measured. ' +
      'Negative means the best compound already beats anything the additive model can build from this ' +
      'series\' parts — what a cliff looks like. Near zero means no further ADDITIVE gain, not that ' +
      'the series is exhausted.';
    let label = `${this.kit.formatEffect(reach ?? 0)} vs best made`;
    if (reach === null)
      ({label, tip} = this.noGainReason(stat));
    const badge = this.kit.badge(label, tip, reach === null);
    if (stat.bestVirtual !== null) {
      badge.classList.add('chem-sar-sum-click');
      badge.onclick = (e: MouseEvent) => {
        e.stopPropagation();
        host.revealCell(stat.matrix, stat.bestVirtual!.ri, stat.bestVirtual!.ci);
      };
    }
    // The one thing the number itself cannot say: whether +0.4 is four of this series' own prediction
    // errors or a quarter of one.
    const rmse = stat.matrix.confidence?.rmse ?? null;
    if (reach === null || rmse === null || !host.activityIsLog)
      return badge;
    return ui.divH([badge, this.tickRun(reach, rmse)], 'chem-sar-sum-gain');
  }

  private startExpandParts(stat: SeriesStat): HTMLElement[] {
    const host = this.kit.host;
    const {matrix} = stat;
    const parts: HTMLElement[] = [];

    if (stat.colRange !== null && stat.rowRange !== null) {
      const story = stat.colRange > stat.rowRange ? 'this is a column story' : 'this is a row story';
      const spread = ui.divText(`Spread: substituent choice spans ${stat.colRange.toFixed(2)} · core ` +
        `choice spans ${stat.rowRange.toFixed(2)} — ${story}`, 'chem-sar-cp-hint');
      ui.tooltip.bind(spread, () => 'Ranges of this series\' own fitted effects, over substituents ' +
        `measured on at least ${MIN_SUPPORT} cores and cores measured at at least ${MIN_SUPPORT} ` +
        'substituents — one floor on both sides, or the looser side of the comparison would come out ' +
        'wider on noise alone. Both are centred on the same fit, so they compare to each other — and ' +
        'to no other series\' numbers.');
      parts.push(spread);
    }

    if (stat.bestRow !== null) {
      const row = matrix.rows[stat.bestRow.ri];
      const folded = Object.values(row.foldedValues);
      const oneStructure = matrix.rows.length > 1 &&
        new Set(matrix.rows.map((r) => r.keySmiles)).size === 1;
      const heading = folded.length === 0 ? 'Best core' :
        oneStructure ? 'Best R-group combination' : 'Best core + fixed R-groups';
      // With R-group columns every row can be labelled "Row N" and draw the same structure, so what
      // distinguishes the winning row is its folded substituents, not its picture.
      const art = folded.length > 0 && oneStructure ?
        ui.divH(folded.slice(0, 3).map((subst) =>
          this.kit.depiction(subst, BENEFIT_MOL_W, BENEFIT_MOL_H)), 'chem-sar-sum-pair') :
        this.kit.depiction(row.keySmiles, CARD_CORE_W, CARD_CORE_H);
      const text = ui.divV([
        ui.divText(`${heading}: ${row.label} · ${this.kit.formatEffect(stat.bestRow.effect)} ` +
          `(n = ${stat.bestRow.n})`, 'chem-sar-sum-prose'),
        ui.divText('Best on average across the substituents it was actually measured at — not the ' +
          'single best compound, and not comparable to another series.', 'chem-sar-cp-hint'),
      ]);
      const block = ui.divH([art, text], 'chem-sar-sum-detail');
      block.classList.add('chem-sar-sum-click');
      const ci = this.bestMeasuredCol(matrix, stat.bestRow.ri);
      if (ci >= 0)
        block.onclick = () => host.revealCell(matrix, stat.bestRow!.ri, ci);
      parts.push(block);
    }

    // The one place on this screen where non-additivity is a finding rather than a trust problem, and
    // where the instruction is the opposite of the leaderboards': hold the pair, vary elsewhere.
    const conf = matrix.confidence;
    const dir = host.higherIsBetter ? 1 : -1;
    // `?? null`, not a strict read: a matrix carried in by a project or layout was serialized by
    // whichever Chem wrote it, and one without these keys yields undefined, which `!== null` passes.
    const outlier = (conf ? (dir > 0 ? conf.hi : conf.lo) : null) ?? null;
    if (conf && outlier !== null && Math.abs(outlier.residual) > conf.rmse &&
      outlier.ri < matrix.rows.length && outlier.ci < matrix.columns.length) {
      const beat = ui.divText(`${matrix.rows[outlier.ri].label} × ` +
        `${shortSmiles(matrix.columns[outlier.ci].substSmiles)} beats its additive expectation by ` +
        `${this.kit.formatEffect(dir * outlier.residual)}, against this series' own ` +
        `±${host.formatActivity(conf.rmse)} — hold that pair and vary elsewhere.`, 'chem-sar-cp-hint');
      beat.classList.add('chem-sar-sum-click');
      ui.tooltip.bind(beat, 'Out-of-sample: the cell was held out and predicted from the rest, so ' +
        'the model could not pull itself toward it. It beats the additive sum in this series\' own ' +
        'units; it does not say the mechanism is understood.');
      beat.onclick = () => host.revealCell(matrix, outlier.ri, outlier.ci);
      parts.push(beat);
    }

    const position = matrix.positions[0] ?? '';
    const reference = matrix.refValues[position];
    const refLine = ui.divText(reference ? `Varies ${position} · reference substituent:` :
      `Varies ${position} · no reference substituent recorded`, 'chem-sar-cp-hint');
    ui.tooltip.bind(refLine, 'The position this series explores, and the most frequently observed ' +
      'substituent at it — the comparator its fitted effects are quoted against. Position labels come ' +
      'from each series\' own decomposition, so R1 here and R1 elsewhere are unrelated.');
    parts.push(reference ?
      ui.divH([refLine, this.kit.depiction(reference, BENEFIT_MOL_W, BENEFIT_MOL_H)], 'chem-sar-sum-detail') :
      refLine);
    return parts;
  }

  /** The most potent measured cell of one row, so an expand lands on a compound. */
  private bestMeasuredCol(matrix: SarMatrix, ri: number): number {
    const dir = this.kit.host.higherIsBetter ? 1 : -1;
    let ci = -1;
    for (let c = 0; c < matrix.columns.length; c++) {
      const cell = matrix.cells[ri][c];
      if (cell.kind === 'real' && cell.value !== null &&
        (ci < 0 || dir * cell.value > dir * matrix.cells[ri][ci].value!))
        ci = c;
    }
    return ci;
  }
}
