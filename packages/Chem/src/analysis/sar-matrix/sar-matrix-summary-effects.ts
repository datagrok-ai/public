/* The Summary's Effects segment: a tab per component with its full fitted ranking, and a last tab
   of what was counted rather than fitted inside each series. Reaches the panel only through its kit. */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {cachedAtomCount} from './sar-matrix-decompose';
import {RoleLevel, RoleSummary} from './sar-matrix-role-fit';
import {bestMeasuredRow, LOSER_ROWS, MIN_SUPPORT, RGroupRow, SeriesStat, SummaryData, SUM_ROWS, SWAP_MIN_PAIRS,
  SWAP_MIN_SERIES, SWAP_ROW_CAP, TRUST_R2} from './sar-matrix-summary-data';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, count} from './sar-matrix-ui-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {STRIP_MIN_FILL, ROLE_SPREAD_TIE, CORES_NOT_COMPARABLE, PANE_EFFECTS, PANE_SERIES, orient, Answer,
  shortSmiles, swapAnchor, roleValueName, roleRanks} from './sar-matrix-summary-common';

/** What this segment reads from the panel: the shared rendering helpers and the state segments share. */
export type SummaryEffectsKit = Pick<SummaryPanel, 'anchored' | 'answerArt' | 'badge' | 'cellLocation' |
  'depiction' | 'effectsTab' | 'expandTopRGroup' | 'flushPaints' | 'formatDelta' | 'formatEffect' | 'hint' |
  'host' | 'leads' | 'openTip' | 'pendingPaints' | 'prose' | 'rankableCores' | 'reason' | 'revealTrust' |
  'roleFitRefusal' | 'scroller' | 'showPane' | 'summaryRow' | 'topCores' | 'trustDot'>;

export class SummaryEffects {
  constructor(private readonly kit: SummaryEffectsKit) {}

  /**
   * One tab per component column, then the per-series evidence.
   *
   * Stacked, this segment is as tall as the decomposition is wide — five component columns is five
   * cards before the two that were already here — and a landing screen that has to be scrolled to be
   * read is not one. The ordering sentence stays above the strip because it is the only finding here
   * that is about every tab at once, and each tab carries its own spread so the comparison between
   * components does not cost a click each.
   */
  build(data: SummaryData): HTMLElement {
    // The R-group card is dropped wherever the fit already ranks the axis role, for the reason the core
    // card is: two rankings of one column, one a count of within-series wins and one an adjusted effect,
    // read as a disagreement rather than as two kinds of evidence.
    const series: HTMLElement[] = this.roleAnswer(data, data.axisRole) === null ?
      [this.kit.anchored(this.buildRGroupCard(data), 'rgroup')] : [];
    // Dropped only where the core column is already ranked by the fit, adjusted for the components it
    // was paired with: two rankings of one column, one confounded and one adjusted, in different units,
    // is a worse screen than either alone. Where the fit declines to rank it, this is the only one left.
    let coreNote: HTMLElement | null = null;
    if (this.roleAnswer(data, data.coreRole) === null) {
      if (data.coresAreSeries)
        series.push(this.kit.anchored(this.buildCoreCard(data, this.kit.topCores(data)), 'cores'));
      else {
        // Where cores are not comparable the card has no ranking it could ever carry, and a card slot
        // spent on one sentence is a slot. It keeps the anchor so the tile that points here still lands.
        coreNote = this.kit.anchored(ui.divText(CORES_NOT_COMPARABLE,
          'chem-sar-cp-hint chem-sar-sum-sub-lead'), 'cores');
      }
    }
    // One card per component, matching the band that links here. A single card carried whichever
    // component the matrix columns enumerated, so clicking the Linker row landed on a card about
    // Warhead — the band's own answer contradicted on arrival.
    const roles = data.roleFit === null ? [] : data.roleFit.roles.map((role) => role.name);
    if (roles.length === 0)
      series.push(this.kit.anchored(this.buildSwapCard(data), 'swaps'));
    else {
      for (const name of roles)
        series.push(this.kit.anchored(this.buildSwapCard(data, name), swapAnchor(name)));
    }

    const roleCards = this.buildRoleCards(data);
    if (roleCards.length === 0) {
      return coreNote === null ? this.effectsScroll(series) :
        ui.divV([coreNote, this.effectsScroll(series)], 'chem-sar-sum-effects');
    }

    const tabs = [...roleCards.map((card) => ({label: card.label, note: card.note, cards: [card.el]})),
      {label: PANE_SERIES, note: 'counted, not fitted' as string | null, cards: series}];
    const open = Math.min(Math.max(this.kit.effectsTab, 0), tabs.length - 1);
    const bar = ui.divH(tabs.map((tab, i) => {
      const parts = [ui.divText(tab.label, 'chem-sar-sum-sub-name')];
      if (tab.note !== null)
        parts.push(ui.divText(tab.note, 'chem-sar-sum-sub-note'));
      const el = ui.divH(parts, 'chem-sar-sum-sub');
      el.classList.toggle('chem-sar-sum-sub-on', i === open);
      el.onclick = () => {
        this.kit.effectsTab = i;
        this.kit.showPane(PANE_EFFECTS);
      };
      return el;
    }), 'chem-sar-sum-sub-bar');

    const lead = this.roleOrdering(data);
    return ui.divV([...(lead === null ? [] : [lead]), bar, this.effectsScroll(tabs[open].cards)],
      'chem-sar-sum-effects');
  }

  private effectsScroll(cards: HTMLElement[]): HTMLElement {
    const grid = ui.div(cards, 'chem-sar-sum-grid');
    // A lone card has no neighbour to sit beside, so left at one track's width it strands the rest of
    // the pane. The width is not decoration here: these row names are SMILES, and at one track three
    // analogs of one scaffold truncate to the same prefix and read as the same structure.
    grid.classList.toggle('chem-sar-sum-grid-one', cards.length === 1);
    const scroll = ui.div([grid], 'chem-sar-sum-scroll');
    this.kit.scroller = scroll;
    return scroll;
  }

  /**
   * MK-D: a signed magnitude against the screen's own model error.
   *
   * `scale` is the largest magnitude in the set this bar is compared against — the card's own rows,
   * except where one fit produced the rows of several cards and they are genuinely on one scale. A bar
   * ending inside the grey band is one the analysis cannot resolve.
   */
  private effectBar(value: number, scale: number, modelError: number | null): HTMLElement {
    const track = ui.div([], 'chem-sar-sum-effbar');
    const half = scale > 0 ? scale : 1;
    // No band at all rather than a zero-width one: a zero-width band reads as "everything here is
    // resolved", which is the opposite of what a missing model error means.
    if (modelError !== null && modelError > 0) {
      const band = ui.div([], 'chem-sar-sum-effband');
      const half2 = Math.min(50, (modelError / half) * 50);
      band.style.left = `${50 - half2}%`;
      band.style.width = `${2 * half2}%`;
      track.appendChild(band);
    }
    track.appendChild(ui.div([], 'chem-sar-sum-effzero'));
    const w = Math.min(50, (Math.abs(value) / half) * 50);
    const bar = ui.div([], `chem-sar-sum-effval chem-sar-sum-eff-${value >= 0 ? 'pos' : 'neg'}`);
    bar.style.left = `${value >= 0 ? 50 : 50 - w}%`;
    bar.style.width = `${w}%`;
    track.appendChild(bar);
    ui.tooltip.bind(track, () => `${this.kit.formatEffect(value)}, against the widest bar it is compared ` +
      `with (${this.kit.formatEffect(half)}).` + (modelError === null ? ' No series has a cross-validated fit, ' +
      'so there is no error band to read it against.' :
      ` The grey band is ± ${this.kit.host.formatActivity(modelError)}, this analysis\' own typical ` +
      'prediction error: a bar ending inside it is a difference this analysis cannot resolve.'));
    return track;
  }

  /** MK-C: one slot per series that tried this group — solid where it won and that series' fit holds,
   *  hatched where it won and the fit does not, outlined where it was tried and did not win. */
  private winTally(row: RGroupRow, placed: string): HTMLElement {
    const shown = Math.min(row.tried, 10);
    const slots: HTMLElement[] = [];
    for (let i = 0; i < shown; i++) {
      const state = i < row.k ? ' chem-sar-sum-slot-on' :
        i < row.m ? ' chem-sar-sum-slot-hatch' : '';
      const slot = ui.div([], `chem-sar-sum-wslot${state}`);
      if (i >= row.k && i < row.m) {
        slot.onclick = (e: MouseEvent) => {
          e.stopPropagation();
          this.kit.revealTrust();
        };
      }
      slots.push(slot);
    }
    const parts: HTMLElement[] = [ui.divH(slots, 'chem-sar-sum-tally')];
    if (row.tried > shown)
      parts.push(ui.divText(`+${row.tried - shown}`, 'chem-sar-sum-tally-more'));
    const box = ui.divH(parts, 'chem-sar-sum-tally-box');
    ui.tooltip.bind(box, () => `${placed} in ${row.k} of the ${row.tried} series that tried it` +
      (row.m > row.k ? `, and in ${row.m - row.k} more whose additive fit does not hold — a group that ` +
      `is ${placed} in additive series and not in non-additive ones is core-dependent. Click a hatched ` +
      'slot to see those series.' : '. Every fit behind those placings holds.'));
    return box;
  }

  /** One role's tile, read off the global fit — the same answer for the axis role and the core role,
   *  since in fragment-columns mode one fit ranks both. Null where that role's card is not ranking. */
  private roleAnswer(data: SummaryData, name: string | null): Answer | null {
    const fit = data.roleFit;
    if (fit === null || this.kit.roleFitRefusal(data) !== null)
      return null;
    const role = fit.roles.find((other) => other.name === name);
    if (role === undefined || !roleRanks(role))
      return null;
    const [top, second] = role.levels;
    // Signed by direction: `levels` is best-first, so on a lower-is-better column the leader's
    // coefficient is the smaller number.
    const dir = this.kit.host.higherIsBetter ? 1 : -1;
    if (!this.kit.leads(dir * top.coef, dir * second.coef, fit.cvRmse)) {
      // A bimodal component has an answer no leaderboard can state: which group, not which value. The
      // top two being inseparable is exactly the case where the split is the finding.
      const split = role.split;
      if (split !== null) {
        return {answer: `two groups: ${count(split.hiCount)} at ${this.kit.formatEffect(split.hiMean)}, ` +
          `${count(split.loCount)} at ${this.kit.formatEffect(split.loMean)}`,
        negative: `against the library mean of ${this.kit.host.formatActivity(fit.mean)}; no single ` +
          `${role.name} leads, but the groups are ${split.gap.toFixed(2)} apart — which group a ` +
          'compound carries matters more than the choice inside it'};
      }
      return {answer: `No single ${role.name} leads`,
        negative: `top two ${this.kit.formatEffect(dir * (top.coef - second.coef))} apart, within the ` +
          `± ${this.kit.host.formatActivity(fit.cvRmse!)} this fit resolves`};
    }
    return {answer: `${this.kit.formatEffect(top.coef)} · over ${count(top.n)} compounds`,
      negative: `against the library mean of ${this.kit.host.formatActivity(fit.mean)}, adjusted for the ` +
        'other components — observational, not a potency', art: this.kit.answerArt(top.value)};
  }

  private card(title: string, subtitle: string, body: HTMLElement[], footer?: HTMLElement): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-sum-card-title'),
      ui.divText(subtitle, 'chem-sar-cp-hint'), ...body];
    if (footer !== undefined)
      parts.push(footer);
    return ui.divV(parts, 'chem-sar-sum-card');
  }

  private buildSwapCard(data: SummaryData, role?: string): HTMLElement {
    const host = this.kit.host;
    const pools = role === undefined ? data.swaps : data.swapsByRole.get(role) ?? [];
    const title = role === undefined ? 'What buys potency inside these series' :
      `Swapping ${role} — what it bought`;
    const subtitle = 'Differences between two compounds that were both made and both measured, alike in ' +
      'everything but one R-group. Nothing here is fitted or predicted: if the pair was never made, it ' +
      'is not on this card.';
    const unit = role === undefined ? 'series' : 'contexts';
    const rows = pools.map((pool) => {
      const {forward, from, to, worst, widest: bestOf, mean, up} = orient(pool);
      const allUp = up === pool.n;
      // Null for a pool grouped off the component columns: the pair was found by matching component
      // values, so there is no one cell of one matrix it can be opened at.
      const best = pool.best;
      const ci = best === null ? -1 : forward ? best.ciTo : best.ciFrom;
      const position = best === null ? '' : best.matrix.columns[ci].position;
      const show = (value: number): string => this.kit.formatDelta(value);
      // Badges, not descriptor tail: the descriptor clips to one line, and these three are exactly the
      // qualifiers that stop the row being over-read — appended last they are the first characters lost.
      const badges: HTMLElement[] = [];
      if (allUp) {
        badges.push(this.kit.badge('all up', 'Every measured pair moved the same way. At 3 pairs that is ' +
          'one-in-four under a coin-flip null.'));
      }
      if (pool.sampled) {
        badges.push(this.kit.badge('sampled row', 'One contributing row carried more than ' +
          `${SWAP_ROW_CAP} measured cells; its most and least potent halves were kept and mid-range ` +
          'pairs dropped.', true));
      }
      // Only on a log scale: the pooled delta is a log ratio there and the model error is in the same
      // units, while on a raw scale the two are different quantities.
      if (host.activityIsLog && data.modelError !== null && Math.abs(worst) < data.modelError) {
        badges.push(this.kit.badge('under resolution', 'The worst case is smaller than this analysis\' own ' +
          'typical prediction error, so the whole row may be noise.', true));
      }
      const pair = ui.divH([
        this.kit.depiction(from, BENEFIT_MOL_W, BENEFIT_MOL_H),
        ui.divText('→', 'chem-sar-sum-arrow'),
        this.kit.depiction(to, BENEFIT_MOL_W, BENEFIT_MOL_H),
      ], 'chem-sar-sum-pair');
      return this.kit.summaryRow({
        depiction: pair,
        name: best === null ? `${role} swap` : `${best.matrix.label} · ${position}`,
        badges,
        desc: `mean ${show(mean)} (${show(worst)} → ${show(bestOf)}) · ${pool.n} pairs · ` +
          `${pool.roots.size} ${unit}`,
        value: show(worst),
        caption: 'worst case',
        valueTip: `This swap was worth at least ${show(worst)} in every one of the ${pool.n} measured ` +
          `pairs, across ${pool.roots.size} ${unit}. Each pair is two measured compounds alike in ` +
          'everything else, so what they do not share cancels exactly and no fit is involved. At 3 ' +
          'pairs "all up" is one-in-four under a coin-flip null.' +
          (pool.sampled ? ` One group carried more than ${SWAP_ROW_CAP} measured compounds; its most ` +
            'and least potent halves were kept and mid-range pairs dropped.' : ''),
        // A pool grouped off the component columns has no one cell to open, so the row selects the
        // compounds that carry the better side instead of going nowhere.
        tip: best === null ? `Select every compound whose ${role} is the one on the right` :
          this.kit.openTip(best.matrix, best.ri, ci),
        onClick: best === null ? () => host.selectRoleValue(role!, to) :
          () => host.revealCell(best.matrix, best.ri, ci, position),
      });
    });
    if (rows.length > 0)
      return this.card(title, subtitle, rows);
    return this.card(title, subtitle, this.kit.reason(data.swapCandidates === 0 && role === undefined ?
      'No row of any series carries two measured substituents, so there is no swap to measure.' :
      `No ${role === undefined ? '' : `${role} `}swap clears ${SWAP_MIN_PAIRS} measured pairs in ` +
      `${SWAP_MIN_SERIES} ${unit}.`));
  }

  private buildRGroupCard(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const title = data.axisRole === null ? 'R-groups — within-series ranking' :
      `${data.axisRole} — within-series ranking`;
    const subtitle = 'Ranked by how often the fitted model places this group first at its position. The ' +
      'number is its margin over that series\' own most-common substituent — a within-series ' +
      'comparison, so two series\' numbers are only loosely comparable.';
    // An additive effect in raw assay units cannot be re-expressed as a fold, which is the one thing
    // that would make it mean something; the swap card refuses the same arithmetic.
    if (!host.activityIsLog) {
      return this.card(title, subtitle, this.kit.reason('The activity is on a raw scale, where a fitted ' +
        'effect is a difference in assay units and says nothing about how large the change is. Set ' +
        'Scaling to lg or −lg, or declare the column higher-is-better if it is already a pIC50.'));
    }
    const slots = data.series.filter((stat) => stat.stripCols !== null);
    // One scale for every bar on this card, so the bars compare down the card — and, deliberately, to
    // nothing on any other card.
    const shown = [...data.rgroups.slice(0, SUM_ROWS), ...data.rgroupsThin, ...data.rgroupLosers];
    const scale = shown.reduce((m, r) => Math.max(m, Math.abs(r.magnitude ?? 0)), 0);
    // With one row the bar is full-length and says nothing.
    const bars = shown.filter((r) => r.magnitude !== null).length > 1 ? scale : 0;
    // One strip scale for the whole card, like the bars: normalised per row, every row's own best
    // difference draws at full height, so a +1.8 row and a +0.05 row terminate alike.
    let stripScale = 0;
    for (const stat of slots) {
      for (const r of shown) {
        const delta = stat.stripCols!.get(r.subst)?.refDelta;
        if (delta != null)
          stripScale = Math.max(stripScale, Math.abs(delta));
      }
    }
    const rows = data.rgroups.slice(0, SUM_ROWS).map((row, i) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {primary: i === 0}));
    const thin = data.rgroupsThin.map((row) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {thin: true}));
    if (rows.length === 0 && thin.length === 0) {
      // Convergence is named because it is the one cause the reader cannot see anywhere else: a fit
      // that stopped short is silently excluded from every pool, and R² says nothing about it.
      return this.card(title, subtitle, this.kit.reason(`No substituent comes first at its position in ` +
        `${SWAP_MIN_SERIES} series whose additive fit holds.` + (data.nonConverged === 0 ? '' :
        ` The additive fit of ${count(data.nonConverged)} of the ${count(this.kit.host.matrices.length)} ` +
        'series did not reach its tolerance, and those are left out of this comparison.')));
    }
    const body: HTMLElement[] = [];
    // Two numbers for one substituent are a tab apart, and nothing else says which question each
    // answers. Only where that second number exists: where the fit declined to rank this column there
    // is no other card to reconcile with.
    if (this.roleAnswer(data, data.axisRole) !== null) {
      body.push(ui.divText('A value can place first in most series and still sit at the library mean: ' +
        'this card counts wins within a series against that series\' own reference, while the ' +
        `${data.axisRole} tab reports one offset pooled over the table.`,
      'chem-sar-sum-prose'));
    }
    if (slots.length > 0 && data.nonConverged > 0) {
      body.push(ui.divText(`${count(data.nonConverged)} further series carry no square and count in ` +
        'no row: their additive fit did not reach its tolerance, so their effects cannot be ordered ' +
        'against another series\'.', 'chem-sar-cp-hint'));
    }
    body.push(...rows);
    if (thin.length > 0) {
      body.push(ui.divText(`Won in only two series — a median over two lineages is the mean of two ` +
        'numbers:', 'chem-sar-cp-hint'));
      body.push(...thin);
    }
    if (data.rgroupLosers.length > 0) {
      const head = data.rgroupLosers[0];
      body.push(ui.divText(`Stop making: ${shortSmiles(head.subst)} came last at its position in ` +
        `${head.k} of the ${head.tried} series it was tried in` +
        (head.magnitude === null ? '.' : `, costing ${this.kit.formatEffect(head.magnitude)} against each ` +
        'series\' own reference.'), 'chem-sar-cp-hint'));
      body.push(...data.rgroupLosers.map((row) =>
        this.rgroupBlock(data, slots, row, bars, stripScale, {loser: true})));
    }
    return this.card(title, subtitle, body, this.sizeVerdict(data));
  }

  /** What one strip slot is. A strip is drawn for any fragment-columns run, but only where a matrix IS
   *  one core is a slot a core — with a series column, or with nothing on the row axis, one slot holds
   *  several cores. */
  private slotNoun(data: SummaryData): string {
    return data.coresAreSeries ? 'core' : 'series';
  }

  /** One leaderboard entry: the row, its per-slot strip where one can be drawn, and an expand that
   *  says what the row does and does not claim. */
  private rgroupBlock(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    stripScale: number, opts: {thin?: boolean, loser?: boolean, primary?: boolean}): HTMLElement {
    const parts: HTMLElement[] = [this.rgroupRow(data, slots, row, scale, opts)];
    if (slots.length > 0)
      parts.push(this.buildStrip(data, slots, row, stripScale));
    const body = ui.divV([]);
    body.style.display = 'none';
    const toggle = ui.divText('What this row claims →', 'chem-sar-cp-hint');
    toggle.classList.add('chem-sar-sum-click');
    let filled = false;
    const open = (): void => {
      if (!filled) {
        filled = true;
        for (const part of this.rgroupExpandParts(data, row))
          body.appendChild(part);
        this.kit.flushPaints();
      }
      body.style.display = '';
    };
    toggle.onclick = () => {
      if (body.style.display === 'none')
        open();
      else
        body.style.display = 'none';
    };
    // The answer tile above asked for its evidence, and the pane it asked for is only built now.
    if (opts.primary && this.kit.expandTopRGroup) {
      this.kit.expandTopRGroup = false;
      open();
    }
    parts.push(toggle, body);
    return ui.divV(parts);
  }

  private rgroupExpandParts(data: SummaryData, row: RGroupRow): HTMLElement[] {
    const slot = this.slotNoun(data);
    // A reading instruction, not a restatement of the marks above it.
    const read = ui.divText(`How to read this: first on most of the ${slot}s it was tried on and never ` +
      `worse than the reference means the group travels — carry it forward. First on one ${slot} while ` +
      `costing potency on the rest means hold it to that ${slot} and vary elsewhere.`,
    'chem-sar-sum-prose');
    const art = ui.divH([this.kit.depiction(row.subst, CARD_CORE_W * 2, CARD_CORE_H * 2), read],
      'chem-sar-sum-detail');
    return [art];
  }

  /**
   * MK-A. One square per slot, in a fixed order that is never re-sorted per row — the row is only
   * readable against its neighbours if the same slot means the same series everywhere.
   *
   * Sign is carried by the direction the fill grows from the mid-rule, and only redundantly by hue:
   * a tint alone renders +0.9 and −0.9 as the same grey. The mid-rule is the legend.
   *
   * `widest` is the card's largest difference, not this row's, so fill heights compare down the card
   * the way the slot order lets the squares compare across it.
   */
  private buildStrip(data: SummaryData, slots: SeriesStat[], row: RGroupRow,
    widest: number): HTMLElement {
    const host = this.kit.host;
    const slot = this.slotNoun(data);
    const squares = slots.map((stat) => {
      const col = stat.stripCols!.get(row.subst);
      if (col === undefined) {
        const none = ui.divText('–', 'chem-sar-sum-sq chem-sar-sum-sq-none');
        ui.tooltip.bind(none, () => `${stat.matrix.label} — never tried on this ${slot}`);
        return none;
      }
      const square = ui.div([], 'chem-sar-sum-sq');
      const delta = col.refDelta;
      if (delta === null) {
        square.classList.add('chem-sar-sum-sq-hatch');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — tried here over ${col.n} measured cells, ` +
          'but this series records no reference substituent to compare against');
      } else if (delta === 0 && stat.matrix.refValues[stat.matrix.columns[col.ci].position] === row.subst) {
        square.classList.add('chem-sar-sum-sq-ref');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — this group IS this series' reference ` +
          `substituent, over ${col.n} measured cells, so there is nothing to compare it against here`);
      } else {
        square.classList.add('chem-sar-sum-sq-val');
        // No difference draws no fill: the floor exists so a small value still reads as a direction,
        // and applied to an exact zero it would draw a win that is not there.
        if (delta !== 0) {
          const fill = ui.div([], `chem-sar-sum-fill chem-sar-sum-${delta > 0 ? 'up' : 'down'}`);
          fill.style.height = `${Math.max(STRIP_MIN_FILL,
            Math.round((widest > 0 ? Math.abs(delta) / widest : 0) * 9))}px`;
          square.appendChild(fill);
        }
        ui.tooltip.bind(square, () => `${stat.matrix.label} · ${this.kit.formatEffect(delta)} against this ` +
          `series' own reference, over ${col.n} measured cells` + (delta === 0 ?
          ' — no difference at all, so nothing grows from the mid-rule.' :
          ` — ${delta > 0 ? 'better than' : 'worse than'} it, so the fill grows ` +
          `${delta > 0 ? 'up' : 'down'} from the mid-rule. Heights compare down this card, against its ` +
          `widest difference (${this.kit.formatEffect(widest)}).`));
      }
      square.classList.add('chem-sar-sum-click');
      square.onclick = () => {
        const ri = bestMeasuredRow(stat.matrix, col.ci, this.kit.host.higherIsBetter ? 1 : -1);
        if (ri >= 0)
          host.revealCell(stat.matrix, ri, col.ci, stat.matrix.columns[col.ci].position);
      };
      return square;
    });
    return ui.divH(squares, 'chem-sar-sum-strip');
  }

  /** One verdict for the whole card, not a number per row: five heavy-atom counts do not add
   *  themselves up in a reader's head, and it is the sum that ends a series at MW 620. */
  private sizeVerdict(data: SummaryData): HTMLElement | undefined {
    const rows = data.rgroups;
    if (rows.length === 0)
      return undefined;
    const line = ui.divText('', 'chem-sar-cp-hint');
    this.kit.pendingPaints.push(() => {
      const rdkit = getRdKitModule();
      let groups = 0;
      let references = 0;
      let n = 0;
      for (const row of rows) {
        const {matrix, ci} = row.win;
        const reference = matrix.refValues[matrix.columns[ci].position];
        const atoms = cachedAtomCount(row.subst, rdkit);
        const refAtoms = reference ? cachedAtomCount(reference, rdkit) : 0;
        if (atoms > 0 && refAtoms > 0) {
          groups += atoms;
          references += refAtoms;
          n++;
        }
      }
      line.innerText = n === 0 ? '' : `Your top ${n} groups average ${Math.round(groups / n)} heavy ` +
        `atoms; the series references they beat average ${Math.round(references / n)}.`;
    });
    ui.tooltip.bind(line, 'Heavy atoms of the fragment, the attachment point included — not ' +
      'molecular weight, and not lipophilicity, which is what actually kills the compound. A smoke ' +
      'alarm, not a property model.');
    return line;
  }

  private rgroupRow(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    opts: {thin?: boolean, loser?: boolean}): HTMLElement {
    const host = this.kit.host;
    const {win} = row;
    const position = win.matrix.columns[win.ci].position;
    const reference = win.matrix.refValues[position] ?? '';
    const atoms = this.kit.badge('', 'Heavy atoms of the fragment, the attachment point included. The ' +
      'model scores potency alone and carries no efficiency term, so an unlabelled leaderboard is a ' +
      'heavy-atom leaderboard wearing a potency label.');
    this.kit.pendingPaints.push(() => {
      const n = cachedAtomCount(row.subst, getRdKitModule());
      atoms.innerText = n ? `${n} atoms` : '—';
    });
    // Attempts, not wins: "best in 5 of 6" is heard as "tried in 6", and a denominator of wins makes
    // every row on the card read as near-perfect.
    const placed = opts.loser ? 'last' : 'best';
    const badges = [this.winTally(row, placed), atoms];
    if (opts.thin)
      badges.push(this.kit.badge('2 series', 'Two lineages only — the number below is the mean of two.', true));
    // Only a single contributing series may have its reference named: each series is measured against
    // its own most-common substituent, and those genuinely differ, so one SMILES cannot stand for the
    // comparator of a median pooled over several.
    const comparator = row.magnitudeSeries > 1 ? 'each series\' own reference' : shortSmiles(reference);
    return this.kit.summaryRow({
      depiction: this.kit.depiction(row.subst, CARD_CORE_W, CARD_CORE_H),
      name: `${win.matrix.label} · ${position}`,
      badges,
      // The strip below carries the per-slot picture where one can be drawn, and a range line beside it
      // would say the same thing twice and worse.
      desc: row.magnitude === null ? 'no reference substituent to compare against' :
        slots.length > 0 ? `vs ${comparator} over ${row.magnitudeSeries} series` :
          `vs ${comparator} · ${this.kit.formatEffect(row.lo)} → ${this.kit.formatEffect(row.hi)} ` +
          `over ${row.magnitudeSeries} series`,
      mark: row.magnitude === null || scale <= 0 ? undefined :
        this.effectBar(row.magnitude, scale, data.modelError),
      value: row.magnitude === null ? '—' : this.kit.formatEffect(row.magnitude),
      caption: 'vs reference',
      valueTip: row.magnitude === null ?
        'No series contributing this placing also carries its reference substituent as a column, so ' +
        'there is no within-series comparison to quote. The row still ranks on how often it placed.' :
        'Median difference against each series\' own most-common substituent' +
        (row.magnitudeSeries > 1 ? ', which differs from series to series' : ` (${reference})`) +
        '. A within-series difference, so it is free of the fact that two series centre their effects ' +
        `over different substituent menus — but only ${row.magnitudeSeries} of the ${row.k} ` +
        'contributing series carry one, and two series\' numbers are still only loosely comparable.',
      tip: this.kit.openTip(win.matrix, win.ri, win.ci),
      onClick: () => host.revealCell(win.matrix, win.ri, win.ci, position),
    });
  }

  /**
   * One card per role column, ordered by that role's noise-corrected spread, or one card naming why
   * nothing is ranked.
   *
   * Absent rather than empty outside fragment-columns mode: a substituent label discovered by
   * fragmentation is local to its own series, so there is no scale on which one pooled offset could be
   * read and the question does not apply at all.
   */
  private buildRoleCards(data: SummaryData): {label: string, note: string | null, el: HTMLElement}[] {
    if (data.axisRole === null)
      return [];
    const refusal = this.kit.roleFitRefusal(data);
    if (refusal !== null) {
      return [{label: 'Components', note: null,
        el: this.kit.anchored(this.card('Component contributions', refusal, this.roleDropped(data)),
          'role-0')}];
    }
    const roles = data.roleFit!.roles;
    // One scale across every role card, against the tab's usual rule: each of these offsets comes from
    // one fit against one library average, so they are on one scale and a shared bar is the honest
    // rendering. Per-matrix effects elsewhere are centred inside their own matrix and are not. Over
    // the rows that will be drawn only — a coefficient on no card would shorten every bar on every
    // other one for a reason the reader cannot see.
    let scale = 0;
    for (const role of roles.filter(roleRanks)) {
      for (const {level} of this.roleRows(role))
        scale = Math.max(scale, Math.abs(level.coef));
    }
    return roles.map((role, i) => ({
      label: role.name,
      note: roleRanks(role) ? role.spread.toFixed(2) : null,
      el: this.kit.anchored(this.buildRoleCard(data, role, scale), `role-${i}`),
    }));
  }

  /** The one finding that is about every component at once, so it sits above the tabs rather than on
   *  any one of them. The order of the roles is stable; the size of the gap between two of them is not,
   *  so it is named as a rank and never quoted as a ratio. Omitted where a role is refusing to rank its
   *  own values, since that refusal is what invalidates the spread ordering. */
  private roleOrdering(data: SummaryData): HTMLElement | null {
    const roles = data.roleFit?.roles;
    if (roles === undefined || this.kit.roleFitRefusal(data) !== null || !roles.every(roleRanks))
      return null;
    // Consecutive roles closer than the tie threshold are joined by an approximation sign rather than
    // an inequality: the order is stable, a gap that small is not, and one chain says so however many
    // component columns there are.
    const chain = roles.map((role, i) => i === 0 ? role.name :
      `${roles[i - 1].spread - role.spread < ROLE_SPREAD_TIE ? '≈' : '>'} ${role.name}`).join(' ');
    const el = this.kit.hint(`Changing ${roles[0].name} moves ${this.kit.host.activityColumnName} most: ${chain}.`,
      'Each component\'s range between its best and worst value, corrected for the fact that a column ' +
      'with more values gets a wider range by chance alone. That correction is what makes columns with ' +
      'different numbers of values comparable. "≈" marks a gap too small to order.');
    el.classList.add('chem-sar-sum-sub-lead');
    return el;
  }

  /** Stated wherever the fit is reported: the connected-block prune drops measured compounds by design,
   *  so the compound count the card prints is not the table. */
  private roleDropped(data: SummaryData): HTMLElement[] {
    const dropped = data.roleFit?.dropped ?? 0;
    if (dropped === 0)
      return [];
    return [ui.divText(`${count(dropped)} compounds share no component value with the rest of the table; ` +
      'only the largest connected block is fitted, because offsets from separate blocks rest on ' +
      'unrelated baselines.', 'chem-sar-sum-prose')];
  }

  /** The rows one role card shows: its best three, then its worst two. Two losers rather than one —
   *  where a role falls into two groups the single lowest value is routinely its least-supported
   *  member, and one row would then stand for a whole cluster by its smallest. */
  private roleRows(role: RoleSummary): {level: RoleLevel, faint: boolean}[] {
    const top = role.levels.slice(0, SUM_ROWS);
    return [...top.map((level) => ({level, faint: false})),
      ...role.levels.slice(Math.max(top.length, role.levels.length - LOSER_ROWS))
        .map((level) => ({level, faint: true}))];
  }

  /**
   * Which component the matrix columns actually enumerate, and the one control that changes it.
   *
   * A card ranking R4 sits beside a grid whose columns are R1, which reads as the two disagreeing. They
   * are answering different questions: the card is the fit over the whole table, the grid enumerates one
   * component at a time.
   */
  private axisSwitchNote(name: string, data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const note = this.kit.hint('Ranked from the fit over the whole table. The SAR Matrix columns enumerate ' +
      `${data.axisRole}, so these values are not the ones across the top there.`,
    `Putting ${name} across the columns means rebuilding: every matrix is reassembled and every fit ` +
      'recomputed, which takes as long as the first run did.');
    const pill = ui.divText(`Put ${name} across the columns — rebuilds`,
      'chem-sar-chip-badge chem-sar-sum-role');
    // What will and will not move, on the control: this ranking is the same fit either way, so a reader
    // who clicks expecting these rows to change waits out a rebuild for a screen that looks identical.
    ui.tooltip.bind(pill, () => `Reassembles every matrix with ${name} across the columns, and pools ` +
      `its measured pairs instead of ${data.axisRole}'s. The offsets on this card do not change — they ` +
      'come from one fit over the whole table, whichever component the columns enumerate.');
    pill.onclick = () => {
      if (!host.computing)
        host.setColumnAxis(name);
    };
    return ui.divH([note, pill], 'chem-sar-sum-chips');
  }

  private buildRoleCard(data: SummaryData, role: RoleSummary, scale: number): HTMLElement {
    const host = this.kit.host;
    const fit = data.roleFit!;
    const title = `${role.name} — offsets from the additive fit`;
    const subtitle = `Additive fit over ${fit.roles.length} component columns. Each offset is against ` +
      `the library mean of ${host.formatActivity(fit.mean)}, adjusted for the other components. ` +
      `R² ${fit.cvR2!.toFixed(2)} ± ${host.formatActivity(fit.cvRmse!)} predicting unseen compounds, ` +
      `n = ${count(fit.compounds)}.`;
    const footer = this.roleFooter();
    // On every card, because each is now the only one its reader can see: a subtitle counting fewer
    // compounds than the tab does, on a card that never explains the difference, is the one omission
    // this note exists to prevent.
    const dropped = this.roleDropped(data);
    if (role.levels.length < 2) {
      return this.card(title, subtitle, [...this.kit.reason(`${role.levels.length === 1 ? 'Only one' : 'No'} ` +
        `${role.name} value carries ${MIN_SUPPORT} or more compounds, so there is nothing to compare.`),
      ...dropped], footer);
    }
    // Prediction is nearly unique even where the split of credit between the roles is not, so this is
    // the only check that this role's share of it means anything.
    if (role.repeat !== null && role.repeat < TRUST_R2) {
      return this.card(title, subtitle, [...this.kit.reason(`${role.name} offsets do not repeat: fitted on ` +
        `each half of the table separately they correlate at r ${role.repeat.toFixed(2)}, below ` +
        `${TRUST_R2}. The other components in this table predict which ${role.name} a compound carries ` +
        'closely enough that the fit cannot separate their contributions on a subset, so these values ' +
        'are not ranked.'), ...dropped], footer);
    }

    const chips: HTMLElement[] = [
      this.kit.badge(`spread ${role.spread.toFixed(2)}`, 'Count-weighted sd of this component\'s offsets, ' +
        'with estimation noise subtracted so columns with different numbers of values compare. ' +
        'Not a range.'),
      this.kit.badge(`${role.levels.length} of ${role.fitted} values ranked`,
        `${role.thin} values seen fewer than ${MIN_SUPPORT} times stay in the fit but carry no readable ` +
        'offset of their own.'),
    ];
    if (role.repeat !== null) {
      chips.push(this.kit.badge(`repeats at r ${role.repeat.toFixed(2)}`,
        'Fitted on each half of the table separately, these offsets correlate this closely. Below ' +
        TRUST_R2 + ' they would not be ranked.'));
    }
    const body: HTMLElement[] = [ui.divH(chips, 'chem-sar-sum-chips')];
    if (role.name !== data.axisRole && host.roleColumns.includes(role.name))
      body.push(this.axisSwitchNote(role.name, data));
    const split = role.split;
    if (split !== null) {
      body.push(ui.divText(`Bimodal: ${split.hiCount} values at ${this.kit.formatEffect(split.hiMean)} ` +
        `(n = ${count(split.hiN)}) and ${split.loCount} at ${this.kit.formatEffect(split.loMean)} ` +
        `(n = ${count(split.loN)}), a gap of ${split.gap.toFixed(2)} against the ` +
        `± ${host.formatActivity(fit.cvRmse!)} this fit resolves.`, 'chem-sar-sum-prose'));
    }

    const best = data.roleBest.get(role.name)!;
    const others = fit.roles.filter((other) => other !== role).map((other) => other.name).join(', ');
    for (const {level, faint} of this.roleRows(role)) {
      const ref = best.get(level.value)!;
      // The axis is the one role pinned in the Vary filter, so its row can land on the position too.
      const position = role.name === data.axisRole ? ref.matrix.columns[ref.ci].position : undefined;
      body.push(this.kit.summaryRow({
        depiction: level.value === '' ? null : this.kit.depiction(level.value, CARD_CORE_W, CARD_CORE_H),
        name: roleValueName(level.value),
        badges: [this.roleCountBadge(role.name, level)],
        desc: `best measured ${host.formatActivity(ref.value)} · ${ref.matrix.label}`,
        mark: this.effectBar(level.coef, scale, fit.cvRmse),
        value: this.kit.formatEffect(level.coef),
        caption: 'offset',
        valueTip: `${this.kit.formatEffect(level.coef)} against the library mean of ` +
          `${host.formatActivity(fit.mean)}, over ${count(level.n)} compounds, adjusted for ${others}. ` +
          'An offset, not a potency.',
        faint,
        tip: role.name === data.coreRole ? `Open ${ref.matrix.label} in the SAR Matrix` :
          this.kit.openTip(ref.matrix, ref.ri, ref.ci),
        onClick: () => role.name === data.coreRole ? host.revealMatrix(ref.matrix) :
          host.revealCell(ref.matrix, ref.ri, ref.ci, position),
      }));
    }
    body.push(...dropped);
    return this.card(title, subtitle, body, footer);
  }

  /** The compound count, which also selects those compounds. An offset is a statement about a subset
   *  of the table, and the only way to check one is to look at the subset it is about. */
  private roleCountBadge(role: string, level: RoleLevel): HTMLElement {
    const badge = this.kit.badge(`${count(level.n)} cpd`, 'Measured compounds carrying this value — click ' +
      'to select them in the table. No ± is printed: one from this count alone would ignore the ' +
      'degrees of freedom the other components cost.');
    badge.classList.add('chem-sar-sum-click');
    badge.onclick = (e: MouseEvent) => {
      e.stopPropagation();
      this.kit.host.selectRoleValue(role, level.value);
    };
    return badge;
  }

  /** The negative clause and the reconciliation, on every role card: each is now the only card its
   *  reader can see, so neither can be carried by a neighbour. */
  private roleFooter(): HTMLElement {
    return ui.divV([
      this.kit.prose('Offset = the mean difference of the compounds carrying this value, adjusted for the ' +
        'other components.',
      'Observational, not causal: it is what the compounds that happen to carry this value did, so it ' +
        'is neither a potency nor a prediction.'),
      this.kit.prose(`Fitted over the whole table. The "${PANE_SERIES}" tab counts instead.`,
        `A count on "${PANE_SERIES}" only sees comparisons that were actually made inside one row — a ` +
        'substitution nobody tried side by side is invisible there, and visible here.'),
    ]);
  }

  private buildCoreCard(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const host = this.kit.host;
    const title = `${data.coreRole ?? 'Cores'} — the average of what was made on each`;
    // Not matrix.scores: the preferred score is a raw column mean, which is exactly the statistic the
    // fitted cards were built to replace, and the potency score reads counts captured before the prune.
    const subtitle = 'The mean activity of each core\'s measured compounds, in the activity column\'s ' +
      'own units. Not corrected for which substituents each core was paired with, so a core only ever ' +
      'tried with good groups scores high — but unlike a fitted offset, these do compare between cores.';
    if (!data.coresAreSeries)
      return this.card(title, subtitle, this.kit.reason(CORES_NOT_COMPARABLE));
    if (cores.length === 0)
      return this.card(title, subtitle, this.kit.reason(`No core holds ${MIN_SUPPORT} measured compounds.`));
    const body: HTMLElement[] = [];
    body.push(...cores.map((stat) => this.kit.summaryRow({
      depiction: this.kit.depiction(stat.matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
      name: stat.matrix.label,
      badges: [this.kit.trustDot(stat.matrix)],
      desc: `${count(stat.cpd)} cpd · best ${host.formatActivity(stat.best!.value)} · ` +
        `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}`,
      value: host.formatActivity(stat.typical),
      caption: 'typical',
      valueTip: 'The fitted mean of this core\'s measured compounds, count-weighted. It is the mean of ' +
        'what was made: a core only ever paired with good caps scores high, and a heavily resampled ' +
        'substituent dominates it.',
      tip: `Open ${stat.matrix.label} in the SAR Matrix`,
      onClick: () => host.revealMatrix(stat.matrix),
    })));
    const dir = host.higherIsBetter ? 1 : -1;
    let holder: SeriesStat | null = null;
    // Inside the ranked tier: a compound on a core cut at another depth is not an alternative to the
    // cores listed here, so naming it would set up a comparison the card refuses to make.
    for (const stat of this.kit.rankableCores(data).filter((s) => s.tier === cores[0].tier)) {
      if (stat.best !== null && (holder === null || dir * stat.best.value > dir * holder.best!.value))
        holder = stat;
    }
    // The carry-forward decision nothing else on the screen states: the best single compound and the
    // best scaffold on average are routinely not the same core.
    if (holder !== null && holder.matrix.id !== cores[0].matrix.id) {
      body.push(ui.divText(`${cores[0].matrix.label} runs highest on average. Your single best compound ` +
        `is on ${holder.matrix.label}, whose typical is ${host.formatActivity(holder.typical)}: one ` +
        'lucky substituent, not a better scaffold.', 'chem-sar-cp-hint'));
    }
    return this.card(title, subtitle, body);
  }
}
