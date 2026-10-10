/* The Summary's Effects segment: a tab per component with its full fitted ranking, and a last tab
   of what was counted rather than fitted inside each series. */
import * as ui from 'datagrok-api/ui';
import {MIN_SUPPORT, SeriesStat, SummaryData, SWAP_MIN_PAIRS, SWAP_MIN_SERIES,
  SWAP_ROW_CAP} from './sar-matrix-summary-data';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, count} from '../sar-matrix-ui-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {PANE_EFFECTS, PANE_SERIES, orient, Answer, swapAnchor, roleRanks} from './sar-matrix-summary-common';
import {SummaryRoleCards} from './sar-matrix-summary-role-cards';
import {SummaryRGroupCard} from './sar-matrix-summary-rgroup-card';

export class SummaryEffects {
  private readonly roleCards: SummaryRoleCards;
  private readonly rgroupCard: SummaryRGroupCard;

  constructor(private readonly kit: SummaryPanel) {
    this.roleCards = new SummaryRoleCards(kit, this);
    this.rgroupCard = new SummaryRGroupCard(kit, this);
  }

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
    // Dropped wherever the fit already ranks the axis role, for the same reason as the core card below.
    const series: HTMLElement[] = this.roleAnswer(data, data.axisRole) === null ?
      [this.kit.anchored(this.rgroupCard.build(data), 'rgroup')] : [];
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
        coreNote = this.kit.anchored(ui.divText('Cores are not comparable across series here — a core is ' +
          'one row of one matrix and recurs only inside its own fold lineage. Each series\' best core is ' +
          'in its Start-here expand.', 'chem-sar-cp-hint chem-sar-sum-sub-lead'), 'cores');
      }
    }
    const roles = data.roleFit === null ? [] : data.roleFit.roles.map((role) => role.name);
    if (roles.length === 0)
      series.push(this.kit.anchored(this.buildSwapCard(data), 'swaps'));
    else {
      for (const name of roles)
        series.push(this.kit.anchored(this.buildSwapCard(data, name), swapAnchor(name)));
    }

    const roleCards = this.roleCards.build(data);
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

    const lead = this.roleCards.roleOrdering(data);
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
   * A signed magnitude against the screen's own model error.
   *
   * `scale` is the largest magnitude in the set this bar is compared against — the card's own rows,
   * except where one fit produced the rows of several cards and they are genuinely on one scale. A bar
   * ending inside the grey band is one the analysis cannot resolve.
   */
  effectBar(value: number, scale: number, modelError: number | null): HTMLElement {
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

  /** One role's tile, read off the global fit — the same answer for the axis role and the core role,
   *  since in fragment-columns mode one fit ranks both. Null where that role's card is not ranking. */
  roleAnswer(data: SummaryData, name: string | null): Answer | null {
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

  card(title: string, subtitle: string, body: HTMLElement[], footer?: HTMLElement): HTMLElement {
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

  private buildCoreCard(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const host = this.kit.host;
    const title = `${data.coreRole ?? 'Cores'} — the average of what was made on each`;
    // Not matrix.scores: the preferred score is a raw column mean, which is exactly the statistic the
    // fitted cards were built to replace, and the potency score reads counts captured before the prune.
    const subtitle = 'The mean activity of each core\'s measured compounds, in the activity column\'s ' +
      'own units. Not corrected for which substituents each core was paired with, so a core only ever ' +
      'tried with good groups scores high — but unlike a fitted offset, these do compare between cores.';
    if (cores.length === 0)
      return this.card(title, subtitle, this.kit.reason(`No core holds ${MIN_SUPPORT} measured compounds.`));
    const body = cores.map((stat) => this.kit.summaryRow({
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
    }));
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
