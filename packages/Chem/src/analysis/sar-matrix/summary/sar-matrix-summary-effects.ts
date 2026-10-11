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

  /** One view per component column, then the per-series evidence, picked one at a time so the section
   *  does not grow with the number of components. */
  build(data: SummaryData): HTMLElement {
    // Dropped wherever the fit already ranks the axis role, for the same reason as the core card below.
    const series: HTMLElement[] = this.roleAnswer(data, data.axisRole) === null ?
      [this.kit.anchored(this.rgroupCard.build(data), 'rgroup')] : [];
    // Dropped where the fit already ranks the core column: two rankings of one column read as a conflict.
    let coreNote: HTMLElement | null = null;
    if (this.roleAnswer(data, data.coreRole) === null) {
      if (data.coresAreSeries)
        series.push(this.kit.anchored(this.buildCoreCard(data, this.kit.topCores(data)), 'cores'));
      else {
        // No ranking is possible here; a note keeps the anchor the core tile points at.
        coreNote = this.kit.anchored(ui.divText('Cores do not compare across series here; each series\' best ' +
          'core is in its Start here expand.', 'chem-sar-cp-hint'), 'cores');
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
    const labels = tabs.map((tab) => tab.note === null ? tab.label : `${tab.label} (${tab.note})`);
    const show = ui.input.choice('Show', {value: labels[open], items: labels, nullable: false,
      tooltipText: 'A card per component, ranked by the additive fit over the whole table, then what was ' +
        'counted rather than fitted inside each series. The number beside a component is its spread.',
      onValueChanged: (value: string) => {
        this.kit.effectsTab = labels.indexOf(value);
        this.kit.showPane(PANE_EFFECTS);
      }});

    const lead = this.roleCards.roleOrdering(data);
    show.root.classList.add('chem-sar-sum-show');
    return ui.divV([...(lead === null ? [] : [lead]), show.root, this.effectsScroll(tabs[open].cards)],
      'chem-sar-sum-effects');
  }

  private effectsScroll(cards: HTMLElement[]): HTMLElement {
    return ui.divV(cards, 'chem-sar-sum-scroll');
  }

  /** One component's tile from the cross-table fit; null where its card does not rank. */
  roleAnswer(data: SummaryData, name: string | null): Answer | null {
    const fit = data.roleFit;
    if (fit === null || this.kit.roleFitRefusal(data) !== null)
      return null;
    const role = fit.roles.find((other) => other.name === name);
    if (role === undefined || !roleRanks(role))
      return null;
    const [top, second] = role.levels;
    // Signed by direction: on a lower-is-better column the leader's coefficient is the smaller.
    const dir = this.kit.host.higherIsBetter ? 1 : -1;
    if (!this.kit.leads(dir * top.coef, dir * second.coef, fit.cvRmse)) {
      // When the top two are inseparable, a two-group split is the answer.
      const split = role.split;
      if (split !== null) {
        return {answer: `two groups: ${count(split.hiCount)} at ${this.kit.formatEffect(split.hiMean)}, ` +
          `${count(split.loCount)} at ${this.kit.formatEffect(split.loMean)}`,
        negative: `no single ${role.name} leads; the groups are ${split.gap.toFixed(2)} apart, so the group ` +
          'matters more than the choice inside it'};
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
    const parts = [ui.divText(title, 'chem-sar-cp-section-title'),
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
    const subtitle = 'Measured pairs alike in everything but one R-group. Nothing here is fitted or predicted.';
    const unit = role === undefined ? 'series' : 'contexts';
    const rows = pools.map((pool) => {
      const {forward, from, to, worst, widest: bestOf, mean, up} = orient(pool);
      const allUp = up === pool.n;
      // Null for a pool grouped off the component columns: no single cell to open.
      const best = pool.best;
      const ci = best === null ? -1 : forward ? best.ciTo : best.ciFrom;
      const position = best === null ? '' : best.matrix.columns[ci].position;
      const show = (value: number): string => this.kit.formatDelta(value);
      // Badges, not descriptor text, which clips to one line.
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
      // Only on a log scale, where the delta and the model error share units.
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
        valueTip: `At least ${show(worst)} in all ${pool.n} measured pairs, across ${pool.roots.size} ${unit}; ` +
          'no fit involved. At 3 pairs "all up" happens one time in four by chance.' +
          (pool.sampled ? ` A group over ${SWAP_ROW_CAP} compounds kept only its extremes.` : ''),
        // No cell to open, so the row selects the compounds on the better side.
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
    // Not matrix.scores, which is a raw column mean read before pruning.
    const subtitle = 'Mean activity of each core\'s measured compounds, not adjusted for its substituents: a ' +
      'core only tried with good groups scores high. Unlike fitted offsets, these compare between cores.';
    if (cores.length === 0)
      return this.card(title, subtitle, this.kit.reason(`No core holds ${MIN_SUPPORT} measured compounds.`));
    const body = cores.map((stat) => this.kit.summaryRow({
      depiction: this.kit.depiction(stat.matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
      name: stat.matrix.label,
      badges: [this.kit.trustBadge(stat.matrix)],
      desc: `${count(stat.cpd)} cpd · best ${host.formatActivity(stat.best!.value)} · ` +
        `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}`,
      value: host.formatActivity(stat.typical),
      caption: 'typical',
      valueTip: 'Count-weighted fitted mean of what was made on this core; good caps lift it.',
      tip: `Open ${stat.matrix.label} in the SAR Matrix`,
      onClick: () => host.revealMatrix(stat.matrix),
    }));
    const dir = host.higherIsBetter ? 1 : -1;
    let holder: SeriesStat | null = null;
    // Inside the ranked tier only: other depths do not compare.
    for (const stat of this.kit.rankableCores(data).filter((s) => s.tier === cores[0].tier)) {
      if (stat.best !== null && (holder === null || dir * stat.best.value > dir * holder.best!.value))
        holder = stat;
    }
    // The best single compound and the best core on average are often not the same core.
    if (holder !== null && holder.matrix.id !== cores[0].matrix.id) {
      body.push(ui.divText(`${cores[0].matrix.label} runs highest on average. Your single best compound ` +
        `is on ${holder.matrix.label}, whose typical is ${host.formatActivity(holder.typical)}: one ` +
        'lucky substituent, not a better scaffold.', 'chem-sar-cp-hint'));
    }
    return this.card(title, subtitle, body);
  }
}
