/* The Overview's Start here list: the series that win one of five reasons to open them — most compounds,
   widest range, best compound, most trusted predictions, best-validated fit — one row per lineage. */
import * as ui from 'datagrok-api/ui';
import {bestMeasured, MIN_SUPPORT, SeriesStat, SummaryData, TRUST_R2} from './sar-matrix-summary-data';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, count, TAB_TRANSFER} from '../sar-matrix-ui-common';
import {GAIN_TICKS, REASON_GLYPHS, REASON_WORDS, shortSmiles} from './sar-matrix-summary-common';
import type {SummaryOverview} from './sar-matrix-summary-overview';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryStartHere {
  constructor(private readonly kit: SummaryPanel, private readonly overview: SummaryOverview) {}

  /**
   * Each row's value: the best unmade compound the model names in that series, against the best one
   * already measured.
   *
   * An additive FLOOR, not a ceiling on the chemistry — a cliff is exactly what this model cannot see,
   * so a gain near zero means "no further additive gain", never "no further potency".
   */
  build(data: SummaryData): HTMLElement[] {
    const host = this.kit.host;
    const parts: HTMLElement[] = [];
    if (data.startHere.length === 0)
      parts.push(...this.kit.reason('No series holds a measured compound.'));
    const lanes = this.overview.lanesFit(data);
    for (const row of data.startHere) {
      const stat = row.primary;
      const {matrix} = stat;
      const rowKey = stat.best !== null ? matrix.rows[stat.best.ri].keySmiles : matrix.rows[0]?.keySmiles;
      const neighbours = stat.best === null ? 0 :
        host.observedNeighbours(matrix, stat.best.ri, stat.best.ci);
      const expand = this.kit.lazyBody(() => this.startExpandParts(stat));
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
        onChevron: expand.toggle,
      }));
      if (!lanes)
        parts.push(ui.divText(row.reasons.map((i) => REASON_WORDS[i]).join(' · '), 'chem-sar-cp-hint'));
      parts.push(expand.body);
    }
    const transfers = ui.divText('', 'chem-sar-cp-hint');
    this.kit.transferLine = transfers;
    this.kit.syncTransferLine();
    transfers.classList.add('chem-sar-sum-click');
    transfers.onclick = () => host.showTab(TAB_TRANSFER);
    parts.push(transfers);
    return parts;
  }

  /** How many of a series' own typical prediction errors a predicted gain is worth. Unfilled
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

  /** The five fixed reasons a series is worth opening, won or not won, always in the same order
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
    const dir = host.higherIsBetter ? 1 : -1;
    const parts: HTMLElement[] = [];

    if (stat.colRange !== null && stat.rowRange !== null) {
      const story = stat.colRange > stat.rowRange ? 'this is a column story' : 'this is a row story';
      parts.push(this.kit.hint(`Spread: substituent choice spans ${stat.colRange.toFixed(2)} · core ` +
        `choice spans ${stat.rowRange.toFixed(2)} — ${story}`, 'Ranges of this series\' own fitted ' +
        `effects, over substituents measured on at least ${MIN_SUPPORT} cores and cores measured at at ` +
        `least ${MIN_SUPPORT} substituents — one floor on both sides, or the looser side of the ` +
        'comparison would come out wider on noise alone. Both are centred on the same fit, so they ' +
        'compare to each other — and to no other series\' numbers.'));
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
      const ci = bestMeasured(matrix.cells[stat.bestRow.ri], dir);
      if (ci >= 0)
        block.onclick = () => host.revealCell(matrix, stat.bestRow!.ri, ci);
      parts.push(block);
    }

    // The one place on this screen where non-additivity is a finding rather than a trust problem, and
    // where the instruction is the opposite of the leaderboards': hold the pair, vary elsewhere.
    const conf = matrix.confidence;
    // `?? null`, not a strict read: a matrix carried in by a project or layout was serialized by
    // whichever Chem wrote it, and one without these keys yields undefined, which `!== null` passes.
    const outlier = (conf ? (dir > 0 ? conf.hi : conf.lo) : null) ?? null;
    if (conf && outlier !== null && Math.abs(outlier.residual) > conf.rmse &&
      outlier.ri < matrix.rows.length && outlier.ci < matrix.columns.length) {
      const beat = this.kit.hint(`${matrix.rows[outlier.ri].label} × ` +
        `${shortSmiles(matrix.columns[outlier.ci].substSmiles)} beats its additive expectation by ` +
        `${this.kit.formatEffect(dir * outlier.residual)}, against this series' own ` +
        `±${host.formatActivity(conf.rmse)} — hold that pair and vary elsewhere.`,
      'Out-of-sample: the cell was held out and predicted from the rest, so the model could not pull ' +
        'itself toward it. It beats the additive sum in this series\' own units; it does not say the ' +
        'mechanism is understood.');
      beat.classList.add('chem-sar-sum-click');
      beat.onclick = () => host.revealCell(matrix, outlier.ri, outlier.ci);
      parts.push(beat);
    }

    const position = matrix.positions[0] ?? '';
    const reference = matrix.refValues[position];
    const refLine = this.kit.hint(reference ? `Varies ${position} · reference substituent:` :
      `Varies ${position} · no reference substituent recorded`, 'The position this series explores, and ' +
      'the most frequently observed substituent at it — the comparator its fitted effects are quoted ' +
      'against. Position labels come from each series\' own decomposition, so R1 here and R1 elsewhere ' +
      'are unrelated.');
    parts.push(reference ?
      ui.divH([refLine, this.kit.depiction(reference, BENEFIT_MOL_W, BENEFIT_MOL_H)], 'chem-sar-sum-detail') :
      refLine);
    return parts;
  }
}
