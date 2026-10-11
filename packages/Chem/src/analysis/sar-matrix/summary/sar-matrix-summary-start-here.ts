/* The Overview's Start here list: the series that win one of five reasons to open them — most compounds,
   widest range, best compound, most trusted predictions, best-validated fit — one row per lineage. */
import * as ui from 'datagrok-api/ui';
import {bestMeasured, MIN_SUPPORT, SeriesStat, SummaryData, TRUST_R2} from './sar-matrix-summary-data';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, count} from '../sar-matrix-ui-common';
import {REASON_GLYPHS, REASON_WORDS, shortSmiles} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryStartHere {
  constructor(private readonly kit: SummaryPanel) {}

  /** Each row's gain is the best unmade compound against the best measured one: an additive floor, so near
   *  zero means no further additive gain, never no further potency. */
  build(data: SummaryData): HTMLElement[] {
    const host = this.kit.host;
    const parts: HTMLElement[] = [];
    if (data.startHere.length === 0)
      parts.push(...this.kit.reason('No series holds a measured compound.'));
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
            'A leaf series: nothing finer under it.' :
            `One cut coarser than the L${stat.tier - 1} series below it.`),
          this.kit.trustBadge(matrix),
          this.gainBadge(stat),
        ],
        mark: this.reasonLane(row.reasons),
        // Cells the decomposition cannot express are left out of the denominator.
        desc: `${count(stat.cpd)} cpd · ${count(stat.realCells)} of ` +
          `${count(stat.totalCells - stat.impossibleCells)} ` +
          `${stat.impossibleCells > 0 ? 'makeable cells' : 'cells'} measured · ` +
          `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}` +
          (stat.realCells ? ` · typical ${host.formatActivity(stat.typical)}` : ''),
        value: stat.best === null ? '—' : host.formatActivity(stat.best.value),
        caption: 'best',
        valueTip: stat.best === null ? 'This series holds no measured compound.' :
          `${host.formatActivity(stat.best.value)} measured; ${neighbours} measured cells share its core ` +
          'or substituent.',
        tip: `Open ${matrix.label} in the SAR Matrix. "Typical" is the fitted mean of what was made; a big ` +
          'gap to "best" suggests one lucky compound.',
        onClick: () => stat.best === null ? host.revealMatrix(matrix) :
          host.revealCell(matrix, stat.best.ri, stat.best.ci),
        onChevron: expand.toggle,
      }));
      parts.push(expand.body);
    }
    return parts;
  }

  /** All five reasons in a fixed order, won or not, so the lanes compare down the list. */
  private reasonLane(reasons: number[]): HTMLElement {
    const won = new Set(reasons);
    const lane = ui.divH(REASON_GLYPHS.map((glyph, i) => ui.divText(won.has(i) ? glyph : '',
      `chem-sar-sum-reason${won.has(i) ? ' chem-sar-sum-reason-on' : ''}`)), 'chem-sar-sum-lane');
    ui.tooltip.bind(lane, () => reasons.map((i) => REASON_WORDS[i]).join(' · '));
    return lane;
  }

  /** Which part of the gate left this series with nothing to aim at; each calls for a different fix. */
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
      return {label: 'fit unchecked', tip: `${count(stat.virtualCells)} predictions, but too few cells to ` +
        'check the fit: unchecked, not wrong.'};
    }
    if (r2 < TRUST_R2) {
      return {label: 'fit fails', tip: `R² ${r2.toFixed(2)}, below ${TRUST_R2}: effects do not add up here, ` +
        'so the predictions are not estimates.'};
    }
    return {label: 'thin support', tip: 'The fit holds, but every prediction rests on fewer than ' +
      MIN_SUPPORT + ' measured compounds on an axis.'};
  }

  private gainBadge(stat: SeriesStat): HTMLElement {
    const host = this.kit.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const reach = stat.bestVirtual !== null && stat.best !== null ?
      dir * (stat.bestVirtual.value - stat.best.value) : null;
    let tip = 'Best prediction vs best measured. Negative suggests a cliff; near zero means no further ' +
      'additive gain, not an exhausted series.';
    const rmse = stat.matrix.confidence?.rmse ?? null;
    if (reach !== null && rmse !== null && rmse > 0 && host.activityIsLog) {
      tip += ` That is ${(reach / rmse).toFixed(1)}× this series' own leave-one-out error; under 1× the ` +
        'model cannot tell it from the best compound already made.';
    }
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
    return badge;
  }

  private startExpandParts(stat: SeriesStat): HTMLElement[] {
    const host = this.kit.host;
    const {matrix} = stat;
    const dir = host.higherIsBetter ? 1 : -1;
    const parts: HTMLElement[] = [];

    if (stat.colRange !== null && stat.rowRange !== null) {
      const story = stat.colRange > stat.rowRange ? 'this is a column story' : 'this is a row story';
      parts.push(this.kit.hint(`Spread: substituent choice spans ${stat.colRange.toFixed(2)} · core ` +
        `choice spans ${stat.rowRange.toFixed(2)} — ${story}`, 'Ranges of this series\' fitted effects, ' +
        'each over values with at least ' + MIN_SUPPORT + ' measurements; comparable only within the series.'));
    }

    if (stat.bestRow !== null) {
      const row = matrix.rows[stat.bestRow.ri];
      const folded = Object.values(row.foldedValues);
      const oneStructure = matrix.rows.length > 1 &&
        new Set(matrix.rows.map((r) => r.keySmiles)).size === 1;
      const heading = folded.length === 0 ? 'Best core' :
        oneStructure ? 'Best R-group combination' : 'Best core + fixed R-groups';
      // Where every row draws the same core, the folded substituents are what tell rows apart.
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
      const block = ui.divH([art, text], 'chem-sar-sum-pair');
      block.classList.add('chem-sar-sum-click');
      const ci = bestMeasured(matrix.cells[stat.bestRow.ri], dir);
      if (ci >= 0)
        block.onclick = () => host.revealCell(matrix, stat.bestRow!.ri, ci);
      parts.push(block);
    }

    // A cell well above its additive prediction: here non-additivity is a finding, not a trust problem.
    const conf = matrix.confidence;
    // `?? null`: a matrix restored from an older project may lack these keys.
    const outlier = (conf ? (dir > 0 ? conf.hi : conf.lo) : null) ?? null;
    if (conf && outlier !== null && Math.abs(outlier.residual) > conf.rmse &&
      outlier.ri < matrix.rows.length && outlier.ci < matrix.columns.length) {
      const beat = this.kit.hint(`${matrix.rows[outlier.ri].label} × ` +
        `${shortSmiles(matrix.columns[outlier.ci].substSmiles)} beats its additive expectation by ` +
        `${this.kit.formatEffect(dir * outlier.residual)}, against this series' own ` +
        `±${host.formatActivity(conf.rmse)} — hold that pair and vary elsewhere.`,
      'Predicted with this cell held out, so the model could not pull toward it.');
      beat.classList.add('chem-sar-sum-click');
      beat.onclick = () => host.revealCell(matrix, outlier.ri, outlier.ci);
      parts.push(beat);
    }

    const position = matrix.positions[0] ?? '';
    const reference = matrix.refValues[position];
    const refLine = this.kit.hint(reference ? `Varies ${position} · reference substituent:` :
      `Varies ${position} · no reference substituent recorded`, 'The most common substituent, the comparator ' +
      'for this series\' effects. R-labels are local to each series.');
    parts.push(reference ?
      ui.divH([refLine, this.kit.depiction(reference, BENEFIT_MOL_W, BENEFIT_MOL_H)], 'chem-sar-sum-pair') :
      refLine);
    return parts;
  }
}
