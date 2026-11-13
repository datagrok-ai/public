import * as ui from 'datagrok-api/ui';
import {MIN_SUPPORT, SeriesStat, SummaryData, TRUST_R2} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import {REASON_GLYPHS, REASON_WORDS} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryStartHere {
  constructor(private readonly kit: SummaryPanel) {}

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
        depiction: this.kit.depiction(rowKey ?? null),
        name: matrix.label,
        badges: [
          this.kit.badge(`L${stat.tier}`, stat.tier === 1 ?
            'A leaf series: nothing finer under it.' :
            `One cut coarser than the L${stat.tier - 1} series below it.`),
          ...this.kit.trustBadges(matrix),
          ...this.gainBadge(stat),
        ],
        mark: this.reasonLane(row.reasons),
        // One line even in a half-width column; the denominator and the compound count are in the tooltip.
        desc: `${count(stat.realCells)} measured · ` +
          `${host.formatActivity(matrix.minActivity)}–${host.formatActivity(matrix.maxActivity)}` +
          (stat.realCells ? ` · typical ${host.formatActivity(stat.typical)}` : ''),
        value: stat.best === null ? '—' : host.formatActivity(stat.best.value),
        caption: 'best',
        valueTip: stat.best === null ? 'This series holds no measured compound.' :
          `${host.formatActivity(stat.best.value)} measured; ${neighbours} measured cells share its core ` +
          'or substituent.',
        // Cells the decomposition cannot express are left out of the denominator.
        tip: (matrix.realCount === stat.realCells ? '' : `${count(matrix.realCount)} compounds; `) +
          `${count(stat.realCells)} of ` +
          `${count(stat.totalCells - stat.impossibleCells)} ${stat.impossibleCells > 0 ? 'makeable cells' : 'cells'} ` +
          `measured. Click to expand; the grid icon opens ${matrix.label} in the SAR Matrix. "Typical" is the fitted ` +
          'mean of what was made; a big gap to "best" suggests one lucky compound.',
        onClick: () => stat.best === null ? host.revealMatrix(matrix) :
          host.revealCell(matrix, stat.best.ri, stat.best.ci),
        onChevron: expand.toggle,
      }));
      parts.push(expand.body);
    }
    return parts;
  }

  /** Every reason in a fixed slot, won or not, so the lanes compare down the list. */
  private reasonLane(reasons: number[]): HTMLElement {
    const won = new Set(reasons);
    const lane = ui.divH(REASON_GLYPHS.map((glyph, i) => ui.divText(won.has(i) ? glyph : '',
      `chem-sar-sum-reason${won.has(i) ? ' chem-sar-sum-reason-on' : ''}`)), 'chem-sar-sum-lane');
    ui.tooltip.bind(lane, () => reasons.map((i) => REASON_WORDS[i]).join(' · '));
    return lane;
  }

  /** Null where the trust badge beside it already says why. */
  private noGainReason(stat: SeriesStat): {label: string, tip: string} | null {
    const r2 = stat.matrix.confidence?.r2 ?? null;
    if (!this.kit.host.predictVirtual) {
      return {label: 'predictions off', tip: '"Predict virtual analogs" is off, so this series has no ' +
        'predicted cells to compare against. Turn it on in the viewer\'s properties.'};
    }
    if (stat.virtualCells === 0) {
      return {label: 'nothing unfilled', tip: 'Every cell this series can express is already measured ' +
        'or cannot be made, so there is no unfilled cell to predict.'};
    }
    if (r2 === null || r2 < TRUST_R2)
      return null;
    return {label: 'thin support', tip: 'The fit holds, but every prediction rests on fewer than ' +
      MIN_SUPPORT + ' measured compounds on an axis.'};
  }

  private gainBadge(stat: SeriesStat): HTMLElement[] {
    const host = this.kit.host;
    if (stat.bestVirtual === null || stat.best === null) {
      const reason = this.noGainReason(stat);
      return reason === null ? [] : [this.kit.badge(reason.label, reason.tip)];
    }
    const reach = this.kit.dir * (stat.bestVirtual.value - stat.best.value);
    const unsure = host.activityIsLog && Math.abs(reach) < (stat.matrix.confidence?.rmse ?? 0);
    return [this.kit.badge(`${this.kit.formatEffect(reach)} predicted`, 'Best trusted prediction for an unmade ' +
      'analog, minus the best compound made.' +
      (unsure ? ' Smaller than this series\' error, so not a sure gain.' : ''))];
  }

  private startExpandParts(stat: SeriesStat): HTMLElement[] {
    if (stat.bestRow === null)
      return [];
    const {matrix} = stat;
    const row = matrix.rows[stat.bestRow.ri];
    const folded = Object.values(row.foldedValues);
    const oneStructure = matrix.rows.length > 1 &&
      new Set(matrix.rows.map((r) => r.keySmiles)).size === 1;
    const heading = folded.length === 0 ? 'Best core' :
      oneStructure ? 'Best R-group combination' : 'Best core + fixed R-groups';
    // Where every row draws the same core, the folded substituents are what tell rows apart.
    const art = folded.length > 0 && oneStructure ?
      ui.divH(folded.slice(0, 3).map((subst) => this.kit.depiction(subst)), 'chem-sar-sum-pair') :
      this.kit.depiction(row.keySmiles);
    const text = ui.divV([
      ui.divText(`${heading}: ${row.label} · ${this.kit.formatEffect(stat.bestRow.effect)} ` +
        `(n = ${stat.bestRow.n})`, 'chem-sar-sum-prose'),
      ui.divText('Best on average over its measured substituents — not the single best compound, and not ' +
        'comparable across series.', 'chem-sar-cp-hint'),
    ]);
    return [ui.divH([art, text], 'chem-sar-sum-pair chem-sar-sum-expand')];
  }
}
