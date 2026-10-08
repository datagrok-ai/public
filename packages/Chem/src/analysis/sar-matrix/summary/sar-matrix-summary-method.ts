/* The scale band pinned above every segment, and the Method segment: what passes the trust gate and
   which series the fit holds in. */
import * as ui from 'datagrok-api/ui';
import {SarMatrix} from '../sar-matrix-types';
import {MIN_SUPPORT, SummaryData, TRUST_R2} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, count} from '../sar-matrix-ui-common';
import {PANE_MAKING, TRUST_LIST_MAX} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryMethod {
  constructor(private readonly kit: SummaryPanel) {}

  build(data: SummaryData): HTMLElement {
    const parts: HTMLElement[] = [this.buildChips(data),
      this.kit.anchored(this.buildTrustSection(data), 'trust')];
    const scroll = ui.divV(parts, 'chem-sar-sum-scroll');
    this.kit.scroller = scroll;
    return scroll;
  }

  /**
   * The two things that invalidate every ranking below, and nothing else — this band does not scroll,
   * so anything pinned here is width the answers never get back.
   */
  scaleBand(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const lines: HTMLElement[] = [];
    const line = (text: string): HTMLElement => ui.divText(text, 'chem-sar-sum-orient-line');

    const transform = host.scalingLabel === 'raw' ? 'untransformed' : `${host.scalingLabel} applied`;
    const direction = host.higherIsBetter ? 'higher is better' : 'lower is better';
    // Log-ness of a raw column is inferred from the declared direction, not read off the data, and
    // every fold and log-unit claim on the tab rests on it. Percent inhibition, ΔTm and ΔG are raw and
    // higher-is-better too, so the inference has to be on screen rather than in a getter.
    // Stated on the hover rather than on the band: it qualifies the fold figures further down, not the
    // column.
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
    // Never printed as 0 when no series has a fit: a floor of zero reads as "everything is resolved".
    if (data.modelError !== null)
      head.push(line(`± ${host.formatActivity(data.modelError)} model error`));
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

  /**
   * The observed range as a rule, spanning exactly lo..hi. It must not extend to zero: on a column that
   * is negative throughout — log solubility, ΔG — the right-hand label would sit at the track's end
   * where the value is zero and not the label's number.
   */
  private rangeRule(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const lo = data.minObserved!;
    const hi = data.maxObserved!;
    const span = hi - lo || 1;
    const track = ui.div([ui.div([], 'chem-sar-sum-rule-bar')], 'chem-sar-sum-rule');
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

  private buildChips(data: SummaryData): HTMLElement {
    const items = [
      this.kit.badge(`${count(data.trustedCells)} predictions pass the trust gate · ` +
        `${count(data.trustedStructures)} distinct structures`,
      // The second constant is concatenated rather than interpolated: the bundler folds two adjacent
      // template operands that each carry a compile-time constant and drops the first one's trailing
      // text, which loses a clause without failing the build.
      `Gate: ${MIN_SUPPORT} supporting compounds on both axes, in a series whose own R² ≥ ` +
      TRUST_R2 + '. The distinct-structure count is the actionable one — it deduplicates the same ' +
      'analog predicted in several tiers.'),
    ];
    if (data.trustedNoStructure > 0) {
      items.push(this.kit.badge(`${count(data.trustedNoStructure)} predictions have a value and no structure`,
        `These pass the same trust gate as the analogs in ${PANE_MAKING}, but the core carries an ` +
        'attachment point none of the picked fragment columns fills, so no structure can be completed ' +
        `over it. Add that column to the R-group columns and they appear in ${PANE_MAKING}.`, true));
    }
    return ui.divH(items, 'chem-sar-sum-chips');
  }

  /**
   * What the two R² on this tab mean and what the gate does with them, then the series at each end of it.
   * The poorly-fitting list opens itself when something elsewhere on the tab links here, since that is
   * the one arrival where the reader came for the list rather than the explanation.
   */
  private buildTrustSection(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const parts: HTMLElement[] = [
      ui.divText('Fit quality', 'chem-sar-sum-card-title'),
      // Two quantities with one name on one screen, and nothing else says they are not the same number.
      this.kit.prose('Two different R² appear on this tab: one for the whole table, one for each series.',
        'Both are scored by predicting compounds that were left out of the fit, so they say how well ' +
        'the model predicts rather than how well it describes what it was shown. In a series, 1 means ' +
        'its substituent effects add perfectly, 0 means the model does no better than simply using that ' +
        'series\' average, and below 0 it does worse than that average.'),
      this.kit.prose(`Trust gate: ${MIN_SUPPORT} measured compounds on both axes, in a series whose own R² ` +
        'is at least ' + TRUST_R2 + '.',
      'Predictions that fail it are still computed and still drawn in the matrix. They are left out of ' +
        'the counts above, and out of Worth making unless nothing at all clears the gate.'),
    ];
    const worst = [...data.lowR2].sort((a, b) =>
      (a.confidence!.r2 - b.confidence!.r2) || (a.id < b.id ? -1 : 1));
    if (worst.length === 0) {
      parts.push(this.kit.prose('Every cross-validated fit holds up.',
        'Substituent effects add in each of them, so their predictions can be read as estimates.'));
    } else {
      parts.push(this.foldable(`Where the additive model does not hold · ${count(worst.length)}`,
        `${count(data.lowR2Virtual)} predicted cells inherit a fit that reproduces its own measured ` +
        'cells poorly.', this.trustList(worst), this.kit.openTrustList));
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

  /** The lead line above the fold stays, since a count is worth reading without opening anything. */
  private foldable(title: string, lead: string, body: HTMLElement, open: boolean): HTMLElement {
    return ui.divV([this.kit.foldHead(title, body, open), ui.divText(lead, 'chem-sar-cp-hint'), body]);
  }

  private trustList(matrices: SarMatrix[]): HTMLElement {
    const host = this.kit.host;
    const list = ui.div(matrices.slice(0, TRUST_LIST_MAX).map((matrix) => {
      const conf = matrix.confidence!;
      return this.kit.summaryRow({
        // Drawn: a series name alone does not show which scaffold the fit is about.
        depiction: this.kit.depiction(matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
        name: matrix.label,
        badges: [this.kit.badge(`R² ${conf.r2.toFixed(2)} ± ${host.formatActivity(conf.rmse)}`,
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
}
