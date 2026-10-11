/* The scale line pinned above every section: the endpoint, its scale and direction, the observed range
   and the model error — the things every ranking under it depends on. */
import * as ui from 'datagrok-api/ui';
import {SummaryData} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryScale {
  constructor(private readonly kit: SummaryPanel) {}

  build(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const transform = host.scalingLabel === 'raw' ? 'untransformed' : `${host.scalingLabel} applied`;
    const parts = [host.activityColumnName || 'Activity', transform,
      host.higherIsBetter ? 'higher is better' : 'lower is better'];
    parts.push(data.minObserved === null ? 'nothing observed' :
      `observed ${host.formatActivity(data.minObserved)} to ${host.formatActivity(data.maxObserved!)}`);
    // Never printed as 0 when no series has a fit: a floor of zero reads as "everything is resolved".
    if (data.modelError !== null)
      parts.push(`± ${host.formatActivity(data.modelError)} model error`);
    // Log-ness of a raw column is inferred from its direction, so the hover says so.
    const assumedLog = host.scalingLabel === 'raw' && host.activityIsLog;
    const line = this.kit.hint(parts.join(' · '), `Every ranking depends on this; ${count(data.measuredCells)} ` +
      'measured cells.' + (!assumedLog ? '' : ' The column is assumed to be a log already.') +
      (data.modelError === null ? '' : ' Model error, not assay error.'));
    line.classList.add('chem-sar-sum-scale');
    // Only after a transform: some raw columns are legitimately negative.
    if (data.maxObserved === null || data.maxObserved >= 0 || host.scalingLabel === 'raw')
      return line;
    return ui.divV([line, ui.divText('Every observed value is negative — check Scaling.',
      'chem-sar-cp-hint chem-sar-chip-partial')]);
  }
}
