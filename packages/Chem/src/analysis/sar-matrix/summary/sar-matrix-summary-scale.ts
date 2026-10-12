import * as ui from 'datagrok-api/ui';
import {SummaryData} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryScale {
  constructor(private readonly kit: SummaryPanel) {}

  build(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const parts = [host.activityColumnName || 'Activity'];
    if (host.scalingLabel !== 'raw')
      parts.push(`${host.scalingLabel} applied`);
    parts.push(host.higherIsBetter ? 'higher is better' : 'lower is better');
    parts.push(data.minObserved === null ? 'nothing observed' :
      `observed ${host.formatActivity(data.minObserved)} to ${host.formatActivity(data.maxObserved!)}`);
    // Log-ness of a raw column is inferred from its direction, so the hover says so.
    const assumedLog = host.scalingLabel === 'raw' && host.activityIsLog;
    const line = this.kit.hint(parts.join(' · '), `Every ranking depends on this; ${count(data.measuredCells)} ` +
      'measured cells.' + (host.scalingLabel === 'raw' ? ' Values are used untransformed.' : '') +
      (!assumedLog ? '' : ' The column is assumed to be a log already.'));
    line.classList.add('chem-sar-sum-scale');
    // Only after a transform: some raw columns are legitimately negative.
    if (data.maxObserved === null || data.maxObserved >= 0 || host.scalingLabel === 'raw')
      return line;
    return ui.divV([line, ui.divText('Every observed value is negative — check Scaling.',
      'chem-sar-cp-hint chem-sar-chip-partial')]);
  }
}
