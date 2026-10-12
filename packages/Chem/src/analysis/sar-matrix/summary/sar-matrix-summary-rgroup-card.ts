import * as ui from 'datagrok-api/ui';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {cachedAtomCount} from '../build/sar-matrix-decompose';
import {RGroupRow, SummaryData, SUM_ROWS, SWAP_MIN_SERIES} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import {shortSmiles} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRGroupCard {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  build(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const title = data.axisRole === null ? 'R-groups — within-series ranking' :
      `${data.axisRole} — within-series ranking`;
    const subtitle = 'Groups that came first at their position in most series that tried them: carry them ' +
      'forward. The number is the margin over each series\' most common substituent.';
    if (!host.activityIsLog) {
      return this.effects.card(title, subtitle, this.kit.reason('The activity is on a raw scale, where a fitted ' +
        'effect is a difference in assay units and says nothing about how large the change is. Set ' +
        'Scaling to lg or −lg, or declare the column higher-is-better if it is already a pIC50.'));
    }
    const rows = data.rgroups.slice(0, SUM_ROWS).map((row) => this.rgroupRow(row, {}));
    const thin = data.rgroupsThin.map((row) => this.rgroupRow(row, {thin: true}));
    if (rows.length === 0 && thin.length === 0) {
      // Non-convergence is named because nothing else on screen shows it.
      return this.effects.card(title, subtitle, this.kit.reason(`No substituent comes first at its position in ` +
        `${SWAP_MIN_SERIES} series whose additive fit holds.` + (data.nonConverged === 0 ? '' :
        ` ${count(data.nonConverged)} series whose fit did not converge are left out.`)));
    }
    const body: HTMLElement[] = [];
    if (this.effects.roleAnswer(data, data.axisRole) !== null) {
      body.push(ui.divText(`This card counts wins within each series; the ${data.axisRole} view reports one ` +
        'effect over the whole table.',
      'chem-sar-sum-prose'));
    }
    if (data.nonConverged > 0) {
      body.push(ui.divText(`${count(data.nonConverged)} further series count in no row: their additive ` +
        'fit did not reach its tolerance, so their effects cannot be ordered against another series\'.',
      'chem-sar-cp-hint'));
    }
    body.push(...rows);
    if (thin.length > 0) {
      body.push(ui.divText(`Won in only two series — a median over two lineages is the mean of two ` +
        'numbers:', 'chem-sar-cp-hint'));
      body.push(...thin);
    }
    if (data.rgroupLosers.length > 0) {
      body.push(ui.divText('Stop making — last at their position in most series that tried them:',
        'chem-sar-cp-hint'));
      body.push(...data.rgroupLosers.map((row) => this.rgroupRow(row, {loser: true})));
    }
    return this.effects.card(title, subtitle, body, this.sizeVerdict(data));
  }

  private tallyBadge(row: RGroupRow, placed: string): HTMLElement {
    return this.kit.badge(`${placed} in ${row.k} of ${row.tried}`,
      `${placed} in ${row.k} of the ${row.tried} series that tried it` + (row.m > row.k ?
        `, plus ${row.m - row.k} whose fit does not hold: likely core-dependent.` : '.'));
  }

  /** One verdict for the card: the total size of the groups matters, not each one's. */
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
      if (n === 0)
        return;
      const won = Math.round(groups / n);
      const beaten = Math.round(references / n);
      line.innerText = won > beaten ? `The winners are bigger than what they beat (${won} vs ${beaten} heavy ` +
        'atoms): part of the gain may be size alone.' :
        `The winners are no bigger than what they beat (${won} vs ${beaten} heavy atoms).`;
    });
    ui.tooltip.bind(line, 'Heavy atoms including the attachment point, averaged over the groups above and the ' +
      'series references they are compared with. Potency alone tends to favour larger groups.');
    return line;
  }

  private rgroupRow(row: RGroupRow, opts: {thin?: boolean, loser?: boolean}): HTMLElement {
    const host = this.kit.host;
    const {win} = row;
    const position = win.matrix.columns[win.ci].position;
    const reference = win.matrix.refValues[position] ?? '';
    const atoms = this.kit.badge('', 'Heavy atoms including the attachment point: potency alone favours size.');
    this.kit.pendingPaints.push(() => {
      const n = cachedAtomCount(row.subst, getRdKitModule());
      atoms.innerText = n ? `${n} atoms` : '—';
    });
    const badges = [this.tallyBadge(row, opts.loser ? 'last' : 'best'), atoms];
    if (opts.thin)
      badges.push(this.kit.badge('2 series', 'Two lineages only — the number below is the mean of two.', true));
    // Only one series' reference can be named; a pooled median has several.
    const comparator = row.magnitudeSeries > 1 ? 'each series\' own reference' : shortSmiles(reference);
    return this.kit.summaryRow({
      depiction: this.kit.depiction(row.subst),
      name: `${win.matrix.label} · ${position}`,
      badges,
      desc: row.magnitude === null ? 'no reference substituent to compare against' :
        `vs ${comparator} · ${this.kit.formatEffect(row.lo)} → ${this.kit.formatEffect(row.hi)} ` +
        `over ${row.magnitudeSeries} series`,
      value: row.magnitude === null ? '—' : this.kit.formatEffect(row.magnitude),
      caption: 'vs reference',
      valueTip: row.magnitude === null ?
        'No contributing series has its reference substituent, so the row ranks on placings alone.' :
        'Median difference against each series\' most common substituent' +
        (row.magnitudeSeries > 1 ? '' : ` (${reference})`) + `, over ${row.magnitudeSeries} of ${row.k} series.`,
      tip: this.kit.openTip(win.matrix, win.ri, win.ci),
      onClick: () => host.revealCell(win.matrix, win.ri, win.ci, position),
    });
  }
}
