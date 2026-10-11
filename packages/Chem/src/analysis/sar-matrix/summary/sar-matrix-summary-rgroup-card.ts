/* The Effects segment's R-group card: the substituents the fitted model most often places first at
   their position, counted over the series that tried them. */
import * as ui from 'datagrok-api/ui';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {cachedAtomCount} from '../build/sar-matrix-decompose';
import {RGroupRow, SummaryData, SUM_ROWS, SWAP_MIN_SERIES} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, count} from '../sar-matrix-ui-common';
import {shortSmiles} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRGroupCard {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  build(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const title = data.axisRole === null ? 'R-groups — within-series ranking' :
      `${data.axisRole} — within-series ranking`;
    const subtitle = 'How often the fit places a group first at its position; the number is its margin over ' +
      'that series\' most common substituent.';
    // An additive effect in raw assay units has no fold to state.
    if (!host.activityIsLog) {
      return this.effects.card(title, subtitle, this.kit.reason('The activity is on a raw scale, where a fitted ' +
        'effect is a difference in assay units and says nothing about how large the change is. Set ' +
        'Scaling to lg or −lg, or declare the column higher-is-better if it is already a pIC50.'));
    }
    const rows = data.rgroups.slice(0, SUM_ROWS).map((row, i) =>
      this.rgroupBlock(data, row, {primary: i === 0}));
    const thin = data.rgroupsThin.map((row) => this.rgroupBlock(data, row, {thin: true}));
    if (rows.length === 0 && thin.length === 0) {
      // Non-convergence is named because nothing else on screen shows it.
      return this.effects.card(title, subtitle, this.kit.reason(`No substituent comes first at its position in ` +
        `${SWAP_MIN_SERIES} series whose additive fit holds.` + (data.nonConverged === 0 ? '' :
        ` ${count(data.nonConverged)} series whose fit did not converge are left out.`)));
    }
    const body: HTMLElement[] = [];
    // Where the cross-table fit also ranks this column, say how the two numbers differ.
    if (this.effects.roleAnswer(data, data.axisRole) !== null) {
      body.push(ui.divText(`This card counts wins within each series; the ${data.axisRole} view reports one ` +
        'offset over the whole table.',
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
      const head = data.rgroupLosers[0];
      body.push(ui.divText(`Stop making: ${shortSmiles(head.subst)} came last at its position in ` +
        `${head.k} of the ${head.tried} series it was tried in` +
        (head.magnitude === null ? '.' : `, costing ${this.kit.formatEffect(head.magnitude)} against each ` +
        'series\' own reference.'), 'chem-sar-cp-hint'));
      body.push(...data.rgroupLosers.map((row) => this.rgroupBlock(data, row, {loser: true})));
    }
    return this.effects.card(title, subtitle, body, this.sizeVerdict(data));
  }

  /** Out of the series that tried the group, not out of its wins. */
  private tallyBadge(row: RGroupRow, placed: string): HTMLElement {
    return this.kit.badge(`${placed} in ${row.k} of ${row.tried}`,
      `${placed} in ${row.k} of the ${row.tried} series that tried it` + (row.m > row.k ?
        `, plus ${row.m - row.k} whose fit does not hold: likely core-dependent.` : '.'));
  }

  /** One leaderboard entry: the row, and an expand that says what the row does and does not claim. */
  private rgroupBlock(data: SummaryData, row: RGroupRow,
    opts: {thin?: boolean, loser?: boolean, primary?: boolean}): HTMLElement {
    const parts: HTMLElement[] = [this.rgroupRow(row, opts)];
    const expand = this.kit.lazyBody(() => [this.rgroupClaim(data, row)]);
    const toggle = ui.divText('What this row claims →', 'chem-sar-cp-hint');
    toggle.classList.add('chem-sar-sum-click');
    toggle.onclick = expand.toggle;
    // The answer tile above asked for its evidence, and the pane it asked for is only built now.
    if (opts.primary && this.kit.expandTopRGroup) {
      this.kit.expandTopRGroup = false;
      expand.open();
    }
    parts.push(toggle, expand.body);
    return ui.divV(parts);
  }

  /** A reading instruction, not a restatement of the marks above it. */
  private rgroupClaim(data: SummaryData, row: RGroupRow): HTMLElement {
    // Only where a matrix IS one core is a series a core; with a series column one series holds several.
    const slot = data.coresAreSeries ? 'core' : 'series';
    const read = ui.divText(`First on most ${slot}s it was tried on: carry it forward. First on one ${slot} ` +
      'only: keep it there and vary elsewhere.',
    'chem-sar-sum-prose');
    return ui.divH([this.kit.depiction(row.subst, CARD_CORE_W * 2, CARD_CORE_H * 2), read], 'chem-sar-sum-pair');
  }

  /** One size verdict for the card: the sum of the groups is what matters, not each one. */
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
    ui.tooltip.bind(line, 'Heavy atoms including the attachment point: a size warning, not a property model.');
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
      depiction: this.kit.depiction(row.subst, CARD_CORE_W, CARD_CORE_H),
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
