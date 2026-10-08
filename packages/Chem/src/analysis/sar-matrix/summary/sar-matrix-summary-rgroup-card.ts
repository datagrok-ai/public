/* The Effects segment's R-group card: the substituents the fitted model most often places first at
   their position, counted over the series that tried them. */
import * as ui from 'datagrok-api/ui';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {cachedAtomCount} from '../build/sar-matrix-decompose';
import {bestMeasured, RGroupRow, SeriesStat, SummaryData, SUM_ROWS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, count, tipText} from '../sar-matrix-ui-common';
import {STRIP_MIN_FILL, shortSmiles} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRGroupCard {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  build(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const title = data.axisRole === null ? 'R-groups — within-series ranking' :
      `${data.axisRole} — within-series ranking`;
    const subtitle = 'Ranked by how often the fitted model places this group first at its position. The ' +
      'number is its margin over that series\' own most-common substituent — a within-series ' +
      'comparison, so two series\' numbers are only loosely comparable.';
    // An additive effect in raw assay units cannot be re-expressed as a fold, which is the one thing
    // that would make it mean something.
    if (!host.activityIsLog) {
      return this.effects.card(title, subtitle, this.kit.reason('The activity is on a raw scale, where a fitted ' +
        'effect is a difference in assay units and says nothing about how large the change is. Set ' +
        'Scaling to lg or −lg, or declare the column higher-is-better if it is already a pIC50.'));
    }
    const slots = data.series.filter((stat) => stat.stripCols !== null);
    // One scale for every bar on this card, so the bars compare down the card — and, deliberately, to
    // nothing on any other card.
    const shown = [...data.rgroups.slice(0, SUM_ROWS), ...data.rgroupsThin, ...data.rgroupLosers];
    const scale = shown.reduce((m, r) => Math.max(m, Math.abs(r.magnitude ?? 0)), 0);
    // With one row the bar is full-length and says nothing.
    const bars = shown.filter((r) => r.magnitude !== null).length > 1 ? scale : 0;
    // One strip scale for the whole card, like the bars: normalised per row, every row's own best
    // difference draws at full height, so a +1.8 row and a +0.05 row terminate alike.
    let stripScale = 0;
    for (const stat of slots) {
      for (const r of shown) {
        const delta = stat.stripCols!.get(r.subst)?.refDelta;
        if (delta != null)
          stripScale = Math.max(stripScale, Math.abs(delta));
      }
    }
    const rows = data.rgroups.slice(0, SUM_ROWS).map((row, i) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {primary: i === 0}));
    const thin = data.rgroupsThin.map((row) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {thin: true}));
    if (rows.length === 0 && thin.length === 0) {
      // Convergence is named because it is the one cause the reader cannot see anywhere else: a fit
      // that stopped short is silently excluded from every pool, and R² says nothing about it.
      return this.effects.card(title, subtitle, this.kit.reason(`No substituent comes first at its position in ` +
        `${SWAP_MIN_SERIES} series whose additive fit holds.` + (data.nonConverged === 0 ? '' :
        ` The additive fit of ${count(data.nonConverged)} of the ${count(this.kit.host.matrices.length)} ` +
        'series did not reach its tolerance, and those are left out of this comparison.')));
    }
    const body: HTMLElement[] = [];
    // Two numbers for one substituent are a tab apart, and nothing else says which question each
    // answers. Only where that second number exists: where the fit declined to rank this column there
    // is no other card to reconcile with.
    if (this.effects.roleAnswer(data, data.axisRole) !== null) {
      body.push(ui.divText('A value can place first in most series and still sit at the library mean: ' +
        'this card counts wins within a series against that series\' own reference, while the ' +
        `${data.axisRole} tab reports one offset pooled over the table.`,
      'chem-sar-sum-prose'));
    }
    if (slots.length > 0 && data.nonConverged > 0) {
      body.push(ui.divText(`${count(data.nonConverged)} further series carry no square and count in ` +
        'no row: their additive fit did not reach its tolerance, so their effects cannot be ordered ' +
        'against another series\'.', 'chem-sar-cp-hint'));
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
      body.push(...data.rgroupLosers.map((row) =>
        this.rgroupBlock(data, slots, row, bars, stripScale, {loser: true})));
    }
    return this.effects.card(title, subtitle, body, this.sizeVerdict(data));
  }

  /** One slot per series that tried this group — solid where it won and that series' fit holds,
   *  hatched where it won and the fit does not, outlined where it was tried and did not win. */
  private winTally(row: RGroupRow, placed: string): HTMLElement {
    const shown = Math.min(row.tried, 10);
    const slots: HTMLElement[] = [];
    for (let i = 0; i < shown; i++) {
      const state = i < row.k ? ' chem-sar-sum-slot-on' :
        i < row.m ? ' chem-sar-sum-slot-hatch' : '';
      const slot = ui.div([], `chem-sar-sum-wslot${state}`);
      if (i >= row.k && i < row.m) {
        slot.onclick = (e: MouseEvent) => {
          e.stopPropagation();
          this.kit.revealTrust();
        };
      }
      slots.push(slot);
    }
    const parts: HTMLElement[] = [ui.divH(slots, 'chem-sar-sum-tally')];
    if (row.tried > shown)
      parts.push(ui.divText(`+${row.tried - shown}`, 'chem-sar-sum-tally-more'));
    const box = ui.divH(parts, 'chem-sar-sum-tally-box');
    ui.tooltip.bind(box, () => `${placed} in ${row.k} of the ${row.tried} series that tried it` +
      (row.m > row.k ? `, and in ${row.m - row.k} more whose additive fit does not hold — a group that ` +
      `is ${placed} in additive series and not in non-additive ones is core-dependent. Click a hatched ` +
      'slot to see those series.' : '. Every fit behind those placings holds.'));
    return box;
  }

  /** What one strip slot is. A strip is drawn for any fragment-columns run, but only where a matrix IS
   *  one core is a slot a core — with a series column, or with nothing on the row axis, one slot holds
   *  several cores. */
  private slotNoun(data: SummaryData): string {
    return data.coresAreSeries ? 'core' : 'series';
  }

  /** One leaderboard entry: the row, its per-slot strip where one can be drawn, and an expand that
   *  says what the row does and does not claim. */
  private rgroupBlock(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    stripScale: number, opts: {thin?: boolean, loser?: boolean, primary?: boolean}): HTMLElement {
    const parts: HTMLElement[] = [this.rgroupRow(data, slots, row, scale, opts)];
    if (slots.length > 0)
      parts.push(this.buildStrip(data, slots, row, stripScale));
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
    const slot = this.slotNoun(data);
    const read = ui.divText(`How to read this: first on most of the ${slot}s it was tried on and never ` +
      `worse than the reference means the group travels — carry it forward. First on one ${slot} while ` +
      `costing potency on the rest means hold it to that ${slot} and vary elsewhere.`,
    'chem-sar-sum-prose');
    return ui.divH([this.kit.depiction(row.subst, CARD_CORE_W * 2, CARD_CORE_H * 2), read], 'chem-sar-sum-detail');
  }

  /**
   * One square per slot, in a fixed order that is never re-sorted per row — the row is only
   * readable against its neighbours if the same slot means the same series everywhere.
   *
   * Sign is carried by the direction the fill grows from the mid-rule, and only redundantly by hue:
   * a tint alone renders +0.9 and −0.9 as the same grey. The mid-rule is the legend.
   *
   * `widest` is the card's largest difference, not this row's, so fill heights compare down the card
   * the way the slot order lets the squares compare across it.
   */
  private buildStrip(data: SummaryData, slots: SeriesStat[], row: RGroupRow,
    widest: number): HTMLElement {
    const host = this.kit.host;
    const slot = this.slotNoun(data);
    const squares = slots.map((stat) => {
      const col = stat.stripCols!.get(row.subst);
      if (col === undefined) {
        return tipText('–', 'chem-sar-sum-sq chem-sar-sum-sq-none',
          `${stat.matrix.label} — never tried on this ${slot}`);
      }
      const square = ui.div([], 'chem-sar-sum-sq');
      const delta = col.refDelta;
      if (delta === null) {
        square.classList.add('chem-sar-sum-sq-hatch');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — tried here over ${col.n} measured cells, ` +
          'but this series records no reference substituent to compare against');
      } else if (delta === 0 && stat.matrix.refValues[stat.matrix.columns[col.ci].position] === row.subst) {
        square.classList.add('chem-sar-sum-sq-ref');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — this group IS this series' reference ` +
          `substituent, over ${col.n} measured cells, so there is nothing to compare it against here`);
      } else {
        square.classList.add('chem-sar-sum-sq-val');
        // No difference draws no fill: the floor exists so a small value still reads as a direction,
        // and applied to an exact zero it would draw a win that is not there.
        if (delta !== 0) {
          const fill = ui.div([], `chem-sar-sum-fill chem-sar-sum-${delta > 0 ? 'up' : 'down'}`);
          fill.style.height = `${Math.max(STRIP_MIN_FILL,
            Math.round((widest > 0 ? Math.abs(delta) / widest : 0) * 9))}px`;
          square.appendChild(fill);
        }
        ui.tooltip.bind(square, () => `${stat.matrix.label} · ${this.kit.formatEffect(delta)} against this ` +
          `series' own reference, over ${col.n} measured cells` + (delta === 0 ?
          ' — no difference at all, so nothing grows from the mid-rule.' :
          ` — ${delta > 0 ? 'better than' : 'worse than'} it, so the fill grows ` +
          `${delta > 0 ? 'up' : 'down'} from the mid-rule. Heights compare down this card, against its ` +
          `widest difference (${this.kit.formatEffect(widest)}).`));
      }
      square.classList.add('chem-sar-sum-click');
      square.onclick = () => {
        const ri = bestMeasured(stat.matrix.cells.map((cells) => cells[col.ci]), host.higherIsBetter ? 1 : -1);
        if (ri >= 0)
          host.revealCell(stat.matrix, ri, col.ci, stat.matrix.columns[col.ci].position);
      };
      return square;
    });
    return ui.divH(squares, 'chem-sar-sum-strip');
  }

  /** One verdict for the whole card, not a number per row: five heavy-atom counts do not add
   *  themselves up in a reader's head, and it is the sum that ends a series at MW 620. */
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
    ui.tooltip.bind(line, 'Heavy atoms of the fragment, the attachment point included — not ' +
      'molecular weight, and not lipophilicity, which is what actually kills the compound. A smoke ' +
      'alarm, not a property model.');
    return line;
  }

  private rgroupRow(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    opts: {thin?: boolean, loser?: boolean}): HTMLElement {
    const host = this.kit.host;
    const {win} = row;
    const position = win.matrix.columns[win.ci].position;
    const reference = win.matrix.refValues[position] ?? '';
    const atoms = this.kit.badge('', 'Heavy atoms of the fragment, the attachment point included. The ' +
      'model scores potency alone and carries no efficiency term, so an unlabelled leaderboard is a ' +
      'heavy-atom leaderboard wearing a potency label.');
    this.kit.pendingPaints.push(() => {
      const n = cachedAtomCount(row.subst, getRdKitModule());
      atoms.innerText = n ? `${n} atoms` : '—';
    });
    // Attempts, not wins: "best in 5 of 6" is heard as "tried in 6", and a denominator of wins makes
    // every row on the card read as near-perfect.
    const placed = opts.loser ? 'last' : 'best';
    const badges = [this.winTally(row, placed), atoms];
    if (opts.thin)
      badges.push(this.kit.badge('2 series', 'Two lineages only — the number below is the mean of two.', true));
    // Only a single contributing series may have its reference named: each series is measured against
    // its own most-common substituent, and those genuinely differ, so one SMILES cannot stand for the
    // comparator of a median pooled over several.
    const comparator = row.magnitudeSeries > 1 ? 'each series\' own reference' : shortSmiles(reference);
    return this.kit.summaryRow({
      depiction: this.kit.depiction(row.subst, CARD_CORE_W, CARD_CORE_H),
      name: `${win.matrix.label} · ${position}`,
      badges,
      // The strip below carries the per-slot picture where one can be drawn, and a range line beside it
      // would say the same thing twice and worse.
      desc: row.magnitude === null ? 'no reference substituent to compare against' :
        slots.length > 0 ? `vs ${comparator} over ${row.magnitudeSeries} series` :
          `vs ${comparator} · ${this.kit.formatEffect(row.lo)} → ${this.kit.formatEffect(row.hi)} ` +
          `over ${row.magnitudeSeries} series`,
      mark: row.magnitude === null || scale <= 0 ? undefined :
        this.effects.effectBar(row.magnitude, scale, data.modelError),
      value: row.magnitude === null ? '—' : this.kit.formatEffect(row.magnitude),
      caption: 'vs reference',
      valueTip: row.magnitude === null ?
        'No series contributing this placing also carries its reference substituent as a column, so ' +
        'there is no within-series comparison to quote. The row still ranks on how often it placed.' :
        'Median difference against each series\' own most-common substituent' +
        (row.magnitudeSeries > 1 ? ', which differs from series to series' : ` (${reference})`) +
        '. A within-series difference, so it is free of the fact that two series centre their effects ' +
        `over different substituent menus — but only ${row.magnitudeSeries} of the ${row.k} ` +
        'contributing series carry one, and two series\' numbers are still only loosely comparable.',
      tip: this.kit.openTip(win.matrix, win.ri, win.ci),
      onClick: () => host.revealCell(win.matrix, win.ri, win.ci, position),
    });
  }
}
