/* The Summary's Worth-making segment: the compounds the table holds but never assayed, and the
   analogs nobody has made, ranked by what each would gain over its own series' best. It owns the
   grid it draws. */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {closeGridQuietly} from '../sar-matrix-types';
import {ANALOG_LIST_MAX, BEST_FIT_MIN_N, MIN_SUPPORT, SummaryData, SummaryRow, supportOf,
  TRUST_R2} from './sar-matrix-summary-data';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {ANALOG_W, CARD_CORE_H, CARD_CORE_W, count, MatrixCellRef, moleculeGrid} from '../sar-matrix-ui-common';

/** One row of the Worth-making grid and which list it came off. */
interface AnalogRow {
  row: SummaryRow;
  evidence: 'ranked' | 'thin' | 'ungated';
}

export class SummaryMaking {
  private analogGrid: DG.Grid | null = null;
  private analogSub: Subscription | null = null;
  private analogSources: MatrixCellRef[] = [];

  constructor(private readonly kit: SummaryPanel) {}

  close(): void {
    this.analogSub?.unsubscribe();
    this.analogSub = null;
    closeGridQuietly(this.analogGrid);
    this.analogGrid = null;
    this.analogSources = [];
  }

  /** Not a scroll pane: the grid needs a definite height to virtualize against. */
  build(data: SummaryData): HTMLElement {
    return ui.divV([this.buildShelfSegment(data), this.buildAnalogBlock(data)], 'chem-sar-sum-making');
  }

  /** Compounds the table holds with no activity value. Not gated on R²: an assay is cheap. */
  private buildShelfSegment(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const body: HTMLElement[] = [];
    body.push(ui.divText('Already in the table, never assayed — test these first',
      'chem-sar-cp-section-title'));
    body.push(this.kit.prose('Rows of your table with no activity value, ranked by prediction. No synthesis needed.',
      'No fit gate, only ≥ 2 measured compounds on both axes: an assay is cheap.'));
    // Filtered inside the pool, so a support-1 extrapolation never tops the list.
    const shelf = data.test.take((row) => supportOf(row) < 2);
    if (shelf.length === 0) {
      body.push(...this.kit.reason(host.predictUnmeasured ?
        'No untested compound rests on two measured compounds on both axes.' :
        '"Predict untested compounds" is off, so untested compounds carry no predicted value to rank by.'));
    }
    for (const row of shelf) {
      const {matrix, ri, ci} = row;
      const cell = matrix.cells[ri][ci];
      const support = cell.support ?? 0;
      body.push(this.kit.summaryRow({
        depiction: this.kit.depiction(cell.smiles, CARD_CORE_W, CARD_CORE_H),
        name: host.cellIdText(cell) ?? '(no id)',
        badges: [
          this.kit.badge(`n=${support}`, 'Measured compounds backing the prediction on the weaker of this ' +
            'core and this substituent.'),
          this.kit.trustBadge(matrix),
        ],
        desc: this.kit.cellLocation(matrix, ri, ci),
        value: `~${host.formatActivity(cell.value!)}`,
        caption: 'predicted',
        cart: {matrix, ri, ci},
        valueTip: 'Predicted, not measured — the compound exists, only the number is estimated. This ' +
          'segment is deliberately not gated on R²: a plate is cheap and a synthesis is not.',
        tip: this.kit.openTip(matrix, ri, ci),
        onClick: () => host.revealCell(matrix, ri, ci),
      }));
    }
    if (shelf.length > 0) {
      body.push(ui.divH([ui.button(shelf.length === 1 ? 'Queue this compound for testing' :
        `Queue these ${shelf.length} for testing`,
      () => host.addCellsToMakeList(shelf, 'Nothing to add.'))], 'chem-sar-sum-chips'));
    }

    return ui.divV(body, 'chem-sar-sum-shelf');
  }

  /** Analogs ranked on gain over their series' best measured compound, in that series' own error. Not
   *  on the predicted value, which only ranks series by what was already made. */
  private buildAnalogBlock(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const body: HTMLElement[] = [
      ui.divText('Worth making — ranked by gain over what this series has already made',
        'chem-sar-cp-section-title'),
      ui.divText('Predicted, not measured. The additive model cannot see a cliff, so a gain is a hypothesis.',
        'chem-sar-cp-hint'),
    ];
    if (data.trustedNoStructure > 0) {
      body.push(ui.divText(`${count(data.trustedNoStructure)} more predictions have no structure: add the ` +
        'R-group column that fills their open attachment point.', 'chem-sar-cp-hint'));
    }

    const ranked = data.analogs.all;
    const thin = data.analogsThin.all;
    const fallback = ranked.length === 0 && thin.length === 0;
    const spare = fallback ? data.analogsAny.all : [];
    if (fallback && spare.length === 0) {
      body.push(...this.kit.reason(!host.predictVirtual ?
        'Predicted analogs are off — turn on "Predict virtual analogs" to see what to make.' :
        'No series predicts anything above the best compound it has already measured.'));
      return ui.divV(body, 'chem-sar-sum-making-block');
    }
    if (fallback) {
      // The gate withholds confidence, not a next structure.
      body.push(ui.divText('Nothing clears the trust gate; these are the largest gains without it. ' +
        this.gateShortfall(data), 'chem-sar-sum-prose chem-sar-chip-partial'));
    }

    const rows: AnalogRow[] = [...ranked.map((row): AnalogRow => ({row, evidence: 'ranked'})),
      ...thin.map((row): AnalogRow => ({row, evidence: 'thin'})),
      ...spare.map((row): AnalogRow => ({row, evidence: 'ungated'}))];
    const structures = new Set(rows.map((r) => r.row.key)).size;
    const series = new Set(rows.map((r) => r.row.matrix.id)).size;
    const dropped = data.withheldBelowError + data.withheldThinSupport + data.withheldFitFails +
      data.withheldUnchecked + data.withheldNotConverged;
    const chips: HTMLElement[] = [
      this.kit.badge(`${count(structures)} analogs · ${count(series)} series`,
        'Molecules, not cells: each structure is listed once.'),
    ];
    if (dropped > 0) {
      chips.push(this.kit.badge(`${count(dropped)} dropped: ${count(data.withheldBelowError)} gain under ` +
        `this series' own error · ${count(data.withheldThinSupport)} thin support`,
      'A gain under the series\' own error cannot be told from what was made; thin support means fewer ' +
      `than ${MIN_SUPPORT} measured compounds on an axis.`, true));
    }
    if (data.withheldFitFails > 0) {
      chips.push(this.kit.badge(`${count(data.withheldFitFails)} in a fit checked and failed`,
        `These series' own R² is below ${TRUST_R2}.`, true));
    }
    if (data.withheldUnchecked + data.withheldNotConverged > 0) {
      chips.push(this.kit.badge(`${count(data.withheldUnchecked)} fit never checked · ` +
        `${count(data.withheldNotConverged)} fit did not converge`,
      'Unchecked: too few cells to verify, not wrong. Not converged: the fit stopped short of its ' +
      'tolerance.', true));
    }
    if (data.alreadyHeld > 0) {
      chips.push(this.kit.badge(`${count(data.alreadyHeld)} already in the table — see the block above`,
        'The structure is one the table already carries, so this is an assay to run rather than a ' +
        'synthesis.'));
    }
    body.push(ui.divH(chips, 'chem-sar-sum-chips'));

    if (data.analogOverflow > 0) {
      body.push(ui.divText(`Top ${count(ANALOG_LIST_MAX)} per list; ${count(data.analogOverflow)} more ` +
        'qualified. Narrow the analysis to reach them.', 'chem-sar-cp-hint'));
    }

    body.push(ui.divH([
      // "These", not "all": the list is capped.
      ui.button(`Add these ${rows.length} to Make list`,
        () => host.addCellsToMakeList(rows.map((r) => r.row), 'Nothing to add.')),
      ui.button('Add selected', () => host.addCellsToMakeList(this.selectedAnalogs(),
        'Select rows in the list first.')),
    ], 'chem-sar-sum-chips'));
    if (thin.length > 0) {
      body.push(ui.divText(ranked.length === 0 ?
        'Thin evidence: every fit rests on fewer than ' + BEST_FIT_MIN_N + ' checkable cells; none is ranked.' :
        `${count(thin.length)} rows marked thin rest on fewer than ${BEST_FIT_MIN_N} checkable cells.`,
      'chem-sar-cp-hint'));
    }
    body.push(this.buildAnalogGrid(data, rows));
    return ui.divV(body, 'chem-sar-sum-making-block');
  }

  /** What the ungated candidates are short of; each cause calls for a different fix. */
  private gateShortfall(data: SummaryData): string {
    if (this.kit.host.matrices.every((matrix) => !matrix.confidence))
      return 'No series has a cross-validated fit, so none of them could be checked.';
    const parts: string[] = [];
    if (data.withheldThinSupport > 0) {
      parts.push(`${count(data.withheldThinSupport)} rest on fewer than ${MIN_SUPPORT} measured ` +
        'compounds on one axis');
    }
    if (data.withheldFitFails > 0)
      parts.push(`${count(data.withheldFitFails)} sit in a series whose fit was checked and failed`);
    if (data.withheldUnchecked > 0)
      parts.push(`${count(data.withheldUnchecked)} sit in a series with too few cells to check`);
    return parts.length === 0 ? '' : `Of what was withheld, ${parts.join('; ')}.`;
  }

  private selectedAnalogs(): MatrixCellRef[] {
    const grid = this.analogGrid;
    if (grid === null)
      return [];
    const selection = grid.dataFrame.selection;
    return this.analogSources.filter((_ref, i) => selection.get(i));
  }

  /** A grid, so the list can be sorted, filtered and sent to the Make list. */
  private buildAnalogGrid(data: SummaryData, rows: AnalogRow[]): HTMLElement {
    const host = this.kit.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const best = new Map(data.series.map((s) => [s.matrix.id, s.best]));
    const items = rows.map(({row, evidence}) => {
      const {matrix, ri, ci} = row;
      const cell = matrix.cells[ri][ci];
      const anchor = best.get(matrix.id) ?? null;
      // Nullable: an ungated row can come off a series whose fit was never cross-validated at all.
      return {matrix, ri, ci, cell, evidence, conf: matrix.confidence ?? null,
        gain: anchor === null ? 0 : dir * (cell.value! - anchor.value)};
    });
    type Item = typeof items[number];
    const {FLOAT, INT, STRING} = DG.COLUMN_TYPE;
    const spec: [string, DG.ColumnType, (item: Item) => string | number, string?][] = [
      ['Analog', STRING, (it) => it.cell.smiles ?? '', 'Not in this dataset, which is not the same as never made.'],
      ['Predicted', FLOAT, (it) => it.cell.value!, 'Fitted core effect plus substituent effect; not measured.'],
      ['Gain', FLOAT, (it) => it.gain, 'Over the best measured compound of the same series.'],
      ['Gain / error', FLOAT, (it) => it.conf !== null && it.conf.rmse > 0 ? it.gain / it.conf.rmse : NaN,
        'The ranking: gain in multiples of the series\' own error. Under 1 it is within noise.'],
      ['n', INT, (it) => it.cell.support ?? 0,
        'Measured compounds on the weaker of this core and this substituent.'],
      ['Series', STRING, (it) => it.matrix.label],
      ['Evidence', STRING, (it) => it.evidence, `"thin": fewer than ${BEST_FIT_MIN_N} checkable cells, so ` +
        'listed, not ranked. "ungated": cleared no gate, shown because nothing else did.'],
    ];
    const frame = DG.DataFrame.fromColumns(spec.map(([name, type, value, description]) => {
      const column = DG.Column.fromList(type, name, items.map(value));
      if (description !== undefined)
        column.setTag(DG.TAGS.DESCRIPTION, description);
      return column;
    }));
    frame.col('Analog')!.semType = DG.SEMTYPE.MOLECULE;

    this.analogSources = items;
    const grid = moleculeGrid(frame, [['Analog', ANALOG_W]]);
    this.analogGrid = grid;
    // A click, not the current row: the grid sets its current row itself while it is built.
    this.analogSub = grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const at = cell.tableRowIndex ?? -1;
      if (at >= 0 && at < items.length)
        host.revealCell(items[at].matrix, items[at].ri, items[at].ci);
    });
    return ui.div([grid.root], 'chem-sar-sum-analog-grid');
  }
}
