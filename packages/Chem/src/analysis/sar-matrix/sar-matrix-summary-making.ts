/* The Summary's Worth-making segment: the compounds the table holds but never assayed, and the
   analogs nobody has made, ranked by what each would gain over its own series' best. It owns the
   grid it draws and reaches the panel only for the host and the rendering helpers every segment
   shares. */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {closeGridQuietly} from './sar-matrix-types';
import {ANALOG_LIST_MAX, BEST_FIT_MIN_N, MIN_SUPPORT, SummaryData, SummaryRow, supportOf,
  TRUST_R2} from './sar-matrix-summary-data';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {ANALOG_W, CARD_CORE_H, CARD_CORE_W, CELL_H, CORE_W, count,
  MatrixCellRef} from './sar-matrix-ui-common';

const ANALOG_COLS = {
  analog: 'Analog', predicted: 'Predicted', gain: 'Gain', interest: 'Gain / error', support: 'n',
  r2: 'R²', rmse: '± error', neighbours: 'Tried around it', series: 'Series', tier: 'Tier',
  evidence: 'Evidence', core: 'Core', fixed: 'Fixed R-groups', rgroup: 'R-group',
};

/** One row of the Worth-making grid and which list it came off: `ranked` cleared every gate, `thin` rests
 *  on a fit with too few cross-validatable cells to rank by, `ungated` cleared no gate at all and is
 *  shown only when the other two lists are empty. */
interface AnalogRow {
  row: SummaryRow;
  evidence: string;
}

/** What this segment needs from the panel. The panel satisfies it as it is. */
export type MakingKit = Pick<SummaryPanel, 'host' | 'prose' | 'reason' | 'summaryRow' | 'depiction' |
  'badge' | 'trustDot' | 'cellLocation' | 'openTip' | 'revealTrust'>;

export class SummaryMaking {
  private analogGrid: DG.Grid | null = null;
  private analogSub: Subscription | null = null;
  private analogSources: MatrixCellRef[] = [];

  constructor(private readonly kit: MakingKit) {}

  close(): void {
    this.analogSub?.unsubscribe();
    this.analogSub = null;
    closeGridQuietly(this.analogGrid);
    this.analogGrid = null;
    this.analogSources = [];
  }

  /** Not a scroll pane: the ranked list is a DG.Grid, which needs a definite height to virtualize
   *  against, and a grid inside a scrolling parent resolves to its own minimum however tall the dock
   *  is. The shelf above it keeps a bounded share and scrolls inside it. */
  build(data: SummaryData): HTMLElement {
    return ui.divV([this.buildShelfSegment(data), this.buildAnalogBlock(data)], 'chem-sar-sum-making');
  }

  /** Compounds the dataset already holds with no activity value. Deliberately not gated on R²: a plate
   *  is cheap and a synthesis is not. */
  private buildShelfSegment(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const body: HTMLElement[] = [];
    body.push(ui.divText('Already in the table, never assayed — test these first',
      'chem-sar-sum-card-title'));
    body.push(this.kit.prose('No synthesis needed: these compounds are rows of your table that carry no ' +
      'value in the activity column, ranked by what the model predicts for them.',
    'They pass no fit gate, only support ≥ 2 measured compounds on both axes — running an assay on a ' +
      'compound you already have is cheap, so the bar it has to clear is lower than for a synthesis.'));
    // Filtered inside the pool's own passes, not afterwards: a pure extrapolation at support 1 would
    // otherwise top the first card a chemist sees, with nothing but a faint value to say so.
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
          this.kit.trustDot(matrix),
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
      const refs: MatrixCellRef[] = shelf.map((row) => ({matrix: row.matrix, ri: row.ri, ci: row.ci}));
      body.push(ui.divH([ui.button(refs.length === 1 ? 'Queue this compound for testing' :
        `Queue these ${refs.length} for testing`,
      () => host.addCellsToMakeList(refs, 'Nothing to add.'))], 'chem-sar-sum-foot'));
    }

    return ui.divV(body, 'chem-sar-sum-shelf');
  }

  /**
   * The analogs the model says are worth making, ranked on the gain each buys over the best
   * compound its own series has already made, in that series' own leave-one-out error.
   *
   * Not on the predicted value. Under an additive model the global maximum is the best row crossed
   * with the best column — one obvious corner per series, which a chemist reads straight off the grid;
   * the intercept it carries is the series' own observed mean, so ranking on it ranks series by what
   * was already made; and a high prediction is a reason to synthesise only if it beats the best
   * compound already measured. The division by the series' own error is a reliability weighting:
   * a +1.2 where predictions are typically off by 1.0 is a worse offer than a +0.7 where they are off
   * by 0.15.
   */
  private buildAnalogBlock(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const body: HTMLElement[] = [
      ui.divText('Worth making — ranked by gain over what this series has already made',
        'chem-sar-sum-card-title'),
      ui.divText('Every number here is a model output: this series\' fitted core effect plus its ' +
        'fitted substituent effect. Nothing here has been measured. The additive model cannot see a ' +
        'cliff, so a large gain is a hypothesis to test and a small one is not a refutation.',
      'chem-sar-cp-hint'),
    ];
    if (data.trustedNoStructure > 0) {
      body.push(ui.divText(`${count(data.trustedNoStructure)} further predictions pass the same gate but ` +
        'have no structure to make: the core carries an attachment point none of the picked fragment ' +
        'columns fills. Add that column and they appear here.', 'chem-sar-cp-hint'));
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
      // The gate withholds confidence, not the ranking. Naming the best candidates and saying what they
      // are short of leaves a next structure on the screen; refusing to name any leaves none.
      body.push(ui.divText(`Nothing clears the trust gate, so these are the largest predicted gains ` +
        `without it — ranked on gain alone, not on gain over error. ` +
        `${this.gateShortfall(data)}`, 'chem-sar-sum-warn'));
    }

    const rows: AnalogRow[] = [...ranked.map((row) => ({row, evidence: 'ranked'})),
      ...thin.map((row) => ({row, evidence: 'thin'})),
      ...spare.map((row) => ({row, evidence: 'ungated'}))];
    const structures = new Set(rows.map((r) => r.row.key)).size;
    const series = new Set(rows.map((r) => r.row.matrix.id)).size;
    const dropped = data.withheldBelowError + data.withheldThinSupport + data.withheldFitFails +
      data.withheldUnchecked + data.withheldNotConverged;
    const chips: HTMLElement[] = [
      this.kit.badge(`${count(structures)} analogs · ${count(series)} series`,
        'Molecules, not cells: one structure is proposed by every tier that folded its core, and only ' +
        'the best-evidenced occurrence is listed.'),
    ];
    if (dropped > 0) {
      chips.push(this.kit.badge(`${count(dropped)} dropped: ${count(data.withheldBelowError)} gain under ` +
        `this series' own error · ${count(data.withheldThinSupport)} thin support`,
      'A prediction whose gain is smaller than its series\' own typical prediction error is one the ' +
      'model cannot distinguish from the best compound that series has already measured. Thin support ' +
      `means fewer than ${MIN_SUPPORT} measured compounds on one of the two axes.`, true));
    }
    // Only the checked-and-failed cause has a list: a fit that was never cross-validated and a fit that
    // stopped short of its tolerance appear on no screen, so they say so rather than pointing at one.
    if (data.withheldFitFails > 0) {
      const drop = this.kit.badge(`${count(data.withheldFitFails)} in a fit checked and failed →`,
        `These series' own R² is below ${TRUST_R2}. Click for the list.`, true);
      drop.onclick = () => this.kit.revealTrust();
      chips.push(drop);
    }
    if (data.withheldUnchecked + data.withheldNotConverged > 0) {
      chips.push(this.kit.badge(`${count(data.withheldUnchecked)} fit never checked · ` +
        `${count(data.withheldNotConverged)} fit did not converge`,
      'Unchecked means the series has too few cross-validatable cells to verify — unverified, not ' +
      'wrong. A fit that did not reach its tolerance is excluded from every pooled comparison, and R² ' +
      'says nothing about it. Neither has a list to open.', true));
    }
    if (data.alreadyHeld > 0) {
      chips.push(this.kit.badge(`${count(data.alreadyHeld)} already in the table — see the block above`,
        'The structure is one the table already carries, so this is an assay to run rather than a ' +
        'synthesis.'));
    }
    body.push(ui.divH(chips, 'chem-sar-sum-chips'));

    if (data.analogOverflow > 0) {
      body.push(ui.divText(`Ranked and thin hold the top ${count(ANALOG_LIST_MAX)} each by gain in ` +
        `their own series' error: ${count(data.analogOverflow)} further structures cleared the same ` +
        'gate and are in neither. Narrow the analysis to reach them.', 'chem-sar-cp-hint'));
    }

    const refs: MatrixCellRef[] = rows.map((r) => ({matrix: r.row.matrix, ri: r.row.ri, ci: r.row.ci}));
    body.push(ui.divH([
      // "These", not "all": the list is capped, and a button promising all of them would be the cap
      // stated as a total for the third time.
      ui.button(`Add these ${refs.length} to Make list`,
        () => host.addCellsToMakeList(refs, 'Nothing to add.')),
      ui.button('Add selected', () => host.addCellsToMakeList(this.selectedAnalogs(),
        'Select rows in the list first.')),
    ], 'chem-sar-sum-foot'));
    if (thin.length > 0) {
      body.push(ui.divText(ranked.length === 0 ?
        'Thin evidence — fewer than ' + BEST_FIT_MIN_N + ' cross-validatable cells behind every one of ' +
        'these series\' fits, so none of them is ranked.' :
        `The ${count(thin.length)} rows marked thin rest on fewer than ${BEST_FIT_MIN_N} ` +
        'cross-validatable cells and are not ranked against the rows above them.', 'chem-sar-cp-hint'));
    }
    body.push(this.buildAnalogGrid(data, rows));
    return ui.divV(body, 'chem-sar-sum-making-block');
  }

  /** What the best ungated candidates are short of, in the order a reader can act on: a fit nobody
   *  checked is a different request from a fit that was checked and failed. */
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

  /** A DG.Grid, because browsing is the point: the chemist sorts by support, filters to one series,
   *  selects fifteen and sends them to the Make list. */
  private buildAnalogGrid(data: SummaryData, rows: AnalogRow[]): HTMLElement {
    const host = this.kit.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const byMatrix = new Map(host.matrices.map((m, i) => [m.id, host.matrixTiers[i]]));
    const best = new Map(data.series.map((s) => [s.matrix.id, s.best]));
    // With fragment columns every row of a series can carry the same core structure, and drawing the
    // same picture on every grid row is worse than drawing none.
    const degenerate = rows.every(({row}) =>
      new Set(row.matrix.rows.map((r) => r.keySmiles)).size === 1);

    const analog: string[] = [];
    const predicted: number[] = [];
    const gains: number[] = [];
    const interest: number[] = [];
    const support: number[] = [];
    const r2: number[] = [];
    const rmse: number[] = [];
    const neighbours: number[] = [];
    const seriesName: string[] = [];
    const tier: string[] = [];
    const evidence: string[] = [];
    const core: string[] = [];
    const rgroup: string[] = [];
    const sources: MatrixCellRef[] = [];

    for (const {row, evidence: kind} of rows) {
      const {matrix, ri, ci} = row;
      const cell = matrix.cells[ri][ci];
      // Nullable: an ungated row can come off a series whose fit was never cross-validated at all.
      const conf = matrix.confidence ?? null;
      const anchor = best.get(matrix.id) ?? null;
      const gain = anchor === null ? 0 : dir * (cell.value! - anchor.value);
      analog.push(cell.smiles ?? '');
      predicted.push(cell.value!);
      gains.push(gain);
      interest.push(conf !== null && conf.rmse > 0 ? gain / conf.rmse : NaN);
      support.push(cell.support ?? 0);
      r2.push(conf?.r2 ?? NaN);
      rmse.push(conf?.rmse ?? NaN);
      // Only for the rows that survived the ranking: this is O(rows + columns) per call.
      neighbours.push(host.observedNeighbours(matrix, ri, ci));
      seriesName.push(matrix.label);
      tier.push(`L${byMatrix.get(matrix.id) ?? 1}`);
      evidence.push(kind);
      core.push(degenerate ? Object.values(matrix.rows[ri].foldedValues).join(' · ') :
        matrix.rows[ri].keySmiles);
      rgroup.push(matrix.columns[ci].substSmiles);
      sources.push({matrix, ri, ci});
    }

    const columns = [
      DG.Column.fromStrings(ANALOG_COLS.analog, analog),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.predicted, predicted),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.gain, gains),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.interest, interest),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, ANALOG_COLS.support, support),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.r2, r2),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.rmse, rmse),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, ANALOG_COLS.neighbours, neighbours),
      DG.Column.fromStrings(ANALOG_COLS.series, seriesName),
      DG.Column.fromStrings(ANALOG_COLS.tier, tier),
      DG.Column.fromStrings(ANALOG_COLS.evidence, evidence),
      DG.Column.fromStrings(degenerate ? ANALOG_COLS.fixed : ANALOG_COLS.core, core),
      DG.Column.fromStrings(ANALOG_COLS.rgroup, rgroup),
    ];
    const descriptions: {[name: string]: string} = {
      [ANALOG_COLS.analog]: 'The assembled structure — a combination this dataset has no row for. It ' +
        'knows nothing of your collection, the catalogue or the literature, so this is "not in this ' +
        'dataset", never "never made".',
      [ANALOG_COLS.predicted]: 'A model output in the host column\'s units: this series\' fitted core ' +
        'effect plus its fitted substituent effect. Nothing here was measured.',
      [ANALOG_COLS.gain]: 'Against the best MEASURED compound of this same series — never against ' +
        'another series, and never against the dataset as a whole.',
      [ANALOG_COLS.interest]: 'The ranking quantity: the gain in multiples of this series\' own ' +
        'leave-one-out error. Under one error the model cannot tell this analog from the compound you ' +
        'already have.',
      [ANALOG_COLS.support]: 'Measured compounds on the weaker of this core and this substituent.',
      [ANALOG_COLS.r2]: 'How well this series\' model predicts a measured compound when that compound ' +
        'is left out of the fit. A property of the model, not of the data.',
      [ANALOG_COLS.rmse]: 'Typical size of this series\' prediction error.',
      [ANALOG_COLS.neighbours]: 'Measured cells sharing this core OR this substituent — how much has ' +
        'already been tried around this combination. Not the size of the series.',
      [ANALOG_COLS.evidence]: `"thin" means fewer than ${BEST_FIT_MIN_N} cross-validatable cells behind ` +
        'this series\' fit, so its error is not a scale to rank by; those rows are listed, not ranked. ' +
        '"ungated" means the row cleared no trust gate and is shown because nothing else did.',
      [degenerate ? ANALOG_COLS.fixed : ANALOG_COLS.core]: degenerate ?
        'The substituents this row holds fixed — every row of these series draws the same core.' :
        'The core this analog is built on.',
    };
    for (const column of columns) {
      const text = descriptions[column.name];
      if (text !== undefined)
        column.setTag(DG.TAGS.DESCRIPTION, text);
    }
    for (const name of [ANALOG_COLS.analog, ANALOG_COLS.rgroup, degenerate ? '' : ANALOG_COLS.core]) {
      const column = columns.find((c) => c.name === name);
      if (column !== undefined)
        column.semType = DG.SEMTYPE.MOLECULE;
    }

    const frame = DG.DataFrame.fromColumns(columns);
    this.analogSources = sources;
    const grid = DG.Viewer.grid(frame);
    this.analogGrid = grid;
    // The grid root is a ui-box, which pins itself to a fixed size and leaves the rest blank.
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    grid.setOptions({rowHeight: CELL_H});
    for (const [name, width] of [[ANALOG_COLS.analog, ANALOG_W], [ANALOG_COLS.core, CORE_W],
      [ANALOG_COLS.rgroup, CARD_CORE_W]] as [string, number][]) {
      const gridCol = grid.col(name);
      if (gridCol)
        gridCol.width = width;
    }
    // A click, not the current row: setting a current row is something the grid does to itself while it
    // is built, and landing on a cell switches the outer tab.
    this.analogSub = grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const at = cell.tableRowIndex ?? -1;
      if (at >= 0 && at < sources.length) {
        const ref = sources[at];
        host.revealCell(ref.matrix, ref.ri, ref.ci);
      }
    });
    return ui.div([grid.root], 'chem-sar-sum-analog-grid');
  }
}
