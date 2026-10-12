import * as ui from 'datagrok-api/ui';

import {BEST_FIT_MIN_N, MIN_LISTED_SUPPORT, MIN_SUPPORT, SummaryData, SummaryRow} from './sar-matrix-summary-data';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {count} from '../sar-matrix-ui-common';

/** Analog rows shown before the rest fold away; the whole ranked list goes to the Make list, which sorts it. */
const MAKING_ROWS = 8;

type Best = {value: number} | null;

export class SummaryMaking {
  constructor(private readonly kit: SummaryPanel) {}

  build(data: SummaryData): HTMLElement {
    return ui.divV([this.buildShelf(data), this.buildAnalogs(data)], 'chem-sar-sum-making');
  }

  private buildShelf(data: SummaryData): HTMLElement {
    const shelf = data.test.take();
    return this.section(`Already in the table, never assayed · ${count(shelf.length)}`,
      'Compounds in your table with no activity value, ranked by their predicted value. Testing them needs ' +
      'no synthesis. Click to open.', () => this.shelfBody(data, shelf));
  }

  private shelfBody(data: SummaryData, shelf: SummaryRow[]): HTMLElement[] {
    const host = this.kit.host;
    const body: HTMLElement[] = [
      this.kit.prose('Rows of your table with no activity value, ranked by prediction. No synthesis needed.',
        `No fit gate, only ≥ ${MIN_LISTED_SUPPORT} measured compounds on both axes: an assay is cheap.`),
    ];
    if (shelf.length === 0) {
      body.push(...this.kit.reason(host.predictUnmeasured ?
        `No untested compound rests on ${MIN_LISTED_SUPPORT} measured compounds on both axes.` :
        '"Predict untested compounds" is off, so untested compounds carry no predicted value to rank by.'));
      return body;
    }
    const bests = this.seriesBests(data);
    for (const row of shelf) {
      const cell = row.matrix.cells[row.ri][row.ci];
      body.push(...this.row(row, bests.get(row.matrix.id) ?? null, host.cellIdText(cell) ?? '(no id)',
        this.kit.cellLocation(row.matrix, row.ri, row.ci), this.kit.trustBadges(row.matrix), '~'));
    }
    body.push(ui.divH([ui.button(shelf.length === 1 ? 'Queue this compound for testing' :
      `Queue these ${shelf.length} for testing`, () => host.addCellsToMakeList(shelf, 'Nothing to add.'))],
    'chem-sar-sum-chips'));
    return body;
  }

  /** Ranked on gain, not on the predicted value, which only ranks series by what was already made. */
  private buildAnalogs(data: SummaryData): HTMLElement {
    const ranked = data.analogs.all;
    // One structure can be predicted in a thin series and in its coarser parent.
    const rankedKeys = new Set(ranked.map((row) => row.key));
    const thin = data.analogsThin.all.filter((row) => !rankedKeys.has(row.key));
    const fallback = ranked.length === 0 && thin.length === 0;
    const rows = [...ranked, ...thin, ...(fallback ? data.analogsAny.all : [])];
    const title = rows.length === 0 ? 'Worth making · none' :
      `Worth making · ${count(rows.length)} analog${rows.length === 1 ? '' : 's'}`;
    const tip = 'Predicted analogs nobody has made, ranked by how much they would beat the best compound of ' +
      'their series. Click to open.';
    return this.section(title, tip,
      () => this.analogBody(data, rows, new Set(thin.map((row) => row.key)), fallback), 'chem-sar-sum-making-block');
  }

  private analogBody(data: SummaryData, rows: SummaryRow[], thinKeys: Set<string>, fallback: boolean):
    HTMLElement[] {
    const host = this.kit.host;
    const body: HTMLElement[] = [ui.divText('Ranked by gain over the best compound each series made. Predicted, ' +
      'not measured: the additive model cannot see a cliff.', 'chem-sar-cp-hint')];
    if (data.trustedNoStructure > 0) {
      body.push(ui.divText(`${count(data.trustedNoStructure)} more predictions have no structure: add the ` +
        'R-group column that fills their open attachment point.', 'chem-sar-cp-hint'));
    }
    if (rows.length === 0) {
      body.push(...this.kit.reason(!host.predictVirtual ?
        'Predicted analogs are off — turn on "Predict virtual analogs" to see what to make.' :
        `No prediction resting on ${MIN_LISTED_SUPPORT} or more measured compounds on both axes beats its ` +
        'series\' best measured compound.'));
      return body;
    }
    if (fallback) {
      body.push(ui.divText('Nothing clears the trust gate; these are the largest gains without it. ' +
        this.gateShortfall(data), 'chem-sar-sum-prose chem-sar-chip-partial'));
    }

    const bests = this.seriesBests(data);
    const draw = (part: SummaryRow[]): HTMLElement[] => part.flatMap((row) => {
      const best = bests.get(row.matrix.id) ?? null;
      return this.row(row, best, this.kit.cellLocation(row.matrix, row.ri, row.ci), this.gainText(row, best),
        // An ungated list says so once above it, so only a thin row among ranked ones is marked.
        thinKeys.has(row.key) ? [this.kit.badge('thin', `The series has fewer than ${BEST_FIT_MIN_N} checkable ` +
          'cells, so this is listed, not ranked.', true)] : [], '');
    });
    body.push(...draw(rows.slice(0, MAKING_ROWS)));
    if (rows.length > MAKING_ROWS)
      body.push(...this.kit.lazyFold(`${count(rows.length - MAKING_ROWS)} more`, () => draw(rows.slice(MAKING_ROWS))));
    body.push(ui.divH([ui.button(rows.length === 1 ? 'Add it to Make list' : `Add all ${rows.length} to Make list`,
      () => host.addCellsToMakeList(rows, 'Nothing to add.'))], 'chem-sar-sum-chips'));
    if (data.analogOverflow > 0) {
      body.push(ui.divText(`${count(data.analogOverflow)} more qualified past these; narrow the analysis to ` +
        'reach them.', 'chem-sar-cp-hint'));
    }
    return body;
  }

  private section(title: string, tip: string, fill: () => HTMLElement[], cls = ''): HTMLElement {
    const [head, body] = this.kit.lazyFold(title, fill);
    ui.tooltip.bind(head, tip);
    return ui.divV([head, body], `chem-sar-sum-card ${cls}`.trim());
  }

  private seriesBests(data: SummaryData): Map<string, Best> {
    return new Map(data.series.map((stat) => [stat.matrix.id, stat.best]));
  }

  private row(row: SummaryRow, best: Best, name: string, desc: string, badges: HTMLElement[],
    approx: string): HTMLElement[] {
    const host = this.kit.host;
    const {matrix, ri, ci} = row;
    const cell = matrix.cells[ri][ci];
    const expand = this.kit.lazyBody(() => this.detail(row, best));
    return [this.kit.summaryRow({
      depiction: this.kit.depiction(cell.smiles),
      name,
      badges: [this.kit.badge(`n=${cell.support ?? 0}`, 'Measured compounds behind the prediction on the weaker ' +
        'of its core and its substituent.'), ...badges],
      desc,
      value: `${approx}${host.formatActivity(cell.value!)}`,
      caption: 'predicted',
      cart: {matrix, ri, ci},
      valueTip: 'Fitted core effect plus substituent effect; not measured.',
      tip: 'Click to expand; the grid icon opens it in the SAR Matrix.',
      onClick: () => host.revealCell(matrix, ri, ci),
      onChevron: expand.toggle,
    }), expand.body];
  }

  private gainText(row: SummaryRow, best: Best): string {
    if (best === null)
      return 'its series has no measured compound';
    const gain = this.kit.dir * (row.matrix.cells[row.ri][row.ci].value! - best.value);
    return `gain ${this.kit.formatEffect(gain)} over the best compound its series made`;
  }

  private detail(row: SummaryRow, best: Best): HTMLElement[] {
    const {matrix, ri, ci} = row;
    const real = (c: {kind: string} | undefined): boolean => c?.kind === 'real';
    const coreN = matrix.cells[ri].filter(real).length;
    const substN = matrix.cells.filter((cells) => real(cells[ci])).length;
    const half = (smiles: string, label: string): HTMLElement => ui.divV([
      this.kit.depiction(smiles), ui.divText(label, 'chem-sar-cp-hint')]);
    const facts = [`${count(coreN)} with this core`, `${count(substN)} with this substituent`];
    if (best !== null)
      facts.push(`series best ${this.kit.host.formatActivity(best.value)}`);
    return [ui.divH([half(matrix.rows[ri].keySmiles, 'core'), ui.divText('+', 'chem-sar-sum-arrow'),
      half(matrix.columns[ci].substSmiles, 'substituent'), ui.divText(facts.join(' · '), 'chem-sar-cp-hint')],
    'chem-sar-sum-pair chem-sar-sum-expand')];
  }

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
}
