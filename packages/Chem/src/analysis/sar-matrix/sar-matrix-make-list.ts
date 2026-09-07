import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {renderMolecule} from '../../rendering/render-molecule';
import {fitAdditiveEffects} from './sar-matrix-assemble';
import {closeGridQuietly, SarMatrix, SarMatrixCell} from './sar-matrix-types';
import {ANALOG_W, CELL_H, CELL_W, CORE_W, MatrixCellRef} from './sar-matrix-ui-common';

/** Spelled out in the Context Panel, where there is room to say what produced the number. */
const FREE_WILSON_METHOD = 'local Free-Wilson (row + column effects)';
/** The make-list column says only which of the two a row is; the panel carries the detail. */
const METHOD_PREDICTED = 'predicted';
const METHOD_MEASURED = 'measured';
/** Whether a make-list row is a compound that exists or one still to be made. `Untested` is the third
 *  case: in hand, never assayed — the row to act on by testing, not by synthesizing. */
const STATUS_MADE = 'Synthesized';
const STATUS_UNTESTED = 'Untested';
const STATUS_VIRTUAL = 'Virtual';
const MAKELIST_NAME = 'SAR make list';
/** The make-list column a row is identified by: the same analog reached twice is one entry. */
export const MAKELIST_STRUCTURE = 'Structure';
/** Fixed caption: what the value is measured on is carried per row, so the column name stays stable
 *  as rows collected under different settings accumulate. */
const MAKELIST_ACTIVITY = 'Activity';
/** Padding the context panel's own box adds around a drawn structure. */
const BOX_CHROME = 12;
const CP_STRUCT_MIN_W = 200;
const CP_STRUCT_MAX_W = 520;
const CP_STRUCT_ASPECT = 0.55;
const PANEL_CHROME = 30;

/** What the panel reads back from the viewer: how to format an activity, which cell is selected, and
 *  the column/scale that number is measured on. */
export interface MakeListHost {
  readonly activityColumnName: string;
  readonly scalingLabel: string;
  readonly selectedCell: MatrixCellRef | null;
  formatActivity(value: number): string;
  showMoleculeContext(smiles: string, build?: () => HTMLElement): void;
}

/** The Make list tab and the Context Panel content for a clicked cell. Both are the same concern:
 *  what a single analog is, why it was predicted, and collecting it for synthesis. */
export class MakeListPanel {
  readonly root = ui.divV([], 'chem-sar-makelist');
  private makeList: DG.DataFrame | null = null;
  private makeListGrid: DG.Grid | null = null;
  private readonly sources: MatrixCellRef[] = [];
  /** Current-cell subscription on the list; torn down with the grid it drives. */
  private currentCellSub: {unsubscribe(): void} | null = null;
  /** Context-panel structure size observer; replaced per click so they do not stack. */
  private cpStructureSub: {unsubscribe(): void} | null = null;
  private cpStructureToken = 0;

  constructor(private readonly host: MakeListHost) {}

  /** Drop the tab grid and the context-panel observer: the panel outlives the viewer, so its
   *  observer must not keep firing into a closed one. */
  release(): void {
    this.releaseMakeListGrid();
    this.cpStructureSub?.unsubscribe();
    this.cpStructureSub = null;
  }

  /** Every assembled virtual cell in the given matrices, as (matrix, row, column) references. */
  private analogCells(matrices: SarMatrix[]): MatrixCellRef[] {
    const out: MatrixCellRef[] = [];
    for (const matrix of matrices) {
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const cell = matrix.cells[ri][ci];
          if (cell.kind === 'virtual' && cell.smiles !== null && cell.value !== null)
            out.push({matrix, ri, ci});
        }
      }
    }
    return out;
  }
  /** A make-list table of matrix cells, observed and predicted alike. `Status` says whether the
   *  compound exists and `Method` where its number came from — the two are independent, because a
   *  compound the set holds but never tested is made AND predicted. `Support (n)` is left empty rather
   *  than zeroed for a measured value, where zero would read as no support. */
  private buildAnalogTable(cells: MatrixCellRef[]): DG.DataFrame {
    const molCol = (name: string, values: string[]): DG.Column => {
      const col = DG.Column.fromStrings(name, values);
      col.semType = DG.SEMTYPE.MOLECULE;
      return col;
    };
    const cell = (c: MatrixCellRef): SarMatrixCell => c.matrix.cells[c.ri][c.ci];
    const isVirtual = (c: MatrixCellRef): boolean => cell(c).kind === 'virtual';
    const isPredicted = (c: MatrixCellRef): boolean => cell(c).kind !== 'real';
    const support = DG.Column.fromList('int', 'Support (n)',
      cells.map((c) => isPredicted(c) ? cell(c).support ?? 0 : null));
    const df = DG.DataFrame.fromColumns([
      molCol(MAKELIST_STRUCTURE, cells.map((c) => cell(c).smiles!)),
      DG.Column.fromStrings('Status', cells.map((c) => isVirtual(c) ? STATUS_VIRTUAL :
        cell(c).kind === 'unmeasured' ? STATUS_UNTESTED : STATUS_MADE)),
      DG.Column.fromList('double', MAKELIST_ACTIVITY, cells.map((c) => cell(c).value!)),
      // Per row, not in the caption: the list outlives the settings it was collected under, so rows
      // gathered before a scaling change would otherwise sit unmarked under the current name.
      DG.Column.fromStrings('Activity basis',
        cells.map(() => `${this.host.activityColumnName || 'activity'} (${this.host.scalingLabel})`)),
      support,
      DG.Column.fromStrings('Series', cells.map((c) => c.matrix.label)),
      molCol('Core', cells.map((c) => c.matrix.rows[c.ri].keySmiles)),
      molCol('Substituent', cells.map((c) => c.matrix.columns[c.ci].substSmiles)),
      DG.Column.fromStrings('Method', cells.map((c) => isPredicted(c) ? METHOD_PREDICTED : METHOD_MEASURED)),
    ]);
    df.name = MAKELIST_NAME;
    return df;
  }
  /** Append analogs, creating the list on first use. A structure already present is skipped — the
   *  same analog is reachable from several matrices — so the count of what actually landed is
   *  returned alongside the total. Deduplication runs over the incoming batch before the table is
   *  built, which keeps `sources` row-for-row with the list. */
  private collectAnalogs(cells: MatrixCellRef[]): {added: number, total: number} {
    const seen = new Set<string>(this.makeList?.col(MAKELIST_STRUCTURE)?.toList() ?? []);
    const fresh = cells.filter((c) => {
      const smiles = c.matrix.cells[c.ri][c.ci].smiles!;
      if (seen.has(smiles))
        return false;
      seen.add(smiles);
      return true;
    });
    if (fresh.length > 0) {
      const table = this.buildAnalogTable(fresh);
      if (this.makeList === null)
        this.makeList = table;
      else {
        const cols = table.columns.toList();
        for (let i = 0; i < table.rowCount; i++)
          this.makeList.rows.addNew(cols.map((c) => c.get(i)));
      }
      this.sources.push(...fresh);
    }
    this.renderMakeList();
    return {added: fresh.length, total: this.makeList?.rowCount ?? 0};
  }
  releaseMakeListGrid(): void {
    this.currentCellSub?.unsubscribe();
    this.currentCellSub = null;
    closeGridQuietly(this.makeListGrid);
    this.makeListGrid = null;
  }

  /** Send the clicked row's compound to the Context Panel, with its SAR context gated in — a
   *  make-list row otherwise loses every trace of the matrix cell that produced it. Clicking a Core
   *  or Substituent cell asks about that fragment, so that is what opens instead. */
  private showCurrentCompound(): void {
    const df = this.makeList;
    if (df === null)
      return;
    const cell = df.currentCell;
    const rowIdx = cell?.rowIndex ?? df.currentRowIdx;
    if (rowIdx < 0 || rowIdx >= this.sources.length)
      return;
    const column = cell?.column ?? null;
    if (column !== null && column.semType === DG.SEMTYPE.MOLECULE && column.name !== MAKELIST_STRUCTURE) {
      const fragment = column.get(rowIdx);
      if (fragment)
        this.host.showMoleculeContext(fragment);
      return;
    }
    const smiles = df.col(MAKELIST_STRUCTURE)!.get(rowIdx);
    const source = this.sources[rowIdx];
    if (smiles)
      this.host.showMoleculeContext(smiles, () => this.buildCellPanel(source.matrix, source.ri, source.ci));
  }
  /** The make-list tab: the collected analogs as a grid, or a note when nothing has been collected. */
  renderMakeList(): void {
    this.releaseMakeListGrid();
    ui.empty(this.root);
    const rows = this.makeList?.rowCount ?? 0;
    if (this.makeList === null || rows === 0) {
      this.root.appendChild(ui.divText('Nothing collected yet. Select any cell in the matrix — ' +
        'predicted or already made — and add it from the matrix header, or right-click ' +
        'the matrix to add a whole series of predictions.',
      'chem-sar-empty-note'));
      return;
    }
    // Handed out as a copy: the tab keeps collecting, and edits to the opened table can't corrupt it.
    const open = ui.button('Add to workspace', () => grok.shell.addTableView(this.makeList!.clone()));
    ui.tooltip.bind(open, 'Add a copy to the workspace as a table, where it can be saved, exported or joined');
    // Two buttons rather than one that reads the selection: clicking a row is how the list drives the
    // Context Panel, so a row is current from the first click onwards and a single button could never
    // be asked to clear everything again.
    const remove = ui.button('Remove', () => {
      const picked = this.pickedRows();
      if (picked.length === 0) {
        grok.shell.info('Click a row in the list first, or use Clear to empty it.');
        return;
      }
      // Highest index first, so the rows still to be dropped keep the positions they were found at.
      for (const i of [...picked].sort((a, b) => b - a)) {
        this.makeList!.rows.removeAt(i);
        this.sources.splice(i, 1);
      }
      if (this.makeList!.rowCount === 0) {
        this.makeList = null;
        this.sources.length = 0;
      }
      this.renderMakeList();
    });
    ui.tooltip.bind(remove, () => {
      const picked = this.pickedRows();
      return picked.length === 0 ? 'Remove — click a row in the list first' :
        picked.length === 1 ? 'Remove the selected compound from the list' :
          `Remove the ${picked.length} selected compounds from the list`;
    });
    const clear = ui.button('Clear', () => {
      this.makeList = null;
      this.sources.length = 0;
      this.renderMakeList();
    });
    ui.tooltip.bind(clear, 'Empty the whole list');
    // Kept away from the collecting actions: it discards the whole list, one slip from undoing them.
    clear.classList.add('chem-sar-bar-end');
    this.root.appendChild(ui.divH([
      ui.divText(`${rows} compound${rows === 1 ? '' : 's'}`, 'chem-sar-main-title'),
      ui.divH([open, remove]),
      clear,
    ], 'chem-sar-main-bar'));
    this.makeListGrid = DG.Viewer.grid(this.makeList);
    // The grid root is a ui-box, which pins itself to a fixed size and leaves the rest of the tab blank.
    this.makeListGrid.root.style.width = '100%';
    this.makeListGrid.root.style.height = '100%';
    // Sized like the matrix cells these came from: the text-oriented defaults leave a molecule column
    // too short and too narrow to read.
    this.makeListGrid.setOptions({rowHeight: CELL_H});
    for (const [name, width] of [[MAKELIST_STRUCTURE, ANALOG_W], ['Core', CORE_W], ['Substituent', CELL_W]] as
      [string, number][]) {
      const gridCol = this.makeListGrid.col(name);
      if (gridCol)
        gridCol.width = width;
    }
    this.currentCellSub = this.makeList.onCurrentCellChanged.subscribe(() => this.showCurrentCompound());
    this.root.appendChild(ui.div([this.makeListGrid.root], 'chem-sar-makelist-grid'));
  }
  /** Bulk: every assembled virtual analog in the given matrices into the make-list. */
  addMatrixAnalogsToMakeList(matrices: SarMatrix[]): void {
    this.addCellsToMakeList(this.analogCells(matrices), 'No assembled virtual analogs to add.');
  }
  /** Bulk: an explicit set of cells — a transfer's predicted analogs, say. Cells with no structure or
   *  no value are dropped rather than added as blanks. */
  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void {
    const usable = cells.filter((c) => {
      const cell = c.matrix.cells[c.ri][c.ci];
      return cell.kind !== 'empty' && cell.smiles !== null && cell.value !== null;
    });
    if (!usable.length) {
      grok.shell.info(emptyMessage);
      return;
    }
    const {added, total} = this.collectAnalogs(usable);
    const skipped = usable.length - added;
    grok.shell.info(`Added ${added} to the Make list` +
      (skipped > 0 ? `, ${skipped} already there` : '') + ` (${total} total).`);
  }
  /** Collect whatever cell is selected in the matrix, so adding one does not depend on the Context
   *  Panel being open at the time. */
  addSelectedToMakeList(): void {
    const cell = this.host.selectedCell;
    if (cell === null) {
      grok.shell.info('Select a cell in the matrix first.');
      return;
    }
    this.addAnalogToMakeList(cell.matrix, cell.ri, cell.ci);
  }
  /** Per-cell action: append one compound to the make-list, predicted or already made. */
  addAnalogToMakeList(matrix: SarMatrix, ri: number, ci: number): void {
    const cell = matrix.cells[ri][ci];
    if (cell.kind === 'empty' || cell.smiles === null || cell.value === null) {
      grok.shell.info('This cell has no structure to add.');
      return;
    }
    const kind = cell.kind === 'virtual' ? 'Virtual analog' :
      cell.kind === 'unmeasured' ? 'Untested compound' : 'Synthesized compound';
    const {added, total} = this.collectAnalogs([{matrix, ri, ci}]);
    grok.shell.info(added > 0 ? `${kind} added to the Make list (${total} total).` :
      `${kind} is already on the Make list (${total} total).`);
  }
  /** Rows the Clear button acts on; empty discards the whole list. A cell click only moves the current
   *  row, so that counts as picking it — otherwise clicking a compound then Clear would wipe everything. */
  private pickedRows(): number[] {
    const df = this.makeList;
    if (df === null)
      return [];
    const selected = df.selection;
    if (selected.trueCount > 0) {
      const out: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if (selected.get(i))
          out.push(i);
      }
      return out;
    }
    return df.currentRowIdx >= 0 && df.currentRowIdx < df.rowCount ? [df.currentRowIdx] : [];
  }

  /** A context-panel section header with an optional provenance badge. */
  private cpSection(title: string, badge?: string): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-cp-section-title')];
    if (badge)
      parts.push(ui.divText(badge, 'chem-sar-cp-badge'));
    return ui.divH(parts, 'chem-sar-cp-section');
  }
  /** A label / value row in the context panel. */
  private cpRow(label: string, value: HTMLElement | string): HTMLElement {
    const v = typeof value === 'string' ? ui.divText(value, 'chem-sar-cp-rv') : value;
    return ui.divH([ui.divText(label, 'chem-sar-cp-rl'), v], 'chem-sar-cp-row');
  }
  /** The cell's structure, drawn to the panel's actual width and redrawn on resize past a few pixels
   *  (each redraw is an RDKit render). */
  private cpStructure(smiles: string): HTMLElement {
    const box = ui.div([], 'chem-sar-cp-structbox');
    let drawnAt = 0;
    const draw = (avail: number): void => {
      const width = Math.min(CP_STRUCT_MAX_W, Math.max(CP_STRUCT_MIN_W, avail));
      if (Math.abs(width - drawnAt) < 8)
        return;
      drawnAt = width;
      ui.empty(box);
      box.appendChild(renderMolecule(smiles,
        {width, height: Math.round(width * CP_STRUCT_ASPECT), popupMenu: false}));
    };
    // A newer panel may claim the single observer slot while this one is still awaiting layout.
    const token = ++this.cpStructureToken;
    this.cpStructureSub?.unsubscribe();
    // Wired up only once the box is in the DOM, else clientWidth is 0 and the first draw is too small.
    ui.tools.waitForElementInDom(box).then(() => {
      if (token !== this.cpStructureToken)
        return;
      // Measure the panel container (what a drag resizes) and the box, smaller wins, so the structure
      // grows with the panel but never gets stuck at a stale box width.
      const panel = box.closest('.panel-content') as HTMLElement | null;
      const fit = (): void => draw(panel === null ? Math.floor(box.clientWidth) - BOX_CHROME :
        Math.min(Math.floor(box.clientWidth) - BOX_CHROME, Math.floor(panel.clientWidth) - PANEL_CHROME));
      this.cpStructureSub?.unsubscribe();
      this.cpStructureSub = DG.debounce(ui.onSizeChanged(panel ?? box), 50).subscribe(fit);
      fit();
    });
    return box;
  }
  /** Where a predicted value came from: the additive model `rowMean + colMean - grandMean` with the
   *  contributing cells and the arithmetic spelled out. */
  private cpPrediction(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    // The effects the prediction was built from, not the row and column means: with holes the margins
    // absorb each other, so deriving them here would print a sum contradicting the number above.
    const {grandMean, rowEffect, colEffect, rowN, colN, observed} =
      fitAdditiveEffects(matrix.cells, matrix.rows.length, matrix.columns.length);
    const coreEffect = rowEffect[rowIdx];
    const substEffect = colEffect[colIdx];
    const signed = (v: number): string => `${v < 0 ? '−' : '+'}${this.host.formatActivity(Math.abs(v))}`;

    const block = ui.divV([]);
    block.appendChild(this.cpRow(`Matrix mean (n = ${observed})`, this.host.formatActivity(grandMean)));
    block.appendChild(this.cpRow(`Core effect (n = ${rowN[rowIdx]})`, signed(coreEffect)));
    block.appendChild(this.cpRow(`Substituent effect (n = ${colN[colIdx]})`, signed(substEffect)));
    block.appendChild(this.cpRow('Sum', `${this.host.formatActivity(grandMean)} ` +
      `${signed(coreEffect)} ${signed(substEffect)} = ` +
      `${this.host.formatActivity(grandMean + coreEffect + substEffect)}`));

    block.appendChild(this.cpReferences(matrix, rowIdx, colIdx));
    return block;
  }
  /** Observed compounds sharing this cell's core or substituent, the cell itself excluded. The matrix
   *  filter does not apply: it narrows what you browse, not what a compound's SAR is. */
  private cpReferences(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    const block = ui.divV([]);
    const self = matrix.cells[rowIdx][colIdx];
    const collect = (pick: (ri: number, ci: number) => boolean): {cell: SarMatrixCell, ri: number, ci: number}[] => {
      const found: {cell: SarMatrixCell, ri: number, ci: number}[] = [];
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const c = matrix.cells[ri][ci];
          if (c !== self && c.kind === 'real' && c.value !== null && pick(ri, ci))
            found.push({cell: c, ri, ci});
        }
      }
      return found;
    };
    const section = (title: string, all: {cell: SarMatrixCell, ri: number, ci: number}[]): void => {
      if (all.length === 0)
        return;
      block.appendChild(this.cpSection(title, `${all.length}`));
      block.appendChild(ui.divH(all.map((e) => this.cpFragment(e.cell.smiles, e.cell.value ?? 0)),
        'chem-sar-cp-decomp'));
    };
    section('Measured with this core', collect((ri) => ri === rowIdx));
    section('Measured with this substituent', collect((_ri, ci) => ci === colIdx));
    return block;
  }
  /** A framed fragment tile: the structure with the compound's value beneath; structure omitted for
   *  an empty SMILES. */
  private cpFragment(smiles: string | null, value?: number): HTMLElement {
    const parts: HTMLElement[] = [];
    if (smiles)
      parts.push(ui.div([renderMolecule(smiles, {width: 78, height: 52, popupMenu: false})], 'chem-sar-cp-frag-box'));
    if (value !== undefined)
      parts.push(ui.divText(this.host.formatActivity(value), 'chem-sar-cp-rv'));
    return ui.divV(parts, 'chem-sar-cp-frag');
  }
  /** SAR context for the selected cell: header, structure, potency, the prediction/references, and a
   *  make-list action. */
  buildCellPanel(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    const cell = matrix.cells[rowIdx][colIdx];
    const row = matrix.rows[rowIdx];
    const column = matrix.columns[colIdx];
    const panel = ui.divV([], 'chem-sar-matrix-cp');

    if (cell.kind === 'empty' || cell.value === null) {
      panel.appendChild(ui.divText('No compound and no prediction for this combination.'));
      return panel;
    }
    const isVirtual = cell.kind === 'virtual';
    // A compound the set holds but never tested carries a predicted number on a real structure, so it
    // reads as a compound everywhere the structure matters and as a prediction everywhere the number
    // does.
    const isUntested = cell.kind === 'unmeasured';
    const isPredicted = isVirtual || isUntested;

    const header = ui.divH([ui.h2(isVirtual ? 'Virtual analog' :
      isUntested ? 'Compound · never tested' : 'Compound')], 'chem-sar-cp-head');
    panel.appendChild(header);

    if (cell.smiles)
      panel.appendChild(this.cpStructure(cell.smiles));

    panel.appendChild(this.cpSection(isPredicted ? 'Predicted potency' : 'Potency'));
    panel.appendChild(this.cpRow(isPredicted ? 'Predicted' : 'Observed',
      ui.divText(this.host.formatActivity(cell.value), 'chem-sar-cp-value')));
    if (isPredicted) {
      panel.appendChild(this.cpRow('Method', FREE_WILSON_METHOD));
      // Call out a single-compound arm: that's where the estimate is really an extrapolation.
      const refs = cell.references ?? 0;
      const weakest = cell.support ?? 0;
      panel.appendChild(this.cpRow('Reference points',
        `n = ${refs}${weakest <= 1 ? ' · one arm has a single compound' : ''}`));
      panel.appendChild(this.cpPrediction(matrix, rowIdx, colIdx));
    } else
      panel.appendChild(this.cpReferences(matrix, rowIdx, colIdx));

    panel.appendChild(this.cpSection('Decomposition', 'R-group'));
    const parts = [this.cpFragment(row.keySmiles)];
    matrix.positions.forEach((position) => {
      const v = position === column.position ? column.substSmiles : matrix.refValues[position];
      if (v)
        parts.push(this.cpFragment(v));
    });
    panel.appendChild(ui.divH(parts, 'chem-sar-cp-decomp'));

    if (cell.smiles) {
      panel.appendChild(this.cpSection('Design action'));
      panel.appendChild(ui.divText(isVirtual ?
        'This core × substituent is not in the dataset. Add it to the Make list for synthesis triage.' :
        isUntested ?
          'You already have this compound — it just has no activity value, so the number above is ' +
          'predicted. Add it to the list to test it rather than synthesize anything.' :
          'This compound is already in the dataset. Add it to the Make list to carry it alongside the ' +
          'predictions as a measured reference — for re-testing, or to ship the whole set as one batch.',
      'chem-sar-cp-hint'));
      // One label for one action, named after the tab it fills. What the compound is for — synthesis
      // or a test — is what the hint above says, not the button.
      const btn = ui.button('Add to Make list', () => this.addAnalogToMakeList(matrix, rowIdx, colIdx));
      btn.classList.add('chem-sar-cp-generate');
      panel.appendChild(btn);
    }
    return panel;
  }
}
