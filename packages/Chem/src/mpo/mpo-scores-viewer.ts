import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';

import {ChemSearchBaseViewer} from '../analysis/chem-search-base-viewer';
import {renderMolecule} from '../rendering/render-molecule';

export class MpoScoreViewer extends ChemSearchBaseViewer {
  mpoColumn: DG.Column;
  curIdx: number = 0;
  mode: 'best' | 'worst';
  interacting = false;

  // A molecule's structure doesn't change when the profile is edited — only its score and
  // whether it ranks top/bottom. So the rendered molecule DOM is reused across renders by row
  // index; only the cheap wrapper (score text + highlight) is rebuilt.
  private molDivCache = new DG.LruCache<number, HTMLElement>(64);
  private molSize?: {width: number; height: number};

  constructor(df: DG.DataFrame, mpoColName: string, mode: 'best' | 'worst' = 'best') {
    const moleculeCol = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
    super('MPO Scores', moleculeCol!);

    const col = df.col(mpoColName);
    if (!col)
      throw new Error(`Column ${mpoColName} not found.`);
    this.mpoColumn = col;
    this.mode = mode;
    this.limit = 4;
    this.root.classList.add('chem-mpo-scores-viewer');
  }

  // Rendering is driven explicitly by the context panel's render() calls, so we skip the
  // base viewer's dataframe-event subscriptions — onMetadataChanged in particular fires on
  // every computeMpo tag write and would re-rank all rows redundantly.
  async onTableAttached(): Promise<void> {
    this.init();
    if (this.dataFrame) {
      this.moleculeColumn = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      this.moleculeColumnName = this.moleculeColumn?.name ?? '';
    }
    await this.render(true);
  }

  private showMessage(msg: string) {
    ui.empty(this.root);
    this.root.append(ui.divText(msg));
  }

  private handleClick(event: MouseEvent, rowIdx: number) {
    if (!this.dataFrame)
      return;

    this.interacting = true;
    try {
      if (event.shiftKey || event.altKey)
        this.dataFrame.selection.set(rowIdx, true);
      else if (event.metaKey)
        this.dataFrame.selection.set(rowIdx, !this.dataFrame.selection.get(rowIdx));
      else {
        this.dataFrame.currentRowIdx = rowIdx;
        this.gridSelect = true;
      }
      this.renderInternal();
    } finally {
      this.interacting = false;
    }
  }

  private getMolDiv(rowIdx: number): HTMLElement {
    this.molSize ??= {width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height};
    return this.molDivCache.getOrCreate(rowIdx, () => renderMolecule(this.moleculeColumn?.get(rowIdx), this.molSize!));
  }

  private createGrid(rowIdx: number): HTMLElement {
    const score = Number(this.mpoColumn.get(rowIdx) ?? 0);
    const molDiv = this.getMolDiv(rowIdx);

    const scoreDiv = ui.divText(score.toFixed(2));
    scoreDiv.style.marginTop = '6px';
    const propsDiv = ui.divV([scoreDiv], {style: {alignItems: 'center'}});

    const grid = ui.div([molDiv, propsDiv], {style: {position: 'relative'}});

    let divClass = 'd4-flex-col';
    if (rowIdx === this.curIdx) {
      grid.style.backgroundColor = '#ddffd9';
      divClass += ' d4-current';
    }
    if (this.dataFrame.selection.get(rowIdx)) {
      grid.style.backgroundColor = divClass === 'd4-flex-col d4-selected' ? '#f8f8df' : '#d3f8bd';
      divClass += ' d4-selected';
    }

    $(grid).addClass(divClass);
    grid.addEventListener('click', (event) => this.handleClick(event, rowIdx));

    return grid;
  }

  async renderInternal(): Promise<void> {
    if (!this.moleculeColumn || !this.mpoColumn) {
      this.showMessage('Molecule or MPO column not found');
      return;
    }

    ui.empty(this.root);

    this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;

    const rowCount = Math.min(this.limit, this.dataFrame.rowCount);
    const scores = this.mpoColumn.getRawData();
    const best = this.mode === 'best';
    // a ranks at least as high as b (higher score for 'best', lower for 'worst')
    const ranksFirst = (a: number, b: number): boolean => best ? scores[a] >= scores[b] : scores[a] <= scores[b];
    // O(n) partial selection of the top `rowCount` indices — avoids sorting all rows.
    const sortedRows: number[] = [];
    for (let i = 0; i < scores.length; i++) {
      if (sortedRows.length === rowCount && ranksFirst(sortedRows[rowCount - 1], i))
        continue;
      let pos = 0;
      while (pos < sortedRows.length && ranksFirst(sortedRows[pos], i))
        pos++;
      sortedRows.splice(pos, 0, i);
      if (sortedRows.length > rowCount)
        sortedRows.pop();
    }

    const grids = sortedRows.map((i) => this.createGrid(i));
    const panel = ui.divH(grids, 'chem-viewer-grid');
    this.root.appendChild(ui.panel([ui.divV([panel])]));
  }
}
