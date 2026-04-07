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

  private createGrid(rowIdx: number): HTMLElement {
    const mol = this.moleculeColumn?.get(rowIdx);
    const score = Number(this.mpoColumn.get(rowIdx) ?? 0);

    const molDiv = renderMolecule(mol, {
      width: this.sizesMap[this.size].width,
      height: this.sizesMap[this.size].height,
    });

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
    const sortedRows = Array.from({length: this.dataFrame.rowCount}, (_, i) => i)
      .sort((a, b) =>
        this.mode === 'best' ?
          Number(this.mpoColumn.get(b)) - Number(this.mpoColumn.get(a)) :
          Number(this.mpoColumn.get(a)) - Number(this.mpoColumn.get(b)),
      )
      .slice(0, rowCount);

    const grids = sortedRows.map((i) => this.createGrid(i));

    const panel = ui.divH(grids, 'chem-viewer-grid');
    this.root.appendChild(ui.panel([ui.divV([panel])]));
  }
}
