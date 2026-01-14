import {ChemSearchBaseViewer} from '../analysis/chem-search-base-viewer';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import {renderMolecule} from '../rendering/render-molecule';

export class MpoScoreViewer extends ChemSearchBaseViewer {
  mpoColumn: DG.Column;
  curIdx: number = 0;
  mode: 'best' | 'worst';

  constructor(df: DG.DataFrame, mpoColName: string, mode: 'best' | 'worst' = 'best') {
    super('MPO Scores', df.getCol('molecule'));
    const col = df.col(mpoColName);
    if (!col) throw new Error(`Column ${mpoColName} not found`);
    this.mpoColumn = col;
    this.mode = mode;
    this.limit = 5;

    this.root.style.display = 'flex';
    this.root.style.flexDirection = 'column';
    this.root.style.width = '100%';
    this.root.style.height = '100%';
  }

  async renderInternal(): Promise<void> {
    if (!this.moleculeColumn || !this.mpoColumn) {
      ui.empty(this.root);
      this.root.append(ui.divText('Molecule or MPO column not found'));
      return;
    }

    ui.empty(this.root);

    const rowCount = Math.min(this.limit, this.dataFrame.rowCount);
    let sortedRows = Array.from({length: this.dataFrame.rowCount}, (_, i) => i);

    if (this.mode === 'best')
      sortedRows.sort((a, b) => this.mpoColumn.get(b) - this.mpoColumn.get(a));
    else
      sortedRows.sort((a, b) => this.mpoColumn.get(a) - this.mpoColumn.get(b));

    sortedRows = sortedRows.slice(0, rowCount);

    const grids: HTMLElement[] = [];
    this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;

    for (const i of sortedRows) {
      const mol = this.moleculeColumn.get(i);
      const score = this.mpoColumn.get(i);

      const molDiv = renderMolecule(mol, {
        width: this.sizesMap[this.size].width,
        height: this.sizesMap[this.size].height,
      });

      const scoreDiv = ui.divText(score.toFixed(2));
      scoreDiv.style.marginTop = '6px';
      const propsDiv = ui.divV([scoreDiv], {style: {alignItems: 'center'}});

      const grid = ui.div([molDiv, propsDiv], {style: {position: 'relative'}});
      const divClass = 'd4-flex-col';
      $(grid).addClass(divClass);
      grids.push(grid);
    }

    const panel = ui.divH(grids, 'chem-viewer-grid');
    this.root.appendChild(ui.panel([ui.divV([panel])]));
  }
}
