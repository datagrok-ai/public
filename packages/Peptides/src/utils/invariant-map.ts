import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import {PeptidesModel} from '../model';
import {isGridCellInvalid} from './misc';

const CELL_SIZE = 20; // 20px cell height and width

export class InvariantMap extends DG.Filter {
  model: PeptidesModel | null = null;
  chosenCells: {[position: string]: string[]} = {};

  constructor() {
    super();
  }

  get caption(): string {return 'Invariant Map';}

  get filterSummary(): string {
    let summary = '';
    for (const [pos, aarList] of Object.entries(this.chosenCells))
      summary += `${pos}: ${aarList}\n`;

    return summary;
  }

  get isFiltering(): boolean {return true && super.isFiltering;}

  get isReadyToApplyFilter(): boolean {return this.model != null;}

  async attach(df: DG.DataFrame): Promise<void> {
    super.attach(df);
    this.model = await PeptidesModel.getInstance(df);
    this.render(true);
  }

  saveState(): any {
    const state = super.saveState();
    state.chosenCells = JSON.stringify(this.chosenCells);
    return state;
  }

  applyState(state: any): void {
    super.applyState(state);
    if (state.chosenCells) {
      this.chosenCells = JSON.parse(state.chosenCells);
      this.render();
    }
  }

  applyFilter(): void {
    this.dataFrame?.filter.init((bitsetIndex) => {
      for (const [position, aarList] of Object.entries(this.chosenCells)) {
        if (aarList.length != 0 && !aarList.includes(this.dataFrame!.get(position, bitsetIndex)))
          return false;
      }
      return true;
    });
  }

  render(initChosenCells: boolean = false): void {
    if (this.model == null)
      return;

    const invariantDf = this.model.statsDf.groupBy([C.COLUMNS_NAMES.MONOMER])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.COUNT, '')
      .aggregate();
    invariantDf.getCol(C.COLUMNS_NAMES.MONOMER).semType = C.SEM_TYPES.MONOMER;
    const orderedColNames = invariantDf.columns.names().sort((a, b) => {
      const aInt = parseInt(a);
      const bInt = parseInt(b);
      if (isNaN(aInt))
        return -1;
      else if (isNaN(bInt))
        return 1;
      return aInt - bInt;
    });

    // Create grid and set properties
    const invariantGrid = invariantDf.plot.grid();
    const gridCols = invariantGrid.columns;
    gridCols.rowHeader!.visible = false;
    gridCols.setOrder(orderedColNames);

    for (let gridColIndex = 0; gridColIndex < gridCols.length; ++gridColIndex) {
      const gridCol = gridCols.byIndex(gridColIndex)!
      if (gridCol.name != C.COLUMNS_NAMES.MONOMER) 
        gridCol.width = CELL_SIZE;
    }

    if (initChosenCells) {
      this.chosenCells = {};
      for (const col of invariantDf.columns)
        if (col.name != C.COLUMNS_NAMES.MONOMER)
          this.chosenCells[col.name] = [];
    }

    invariantGrid.root.addEventListener('click', (ev) => {
      invariantDf.currentRowIdx = -1;
      const gridCell = invariantGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name == C.COLUMNS_NAMES.MONOMER)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = invariantDf.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!);
      const aarList = this.chosenCells[position];
      const aarIndex = aarList.indexOf(aar);

      if (aarIndex != -1)
        aarList.splice(aarIndex, 1);
      else
        aarList.push(aar);

      invariantGrid.invalidate();
      this.applyFilter();
    });

    invariantGrid.onCellRender.subscribe((args) => {
      //FIXME: for some reason it doesn't work when I set right away
      const gridProps = invariantGrid.props;
      gridProps.allowBlockSelection = false;
      gridProps.allowColSelection = false;
      gridProps.allowRowSelection = false;
      gridProps.allowEdit = false;
      gridProps.rowHeight = CELL_SIZE;

      const gc = args.cell;
      const tableColName = gc.tableColumn?.name;
      const tableRowIndex = gc.tableRowIndex;

      if (isGridCellInvalid(gc) || tableColName == C.COLUMNS_NAMES.MONOMER ||
        tableRowIndex == null || tableRowIndex == -1)
        return;

      const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
        tableColName : invariantDf.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
      const currentAAR: string = invariantDf.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);
      const canvasContext = args.g;
      const bound = args.bounds;

      canvasContext.font = '10px Roboto';
      canvasContext.textAlign = 'center';
      canvasContext.textBaseline = 'middle';
      canvasContext.fillStyle = '#000';
      canvasContext.fillText(gc.cell.value, bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

      const aarSelection = this.chosenCells[currentPosition];
      if (aarSelection.includes(currentAAR)) {
        canvasContext.strokeStyle = '#000';
        canvasContext.lineWidth = 1;
        canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
      }
      args.preventDefault();
    });

    const gridHost = ui.box(invariantGrid.root);
    gridHost.style.width = '100%';
    this.root.appendChild(gridHost);
  }
}
