import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import {PeptidesModel} from '../model';
import {isGridCellInvalid} from './misc';

export class InvariantMap extends DG.Filter {
  model!: PeptidesModel;
  chosenCells: Map<string, string[]> = new Map();

  constructor() {
    super();
  }

  get filterSummary(): string {

  }

  async attach(df: DG.DataFrame): Promise<void> {
    super.attach(df);
    this.model = await PeptidesModel.getInstance(df);
    this.render();
  }

  applyFilter(): void {

  }

  render(): void {
    const invariantDf = this.model.statsDf.groupBy([C.COLUMNS_NAMES.MONOMER])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.COUNT, '')
      .aggregate();
    const invariantGrid = invariantDf.plot.grid();
    this.chosenCells = new Map();
    for (const col of invariantDf.columns)
      if (col.name != C.COLUMNS_NAMES.MONOMER)
        this.chosenCells.set(col.name, []);

    invariantGrid.root.addEventListener('click', (ev) => {
      const gridCell = invariantGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name == C.COLUMNS_NAMES.MONOMER)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = invariantDf.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!);
      const aarList = this.chosenCells.get(position)!;
      const aarIndex = aarList.indexOf(aar);

      if (aarIndex == -1)
        aarList.splice(aarIndex, 1);
      else
        aarList.push(position);

      invariantGrid.invalidate();
    });

    invariantGrid.onCellRender.subscribe((args) => {
      const gc = args.cell;
      const tableColName = gc.tableColumn?.name;
      const tableRowIndex = gc.tableRowIndex;
      if (isGridCellInvalid(gc) || tableColName == C.COLUMNS_NAMES.MONOMER || !tableRowIndex || tableRowIndex == -1)
        return;

      const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
        tableColName : invariantDf.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
      const currentAAR: string = invariantDf.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);

      const aarSelection = this.chosenCells.get(currentPosition)!;
      if (aarSelection && aarSelection.includes(currentAAR)) {
        const canvasContext = args.g;
        const bound = args.bounds;
        canvasContext.strokeStyle = '#000';
        canvasContext.lineWidth = 1;
        canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
      }
    });

    this.root.appendChild(invariantGrid.root);
  }
}
