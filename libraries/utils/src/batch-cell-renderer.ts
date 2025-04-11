import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { LruCache } from "datagrok-api/dg";
import { emit } from 'process';


export abstract class BatchCellRenderer<T> extends DG.GridCellRenderer {
  currentTimeout: any = null;
  isRunning: boolean = false;
  loadedValues: Map<string, T> = new Map<string, T>();
  cellsToLoad: DG.GridCell[] = [];
  cellsToLoadOnTimeoutComplete: DG.GridCell[] = [];
  cache: LruCache<string, T>;

  constructor() {
    super();
    this.cache = new LruCache<string, T>();
  }

  abstract loadData(keys: string[]): Promise<Map<string, T>>;

  async retrieveBatch(gridCell: DG.GridCell) {
    let found = this.cellsToLoad.filter((e) => {
      return e.grid === gridCell.grid && e.value === gridCell.value && e.gridColumn?.name === gridCell.gridColumn?.name && e.tableRowIndex === gridCell.tableRowIndex;
    }) ?? [];

    if (found.length > 0)
      return;

    if (this.isRunning) {
      this.cellsToLoadOnTimeoutComplete.push(gridCell);
      return;
    }

    this.cellsToLoad.push(gridCell);
    if (this.currentTimeout)
      clearTimeout(this.currentTimeout);

    this.currentTimeout = setTimeout(async () => {
      this.isRunning = true;
      this.loadedValues = await this.loadData(this.cellsToLoad.map((e) => e.cell.valueString));

      for (const [key, value] of this.loadedValues) {
        this.cache.set(key, value);
      }

      this.cellsToLoad.forEach((x) => {
        if (!this.loadedValues.has(x.cell.valueString))
          this.cache.set(x.cell.valueString, null);
        x.grid.invalidate();
      })

      this.loadedValues.clear();
      this.cellsToLoad = [];
      this.currentTimeout = null;
      this.isRunning = false;
      if (this.cellsToLoadOnTimeoutComplete.length > 0)
        this.cellsToLoadOnTimeoutComplete.forEach((cell) => { this.retrieveBatch(cell) });
      this.cellsToLoadOnTimeoutComplete = [];
    });
  }
}