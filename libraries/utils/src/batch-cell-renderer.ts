import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { LruCache } from "datagrok-api/dg";
import { emit } from 'process';


export abstract class BatchCellRenderer<T> extends DG.GridCellRenderer {
  currentTimeout: any = null;
  isRunning: boolean = false;
  cellsToLoad: { value: DG.GridCell, renderFuncs: Map<HTMLCanvasElement, () => (void)> }[] = [];
  cellsToLoadOnTimeoutComplete: { value: DG.GridCell, renderFunc: () => (void), canvas: HTMLCanvasElement }[] = [];
  cache: LruCache<string, T>;
  delay: number = 100;

  constructor() {
    super();
    this.cache = new LruCache<string, T>();
  }

  abstract loadData(keys: string[]): Promise<Map<string, T>>;

  /** 
    gridCell - cell that has to be loaded
    funcToRender - closure function to render cell
    canvas - cell's canvas
  **/
  async retrieveBatch(gridCell: DG.GridCell, funcToRender: () => void, canvas: HTMLCanvasElement) {
    const found = this.cellsToLoad.find((e) => {
      return e.value.value === gridCell.value && e.value.gridColumn === gridCell.gridColumn && e.value.tableRowIndex === gridCell.tableRowIndex;;
    });

    if (found) {
      found.renderFuncs.set(canvas, funcToRender);
      return;
    }

    if (this.isRunning) {
      this.cellsToLoadOnTimeoutComplete.push({ value: gridCell, renderFunc: funcToRender, canvas: canvas });
      return;
    }

    this.cellsToLoad.push({ value: gridCell, renderFuncs: new Map().set(canvas, funcToRender)});
    if (this.currentTimeout)
      clearTimeout(this.currentTimeout);

    this.currentTimeout = setTimeout(async () => {
      this.isRunning = true;
      
      const loadedValues = await this.loadData(Array.from(new Set(this.cellsToLoad.map((e) => e.value.cell.valueString))));

      for (const [key, value] of loadedValues) {
        this.cache.set(key, value);
      }

      this.cellsToLoad.forEach((x) => {
        if (!loadedValues.has(x.value.cell.valueString))
          this.cache.set(x.value.cell.valueString, null);
        for(const func of x.renderFuncs.values()){
          func();
        }
      })

      this.cellsToLoad = [];
      this.currentTimeout = null;
      this.isRunning = false;

      if (this.cellsToLoadOnTimeoutComplete.length > 0)
        this.cellsToLoadOnTimeoutComplete.forEach((cell) => { this.retrieveBatch(cell.value, cell.renderFunc, cell.canvas) });

      this.cellsToLoadOnTimeoutComplete = [];
    }, this.delay);
  }
}