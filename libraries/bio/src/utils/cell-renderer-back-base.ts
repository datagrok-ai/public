import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ILogger} from './logger';
import {Unsubscribable} from 'rxjs';
import {CellRendererBackAsyncBase} from './cell-renderer-async-base';

type GridCellRendererTemp<TBack> = {
  rendererBack: TBack;
}

export function getGridCellRendererBack<TValue, TBack extends CellRendererBackBase<TValue>>(
  gridCell: DG.GridCell
): [DG.GridColumn | null, DG.Column<TValue>, GridCellRendererTemp<TBack>] {
  let temp: GridCellRendererTemp<TBack> | null = null;

  let gridCol: DG.GridColumn | null = null;
  try { gridCol = gridCell.gridColumn;} catch { gridCol = null;}
  if (gridCol) temp = gridCol.temp as GridCellRendererTemp<TBack>;

  // let tableCol: DG.Column<TValue> | null = null;
  // try { tableCol = gridCell.cell.column; } catch { tableCol = null;}
  // if (tableCol) temp = tableCol.temp as GridCellRendererTemp<TBack>;
  const tableCol: DG.Column<TValue> = gridCell.cell.column;
  temp = temp ?? tableCol.temp;

  if (!temp)
    throw new Error(`Grid cell renderer back store (GridColumn or Column) not found.`);

  return [gridCol, tableCol, temp];
}

export abstract class CellRendererBackBase<TValue> {
  protected subs: Unsubscribable[] = [];
  protected destroyed: boolean = false;

  protected constructor(
    protected readonly gridCol: DG.GridColumn | null,
    protected readonly tableCol: DG.Column<TValue>,
    protected readonly logger: ILogger,
  ) {
    this.reset();
    if (this.tableCol && this.tableCol.dataFrame) {
      this.subs.push(this.tableCol.dataFrame.onDataChanged.subscribe(() => {
        try { this.reset(); } catch (err) { this.logger.error(err); }
      }));
    }

    if (this.tableCol) {
      this.subs.push(grok.events.onTableRemoved.subscribe((eventData) => {
        try {
          const eventDf: DG.DataFrame = eventData.args.dataFrame;
          if (this.tableCol?.dataFrame.id === eventDf.id && !this.destroyed)
            this.destroy();
        } catch (err) { this.logger.error(err); }
      }));
    }

    if (this.gridCol) {
      this.subs.push(grok.events.onViewRemoving.subscribe((eventData) => {
        try {
          const eventView = eventData.args.view;
          if (
            this.gridCol?.dart &&
            this.gridCol.grid && this.gridCol?.grid?.dart &&
            this.gridCol?.grid?.view?.id === eventView.id && !this.destroyed
          ) this.destroy();
        } catch (err: any) { this.logger.error(err); }
      }));
    }
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = ++CellRendererBackBase.viewerCounter;

  protected toLog(): string {
    return `${this.constructor.name}<${this.viewerId}>`;
  }

  protected destroy(): void {
    for (const sub of this.subs)
      sub.unsubscribe();
    this.destroyed = true;
  }

  protected invalidateGrid(): void {
    if (this.gridCol && this.gridCol.dart) this.gridCol.grid?.invalidate();
  }

  protected abstract reset(): void;
}
