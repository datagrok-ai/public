/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {fromEvent, Subject} from 'rxjs';
import {take, takeUntil} from 'rxjs/operators';
import {PlateWidget} from './plate-widget';

const DRAG_THRESHOLD = 5;

export interface WellSelectionEvent {
  wells: Array<{row: number, col: number, dataIndex: number}>;
  isDrag: boolean;
}

export type SelectionMode = 'select' | 'passthrough';

export class PlateSelectionController {
  private _enabled: boolean = false;
  private _mode: SelectionMode = 'select';
  private _isDragging: boolean = false;
  private _selectionRect: DG.Rect | null = null;
  private _dragStartPoint: DG.Point | null = null;

  private selectionCompleteSubject = new Subject<WellSelectionEvent>();
  public get onSelectionComplete() { return this.selectionCompleteSubject.asObservable(); }

  private _onDestroy = new Subject<void>();

  constructor(
    private plateWidget: PlateWidget,
    private canvas: HTMLCanvasElement
  ) {}

  get enabled(): boolean { return this._enabled; }

  setMode(mode: SelectionMode): void {
    this._mode = mode;
  }

  enable(): void {
    if (this._enabled) return;
    this._enabled = true;
    this.initSelectionEvents();
  }

  disable(): void {
    if (!this._enabled) return;
    this._enabled = false;
    this.clearSelection();
    this._onDestroy.next();
  }

  clearSelection(): void {
    this.plateWidget.plate.data.selection.setAll(false, true);
    this.clearSelectionRect();
  }

  private initSelectionEvents(): void {
    const interactionElement = this.canvas;
    if (!interactionElement) return;

    fromEvent<MouseEvent>(interactionElement, 'mousedown')
      .pipe(takeUntil(this._onDestroy))
      .subscribe((e: MouseEvent) => {
        if (!this._enabled || e.button !== 0) return;

        this._dragStartPoint = new DG.Point(e.offsetX, e.offsetY);
        this._isDragging = false; // Reset drag state
        this._selectionRect = null;

        const mouseMoveStream = fromEvent<MouseEvent>(document, 'mousemove').pipe(takeUntil(this._onDestroy));
        const mouseUpStream = fromEvent<MouseEvent>(document, 'mouseup').pipe(takeUntil(this._onDestroy), take(1));

        mouseMoveStream.pipe(takeUntil(mouseUpStream)).subscribe((move_e: MouseEvent) => {
          if (!this._dragStartPoint) return;

          const dx = Math.abs(move_e.offsetX - this._dragStartPoint.x);
          const dy = Math.abs(move_e.offsetY - this._dragStartPoint.y);

          if (!this._isDragging && (dx > DRAG_THRESHOLD || dy > DRAG_THRESHOLD)) {
            this._isDragging = true;
            this._selectionRect = new DG.Rect(this._dragStartPoint.x, this._dragStartPoint.y, 0, 0);
          }

          if (this._isDragging && this._selectionRect) {
            move_e.preventDefault();
            move_e.stopPropagation();

            const canvasBounds = interactionElement.getBoundingClientRect();
            const currentX = move_e.clientX - canvasBounds.left;
            const currentY = move_e.clientY - canvasBounds.top;

            this._selectionRect.width = currentX - this._selectionRect.x;
            this._selectionRect.height = currentY - this._selectionRect.y;

            this.drawSelectionRect();
          }
        });

        mouseUpStream.subscribe((up_e: MouseEvent) => {
          if (this._isDragging) {
            this.finalizeDragSelection();
            const selectedWells = this.getSelectedWells();

            if (selectedWells.length > 0) {
              this.selectionCompleteSubject.next({
                wells: selectedWells,
                isDrag: true
              });
            }
          } else {
            const well = this.plateWidget.hitTest(up_e.offsetX, up_e.offsetY);

            if (well) {
              this.plateWidget.fireWellClick(well.row, well.col);
            } else {
            }
          }

          // Reset state
          this._isDragging = false;
          this.clearSelectionRect();
          this._selectionRect = null;
          this._dragStartPoint = null;
          this.plateWidget.grid.invalidate();
        });
      });
  }

  private finalizeDragSelection(): void {
    if (!this._selectionRect || !this.plateWidget.grid) return;

    const selection = this.plateWidget.plate.data.selection;
    const plate = this.plateWidget.plate;
    const r = this._selectionRect;

    const normalizedRect = new DG.Rect(
      r.width < 0 ? r.x + r.width : r.x,
      r.height < 0 ? r.y + r.height : r.y,
      Math.abs(r.width),
      Math.abs(r.height)
    );

    for (let row = 0; row < plate.rows; row++) {
      for (let col = 0; col < plate.cols; col++) {
        const gridCol = this.plateWidget.grid.columns.byIndex(col + 1);
        if (gridCol) {
          const cell = this.plateWidget.grid.cell(gridCol.name, row);
          if (cell && normalizedRect.contains(cell.bounds.midX, cell.bounds.midY))
            selection.set(plate._idx(row, col), true, false);
        }
      }
    }

    selection.fireChanged();
  }

  private getSelectedWells(): Array<{row: number, col: number, dataIndex: number}> {
    const selection = this.plateWidget.plate.data.selection;
    const wells: Array<{row: number, col: number, dataIndex: number}> = [];

    for (const dataIndex of selection.getSelectedIndexes()) {
      const [row, col] = this.plateWidget.plate.rowIndexToExcel(dataIndex);
      wells.push({row, col, dataIndex});
    }

    return wells;
  }

  private drawSelectionRect(): void {
    if (!this.canvas || !this._selectionRect) return;

    const g = this.canvas.getContext('2d')!;
    g.clearRect(0, 0, this.canvas.width, this.canvas.height);

    Object.assign(g, {
      strokeStyle: 'rgba(0, 128, 255, 0.7)',
      fillStyle: 'rgba(0, 128, 255, 0.2)',
      lineWidth: 1
    });

    g.strokeRect(
      this._selectionRect.x,
      this._selectionRect.y,
      this._selectionRect.width,
      this._selectionRect.height
    );
    g.fillRect(
      this._selectionRect.x,
      this._selectionRect.y,
      this._selectionRect.width,
      this._selectionRect.height
    );
  }

  private clearSelectionRect(): void {
    if (!this.canvas) return;
    this.canvas.getContext('2d')!.clearRect(0, 0, this.canvas.width, this.canvas.height);
  }

  destroy(): void {
    this.disable();
    this._onDestroy.complete();
    this.selectionCompleteSubject.complete();
  }
}
