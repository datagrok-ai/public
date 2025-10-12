/* eslint-disable prefer-const */
/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {mapFromRow, safeLog, toExcelPosition} from './utils';
import {Plate, PLATE_OUTLIER_WELL_NAME} from './plate';
import {fromEvent, Subject, Subscription} from 'rxjs';
import {filter, take, takeUntil} from 'rxjs/operators';
import './plate-widget.css';
import * as grok from 'datagrok-api/grok';
import {IPlateWellValidator, plateWellValidators} from './plate-well-validators';

const colorPalette = [
  DG.Color.fromHtml('#e41a1c'), DG.Color.fromHtml('#377eb8'), DG.Color.fromHtml('#4daf4a'),
  DG.Color.fromHtml('#984ea3'), DG.Color.fromHtml('#ff7f00'), DG.Color.fromHtml('#f781bf'),
  DG.Color.fromHtml('#a65628'), DG.Color.fromHtml('#999999'), DG.Color.fromHtml('#66c2a5'),
  DG.Color.fromHtml('#fc8d62'), DG.Color.fromHtml('#8da0cb'), DG.Color.fromHtml('#e78ac3'),
];

export type AnalysisOptions = {
  roleName: string,
  concentrationName: string,
  valueName: string
  controlColumns: string[],
  normalize: boolean,
  autoFilterOutliers: boolean,
  submitAction?: (plate: Plate, curvesDf: DG.DataFrame) => void,
  categorizeFormula: string,
  statisticsColumns: string[],
  plateStatistics?: { [key: string]: (plate: Plate) => string | number }
}

const colorScheme = [DG.Color.white, DG.Color.gray];
export const dimensions = new Map([
  [96, {rows: 8, cols: 12}],
  [384, {rows: 16, cols: 24}],
  [1536, {rows: 32, cols: 48}],
]);

export interface IAnalysisWidgetCoordinator {
  onPlateDataChanged(changeType: 'outlier' | 'data' | 'layer', details: any): void;
  refreshAnalysisView(): void;
}

type InteractionMode = 'default' | 'outlier';
const DRAG_THRESHOLD = 5;


export class PlateWidget extends DG.Widget {
  _colorColumn?: DG.Column | undefined;
  _plate: Plate = new Plate(8, 12);
  grid: DG.Grid = DG.Viewer.heatMap(DG.DataFrame.create());
  colorColumnName: string = '';
  detailsDiv?: HTMLElement;
  wellDetailsDiv?: HTMLElement;
  plateDetailsDiv?: HTMLElement;
  plateActionsDiv?: HTMLElement;
  tabs: DG.TabControl = DG.TabControl.create(true);
  tabsContainer: HTMLElement = ui.divH([], 'assay-plates--plate-widget__tabs-container');
  _editable: boolean = false;
  mapFromRowFunc: (row: DG.Row) => Record<string, any> = mapFromRow;
  grids: Map<string, DG.Grid> = new Map();
  wellValidators: IPlateWellValidator[] = plateWellValidators;
  wellValidationErrors: Map<string, string[]> = new Map();

  private hoveredCell: { row: number, col: number } | null = null;
  private _isDragging: boolean = false;
  private _selectionRect: DG.Rect | null = null;
  private _canvas: HTMLCanvasElement | null = null;
  private _onDestroy = new Subject<void>();
  private analysisCoordinators: IAnalysisWidgetCoordinator[] = [];
  private outlierSubscription?: Subscription;
  private interactionMode: InteractionMode = 'default';
  private _dragStartPoint: DG.Point | null = null;

  // private roleColorMap: Map<string, number> = new Map();
  // private nextColorIndex = 0;


  get editable() { return this._editable; }
  set editable(x: boolean) {
    this._editable = x;
    this.syncGrids();
  }

  get plate(): Plate { return this._plate; }
  set plate(p: Plate) {
    this._plate = p;
    if (this.outlierSubscription)
      this.outlierSubscription.unsubscribe();

    this.outlierSubscription = p.onOutlierChanged.subscribe((change) => {
      this.grid.invalidate();
      this.notifyAnalysisCoordinators('outlier', change);
    });

    // this.roleColorMap.clear();
    // this.nextColorIndex = 0;
    this.refresh();
  }

  private isPointInWell(x: number, y: number, gridRow: number, gridCol: number): boolean {
    if (gridRow < 0 || gridCol <= 0) return false;
    const cell = this.grid.cell(this.grid.columns.byIndex(gridCol)!.name, gridRow);
    if (!cell) return false;
    const bounds = cell.bounds;
    const centerX = bounds.midX;
    const centerY = bounds.midY;
    const radius = Math.min(bounds.height / 2, bounds.width / 2) * 0.8;
    const distance = Math.sqrt(Math.pow(x - centerX, 2) + Math.pow(y - centerY, 2));
    return distance <= radius;
  }

  public addAnalysisTab(name: string, content: HTMLElement): void {
    const pane = this.tabs.addPane(name, () => content);
    this.tabs.currentPane = pane;
  }
  private hitTest(canvasX: number, canvasY: number): { row: number, col: number } | null {
    for (let row = 0; row < this.plate.rows; row++) {
      for (let col = 0; col < this.plate.cols; col++) {
        if (this.isPointInWell(canvasX, canvasY, row, col + 1))
          return {row, col};
      }
    }
    return null;
  }

  constructor() {
    super(ui.div([], 'assay-plates--plate-widget'));
    this.tabs.root.classList.add('assay-plates--plate-widget__tabs');

    this.tabsContainer.appendChild(this.tabs.root);

    this.tabs.onTabChanged.subscribe(() => {
      const isOutlierMode = this.tabs.currentPane.name.includes('Outliers');
      this.interactionMode = isOutlierMode ? 'outlier' : 'default';
      this.grid.root.classList.toggle('outlier-mode', isOutlierMode);
      this.plate.data.selection.setAll(false, true);
      this.grid.invalidate();
    });

    const mainContainer = ui.divV([this.tabsContainer], 'assay-plates--plate-widget__main-content');
    this.root.appendChild(mainContainer);

    this.setupGrid();

    ui.tools.waitForElementInDom(this.grid.root).then(() => {
      this._canvas = this.grid.overlay;
      this.initSelectionEvents();
    });
  }

  private setupGrid() {
    Object.assign(this.grid.props, {
      allowRowSelection: false, allowColSelection: false, allowBlockSelection: false,
      showCurrentRowIndicator: false, showMouseOverRowIndicator: false, showCurrentCellOutline: false,
      showHeatmapScrollbars: false, colHeaderHeight: 30, allowColHeaderResizing: false,
      allowColResizing: false, allowRowResizing: false,
    });
    this.grid.onCellRender.subscribe((args) => this.renderCell(args, true));
    this.grid.onCellRendered.subscribe((args) => {
      const cell = args.cell;
      if (cell.gridRow !== null && cell.gridRow >= 0 && cell.gridColumn.idx > 0) {
        try {
          const pos = toExcelPosition(cell.gridRow, cell.gridColumn.idx - 1);
          const errors = this.wellValidationErrors.get(pos);
          if (errors && errors.length > 0)
            DG.Paint.marker(args.g, DG.MARKER_TYPE.CROSS_BORDER, cell.bounds.midX, cell.bounds.midY, DG.Color.red, 10);
        } catch (e) { console.error(e); }
      }
    });

    this.setupHoverEvents(this.grid);

    this.grid.onCellClick.subscribe((gc: DG.GridCell) => {
      if (!gc.isTableCell) return;
      if (this.interactionMode === 'outlier' && this.editable) {
        const row = gc.gridRow!;
        const col = gc.gridColumn.idx - 1;
        const currentState = this.plate.isOutlier(row, col);

        this.plate.markOutlier(row, col, !currentState);
        this.grid.invalidate();
      }
    });
  }
  private setupHoverEvents(grid: DG.Grid) {
    grid.onCellMouseEnter.subscribe((gc: DG.GridCell) => {
      if (gc.isTableCell) {
        gc.grid.root.style.cursor = 'pointer';
        this.hoveredCell = {row: gc.gridRow!, col: gc.gridColumn.idx - 1};
        grid.invalidate();
      }
    });

    grid.onCellMouseLeave.subscribe((gc: DG.GridCell) => {
      gc.grid.root.style.cursor = 'default';
      this.hoveredCell = null;
      grid.invalidate();
    });
  }

  registerAnalysisCoordinator(coordinator: IAnalysisWidgetCoordinator): void {
    this.analysisCoordinators.push(coordinator);
  }

  private notifyAnalysisCoordinators(changeType: 'outlier' | 'data' | 'layer', details: any): void {
    this.analysisCoordinators.forEach((coord) => coord.onPlateDataChanged(changeType, details));
  }

  private initSelectionEvents() {
    const interactionElement = this.grid.overlay;
    if (!interactionElement) return;

    fromEvent<MouseEvent>(interactionElement, 'mousedown').pipe(takeUntil(this._onDestroy)).subscribe((e: MouseEvent) => {
      if (!this.editable || e.button !== 0) return;
      this._dragStartPoint = new DG.Point(e.offsetX, e.offsetY);
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
          const selection = this.plate.data.selection;
          if (selection.trueCount > 0) {
            if (this.interactionMode === 'default') {
              // this.showRoleAssignmentPopup();
            } else if (this.interactionMode === 'outlier') {
              for (const i of selection.getSelectedIndexes()) {
                const [row, col] = this.plate.rowIndexToExcel(i);

                this.plate.markOutlier(row, col, true);
              }
              selection.setAll(false, true);
            }
          }
        } else if (this._dragStartPoint) {
          const well = this.hitTest(up_e.offsetX, up_e.offsetY);
          if (well) {
            if (this.interactionMode === 'default') {
              const dataIndex = this.plate._idx(well.row, well.col);
              const selection = this.plate.data.selection;
              selection.set(dataIndex, !selection.get(dataIndex), true);
            } else if (this.interactionMode === 'outlier') {
              const currentState = this.plate.isOutlier(well.row, well.col);

              this.plate.markOutlier(well.row, well.col, !currentState);
            }
          }
        }

        this._isDragging = false;
        this.clearSelectionRect();
        this._selectionRect = null;
        this._dragStartPoint = null;
        this.grid.invalidate();
      });
    });
  }

  private finalizeDragSelection() {
    if (!this._selectionRect || !this.grid) return;
    const selection = this.plate.data.selection;
    const r = this._selectionRect;
    const normalizedRect = new DG.Rect(
      r.width < 0 ? r.x + r.width : r.x,
      r.height < 0 ? r.y + r.height : r.y,
      Math.abs(r.width),
      Math.abs(r.height)
    );
    for (let row = 0; row < this.plate.rows; row++) {
      for (let col = 0; col < this.plate.cols; col++) {
        const gridCol = this.grid.columns.byIndex(col + 1);
        if (gridCol) {
          const cell = this.grid.cell(gridCol.name, row);
          if (cell && normalizedRect.contains(cell.bounds.midX, cell.bounds.midY))
            selection.set(this.plate._idx(row, col), true, false);
        }
      }
    }
    selection.fireChanged();
  }

  // private showRoleAssignmentPopup(): void {
  //   const selection = this.plate.data.selection;
  //   if (selection.trueCount === 0) return;
  //   const roleInput = ui.input.string('Role', {value: 'Sample'});
  //   let popup: Element;
  //   let anchorDiv: HTMLElement;

  //   const assignAction = () => {
  //     const roleName = roleInput.value.trim();
  //     if (!roleName) {
  //       grok.shell.warning('Role name cannot be empty.');
  //       return;
  //     }
  //     const existingRoleCols = this.plate.getLayersByType(LayerType.LAYOUT);
  //     const selectedIndexes = selection.getSelectedIndexes();

  //     for (const colName of existingRoleCols) {
  //       const col = this.plate.data.col(colName);
  //       if (col) {
  //         for (const i of selectedIndexes)
  //           col.set(i, false, false);
  //       }
  //     }

  //     let newRoleCol = this.plate.data.col(roleName);
  //     if (newRoleCol === null) {
  //       newRoleCol = this.plate.data.columns.addNewBool(roleName);
  //       this.plate.registerLayer(roleName, LayerType.LAYOUT, 'user-assignment');
  //     }
  //     for (const i of selectedIndexes)
  //       newRoleCol.set(i, true, false);

  //     this.plate.data.fireValuesChanged();
  //     // this.updateRoleSummary();
  //     this.refresh();
  //     selection.setAll(false, true);
  //     this.grid.invalidate();

  //     if (popup && anchorDiv) {
  //       popup.remove();
  //       anchorDiv.remove();
  //     }
  //   };

  //   const popupContent = ui.divV([
  //     ui.h3(`${selection.trueCount} wells selected`),
  //     roleInput,
  //     ui.button('ASSIGN', assignAction)
  //   ], {style: {padding: '10px'}});

  //   roleInput.input.addEventListener('keydown', (e) => {
  //     if (e.key === 'Enter') {
  //       e.preventDefault();
  //       assignAction();
  //     }
  //   });

  //   const selectedIndices = selection.getSelectedIndexes();
  //   let minRow = this.plate.rows; let maxRow = -1; let minCol = this.plate.cols; let maxCol = -1;
  //   for (const idx of selectedIndices) {
  //     const row = Math.floor(idx / this.plate.cols);
  //     const col = idx % this.plate.cols;
  //     if (row < minRow) minRow = row;
  //     if (row > maxRow) maxRow = row;
  //     if (col < minCol) minCol = col;
  //     if (col > maxCol) maxCol = col;
  //   }

  //   const firstCell = this.grid.cell(this.grid.columns.byIndex(minCol + 1)!.name, minRow);
  //   const lastCell = this.grid.cell(this.grid.columns.byIndex(maxCol + 1)!.name, maxRow);
  //   const selectionBounds = firstCell.bounds.union(lastCell.bounds);
  //   const canvasBounds = this.grid.overlay.getBoundingClientRect();

  //   anchorDiv = ui.div('', {
  //     style: {
  //       position: 'absolute',
  //       left: `${canvasBounds.left + selectionBounds.x + selectionBounds.width / 2}px`,
  //       top: `${canvasBounds.top + selectionBounds.y}px`,
  //       width: '0px', height: '0px',
  //     }
  //   });
  //   document.body.appendChild(anchorDiv);
  //   popup = ui.showPopup(popupContent, anchorDiv);
  // }

  private drawSelectionRect() {
    if (!this._canvas || !this._selectionRect) return;
    const g = this._canvas.getContext('2d')!;
    g.clearRect(0, 0, this._canvas.width, this._canvas.height);
    Object.assign(g, {strokeStyle: 'rgba(0, 128, 255, 0.7)', fillStyle: 'rgba(0, 128, 255, 0.2)', lineWidth: 1});
    g.strokeRect(this._selectionRect.x, this._selectionRect.y, this._selectionRect.width, this._selectionRect.height);
    g.fillRect(this._selectionRect.x, this._selectionRect.y, this._selectionRect.width, this._selectionRect.height);
  }

  private clearSelectionRect() {
    if (!this._canvas) return;
    this._canvas.getContext('2d')!.clearRect(0, 0, this._canvas.width, this._canvas.height);
  }

  static fromPlate(plate: Plate, useSimpleView: boolean = false) {
    if (useSimpleView) {
      const pw = new PlateWidget();
      pw.plate = plate;
      return pw;
    }
    return PlateWidget.detailedView(plate);
  }

  static detailedView(plate: Plate): PlateWidget {
    const pw = new PlateWidget();
    pw.plate = plate;

    ui.empty(pw.root);

    pw.detailsDiv = ui.divV([], 'assay-plates--plate-widget__details');
    pw.wellDetailsDiv = ui.div();
    pw.detailsDiv.appendChild(pw.wellDetailsDiv);

    const mainContainer = ui.divH([
      pw.tabs.root,
      pw.detailsDiv,
    ], 'assay-plates--plate-widget__main-container');

    pw.root.appendChild(mainContainer);

    pw.grid.onCurrentCellChanged
      .pipe(filter((gc) => gc.isTableCell && gc.gridRow !== null))
      .subscribe((gc) => {
        if (pw.wellDetailsDiv) {
          ui.empty(pw.wellDetailsDiv);
          const map = pw.mapFromRowFunc(pw.plate.data.rows.get(plate._idx(gc.gridRow!, gc.gridColumn.idx - 1)));
          pw.wellDetailsDiv.appendChild(ui.tableFromMap(map));
        }
      });
    return pw;
  }

  initPlateGrid(grid: DG.Grid, isSummary: boolean = false) {
    Object.assign(grid.props, {
      showHeatmapScrollbars: false, allowColReordering: false, allowRowReordering: false,
      allowSorting: false, allowEdit: !isSummary, showRowGridlines: false, showColumnGridlines: false,
      heatmapColors: false, colHeaderHeight: 25, allowRowSelection: false, allowColSelection: false,
      allowBlockSelection: false, showCurrentRowIndicator: false, showMouseOverRowIndicator: false,
      showCurrentCellOutline: false, allowColHeaderResizing: false, allowColResizing: false,
      allowRowResizing: false, selectedRowsColor: 0x00000000, mouseOverRowColor: 0x00000000,
      currentRowColor: 0x00000000,
    });
    grid.onCellRender.subscribe((args) => this.renderCell(args, isSummary));
    this.setupHoverEvents(grid);
  }

  // private getDisplayName(layerName: string, layerType: LayerType): string {
  //   return layerName;
  // }

  // private _getRoleColor(roleName: string): number {
  //   if (!this.roleColorMap.has(roleName)) {
  //     this.roleColorMap.set(roleName, colorPalette[this.nextColorIndex % colorPalette.length]);
  //     this.nextColorIndex++;
  //   }
  //   return this.roleColorMap.get(roleName)!;
  // }

  refresh() {
    const currentPaneName = this.tabs.currentPane?.name;
    this.tabs.clear();
    this.grids.clear();
    this.tabs.addPane('Summary', () => this.grid.root);
    this.tabs.addPane(`X Outliers X`, () => this.grid.root);

    const allColumns = this.plate.data.columns.toList()
      .filter((c) => c.name !== PLATE_OUTLIER_WELL_NAME);

    for (const col of allColumns)
      this.tabs.addPane(col.name, () => this.createLayerGrid(col.name, false));


    if (currentPaneName)
      this.tabs.currentPane = this.tabs.panes.find((p) => p.name === currentPaneName) ?? this.tabs.getPane('Summary');
    else if (this.tabs.panes.length > 0)
      this.tabs.currentPane = this.tabs.getPane('Summary');

    const t = this.plate.data;
    const nonLayoutCols = t.columns.toList();
    this._colorColumn = nonLayoutCols.find((c) => c.semType === 'Activity' || c.name.toLowerCase() === 'activity') ??
                      nonLayoutCols.find((c) => c.semType === 'Concentration' || c.name.toLowerCase().includes('concentration')) ??
                      nonLayoutCols.find((c) => c.type === DG.TYPE.FLOAT);

    this.grid.dataFrame = DG.DataFrame.create(this.plate.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= this.plate.cols; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    this.grid.invalidate();
  }


  private createLayerGrid(layer: string, isRole: boolean): HTMLElement {
    const df = this.plate.toGridDataFrame(layer);
    const grid = DG.Viewer.grid(df);
    grid.props.showRowHeader = true;
    grid.props.allowEdit = !isRole;
    return grid.root;
  }

  private configureOutlierGrid(grid: DG.Grid, layerName: string): void {
    grid.onCellRender.subscribe((args) => {
      if (args.cell.isTableCell) {
        const dataRow = this.plate._idx(args.cell.gridRow!, args.cell.gridColumn.idx - 1);
        const outlierCol = this.plate.data.col(layerName);
        if (outlierCol?.get(dataRow)) {
          const g = args.g;
          const bounds = args.bounds;
          g.fillStyle = 'rgba(255, 0, 0, 0.3)';
          g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
          g.fillStyle = 'red';
          g.font = '12px Arial';
          g.textAlign = 'center';
          g.textBaseline = 'middle';
          g.fillText('OUT', bounds.midX, bounds.midY);
        }
      }
    });
  }

  renderCell(args: DG.GridCellRenderArgs, summary: boolean = false) {
    const gc = args.cell;
    const g = args.g;
    const x = args.bounds.x;
    const y = args.bounds.y;
    const w = args.bounds.width;
    const h = args.bounds.height;

    if (summary) {
      g.save();
      g.fillStyle = 'white';
      g.fillRect(x, y, w, h);
      g.restore();
    }

    g.textAlign = 'center';
    g.textBaseline = 'middle';
    g.font = `${Math.ceil(Math.min(...[16, w - 1, h - 1]))}px Roboto, Roboto Local`;

    const hasRowHeader = gc.grid.columns.byIndex(0)?.cellType === 'row header';

    if (gc.isColHeader && gc.gridColumn.idx > (hasRowHeader ? 1 : 0)) {
      g.fillStyle = 'grey';
      g.fillText('' + (gc.gridColumn.idx - (hasRowHeader ? 1 : 0)), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (gc.isColHeader) {
      args.preventDefault();
      return;
    } else if ((gc.gridColumn.name == '0' || gc.gridColumn.idx == 0) && gc.gridRow !== null && gc.gridRow >= 0) {
      g.fillStyle = 'grey';
      const prefix = gc.gridRow > 25 ? 'A' : '';
      g.fillText(prefix + String.fromCharCode(65 + gc.gridRow % 26), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (summary && h > 0 && gc.gridRow !== null && gc.gridRow >= 0 && gc.gridColumn.idx > 0) {
      const dataRow = this._plate._idx(gc.gridRow, gc.gridColumn.idx - 1);
      const isHovered = this.hoveredCell && this.hoveredCell.row === gc.gridRow && this.hoveredCell.col === gc.gridColumn.idx - 1;

      g.beginPath();
      const r = Math.min(h / 2, w / 2) * 0.8;
      g.ellipse(x + w / 2, y + h / 2, r, r, 0, 0, 2 * Math.PI);

      if (this._colorColumn && !this._colorColumn.isNone(dataRow)) {
        const color = this.getColor(dataRow);
        g.fillStyle = DG.Color.toHtml(color);
      } else {
        g.fillStyle = DG.Color.toHtml(DG.Color.lightGray);
      }
      g.fill();

      const isSelected = this.plate.data.selection.get(dataRow) && this.interactionMode === 'default';
      if (isHovered || isSelected) {
        g.shadowColor = 'rgba(40, 255, 140, 0.9)';
        g.shadowBlur = 10;
        g.strokeStyle = 'rgba(40, 255, 140, 0.9)';
        g.lineWidth = isHovered ? 4 : 3;
      } else {
        g.strokeStyle = 'grey';
        g.lineWidth = 1;
      }
      g.stroke();
      g.shadowBlur = 0;

      const outlierCol = this.plate.data.col(PLATE_OUTLIER_WELL_NAME);
      if (outlierCol?.get(dataRow)) {
        g.strokeStyle = 'rgba(255, 0, 0, 0.8)';
        g.lineWidth = 3;
        g.lineCap = 'round';
        const crossSize = r * 0.8;
        const centerX = x + w / 2;
        const centerY = y + h / 2;
        g.beginPath();
        g.moveTo(centerX - crossSize, centerY - crossSize);
        g.lineTo(centerX + crossSize, centerY + crossSize);
        g.moveTo(centerX + crossSize, centerY - crossSize);
        g.lineTo(centerX - crossSize, centerY + crossSize);
        g.stroke();
      }
    }

    if (summary)
      args.preventDefault();
  }

  getColor(dataRow: number) {
    if (!this._colorColumn) return DG.Color.white;
    const isLog = this._colorColumn.name?.toLowerCase()?.includes('conc');
    const val = this._colorColumn.get(dataRow)!;
    const stats = this._colorColumn.stats;
    const reducedVal = isLog ? safeLog(val) : val;
    const min = isLog ? safeLog(stats.min) : stats.min;
    const max = isLog ? safeLog(stats.max) : stats.max;
    return DG.Color.scaleColor(reducedVal, min, max, undefined, colorScheme);
  }

  syncGrids() {
    for (const grid of this.grids.values())
      grid.props.allowEdit = this.editable;
  }

  detach() {
    this.outlierSubscription?.unsubscribe();
    this._onDestroy.next();
    this._onDestroy.complete();
    super.detach();
  }
}
