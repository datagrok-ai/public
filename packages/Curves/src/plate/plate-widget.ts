/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {MARKER_TYPE, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {div} from 'datagrok-api/ui';
import {mapFromRow, safeLog, toExcelPosition} from './utils';
import {LayerType, Plate, PLATE_OUTLIER_WELL_NAME} from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {IPlateWellValidator, plateWellValidators} from './plate-well-validators';
import {fromEvent, Subject, Subscription} from 'rxjs';
import {debounceTime, filter, takeUntil, take} from 'rxjs/operators';
import './plate-widget.css';

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
    plateStatistics?: {[key: string]: (plate: Plate) => string | number}
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
// NEW: Constant to define the minimum pixel distance to initiate a drag action
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
  tabsContainer: HTMLElement = ui.divH([], 'plate-widget__tabs-container');
  _editable: boolean = false;
  mapFromRowFunc: (row: DG.Row) => Record<string, any> = mapFromRow;
  grids: Map<string, DG.Grid> = new Map();
  wellValidators: IPlateWellValidator[] = plateWellValidators;
  wellValidationErrors: Map<string, string[]> = new Map();

  // MODIFIED: State variables for the new interaction logic
  private _isDragging: boolean = false;
  private _selectionRect: DG.Rect | null = null;
  private _canvas: HTMLCanvasElement | null = null;
  private _onDestroy = new Subject<void>();
  private analysisCoordinators: IAnalysisWidgetCoordinator[] = [];
  private outlierSubscription?: Subscription;
  private interactionMode: InteractionMode = 'default';
  private _dragStartPoint: DG.Point | null = null;


  roleSummaryDiv: HTMLElement = ui.divH([], 'plate-widget__role-summary');
  get editable() { return this._editable; }
  set editable(x: boolean) {
    this._editable = x;
    this.syncGrids();
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

  // NEW: Hit-testing function to find a well from canvas coordinates
  private hitTest(canvasX: number, canvasY: number): { row: number, col: number } | null {
    for (let row = 0; row < this.plate.rows; row++) {
      for (let col = 0; col < this.plate.cols; col++) {
        // Check if the point is inside the circular area of the well
        if (this.isPointInWell(canvasX, canvasY, row, col + 1))
          return {row, col};
      }
    }
    return null;
  }


  constructor() {
    super(ui.div([], 'plate-widget'));

    this.tabs.root.classList.add('plate-widget__tabs');
    this.tabsContainer.appendChild(this.tabs.root);

    this.tabs.onTabChanged.subscribe(() => {
      const isOutlierMode = this.tabs.currentPane.name.includes('Outliers');
      this.interactionMode = isOutlierMode ? 'outlier' : 'default';
      this.grid.root.classList.toggle('outlier-mode', isOutlierMode);
      this.plate.data.selection.setAll(false, true);
      this.grid.invalidate();
    });

    const mainContainer = ui.divV([this.tabsContainer, this.roleSummaryDiv]);
    mainContainer.style.flexGrow = '1';
    mainContainer.style.width = '100%';
    this.root.appendChild(mainContainer);


    this.grid.props.allowRowSelection = false;
    this.grid.props.allowColSelection = false;
    this.grid.props.allowBlockSelection = false;
    this.grid.props.showCurrentRowIndicator = false;
    this.grid.props.showMouseOverRowIndicator = false;
    this.grid.props.showCurrentCellOutline = false;
    this.grid.props.showHeatmapScrollbars = false;
    this.grid.props.colHeaderHeight = 30;
    this.grid.props.allowColHeaderResizing = false;
    this.grid.props.allowColResizing = false;
    this.grid.props.allowRowResizing = false;

    this.grid.onCellRender.subscribe((args) => this.renderCell(args, true));
    this.grid.onCellRendered.subscribe((args) => {
      const cell = args.cell;
      if (cell.gridRow >= 0 && cell.gridColumn.idx > 0) {
        try {
          const pos = toExcelPosition(cell.gridRow, cell.gridColumn.idx - 1);
          const errors = this.wellValidationErrors.get(pos);
          if (errors && errors.length > 0)
            DG.Paint.marker(args.g, DG.MARKER_TYPE.CROSS_BORDER, cell.bounds.midX, cell.bounds.midY, DG.Color.red, 10);
        } catch (e) {
          console.error(e);
        }
      }
    });

    let hoveredCell: {row: number, col: number} | null = null;

    this.grid.onCellMouseEnter.subscribe((gc: DG.GridCell) => {
      if (gc.isTableCell) {
        hoveredCell = {row: gc.gridRow, col: gc.gridColumn.idx - 1};
        this.grid.invalidate();
      }
    });

    this.grid.onCellMouseLeave.subscribe((gc: DG.GridCell) => {
      hoveredCell = null;
      this.grid.invalidate();
    });

    // MODIFIED: This is now a fallback and can be removed, as our new system handles clicks.
    // We'll leave it for now in case other parts of the system rely on it.
    this.grid.onCellClick.subscribe((gc: DG.GridCell) => {
      if (!gc.isTableCell) return;
      if (this.interactionMode === 'outlier' && this.editable) {
        const row = gc.gridRow;
        const col = gc.gridColumn.idx - 1;
        const currentState = this.plate.isOutlier(row, col);
        this.plate.markOutlierWithSource(row, col, !currentState, 'summary-click-toggle');
        this.grid.invalidate();
      }
    });

    (this as any)._hoveredCell = () => hoveredCell;


    ui.tools.waitForElementInDom(this.grid.root).then(() => {
      this._canvas = this.grid.overlay;
      this.initSelectionEvents();
    });
  }


  registerAnalysisCoordinator(coordinator: IAnalysisWidgetCoordinator): void {
    this.analysisCoordinators.push(coordinator);
  }

  private notifyAnalysisCoordinators(changeType: 'outlier' | 'data' | 'layer', details: any): void {
    this.analysisCoordinators.forEach((coord) => coord.onPlateDataChanged(changeType, details));
  }


  private refreshTabs(): void {
    const currentPaneName = this.tabs.currentPane?.name;

    this.tabs.clear();
    this.grids.clear();

    this.tabs.addPane('Summary', () => this.grid.root);

    const layerInfo: Record<LayerType, {icon: string, layers: string[]}> = {
      [LayerType.ORIGINAL]: {icon: '', layers: this.plate.getLayersByType(LayerType.ORIGINAL)},
      [LayerType.LAYOUT]: {icon: '️(Layout)', layers: this.plate.getLayersByType(LayerType.LAYOUT)},
      [LayerType.DERIVED]: {icon: '(Derived)', layers: this.plate.getLayersByType(LayerType.DERIVED)},
      [LayerType.OUTLIER]: {icon: '🚫', layers: []},
    };

    const allRegisteredLayers = new Set<string>();
    for (const type of [LayerType.ORIGINAL, LayerType.LAYOUT, LayerType.DERIVED]) {
      for (const layerName of layerInfo[type].layers) {
        allRegisteredLayers.add(layerName);
        const paneName = `${layerInfo[type].icon} ${this.getDisplayName(layerName, type)}`;
        this.tabs.addPane(paneName, () => this.createLayerGrid(layerName));
      }
    }

    this.tabs.addPane(`🚫 Outliers`, () => this.grid.root);

    for (const layerName of this.plate.getLayerNames()) {
      if (!allRegisteredLayers.has(layerName) && layerName !== PLATE_OUTLIER_WELL_NAME)
        this.tabs.addPane(`❔ ${layerName}`, () => this.createLayerGrid(layerName));
    }

    if (currentPaneName) {
      const restoredPane = this.tabs.panes.find((p) => p.name === currentPaneName);
      if (restoredPane)
        this.tabs.currentPane = restoredPane;
    }
  }

  private layerExists(layerName: string): boolean {
    return this.tabs.panes.some((pane) =>
      pane.name.includes(layerName) ||
            (layerName === PLATE_OUTLIER_WELL_NAME && pane.name.includes('Outliers'))
    );
  }

  // MODIFIED: Complete overhaul of the selection logic.
  private initSelectionEvents() {
    const interactionElement = this.grid.overlay;
    if (!interactionElement) return;

    fromEvent<MouseEvent>(interactionElement, 'mousedown')
      .pipe(takeUntil(this._onDestroy))
      .subscribe((e: MouseEvent) => {
        if (!this.editable || e.button !== 0) return;

        this._dragStartPoint = new DG.Point(e.offsetX, e.offsetY);

        const mouseMoveStream = fromEvent<MouseEvent>(document, 'mousemove').pipe(takeUntil(this._onDestroy));
        const mouseUpStream = fromEvent<MouseEvent>(document, 'mouseup').pipe(takeUntil(this._onDestroy));

        mouseMoveStream
          .pipe(takeUntil(mouseUpStream))
          .subscribe((move_e: MouseEvent) => {
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

              // MODIFIED: Only draw the rectangle, do not update the selection state here.
              this.drawSelectionRect();
            }
          });

        mouseUpStream.pipe(take(1)).subscribe((up_e: MouseEvent) => {
          if (this._isDragging) {
            // MODIFIED: Finalize the selection state here, on mouseup.
            this.finalizeDragSelection();
            const selection = this.plate.data.selection;
            if (selection.trueCount > 0) {
              if (this.interactionMode === 'default') {
                this.showRoleAssignmentPopup();
              } else if (this.interactionMode === 'outlier') {
                for (const i of selection.getSelectedIndexes()) {
                  const [row, col] = this.plate.rowIndexToExcel(i);
                  this.plate.markOutlierWithSource(row, col, true, 'summary-drag-set');
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
                this.plate.markOutlierWithSource(well.row, well.col, !currentState, 'summary-click-toggle');
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
  updateRoleSummary() {
    ui.empty(this.roleSummaryDiv);
    const roleColumn = this.plate.data.col('Role');
    if (!roleColumn || roleColumn.stats.valueCount === 0) return;

    const countsDf = this.plate.data.groupBy(['Role']).count().aggregate();
    const roleColumnColors = roleColumn.meta.colors;
    const originalRolesList = roleColumn.toList();

    for (const row of countsDf.rows) {
      const role = row.Role;
      const count = row.count;
      const firstIndex = originalRolesList.indexOf(role);
      const color = (firstIndex !== -1) ?
        DG.Color.toHtml(roleColumnColors.getColor(firstIndex)) :
        DG.Color.toHtml(DG.Color.lightGray);

      const legendItem = ui.divH([
        ui.div('', {style: {width: '12px', height: '12px', borderRadius: '3px', backgroundColor: color, border: '1px solid var(--grey-3)'}}),
        ui.divText(`${role} (${count} wells)`),
      ], 'role-summary__item');
      this.roleSummaryDiv.appendChild(legendItem);
    }
  }

  private finalizeDragSelection() {
    if (!this._selectionRect || !this.grid || this.grid.columns.length === 0)
      return;

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
          if (cell) {
            const cellBounds = cell.bounds;
            if (normalizedRect.contains(cellBounds.midX, cellBounds.midY)) {
              const dataIndex = this.plate._idx(row, col);
              selection.set(dataIndex, true, false);
            }
          }
        }
      }
    }
    selection.fireChanged();
  }

  private showRoleAssignmentPopup(): void {
    const selection = this.plate.data.selection;
    if (selection.trueCount === 0) return;

    const roles = ['Control', 'Buffer', 'Assay Reagent', 'Sample'];
    const roleInput = ui.input.multiChoice<string>('Roles', {items: roles, value: []});

    const popupContent = ui.divV([
      ui.h3(`${selection.trueCount} wells selected`),
      roleInput,
      ui.button('Assign', () => {
        const selectedRoles = roleInput.value;
        if (!selectedRoles || selectedRoles.length === 0) return;

        let roleCol = this.plate.data.col('Role');
        if (roleCol === null) {
          roleCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Role', this.plate.data.rowCount);
          this.plate.data.columns.add(roleCol);
          this.plate.registerLayer('Role', LayerType.LAYOUT, 'user-assignment');
        }

        const selectedIndexes = selection.getSelectedIndexes();
        for (const i of selectedIndexes)
          (roleCol as DG.Column<string>).set(i, selectedRoles.join(','));

        this._colorColumn = roleCol;
        if (this._colorColumn.isCategorical && this._colorColumn.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
          this._colorColumn.meta.colors.setCategorical();


        this.updateRoleSummary();
        this.refresh();
        selection.setAll(false, true); // Clear selection after assigning
        this.grid.invalidate();

        popup.remove();
        anchorDiv.remove();
      })
    ], 'd4-menu-item-container');

    popupContent.style.padding = '10px';

    const selectedIndices = selection.getSelectedIndexes();
    let minRow = this.plate.rows; let maxRow = -1; let minCol = this.plate.cols; let maxCol = -1;

    for (const idx of selectedIndices) {
      const row = Math.floor(idx / this.plate.cols);
      const col = idx % this.plate.cols;
      if (row < minRow) minRow = row;
      if (row > maxRow) maxRow = row;
      if (col < minCol) minCol = col;
      if (col > maxCol) maxCol = col;
    }

    const firstCell = this.grid.cell(this.grid.columns.byIndex(minCol + 1)!.name, minRow);
    const lastCell = this.grid.cell(this.grid.columns.byIndex(maxCol + 1)!.name, maxRow);
    const selectionBounds = firstCell.bounds.union(lastCell.bounds);

    const canvasBounds = this.grid.overlay.getBoundingClientRect();

    const anchorDiv = ui.div('', {
      style: {
        position: 'absolute',
        left: `${canvasBounds.left + selectionBounds.x + selectionBounds.width / 2}px`,
        top: `${canvasBounds.top + selectionBounds.y}px`,
        width: '0px',
        height: '0px',
      }
    });
    document.body.appendChild(anchorDiv);

    const popup = ui.showPopup(popupContent, anchorDiv);

    const closePopup = () => {
      if (document.body.contains(popup)) popup.remove();
      if (document.body.contains(anchorDiv)) anchorDiv.remove();
      document.removeEventListener('mousedown', onMouseDown, true);
    };

    const onMouseDown = (e: MouseEvent) => {
      if (!popup.contains(e.target as Node)) closePopup();
    };

    setTimeout(() => document.addEventListener('mousedown', onMouseDown, true), 0);
  }

  private drawSelectionRect() {
    if (!this._canvas || !this._selectionRect) return;
    const g = this._canvas.getContext('2d')!;
    g.clearRect(0, 0, this._canvas.width, this._canvas.height);
    g.strokeStyle = 'rgba(0, 128, 255, 0.7)';
    g.fillStyle = 'rgba(0, 128, 255, 0.2)';
    g.lineWidth = 1;
    g.strokeRect(this._selectionRect.x, this._selectionRect.y, this._selectionRect.width, this._selectionRect.height);
    g.fillRect(this._selectionRect.x, this._selectionRect.y, this._selectionRect.width, this._selectionRect.height);
  }

  private clearSelectionRect() {
    if (!this._canvas) return;
    const g = this._canvas.getContext('2d')!;
    g.clearRect(0, 0, this._canvas.width, this._canvas.height);
  }

  private updateSelection() {
    if (!this._selectionRect || !this.grid || this.grid.columns.length === 0)
      return;

    const selection = this.plate.data.selection;
    // MODIFIED: The line `selection.setAll(false, false)` is REMOVED from here.

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
          if (cell) {
            const cellBounds = cell.bounds;
            if (normalizedRect.contains(cellBounds.midX, cellBounds.midY)) {
              const dataIndex = this.plate._idx(row, col);
              // This now ADDs to the selection.
              selection.set(dataIndex, true, false);
            }
          }
        }
      }
    }
    selection.fireChanged();
    this.grid.invalidate();
  }


  static fromPlate(plate: Plate, useSimpleView: boolean = false) {
    if (useSimpleView) {
      // Simple view without the detailed layout
      const pw = new PlateWidget();
      pw.plate = plate;
      return pw;
    }
    return PlateWidget.detailedView(plate);
  }

  static detailedView(plate: Plate): PlateWidget {
    const pw = new PlateWidget();
    pw.plate = plate;

    pw.roleSummaryDiv.remove();

    pw.detailsDiv = ui.divV([], 'plate-widget__details');
    pw.wellDetailsDiv = ui.div();
    pw.detailsDiv.appendChild(pw.wellDetailsDiv);

    const gridAndSummaryWrapper = ui.divV([
      pw.tabs.root,
      pw.roleSummaryDiv,
    ], {style: {flexGrow: '1', display: 'flex', flexDirection: 'column'}});

    const mainContainer = ui.divH([
      gridAndSummaryWrapper,
      pw.detailsDiv,
    ], 'plate-widget__main-container');

    ui.empty(pw.tabsContainer);
    pw.tabsContainer.appendChild(mainContainer);

    pw.grid.onCurrentCellChanged
      .pipe(filter((gc) => gc.gridRow >= 0 && gc.gridColumn.idx > 0))
      .subscribe((gc) => {
        if (pw.wellDetailsDiv) {
          ui.empty(pw.wellDetailsDiv);
          const map = pw.mapFromRowFunc(pw.plate.data.rows.get(plate._idx(gc.gridRow, gc.gridColumn.idx - 1)));
          pw.wellDetailsDiv.appendChild(ui.tableFromMap(map));
        }
      });

    return pw;
  }

  initPlateGrid(grid: DG.Grid, isSummary: boolean = false) {
    // Core display settings
    grid.props.showHeatmapScrollbars = false;
    grid.props.allowColReordering = false;
    grid.props.allowRowReordering = false;
    grid.props.allowSorting = false;
    grid.props.allowEdit = !isSummary;
    grid.props.showRowGridlines = false;
    grid.props.showColumnGridlines = false;
    grid.props.heatmapColors = false;
    grid.props.colHeaderHeight = 25;

    // Disable selection visuals - using actual properties from IGridSettings
    grid.props.allowRowSelection = false;
    grid.props.allowColSelection = false;
    grid.props.allowBlockSelection = false;
    grid.props.showCurrentRowIndicator = false;
    grid.props.showMouseOverRowIndicator = false;
    grid.props.showCurrentCellOutline = false; // This should help with the black box

    // Disable resizing - correct properties
    grid.props.allowColHeaderResizing = false;
    grid.props.allowColResizing = false;
    grid.props.allowRowResizing = false;

    // Make selection colors transparent (using 0 for transparent)
    grid.props.selectedRowsColor = 0x00000000; // Fully transparent
    grid.props.mouseOverRowColor = 0x00000000; // Fully transparent
    grid.props.currentRowColor = 0x00000000; // Fully transparent

    grid.onCellRender.subscribe((args) => this.renderCell(args, isSummary));
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

    this.refresh();
  }


  private getDisplayName(layerName: string, layerType: LayerType): string {
    return layerName;
  }


  refresh() {
    this.tabs.clear();
    this.grids.clear();
    this.tabs.addPane('Summary', () => this.grid.root);
    this.tabs.addPane(`🚫 Outliers`, () => this.grid.root);

    const layerInfo: Record<LayerType, {icon: string, layers: string[]}> = {
      [LayerType.ORIGINAL]: {icon: '', layers: this.plate.getLayersByType(LayerType.ORIGINAL)},
      [LayerType.LAYOUT]: {icon: '️(Layout)', layers: this.plate.getLayersByType(LayerType.LAYOUT)},
      [LayerType.DERIVED]: {icon: '(Derived)', layers: this.plate.getLayersByType(LayerType.DERIVED)},
      [LayerType.OUTLIER]: {icon: '🚫', layers: []},
    };

    const allRegisteredLayers = new Set<string>();

    for (const type of [LayerType.ORIGINAL, LayerType.LAYOUT, LayerType.DERIVED]) {
      for (const layerName of layerInfo[type].layers) {
        allRegisteredLayers.add(layerName);
        const paneName = `${layerInfo[type].icon} ${this.getDisplayName(layerName, type)}`;
        this.tabs.addPane(paneName, () => this.createLayerGrid(layerName));
      }
    }

    for (const layerName of this.plate.getLayerNames()) {
      if (!allRegisteredLayers.has(layerName) && layerName !== PLATE_OUTLIER_WELL_NAME)
        this.tabs.addPane(`❔ ${layerName}`, () => this.createLayerGrid(layerName));
    }

    ui.tools.waitForElementInDom(this.tabs.root).then(() => {
      if (this.tabs.panes.length > 0)
        this.tabs.currentPane = this.tabs.getPane('Summary');
    });


    this.syncGrids();

    const t = this.plate.data;
    this._colorColumn =
            t.columns.firstWhere((c) => c.semType === 'Activity') ??
            t.columns.firstWhere((c) => c.semType === 'Concentration') ??
            t.columns.firstWhere((c) => c.name.toLowerCase() === 'activity') ??
            t.columns.firstWhere((c) => c.name.toLowerCase().includes('concentration')) ??
            t.columns.firstWhere((c) => c.type === DG.TYPE.FLOAT && !['row', 'col'].includes(c.name.toLowerCase()));


    this.grid.dataFrame = DG.DataFrame.create(this.plate.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= this.plate.cols; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    if (this.grid && this.grid.root)
      this.grid.root.style.width = '100%';

    this.grid.invalidate();
  }

  private createLayerGrid(layer: string): HTMLElement {
    const df = this.plate.toGridDataFrame(layer);
    const grid = DG.Viewer.heatMap(df);
    grid.columns.add({gridColumnName: '0', cellType: 'string', index: 1});

    df.onValuesChanged.pipe(debounceTime(1000)).subscribe(() => {
      const p = this.plate;
      for (let i = 0; i < df.rowCount; i++) {
        for (let j = 0; j < df.columns.length - 1; j++) {
          const newValue = df.get(`${j + 1}`, i);
          const plateCol = p.getColumn(layer);
          if (plateCol)
            plateCol.set(p._idx(i, j), newValue);
        }
      }
    });

    this.initPlateGrid(grid, false);
    this.grids.set(layer, grid);
    return grid.root;
  }
  private configureOutlierGrid(grid: DG.Grid, layerName: string): void {
    grid.onCellRender.subscribe((args) => {
      if (args.cell.gridColumn.idx > 0 && args.cell.gridRow >= 0) {
        const dataRow = this.plate._idx(args.cell.gridRow, args.cell.gridColumn.idx - 1);
        const outlierCol = this.plate.data.col(layerName);

        if (outlierCol && outlierCol.get(dataRow)) {
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
    const isColoredByConc = this._colorColumn?.name?.toLowerCase()?.includes('conc');

    const hasRowHeader = gc.grid.columns.byIndex(0)?.cellType === 'row header';

    if (gc.isColHeader && gc.gridColumn.idx > (hasRowHeader ? 1 : 0)) {
      g.fillStyle = 'grey';
      g.fillText('' + (gc.gridColumn.idx - (hasRowHeader ? 1 : 0)), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (gc.isColHeader) {
      args.preventDefault();
      return;
    } else if ((gc.gridColumn.name == '0' || gc.gridColumn.idx == 0) && gc.gridRow >= 0) {
      g.fillStyle = 'grey';
      const prefix = gc.gridRow > 25 ? 'A' : '';
      g.fillText(prefix + String.fromCharCode(65 + gc.gridRow % 26), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (summary && h > 0 && gc.gridRow >= 0 && gc.gridColumn.idx > 0) {
      const dataRow = this._plate._idx(gc.gridRow, gc.gridColumn.idx - 1);

      const hoveredCell = (this as any)._hoveredCell ? (this as any)._hoveredCell() : null;
      const isHovered = hoveredCell &&
                                hoveredCell.row === gc.gridRow &&
                                hoveredCell.col === gc.gridColumn.idx - 1;

      g.beginPath();
      const r = Math.min(h / 2, w / 2) * 0.8;
      g.ellipse(x + w / 2, y + h / 2, r, r, 0, 0, 2 * Math.PI);

      if (this._colorColumn) {
        if (this._colorColumn.isCategorical && this._colorColumn.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
          this._colorColumn.meta.colors.setCategorical();
        const color = this._colorColumn.isNone(dataRow) ? DG.Color.white : this._colorColumn.isNumerical ?
          this.getColor(dataRow, isColoredByConc) :
          this._colorColumn.meta.colors.getColor(dataRow);
        g.fillStyle = DG.Color.toHtml(color);
        g.fill();
      }

      // Draw hover effect
      if (isHovered && this.interactionMode !== 'outlier') {
        g.strokeStyle = 'rgba(0, 128, 255, 0.8)';
        g.lineWidth = 3;
      } else if (this.plate.data.selection.get(dataRow) && this.interactionMode === 'default') {
        // Selection ring
        g.strokeStyle = 'rgba(0, 128, 255, 0.5)';
        g.lineWidth = 2;
      } else {
        g.strokeStyle = 'grey';
        g.lineWidth = 1;
      }
      g.stroke();

      const outlierCol = this.plate.data.col(PLATE_OUTLIER_WELL_NAME);
      if (outlierCol?.get(dataRow))
        DG.Paint.marker(g, MARKER_TYPE.CROSS_X_BORDER, x + w / 2, y + h / 2, DG.Color.red, r * 2);
    }

    if (summary)
      args.preventDefault();
  }

  getColor(dataRow: number, isLog?: boolean) {
    const val = this._colorColumn!.get(dataRow)!;
    const reducedVal = isLog ? safeLog((val - this._colorColumn!.min) * 1e9) : val;
    const min = (isLog ? safeLog(Math.max(this._colorColumn!.min, 1)) : this._colorColumn!.min);
    const max = isLog ? safeLog((this._colorColumn!.max - this._colorColumn!.min) * 1e9) : this._colorColumn!.max;
    return DG.Color.scaleColor(reducedVal, min, max, undefined, colorScheme);
  }

  syncGrids() {
    for (const grid of this.grids.values())
      grid.props.allowEdit = this.editable;
  }

  detach() {
    if (this.outlierSubscription)
      this.outlierSubscription.unsubscribe();

    this._onDestroy.next();
    this._onDestroy.complete();
    super.detach();
  }
}
