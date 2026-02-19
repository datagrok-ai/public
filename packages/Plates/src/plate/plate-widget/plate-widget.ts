/* eslint-disable prefer-const */
/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {mapFromRow, toExcelPosition} from './../utils';
import {Plate, PLATE_OUTLIER_WELL_NAME} from './../plate';
import {fromEvent, Subject, Subscription} from 'rxjs';
import {filter, takeUntil} from 'rxjs/operators';
import './plate-widget.css';
import {IPlateWellValidator, plateWellValidators} from './../plate-well-validators';

const colorScheme = [DG.Color.white, DG.Color.gray];

export interface IAnalysisWidgetCoordinator {
  onPlateDataChanged(changeType: 'outlier' | 'data' | 'layer', details: any): void;
  refreshAnalysisView(): void;
}

export interface WellClickEvent {
  row: number;
  col: number;
  dataIndex: number;
}

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
  private _positiveMinCache: number | null = null;
  private _positiveMinCacheCol: DG.Column | undefined = undefined;
  private _onDestroy = new Subject<void>();
  private analysisCoordinators: IAnalysisWidgetCoordinator[] = [];
  private outlierSubscription?: Subscription;

  private wellClickSubject = new Subject<WellClickEvent>();
  public get onWellClick() { return this.wellClickSubject.asObservable(); }

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

    this.refresh();
  }

  public hitTest(canvasX: number, canvasY: number): { row: number, col: number } | null {
    const gridCell = this.grid.hitTest(canvasX, canvasY);

    if (!gridCell || gridCell.gridRow == null || gridCell.gridColumn?.idx == null || gridCell.gridColumn.idx === 0)
      return null;


    const gridRow = gridCell.gridRow;
    const gridCol = gridCell.gridColumn.idx;

    const bounds = gridCell.bounds;
    const centerX = bounds.midX;
    const centerY = bounds.midY;
    const radius = Math.min(bounds.height / 2, bounds.width / 2) * 0.8;
    const distance = Math.sqrt(Math.pow(canvasX - centerX, 2) + Math.pow(canvasY - centerY, 2));

    if (distance <= radius) {
      return {
        row: gridRow,
        col: gridCol - 1
      };
    }

    return null;
  }
  private isPointInWell(x: number, y: number, gridRow: number, gridCol: number): boolean {
    if (gridRow < 0 || gridCol <= 0)
      return false;


    const column = this.grid.columns.byIndex(gridCol);
    if (!column)
      return false;


    const cell = this.grid.cell(column.name, gridRow);
    if (!cell)
      return false;


    const bounds = cell.bounds;

    const centerX = bounds.midX;
    const centerY = bounds.midY;
    const radius = Math.min(bounds.height / 2, bounds.width / 2) * 0.8;
    const distance = Math.sqrt(Math.pow(x - centerX, 2) + Math.pow(y - centerY, 2));


    return distance <= radius;
  }

  public addAnalysisTab(name: string, content: HTMLElement, makeCurrent = true): void {
    const pane = this.tabs.addPane(name, () => content);
    if (makeCurrent) this.tabs.currentPane = pane;
  }

  constructor() {
    super(ui.div([], 'assay-plates--plate-widget'));
    this.tabs.root.classList.add('assay-plates--plate-widget__tabs');
    this.tabsContainer.appendChild(this.tabs.root);

    const mainContainer = ui.divV([this.tabsContainer], 'assay-plates--plate-widget__main-content');
    this.root.appendChild(mainContainer);

    this.setupGrid();
  }
  public fireWellClick(row: number, col: number): void {
    const dataIndex = this.plate._idx(row, col);
    this.wellClickSubject.next({row, col, dataIndex});
  }


  private setupGrid() {
    this.grid.props.allowColReordering=false;
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
    fromEvent<MouseEvent>(this.grid.canvas, 'click')
      .pipe(takeUntil(this._onDestroy))
      .subscribe((event: MouseEvent) => {
        const gridCell = this.grid.hitTest(event.offsetX, event.offsetY);

        if (!gridCell || !gridCell.isTableCell || gridCell.gridColumn.idx === 0)
          return;

        const gridRow = gridCell.gridRow;
        const gridColIdx = gridCell.gridColumn.idx; // 1-based index

        if (this.isPointInWell(event.offsetX, event.offsetY, gridRow, gridColIdx)) {
          const wellCol = gridColIdx - 1; /// Convert to 0-based
          const dataIndex = this.plate._idx(gridRow, wellCol);

          this.wellClickSubject.next({
            row: gridRow,
            col: wellCol,
            dataIndex: dataIndex,
          });
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

    // ui.empty(pw.root);

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
    grid.props.showHeatmapScrollbars = false;
    grid.props.allowColReordering = false;
    grid.props.allowRowReordering = false;
    grid.props.allowSorting = false;
    grid.props.allowEdit = !isSummary;
    grid.props.showRowGridlines = false;
    grid.props.showColumnGridlines = false;
    grid.props.heatmapColors = false;
    grid.props.colHeaderHeight = 25;
    grid.props.allowRowSelection = false;
    grid.props.allowColSelection = false;
    grid.props.allowBlockSelection = false;
    grid.props.showCurrentRowIndicator = false;
    grid.props.showMouseOverRowIndicator = false;
    grid.props.showCurrentCellOutline = false;
    grid.props.allowColHeaderResizing = false;
    grid.props.allowColResizing = false;
    grid.props.allowRowResizing = false;
    grid.props.selectedRowsColor = 0x00000000;
    grid.props.mouseOverRowColor = 0x00000000;
    grid.props.currentRowColor = 0x00000000;
    grid.onCellRender.subscribe((args) => this.renderCell(args, isSummary));
    this.setupHoverEvents(grid);
  }

  refresh() {
    this._positiveMinCache = null;
    this._positiveMinCacheCol = undefined;
    const currentPaneName = this.tabs.currentPane?.name;
    this.tabs.clear();
    this.grids.clear();
    this.tabs.addPane('Summary', () => {
      // this.grid.invalidate();
      return this.grid.root;
    });
    // this.tabs.addPane(`X Outliers X`, () => this.grid.root);
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

      const isSelected = this.plate.data.selection.get(dataRow);
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

      // Render outlier marker
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

  /** Returns the smallest positive value in the color column (cached). */
  private getPositiveMin(): number | null {
    if (this._positiveMinCacheCol === this._colorColumn && this._positiveMinCache !== null)
      return this._positiveMinCache;

    if (!this._colorColumn) return null;

    let min = Infinity;
    for (let i = 0; i < this._colorColumn.length; i++) {
      if (!this._colorColumn.isNone(i)) {
        const v = this._colorColumn.get(i);
        if (v > 0 && v < min) min = v;
      }
    }
    this._positiveMinCache = min === Infinity ? null : min;
    this._positiveMinCacheCol = this._colorColumn;
    return this._positiveMinCache;
  }

  getColor(dataRow: number) {
    if (!this._colorColumn) return DG.Color.white;
    const isLog = this._colorColumn.name?.toLowerCase()?.includes('conc');
    const val = this._colorColumn.get(dataRow)!;
    const stats = this._colorColumn.stats;

    if (isLog) {
      // Non-positive concentrations (e.g. control wells) get the base color
      if (val <= 0) return colorScheme[0];

      const reducedVal = Math.log10(val);
      const positiveMin = this.getPositiveMin();
      const logMin = positiveMin !== null ? Math.log10(positiveMin) : reducedVal;
      const logMax = stats.max > 0 ? Math.log10(stats.max) : reducedVal;

      return DG.Color.scaleColor(reducedVal, Math.min(logMin, logMax), Math.max(logMin, logMax), undefined, colorScheme);
    }

    return DG.Color.scaleColor(val, stats.min, stats.max, undefined, colorScheme);
  }

  syncGrids() {
    for (const grid of this.grids.values())
      grid.props.allowEdit = this.editable;
  }

  detach() {
    this.outlierSubscription?.unsubscribe();
    this._onDestroy.next();
    this._onDestroy.complete();
    this.wellClickSubject.complete();
    super.detach();
  }
}
