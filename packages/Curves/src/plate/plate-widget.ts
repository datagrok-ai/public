/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {MARKER_TYPE, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {div} from 'datagrok-api/ui';
import {mapFromRow, safeLog, toExcelPosition} from './utils';
import {Plate, PLATE_OUTLIER_WELL_NAME} from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {IPlateWellValidator, plateWellValidators} from './plate-well-validators';
import {fromEvent, Subject} from 'rxjs';
import {debounceTime, filter, takeUntil, take} from 'rxjs/operators';

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

  private _isDragging: boolean = false;
  private _selectionRect: DG.Rect | null = null;
  private _canvas: HTMLCanvasElement | null = null;
  private _onDestroy = new Subject<void>();

  // ADDED: To track which role tab is currently active for filtering
  private _activeRoleFilter: string | null = null;

  roleSummaryDiv: HTMLElement = ui.divH([], 'plate-widget__role-summary');

  get editable() { return this._editable; }
  set editable(x: boolean) {
    this._editable = x;
    this.syncGrids();
  }

  constructor() {
    super(ui.div([], 'plate-widget'));

    this.tabs.root.classList.add('plate-widget__tabs');
    this.tabsContainer.appendChild(this.tabs.root);
    const mainContainer = ui.divV([this.tabsContainer, this.roleSummaryDiv]);
    mainContainer.style.flexGrow = '1';
    mainContainer.style.width = '100%';
    this.root.appendChild(mainContainer);

    this.grid.props.allowRowSelection = false;
    this.grid.props.allowColSelection = false;
    this.grid.props.allowBlockSelection = false;
    this.grid.props.showCurrentRowIndicator = false;
    this.grid.props.showMouseOverRowIndicator = false;

    this.grid.props.showHeatmapScrollbars = false;
    this.grid.props.colHeaderHeight = 30;

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

    ui.tools.waitForElementInDom(this.grid.root).then(() => {
      this._canvas = this.grid.overlay;
      this.initSelectionEvents();
    });
  }

  private initSelectionEvents() {
    const interactionElement = this.grid.overlay;
    if (!interactionElement) return;

    fromEvent<MouseEvent>(interactionElement, 'mousedown')
      .pipe(takeUntil(this._onDestroy))
      .subscribe((e: MouseEvent) => {
        if (!this.editable) return;
        e.preventDefault();
        e.stopPropagation();

        this._isDragging = true;
        this._selectionRect = new DG.Rect(e.offsetX, e.offsetY, 0, 0);

        const mouseMoveStream = fromEvent<MouseEvent>(document, 'mousemove');
        const mouseUpStream = fromEvent<MouseEvent>(document, 'mouseup');

        mouseMoveStream
          .pipe(takeUntil(mouseUpStream))
          .subscribe((move_e: MouseEvent) => {
            if (!this._isDragging || !this._selectionRect) return;
            const canvasBounds = interactionElement.getBoundingClientRect();
            const currentX = move_e.clientX - canvasBounds.left;
            const currentY = move_e.clientY - canvasBounds.top;
            this._selectionRect.width = currentX - this._selectionRect.x;
            this._selectionRect.height = currentY - this._selectionRect.y;
            this.drawSelectionRect();
            this.updateSelection();
          });

        mouseUpStream.pipe(take(1)).subscribe(() => {
          this._isDragging = false;
          this.clearSelectionRect();

          if (this.plate.data.selection.trueCount > 0)
            this.showRoleAssignmentPopup();

          this._selectionRect = null;
        });
      });
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
        }

        const selectedIndexes = selection.getSelectedIndexes();

        for (const i of selectedIndexes)
          (roleCol as DG.Column<string>).set(i, selectedRoles.join(','));

        this._colorColumn = roleCol;
        if (this._colorColumn.isCategorical && this._colorColumn.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
          this._colorColumn.meta.colors.setCategorical();

        this.grid.invalidate();
        this.updateRoleSummary();
        selection.setAll(false, true);
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
    selection.setAll(false, false);

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
    this.grid.invalidate();
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
    grid.props.showHeatmapScrollbars = false;
    grid.props.allowColReordering = false;
    grid.props.allowRowReordering = false;
    grid.props.allowSorting = false;
    grid.props.allowEdit = !isSummary;
    grid.props.showRowGridlines = true;
    grid.props.showColumnGridlines = true;
    grid.props.heatmapColors = false;
    grid.props.colHeaderHeight = 25;
    grid.props.allowRowSelection = false;
    grid.props.allowColSelection = false;
    grid.props.allowBlockSelection = false;
    grid.props.showCurrentRowIndicator = false;
    grid.props.showMouseOverRowIndicator = false;

    grid.onCellRender.subscribe((args) => this.renderCell(args, isSummary));
  }

  get plate(): Plate { return this._plate; }
  set plate(p: Plate) {
    console.log(`[DEBUG] 4. PlateWidget: Plate setter called with plate barcode: ${p.barcode}. Refreshing widget.`);
    this._plate = p;
    this.refresh();
  }

  refresh() {
    this.tabs.clear();
    this.grids.clear(); // Clear old grid references

    // --- Main Summary Tab ---
    this.tabs.addPane('Summary', () => this.grid.root);

    // --- Role-based Tabs ---
    const roleCol = this.plate.data.col('Role');
    if (roleCol) {
      this._colorColumn = roleCol;
      // Ensure the 'Role' column is treated as categorical for coloring
      if (roleCol.isCategorical && roleCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
        roleCol.meta.colors.setCategorical();

      for (const role of roleCol.categories) {
        // Add a tab for each role. The content can be empty because the main grid is always visible.
        if (role) this.tabs.addPane(role, () => ui.div());
      }
    }

    // --- Tab Change Logic ---
    this.tabs.onTabChanged.subscribe((tab) => {
      const roleCategories = roleCol?.categories ?? [];
      if (roleCategories.includes(tab.name)) {
        this._activeRoleFilter = tab.name;
      } else {
        // For 'Summary' or any other tab, clear the filter
        this._activeRoleFilter = null;
      }
      this.grid.invalidate(); // Re-render the grid with the new filter
    });


    // Ensure the tabs container takes full width
    if (this.tabsContainer) {
      this.tabsContainer.style.width = '100%';
      this.tabsContainer.style.display = 'flex';
      this.tabsContainer.style.flexDirection = 'column';
    }

    // Set Summary as the default tab
    ui.tools.waitForElementInDom(this.tabs.root).then(() => {
      if (this.tabs.panes.length > 0)
        this.tabs.currentPane = this.tabs.getPane('Summary');
    });

    this.syncGrids();

    this.grid.dataFrame = DG.DataFrame.create(this.plate.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= this.plate.cols; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    if (this.grid && this.grid.root)
      this.grid.root.style.width = '100%';

    this.grid.invalidate();
    this.updateRoleSummary();
  }


  renderCell(args:DG.GridCellRenderArgs, summary: boolean = false) {
    const gc = args.cell;
    args.g.fillStyle = 'grey';
    args.g.strokeStyle = 'grey';
    args.g.lineWidth = 1;
    const g = args.g;
    const x = args.bounds.x; const y = args.bounds.y; const w = args.bounds.width; const h = args.bounds.height;
    g.textAlign = 'center';
    g.textBaseline = 'middle';
    g.font = `${Math.ceil(Math.min(...[16, w - 1, h - 1]))}px  Roboto, Roboto Local`;

    const hasRowHeader = gc.grid.columns.byIndex(0)?.cellType === 'row header';
    if (gc.isColHeader && gc.gridColumn.idx > (hasRowHeader ? 1 : 0)) {
      g.fillText('' + (gc.gridColumn.idx - (hasRowHeader ? 1 : 0)), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (gc.isColHeader) {
      args.preventDefault();
      return;
    } else if ((gc.gridColumn.name == '0' || gc.gridColumn.idx == 0) && gc.gridRow >= 0) {
      const prefix = gc.gridRow > 25 ? 'A' : '';
      g.fillText(prefix + String.fromCharCode(65 + gc.gridRow % 26), x + w / 2, y + h / 2);
      args.preventDefault();
    } else if (summary && h > 0 && gc.gridRow >= 0 && gc.gridColumn.idx > 0) {
      const dataRow = this._plate._idx(gc.gridRow, gc.gridColumn.idx - 1);

      if (this.plate.data.selection.get(dataRow)) {
        g.fillStyle = 'rgba(0, 128, 255, 0.2)';
        g.fillRect(x, y, w, h);
      }

      g.beginPath();
      const r = Math.min(h / 2, w / 2) * 0.8;
      g.ellipse(x + w / 2, y + h / 2, r, r, 0, 0, 2 * Math.PI);

      if (this._colorColumn && this._colorColumn.name === 'Role') {
        const cellRole = this._colorColumn.get(dataRow);
        let color = DG.Color.lightGray;

        if (cellRole) {
          color = this._colorColumn.meta.colors.getColor(dataRow);
          // If a role filter is active, fade out non-matching wells
          if (this._activeRoleFilter && cellRole !== this._activeRoleFilter)
            color = DG.Color.lightLightGray;
        }

        g.fillStyle = DG.Color.toHtml(color);
        g.fill();
        g.strokeStyle = DG.Color.toHtml(DG.Color.gray);
        g.stroke();
      } else {
        // Default rendering if not coloring by Role
        g.strokeStyle = DG.Color.toHtml(DG.Color.gray);
        g.stroke();
      }

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
    this._onDestroy.next();
    this._onDestroy.complete();
    super.detach();
  }
}

