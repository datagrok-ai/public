import * as DG from 'datagrok-api/dg';
import {MARKER_TYPE, TYPE} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {div} from 'datagrok-api/ui';
import {mapFromRow, safeLog, toExcelPosition} from './utils';
import {Plate, PLATE_OUTLIER_WELL_NAME} from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {IPlateWellValidator, plateWellValidators} from './plate-well-validators';
import {debounceTime, filter} from 'rxjs/operators';


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

/** Visualizes multiple layers of {@link Plate} in tabbed interface and lets you edit it */
export class PlateWidget extends DG.Widget {
  _colorColumn?: DG.Column | undefined;
  _plate: Plate = new Plate(8, 12);
  grid: DG.Grid = DG.Viewer.heatMap(DG.DataFrame.create());
  colorColumnName: string = '';
  colorSelector?: DG.InputBase;
  detailsDiv?: HTMLElement;
  wellDetailsDiv?: HTMLElement;
  plateDetailsDiv?: HTMLElement;
  plateActionsDiv?: HTMLElement;
  tabs: DG.TabControl = DG.TabControl.create();
  _editable: boolean = false;
  mapFromRowFunc: (row: DG.Row) => Record<string, any> = mapFromRow;
  grids: Map<string, DG.Grid> = new Map();
  wellValidators: IPlateWellValidator[] = plateWellValidators;
  wellValidationErrors: Map<string, string[]> = new Map();

  get editable() { return this._editable; }
  set editable(x: boolean) {
    this._editable = x;
    this.syncGrids();
  }

  constructor() {
    super(ui.div([], 'curves-plate-widget'));

    this.tabs.root.style.width = '100%';
    this.tabs.root.style.height = '100%';
    this.tabs.addPane('Summary', () => this.grid.root);
    this.root.appendChild(this.tabs.root);

    this.grid.props.showHeatmapScrollbars = false;
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
  }

  static fromPlate(plate: Plate) {
    return PlateWidget.detailedView(plate);
  }

  static detailedView(plate: Plate): PlateWidget {
    const pw = new PlateWidget();
    pw.grid.root.style.width = '100%';
    pw.plate = plate;
    pw.colorSelector = ui.input.column('Color by', {
      table: plate.data,
      value: pw._colorColumn,
      onValueChanged: (v) => { pw._colorColumn = v; pw.grid.invalidate(); }
    });
    pw.detailsDiv = ui.divV([]);
    pw.wellDetailsDiv = div();
    pw.detailsDiv.appendChild(pw.wellDetailsDiv!);
    pw.grid.onCurrentCellChanged.pipe(filter(gc => gc.gridRow >= 0 && gc.gridColumn.idx > 0)).subscribe((gc) => {
      ui.empty(pw.wellDetailsDiv!);
      const map = pw.mapFromRowFunc(pw.plate.data.rows.get(plate._idx(gc.gridRow, gc.gridColumn.idx - 1)))
      pw.wellDetailsDiv!.appendChild(ui.tableFromMap(map));
    });

    pw.root.prepend(pw.colorSelector.root);
    pw.root.append(pw.detailsDiv);
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

    grid.onCellRender.subscribe((args) => this.renderCell(args, isSummary));
  }

  get plate(): Plate { return this._plate; }
  set plate(p: Plate) {
    this._plate = p;
    this.refresh();
  }

  /** Re-renders the UI */
  refresh() {
    this.tabs.clear();
    this.tabs.addPane('Summary', () => this.grid.root);
    for (const layer of this.plate.getLayerNames()) {
      this.tabs.addPane(layer, () => {
        const df = this.plate.toGridDataFrame(layer);
        const grid = DG.Viewer.heatMap(df);
        grid.columns.add({gridColumnName: '0', cellType: 'string', index: 1});   // to show row numbers in Excel format - will be done in initPlateGrid
        df.onValuesChanged.pipe(debounceTime(1000)).subscribe(() => {
          const p = this.plate;
          for (let i = 0; i < df.rowCount; i++)
            for (let j = 0; j < df.columns.length - 1; j++)
              p.data.col(layer)!.set(p._idx(i, j), df.get(`${j + 1}`, i));
          this.wellValidationErrors = p.validateWells(this.wellValidators);
        });
        this.initPlateGrid(grid, false);
        this.grids.set(layer, grid);
        return grid.root;
      });
    }
    this.tabs.currentPane = this.tabs.getPane('Summary');

    this.syncGrids();

    const t = this.plate.data;
    this._colorColumn = t.columns.firstWhere((col) => col.name == 'activity' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => (col.name == 'concentration' || col.name == 'concentrations') && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.name != 'row' && col.name != 'col' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.type == TYPE.FLOAT)
      ?? (t.columns.length > 0 ? t.columns.byIndex(0) : undefined);

    // row header + all columns
    this.grid.dataFrame = DG.DataFrame.create(this.plate.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= this.plate.cols; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    this.grid.invalidate();
  }

  renderCell(args:DG.GridCellRenderArgs, summary: boolean = false) {
    const gc = args.cell;
    args.g.fillStyle = 'grey'; //(args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
    args.g.strokeStyle = 'grey';
    args.g.lineWidth = 1;
    let g = args.g;
    let x = args.bounds.x, y = args.bounds.y, w = args.bounds.width, h = args.bounds.height;
    g.textAlign = 'center';
    g.textBaseline = 'middle';
    g.font = `${Math.ceil(Math.min(...[16, w - 1, h - 1]))}px  Roboto, Roboto Local`;
    const isColoredByConc = this._colorColumn?.name?.toLowerCase()?.includes('conc');

    // column header
    if (gc.isColHeader && gc.gridColumn.idx > 0)
      g.fillText('' + gc.gridColumn.idx, x + w / 2, y + h / 2);
    // row header
    else if ((gc.gridColumn.name == '0' || gc.gridColumn.idx == 0) && gc.gridRow >= 0) {
      const prefix = gc.gridRow > 25 ? 'A' : '';
      g.fillText(prefix + String.fromCharCode(65 + gc.gridRow % 26), x + w / 2, y + h / 2);
    }
    else if (summary && h > 0 && gc.gridRow >= 0 && gc.gridColumn.idx > 0) {
      const dataRow  = this._plate._idx(gc.gridRow, gc.gridColumn.idx - 1);

      g.beginPath();
      const r = Math.min(h / 2, w / 2) * 0.8;
      g.ellipse(x + w / 2, y + h / 2, r, r, 0, 0, 2 * Math.PI);
      if (this._colorColumn) {
        if (this._colorColumn.isCategorical && this._colorColumn.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
          this._colorColumn.meta.colors.setCategorical();
        const color = this._colorColumn.isNone(dataRow) ? DG.Color.white : this._colorColumn.isNumerical
          ? this.getColor(dataRow, isColoredByConc)
          : this._colorColumn.meta.colors.getColor(dataRow);
        g.fillStyle = DG.Color.toHtml(color);
        g.fill();
        const outlierCol = this.plate.data.col(PLATE_OUTLIER_WELL_NAME);
        if (outlierCol?.get(dataRow))
          DG.Paint.marker(g, MARKER_TYPE.CROSS_X_BORDER, x + w / 2, y + w / 2, DG.Color.red, r);
      }

      g.stroke();
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

  /** Applies plate widget options to all grids */
  syncGrids() {
    for (const grid of this.grids.values())
      grid.props.allowEdit = this.editable;
  }
}
