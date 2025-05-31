import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {TYPE} from "datagrok-api/dg";
import {div} from "datagrok-api/ui";
import {safeLog, mapFromRow, toExcelPosition } from './utils';
import {IPlateWellFilter, Plate, PLATE_OUTLIER_WELL_NAME, randomizeTableId} from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {FIT_FUNCTION_4PL_REGRESSION, FIT_FUNCTION_SIGMOID, FitMarkerType, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import { FitConstants } from '../fit/const';
import {FitCellOutlierToggleArgs, setOutlier} from '../fit/fit-renderer';
import { _package } from '../package';
import {savePlate} from "../plates/plates-crud";
import { IPlateWellValidator, plateWellValidators } from './plate-well-validators';
import { debounceTime } from 'rxjs/operators';


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
  _plateData: DG.DataFrame = DG.DataFrame.create();
  _colorColumn?: DG.Column;
  grid: DG.Grid = DG.Viewer.heatMap(this._plateData);
  _posToRow: Map<String, number> = new Map();
  rows = 0;
  cols = 0;
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
    pw.grid.onCurrentCellChanged.subscribe((gc) => {
      ui.empty(pw.wellDetailsDiv!);
      const row = pw.dataRow(gc);
      if (row != null && row >= 0)
        pw.wellDetailsDiv!.appendChild(ui.tableFromMap(pw.mapFromRowFunc(pw.plateData.rows.get(row))));
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

  get plate(): Plate {
    const plate = new Plate(this.rows, this.cols);
    for (const column of this.plateData.columns)
      if (column.name != 'row' && column.name != 'col') {
        plate.data.columns.addNew(column.name, column.type).init(i => {
          const [row, col] = plate.rowIndexToExcel(i);
          const dataIdx = this._posToRow.get(`${row}:${col}`);
          return dataIdx != null ? column.get(dataIdx) : null;
        });
      }
    return plate;
  }
  set plate(p: Plate) {
    this.plateData = p.data;
  }

  get plateData() { return this._plateData; }
  set plateData(t: DG.DataFrame) {
    this._plateData = t;
    this._colorColumn = t.columns.firstWhere((col) => col.name == 'activity' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => (col.name == 'concentration' || col.name == 'concentrations') && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.name != 'row' && col.name != 'col' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.type == TYPE.FLOAT)
      ?? (t.columns.length > 0 ? t.columns.byIndex(0) : undefined);

    let rowCol: DG.Column<number> = this._plateData.col('row')!;
    let colCol: DG.Column<number> = this._plateData.col('col')!;
    this.rows = rowCol?.stats?.max ?? dimensions.get(t.rowCount)?.rows!;
    this.cols = colCol?.stats?.max ?? dimensions.get(t.rowCount)?.cols!;

    if (this.rows == null || this.cols == null)
      throw 'Row/col columns not found, and dataframe length is not of the recognized sizes (96, 384, 1536)';

    for (let i = 0; i < this._plateData.rowCount; i++)
      if (rowCol && colCol)
        this._posToRow.set(`${rowCol.get(i)}:${colCol.get(i)}`, i);
      else
        this._posToRow.set(`${Math.floor(i / this.cols)}:${i % this.cols}`, i);

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
    this.syncGrids();
    this.tabs.currentPane = this.tabs.getPane('Summary');

    // row header + all columns
    this.grid.dataFrame = DG.DataFrame.create(this.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= this.cols; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    this.grid.invalidate();
  }

  dataRow(gc: DG.GridCell): number | undefined {
    // we do not increment gridColumn index by 1 because there is a row number column with index 0 which we do not count
    // so data columns start with 1 anyway
    return this._posToRow.get(`${gc.gridRow}:${gc.gridColumn.idx - 1}`);
  }

  renderCell(args:DG.GridCellRenderArgs, summary: boolean = false) {
    const gc = args.cell;
    args.g.fillStyle = 'grey'; //(args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
    args.g.strokeStyle = 'grey';
    args.g.lineWidth = 1;
    let g = args.g;
    let x = args.bounds.x, y = args.bounds.y, w = args.bounds.width, h = args.bounds.height;
    const dataRow  = this.dataRow(gc);
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
    else if (summary && h > 0 && dataRow != null) {
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
        const outlierCol = this.plateData.col(PLATE_OUTLIER_WELL_NAME);
        if (outlierCol?.get(dataRow)) {
          g.strokeStyle = DG.Color.toHtml(DG.Color.red);
          g.lineWidth = 2;
          g.moveTo(x + w / 2 - r, y + h / 2 - r);
          g.lineTo(x + w / 2 + r, y + h / 2 + r);
          g.moveTo(x + w / 2 + r, y + h / 2 - r);
          g.lineTo(x + w / 2 - r, y + h / 2 + r);
          g.stroke();
        }
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
