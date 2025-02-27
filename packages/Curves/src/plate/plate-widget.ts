import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {TYPE} from "datagrok-api/dg";
import {div} from "datagrok-api/ui";
import {safeLog, tableFromRow } from './utils';
import { IPlateWellFilter, Plate, PLATE_OUTLIER_WELL_NAME } from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {FIT_FUNCTION_4PL_REGRESSION, FitMarkerType, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import { FitConstants } from '../fit/const';
import { FitCellOutlierToggleArgs } from '../fit/fit-renderer';


type AnalysisOptions = {
  roleName?: string,
  concentrationName?: string,
  valueName?: string
  controlColumns?: string[],
  normalize?: boolean,
  autoFilterOutliers?: boolean
}


const colorScheme = [DG.Color.white, DG.Color.gray];
const dimensions = new Map([
  [96, {rows: 8, cols: 12}],
  [384, {rows: 16, cols: 24}],
  [1536, {rows: 32, cols: 48}],
]);

export class PlateWidget extends DG.Widget {

  _plateData: DG.DataFrame = DG.DataFrame.create();
  _colorColumn?: DG.Column<number>;
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

  constructor() {
    super(ui.div([], 'curves-plate-widget'));
    this.root.appendChild(this.grid.root);

    this.subs.push(this.grid.onAfterDrawContent.subscribe(() => {
      this.grid.root.querySelectorAll('.d4-range-selector').forEach((el) => (el as HTMLElement).style.display = 'none');
    }));

    this.grid.onCellRender.subscribe((args) => this.renderCell(args));
  }

  static detailedView(table: DG.DataFrame): PlateWidget {
    const pw = new PlateWidget();
    pw.grid.root.style.width = '100%';
    pw.plateData = table;
    pw.colorSelector = ui.input.column('Color by', {
      table: table,
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
        pw.wellDetailsDiv!.appendChild(ui.tableFromMap(tableFromRow(pw.plateData.rows.get(row))));
    });

    pw.root.prepend(pw.colorSelector.root);
    pw.root.append(pw.detailsDiv);
    return pw;
  }

  static analysisView(plate: Plate, options?: AnalysisOptions): PlateWidget {
    const pw = PlateWidget.detailedView(plate.data);
    const gridRoot = pw.grid.root;
    const detailsRoot = pw.detailsDiv;
    pw.plateActionsDiv = ui.div();
    pw.plateDetailsDiv = ui.div();
    detailsRoot!.prepend(pw.plateDetailsDiv);
    detailsRoot!.append(pw.plateActionsDiv);



    if (!gridRoot || !detailsRoot)
      return pw;
    // place them horizontally
    const container = ui.divH([gridRoot, detailsRoot]);
    pw.root.appendChild(container);
    detailsRoot.style.flex = '0 0 min(250px, 25%)';
    gridRoot.style.flexGrow = '1';
    gridRoot.style.removeProperty('width');

    const defaultOptions: Required<AnalysisOptions> = {roleName: 'layout', concentrationName: 'concentration', valueName: 'readout', normalize: true, controlColumns: ['High Control', 'Low Control'], autoFilterOutliers: true};
    const actOptions: Required<AnalysisOptions> & {normalizedColName?: string} = {...defaultOptions, ...options};
    const drFilterOptions: IPlateWellFilter = {exclude: {[actOptions.roleName]: actOptions.controlColumns}};
    let normed = false;
    if (actOptions.normalize && actOptions.controlColumns.length === 2) {
      const [lMean, hMean] = actOptions.controlColumns.map((colName) => jStat.mean(plate.fieldValues(actOptions.valueName, {match: {[actOptions.roleName]: colName}}))).sort((a,b) => a - b);
      actOptions.normalizedColName = plate.normalize(actOptions.valueName, (v) => (hMean - v) / (hMean - lMean) * 100).name;
      normed = true;
      // mark outliers that go outside the bounds

      if (actOptions.autoFilterOutliers) {
        const values = plate.values([actOptions.normalizedColName], drFilterOptions);
        values.forEach((v) => {
          if (v[actOptions.normalizedColName!] > 106 || v[actOptions.normalizedColName!] < -6) {
            plate._markOutlier(v.innerDfRow, true);
          }
        })
      }
    }

    const series = plate.doseResponseSeries({...drFilterOptions, value: actOptions.valueName, concentration: actOptions.concentrationName, groupBy: actOptions.roleName});

    const seriesVals = Object.entries(series);

    const minXYOpts = normed ? {minY: 0, maxY: 100} : {};

    const roleCol = DG.Column.string(actOptions.roleName, seriesVals.length);
    const curveCol = DG.Column.string('Curve', seriesVals.length);
    curveCol.init((i) => JSON.stringify(
      {
                  "chartOptions": {
                    "xAxisName": actOptions.concentrationName,
                    "yAxisName": actOptions.valueName,
                    "logX": true,
                    "title": `${seriesVals[i][0]}`,
                    'clickToToggle': true,
                  },
                series: [{...seriesVals[i][1], fit: undefined, fitFunction: FIT_FUNCTION_4PL_REGRESSION, clickToToggle: true}]
      }

      ));
    roleCol.init((i) => seriesVals[i][0]);

    const df = DG.DataFrame.fromColumns([roleCol, curveCol]);
    df.name = pw.plateData.name;
    curveCol.semType ='fit';
    const curvesGrid = df.plot.grid();

    const rSquaredFuncParams = {df: df, colName: curveCol.name, propName: 'rSquared', seriesName: 'series 0', seriesNumber: 0, newColName: 'rSquared'};
    const slopeFuncParams = {df: df, colName: curveCol.name, propName: 'slope', seriesName: 'series 0', seriesNumber: 0, newColName: 'slope'};
    const curveMinFuncParams = {df: df, colName: curveCol.name, propName: 'bottom', seriesName: 'series 0', seriesNumber: 0, newColName: 'min'};
    const curveMaxFuncParams = {df: df, colName: curveCol.name, propName: 'top', seriesName: 'series 0', seriesNumber: 0, newColName: 'max'};
    DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(rSquaredFuncParams).callSync({processed: false});
    DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(slopeFuncParams).callSync({processed: false});
    DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(curveMinFuncParams).callSync({processed: false});
    DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(curveMaxFuncParams).callSync({processed: false});

    curvesGrid.root.style.width = '100%';
    pw.root.style.display = 'flex';
    pw.root.style.flexDirection = 'column';
    pw.root.style.height = '100%';

    pw.root.appendChild(curvesGrid.root);

    // when selecting a cell on plate, go to appropriate curve and mark the corresponding point with different marker
    // remember the previous change, so that we can revert it
    let prevSelection: {seriesIndex: number, pointIndex: number, markerType: FitMarkerType, markerSize: number, markerColor: string} | null = null;

    function clearPreviousSelection() {
      try {
        if (!prevSelection || (prevSelection.seriesIndex ?? -1) < 0) {
          return;
        }
        const series = curveCol.get(prevSelection.seriesIndex);
        if (!series)
          return;
        const parsed = JSON.parse(series);
        const points: IFitPoint[] = parsed.series[0]?.points;
        if (points && points.length && points.length > prevSelection.pointIndex) {
          points[prevSelection.pointIndex].marker = prevSelection.markerType;
          points[prevSelection.pointIndex].size = prevSelection.markerSize;
          points[prevSelection.pointIndex].color = prevSelection.markerColor;
          curveCol.set(prevSelection.seriesIndex, JSON.stringify(parsed));
        }
      } catch (e) {
        console.error(e);
      } finally {
        prevSelection = null;
      }

    }
    // highlight selected point in the curve
    pw.subs.push(pw.grid.onCurrentCellChanged.subscribe((gc) => {
      clearPreviousSelection();
      if (gc?.gridRow == null || gc?.gridRow == -1 || !gc?.gridColumn || gc?.gridColumn.idx == 0)
        return;
      const row = pw.dataRow(gc);
      if (row == undefined || row < 0)
        return;
      const catValue = pw.plateData.get(actOptions.roleName, row)?.toLowerCase();
      if (!catValue)
        return;
      const seriesIndex = seriesVals.findIndex(([serName, _]) => serName?.toLowerCase() === catValue);
      if (seriesIndex < 0)
        return;
      const conscentration: number = pw.plateData.get(actOptions.concentrationName, row);
      const value: number = pw.plateData.get(actOptions.valueName, row);

      const pointInSeriesIndex: number = seriesVals[seriesIndex][1].points.findIndex((p) => p.x === conscentration && p.y === value);
      if (pointInSeriesIndex < 0)
        return;

      prevSelection = {seriesIndex, pointIndex: pointInSeriesIndex, markerType: seriesVals[seriesIndex][1].points[pointInSeriesIndex].marker ?? DG.MARKER_TYPE.CIRCLE,
         markerSize: seriesVals[seriesIndex][1].points[pointInSeriesIndex].size ?? FitConstants.POINT_PX_SIZE, markerColor: seriesVals[seriesIndex][1].points[pointInSeriesIndex].color ?? DG.Color.toHtml(DG.Color.getCategoricalColor(0))};

      const curveJSON = JSON.parse(curveCol.get(seriesIndex)!);
      const points: IFitPoint[] = curveJSON.series[0]?.points;
      if (points && points.length && points.length > pointInSeriesIndex) {
        points[pointInSeriesIndex].marker = DG.MARKER_TYPE.SQUARE;
        points[pointInSeriesIndex].size = FitConstants.POINT_PX_SIZE * 2;
        points[pointInSeriesIndex].color = DG.Color.toHtml(DG.Color.getCategoricalColor(1));
        curveCol.set(seriesIndex, JSON.stringify(curveJSON));
      }

      // scroll in grid
      const curveGridRow = curvesGrid.tableRowToGrid(seriesIndex);
      curvesGrid.scrollToCell(curveCol.name, curveGridRow);
    }));

    // mark outliers in the original plate if it is switched from curve manually

    pw.subs.push(grok.events.onCustomEvent('fit-cell-outlier-toggle').subscribe((args: FitCellOutlierToggleArgs) => {
      if (!args || !args.gridCell || !args.series || args.pointIdx == null || args.gridCell.cell.column !== curveCol)
        return;
      const point: IFitPoint = args.series.points[args.pointIdx];
      if (point.meta !== null) {
        plate._markOutlier(point.meta, !!point.outlier);
        pw.grid.invalidate();
      }


    }));


    setTimeout(() => {
      curvesGrid.col(curveCol.name)!.width = 400;
      curvesGrid.props.rowHeight = 200;
    }, 300);

    return pw;
  }

  get plateData() { return this._plateData; }
  set plateData(t: DG.DataFrame) {
    this._plateData = t;
    this._colorColumn = t.columns.firstWhere((col) => col.name == 'activity' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.name == 'activity' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.name != 'row' && col.name != 'col' && col.type == TYPE.FLOAT)
      ?? t.columns.firstWhere((col) => col.type == TYPE.FLOAT)
      ?? t.columns.byIndex(0);

    let rowCol: DG.Column<number> = this._plateData.col('row')!;
    let colCol: DG.Column<number> = this._plateData.col('col')!;
    this.rows = rowCol?.stats?.max ?? dimensions.get(t.rowCount)?.rows;
    this.cols = colCol?.stats?.max ?? dimensions.get(t.rowCount)?.cols;

    if (this.rows == null || this.cols == null)
      throw 'Row/col columns not found, and dataframe length is not of the recognized sizes (92, 384, 1526)';

    for (let i = 0; i < this._plateData.rowCount; i++)
      if (rowCol && colCol)
        this._posToRow.set(`${rowCol.get(i)}:${colCol.get(i)}`, i);
      else
        this._posToRow.set(`${Math.floor(i / this.cols)}:${i % this.cols + 1}`, i);

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
    return this._posToRow.get(`${gc.gridRow}:${gc.gridColumn.idx}`);
  }

  renderCell(args:DG.GridCellRenderArgs) {
    const gc = args.cell;
    args.g.fillStyle = 'grey'; //(args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
    args.g.strokeStyle = 'grey';
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
    else if (gc.gridColumn.idx == 0 && gc.gridRow >= 0)
      g.fillText(String.fromCharCode(65 + gc.gridRow), x + w / 2, y + h / 2);
    else if (h > 0 && dataRow != null) {
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
    args.preventDefault();
  }

  getColor(dataRow: number, isLog?: boolean) {
    const val = this._colorColumn!.get(dataRow)!;
    const reducedVal = isLog ? safeLog((val - this._colorColumn!.min) * 1e9) : val;
    const min = (isLog ? safeLog(Math.max(this._colorColumn!.min, 1)) : this._colorColumn!.min);
    const max = isLog ? safeLog((this._colorColumn!.max - this._colorColumn!.min) * 1e9) : this._colorColumn!.max;
    return DG.Color.scaleColor(reducedVal, min, max, undefined, colorScheme);
  }
}
