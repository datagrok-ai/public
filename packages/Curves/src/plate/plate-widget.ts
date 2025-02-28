import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {TYPE} from "datagrok-api/dg";
import {div} from "datagrok-api/ui";
import {safeLog, mapFromRow } from './utils';
import { IPlateWellFilter, Plate, PLATE_OUTLIER_WELL_NAME } from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {FIT_FUNCTION_4PL_REGRESSION, FIT_FUNCTION_SIGMOID, FitMarkerType, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import { FitConstants } from '../fit/const';
import { FitCellOutlierToggleArgs, setOutlier } from '../fit/fit-renderer';
import { _package } from '../package';


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
const dimensions = new Map([
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
  mapFromRowFunc: (row: DG.Row) => Record<string, any> = mapFromRow; 

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
        pw.wellDetailsDiv!.appendChild(ui.tableFromMap(pw.mapFromRowFunc(pw.plateData.rows.get(row))));
    });

    pw.root.prepend(pw.colorSelector.root);
    pw.root.append(pw.detailsDiv);
    return pw;
  }

  static analysisView(plate: Plate, options?: Partial<AnalysisOptions>): PlateWidget {
    const pw = PlateWidget.detailedView(plate.data);
    const gridRoot = pw.grid.root;
    const detailsRoot = pw.detailsDiv;
    pw.plateActionsDiv = ui.div();
    pw.plateDetailsDiv = ui.div();
    detailsRoot!.prepend(pw.plateDetailsDiv);
    detailsRoot!.append(pw.plateActionsDiv);
    const defaultOptions: AnalysisOptions = {
      roleName: 'layout', concentrationName: 'concentration', valueName: 'readout',
      normalize: true, controlColumns: ['High Control', 'Low Control'], autoFilterOutliers: true,
      categorizeFormula: '${rSquared} > 0.8 && ${Hill} > 0.25 && ${Max} > 80 && ${Max} < 120', statisticsColumns: ['Min', 'Max', 'rSquared', 'Hill', 'IC50', 'AUC'],
    };

    // system of aliases for the column names
    const aliases = {
      concentrationName: ['conc', 'dose', 'conc.', 'dilution', 'concentrations'],
      valueName: ['response', 'value', 'signal', 'raw data', 'raw signal', 'raw', 'response', 'values'],
      roleName: ['role', 'plate layout', 'plate positions', 'group', 'compound', 'well']
    };
    Object.entries(aliases).forEach(([key , value]) => {
      const colName = defaultOptions[key as keyof typeof defaultOptions] as string;
      if (plate.data.columns.contains(colName))
        return;
      const alias = value.find((v) => plate.data.columns.contains(v));
      if (alias)
        (defaultOptions[key as keyof typeof defaultOptions] as string) = alias;
    })

    // find control columns
    if (plate.data.columns.contains(defaultOptions.roleName) && !defaultOptions.controlColumns.every((cc) => plate.data.col(defaultOptions.roleName)!.categories.includes(cc))) {
      //  if control columns are not provided, try to find them in the data
      const controlCols = plate.data.col(defaultOptions.roleName)!.categories.filter((cat) => cat.toLowerCase().includes('control'));
      // take maximum of 2 control columns
      defaultOptions.controlColumns = controlCols.slice(0, 2);
    }

    const actOptions: AnalysisOptions & {normalizedColName?: string} = {...defaultOptions, ...options};

    const statisticsAliases = {
      'rSquared': ['rsquared', 'r2', 'r squared', 'r^2'],
      'slope': ['slope', 'hill', 'steepness', 'hill slope'],
      'bottom': ['bottom', 'min', 'minimum', 'miny', 'min y'],
      'top': ['top', 'max', 'maximum', 'maxy', 'max y'], 
      'interceptX': ['interceptx', 'intercept x','ic50', 'ic 50', 'ic-50', 'ic_50', 'ic 50 value', 'ic50 value', 'ic50value', 'ec50', 'ec 50', 'ec-50', 'ec_50', 'ec 50 value', 'ec50 value', 'ec50value'],
      'auc': ['auc', 'area under the curve', 'area under curve', 'area'],
    };

    const actualStatNames: Record<string, string> = {}; // will hold actual statistic name to its alias column name
    actOptions.statisticsColumns.forEach((stat) => {
      const alias = Object.entries(statisticsAliases).find(([_, aliases]) => aliases.includes(stat.toLowerCase()));
      if (alias)
        actualStatNames[alias[0]] = stat;
    });


    if (!gridRoot || !detailsRoot)
      return pw;
    // place them horizontally
    const container = ui.divH([gridRoot, detailsRoot]);
    pw.root.appendChild(container);
    detailsRoot.style.flex = '0 0 min(250px, 25%)';
    gridRoot.style.flexGrow = '1';
    gridRoot.style.removeProperty('width');

    const drFilterOptions: IPlateWellFilter = {exclude: {[actOptions.roleName]: actOptions.controlColumns}};
    let normed = false;
    if (actOptions.normalize && actOptions.controlColumns.length === 2) {

      const [lStats, hStats] = actOptions.controlColumns.map((colName) => {
        return plate.getStatistics(actOptions.valueName, ['mean', 'std'], {match: {[actOptions.roleName]: colName}});
      }).sort((a,b) => a.mean - b.mean);
  
      defaultOptions.plateStatistics = {
        'Z Prime': (_plate) => 1 - (3 * (hStats.std + lStats.std) / Math.abs(hStats.mean - lStats.mean)),
        'Signal to background': (_plate) => hStats.mean / lStats.mean,
      }
      if (!actOptions.plateStatistics)
        actOptions.plateStatistics = defaultOptions.plateStatistics;
      actOptions.normalizedColName = plate.normalize(actOptions.valueName, (v) => (hStats.mean - v) / (hStats.mean - lStats.mean) * 100).name;
      normed = true;
      const plateStatMap: Record<string, string> = {};
      Object.entries(actOptions.plateStatistics).forEach(([statName, statFunc]) => {
        const value = statFunc(plate);
        plateStatMap[statName] = typeof value === 'number' ? DG.format(value, '#0.0000') : value;
      });
      const statTable = ui.tableFromMap(plateStatMap);
      pw.plateDetailsDiv!.appendChild(statTable);
      // mark outliers that go outside the bounds

      if (actOptions.autoFilterOutliers) {
        plate.markOutliersWhere(actOptions.normalizedColName!, (v) => v > 106 || v < -6, drFilterOptions);
      }
    }

    const series = plate.doseResponseSeries({...drFilterOptions, value: normed ? actOptions.normalizedColName! : actOptions.valueName , concentration: actOptions.concentrationName, groupBy: actOptions.roleName});
    const seriesVals = Object.entries(series);
    const minMax = normed ? {minY: -10, maxY: 110} : 0;
    const roleCol = DG.Column.string(actOptions.roleName, seriesVals.length);
    const curveCol = DG.Column.string('Curve', seriesVals.length);
    curveCol.init((i) => JSON.stringify(
      {
                  "chartOptions": {
                    "xAxisName": actOptions.concentrationName,
                    "yAxisName": normed ? `norm(${actOptions.valueName})` : actOptions.valueName,
                    "logX": true,
                    "title": `${seriesVals[i][0]}`,
                    'clickToToggle': true,
                    ...minMax
                  },
                  // TODO: change to 4PL regression once fixed for normed data
                series: [{...seriesVals[i][1], fit: undefined, fitFunction: FIT_FUNCTION_4PL_REGRESSION, clickToToggle: true, droplines: ['IC50'], name: seriesVals[i][0]}]
      }

      ));
    roleCol.init((i) => seriesVals[i][0]);

    const df = DG.DataFrame.fromColumns([roleCol, curveCol]);
    df.name = pw.plateData.name;
    df.id = `${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}`;

    curveCol.semType ='fit';
    const curvesGrid = df.plot.grid();

    // add statistics columns
    Object.entries(actualStatNames).forEach(([statName, alias]) => {
      const params = {df: df, colName: curveCol.name, propName: statName, seriesName: 'series 0', seriesNumber: 0, newColName: alias};
      DG.Func.find({name: 'addStatisticsColumn'})[0].prepare(params).callSync({processed: false});
    });
    if (actualStatNames['interceptX'])
      df.col(actualStatNames['interceptX']) && (df.col(actualStatNames['interceptX'])!.meta.format = 'scientific');

    //categorize the curves
    
    const criteriaCol = df.columns.addNewString(df.columns.getUnusedName('Criteria'));
    criteriaCol.applyFormula(`if(${actOptions.categorizeFormula}, "Qualified", "Fails Criteria")`, 'string');
    criteriaCol.meta.colors.setCategorical({ "Fails Criteria":4294922560, "Qualified":4283477800 });

    curvesGrid.root.style.width = '100%';
    pw.root.style.display = 'flex';
    pw.root.style.flexDirection = 'column';
    pw.root.style.height = '100%';

    pw.root.appendChild(curvesGrid.root);

    // when selecting a cell on plate, go to appropriate curve and mark the corresponding point with different marker
    // remember the previous change, so that we can revert it
    let prevSelection: {seriesIndex: number, pointIndex: number, markerType: FitMarkerType, markerSize: number, markerColor: string, curvesGridCell?: DG.GridCell} | null = null;

    pw.mapFromRowFunc = (row) => mapFromRow(row, (rowIdx, checkBoxState) => {
      plate._markOutlier(rowIdx, checkBoxState);
      if (prevSelection && prevSelection.curvesGridCell)
      setOutlier(prevSelection.curvesGridCell, {x: 0, y: 0, outlier: !checkBoxState}, 0, prevSelection.pointIndex);
      pw.grid.invalidate();
    });

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
      const value: number = pw.plateData.get(normed ? actOptions.normalizedColName! : actOptions.valueName, row);

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
        points[pointInSeriesIndex].color = DG.Color.toHtml(DG.Color.green);
        curveCol.set(seriesIndex, JSON.stringify(curveJSON));
      }

      // scroll in grid
      const curveGridRow = curvesGrid.tableRowToGrid(seriesIndex);
      prevSelection.curvesGridCell = curvesGrid.cell(curveCol.name, curveGridRow);
      prevSelection.curvesGridCell && curvesGrid.scrollToCell(curveCol.name, curveGridRow);
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


    const s = DG.debounce(curvesGrid.onAfterDrawContent, 300).subscribe(() => {
      s.unsubscribe();
      curvesGrid.col(curveCol.name) && (curvesGrid.col(curveCol.name)!.width = 400);
      curvesGrid.props.rowHeight = 200;
    });

    if (actOptions.submitAction) {
      const btn = ui.button('Save to ELN', () => actOptions.submitAction!(plate, df));
      pw.plateActionsDiv!.appendChild(btn);
    }
    return pw;
  }

  get plateData() { return this._plateData; }
  set plateData(t: DG.DataFrame) {
    this._plateData = t;
    this._colorColumn = t.columns.firstWhere((col) => col.name == 'activity' && col.type == TYPE.FLOAT)

      ?? t.columns.firstWhere((col) => (col.name == 'concentration' || col.name == 'concentrations') && col.type == TYPE.FLOAT)
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


export async function getPlatesFolderPreview(files: DG.FileInfo[]): Promise<DG.Widget | DG.ViewBase | undefined> {

    const csvFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.csv'));
    let csvView: DG.Widget | undefined = undefined;
    if (csvFiles.length > 2) {

      const plate = Plate.fromPlates(await Promise.all(csvFiles.map(async (f) => await Plate.fromCsvTableFile(f.fullPath, f.name.toLowerCase().substring(0, f.name.length - 4)))));
      csvView = PlateWidget.analysisView(plate, {submitAction: () => {grok.shell.info('Plate Submitted')}});
      plate.data.name = `${csvFiles.map((file) => file.name.slice(0, file.name.length - 4)).join('_')}.csv`;
      if (csvFiles.length === files.length)
        return csvView;
    }
    const multiView = new DG.MultiView({viewFactories: {}});

    if (csvView) {
      multiView.addView('Plate 1', () => DG.View.fromRoot(csvView.root), true);
    }

    const xlsxFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.xlsx') && f?.name?.toLowerCase().includes('plate'));

    if (xlsxFiles.length == 0)
      return csvView;

    for (const xlsxFile of xlsxFiles) {
      try {
        const plate = await Plate.fromExcelFileInfo(xlsxFile);
        const pw = PlateWidget.analysisView(plate, {submitAction: () => {grok.shell.info('Plate Submitted')}});
        const v = DG.View.fromRoot(pw.root);
        v.name = xlsxFile.name.substring(0, xlsxFile.name.length - 5);
        multiView.addView(v.name, () => v, true);
      } catch (e) {
        _package.logger.error(e);
      }
    }
    setTimeout(() => {
      const header: HTMLElement | null = multiView.root.querySelector('.d4-tab-header-stripe');
      if (!header)
        return;
      header.style.overflow = 'scroll';

      const headerItems = header.querySelectorAll('.d4-tab-header');
      headerItems?.forEach((el) => (el as HTMLElement).style.whiteSpace = 'nowrap');
    }, 300)

    return multiView;
}
