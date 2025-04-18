import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import wu from 'wu';
import {
  assure,
  firstWhere,
  getMaxPosition,
  jstatStatistics,
  JSTATStatistics,
  numToExcel,
  parseExcelPosition, toStandardSize
} from "./utils";
import type ExcelJS from 'exceljs';
import {findPlatePositions, getPlateFromSheet} from "./excel-plates";
import {FitFunctionType, FitSeries} from '@datagrok-libraries/statistics/src/fit/new-fit-API';
import {AnalysisOptions, PlateWidget} from './plate-widget';
import {inspectCurve} from '../fit/fit-renderer';


/** Represents a well in the experimental plate */
export interface PlateWell {
  row: number;
  col: number;
  volume?: number;
  concentration?: number;
  [key: string]: any;
}


interface IPlateCsvImportOptions {
  rowColumn?: string;
  colColumn?: string;
  positionColumn?: string;
  volumeColumn?: string;
  concentrationColumn?: string;
  field?: string;   // used as a field name when a 12x8 table is passed; if not passed, defaults to 'value'
}


export interface IPlateWellFilter {
  includeEmpty?: boolean;
  match?: {[key: string]: any};
  exclude?: {[key: string]: any};
}


interface ISeriesData {
  x: number[];
  y: number[];
}

// export class PlateLayer<TData> {
//   plate: Plate;
//   data: Float64Array ArrayLike<number> DG.Column<TData>;
//   role: 'concentration' | 'volume';
//
//   /** Array of non-empty values */
//   values(filter?: IPlateWellFilter): Array<any> {
//     const col = this.data.columns.byName(data.name);
//     assure(col != null, `Field does not exist: ${field}`);
//
//     const result = [];
//     for (let i = 0; i < this.rows * this.cols; i++)
//       if ((filter?.includeEmpty ?? false) || !col.isNone(i))
//         result.push(col.isNone(i) ? null : col.get(i));
//
//     return result;
//   }
// }

export const PLATE_OUTLIER_WELL_NAME = 'Outlier';

export function randomizeTableId() {
  return `${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}`;
}

/** Represents experimental plate (typically 96-, 384-, or 1536-well assay plates) */
export class Plate {
  data: DG.DataFrame;
  //layers: PlateLayer[];
  rows: number = 8;
  cols: number = 12;

  constructor(rows: number, cols: number) {
    this.data = DG.DataFrame.create(rows * cols);
    this.rows = rows;
    this.cols = cols;
  }

  /** Creates an empty plate of the standard size to fit the specified positions. t*/
  static autoSize(positions: Iterable<[number, number]>): Plate {
    const [rows, cols] = toStandardSize(getMaxPosition(positions));
    return new Plate(rows, cols);
  }

  _idx(row: number, col: number): number {
    return row * this.cols + col;
  }

  //well(row: number, col: number): PlateWell { return new PlateWell(); }

  _markOutlier(row: number, flag: boolean = true) {
    const outlierCol = this.data.columns.getOrCreate(PLATE_OUTLIER_WELL_NAME, DG.TYPE.BOOL);
    outlierCol.set(row, flag);
  }

  markOutlier(row: number, col: number, flag: boolean = true) {
    this._markOutlier(this._idx(row, col), flag);
  }

  _isOutlier(row: number): boolean {
    const outlierCol = this.data.columns.getOrCreate(PLATE_OUTLIER_WELL_NAME, DG.TYPE.BOOL);
    return outlierCol.get(row);
  }
  isOutlier(row: number, col: number): boolean {
    return this._isOutlier(this._idx(row, col));
  }

  markOutliersWhere(field: string, valueFunc: (fieldValue: any) => boolean, filter?: IPlateWellFilter) {
    this.values([field], filter).forEach((v) => {
      if(valueFunc(v[field]))
        this._markOutlier(v.innerDfRow, true);
    })
  }

  static fromTable(table: DG.DataFrame, field?: string): Plate {
    const rows = table.rowCount;
    // remove row letters
    function containsRowLetters(c: DG.Column): boolean {
      return c.type == DG.TYPE.STRING && DG.range(rows).every((i) => numToExcel(i) == c.get(i).toUpperCase());
    }

    const dataColumns = wu(table.columns).filter((c) => !containsRowLetters(c)).toArray();
    const cols = dataColumns.length;
    const plate = new Plate(rows, cols);

    const nonEmptyCols = dataColumns.filter((c) => c.stats.missingValueCount != rows);
    const toString = (nonEmptyCols.length > 0 && nonEmptyCols.findIndex((c) => c.type != nonEmptyCols[0].type) != -1);
    const type = toString
      ? DG.TYPE.STRING
      : firstWhere(dataColumns, c => c.stats.missingValueCount != rows)?.type ?? DG.TYPE.STRING;

    const dataColumn = plate.data.columns.addNew(table.name ?? field, type);

    for (let colIdx = 0; colIdx < dataColumns.length; colIdx ++)
      if (dataColumns[colIdx].stats.missingValueCount != rows)
        for (let i = 0; i < rows; i++)
          dataColumn.set(plate._idx(i, colIdx), type != DG.TYPE.STRING ? dataColumns[colIdx].get(i) : dataColumns[colIdx].getString(i), false);

    return plate;
  }


  /** Constructs a plate from a dataframe where each row corresponds to a well.
   * Automatically detects position column (has to be in Excel notation). */
  static fromTableByRow(table: DG.DataFrame, field?: string): Plate {
    const posCol = wu(table.columns)
      .find((c) => c.type == DG.TYPE.STRING && /^([A-Za-z]+)(\d+)$/.test(c.get(0)));
    if (!posCol)
      throw 'Column with well positions not found';

    const [rows, cols] = toStandardSize(getMaxPosition(wu(posCol.values()).map(parseExcelPosition)));
    const plate = new Plate(rows, cols);

    for (const col of table.columns) {
      if (col.isEmpty || col == posCol)
        continue;

      const plateCol = plate.data.columns.addNew(col.name, col.type);
      for (let r = 0; r < col.length; r++) {
        const [wellRow, wellCol] = parseExcelPosition(posCol.get(r));
        plateCol.set(wellRow * cols + wellCol, col.get(r), false);
      }
    }

    return plate;
  }

  static async fromCsvTableFile(csvPath: string, field: string, options?: DG.CsvImportOptions): Promise<Plate> {
    return this.fromTable(await grok.dapi.files.readCsv(csvPath, options), field);
  }

  // getClosestRoleName(fileName: string): string {
  //   const roleRoleNames = ['role', 'compound', 'compound name'];
  //   const concentrationRoleNames = ['concentration', 'conc', 'dilution'];
  //   const valueRoleNames = ['value', 'readout', 'response', 'activity', 'signal', 'absorbance', 'fluorescence', 'luminescence', 'intensity', 'od', 'optical density'];
  //   const foundRoleIndex = [roleRoleNames, concentrationRoleNames, valueRoleNames];

  // }

  static fromCsvTable(csv: string, field: string, options?: DG.CsvImportOptions): Plate {
    return this.fromTable(DG.DataFrame.fromCsv(csv, options), field);
  }

  static fromCsv(csv: string, options?: IPlateCsvImportOptions): Plate {
    const df = DG.DataFrame.fromCsv(csv);

    // we do not need first row with row letters
    if (DG.range(df.rowCount).every((i) => df.col(0)!.get(i) == String.fromCharCode(65 + i)))
      df.columns.remove(df.columns.byIndex(0).name);

    // sometimes there is an unnecessary last column in the end
    if (df.columns.length == 13 || df.columns.length == 25 || df.columns.length == 49 && df.col(df.columns.length - 1)?.isEmpty)
      df.columns.remove(df.columns.length - 1);

    // TODO: tall format?
    // if (df.columns.contains('pos'))
    //   return new Plate(0, 0);

    return this.fromTable(df, options?.field ?? 'value');
  }

  /** Merges the attributes from {@link plates} into one plate. */
  static fromPlates(plates: Plate[], name?: string) {
    assure(plates.length > 0, 'Array is empty.');
    assure(plates.every(p => p.rows == plates[0].rows && p.cols == plates[0].cols), 'Plate dimensions differ.');

    const result = plates.reduce((p1, p2) => p1 ? p1.merge(p2) : p2.clone());
    if (name != null)
      result.data.name = name;
    return result;
  }

  static async fromExcelFileInfo(fi: DG.FileInfo): Promise<Plate> {
    return fi.data && fi.data.length ? await Plate.fromExcel(fi.data, fi.friendlyName) : await Plate.fromExcelPath(fi.fullPath, fi.friendlyName);
  }

  static async fromExcelPath(excelPath: string, name?: string): Promise<Plate> {
    const content = await grok.dapi.files.readAsBytes(excelPath);
    return this.fromExcel(content, name);
  }

  /** Constructs the plates from the specified Excel file.
   * Plate positions are detected automatically. */
  static async fromExcel(excelBytes: Uint8Array, name?: string): Promise<Plate> {
    await DG.Utils.loadJsCss(['/js/common/exceljs.min.js']);

    //@ts-ignore
    const loadedExcelJS = window.ExcelJS as ExcelJS;
    const workbook = new loadedExcelJS.Workbook() as ExcelJS.Workbook;
    // @ts-ignore reason: some typescript error with araryLike and arrayBufferLike
    const wb= await workbook.xlsx.load(excelBytes);

    const platePositions = findPlatePositions(wb);
    if (platePositions.length == 0) throw 'Plates not found in "${excelPath}"';
    const p0 = platePositions[0];
    if (!platePositions.every((pc) => pc.rows == p0.rows || pc.cols == p0.cols))
      throw `Plate sizes differ in "${name}"`;

    return Plate.fromPlates(platePositions.map(p => getPlateFromSheet(p)), name);
  }

  /** Generates [count] increasing integer numbers, starting with 0. */
  get wells(): IterableIterator<PlateWell> {
    const self = this; // Capture this context
    return (function* () {
      for (let row = 0; row < self.rows; row++) {
        for (let col = 0; col < self.cols; col++) {
          const well: PlateWell = {
            row: row,
            col: col,
          };
          for (const field of self.data.columns) {
            if (!field.isNone(self._idx(row, col))) {
              well[field.name] = field.get(self._idx(row, col));
            }
          }
          yield well;
        }
      }
    })();
  }

  /** Returns the specified field value at the specified positions (0-based).
   * The following examples all return volume for the C4 well:
   * - plate.get('volume', 2, 3)
   * - plate.get('volume', 'C4')
   * */
  get(field: string, rowIdxOrPos: number | string, colIdx?: number): any {
    assure(this.data.columns.byName(field) != null, `Field does not exist: ${field}`);
    if (typeof rowIdxOrPos === 'number' && !colIdx)
      throw 'Column not defined';

    const [row, col] = (typeof rowIdxOrPos === 'string') ? parseExcelPosition(rowIdxOrPos) : [rowIdxOrPos, colIdx!];
    return this.data.columns.byName(field).get(this._idx(row, col));
  }

  /** Changes all numerical values to the results of the specified normalization function */
  normalize(field: string, f: (value: number) => number, inplace: boolean = false) {
    const originalCol = this.data.getCol(field);
    const col = inplace ? originalCol : this.data.columns.addNewFloat(this.data.columns.getUnusedName(`${field}_normalized`));
    col.init((i) => originalCol.isNone(i) ? null : f(originalCol.get(i)));
    return col;
  }

  matches(i: number, filter: IPlateWellFilter): boolean {
    const cols = this.data.columns;
    // we allow both any and array of any things in matches object, so we need to convert it to array to have one api
    const arMatches = filter.match ? Object.entries(filter.match).reduce((acc, [key, value]) => {acc[key] = Array.isArray(value) ? value : [value]; return acc;}, {} as Record<string, any[]>) : null;
    const arExclude = filter.exclude ? Object.entries(filter.exclude).reduce((acc, [key, value]) => {acc[key] = Array.isArray(value) ? value : [value]; return acc;}, {} as Record<string, any[]>) : null;
    return !filter || (!arMatches && !arExclude) ||
        (
        Object.keys(arMatches ?? {}).filter((m) => cols.contains(m)).every((key) => arMatches![key].some((m) => this.data.columns.byName(key).get(i) === m)) &&
        Object.keys(arExclude ?? {}).filter((m) => cols.contains(m)).every((key) => arExclude![key].every((e) => this.data.columns.byName(key).get(i) !== e)));
  }

  count(filter?: IPlateWellFilter): number {
    if (!filter)
      return this.rows * this.cols;
    let c = 0;
    for (let i = 0; i < this.rows * this.cols; i++)
      if (this.matches(i, filter))
        c++;
    return c;
  }

  /** Array of non-empty values */
  values(fields: string[], filter?: IPlateWellFilter): Array<Record<string, any> & {innerDfRow: number}> {
    const cols = fields.map((f) => this.data.columns.byName(f));
    assure(cols.every((c) => c != null), `Field does not exist: ${fields.find((_, i) => cols[i] == null)}`);
    const colsObj: Record<string, DG.Column> = {};
    for (let i = 0; i < fields.length; i++)
      colsObj[fields[i]] = cols[i];
    const result: (Record<string, any> & {innerDfRow: number}) []  = [];
    for (let i = 0; i < this.rows * this.cols; i++) {
      if (((filter?.includeEmpty ?? false) || cols.every((col) => !col.isNone(i))) && (!filter || this.matches(i, filter))) {
        const res = fields.reduce((acc, f) => { acc[f] = colsObj[f].isNone(i) ? null : colsObj[f].get(i); return acc; }, {} as Record<string, any> & {innerDfRow: number});
        res.innerDfRow = i;
        result.push(res);
      }
    }
    return result;
  }

  /// Returns specified field statistics for the specified filter and specified field name
  getStatistics<T extends JSTATStatistics[]>(field: string, statistics: T, filter?: IPlateWellFilter): Record<T[number], number> {
    const values = this.fieldValues(field, filter);
    const result: Record<JSTATStatistics, number> = {} as Record<JSTATStatistics, number>;
    for (const stat of statistics) {
      result[stat] = jstatStatistics[stat](values);
    }
    return result;
  }

  fieldValues(field: string, filter?: IPlateWellFilter): Array<any> {
    return this.values([field], filter).map((v) => v[field]);
  }

  static inspectSeriesByName(records: Record<string, FitSeries>, seriesName: string, fitFunctionName: FitFunctionType): void {
    if (!seriesName || !fitFunctionName)
      return;
    Plate.inspectSeries(records[seriesName], fitFunctionName);
  }

  static inspectSeries(series: FitSeries, fitFunctionName: FitFunctionType) {
    const seriesName = series?.name ?? 'Series';
    if (!series || !fitFunctionName)
      return;
    const curveCol = DG.Column.string('Curve', 1);
    curveCol.set(0, JSON.stringify({
      chartOptions: {
        logX: true,
        title: seriesName,
      },
      series: [{...series, fit: undefined, fitFunction: fitFunctionName, clickToToggle: true, droplines: ['IC50'], name: seriesName}]
    }), false);
    const df = DG.DataFrame.fromColumns([curveCol]);
    df.name = seriesName;
    df.id = randomizeTableId();
    curveCol.semType = 'fit';
    curveCol.tags[DG.Tags.CellRenderer] = 'fit';
    const grid = DG.Viewer.grid(df);
    const gridCell = grid.cell('Curve', 0);
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = gridCell;
    inspectCurve(gridCell, {width: 480, height: 370}, true);
  }

  doseResponseSeries(options?: IPlateWellFilter & { concentration?: string; value?: string, groupBy?: string}): Record<string, FitSeries> {
    const valueOptions = {includeEmpty: options?.includeEmpty ?? false, exclude: options?.exclude ?? {'role': ['High Control', 'Low Control']}}
    const concKey = options?.concentration ?? 'concentration';
    const valueKey = options?.value ?? 'value';
    const values = this.values([concKey, valueKey, ...(options?.groupBy ? [options.groupBy] : [])], valueOptions);

    const series: Record<string, ISeriesData & Record<string, any>> = {};
    for (const v of values) {
      const group = options?.groupBy ? v[options.groupBy] : '0';
      if (!series[group])
        series[group] = {x: [], y: [], meta: [], outlier: []};
      series[group].x.push(v[concKey]);
      series[group].y.push(v[valueKey]);
      series[group].meta.push(v.innerDfRow);
      series[group].outlier.push(this._isOutlier(v.innerDfRow));
    }
    return Object.fromEntries(Object.entries(series).map(([k, v]) => {
      const fitSeries = new FitSeries(v.x.map((_, i) => ({x: v.x[i], y: v.y[i], outlier: v.outlier[i], meta: v.meta[i]})).sort((a, b) => a.x - b.x));
      fitSeries.name = k;
      return [k, fitSeries];
    }));
  }

  getAnalysisDialog(options: AnalysisOptions) {
    ui.dialog('Plate Analysis')
    .add(PlateWidget.analysisView(this, options))
    .showModal(true);
  }

  getAnalysisView(options: AnalysisOptions) {
    const view = DG.View.fromRoot(PlateWidget.analysisView(this, options).root);
    view.name = 'Plate Analysis';
    return grok.shell.addView(view);
  }

  /** Adds data from the other plate to this plate. Typically, you would apply plate layout */
  merge(plate: Plate) {
    for (const col of plate.data.columns) {
      this.data.columns.add(col.clone());
    }
    return this;
  }

  /** Deep cloning */
  clone(): Plate {
    const cloned = new Plate(this.rows, this.cols);
    cloned.data = this.data.clone();
    return cloned;
  }

  print() {
    for (const fieldCol of this.data.columns) {
      console.log();
      console.log(fieldCol.name);
      console.log(DG.range(this.cols + 1).map((col) => col == 0 ? '' : `${col}`).toArray().join('\t'));
      for (let row = 0; row < this.rows; row++) {
        console.log(numToExcel(row) + '\t' +
          DG.range(this.cols).map((col) => fieldCol.getString(this._idx(row, col))).toArray().join(',\t'));
      }
    }
  }
}





