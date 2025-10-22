/* eslint-disable max-len */
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
  parseExcelPosition, standardPlateSizes, toExcelPosition, toStandardSize
} from './utils';
import type ExcelJS from 'exceljs';
import {findPlatePositions, getPlateFromSheet} from './excel-plates';
import {plateDbColumn, allProperties, plateTypes} from '../plates/plates-crud';
import {Subject} from 'rxjs';
import { IPlateWellValidator } from './plate-well-validators';

/** Represents a well in the experimental plate */
export interface PlateWell {
  row: number;
  col: number;
  tags?: string[]; // New addition for future use
  [key: string]: any;
}

interface IPlateCsvImportOptions {
  rowColumn?: string;
  colColumn?: string;
  positionColumn?: string;
  volumeColumn?: string;
  concentrationColumn?: string;
  field?: string; // used as a field name when a 12x8 table is passed; if not passed, defaults to 'value'
}

export interface IPlateWellFilter {
  includeEmpty?: boolean;
  match?: {[key: string]: any};
  exclude?: {[key: string]: any};
}

export const PLATE_OUTLIER_WELL_NAME = 'Outlier';

/** Represents experimental plate (typically 96-, 384-, or 1536-well assay plates) */
export class Plate {
  id?: number;
  plateTemplateId?: number;
  plateTypeId?: number;
  barcode?: string;

  data: DG.DataFrame;
  details: {[index: string]: any} = {}; // dynamic plate properties, stored in db: plates.plate_details
  rows: number = 8;
  cols: number = 12;


  private outlierChangeSubject = new Subject<{row: number, col: number, isOutlier: boolean, source: string}>();
  public get onOutlierChanged() { return this.outlierChangeSubject.asObservable(); }

  private dataChangeSubject = new Subject<void>();
  public get onDataChanged() { return this.dataChangeSubject.asObservable(); }

  public fireDataChange(): void {
    this.dataChangeSubject.next();
  }


  constructor(rows: number, cols: number) {
    this.data = DG.DataFrame.create(rows * cols);
    this.rows = rows;
    this.cols = cols;
    this.plateTypeId = plateTypes.find((pt) => pt.rows == this.rows && pt.cols == this.cols)?.id;
    this.barcode = Math.round(Math.random() * 10000).toString().padStart(10, '0');
  }

  static autoSize(positions: Iterable<[number, number]>): Plate {
    const [rows, cols] = toStandardSize(getMaxPosition(positions));
    return new Plate(rows, cols);
  }

  _idx(row: number, col: number): number {
    return row * this.cols + col;
  }

  rowIndexToExcel(dataFrameRow: number): [row: number, col: number] {
    return Plate._idxToPos(dataFrameRow, this.cols);
  }

  static _idxToPos(dataFrameRow: number, cols: number): [row: number, col: number] {
    return [Math.floor(dataFrameRow / cols), dataFrameRow % cols];
  }

  private _markOutlier(dataframeRow: number, flag: boolean = true) {
    const outlierCol = this.data.columns.getOrCreate(PLATE_OUTLIER_WELL_NAME, DG.TYPE.BOOL);
    outlierCol.set(dataframeRow, flag);
    const [row, col] = this.rowIndexToExcel(dataframeRow);
    this.outlierChangeSubject.next({row, col, isOutlier: flag, source: ''});
  }

  markOutlier(row: number, col: number, flag: boolean = true) {
    this._markOutlier(this._idx(row, col), flag);
  }

  _isOutlier(dataframeRow: number): boolean {
    const outlierCol = this.data.col(PLATE_OUTLIER_WELL_NAME);
    return outlierCol ? outlierCol.get(dataframeRow) : false;
  }

  isOutlier(row: number, col: number): boolean {
    return this._isOutlier(this._idx(row, col));
  }


  markOutliersWhere(field: string, valueFunc: (fieldValue: any) => boolean, filter?: IPlateWellFilter) {
    this.values([field], filter).forEach((v) => {
      if (valueFunc(v[field]))
        this._markOutlier(v.innerDfRow, true);
    });
  }

  /**
   * Get a column by its name or any of its aliases
   */
  getColumn(name: string): DG.Column | null {
    if (this.data.columns.contains(name))
      return this.data.col(name);
    return null;
  }


  static fromGridDataFrame(table: DG.DataFrame, field?: string): Plate {
    const rows = table.rowCount;
    function containsRowLetters(c: DG.Column): boolean {
      return c.type == DG.TYPE.STRING && DG.range(rows).every((i) => numToExcel(i) == c.get(i).toUpperCase());
    }
    const dataColumns = wu(table.columns).filter((c) => !containsRowLetters(c)).toArray();
    const cols = dataColumns.length;
    const plate = new Plate(rows, cols);
    const nonEmptyCols = dataColumns.filter((c) => c.stats.missingValueCount != rows);
    const toString = (nonEmptyCols.length > 0 && nonEmptyCols.findIndex((c) => c.type != nonEmptyCols[0].type) != -1);
    const type = toString ?
      DG.TYPE.STRING :
      firstWhere(dataColumns, (c) => c.stats.missingValueCount != rows)?.type ?? DG.TYPE.STRING;
    const dataColumn = plate.data.columns.addNew(table.name ?? field, type);
    for (let colIdx = 0; colIdx < dataColumns.length; colIdx ++) {
      if (dataColumns[colIdx].stats.missingValueCount != rows) {
        for (let i = 0; i < rows; i++)
          dataColumn.set(plate._idx(i, colIdx), type != DG.TYPE.STRING ? dataColumns[colIdx].get(i) : dataColumns[colIdx].getString(i), false);
      }
    }
    return plate;
  }
  toGridDataFrame(layer: string): DG.DataFrame {
    const df = DG.DataFrame.create(this.rows, layer);
    const col = this.getColumn(layer) || this.data.columns.byName(layer);

    for (let i = 0; i < this.cols; i++)
      df.columns.addNew(`${i + 1}`, col.type).init((r) => col.get(this._idx(r, i)));

    return df;
  }

  static fromTableByRow(table: DG.DataFrame, options?: {posColName?: string, rowColName?: string, colColName?: string}): Plate {
    const posCol = wu(table.columns).find((c) =>
      (!options?.posColName || c.name.toLowerCase() == options.posColName) &&
    c.type == DG.TYPE.STRING && /^([A-Za-z]+)(\d+)$/.test(c.get(0)));
    const rowColName = options?.rowColName ?? 'row';
    const colColName = options?.colColName ?? 'col';
    const rowCol = table.col(rowColName)?.type === DG.TYPE.INT ? table.col(rowColName) : null;
    const colCol = table.col(colColName)?.type === DG.TYPE.INT ? table.col(colColName) : null;
    const rowToPos: ((i: number) => [row: number, col: number]) | null = posCol ?
      ((i) => parseExcelPosition(posCol!.get(i))) :
      rowCol && colCol ? ((i) => [rowCol!.get(i), colCol?.get(i)]) :
        standardPlateSizes[table.rowCount] ? ((i) => Plate._idxToPos(i, standardPlateSizes[table.rowCount][1])) :
          null;
    if (rowToPos == null)
      throw new Error('Columns with well positions not identified');
    const positions = DG.range(table.rowCount).map(rowToPos);
    const [rows, cols] = toStandardSize(getMaxPosition(positions));
    const plate = new Plate(rows, cols);
    for (const col of table.columns) {
      if (col.isEmpty || col == posCol || col == rowCol || col == colCol)
        continue;
      const plateCol = plate.data.columns.addNew(col.name, col.type);
      for (let r = 0; r < col.length; r++) {
        const [wellRow, wellCol] = rowToPos!(r);
        plateCol.set(plate._idx(wellRow, wellCol), col.get(r), false);
      }
    }
    return plate;
  }

  static fromDbDataFrame(df: DG.DataFrame): Plate {
    const plate = new Plate(df.col('row')!.stats.max + 1, df.col('col')!.stats.max + 1);
    const pidCol = df.col('property_id');
    for (let propBlock = 0; propBlock < df.rowCount / (plate.rows * plate.cols); propBlock++) {
      const start = propBlock * plate.rows * plate.cols;
      const pid = pidCol?.get(start);
      const property = allProperties.find((p) => p.id == pid)!;
      const valueColumn = df.col(plateDbColumn[property.type])!;
      //@ts-ignore
      const newCol = plate.data.columns.addNew(property.name, property.type)
        .init((i) => valueColumn.get(start + i));
    }

    return plate;
  }

  static async fromCsvTableFile(csvPath: string, field: string, options?: DG.CsvImportOptions): Promise<Plate> {
    return this.fromGridDataFrame(await grok.dapi.files.readCsv(csvPath, options), field);
  }

  static fromCsvTable(csv: string, field: string, options?: DG.CsvImportOptions): Plate {
    return this.fromGridDataFrame(DG.DataFrame.fromCsv(csv, options), field);
  }

  static fromCsv(csv: string, options?: IPlateCsvImportOptions): Plate {
    const df = DG.DataFrame.fromCsv(csv);
    if (DG.range(df.rowCount).every((i) => df.col(0)!.get(i) == String.fromCharCode(65 + i)))
      df.columns.remove(df.columns.byIndex(0).name);
    if (df.columns.length == 13 || df.columns.length == 25 || df.columns.length == 49 && df.col(df.columns.length - 1)?.isEmpty)
      df.columns.remove(df.columns.length - 1);
    return this.fromGridDataFrame(df, options?.field ?? 'value');
  }

  static fromPlates(plates: Plate[], name?: string) {
    assure(plates.length > 0, 'Array is empty.');
    assure(plates.every((p) => p.rows == plates[0].rows && p.cols == plates[0].cols), 'Plate dimensions differ.');
    const result = plates.reduce((p1, p2) => p1 ? p1.merge(p2) : p2.clone());
    if (name != null)
      result.data.name = name;
    return result;
  }

  static async fromExcelFile(fi: DG.FileInfo): Promise<Plate> {
    return fi.data && fi.data.length ? await Plate.fromExcel(fi.data, fi.friendlyName) : await Plate.fromExcelPath(fi.fullPath, fi.friendlyName);
  }

  static async fromExcelPath(excelPath: string, name?: string): Promise<Plate> {
    const content = await grok.dapi.files.readAsBytes(excelPath);
    return this.fromExcel(content, name);
  }

  static async fromExcel(excelBytes: Uint8Array, name?: string): Promise<Plate> {
    await DG.Utils.loadJsCss(['/js/common/exceljs.min.js']);
    //@ts-ignore
    const loadedExcelJS = window.ExcelJS as ExcelJS;
    const workbook = new loadedExcelJS.Workbook() as ExcelJS.Workbook;
    //@ts-ignore
    const wb = await workbook.xlsx.load(excelBytes);

    const platePositions = findPlatePositions(wb);
    if (platePositions.length == 0) throw new Error('Plates not found in excel file');
    const p0 = platePositions[0];
    if (!platePositions.every((pc) => pc.rows == p0.rows || pc.cols == p0.cols))
      throw new Error(`Plate sizes differ in "${name}"`);
    return Plate.fromPlates(platePositions.map((p) => getPlateFromSheet(p)), name);
  }

  static demo(): Plate {
    const plate = Plate.fromTableByRow(grok.data.demo.wells(96));
    plate.barcode = 'demo-' + Math.round(Math.random() * 10000).toString().padStart(6, '0');
    return plate;
  }

  get wells(): wu.WuIterable<PlateWell> {
    const self = this;
    return wu((function* () {
      for (let row = 0; row < self.rows; row++) {
        for (let col = 0; col < self.cols; col++) {
          const well: PlateWell = {row: row, col: col};
          for (const field of self.data.columns) {
            if (field.name.toLowerCase() != 'row' && field.name.toLowerCase() != 'col' && !field.isNone(self._idx(row, col)))
              well[field.name] = field.get(self._idx(row, col));
          }
          yield well;
        }
      }
    })());
  }

  get(field: string, rowIdxOrPos: number | string, colIdx?: number): any {
    const column = this.getColumn(field);
    assure(column != null, `Field does not exist: ${field}`);
    if (typeof rowIdxOrPos === 'number' && colIdx == null)
      throw new Error('Column not defined');
    const [row, col] = (typeof rowIdxOrPos === 'string') ? parseExcelPosition(rowIdxOrPos) : [rowIdxOrPos, colIdx!];
    return column!.get(this._idx(row, col));
  }

  normalize(field: string, f: (value: number) => number, inplace: boolean = false) {
    const originalCol = this.getColumn(field) || this.data.getCol(field);
    const col = inplace ? originalCol : this.data.columns.addNewFloat(this.data.columns.getUnusedName(`${field}_normalized`));
    col.init((i) => originalCol.isNone(i) ? null : f(originalCol.get(i)));
    return col;
  }

  matches(i: number, filter: IPlateWellFilter): boolean {
    const cols = this.data.columns;
    const arMatches = filter.match ? Object.entries(filter.match).reduce((acc, [key, value]) => {
      acc[key] = Array.isArray(value) ? value : [value]; return acc;
    }, {} as Record<string, any[]>) : null;
    const arExclude = filter.exclude ? Object.entries(filter.exclude).reduce((acc, [key, value]) => {
      acc[key] = Array.isArray(value) ? value : [value]; return acc;
    }, {} as Record<string, any[]>) : null;
    return !filter || (!arMatches && !arExclude) ||
      (
        Object.keys(arMatches ?? {}).filter((m) => cols.contains(m)).every((key) => arMatches![key].some((m) => this.data.columns.byName(key).get(i) === m)) &&
        Object.keys(arExclude ?? {}).filter((m) => cols.contains(m)).every((key) => arExclude![key].every((e) => this.data.columns.byName(key).get(i) !== e)));
  }

  count(filter?: IPlateWellFilter): number {
    if (!filter)
      return this.rows * this.cols;
    let c = 0;
    for (let i = 0; i < this.rows * this.cols; i++) {
      if (this.matches(i, filter))
        c++;
    }
    return c;
  }

  values(fields: string[], filter?: IPlateWellFilter): Array<Record<string, any> & {innerDfRow: number}> {
    const cols = fields.map((f) => this.getColumn(f) || this.data.columns.byName(f));
    assure(cols.every((c) => c != null), `Field does not exist: ${fields.find((_, i) => cols[i] == null)}`);
    const colsObj: Record<string, DG.Column> = {};
    for (let i = 0; i < fields.length; i++)
      colsObj[fields[i]] = cols[i];
    const result: (Record<string, any> & {innerDfRow: number}) [] = [];
    for (let i = 0; i < this.rows * this.cols; i++) {
      if (((filter?.includeEmpty ?? false) || cols.every((col) => !col.isNone(i))) && (!filter || this.matches(i, filter))) {
        const res = fields.reduce((acc, f) => {
          acc[f] = colsObj[f].isNone(i) ? null : colsObj[f].get(i); return acc;
        }, {} as Record<string, any> & {innerDfRow: number});
        res.innerDfRow = i;
        result.push(res);
      }
    }
    return result;
  }

  getStatistics<T extends JSTATStatistics[]>(field: string, statistics: T, filter?: IPlateWellFilter): Record<T[number], number> {
    const values = this.fieldValues(field, filter);
    const result: Record<JSTATStatistics, number> = {} as Record<JSTATStatistics, number>;
    for (const stat of statistics)
      result[stat] = jstatStatistics[stat](values);
    return result;
  }

  fieldValues(field: string, filter?: IPlateWellFilter): Array<any> {
    return this.values([field], filter).map((v) => v[field]);
  }

  getLayerNames(): string[] {
    return this.data.columns.names();
  }

  merge(plate: Plate) {
    for (const col of plate.data.columns)
      this.data.columns.add(col.clone());
    return this;
  }

  clone(): Plate {
    const cloned = new Plate(this.rows, this.cols);
    cloned.data = this.data.clone();
    cloned.details = {...this.details};
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

  validateWells(validators: IPlateWellValidator[]): Map<string, string[]> {
    const result = new Map<string, string[]>();
    for (const validator of validators) {
      for (let row = 0; row < this.rows; row++) {
        for (let col = 0; col < this.cols; col++) {
          const error = validator.validate(this, row, col);
          if (error) {
            const errors = result.get(`${numToExcel(row)}${col}`) ?? [];
            errors.push(error);
            result.set(`${toExcelPosition(row, col)}`, errors);
          }
        }
      }
    }
    return result;
  }
}
