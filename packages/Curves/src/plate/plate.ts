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
import {FitFunctionType, FitSeries} from '@datagrok-libraries/statistics/src/fit/new-fit-API';
import {AnalysisOptions} from './plate-widget';
import {inspectCurve} from '../fit/fit-renderer';
import {plateDbColumn, wellProperties, plateTypes} from '../plates/plates-crud';
import {PlateDrcAnalysis} from './plate-drc-analysis';
import {IPlateWellValidator} from './plate-well-validators';
import {Subject} from 'rxjs';

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

interface ISeriesData {
  x: number[];
  y: number[];
}

export enum LayerType {
  ORIGINAL = 'original',
  DERIVED = 'derived',
  LAYOUT = 'layout',
  OUTLIER = 'outlier'
}


export interface LayerMetadata {
  type: LayerType;
  source?: string;
  createdAt: Date;
}

export const PLATE_OUTLIER_WELL_NAME = 'Outlier';

export function randomizeTableId() {
  return `${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}-${Math.random()}`;
}

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

  private layerRegistry: Map<string, LayerMetadata> = new Map();
  private scopedAliases: Map<string, Map<string, string>> = new Map();

  private outlierChangeSubject = new Subject<{row: number, col: number, isOutlier: boolean, source: string}>();
  public get onOutlierChanged() { return this.outlierChangeSubject.asObservable(); }


  private changeLog: Array<{
    timestamp: Date;
    action: string;
    details: any;
  }> = [];

  constructor(rows: number, cols: number) {
    this.data = DG.DataFrame.create(rows * cols);
    this.rows = rows;
    this.cols = cols;
    this.plateTypeId = plateTypes.find((pt) => pt.rows == this.rows && pt.cols == this.cols)?.id;
    this.barcode = Math.round(Math.random() * 10000).toString().padStart(10, '0');
    this.logChange('plate-created', {rows, cols});
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

  _markOutlier(row: number, flag: boolean = true) {
    const outlierCol = this.data.columns.getOrCreate(PLATE_OUTLIER_WELL_NAME, DG.TYPE.BOOL);

    if (!this.layerRegistry.has(PLATE_OUTLIER_WELL_NAME)) {
      this.registerLayer(PLATE_OUTLIER_WELL_NAME, LayerType.OUTLIER, 'user-interaction');
      outlierCol.setTag('outlier.sources', '{}');
    }

    outlierCol.set(row, flag);
    const source = 'user-interaction';
    this._updateOutlierSource(outlierCol, row, flag, source);
    this._markOutlierWithSource(row, flag, source);
  }

  markOutlierWithSource(row: number, col: number, flag: boolean = true, source: string = 'user-interaction') {
    const dataRow = this._idx(row, col);
    this._markOutlierWithSource(dataRow, flag, source);
  }
  _markOutlierWithSource(row: number, flag: boolean = true, source: string = 'user-interaction') {
    const outlierCol = this.data.columns.getOrCreate(PLATE_OUTLIER_WELL_NAME, DG.TYPE.BOOL);


    // Register the outlier column as an OUTLIER layer when first created
    if (!this.layerRegistry.has(PLATE_OUTLIER_WELL_NAME)) {
      this.registerLayer(PLATE_OUTLIER_WELL_NAME, LayerType.OUTLIER, 'user-interaction');
      outlierCol.setTag('outlier.sources', '{}');
    }

    outlierCol.set(row, flag);
    this._updateOutlierSource(outlierCol, row, flag, source);
    // Emitting change event
    const [plateRow, plateCol] = this.rowIndexToExcel(row);
    this.outlierChangeSubject.next({row: plateRow, col: plateCol, isOutlier: flag, source});
  }

  // Helper method to manage source tracking via column tags
  private _updateOutlierSource(outlierCol: DG.Column, row: number, flag: boolean, source: string) {
    try {
      const sourcesJson = outlierCol.getTag('outlier.sources') || '{}';
      const sources: Record<number, string[]> = JSON.parse(sourcesJson);

      if (!sources[row])
        sources[row] = [];


      if (flag && !sources[row].includes(source)) {
        sources[row].push(source);
      } else if (!flag) {
      // Clear sources when unmarking outlier
        sources[row] = [];
      }

      outlierCol.setTag('outlier.sources', JSON.stringify(sources));
    } catch (e) {
      console.warn('Failed to update outlier sources:', e);
      // Fallback: just ensure the tag exists
      outlierCol.setTag('outlier.sources', '{}');
    }
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

  markOutliersWhere(field: string, valueFunc: (fieldValue: any) => boolean, filter?: IPlateWellFilter, source: string = 'auto-filter') {
    this.values([field], filter).forEach((v) => {
      if (valueFunc(v[field]))
        this._markOutlierWithSource(v.innerDfRow, true, source);
    });
  }

  /**
   * Get a column by its name or any of its aliases
   */
  getColumn(nameOrAlias: string): DG.Column | null {
  // First check direct column name
    if (this.data.columns.contains(nameOrAlias))
      return this.data.col(nameOrAlias);

    // Then search all scopes for this alias
    for (const scopeMap of this.scopedAliases.values()) {
      if (scopeMap.has(nameOrAlias))
        return this.data.col(scopeMap.get(nameOrAlias)!);
    }

    return null;
  }
  addScopedAlias(scope: string, originalColumnName: string, alias: string): void {
    if (!this.data.columns.contains(originalColumnName))
      throw new Error(`Column ${originalColumnName} does not exist`);

    if (!this.scopedAliases.has(scope))
      this.scopedAliases.set(scope, new Map());
    // Remove this alias from any other column in this scope first
    this.removeScopedAlias(scope, alias);

  this.scopedAliases.get(scope)!.set(alias, originalColumnName);
  this.logChange('scoped-alias-added', {originalColumnName, scope, alias});
  }

  removeScopedAlias(scope: string, alias: string): void {
    const scopeMap = this.scopedAliases.get(scope);
    if (scopeMap && scopeMap.has(alias)) {
      const originalColumnName = scopeMap.get(alias)!;
      scopeMap.delete(alias);
      this.logChange('scoped-alias-removed', {originalColumnName, scope, alias});
    }
  }

  getScopedAliases(scope: string): Map<string, string> {
    return this.scopedAliases.get(scope) || new Map();
  }

  getScopedColumn(scope: string, nameOrAlias: string): DG.Column | null {
  // First try direct column name
    if (this.data.columns.contains(nameOrAlias))
      return this.data.col(nameOrAlias);

    // Then check scoped alias
    const scopeMap = this.scopedAliases.get(scope);
    if (scopeMap && scopeMap.has(nameOrAlias)) {
      const originalColumnName = scopeMap.get(nameOrAlias)!;
      return this.data.col(originalColumnName);
    }

    return null;
  }

  /**
   * Register a column as a layer with metadata
   */
  registerLayer(columnName: string, type: LayerType, source?: string): void {
    if (!this.data.columns.contains(columnName))
      throw new Error(`Column ${columnName} does not exist`);


    if (!this.layerRegistry.has(columnName)) {
      this.layerRegistry.set(columnName, {
        type,
        source: source || 'unknown',
        createdAt: new Date(),
      });

      this.logChange('layer-registered', {columnName, type, source});
    }
  }


  getLayerType(columnName: string): LayerType | undefined {
    return this.layerRegistry.get(columnName)?.type;
  }
  getOutlierSources(row: number, col: number): string[] {
    const dataRow = this._idx(row, col);
    const outlierCol = this.data.col(PLATE_OUTLIER_WELL_NAME);
    if (!outlierCol) return [];

    try {
      const sourcesJson = outlierCol.getTag('outlier.sources') || '{}';
      const sources: Record<number, string[]> = JSON.parse(sourcesJson);
      return sources[dataRow] || [];
    } catch (e) {
      console.warn('Failed to parse outlier sources:', e);
      return [];
    }
  }
  getAllOutlierSources(): Record<number, string[]> {
    const outlierCol = this.data.col(PLATE_OUTLIER_WELL_NAME);
    if (!outlierCol) return {};

    try {
      const sourcesJson = outlierCol.getTag('outlier.sources') || '{}';
      return JSON.parse(sourcesJson);
    } catch (e) {
      console.warn('Failed to parse outlier sources:', e);
      return {};
    }
  }

  getLayersByType(type: LayerType): string[] {
    const layers: string[] = [];
    for (const [name, metadata] of this.layerRegistry) {
      if (metadata.type === type)
        layers.push(name);
    }
    return layers;
  }

  private logChange(action: string, details: any): void {
    this.changeLog.push({
      timestamp: new Date(),
      action,
      details
    });
  }

  getChangeLog(): ReadonlyArray<any> {
    return this.changeLog;
  }

  // --- EXISTING METHODS WITH ALIAS SUPPORT ---

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

    // Register imported column as ORIGINAL
    plate.registerLayer(dataColumn.name, LayerType.ORIGINAL, 'grid-import');

    return plate;
  }

  toGridDataFrame(layer: string): DG.DataFrame {
    const df = DG.DataFrame.create(this.rows, layer);
    const col = this.getColumn(layer) || this.data.columns.byName(layer);

    for (let i = 0; i < this.cols; i++)
      df.columns.addNew(`${i + 1}`, col.type).init((r) => col.get(this._idx(r, i)));

    return df;
  }

  static fromTableByRow(
    table: DG.DataFrame, options?: {
      posColName?: string,
      rowColName?: string,
      colColName?: string
    }): Plate {
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

      // Register imported columns as ORIGINAL
      plate.registerLayer(col.name, LayerType.ORIGINAL, 'table-by-row');
    }
    return plate;
  }

  static fromDbDataFrame(df: DG.DataFrame): Plate {
    const plate = new Plate(df.col('row')!.stats.max + 1, df.col('col')!.stats.max + 1);
    const pidCol = df.col('property_id');
    for (let propBlock = 0; propBlock < df.rowCount / (plate.rows * plate.cols); propBlock++) {
      const start = propBlock * plate.rows * plate.cols;
      const pid = pidCol?.get(start);
      const property = wellProperties.find((p) => p.id == pid)!;
      const valueColumn = df.col(plateDbColumn[property.type])!;
      //@ts-ignore
      const newCol = plate.data.columns.addNew(property.name, property.type)
        .init((i) => valueColumn.get(start + i));

      // Register database columns as ORIGINAL
      plate.registerLayer(property.name, LayerType.ORIGINAL, 'database');
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

    // Register normalized column as DERIVED
    if (!inplace)
      this.registerLayer(col.name, LayerType.DERIVED, 'normalization');


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
      series: [{...series, fit: undefined, fitFunction: fitFunctionName, clickToToggle: true, droplines: ['IC50'], name: seriesName}],
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
    const valueOptions = {includeEmpty: options?.includeEmpty ?? false, exclude: options?.exclude ?? {'role': ['High Control', 'Low Control']}};
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

      const isOutlier = this._isOutlier(v.innerDfRow);
      console.log(`[DEBUG] doseResponseSeries: Row ${v.innerDfRow} isOutlier=${isOutlier}`);
      series[group].outlier.push(isOutlier);
    }

    console.log('--- Debugging doseResponseSeries ---');
    console.log(`Grouping by column: '${options?.groupBy ?? 'undefined'}' on plate with barcode: '${this.barcode}'`);
    const seriesPointCounts = Object.fromEntries(
      Object.entries(series).map(([key, value]) => [key, value.x.length])
    );
    console.log('Found series and their point counts:', seriesPointCounts);
    console.log('Full series data for inspection:', series);

    return Object.fromEntries(Object.entries(series).map(([k, v]) => {
      const fitSeries = new FitSeries(v.x.map((_, i) => ({x: v.x[i], y: v.y[i], outlier: v.outlier[i], meta: v.meta[i]})).sort((a, b) => a.x - b.x));
      fitSeries.name = k;
      return [k, fitSeries];
    }));
  }

  merge(plate: Plate) {
    for (const col of plate.data.columns) {
      this.data.columns.add(col.clone());
      // Register merged columns, preserving their type if known
      const sourceType = plate.getLayerType(col.name);
      if (sourceType)
        this.registerLayer(col.name, sourceType, 'merged');
    }
    return this;
  }

  clone(): Plate {
    const cloned = new Plate(this.rows, this.cols);
    cloned.data = this.data.clone();
    cloned.layerRegistry = new Map(this.layerRegistry);
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

//   validateWells(validators: IPlateWellValidator[]): Map<string, string[]> {
//     const result = new Map<string, string[]>();
//     for (const validator of validators) {
//       for (let row = 0; row < this.rows; row++) {
//         for (let col = 0; col < this.cols; col++) {
//           const error = validator.validate(this, row, col);
//           if (error) {
//             const errors = result.get(`${numToExcel(row)}${col}`) ?? [];
//             errors.push(error);
//             result.set(`${toExcelPosition(row, col)}`, errors);
//           }
//         }
//       }
//     }
//     return result;
//   }
}
