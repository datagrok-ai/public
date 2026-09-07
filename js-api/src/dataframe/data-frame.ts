/**
 * DataFrame - the core tabular data structure.
 * @module dataframe/data-frame
 */

import * as rxjs from 'rxjs';
import {
  TAGS,
  JoinType,
  ColumnType,
  CsvImportOptions,
  ViewerType,
  JOIN_TYPE
} from "../const";
import {__obs, EventData, MapChangeArgs, CellRangeArgs, RowChangeArgs, ColumnChangeArgs, CellChangeArgs} from "../events";
import {toDart, toJs} from "../wrappers";
import {MapBag, MapProxy} from "../proxies";
import {_toJson} from "../utils_convert";
import {Observable} from "rxjs";
import type {Widget} from '../widgets';
import {IDartApi} from "../api/grok_api.g";
import type {Property, TableInfo} from "../entities";
import type {Grid, FormViewer} from "../grid";
import type {BarChartViewer, BoxPlot, HistogramViewer, LineChartViewer, NetworkDiagramViewer, ScatterPlotViewer, TileViewer, Viewer, ViewerOptions} from "../viewer";
import type {IBarChartSettings, IBoxPlotSettings, IFormSettings, IGridSettings, IHistogramSettings, ILineChartSettings, INetworkDiagramSettings, IScatterPlotSettings, ITileViewerSettings} from "../interfaces/d4";
import {BitSet} from "./bit-set";
import {ColumnList} from "./column-list";
import {Row, RowList, Cell} from "./row";
import {Column} from "./column";
import {GroupByBuilder} from "./stats";
import type {CloneOptions, CsvExportOptions, ColumnId, GroupsDescription, JoinOptions} from "./types";
import {DataFrameFormulaLinesHelper, DataFrameAnnotationRegionsHelper} from "./formula-helpers";

declare let grok: any;
declare let DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 *
 * In the API, the terms "Table" and "DataFrame" are used interchangeably.
 *
 * Usage samples: {@link https://public.datagrok.ai/js/samples/data-frame/modification/manipulate}
 * Usage details: {@link https://datagrok.ai/help/develop/advanced/data-frame}
 * Implementation details: {@link https://datagrok.ai/help/develop/admin/architecture#in-memory-database}
 */
export class DataFrame {
  public readonly dart: any;
  /** Columns of the table. */
  public columns: ColumnList;
  /** Rows of the table. */
  public rows: RowList;
  /** Filter mask: rows currently passing all filters. */
  public filter: BitSet;
  /** Auxiliary data that is not persisted (a {@link MapBag}: indexed access plus the map methods; typed `any` so it can be cast to a package's own shape). */
  public temp: any;
  /** Metadata as string key-value pairs; persisted with the table. */
  public tags: MapBag<string>;
  public _meta: DataFrameMetaHelper | undefined;
  private _plot: DataFramePlotHelper | undefined;
  private _dialogs: DataFrameDialogHelper | undefined;

  constructor(dart: any) {
    this.dart = dart;
    this.columns = toJs(api.grok_DataFrame_Columns(this.dart));
    this.rows = new RowList(this, api.grok_DataFrame_Rows(this.dart));
    this.filter = new BitSet(api.grok_DataFrame_Get_Filter(this.dart));
    this.temp = new MapProxy(api.grok_DataFrame_Get_Temp(this.dart), 'temp');
    this.tags = new MapProxy(api.grok_DataFrame_Get_Tags(this.dart), 'tags', 'string');
  }

  /** Creates a {@link DataFrame} with the specified number of rows and no columns. */
  static create(rowCount: number = 0, name?: string): DataFrame {
    const df = new DataFrame(api.grok_DataFrame(rowCount));
    if (name)
      df.name = name;
    return df;
  }

  /** Deserializes a table from its binary (d42) form; see {@link toByteArray}. */
  static fromByteArray(byteArray: Uint8Array): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromByteArray(byteArray));
  }

  /** Creates a {@link DataFrame} from the specified columns. All columns should be of the same length. */
  static fromColumns(columns: Column[]): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromColumns(columns.map((c) => c.dart)));
  }

  /** Creates a {@link DataFrame} from the specified properties with the specified row count. */
  static fromProperties(properties: Property[], rows: number = 0): DataFrame {
    let df = DataFrame.create(rows);
    for (let p of properties)
      df.columns.addNew(p.name, <ColumnType>p.propertyType);
    return df;
  }

  /** Creates a [DataFrame] from a list of objects by using object keys as column names,
   * and object values as values.
   *
   * NOTE: The implementation converts the values to strings first and then parses it,
   * so do not use this method in performance-critical paths (for instance when the
   * number of objects could be big), consider using {@link fromColumns} instead.
   *
   * @param list - List of objects.
   * {@link https://public.datagrok.ai/js/samples/data-frame/construction/create-from-objects} */
  static fromObjects(list: object[]): DataFrame | undefined {
    let table = DataFrame.create(list.length);
    if (list.length === 0)
      return;

    let names = Object.keys(list[0]);
    for (let name of names) {
      let strings = list.map((x: any) => {
        let value = x[<any>name];
        return value === null ? '':  `${value}`;
      });
      table.columns.add(Column.fromStrings(name, strings));
    }

    return table;
  }

  /** Constructs {@link DataFrame} from a comma-separated values string
   * @param csv - The content of the comma-separated values file.
   * {@link https://public.datagrok.ai/js/samples/data-frame/construction/create-from-csv} */
  static fromCsv(csv: string, options?: CsvImportOptions): DataFrame {
    return grok.data.parseCsv(csv, options);
  }

  /** Constructs {@link DataFrame} from the specified JSON string.
   * @param json - JSON document.
   * {@link https://public.datagrok.ai/js/samples/data-frame/construction/create-from-json} */
  static fromJson(json: string): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromJson(json));
  }

  /** A helper to conveniently access certain metadata properties stored in {@link tags} */
  get meta(): DataFrameMetaHelper {
    if (this._meta == undefined)
      this._meta = new DataFrameMetaHelper(this);
    return this._meta;
  }

  /** A helper for creating plots for this dataframe */
  get plot(): DataFramePlotHelper {
    if (this._plot == undefined)
      this._plot = new DataFramePlotHelper(this);
    return this._plot;
  }

  /** Standard dialogs for this table. */
  get dialogs(): DataFrameDialogHelper {
    if (this._dialogs == undefined)
      this._dialogs = new DataFrameDialogHelper(this);
    return this._dialogs;
  }

  /** Returns number of rows in the table. */
  get rowCount(): number {
    return api.grok_DataFrame_RowCount(this.dart);
  }

  /** Returns a {@link BitSet} with selected rows. */
  get selection(): BitSet {
    return new BitSet(api.grok_DataFrame_Get_Selection(this.dart));
  }

  /** Name of the dataframe.  */
  get name(): string { return api.grok_DataFrame_Get_Name(this.dart); }
  set name(s: string) { api.grok_DataFrame_Set_Name(this.dart, s); }

  /** Returns the value of the specified tag, or null if it does not exist. */
  getTag(tag: string): string | null {
    return api.grok_DataFrame_Get_Tag(this.dart, tag);
  }

  /** Sets a tag to the specified value.
   * @param tag - Key.
   * @param value - Value. */
  setTag(tag: string, value: string): DataFrame {
    if (typeof(value) !== 'string')
      throw new Error(`Tags must be strings, passed '${typeof(value)}'`);

    api.grok_DataFrame_Set_Tag(this.dart, tag, value);
    return this;
  }

  /** Returns i-th row.
   * _NOTE_: Do not use in performance-critical paths, consider accessing values via the {@link Column} instance. */
  row(rowIndex: number): Row {
    return new Row(this, rowIndex);
  }

  /** Returns idx-th value of the specified columns.
   * @param name - Column name.
   * @param idx - Row index. */
  get(name: string, idx: number): any {
    return this.getCol(name).get(idx);
  }

  /** Sets idx-th value of the specified columns.
   * @param name - Column name.
   * @param idx - Row index.
   * @param value - Value.  */
  set(name: string, idx: number, value: any): void {
    this.getCol(name).set(idx, value);
  }

  /** Returns a {@link Column} with the specified name.
   * @param nameOrIndex - Column name. */
  col(nameOrIndex: string | number): Column | null {
    if (nameOrIndex == null)
      return null;
    if (typeof nameOrIndex === 'string')
      return toJs(api.grok_DataFrame_ColumnByName(this.dart, nameOrIndex));
    else
      return this.columns.byIndex(nameOrIndex);
  }

  /** Returns a {@link Cell} with the specified row and column.
   * @param idx - Row index.
   * @param name - Column name. */
  cell(idx: number, name: string): Cell {
    return new Cell(api.grok_DataFrame_Cell(this.dart, idx, name));
  }


  /** Same as {@link col}, but throws Error if column is not found
   * @param name - Column name. */
  getCol(name: string): Column {
    let c = this.col(name);
    if (c === null)
      throw new Error(`No such column: ${name}`);
    return c;
  }

  /** Exports the content to comma-separated-values format.
   * @param options - options for the export
   * @param grid - if specified, takes visible columns, column and row order from the grid. */
  toCsv(options?: CsvExportOptions, grid?: Grid): string {
    return api.grok_DataFrame_ToCsv(this.dart, options, grid?.dart);
  }

  /** Exports the content to comma-separated-values format asynchronously
   * with converting the molblock columns to smiles if specified.
   * @param options - options for the export
   * @param grid - if specified, takes visible columns, column and row order from the grid. */
  async toCsvEx(options?: CsvExportOptions, grid?: Grid): Promise<string> {
    return api.grok_DataFrame_ToCsvEx(this.dart, options, grid?.dart);
  }

  /** Same as {@link toCsvEx}: the asynchronous export, which also converts molblock columns to SMILES when requested. */
  toCsvAsync(options?: CsvExportOptions, grid?: Grid): Promise<string> {
    return this.toCsvEx(options, grid);
  }

  /** Converts the contents to array of objects, with column names as keys.
   * Keep in mind that the internal DataFrame format is far more efficient than JSON, so
   * use it only as a convenience for working with relatively small datasets. */
  toJson(): any[] {
    const rows = this.rowCount;
    const result: any[] = Array.from({ length: rows }, () => ({}));
    for (const col of this.columns)
      for (let i = 0; i < rows; i++)
        if (!col.isNone(i))
          result[i][col.name] = col.get(i);

    return result;
  }

  /** Exports dataframe to binary */
  toByteArray(): Uint8Array {
    return api.grok_DataFrame_ToByteArray(this.dart);
  }

  /** Exports dataframe to Parquet format.
   *  @param compress - If true, applies GZIP compression. */
  toParquet(compress: boolean = false): Uint8Array {
    return api.grok_DataFrame_ToParquet(this.dart, compress);
  }

  /** Exports dataframe to Arrow IPC format. */
  toArrow(): Uint8Array {
    return api.grok_DataFrame_ToArrow(this.dart);
  }

  /** Creates a new dataframe from the specified rows and columns (all, when omitted).
   * @param rows - Rows to include.
   * @param columns - Names of the columns to include.
   * @param saveSelection - Whether the selection is copied.
   * @param saveTags - Whether the tags are copied (default). */
  clone(options?: CloneOptions): DataFrame;
  clone(rowMask?: BitSet | null, columnIds?: string[] | null, saveSelection?: boolean, saveTags?: boolean): DataFrame;
  clone(rowMask: BitSet | CloneOptions | null = null, columnIds: string[] | null = null, saveSelection: boolean = false, saveTags: boolean = true): DataFrame {
    const o: CloneOptions = rowMask === null || rowMask instanceof BitSet ? {rows: rowMask, columns: columnIds, saveSelection, saveTags} : rowMask;
    return new DataFrame(api.grok_DataFrame_Clone(this.dart, toDart(o.rows ?? null), o.columns ?? null, o.saveSelection ?? false, o.saveTags ?? true));
  }

  /** Current row.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get currentRow(): Row { return new Row(this, api.grok_DataFrame_Get_CurrentRowIdx(this.dart)); }
  set currentRow(row: Row) { api.grok_DataFrame_Set_CurrentRowIdx(this.dart, row.idx); }

  /** Index of the current row. */
  get currentRowIdx(): number { return api.grok_DataFrame_Get_CurrentRowIdx(this.dart); }
  set currentRowIdx(idx) { api.grok_DataFrame_Set_CurrentRowIdx(this.dart, idx); }

  /** Index of the mouse-over row. */
  get mouseOverRowIdx(): number { return api.grok_DataFrame_Get_MouseOverRowIdx(this.dart); }
  set mouseOverRowIdx(idx) { api.grok_DataFrame_Set_MouseOverRowIdx(this.dart, idx); }

  /** Current column.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get currentCol(): Column {
    return toJs(api.grok_DataFrame_Get_CurrentCol(this.dart));
  }

  /** Current column.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  set currentCol(col: Column) {
    api.grok_DataFrame_Set_CurrentCol(this.dart, col.dart);
  }

  /** Current cell.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get currentCell(): Cell {
    return new Cell(api.grok_DataFrame_Get_CurrentCell(this.dart));
  }

  /** Current cell.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  set currentCell(cell: Cell) {
    api.grok_DataFrame_Set_CurrentCell(this.dart, cell.dart);
  }

  /** Converts a column with the specified name to [newType],
   * removes the original column from its dataframe and adds the new column to it.
   * @deprecated Use {@link changeColumnsType} instead.
   * @param newType - @see {@link COLUMN_TYPE}
   * @param format - number format */
  changeColumnType(column: string | Column, newType: ColumnType, format: string | null = null): Column {
    return this.changeColumnsType([column], newType, format)[0];
  }

  /** Converts columns with the specified names to [newType],
   * removes the original columns from the dataframe and adds the new columns to it.
   * @param newType - @see {@link COLUMN_TYPE}
   * @param format - number format */
  changeColumnsType(columns: (string | Column)[], newType: ColumnType, format: string | null = null): Column[] {
    return toJs(api.grok_DataFrame_ChangeColumnsType(this.dart, columns.map(toDart), newType, format));
  }

  /**
   * Returns [Int32Array] that contains sorted order, or null for unsorted (original) order.
   * See also Column.getSortedOrder.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/sorting/sorted-order}
   * @param sortByColumnIds - Collection of [Column]s to use as keys for sorting.
   * @param sortOrders - List of sort orders for [sortByCols]. True == ascending.
   * @param rowMask - Mask of the rows to sort. Result array will contain [rowIndexes.length] elements. */
  getSortedOrder(sortByColumnIds: ColumnId[], sortOrders: boolean[] | null = null, rowMask: BitSet | null = null): Int32Array {
    return api.grok_DataFrame_GetSortedOrder(this.dart, sortByColumnIds.map(toDart), sortOrders, toDart(rowMask));
  }

  /** Begins building a query, using the specified columns as keys.
   * @param columnNames - Names of the columns to be used as keys. */
  groupBy(columnNames: string[] = []): GroupByBuilder {
    return new GroupByBuilder(api.grok_DataFrame_GroupBy(this.dart, columnNames));
  }

  /**
   * Unpivots the table (converts from 'wide' representation with many columns to 'tall and skinny').
   * @param copyColumnNames - columns to copy
   * @param mergeColumnNames - columns to merge. Column name will become a value in the [categoryColumnName] column,
   *                                      and column value will become a value in the [valueColumnName] column. */
  unpivot(copyColumnNames: string[], mergeColumnNames: string[], categoryColumnName: string = 'Category', valueColumnName: string = 'Value'): DataFrame {
    return new DataFrame(api.grok_DataFrame_Unpivot(this.dart, copyColumnNames, mergeColumnNames, categoryColumnName, valueColumnName));
  }

  /**
   * Merges two tables by the specified key columns.
   * @param t2 - a table to join
   * @param keyColumns1 - key column names from the first table
   * @param keyColumns2 - key column names from the second table
   * @param valueColumns1 - column names to copy from the first table.
   * Pass null to add all columns, an empty array [] to not add any columns, or an array with column names to add them specifically.
   * @param valueColumns2 - column names to copy from the second table
   * @param joinType - inner, outer, left, or right. See [DG.JOIN_TYPE]
   * @param inPlace - merges content in-place into the source table
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/join-link/join-tables} */
  join(t2: DataFrame, options: JoinOptions): DataFrame;
  join(t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1?: string[] | null, valueColumns2?: string[] | null, joinType?: JoinType, inPlace?: boolean): DataFrame;
  join(t2: DataFrame, keyColumns1: string[] | JoinOptions, keyColumns2: string[] = [], valueColumns1: string[] | null = null, valueColumns2: string[] | null = null, joinType: JoinType = JOIN_TYPE.INNER, inPlace: boolean = false): DataFrame {
    const o: JoinOptions = Array.isArray(keyColumns1) ? {keys: keyColumns1, keys2: keyColumns2, columns: valueColumns1, columns2: valueColumns2, type: joinType, inPlace} : keyColumns1;
    return new DataFrame(api.grok_JoinTables(this.dart, t2.dart, o.keys, o.keys2 ?? o.keys, o.columns ?? null, o.columns2 ?? null, o.type ?? JOIN_TYPE.INNER, o.inPlace ?? false));
  }

  /** Clears all active filters and unsets the filter bitset. */
  resetFilter(): void {
    api.grok_DataFrame_ResetFilter(this.dart);
  }

  /**
   * Appends two tables ('union' in SQL).
   * @param inPlace - whether to create a new table, or modify 'this' one. */
  append(t2: DataFrame, inPlace: boolean = false, columnsToAppend: string[] | null = null): DataFrame {
    return new DataFrame(api.grok_DataFrame_Append(this.dart, t2.dart, inPlace, columnsToAppend));
  }

  /** Appends the rows of [t] in place, adding its missing columns. */
  appendMerge(t: DataFrame): void {
    api.grok_DataFrame_Append_Merge(this.dart, t.dart);
  }

  _event(event: string): Observable<any> {
    return __obs(event, this.dart);
  }

  /** Observes the table event with the specified id. */
  onEvent(event: string): Observable<any> {
    return __obs(event, this.dart);
  }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onValuesChanged(): Observable<EventData<CellRangeArgs>> { return this._event('ddt-values-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentRowChanged(): Observable<EventData<RowChangeArgs>> { return this._event('ddt-current-row-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMouseOverRowChanged(): Observable<EventData<RowChangeArgs>> { return this._event('ddt-mouse-over-row-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentColChanged(): Observable<EventData<ColumnChangeArgs>> { return this._event('ddt-current-col-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMouseOverColChanged(): Observable<EventData<ColumnChangeArgs>> { return this._event('ddt-mouse-over-col-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentCellChanged(): Observable<EventData<CellChangeArgs>> { return this._event('ddt-current-cell-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMouseOverRowGroupChanged(): Observable<any> { return this._event('ddt-mouse-over-row-group-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onNameChanged(): Observable<any> { return this._event('ddt-table-name-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMetadataChanged(): Observable<EventData<MapChangeArgs<string, string>>> { return this._event('ddt-table-metadata-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onColumnNameChanged(): Observable<any> { return this._event('ddt-column-name-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onColumnSelectionChanged(): Observable<any> { return this._event('ddt-column-selection-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onColumnsChanged(): Observable<any> { return this._event('ddt-columns-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onColumnsAdded(): Observable<any> { return this._event('ddt-columns-added'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onColumnsRemoved(): Observable<any> { return this._event('ddt-columns-removed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onRowsAdded(): Observable<any> { return this._event('ddt-rows-added'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onRowsRemoved(): Observable<any> { return this._event('ddt-rows-removed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onRowsFiltered(): Observable<any> { return this._event('ddt-rows-filtered'); }


  get onRowsFiltering(): Observable<any> { return this._event('ddt-rows-filtering'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/advanced/semantic-type-detection} */
  get onSemanticTypeDetecting(): Observable<any> { return this._event('ddt-semantic-type-detecting'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/advanced/semantic-type-detection} */
  get onSemanticTypeDetected(): Observable<any> { return this._event('ddt-semantic-type-detected'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onDataChanged(): Observable<any> {
    return rxjs.merge(this.onValuesChanged, this.onColumnsAdded,
      this.onColumnsRemoved, this.onRowsAdded, this.onRowsRemoved);
  }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onSelectionChanged(): Observable<any> { return this.selection.onChanged; }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onFilterChanged(): Observable<any> { return this.filter.onChanged; }

  /** Fires {@link onValuesChanged} after values were set without notification. */
  fireValuesChanged(): void {
    api.grok_DataFrame_FireValuesChanged(this.dart);
  }


  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }

  /** Id of the dataframe. */
  get id(): string {
    return this.tags[TAGS.ID];
  }

  set id(id: string) {
    this.tags[TAGS.ID] = id;
  }

  /** Point density of the two columns on an [xBins] × [yBins] grid. */
  getDensity(xBins: number, yBins: number, xColName: string, yColName: string): Int32Array {
    return api.grok_MathActions_GetDensity(this.dart, xBins, yBins, xColName, yColName);
  }

  /** Server-side metadata of this table. */
  getTableInfo(): TableInfo {
    return toJs(api.grok_DataFrame_Get_TableInfo(this.dart));
  }

  _exportReopen(): DataFrame {
    return toJs(api.grok_DataFrame_ExportAndReopen(this.dart));
  }
}


export class DataFrameMetaHelper {
  df: DataFrame;

  readonly formulaLines: DataFrameFormulaLinesHelper;
  readonly annotationRegions: DataFrameAnnotationRegionsHelper;

  async detectSemanticTypes() {
    await grok.data.detectSemanticTypes(this.df);
  }

  constructor(df: DataFrame) {
    this.df = df;
    this.formulaLines = new DataFrameFormulaLinesHelper(this.df);
    this.annotationRegions = new DataFrameAnnotationRegionsHelper(this.df);
  }

  /** This data will be picked up by {@link Grid} to construct groups. */
  setGroups(groups: GroupsDescription | null): void {
    if (groups)
      this.df.tags['.columnGroups'] = JSON.stringify(groups);
    else
      delete this.df.tags['.columnGroups'];
  }
}


export class DataFramePlotHelper {
  private readonly df: DataFrame;
  constructor(df: DataFrame) {
    this.df = df;
  }

  fromType(viewerType: ViewerType, options: object | null = null): Promise<Widget> {
    return toJs(api.grok_Viewer_FromType_Async(viewerType, this.df.dart, _toJson(options)));
  }

  scatter(options?: null | ViewerOptions<IScatterPlotSettings>): ScatterPlotViewer { return DG.Viewer.scatterPlot(this.df, options); }
  grid(options?: null | ViewerOptions<IGridSettings>): Grid { return DG.Viewer.grid(this.df, options); }
  tile(options?: null | ViewerOptions<ITileViewerSettings>): TileViewer { return DG.Viewer.tile(this.df, options); }
  form(options?: null | ViewerOptions<IFormSettings>): FormViewer { return DG.Viewer.form(this.df, options); }
  histogram(options?: null | ViewerOptions<IHistogramSettings>): HistogramViewer { return DG.Viewer.histogram(this.df, options); }
  bar(options?: null | ViewerOptions<IBarChartSettings>): BarChartViewer { return DG.Viewer.barChart(this.df, options); }
  heatMap(options?: null | ViewerOptions<IGridSettings>): Grid { return DG.Viewer.heatMap(this.df, options); }
  box(options?: null | ViewerOptions<IBoxPlotSettings>): BoxPlot { return DG.Viewer.boxPlot(this.df, options); }
  line(options?: null | ViewerOptions<ILineChartSettings>): LineChartViewer { return DG.Viewer.lineChart(this.df, options); }
  network(options?: null | ViewerOptions<INetworkDiagramSettings>): NetworkDiagramViewer { return DG.Viewer.network(this.df, options); }

}

export class DataFrameDialogHelper {
  private readonly df: DataFrame;
  constructor(df: DataFrame) {
    this.df = df;
  }

  async addNewColumn(): Promise<void> { (await grok.functions.eval('AddNewColumn')).prepare().edit(); }
}
