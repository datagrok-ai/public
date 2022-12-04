import * as rxjs from 'rxjs';
import {
  AGG,
  TYPE,
  TAGS,
  JoinType,
  SemType,
  ColumnType,
  SimilarityMetric,
  AggregationType,
  CsvImportOptions,
  IndexPredicate, FLOAT_NULL, ViewerType, ColorCodingType, MarkerCodingType, ColorType, ColumnAggregationType
} from "./const";
import {__obs, EventData, MapChangeArgs, observeStream} from "./events";
import {toDart, toJs} from "./wrappers";
import {SIMILARITY_METRIC} from "./const";
import {MapProxy, _getIterator, _toIterable, _toJson, DartList} from "./utils";
import {Observable}  from "rxjs";
import {filter} from "rxjs/operators";
import {Widget} from "./widgets";
import {Grid} from "./grid";
import {FilterState, ScatterPlotViewer, Viewer} from "./viewer";
import {Property, TableInfo} from "./entities";
import {FormulaLinesHelper} from "./helpers";
import dayjs from "dayjs";

declare let grok: any;
declare let DG: any;
let api = <any>window;
type RowPredicate = (row: Row) => boolean;
type Comparer = (a: any, b: any) => number;
type IndexSetter = (index: number, value: any) => void;
type ColumnId = number | string | Column;


/** Column CSV export options */
export interface ColumnCsvExportOptions {
  /** Custom column format to be used */
  format: string;

  /** Additional options */
  [index: string]: any;
}


/** Column name -> options */
export interface ColumnsCsvExportOptions {
  [index: string]: ColumnCsvExportOptions;
}


/** Csv export options to be used in {@link DataFrame.toCsv} */
export interface CsvExportOptions {

  /** Field delimiter; comma if not specified */
  delimiter?: string;

  /** New line character; \n if not specified */
  newLine?: string;

  /** Textual representation of the missing value. Empty string if not specified. */
  missingValue?: string;

  /** Whether the header row containing column names is included. True if not specified. */
  includeHeader?: boolean;

  /** Whether only selected columns are included. False if not specified. */
  selectedColumnsOnly?: boolean;

  /** Whether only filtered rows are included. Will be combined with [selectedRowsOnly]. */
  filteredRowsOnly?: boolean;

  /** Whether only selected rows are included. Will be combined with [filteredRowsOnly]. */
  selectedRowsOnly?: boolean;

  /** Column order */
  columns?: string[];

  /** Expands qualified numbers into two columns: `sign(column)` and `column` */
  qualifierAsColumn?: boolean;

  /** Column-specific formats (column name -> format).
      For format examples, see [dateTimeFormatters]. */
  columnOptions?: ColumnsCsvExportOptions;
}


/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 *
 * In the API, the terms "Table" and "DataFrame" are used interchangeably.
 *
 * Usage samples: {@link https://public.datagrok.ai/js/samples/data-frame/manipulate}
 * Usage details: {@link https://datagrok.ai/help/develop/advanced/data-frame}
 * Implementation details: {@link https://datagrok.ai/help/develop/admin/architecture#in-memory-database}
 */
export class DataFrame {
  public readonly dart: any;
  public columns: ColumnList;
  public rows: RowList;
  public filter: BitSet;
  public temp: any;
  public tags: any;
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
  static create(rowCount: number = 0): DataFrame {
    return new DataFrame(api.grok_DataFrame(rowCount));
  }

  static fromByteArray(byteArray: Uint8Array): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromByteArray(byteArray));
  }

  /** Creates a {@link DataFrame} from the specified columns. All columns should be of the same length. */
  static fromColumns(columns: Column[]): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromColumns(columns.map((c) => c.dart)));
  }

  /** Creates a {@link DataFrame} from the specified properties with the specified row count. */
  static fromProperties(properties: Property[], rows: number = 0) {
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
   * @param {object[]} list - List of objects.
   * @returns {DataFrame}
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/construction/create-from-objects}
   * */
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
   * @param {string} csv - The content of the comma-separated values file.
   * @param {CsvImportOptions} options
   * @returns {DataFrame}
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/construction/create-from-csv}
   * */
  static fromCsv(csv: string, options?: CsvImportOptions): DataFrame {
    return grok.data.parseCsv(csv, options);
  }

  /** Constructs {@link DataFrame} from the specified JSON string.
   * @param {string} json - JSON document.
   * @returns {DataFrame}
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/construction/create-from-json}
   * */
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

  /** Returns the value of the specified tag, or null if it does not exist.
   * @returns {string} */
  getTag(tag: string): string | null {
    return api.grok_DataFrame_Get_Tag(this.dart, tag);
  }

  /** Sets a tag to the specified value.
   * @param {string} tag - Key.
   * @param {string} value - Value. */
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
   * @param {string} name - Column name.
   * @param {number} idx - Row index. */
  get(name: string, idx: number): any {
    return this.getCol(name).get(idx);
  }

  /** Sets idx-th value of the specified columns.
   * @param {string} name - Column name.
   * @param {number} idx - Row index.
   * @param value - Value. */
  set(name: string, idx: number, value: any): void {
    this.getCol(name).set(idx, value);
  }

  /** Returns a {@link Column} with the specified name.
   * @param {string} name - Column name.
   * @returns {Column} */
  col(name: string): Column | null {
    return toJs(api.grok_DataFrame_ColumnByName(this.dart, name));
  }

  /** Returns a {@link Cell} with the specified row and column.
   * @param {number} idx - Row index.
   * @param {string} name - Column name.
   * @returns {Cell} */
  cell(idx: number, name: string): Cell {
    return new Cell(api.grok_DataFrame_Cell(this.dart, idx, name));
  }

  /** Same as {@link col}, but throws Error if column is not found
   * @param {string} name - Column name.
   * @returns {Column} */
  getCol(name: string): Column {
    let c = this.col(name);
    if (c === null)
      throw new Error(`No such column: ${name}`);
    return c;
  }

  /** Exports the content to comma-separated-values format. */
  toCsv(options?: CsvExportOptions): string {
    return api.grok_DataFrame_ToCsv(this.dart, options);
  }

  /** Exports dataframe to binary */
  toByteArray(): Uint8Array {
    return api.grok_DataFrame_ToByteArray(this.dart);
  }

  /** Creates a new dataframe from the specified row mask and a list of columns.
   * @param {BitSet} rowMask - Rows to include.
   * @param {string[]} columnIds - Columns to include.
   * @param {boolean} saveSelection - Whether selection should be saved. */
  clone(rowMask: BitSet | null = null, columnIds: string[] | null = null, saveSelection: boolean = false): DataFrame {
    return new DataFrame(api.grok_DataFrame_Clone(this.dart, toDart(rowMask), columnIds, saveSelection));
  }

  /** Current row.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get currentRow(): Row {
    return new Row(this, api.grok_DataFrame_Get_CurrentRowIdx(this.dart));
  }

  /** Current row.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  set currentRow(idx) {
    api.grok_DataFrame_Set_CurrentRowIdx(this.dart, idx);
  }

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
   * @param {string|Column} column
   * @param {string} newType
   * @param {string=} format
   * @returns {Column} */
  changeColumnType(column: string | Column, newType: ColumnType, format: string | null = null): Column {
    return toJs(api.grok_DataFrame_ChangeColumnType(this.dart, toDart(column), newType, format));
  }

  /**
   * Returns [Int32Array] that contains sorted order, or null for unsorted (original) order.
   * See also Column.getSortedOrder.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/sorting/sorted-order}
   * @param {Object[]} sortByColumnIds - Collection of [Column]s to use as keys for sorting.
   * @param {boolean[]} sortOrders - List of sort orders for [sortByCols]. True == ascending.
   * @param {BitSet} rowMask - Mask of the rows to sort. Result array will contain [rowIndexes.length] elements.
   * @returns {Int32Array}
   * */
  getSortedOrder(sortByColumnIds: ColumnId[], sortOrders: boolean[] | null = null, rowMask: BitSet | null = null): Int32Array {
    return api.grok_DataFrame_GetSortedOrder(this.dart, sortByColumnIds.map(toDart), sortOrders, toDart(rowMask));
  }

  /** Begins building a query, using the specified columns as keys.
   * @param {string[]} columnNames - Names of the columns to be used as keys.
   * @returns {GroupByBuilder}
   *  */
  groupBy(columnNames: string[] = []): GroupByBuilder {
    return new GroupByBuilder(api.grok_DataFrame_GroupBy(this.dart, columnNames));
  }

  /**
   * Unpivots the table (converts from 'wide' representation with many columns to 'tall and skinny').
   * @param {String[]} copyColumnNames - columns to copy
   * @param {String[]} mergeColumnNames - columns to merge. Column name will become a value in the [categoryColumnName] column,
   *                                      and column value will become a value in the [valueColumnName] column.
   * @param {String} categoryColumnName
   * @param {String} valueColumnName
   * */
  unpivot(copyColumnNames: string[], mergeColumnNames: string[], categoryColumnName: string = 'Category', valueColumnName: string = 'Value'): DataFrame {
    return new DataFrame(api.grok_DataFrame_Unpivot(this.dart, copyColumnNames, mergeColumnNames, categoryColumnName, valueColumnName));
  }

  /**
   * Merges two tables by the specified key columns.
   * @param {DataFrame} t2
   * @param {string[]} keyColumns1
   * @param {string[]} keyColumns2
   * @param {string[]} valueColumns1
   * @param {string[]} valueColumns2
   * @param {JoinType} joinType
   * @param {boolean} inPlace - merges content in-place into the source table
   * @returns {DataFrame}
   * */
  join(t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[], joinType: JoinType, inPlace: boolean): DataFrame {
    return new DataFrame(api.grok_JoinTables(this.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }

  /**
   * Appends two tables ('union' in SQL).
   * @param {DataFrame} t2
   * @param {boolean} inPlace - whether to create a new table, or modify 'this' one.
   * @param {String[]} columnsToAppend
   * @returns {DataFrame}
   * */
  append(t2: DataFrame, inPlace: boolean = false, columnsToAppend: string[] | null = null): DataFrame {
    return new DataFrame(api.grok_DataFrame_Append(this.dart, t2.dart, inPlace, columnsToAppend));
  }

  _event(event: string): Observable<any> {
    return __obs(event, this.dart);
  }

  onEvent(event: string): Observable<any> {
    return __obs(event, this.dart);
  }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onValuesChanged(): Observable<any> { return this._event('ddt-values-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentRowChanged(): Observable<any> { return this._event('ddt-current-row-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMouseOverRowChanged(): Observable<any> { return this._event('ddt-mouse-over-row-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentColChanged(): Observable<any> { return this._event('ddt-current-col-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onMouseOverColChanged(): Observable<any> { return this._event('ddt-mouse-over-col-changed'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/current-elements} */
  get onCurrentCellChanged(): Observable<any> { return this._event('ddt-current-cell-changed'); }

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

  /** @returns {Observable} */
  get onRowsFiltering(): Observable<any> { return this._event('ddt-rows-filtering'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/advanced/semantic-type-detection} */
  get onSemanticTypeDetecting(): Observable<any> { return this._event('ddt-semantic-type-detecting'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/advanced/semantic-type-detection} */
  get onSemanticTypeDetected(): Observable<any> { return this._event('ddt-semantic-type-detected'); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onDataChanged(): Observable<any> {
    return rxjs.concat(this.onValuesChanged, this.onColumnsAdded,
      this.onColumnsRemoved, this.onRowsAdded, this.onRowsRemoved);
  }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onSelectionChanged(): Observable<any> { return this.selection.onChanged; }

  /** Sample: {@link https://public.datagrok.ai/js/samples/data-frame/events/events} */
  get onFilterChanged(): Observable<any> { return this.filter.onChanged; }

  fireValuesChanged(): void {
    api.grok_DataFrame_FireValuesChanged(this.dart);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }

  /** Id of the dataframe.
   * @returns {string}*/
  get id(): string {
    return this.tags[TAGS.ID];
  }

  set id(id: string) {
    this.tags[TAGS.ID] = id;
  }

  getDensity(xBins: number, yBins: number, xColName: string, yColName: string): Int32Array {
    return api.grok_MathActions_GetDensity(this.dart, xBins, yBins, xColName, yColName);
  }

  getTableInfo(): TableInfo {
    return toJs(api.grok_DataFrame_Get_TableInfo(this.dart));
  }
}

/** Represents a row. Allows for quick property access like "row.height". */
export class Row {
  table: DataFrame;
  readonly idx: number;

  /**
   * Creates a {@link Row} from the specified {@link DataFrame} and row index.
   * @param {DataFrame} table
   * @param {number} idx */
  constructor(table: DataFrame, idx: number) {

    /** @member {DataFrame} */
    this.table = table;
    /** @member {number} */
    this.idx = idx;

    return new Proxy(this, {
      set(target, name: string, value) {
        if (target.hasOwnProperty(name)) {
          Object.entries(target)[<any>name] = value;
          return true;
        }
        target.table.set(name, target.idx, value);
        return true;
      },
      get(target: any, name) {
        if (name == 'cells' || name == 'get' || target.hasOwnProperty(name))
          return target[<any>name];
        return target.table.get(name, target.idx);
      }
    });
  }

  get cells(): Iterable<Cell> { return _toIterable(api.grok_Row_Get_Cells(this.table.dart, this.idx)); }

  /** Returns this row's value for the specified column
   * @param {string} columnName
   * @returns {Object} */
  [columnName: string]: any

  /** Returns this row's value for the specified column
   * @param {string} columnName
   * @returns {Object} */
  get(columnName: string): any {
    return this.table.getCol(columnName).get(this.idx);
  }

  toDart(): any { return api.grok_Row(this.table.dart, this.idx); }
}

type KnownColumnTags = 'format' | 'colors';

/** Strongly-typed column.
 * Use {@link get} and {@link set} to access elements by index.
 * */
export class Column<T = any> {
  public dart: any;
  public temp: any;
  public tags: any;
  private _dialogs: ColumnDialogHelper | undefined;
  private _meta: ColumnMetaHelper | undefined;

  constructor(dart: any) {
    this.dart = dart;
    this.temp = new MapProxy(api.grok_Column_Get_Temp(this.dart));
    this.tags = new MapProxy(api.grok_Column_Get_Tags(this.dart));
  }

  get meta(): ColumnMetaHelper {
    if (this._meta == undefined)
      this._meta = new ColumnMetaHelper(this);
    return this._meta;
  }

  /** Creates a {@link Column} from the list of string values
   * Please note that method performs type promotion if all listed values are numeric
   *
   * @param {string} name - Column name
   * @param {Array} list - List of column values
   *
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/construction/create-from-columns}
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/construction/create-from-arrays}
   *
   */
  static fromStrings(name: string, list: string[]): Column {
    return toJs(api.grok_Column_FromStrings(name, list));
  }

  /** Creates a {@link Column} with explicitly specified type
   *
   * @param type - column type code {@link COLUMN_TYPE}
   * @param {string} name - Column name
   * @param length - Column length (should match row count of the data frame )
   *
   * {@link DataFrame.create}
   * {@see COLUMN_TYPE}
   */
  static fromType(type: ColumnType, name?: string | null, length: number = 0): Column {
    return toJs(api.grok_Column_FromType(type, name, length));
  }

  /** [array] will be not be copied and will be used as column's storage */
  static fromInt32Array(name: string, array: Int32Array, length: number | null = null): Column<number> {
    return toJs(api.grok_Column_FromInt32Array(name, array, length));
  }

  /** [array] will be not be copied and will be used as column's storage */
  static fromFloat32Array(name: string, array: Float32Array, length: number | null = null): Column<number> {
    return toJs(api.grok_Column_FromFloat32Array(name, array, length));
  }

  /**
   * Creates a {@link Column} from the list of values.
   * @param {string} type
   * @param {string} name
   * @param {object[]} list
   * @returns {Column} */
  static fromList(type: ColumnType, name: string, list: any[]): Column {
    return toJs(api.grok_Column_FromList(type, name, list));
  }

  /** Creates a {Column} from the bitset.
   * @param {string} name
   * @param {BitSet} bitset
   * @returns {Column} */
  static fromBitSet(name: string, bitset: BitSet): Column<boolean> {
    return toJs(api.grok_Column_FromBitSet(name, bitset.dart));
  }

  /** Creates an integer column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static int(name: string, length: number = 0): Column<number> {
    return Column.fromType(TYPE.INT, name, length);
  }

  /** Creates a floating point column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static float(name: string, length: number = 0): Column<number> {
    return Column.fromType(TYPE.FLOAT, name, length);
  }

  /** Creates a string column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static string(name: string, length: number = 0): Column<string> {
    return Column.fromType(TYPE.STRING, name, length);
  }

  /** Creates a boolean column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static bool(name: string, length: number = 0): Column<boolean> {
    return Column.fromType(TYPE.BOOL, name, length);
  }

  /** Creates a datetime column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static dateTime(name: string, length: number = 0): Column {
    return Column.fromType(TYPE.DATE_TIME, name, length);
  }

  /** Creates a column containing dataframes with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static dataFrame(name: string, length: number = 0): Column<DataFrame> {
    return Column.fromType(TYPE.DATA_FRAME, name, length);
  }

  /** Creates a qualified number column with the specified name and length.
   *  Initialized values with [values], if it is specified; strips out the qualifier
   *  part if [exact] is true.
   *
   * @param {string} name
   * @param {number} length
   * @param {number[]} values
   * @param {boolean} exact - if true, strips out qualifier from [values].
   * */
  static qnum(name: string, length: number = 0, values: number[] = [], exact: boolean = true): Column<number> {
    let col = Column.fromType(TYPE.QNUM, name, length);
    if (values !== null) {
      let buffer = col.getRawData();
      for (let i = 0; i < length; i++) {
        let val = values[i] === undefined || values[i] === null ? FLOAT_NULL : values[i];
        buffer[i] = exact ? Qnum.exact(val) : val;
      }
      col.setRawData(buffer);
    }
    return col;
  }

  /** Column data type.
   * @type {string} */
  get type(): ColumnType {
    return api.grok_Column_Get_Type(this.dart);
  }

  /** Is this column virtual
   * @type {boolean} */
  get isVirtual(): boolean {
    return api.grok_Column_IsVirtual(this.dart);
  }

  /** Number of elements
   * @type {number} */
  get length(): number {
    return api.grok_Column_Get_Length(this.dart);
  }

  /** Parent table
   * @type {DataFrame} */
  get dataFrame(): DataFrame {
    return toJs(api.grok_Column_Get_DataFrame(this.dart));
  }

  /** Semantic type
   * @type {string} */
  get semType(): SemType {
    return api.grok_Column_Get_SemType(this.dart);
  }

  set semType(s: SemType) {
    api.grok_Column_Set_SemType(this.dart, s);
  }

  /** Layout column ID
   @type {string} */
  get layoutColumnId(): string {
    return api.grok_Column_Get_LayoutColumnId(this.dart);
  }

  set layoutColumnId(s: string) {
    api.grok_Column_Set_LayoutColumnId(this.dart, s);
  }

  /** @type {string} */
  get name(): string {
    return api.grok_Column_Get_Name(this.dart);
  }

  set name(s: string) {
    api.grok_Column_Set_Name(this.dart, s);
  }

  /** Version of the column. Increases each time the column was changed
   * @returns {number} */
  get version(): number {
    return api.grok_Column_Get_Version(this.dart);
  }

  get dialogs(): ColumnDialogHelper {
    if (this._dialogs == undefined)
      this._dialogs = new ColumnDialogHelper(this);
    return this._dialogs;
  }

  // Obsolete. Recommended method is "meta.colors".
  get colors(): ColumnColorHelper {
    return this.meta.colors;
  }

  // Obsolete. Recommended method is "meta.markers".
  get markers(): ColumnMarkerHelper {
    return this.meta.markers;
  }

  /**
   * Initializes all values in the column to [columnInitializer].
   * @param {string | number | boolean | Function} valueInitializer value, or a function that returns value by index
   * @returns {Column}
   * */
  init(valueInitializer: string | number | boolean | ((ind: number) => any)): Column {
    let type = typeof valueInitializer;
    if (type === 'function')
      api.grok_Column_Init(this.dart, valueInitializer);
    else if (type === 'number' || type === 'string' || type === 'boolean')
      api.grok_Column_SetAllValues(this.dart, valueInitializer);
    return this;
  }

  /** FOR EXPERT USE ONLY!
   *
   * Returns the raw buffer containing data.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/performance/access}
   * Return type depends on the column type:
   * {Int32Array} for ints, {@link INT_NULL} represents null.
   * {Float32Array} for floats, {@link FLOAT_NULL} represents null.
   * {Float64Array} for qnums, {@link FLOAT_NULL} represents null.
   * {Float64Array} for datetime, in microseconds since epoch, {@link FLOAT_NULL} represents null.
   * {Int32Array} for strings indexes of {@link categories}.
   * {Uint32Array} bit array.
   * @returns {Array} */
  getRawData(): Int32Array | Float32Array | Float64Array | Uint32Array {
    // a hack that extracts the real underlying array from the Dart Column
    const handle = api.grok_Column_GetRawData(this.dart);
    const TypedArray = Object.getPrototypeOf(Uint8Array);
    for (const k of Object.keys(handle)) {
      const v = handle[k];
      if (v instanceof TypedArray)
        return v;
    }
    return api.grok_Column_GetRawDataDartium(this.dart);
  }

  /** Linearly maps idx-th value to the [0,1] interval, where 0 represents column minimum, and 1 represents maximum. */
  scale(idx: number): number {
    return api.grok_Column_Scale(this.dart, idx);
  }

  setRawData(rawData: Int32Array | Float32Array | Float64Array | Uint32Array, notify: boolean = true): void {
    api.grok_Column_SetRawData(this.dart, rawData, notify);
  }

  /** Gets i-th value
   * @param {number} row - row index
   * @returns {object} - or null if isNone(i) */
  get(row: number): T | null {
    return api.grok_Column_GetValue(this.dart, row);
  }

  /** Returns i-th value as string, taking into account value format defined for the column.
   *  An empty string is returned if there is no value.
   * @param {number} i
   * @returns {string} */
  getString(i: number): string {
    return api.grok_Column_GetString(this.dart, i);
  }

  /** Attempts to set i-th value by converting a provided string to the corresponding strongly-typed value.
   *  Returns true if text was successfully parsed and set, otherwise false.
   *  Examples: dateColumn.setString('April 1, 2020');
   *            intColumn.setString('42');
   * @param {number} i
   * @param {string} str
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired
   * @returns {boolean} */
  setString(i: number, str: string, notify: boolean = true): boolean {
    return api.grok_Column_SetString(this.dart, i, str, notify);
  }

  /**
   * Sets [i]-th value to [x], and optionally notifies the dataframe about this change.
   * @param {number} i
   * @param value
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired
   */
  set(i: number, value: T | null, notify: boolean = true): void {
    api.grok_Column_SetValue(this.dart, i, toDart(value), notify);
  }

  /** Returns whether i-th value is missing.
   * @param {number} i - Row index.
   * @returns {boolean} */
  isNone(i: number): boolean {
    return api.grok_Column_IsNone(this.dart, i);
  }

  /** Gets the value of the specified tag.
   *  @param {string} tag
   *  @returns {string} */
  getTag(tag: string): string {
    return api.grok_Column_Get_Tag(this.dart, tag);
  }

  /** Sets a tag to the specified value.
   * @param {string} tag - Key.
   * @param {string} value - Value.
   * @returns {Column}. * */
  setTag(tag: string, value: string): Column {
    api.grok_Column_Set_Tag(this.dart, tag, value);
    return this;
  }

  /** Compacts the internal column representation.
   *  Currently, it only affects string columns where values were modified. */
  compact() {
    return api.grok_Column_Compact(this.dart);
  }

  /** Copies column values to an array.
   * Avoid using this method for performance-critical routines; consider using {@link getRawData} */
  toList(): Array<any> {
    return api.grok_Column_ToList(this.dart);
  }

  /** Returns all unique strings in a sorted order. Applicable to string column only. */
  get categories(): string[] {
    return api.grok_Column_Categories(this.dart);
  }

  /** Returns i-th category. Applicable to string column only. */
  getCategory(categoryIndex: number): string {
    return api.grok_Column_GetCategory(this.dart, categoryIndex);
  }

  /** Sets order of categories */
  setCategoryOrder(order: string[]) {
    api.grok_Column_SetCategoryOrder(this.dart, order);
  }

  /** Gets order of categories */
  getCategoryOrder(): string[] { return api.grok_Column_GetCategoryOrder(this.dart); }

  /** Returns an array of indexes sorted using [valueComparer]. */
  getSortedOrder(): Int32Array { return api.grok_Column_GetSortedOrder(this.dart); }

  /** Value comparison function to be used for sorting. Null means default sorting.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/sorting/custom-comparer}  */
  get valueComparer(): Comparer | null { return api.grok_Column_Get_ValueComparer(this.dart); }
  set valueComparer( cmp: Comparer | null) { api.grok_Column_Set_ValueComparer(this.dart, cmp); }

  /** Column's minimum value. The result is cached. */
  get min(): number { return api.grok_Column_Min(this.dart); }

  /** Column's maximum value. The result is cached. */
  get max(): number { return api.grok_Column_Max(this.dart); }

  /** Checks whether the column passes the specified [filter].
   * [filter] can be either specific data [type] such as 'int' or 'string',
   * more broadly - 'numerical', or 'categorical', or null for any columns. */
  matches(filter: ColumnType | 'numerical' | 'categorical' | null): boolean {
    return api.grok_Column_Matches(this.dart, filter);
  }

  /** Basic descriptive statistics. The result is cached. */
  get stats(): Stats { return Stats.fromColumn(this); }

  /** An iterator over all values in this column. */
  * values() {
    for (let i = 0; i < this.length; i++) {
      yield this.get(i);
    }
  }

  /** Applies the specified formula to a calculated column.
   * Returns a new column object, if applied successfully, and null otherwise.*/
  async applyFormula(formula: string, type: ColumnType | null = null, treatAsString: boolean = false): Promise<Column | null> {
    if (!(this.name && this.dataFrame?.columns.contains(this.name)))
      return null;
    if (type == null)
      type = this.type;
    return new Promise((resolve, reject) => api.grok_Column_ApplyFormula(
      this.dart, formula, type, treatAsString, (c: any) => resolve(toJs(c)), (e: any) => reject(e)));
  }

  /** Creates and returns a new column by converting [column] to the specified [newType]. */
  convertTo(newType: string, format: string | null = null): Column {
    return toJs(api.grok_Column_ConvertTo(this.dart, newType, format));
  }

  /** @returns {string} - string representation of this column */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }

  /** Aggregates a column using the [type] function, which corresponds to
   *  - [DG.AGG] for int, bigint, float, qnum columns,
   *  - [DG.STR_AGG] and [DG.STAT_COUNTS] for string columns,
   *  - [DG.STAT_COUNTS], [DG.AGG.MIN], [DG.AGG.MAX], [DG.AGG.AVG] for datetime columns,
   *  - [DG.AGG.TOTAL_COUNT] and [DG.AGG.MISSING_VALUE_COUNT] for virtual columns.
  */
  aggregate(type: ColumnAggregationType): any {
    return api.grok_Column_Aggregate(this.dart, type);
  }
}
export class BigIntColumn extends Column<BigInt> {
  /**
   * Gets [i]-th value.
   */
  get(row: number): BigInt | null {
    let v = api.grok_BigIntColumn_GetValue(this.dart, row);
    if (v == null)
      return null;
    return BigInt(v);
  }

  /**
   * Sets [i]-th value to [x], and optionally notifies the dataframe about this change.
   */
  set(i: number, value: BigInt | null, notify: boolean = true): void {
    api.grok_BigIntColumn_SetValue(this.dart, i, value?.toString(), notify);
  }
}

export class DateTimeColumn extends Column<dayjs.Dayjs> {
  /**
   * Gets [i]-th value.
   */
  get(row: number): dayjs.Dayjs | null {
    let v = api.grok_DateTimeColumn_GetValue(this.dart, row);
    if (v == null)
      return null;
    return dayjs(v);
  }

  /**
   * Sets [i]-th value to [x], and optionally notifies the dataframe about this change.
   */
  set(i: number, value: dayjs.Dayjs | null, notify: boolean = true): void {
    // @ts-ignore
    if (value == '')
      value = null;
    if (!(dayjs.isDayjs(value) || value == null))
      value = dayjs(value);
    api.grok_DateTimeColumn_SetValue(this.dart, i, value?.valueOf(), notify);
  }
}


export class ObjectColumn extends Column<any> {
  /**
   * Gets [i]-th value.
   */
  get(row: number): any | null {
    return DG.toJs(api.grok_Column_GetValue(this.dart, row));
  }
}

/** Columns in a [DataFrame]. */
export class ColumnList {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Number of elements in the column. */
  get length(): number { return api.grok_ColumnList_Length(this.dart); }

  /** Column with the corresponding name (case-insensitive). */
  byName(name: string): Column { return toJs(api.grok_ColumnList_ByName(this.dart, name)); }

  /** Maps names to columns. */
  byNames(names: string[]): Column[] { return names.map(name => this.byName(name)); }

  /** Column by index */
  byIndex(index: number): Column { return toJs(api.grok_ColumnList_ByIndex(this.dart, index)); }

  /** First column of [semType], or null. */
  bySemType(semType: SemType): Column | null {
    let col = api.grok_ColumnList_BySemType(this.dart, semType);
    return col == null ? null : col;
  }

  /** All columns of [semType], or empty list. */
  bySemTypeAll(semType: SemType): Column[] {
    return api.grok_ColumnList_BySemTypeAll(this.dart, semType);
  }

  /** Finds columns by the corresponding semTypes, or null, if any of the sem types could not be found. */
  bySemTypesExact(semTypes: SemType[]): Column[] | null {
    let columns = <any>[];
    for (let semType of semTypes) {
      let col = this.bySemType(semType);
      if (col == null)
        return null;
      columns.push(col);
    }
    return columns;
  }

  /** Finds columns by specified tags and values: {'tag': 'value'}.
   * Pass null or undefined to match any value of a tag.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/find-columns} */
  byTags(tags: object): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_ByTags(this.dart, tags));
  }

  /** Finds categorical columns.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/find-columns} */
  get categorical(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_Categorical(this.dart));
  }

  /** Finds numerical columns.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/find-columns} */
  get numerical(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_Numerical(this.dart));
  }

  /** Array containing column names. */
  names(): string[] { return api.grok_ColumnList_Names(this.dart); }

  /** Creates an array of columns. */
  toList(): Column[] {
    return this.names().map((name: string) => this.byName(name));
  }

  /** Adds a column, and optionally notifies the parent dataframe.
   * @param {Column} column
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired
   * @returns {Column} */
  add(column: Column, notify: boolean = true): Column {
    api.grok_ColumnList_Add(this.dart, column.dart, notify);
    return column;
  }

  getOrCreate(name: string, type: ColumnType, length: number): Column {
    return this.contains(name) ?
    this.byName(name) :
    this.add(DG.Column.fromType(type, name, length));
  }

  /** Inserts a column, and optionally notifies the parent dataframe.
   * @param {Column} column
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired
   * @param {int} index
   * @returns {Column} */
  insert(column: Column, index: number | null = null, notify: boolean = true): Column {
    api.grok_ColumnList_Insert(this.dart, column.dart, index, notify);
    return column;
  }

  /** Adds an empty column of the specified type.
   * @param {string} name
   * @param {ColumnType} type
   * @returns {Column} */
  addNew(name: string, type: ColumnType): Column {
    return toJs(api.grok_ColumnList_AddNew(this.dart, name, type));
  }

  /** Adds calculated column.
   * @param {string} name
   * @param {string} expression
   * @param {ColumnType} type
   * @param {bool} treatAsString - if true, [expression] is not evaluated as formula and is treated as a regular string value instead
   * @returns {Column} */
  addNewCalculated(name: string, expression: string, type: ColumnType | 'auto' = 'auto', treatAsString: boolean = false): Promise<Column> {
    return new Promise((resolve, reject) => api.grok_ColumnList_AddNewCalculated(this.dart, name, expression, type, treatAsString, (c: any) => resolve(toJs(c)), (e: any) => reject(e)));
  }

  _getNewCalculated(name: string, expression: string, type: ColumnType | 'auto' = 'auto', treatAsString: boolean = false): Promise<Column> {
    return new Promise((resolve, reject) => api.grok_ColumnList_GetNewCalculated(this.dart, name, expression, type, treatAsString, (c: any) => resolve(toJs(c)), (e: any) => reject(e)));
  }

  /** Creates and adds a string column. */
  addNewString(name: string): Column<string> { return this.addNew(name, TYPE.STRING); }

  /** Creates and adds an integer column
   *  {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns} */
  addNewInt(name: string): Column<number> { return this.addNew(name, TYPE.INT); }

  /** Creates and adds a float column */
  addNewFloat(name: string): Column<number> { return this.addNew(name, TYPE.FLOAT); }

  /** Creates and adds a qualified number column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewQnum(name: string): Column<number> { return this.addNew(name, TYPE.QNUM); }

  /** Creates and adds a datetime column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewDateTime(name: string): Column { return this.addNew(name, TYPE.DATE_TIME); }

  /** Creates and adds a boolean column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewBool(name: string): Column<boolean> { return this.addNew(name, TYPE.BOOL); }

  /** Creates and adds a boolean column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewBytes(name: string): Column<Uint8Array> { return this.addNew(name, TYPE.BYTE_ARRAY); }

  /** Creates and adds a virtual column.
   * @param {string} name
   * @param getValue - value constructor function that accepts int index and returns value
   * @param setValue - function that gets invoked when a column cell value is set
   * @param {String} type - column type
   * @returns {Column}
   *
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/advanced/virtual-int-column}
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/advanced/virtual-columns}
   * */
  addNewVirtual(
      name: string,
      getValue: (ind: number) => any, type = TYPE.OBJECT,
      setValue: IndexSetter | null = null): Column {
    return toJs(api.grok_ColumnList_AddNewVirtual(this.dart, name, getValue, setValue, type));
  }

  /** Removes column by name (case-insensitive).*/
  remove(column: string, notify: boolean = true): ColumnList {
    api.grok_ColumnList_Remove(this.dart, column, notify);
    return this;
  }

  /** Checks whether this list contains a column with the specified name. The check is case-insensitive.
   * @returns {boolean} */
  contains(columnName: string): boolean {
    return api.grok_ColumnList_Contains(this.dart, columnName);
  }

  /** Replaces the column with the new column.
   * @param {Column} columnToReplace
   * @param {Column} newColumn
   * */
  replace(columnToReplace: Column | string, newColumn: Column): Column {
    return toJs(api.grok_ColumnList_Replace(this.dart, (typeof columnToReplace === 'string') ? columnToReplace:  columnToReplace.dart, newColumn.dart));
  }

  /** Returns a name that does not exist in column list.
   * If column list does not contain column with [name], returns [name].
   * Otherwise, tries [choices], and if the names are taken already, returns a string in a form of 'name (i)'.
   * */
  getUnusedName(name: string, choices?: string[]): string {
    return api.grok_ColumnList_GetUnusedName(this.dart, name, choices);
  }

  /** Iterates over all columns. */
  [Symbol.iterator]() : IterableIterator<Column> {
    return _getIterator(this.dart) as IterableIterator<Column>;
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}


/**
 * Row matcher.
 * */
export class RowMatcher {
  private readonly dart: any;

  constructor(dart: any) { this.dart = dart; }

  select() { api.grok_RowMatcher_Select(this.dart); }
  filter() { api.grok_RowMatcher_Filter(this.dart); }
  highlight() { api.grok_RowMatcher_Highlight(this.dart); }

  toDataFrame() { return toJs(api.grok_RowMatcher_ToDataFrame(this.dart)); }
}

/**
 * Value matcher.
 * See usage example: {@link https://public.datagrok.ai/js/samples/data-frame/value-matching/value-matcher}
 * */
export class ValueMatcher {
  private readonly dart: any;

  constructor(dart: any) { this.dart = dart; }

  static forColumn(column: Column, pattern: string) {
    return new ValueMatcher(api.grok_ValueMatcher_ForColumn(column.dart, pattern));
  }

  /** @returns {ValueMatcher} */
  static numerical(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_Numerical(pattern)); }

  /** @returns {ValueMatcher} */
  static string(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_String(pattern)); }

  /** @returns {ValueMatcher} */
  static dateTime(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_DateTime(pattern)); }

  /** @returns {ValueMatcher} */
  static bool(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_BoolMatcher(pattern)); }

  get pattern() { return api.grok_ValueMatcher_Get_Pattern(this.dart); }
  get operator() { return api.grok_ValueMatcher_Get_Operator(this.dart); }

  match(value: any) { return api.grok_ValueMatcher_Match(this.dart, value); }
  validate(value: any) { return api.grok_ValueMatcher_Validate(this.dart, value); }
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
  private readonly dart: any;
  private readonly table: any;

  constructor(table: DataFrame, dart: any) {
    /** @member {DataFrame} */
    this.table = table;
    this.dart = dart;
  }

  /** List of textual descriptions of currently applied filters */
  get filters(): DartList<string> { return DartList.fromDart(api.grok_RowList_Get_Filters(this.dart)); }

  /** Removes specified rows
   * @param {number} idx
   * @param {number} [count=1] - Number of rows to remove.
   * @param notify - Whether a change notification should be fired. */
  removeAt(idx: number, count: number = 1, notify: boolean = true): void {
    api.grok_RowList_RemoveAt(this.dart, idx, count, notify);
  }

  /** Inserts empty rows at the specified position
   * @param {number} idx
   * @param {number} [count=1] - Number of rows to insert.
   * @param notify - Whether a change notification should be fired. */
  insertAt(idx: number, count: number = 1, notify: boolean = true): void {
    api.grok_RowList_InsertAt(this.dart, idx, count, notify);
  }

  /** Appends a new row with the specified values
   * @param values - List of values (length and types should match columns)
   * @param notify - Whether a change notification should be fired.
   * @returns {Row} */
  addNew(values: any[] | null = null, notify: boolean = true): Row {
    return new Row(this.table, api.grok_RowList_AddNew(this.dart, values, notify));
  }

  /** Iterates over all rows.
   * @returns {Iterable.<Row>}
   * */
  * [Symbol.iterator]() {
    for (let i = 0; i < this.table.rowCount; i++)
      yield new Row(this.table, i);
  }

  /** Sets values for the specified row.
   * @param {number} idx - Row index.
   * @param values - List of values (length and types should match columns)
   * @param notify - Raise onDataChanged event */
  setValues(idx: number, values: any[], notify: boolean = true): void {
    api.grok_RowList_SetValues(this.dart, idx, values, notify);
  }

  /**
   * Creates a query matcher.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/row-matching/patterns}
   * @param {String|Object} query
   * @returns {RowMatcher}
   * */
  match(query: string | object): RowMatcher {
    return toJs(api.grok_RowList_Match(this.dart, query));
  }

  /**
   * @param {BitSet} bitset
   * @param {RowPredicate} rowPredicate
   * */
  _applyPredicate(bitset: BitSet, rowPredicate: RowPredicate): void {
    for (let row of this) {
      bitset.set(row.idx, rowPredicate(row), false);
    }
  }

  /** Selects rows by predicate. See {@link DataFrame.selection}
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/row-matching/select-rows}
   * */
  select(rowPredicate: RowPredicate): void {
    this._applyPredicate(this.table.selection, rowPredicate);
    this.table.selection.fireChanged();
  }

  /** Filters rows by predicate. See {@link DataFrame.filter}
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/row-matching/select-rows}
   * */
  filter(rowPredicate: RowPredicate): void {
    this._applyPredicate(this.table.filter, rowPredicate);
    this.table.filter.fireChanged();
  }

  /** Viewers that filter rows should subscribe to DataFrame.onRowsFiltering event.
   * When filtering conditions are changed, viewers should call requestFilter(). */
  requestFilter(): void {
    api.grok_RowList_RequestFilter(this.dart);
  }

  /** Adds a filter state. This should be done in the onRowsFiltering handler.
   * This is needed for filter synchronization. */
  addFilterState(state: FilterState): void {
    api.grok_RowList_AddFilterState(this.dart, state);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}

export class RowGroup {
  private readonly dart: any;

  constructor(dart: any) { this.dart = dart; }

  get dataFrame(): DataFrame {
    return toJs(api.grok_RowGroup_Get_DataFrame(this.dart));
  }
}

/** Represents a table cell. */
export class Cell {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Corresponding table.
   * @returns {DataFrame} */
  get dataFrame(): DataFrame {
    return new DataFrame(api.grok_Cell_Get_DataFrame(this.dart));
  }

  /** Corresponding row.
   * @returns {Row} */
  get row(): Row {
    return new Row(this.dataFrame, this.rowIndex);
  }

  /** Index of the corresponding row.
   * @returns {number} */
  get rowIndex(): number {
    return api.grok_Cell_Get_RowIndex(this.dart);
  }

  /** Corresponding column.
   * @returns {Column} */
  get column(): Column {
    return toJs(api.grok_Cell_Get_Column(this.dart));
  }

  /** Cell value.
   * @returns {*} */
  get value(): any { return api.grok_Cell_Get_Value(this.dart); }
  set value(x: any) { api.grok_Cell_Set_Value(this.dart, x); }

  /** Whether the cell is empty */
  isNone(): boolean { return this.column.isNone(this.rowIndex); }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}

/**
 * Efficient bit storage and manipulation.
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 */
export class BitSet {
  public dart: any;

  /** Creates a {BitSet} from the specified Dart object. */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Creates a {BitSet} from the string representing the bitset.
   * @param {string} zerosOnes - A string containing '1' and '0'.
   * @returns {BitSet} */
  static fromString(zerosOnes: string): BitSet {
    return new BitSet(api.grok_BitSet_FromString(zerosOnes));
  }

  /** Creates a {BitSet} from the ArrayBuffer representing the bitset.
   * @param {ArrayBuffer} buffer - An array containing 1 and 0.
   * @param {Number} bitLength - count of bits.
   * @returns {BitSet} */
  static fromBytes(buffer: ArrayBuffer, bitLength: number): BitSet {
    if (bitLength == null || !Number.isInteger(bitLength) || bitLength < 0)
      bitLength = buffer.byteLength * 8;
    return new BitSet(api.grok_BitSet_FromBytes(buffer, bitLength));
  }

  /** Creates a {BitSet} of the specified length with all bits set to false.
   * @param {number} length - Number of bits.
   * @param {Function} f - when specified, Sets all bits by setting i-th bit to the results of f(i)
   * @returns {BitSet} */
  static create(length: number, f?: IndexPredicate | null): BitSet {
    let bitset = new BitSet(api.grok_BitSet(length));
    if (f != null)
      bitset.init(f);
    return bitset;
  }

  /** Returns the underlying storage. Be careful with the
   * direct manipulations, as some statistics (set count, etc) are cached. */
  getBuffer(): Int32Array {
    return api.grok_BitSet_Get_Buffer(this.dart);
  }

  toBinaryString(): string {
    return api.grok_BitSet_ToBinaryString(this.dart);
  }

  /** Number of bits in a bitset
   * @type {number} */
  get length(): number {
    return api.grok_BitSet_Get_Length(this.dart);
  }

  /** Number of set bits */
  get trueCount(): number {
    return api.grok_BitSet_Get_TrueCount(this.dart);
  }

  /** Number of unset bits */
  get falseCount(): number {
    return api.grok_BitSet_Get_FalseCount(this.dart);
  }

  /** Whether any bits are set to true. */
  get anyTrue(): boolean { return this.trueCount > 0; }

  /** Whether any bits are set to false. */
  get anyFalse(): boolean { return this.falseCount > 0; }

  /** Clones a bitset
   *  @returns {BitSet} */
  clone(): BitSet {
    return new BitSet(api.grok_BitSet_Clone(this.dart));
  }

  /** Inverts a bitset.
   * @returns {BitSet} */
  invert(notify: boolean = true): BitSet {
    api.grok_BitSet_Invert(this.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise AND operation against the
   *  specified bitset. Returns this. */
  and(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_And(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise OR operation against the
   *  specified bitset. Returns this. */
  or(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_Or(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise XOR operation against the
   *  specified bitset. Returns this. */
  xor(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_Xor(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise AND_NOT operation against the
   *  specified bitset. Returns this. */
  andNot(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_AndNot(this.dart, other.dart, notify);
    return this;
  }

  /** Sets all bits to x
   * @param {boolean} x
   * @param {boolean} notify - whether BitSet's `changed` event should be fired
   * @returns {BitSet} */
  setAll(x: boolean, notify: boolean = true): BitSet {
    api.grok_BitSet_SetAll(this.dart, x, notify);
    return this;
  }

  /** Finds the first index of value x, going forward from i-th position.
   * @param {number} i - index
   * @param {boolean} x
   * @returns {number} */
  findNext(i: number, x: boolean): number {
    return api.grok_BitSet_FindNext(this.dart, i, x);
  }

  /** Finds the first index of value x, going forward from i-th position, or -1 if not found.
   * @param {number} i - Index to start searching from.
   * @param {boolean} x - Value to search for.
   * @returns {number} */
  findPrev(i: number, x: boolean): number {
    return api.grok_BitSet_FindPrev(this.dart, i, x);
  }

  /** Gets i-th bit
   * @param {number} i
   * @returns {boolean} */
  get(i: number): boolean {
    return api.grok_BitSet_GetBit(this.dart, i);
  }

  /** Sets i-th bit to x
   * @param {number} i
   * @param {boolean} x
   * @param {boolean} notify - whether BitSet's `changed` event should be fired
   * */
  set(i: number, x: boolean, notify: boolean = true): void {
    api.grok_BitSet_SetBit(this.dart, i, x, notify);
  }

  /** Sets [i]-th bit to [value], does not check bounds */
  // setFast(i: number, value: boolean): void {
  //   let buf = api.grok_BitSet_GetBuffer(this.dart);
  //   let idx = (i | 0) / 0x20;
  //
  //   if (value)
  //     buf[idx] |= 1 << (i & 0x1f);
  //   else
  //     buf[idx] &= ~(1 << (i & 0x1f));
  // }

  /** Sets all bits by setting i-th bit to the results of f(i)
   * @param {Function} f - function that accepts bit index and returns bit value
   * @param {boolean} notify - whether BitSet's `changed` event should be fired
   * @returns {BitSet} - this
   * */
  init(f: IndexPredicate, notify: boolean = true): BitSet {
    let buf = api.grok_BitSet_Get_Buffer(this.dart);
    let length = this.length;

    for (let i = 0; i < length; i++)
      buf[i] = 0;

    for (let i = 0; i < length; i++) {
      let idx = (i / 0x20) | 0;
      if (f(i))
        buf[idx] |= 1 << (i & 0x1f);
    }

    api.grok_BitSet_SetBuffer(this.dart, buf, notify);
    return this;
  }

  /** Indexes of all set bits. The result is cached.
   *  @returns {Int32Array} */
  getSelectedIndexes(): Int32Array {
    return api.grok_BitSet_GetSelectedIndexes(this.dart);
  }

  /** Copies the content from the other {BitSet}.
   * @param {BitSet} b - BitSet to copy from.
   * @returns {BitSet} */
  copyFrom(b: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_CopyFrom(this.dart, b.dart, notify);
    return this;
  }

  fireChanged(): void {
    api.grok_BitSet_FireChanged(this.dart);
  }

  /** @returns {Observable} - fires when the bitset gets changed. */
  get onChanged(): Observable<any> {
    return observeStream(api.grok_BitSet_Changed(this.dart));
  }

  /** Finds the value of similarity between two BitSets.
   * @param {BitSet} b - second BitSet.
   * @param {SimilarityMetric} metric - similarity metric to use.
   * @returns {number} */
  similarityTo(b: BitSet, metric: SimilarityMetric = SIMILARITY_METRIC.TANIMOTO): number {
    return api.grok_BitSet_SimilarityTo(this.dart, b.dart, metric);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }

  handleClick(rowIndexPredicate: IndexPredicate, mouseEvent: MouseEvent, modifiedSelectOnly: boolean = false) {
    api.grok_Utils_SelectRowsWhere(this.dart, rowIndexPredicate, mouseEvent.ctrlKey, mouseEvent.shiftKey, mouseEvent.metaKey, modifiedSelectOnly);
  }
}


/** Represents basic descriptive statistics calculated for a {Column}.
 *  See samples: {@link https://public.datagrok.ai/js/samples/data-frame/stats} */
export class Stats {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Calculates statistics for the specified column.
   * @param {Column} col
   * @param {BitSet} mask
   * @returns {Stats} */
  static fromColumn(col: Column, mask: BitSet | null = null): Stats {
    return new Stats(api.grok_Stats_FromColumn(col.dart, toDart(mask)));
  }

  /** Total number of values (including missing values). */
  get totalCount(): number {
    return api.grok_Stats_Get_TotalCount(this.dart);
  }

  /** Number of missing (empty) values. */
  get missingValueCount(): number {
    return api.grok_Stats_Get_MissingValueCount(this.dart);
  }

  /** Number of unique values. */
  get uniqueCount(): number {
    return api.grok_Stats_Get_UniqueCount(this.dart);
  }

  /** Number of non-empty values. */
  get valueCount(): number {
    return api.grok_Stats_Get_ValueCount(this.dart);
  }

  /** @returns {number} - minimum */
  get min(): number {
    return api.grok_Stats_Get_Min(this.dart);
  }

  /** @returns {number} - maximum */
  get max(): number {
    return api.grok_Stats_Get_Max(this.dart);
  }

  /** @returns {number} - sum */
  get sum(): number {
    return api.grok_Stats_Get_Sum(this.dart);
  }

  /** @returns {number} - average */
  get avg(): number {
    return api.grok_Stats_Get_Avg(this.dart);
  }

  /** @returns {number} - standard deviation */
  get stdev(): number {
    return api.grok_Stats_Get_Stdev(this.dart);
  }

  /** @returns {number} - variance */
  get variance(): number {
    return api.grok_Stats_Get_Variance(this.dart);
  }

  /** @returns {number} - skewness */
  get skew(): number {
    return api.grok_Stats_Get_Skew(this.dart);
  }

  /** @returns {number} - kurtosis */
  get kurt(): number {
    return api.grok_Stats_Get_Kurt(this.dart);
  }

  /** @returns {number} - median value */
  get med(): number {
    return api.grok_Stats_Get_Med(this.dart);
  }

  /** @returns {number} - first quartile */
  get q1(): number {
    return api.grok_Stats_Get_Q1(this.dart);
  }

  /** @returns {number} - second quartile */
  get q2(): number {
    return api.grok_Stats_Get_Q2(this.dart);
  }

  /** @returns {number} - third quartile */
  get q3(): number {
    return api.grok_Stats_Get_Q3(this.dart);
  }

  /** Pearson correlation
   * @param {Column} otherColumn
   * @returns {number} */
  corr(otherColumn: Column): number { return api.grok_Stats_Corr(this.dart, otherColumn.dart); }

  /** Spearman correlation
   * @param {Column} otherColumn
   * @returns {number} */
  spearmanCorr(otherColumn: Column): number { return api.grok_Stats_SpearmanCorr(this.dart, otherColumn.dart); }

  /** Returns distributions of [valueColumn] for each category in [catColumn]. */
  static histogramsByCategories(valueColumn: Column, catColumn: Column): Int32Array[] {
    return api.grok_Stats_HistogramsByCategories(valueColumn.dart, catColumn.dart);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}

/**
 * Fluid API for building an aggregation query against a {@link DataFrame}.
 * Build a query by calling the following methods: {@link key}, {@link pivot}, {@link count},
 * {@link uniqueCount}, {@link missingValueCount}, {@link valueCount}, {@link min}, {@link max}, {@link sum},
 * {@link avg}, {@link stdev}, {@link variance}, {@link q1}, {@link q2}, {@link q3}.
 *
 * When the query is constructed, execute it by calling {@link aggregate}, which will
 * produce a {@link DataFrame}.
 *
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 *
 * @example
 * let avgAgesByRaceAndSex = demographicsTable
 *   .groupBy(['race', 'sex'])
 *   .avg('age')
 *   .aggregate();
 */
export class GroupByBuilder {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Performs the aggregation
   *  @returns {DataFrame} */
  aggregate(): DataFrame {
    return new DataFrame(api.grok_GroupByBuilder_Aggregate(this.dart));
  }

  /**
   * Performs the aggregation
   * @param {AggregationType} agg - Aggregation type.
   * @param {string} colName - Column name.
   * @param {string} resultColName - Name of the resulting column. Default value is agg(colName).
   * @returns {GroupByBuilder}
   * */
  add(agg: AggregationType, colName?: string | null, resultColName?: string | null): GroupByBuilder {
    api.grok_GroupByBuilder_Add(this.dart, agg, colName, resultColName);
    return this;
  }

  /** Adds a key column to group values on.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  key(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.KEY, srcColName, resultColName);
  }

  /** Adds a column to pivot values on.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  pivot(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.PIVOT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts rows, including these will null values.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  count(resultColName: string = 'count'): GroupByBuilder {
    return this.add(AGG.TOTAL_COUNT, null, resultColName);
  }

  /** Adds an aggregation that counts number of unique values in the specified column.
   * See also {@link count}, {@link valueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  uniqueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.UNIQUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts number of missing values in the specified column.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  missingValueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MISSING_VALUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that counts rows, including these will null values.
   * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  valueCount(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.VALUE_COUNT, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates minimum value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  min(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MIN, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates maximum value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  max(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MAX, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates sum of the values for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  sum(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.SUM, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates median value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  med(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.MED, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates average value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  avg(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.AVG, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates standard deviation for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  stdev(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.STDEV, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates varians for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  variance(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.VARIANCE, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates first quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q1(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q1, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates second quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q2(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q2, srcColName, resultColName);
  }

  /** Adds an aggregation that calculates third quartile for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  q3(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.Q3, srcColName, resultColName);
  }

  /** Adds an aggregation that takes first value for the specified column.
   * Call {@link aggregate} when the query is constructed.
   * @param {string} srcColName - column name in the source table
   * @param {string} [resultColName] - column name in the resulting DataFrame
   * @returns {GroupByBuilder} */
  first(srcColName: string, resultColName: string | null = null): GroupByBuilder {
    return this.add(AGG.FIRST, srcColName, resultColName);
  }

  /** Gets groups of DataFrames
   * @returns {Map} - where keys are stings in format 'columnName=value' and values are DataFrames */
  getGroups(): Map<string, DataFrame> {
    return api.grok_GroupByBuilder_GetGroups(this.dart);
  }

  /**
   * Specifies the filter for the source rows.
   * @input {String|Object} pattern
   * @returns {GroupByBuilder}
   **/
  where(pattern: string | object): GroupByBuilder {
    api.grok_GroupByBuilder_Where(this.dart, pattern);
    return this;
  }

  /**
   * @param {BitSet} bitset
   * @returns {GroupByBuilder} */
  whereRowMask(bitset: BitSet): GroupByBuilder {
    if (bitset != null)
      api.grok_GroupByBuilder_WhereBitSet(this.dart, bitset.dart);
    return this;
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}

export const QNUM_LESS = 1;
export const QNUM_EXACT = 2;
export const QNUM_GREATER = 3;

let _qnumBuf = new DataView(new ArrayBuffer(8));

/**
 *  A set of static methods for working with qualified numbers.
 *  The internal representation of a qualified number is a regular double precision floating point
 *  number (IEEE 754), except the two least significant bits in mantissa are reserved
 *  for holding the qualifier ([LESS], [EXACT], [GREATER]).
 *
 *  The advantage of that representation is that the standard arithmetic operations could be
 *  performed directly on the number, without unpacking it. This is especially important for batch
 *  operations such as aggregation or sorting. While there is a loss of precision, it is rather
 *  insignificant (50 bits for storing mantissa instead of 52), which makes perfect sense
 *  considering that qualified numbers represent imprecise measurements.
 *
 *  Use [create], [getValue], and [getQ] methods for packing/unpacking.
 * */
export class Qnum {
  /**
   * Extracts the qualifier ({@link QNUM_LESS}, {@link QNUM_EXACT}, {@link QNUM_GREATER}).
   * See also {@link getValue}
   * @param {number} x
   * @returns {number}
   * */
  static getQ(x: number): number {
    _qnumBuf.setFloat64(0, x);
    return _qnumBuf.getInt8(7) & 0x03;
  }

  /**
   * Extracts the value from x, stripping the qualifier .
   * See also {@link getQ}
   * @param {number} x
   * @returns {number}
   * */
  static getValue(x: number): number {
    _qnumBuf.setFloat64(0, x);
    let last = _qnumBuf.getInt8(7) & 0xFC;
    _qnumBuf.setInt8(7, last);
    return _qnumBuf.getFloat64(0);
  }

  /**
   * Creates a QNum value out of the [value] and qualifier [q].
   * @param {number} value
   * @param {number} q
   * @returns {number}
   * */
  static create(value: number, q: number = QNUM_EXACT): number {
    _qnumBuf.setFloat64(0, value);
    let last = _qnumBuf.getInt8(7);
    _qnumBuf.setInt8(7, (last & 0xFC) | q);
    return _qnumBuf.getFloat64(0);
  }

  static exact(x: number): number {
    return Qnum.create(x, QNUM_EXACT)
  };

  static less(x: number): number {
    return Qnum.create(x, QNUM_LESS)
  };

  static greater(x: number): number {
    return Qnum.create(x, QNUM_GREATER);
  }

  /**
   * Parses a string into a qualified number.
   * @param {string} s
   * @returns {number}
   * */
  static parse(s: string): number {
    return api.grok_Qnum_Parse(s);
  }

  /**
   * Converts a qualified number to a string representation.
   * @param {number} x
   * @returns {string}
   * */
  static toString(x: number): string {
    return api.grok_Qnum_ToString(x);
  }

  /**
   * Returns the string representation of the qualifier.
   * @param {number} x
   * @returns {string}
   * */
  static qualifier(x: number): string {
    return api.grok_Qnum_Qualifier(x);
  }
}

export class DataFrameMetaHelper {
  private readonly _df: DataFrame;

  readonly formulaLines: DataFrameFormulaLinesHelper;

  async detectSemanticTypes() {
    await grok.data.detectSemanticTypes(this._df);
  }

  constructor(df: DataFrame) {
    this._df = df;
    this.formulaLines = new DataFrameFormulaLinesHelper(this._df);
  }
}

export class DataFrameFormulaLinesHelper extends FormulaLinesHelper {
  readonly df: DataFrame;

  get storage(): string { return this.df.getTag(DG.TAGS.FORMULA_LINES) ?? ''; }
  set storage(value: string) { this.df.setTag(DG.TAGS.FORMULA_LINES, value); }

  constructor(df: DataFrame) {
    super();
    this.df = df;
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
  scatter(options: object | null = null): ScatterPlotViewer { return DG.Viewer.scatterPlot(this.df, options); }
  grid(options: object | null = null): Grid { return DG.Viewer.grid(this.df, options); }
  histogram(options: object | null = null): Viewer { return DG.Viewer.histogram(this.df, options); }
  bar(options: object | null = null): Viewer { return DG.Viewer.barChart(this.df, options); }
  heatMap(options: object | null = null): Viewer { return DG.Viewer.heatMap(this.df, options); }
  box(options: object | null = null): Viewer { return DG.Viewer.boxPlot(this.df, options); }
  line(options: object | null = null): Viewer { return DG.Viewer.lineChart(this.df, options); }
  network(options: object | null = null): Viewer { return DG.Viewer.network(this.df, options); }

}

export class DataFrameDialogHelper {
  private readonly df: DataFrame;
  constructor(df: DataFrame) {
    this.df = df;
  }

  async addNewColumn(): Promise<void> { (await grok.functions.eval('AddNewColumn')).prepare().edit(); }
}

export class ColumnDialogHelper {
  private readonly column: Column;
  constructor(column: Column) {
    this.column = column;
  }

  /** Opens an editor dialog with preview for a calculated column. */
  editFormula(): void {
    let formula = this.column.getTag('formula');
    // let df = this.column.dataFrame;
    if (formula == null)
      formula = '';
    if (!(this.column.name && this.column.dataFrame?.columns.contains(this.column.name)))
      return;
    let params = { table: this.column.dataFrame, expression: formula, name: this.column.name, type: this.column.type };
    let call = DG.Func.byName('AddNewColumn').prepare(params);
    call.aux['addColumn'] = false;
    call.edit();
    let sub = grok.functions.onAfterRunAction.pipe(filter(c => c == call)).subscribe(() => {
      let newCol = call.getOutputParamValue();
      for (let [key, value] of this.column.tags) if (key !== DG.TAGS.FORMULA) newCol.setTag(key, value);
      let name = newCol.name;
      newCol = this.column.dataFrame.columns.replace(this.column, newCol);
      newCol.name = name;
      api.grok_Column_UpdateDependentColumns(newCol.dart);
      sub.unsubscribe();
    });
  }
}

export class ColumnColorHelper {
  private readonly column: Column;

  constructor(column: Column) {
    this.column = column;
  }

  getType(): ColorCodingType {
    if (this.column.tags.has(DG.TAGS.COLOR_CODING_TYPE))
      return this.column.tags[DG.TAGS.COLOR_CODING_TYPE];
    else if (this.column.tags.has(DG.TAGS.COLOR_CODING_CATEGORICAL))
      return DG.COLOR_CODING_TYPE.CATEGORICAL;
    return DG.COLOR_CODING_TYPE.OFF;
  }

  setLinear(range: ColorType[] | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.LINEAR;
    if (range != null)
      this.column.tags[DG.TAGS.COLOR_CODING_LINEAR] = JSON.stringify(range);
  }

  setCategorical(colorMap: {} | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.CATEGORICAL;
    if (colorMap != null)
      this.column.tags[DG.TAGS.COLOR_CODING_CATEGORICAL] = JSON.stringify(colorMap);
  }

  setConditional(rules: {[index: string]: number | string}  | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.CONDITIONAL;
    if (rules != null) {
      for (let [rule, color] of Object.entries(rules)) {
        if (typeof color === 'number')
          rules[rule] = DG.Color.toHtml(color);
      }
      this.column.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = JSON.stringify(rules);
    }
  }

  getColor(i: number): number {
    return api.grok_Column_GetColor(this.column.dart, i);
  }

  getColors(): Uint32Array {
    return api.grok_Column_GetColors(this.column.dart);
  }
}

export class ColumnMarkerHelper {
  private readonly column: Column;

  constructor(column: Column) {
    this.column = column;
  }

  assign(category: string, marker: MarkerCodingType): ColumnMarkerHelper {
    this.setMarkerCoding(category, marker);
    return this;
  }

  default(marker: MarkerCodingType): ColumnMarkerHelper {
    return this.assign('~DEFAULT', marker);
  }

  reset(): ColumnMarkerHelper {
    this.column.tags[DG.TAGS.MARKER_CODING] = '{}';
    return this;
  }

  setMarkerCoding(category: string, marker: string): void {
    let jsonTxt: string | null = this.column.getTag(DG.TAGS.MARKER_CODING);
    let jsonMap: {[key: string]: string} = jsonTxt ? JSON.parse(jsonTxt) : {};
    jsonMap[category] = marker;
    jsonTxt = JSON.stringify(jsonMap);
    this.column.setTag(DG.TAGS.MARKER_CODING, jsonTxt);
  }
}

export class ColumnMetaHelper {
  private readonly column: Column;
  private _colors: ColumnColorHelper | undefined;
  private _markers: ColumnMarkerHelper | undefined;

  constructor(column: Column) {
    this.column = column;
  }

  get colors(): ColumnColorHelper {
    if (this._colors == undefined)
      this._colors = new ColumnColorHelper(this.column);
    return this._colors;
  }

  get markers(): ColumnMarkerHelper {
    if (this._markers == undefined)
      this._markers = new ColumnMarkerHelper(this.column);
    return this._markers;
  }

  /** Returns the format of the dataframe column. See also [GridColumn.format] */
  get format(): string | null {
    return this.column.getTag(TAGS.FORMAT) ?? api.grok_Column_GetAutoFormat(this.column.dart);
  }
}
