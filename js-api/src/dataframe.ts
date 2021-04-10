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
  IndexPredicate
} from "./const";
import {__obs, observeStream} from "./events";
import {toDart, toJs} from "./wrappers";
import {SIMILARITY_METRIC} from "./const";
import {_getIterator, _toIterable} from "./utils";
import {Observable}  from "rxjs";

declare let grok: any;
declare let DG: any;
let api = <any>window;
type RowPredicate = (row: Row) => boolean;

/**
 * Finds the item by its unique id.
 *
 * @typedef {function(Row): boolean} RowPredicate
 * @typedef {function(Column): boolean} ColumnPredicate
 */

/**
 * Proxies a Dart Map, API-compliant to ES 2015+
 */

const MapProxy = new Proxy(class {
    d: any;
    constructor(d: any) {
      this.d = d;
    }
    keys(): Iterable<any> {
      return _toIterable(api.grok_Map_Keys(this.d));
    }
    values() {
      return _toIterable(api.grok_Map_Values(this.d));
    }
    * [Symbol.iterator] () {
      for (let key of this.keys()) {
        const value = DG.toJs(api.grok_Map_Get(this.d, key));
        yield [key, value];
      }
    }
    entries() {
      return this;
    }
    forEach(callback: (key: string, value: any) => void) {
      for (const [key, value] of this) {
        callback(key, value);
      }
    }
    delete(key: string) {
      return DG.toJs(api.grok_Map_Delete(this.d, DG.toDart(key)));
    }
    get(key: any) {
      return DG.toJs(api.grok_Map_Get(this.d, key));
    }
    has(key: string) {
      return api.grok_Map_Has(this.d, DG.toDart(key));
    }
    set(key: string, value: any) {
      api.grok_Map_Set(this.d, key, DG.toDart(value));
      return this;
    }
    clear() {
      api.grok_Map_Clear(this.d);
    }
    size() {
      return DG.toJs(api.grok_Map_Size(this.d));
    }
  }, {
    construct(target, args) {
      // @ts-ignore
      return new Proxy(new target(...args), {
        get: function (target: any, prop) {
          const val = target[prop];
          if (typeof val === 'function') {
            return function (...args :string[]) {
              return val.apply(target, args);
            };
          } /* else if (typeof val === 'number') {
            // for the case of a getter like name
            return val;
          } */ else {
            return DG.toJs(api.grok_Map_Get(target.d, prop));
          }
        },
        set: function (target, prop, value) {
          api.grok_Map_Set(target.d, prop, DG.toDart(value));
          return true;
        },
        deleteProperty: function (target, prop) {
          api.grok_Map_Delete(target.d, DG.toDart(prop));
          return true;
        },
        has: function (target, prop) {
          return api.grok_Map_Has(target.d, DG.toDart(prop));
        },
        getOwnPropertyDescriptor(target, prop) {
          return {
            enumerable: true,
            configurable: true
          };
        },
        ownKeys: function (target) {
          return Array.from(target.keys());
        }
      });
    }
  }
);

/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 *
 * In the API, the terms "Table" and "DataFrame" are used interchangeably.
 * See usage samples: https://public.datagrok.ai/js/samples/data-frame/manipulate
 */
export class DataFrame {
  public readonly d: any;
  public columns: ColumnList;
  public rows: RowList;
  public filter: BitSet;
  public temp: any;
  public tags: any;

  constructor(d: any) {
    this.d = d;
    this.columns = toJs(api.grok_DataFrame_Columns(this.d));
    this.rows = new RowList(this, api.grok_DataFrame_Rows(this.d));
    this.filter = new BitSet(api.grok_DataFrame_Get_Filter(this.d));

    this.temp = new MapProxy(api.grok_DataFrame_Get_Temp(this.d));
    this.tags = new MapProxy(api.grok_DataFrame_Get_Tags(this.d));

    // return new Proxy(this, {
    //     get(target, name) {
    //         if (target.hasOwnProperty(name))
    //             return target[name];
    //         return target.table.get(name, target.idx);
    //     }
    // });
  }

  /** Creates a {@link DataFrame} with the specified number of rows and no columns.
   * @param {number} rowCount
   * @returns {DataFrame} */
  static create(rowCount: number = 0): DataFrame {
    return new DataFrame(api.grok_DataFrame(rowCount));
  }

  /** Creates a {@link DataFrame} from the specified columns. All columns should be of the same length.
   * @param {Column[]} columns
   * @returns {DataFrame} */
  static fromColumns(columns: Column[]): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromColumns(columns.map((c) => c.d)));
  }

  /** Constructs {@link DataFrame} from a comma-separated values string
   * @param {string} csv - The content of the comma-separated values file.
   * @param {CsvImportOptions} options
   * @returns {DataFrame} */
  static fromCsv(csv: string, options?: CsvImportOptions): DataFrame {
    return grok.data.parseCsv(csv, options);
  }

  /** Constructs {@link DataFrame} from the specified JSON string.
   * @param {string} json - JSON document.
   * @returns {DataFrame} */
  static fromJson(json: string): DataFrame {
    return new DataFrame(api.grok_DataFrame_FromJson(json));
  }

  /** Returns number of rows in the table.
   * @returns {number} */
  get rowCount(): number {
    return api.grok_DataFrame_RowCount(this.d);
  }

  /** Returns a {@link BitSet} with selected rows.
   * @returns {BitSet} */
  get selection(): BitSet {
    return new BitSet(api.grok_DataFrame_Get_Selection(this.d));
  }

  /** Name of the dataframe.
   * @returns {string}*/
  get name(): string {
    return api.grok_DataFrame_Get_Name(this.d);
  }

  set name(s: string) {
    api.grok_DataFrame_Set_Name(this.d, s);
  }

  /** Returns the value of the specified tag, or null if it does not exist.
   * @returns {string} */
  getTag(tag: string): string | null {
    return api.grok_DataFrame_Get_Tag(this.d, tag);
  }

  /** Sets a tag to the specified value.
   * @param {string} tag - Key.
   * @param {string} value - Value. */
  setTag(tag: string, value: string): void {
    api.grok_DataFrame_Set_Tag(this.d, tag, value);
  }

  /** Returns i-th row.
   * @param {number} i - Row index.
   * @returns {Row} */
  row(i: number): Row {
    return new Row(this, i);
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
    return toJs(api.grok_DataFrame_ColumnByName(this.d, name));
  }

  /** Returns a {@link Cell} with the specified name.
   * @param {number} idx - Row index.
   * @param {string} name - Column name.
   * @returns {Cell} */
  cell(idx: number, name: string): Cell {
    return new Cell(api.grok_DataFrame_Cell(this.d, idx, name));
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

  /** Exports the content to comma-separated-values format.
   * @returns {string} */
  toCsv(): string {
    return api.grok_DataFrame_ToCsv(this.d);
  }

  /** Creates a new dataframe from the specified row mask and a list of columns.
   * @param {BitSet} rowMask - Rows to include.
   * @param {string[]} columnIds - Columns to include.
   * @param {boolean} saveSelection - Whether selection should be saved. */
  clone(rowMask: BitSet | null = null, columnIds: string[] | null = null, saveSelection: boolean = false): DataFrame {
    return new DataFrame(api.grok_DataFrame_Clone(this.d, toDart(rowMask), columnIds, saveSelection));
  }

  /** Current row
   * @type {Row} */
  get currentRow(): Row {
    return new Row(this, api.grok_DataFrame_Get_CurrentRowIdx(this.d));
  }

  set currentRow(idx) {
    api.grok_DataFrame_Set_CurrentRowIdx(this.d, idx);
  }

  /** Current column
   * @type {Column} */
  get currentCol(): Column {
    return toJs(api.grok_DataFrame_Get_CurrentCol(this.d));
  }

  set currentCol(col: Column) {
    api.grok_DataFrame_Set_CurrentCol(this.d, col.d);
  }

  /** Current cell
   * @type {Cell} */
  get currentCell(): Cell {
    return new Cell(api.grok_DataFrame_Get_CurrentCell(this.d));
  }

  set currentCell(cell: Cell) {
    api.grok_DataFrame_Set_CurrentCell(this.d, cell.d);
  }

  /** Creates a [DataFrame] from list of objects by using object keys as column names,
   * and object values as values.
   * @param {object[]} list - List of objects. **/
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

  /** Converts a column with the specified name to [newType],
   * removes the original column from its dataframe and adds the new column to it.
   * @param {string|Column} column
   * @param {string} newType
   * @param {string=} format
   * @returns {Column} */
  changeColumnType(column: string | Column, newType: ColumnType, format: string | null = null): Column {
    return toJs(api.grok_DataFrame_ChangeColumnType(this.d, toDart(column), newType, format));
  }

  /**
   * Returns [Int32Array] that contains sorted order, or null for unsorted (original) order.
   * @param {Object[]} sortByColumnIds - Collection of [Column]s to use as keys for sorting.
   * @param {boolean[]} sortOrders - List of sort orders for [sortByCols]. True == ascending.
   * @param {BitSet} rowMask - Mask of the rows to sort. Result array will contain [rowIndexes.length] elements.
   * @returns {Int32Array}
   * */
  getSortedOrder(sortByColumnIds: object[], sortOrders: boolean[] | null = null, rowMask: BitSet | null = null): Int32Array {
    return api.grok_DataFrame_GetSortedOrder(this.d, sortByColumnIds.map(toDart), sortOrders, toDart(rowMask));
  }

  /** Begins building a query, using the specified columns as keys.
   * @param {string[]} columnNames - Names of the columns to be used as keys.
   * @returns {GroupByBuilder}
   *  */
  groupBy(columnNames: string[] = []): GroupByBuilder {
    return new GroupByBuilder(api.grok_DataFrame_GroupBy(this.d, columnNames));
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
    return new DataFrame(api.grok_DataFrame_Unpivot(this.d, copyColumnNames, mergeColumnNames, categoryColumnName, valueColumnName));
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
    return new DataFrame(api.grok_JoinTables(this.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }

  /**
   * Appends two tables ('union' in SQL).
   * @param {DataFrame} t2
   * @param {boolean} inPlace - whether to create a new table, or modify 'this' one.
   * @param {String[]} columnsToAppend
   * @returns {DataFrame}
   * */
  append(t2: DataFrame, inPlace: boolean = false, columnsToAppend: string[] | null = null): DataFrame {
    return new DataFrame(api.grok_DataFrame_Append(this.d, t2.d, inPlace, columnsToAppend));
  }

  /** @returns {Observable} */ _event(event: string): Observable<any> {
    return __obs(event, this.d);
  }

  /** @returns {Observable} */ get onValuesChanged(): Observable<any> {
    return this._event('ddt-values-changed');
  }

  /** @returns {Observable} */ get onCurrentRowChanged(): Observable<any> {
    return this._event('ddt-current-row-changed');
  }

  /** @returns {Observable} */ get onMouseOverRowChanged(): Observable<any> {
    return this._event('ddt-mouse-over-row-changed');
  }

  /** @returns {Observable} */ get onCurrentColChanged(): Observable<any> {
    return this._event('ddt-current-col-changed');
  }

  /** @returns {Observable} */ get onMouseOverColChanged(): Observable<any> {
    return this._event('ddt-mouse-over-col-changed');
  }

  /** @returns {Observable} */ get onCurrentCellChanged(): Observable<any> {
    return this._event('ddt-current-cell-changed');
  }

  /** @returns {Observable} */ get onMouseOverRowGroupChanged(): Observable<any> {
    return this._event('ddt-mouse-over-row-group-changed');
  }

  /** @returns {Observable} */ get onNameChanged(): Observable<any> {
    return this._event('ddt-table-name-changed');
  }

  /** @returns {Observable} */ get onMetadataChanged(): Observable<any> {
    return this._event('ddt-table-metadata-changed');
  }

  /** @returns {Observable} */ get onColumnNameChanged(): Observable<any> {
    return this._event('ddt-table-column-name-changed');
  }

  /** @returns {Observable} */ get onColumnSelectionChanged(): Observable<any> {
    return this._event('ddt-column-selection-changed');
  }

  /** @returns {Observable} */ get onColumnsChanged(): Observable<any> {
    return this._event('ddt-columns-changed');
  }

  /** @returns {Observable} */ get onColumnsAdded(): Observable<any> {
    return this._event('ddt-columns-added');
  }

  /** @returns {Observable} */ get onColumnsRemoved(): Observable<any> {
    return this._event('ddt-columns-removed');
  }

  /** @returns {Observable} */ get onRowsAdded(): Observable<any> {
    return this._event('ddt-rows-added');
  }

  /** @returns {Observable} */ get onRowsRemoved(): Observable<any> {
    return this._event('ddt-rows-removed');
  }

  /** @returns {Observable} */ get onRowsFiltering(): Observable<any> {
    return this._event('ddt-rows-filtering');
  }

  /** @returns {Observable} */ get onSemanticTypeDetecting(): Observable<any> {
    return this._event('ddt-semantic-type-detecting');
  }

  /** @returns {Observable} */ get onSemanticTypeDetected(): Observable<any> {
    return this._event('ddt-semantic-type-detected');
  }

  /** @returns {Observable} */ get onDataChanged(): Observable<any> {
    return rxjs.concat(this.onValuesChanged, this.onColumnsAdded,
      this.onColumnsRemoved, this.onRowsAdded, this.onRowsRemoved);
  }

  /** @returns {Observable} */ get onSelectionChanged(): Observable<any> {
    return this.selection.onChanged;
  }

  /** @returns {Observable} */ get onFilterChanged(): Observable<any> {
    return this.filter.onChanged;
  }

  fireValuesChanged(): void {
    api.grok_DataFrame_FireValuesChanged(this.d);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
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
    return api.grok_MathActions_GetDensity(this.d, xBins, yBins, xColName, yColName);
  }
}

/** Represents a row. Allows for quick property access like "row.height". */
export class Row {
  private table: DataFrame;
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
        if (target.hasOwnProperty(name))
          return target[<any>name];
        return target.table.get(name, target.idx);
      }
    });
  }

  /** Returns this row's value for the specified column
   * @param {string} columnName
   * @returns {Object} */
  get(columnName: string): any {
    return this.table.getCol(columnName).get(this.idx);
  }
}

/** Strongly-typed column. */
export class Column {
  public d: any;
  private temp: any;
  private tags: any;

  constructor(d: any) {
    this.d = d;
    this.temp = new MapProxy(api.grok_Column_Get_Temp(this.d));
    this.tags = new MapProxy(api.grok_Column_Get_Tags(this.d));
    //
    // return new Proxy(this, {
    //     get(target, x) {
    //         if (typeof x === 'number')
    //             return target.get(x);
    //         if (target.hasOwnProperty(x))
    //             return target[x];
    //     }
    // });
  }

  static fromStrings(name: string, list: string[]): Column {
    return toJs(api.grok_Column_FromStrings(name, list));
  }

  static fromType(type: ColumnType, name?: string | null, length: number = 0): Column {
    return toJs(api.grok_Column_FromType(type, name, length));
  }

  /** [array] will be not be copied and will be used as column's storage */
  static fromInt32Array(name: string, array: Int32Array, length: number | null = null): Column {
    return toJs(api.grok_Column_FromInt32Array(name, array, length));
  }

  /** [array] will be not be copied and will be used as column's storage */
  static fromFloat32Array(name: string, array: Float32Array, length: number | null = null): Column {
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
  static fromBitSet(name: string, bitset: BitSet): Column {
    return toJs(api.grok_Column_FromBitSet(name, bitset.d));
  }

  /** Creates an integer column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static int(name: string, length: number = 0): Column {
    return Column.fromType(TYPE.INT, name, length);
  }

  /** Creates a floating point column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static float(name: string, length: number = 0): Column {
    return Column.fromType(TYPE.FLOAT, name, length);
  }

  /** Creates a string column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static string(name: string, length: number = 0): Column {
    return Column.fromType(TYPE.STRING, name, length);
  }

  /** Creates a boolean column with the specified name and length.
   * @param {string} name
   * @param {number} length
   * @returns {Column} */
  static bool(name: string, length: number = 0): Column {
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
  static dataFrame(name: string, length: number = 0): Column {
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
  static qnum(name: string, length: number = 0, values: number[] = [], exact: boolean = true): Column {
    let col = Column.fromType(TYPE.QNUM, name, length);
    if (values !== null) {
      let buffer = col.getRawData();
      for (let i = 0; i < length; i++)
        buffer[i] = exact ? Qnum.exact(values[i]):  values[i];
      col.setRawData(buffer);
    }
    return col;
  }

  /** Column data type.
   * @type {string} */
  get type(): ColumnType {
    return api.grok_Column_Get_Type(this.d);
  }

  /** Number of elements
   * @type {number} */
  get length(): number {
    return api.grok_Column_Get_Length(this.d);
  }

  /** Parent table
   * @type {DataFrame} */
  get dataFrame(): DataFrame {
    return toJs(api.grok_Column_Get_DataFrame(this.d));
  }

  /** Semantic type
   * @type {string} */
  get semType(): SemType {
    return api.grok_Column_Get_SemType(this.d);
  }

  set semType(s: SemType) {
    api.grok_Column_Set_SemType(this.d, s);
  }

  /** Layout column ID
   @type {string} */
  get layoutColumnId(): string {
    return api.grok_Column_Get_LayoutColumnId(this.d);
  }

  set layoutColumnId(s: string) {
    api.grok_Column_Set_LayoutColumnId(this.d, s);
  }

  /** @type {string} */
  get name(): string {
    return api.grok_Column_Get_Name(this.d);
  }

  set name(s: string) {
    api.grok_Column_Set_Name(this.d, s);
  }

  /**
   * Initializes all values in the column to [columnInitializer].
   * @param {string | number | Function} valueInitializer
   * @returns {Column}
   * */
  init(valueInitializer: string | number | ((ind: number) => any)): Column {
    let type = typeof valueInitializer;
    if (type === 'function')
      api.grok_Column_Init(this.d, valueInitializer);
    else if (type === 'number' || type === 'string')
      api.grok_Column_SetAllValues(this.d, valueInitializer);
    return this;
  }

  /** Returns the raw buffer containing data. Return type depends on the column type:
   * {Int32Array} for ints, {@link INT_NULL} represents null.
   * {Float32Array} for floats, {@link FLOAT_NULL} represents null.
   * {Float64Array} for qnums, {@link FLOAT_NULL} represents null.
   * {Float64Array} for datetime, in microseconds since epoch, {@link FLOAT_NULL} represents null.
   * {Int32Array} for strings indexes of {@link categories}.
   * {Uint32Array} bit array.
   * @returns {Array} */
  getRawData(): Int32Array | Float32Array | Float64Array | Uint32Array {
    return api.grok_Column_GetRawData(this.d);
  }

  setRawData(rawData: Int32Array | Float32Array | Float64Array | Uint32Array, notify: boolean = true): void {
    api.grok_Column_SetRawData(this.d, rawData, notify);
  }

  /** Gets i-th value
   * @param {number} row - row index
   * @returns {object} - or null if isNone(i) */
  get(row: number): any {
    return api.grok_Column_GetValue(this.d, row);
  }

  // /** Returns i-th value as integer. Throws exception if column's type is not TYPE_INT or TYPE_DATE_TIME.
  //  * @param {number} i
  //  * @returns {number} */
  // getInt(i);
  //
  // /** Returns i-th value as a integer. Works for IntColumns only.
  //  * @param {number} i
  //  * @returns {number} */
  // getFloat(i)
  //
  // /** Returns i-th value as boolean. Works for IntColumns only.
  //  * @param {number} i
  //  * @returns {number} */
  // getBool(i)
  //
  // /** Returns i-th value as integer. Works for IntColumns only.
  //  * @param {number} i
  //  * @returns {number} */
  // getDateTime(i)

  /** Returns i-th value as string, taking into account value format defined for the column.
   *  An empty string is returned if there is no value.
   * @param {number} i
   * @returns {string} */
  getString(i: number): string {
    return api.grok_Column_GetString(this.d, i);
  }

  /** Attempts to set i-th value by converting a provided string to the corresponding strongly-typed value.
   *  Returns true if text was successfully parsed and set, otherwise false.
   *  Examples: dateColumn.setString('April 1, 2020');
   *            intColumn.setString('42');
   * @param {number} i
   * @param {string} str
   * @param {boolean} notify
   * @returns {boolean} */
  setString(i: number, str: string, notify: boolean = true): void {
    api.grok_Column_SetString(this.d, i, str, notify);
  }

  /**
   * Sets [i]-th value to [x]
   * @param {number} i
   * @param x
   * @param {boolean} notify
   */
  set(i: number, x: any, notify: boolean = true): void {
    api.grok_Column_SetValue(this.d, i, toDart(x), notify);
  }

  /** Returns whether i-th value is missing.
   * @param {number} i - Row index.
   * @returns {boolean} */
  isNone(i: number): boolean {
    return api.grok_Column_IsNone(this.d, i);
  }

  /** Gets the value of the specified tag.
   *  @param {string} tag
   *  @returns {string} */
  getTag(tag: string): string {
    return api.grok_Column_Get_Tag(this.d, tag);
  }

  /** Sets a tag to the specified value.
   * @param {string} tag - Key.
   * @param {string} value - Value.
   * @returns {Column}. * */
  setTag(tag: string, value: string): Column {
    api.grok_Column_Set_Tag(this.d, tag, value);
    return this;
  }

  /** Compacts the internal column representation.
   *  Currently, it only affects string columns where values were modified. */
  compact() {
    return api.grok_Column_Compact(this.d);
  }

  /** Copies column values to an array.
   *  @returns {Array} */
  toList(): Array<any> {
    return api.grok_Column_ToList(this.d);
  }

  /** Returns all unique strings in a sorted order. Applicable to string column only.
   * @returns {string[]} */
  get categories(): string[] {
    return api.grok_Column_Categories(this.d);
  }

  /** Sets order of categories
   * @param {string[]} order */
  setCategoryOrder(order: string[]) {
    api.grok_Column_SetCategoryOrder(this.d, order);
  }

  /** Gets order of categories
   * @returns string[] */
  getCategoryOrder() {
    return api.grok_Column_GetCategoryOrder(this.d);
  }

  /** Column's minimum value. The result is cached.
   * @returns {number} */
  get min(): number {
    return api.grok_Column_Min(this.d);
  }

  /** Column's maximum value. The result is cached.
   * @returns {number} */
  get max(): number {
    return api.grok_Column_Max(this.d);
  }

  /** Checks whether the column passes the specified [filter].
   * [filter] can be either specific data [type] such as 'int' or 'string', more broadly - 'numerical', or 'categorical', or null for any columns.
   * @returns {boolean} */
  matches(filter: ColumnType | 'numerical' | 'categorical' | null): boolean {
    return api.grok_Column_Matches(this.d, filter);
  }

  /** Basic descriptive statistics. The result is cached.
   * @returns {Stats} */
  get stats(): Stats {
    return Stats.fromColumn(this);
  }

  /** An iterator over all values in this column. */
  * values() {
    for (let i = 0; i < this.length; i++) {
      yield this.get(i);
    }
  }

  /** Creates and returns a new column by converting [column] to the specified [newType].
   *  @param {string} newType
   *  @param {string} format
   *  @returns {Column} */
  convertTo(newType: string, format: string | null = null): Column {
    return toJs(api.grok_Column_ConvertTo(this.d, newType, format));
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
  }
}

/** Columns in a [DataFrame]. */
export class ColumnList {
  private readonly d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Number of elements in the column.
   * @returns {number} */
  get length(): number {
    return api.grok_ColumnList_Length(this.d);
  }

  /** Column with the corresponding name (case-insensitive).
   * @param {string} name - Column name
   * @returns {Column} */
  byName(name: string): Column {
    return toJs(api.grok_ColumnList_ByName(this.d, name));
  }

  /** Maps names to columns.
   * @param {string[]} names - Column names
   * @returns {Column[]} */
  byNames(names: string[]): Column[] {
    return names.map(name => this.byName(name));
  }

  /** Column by index.
   * @param {number} index - Index of the column.
   * @returns {Column} */
  byIndex(index: number): Column {
    return toJs(api.grok_ColumnList_ByIndex(this.d, index));
  }

  /** First column of [semType], or null.
   * @returns {Column} */
  bySemType(semType: SemType): Column | null {
    let col = api.grok_ColumnList_BySemType(this.d, semType);
    return col == null ? null : toJs(col);
  }

  /** Finds columns by the corresponding semTypes, or null, if any of the sem types could not be found.
   * @returns {Column[]} */
  bySemTypesExact(semTypes: SemType[]): Column[] | null {
    let columns = [];
    for (let semType of semTypes) {
      let col = this.bySemType(semType);
      if (col == null)
        return null;
      columns.push(col);
    }
    return columns;
  }

  /** @returns {Iterable.<Column>} */
  byTags(tags: string[]): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_ByTags(this.d, tags));
  }

  /** @returns {Iterable.<Column>} */
  get categorical(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_Categorical(this.d));
  }

  /** @returns {Iterable.<Column>} */
  get numerical(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_Numerical(this.d));
  }

  /** Array containing column names.
   * @returns {string[]} */
  names(): string[] {
    return api.grok_ColumnList_Names(this.d);
  }

  /** Creates an array of columns.
   * @returns {Column[]} */
  toList(): Column[] {
    return this.names().map((name: string) => this.byName(name));
  }

  /** Adds a column, and optionally notifies the parent dataframe.
   * @param {Column} column
   * @param {boolean} notify
   * @returns {Column} */
  add(column: Column, notify: boolean = true): Column {
    api.grok_ColumnList_Add(this.d, column.d, notify);
    return column;
  }

  /** Inserts a column, and optionally notifies the parent dataframe.
   * @param {Column} column
   * @param {boolean} notify
   * @param {int} index
   * @returns {Column} */
  insert(column: Column, notify: boolean = true, index: number): Column {
    api.grok_ColumnList_Insert(this.d, column.d, notify, index);
    return column;
  }

  /** Adds an empty column of the specified type.
   * @param {string} name
   * @param {ColumnType} type
   * @returns {Column} */
  addNew(name: string, type: ColumnType): Column {
    return toJs(api.grok_ColumnList_AddNew(this.d, name, type));
  }

  /** Adds calculated column.
   * @param {string} name
   * @param {string} expression
   * @param {ColumnType} type
   * @param {bool} treatAsString
   * @returns {Column} */
  addNewCalculated(name: string, expression: string, type: ColumnType | null, treatAsString: boolean | null): Promise<Column> {
    return new Promise((resolve, reject) => toJs(api.grok_ColumnList_AddNewCalculated(this.d, name, expression, type, treatAsString, (c: any) => resolve(c), (e: any) => reject(e))));
  }

  /** Adds a string column
   * @param {string} name
   * @returns {Column} */
  addNewString(name: string): Column { return this.addNew(name, TYPE.STRING); }

  /** Adds a new integer column
   * @param {string} name
   * @returns {Column} */
  addNewInt(name: string): Column { return this.addNew(name, TYPE.INT); }

  /** Adds a new float column
   * @param {string} name
   * @returns {Column} */
  addNewFloat(name: string): Column { return this.addNew(name, TYPE.FLOAT); }

  /** Adds a new qualified number column
   * @param {string} name
   * @returns {Column} */
  addNewQnum(name: string): Column { return this.addNew(name, TYPE.QNUM); }

  /** Adds a new datetime column
   * @param {string} name
   * @returns {Column} */
  addNewDateTime(name: string): Column { return this.addNew(name, TYPE.DATE_TIME); }

  /** Adds a new boolean column
   * @param {string} name
   * @returns {Column} */
  addNewBool(name: string): Column { return this.addNew(name, TYPE.BOOL); }

  /** Adds a virtual column.
   * @param {string} name
   * @param {Function} getValue - value constructor function that accepts int index and returns value
   * @param {String} type - column type
   * @returns {Column} */
  addNewVirtual(name: string, getValue: (ind: number) => any, type = TYPE.OBJECT): Column {
    return toJs(api.grok_ColumnList_AddNewVirtual(this.d, name, getValue, type));
  }

  /** Removes column by name (case-insensitive).
   * @param {string} column
   * @returns {ColumnList} */
  remove(column: string): ColumnList {
    api.grok_ColumnList_Remove(this.d, column);
    return this;
  }

  /** Checks whether this list contains a column with the specified name. The check is case-insensitive.
   * @returns {boolean} */
  contains(columnName: string): boolean {
    return api.grok_ColumnList_Contains(this.d, columnName);
  }

  /** Replaces the column with the new column.
   * @param {Column} columnToReplace
   * @param {Column} newColumn */
  replace(columnToReplace: Column | string, newColumn: Column): void {
    api.grok_ColumnList_Replace(this.d, (typeof columnToReplace === 'string') ? columnToReplace:  columnToReplace.d, newColumn.d);
  }

  /** Iterates over all columns.
   * @returns {Iterable.<Column>}
   * */
  [Symbol.iterator]() {
    return _getIterator(this.d);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
  }
}


/**
 * Row matcher.
 * */
export class RowMatcher {
  private readonly d: any;

  constructor(d: any) { this.d = d; }

  select() { api.grok_RowMatcher_Select(this.d); }
  filter() { api.grok_RowMatcher_Filter(this.d); }
  highlight() { api.grok_RowMatcher_Highlight(this.d); }

  toDataFrame() { return toJs(api.grok_RowMatcher_ToDataFrame(this.d)); }
}

/**
 * Value matcher.
 * */
export class ValueMatcher {
  private readonly d: any;

  constructor(d: any) { this.d = d; }

  static forColumn(column: Column, pattern: string) {
    return new ValueMatcher(api.grok_ValueMatcher_ForColumn(column.d, pattern));
  }

  /** @returns {ValueMatcher} */
  static numerical(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_Numerical(pattern)); }

  /** @returns {ValueMatcher} */
  static string(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_String(pattern)); }

  /** @returns {ValueMatcher} */
  static dateTime(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_DateTime(pattern)); }

  /** @returns {ValueMatcher} */
  static bool(pattern: string) { return new ValueMatcher(api.grok_ValueMatcher_BoolMatcher(pattern)); }

  get pattern() { return api.grok_ValueMatcher_Get_Pattern(this.d); }
  get operator() { return api.grok_ValueMatcher_Get_Operator(this.d); }

  match(value: any) { return api.grok_ValueMatcher_Match(this.d, value); }
  validate(value: any) { return api.grok_ValueMatcher_Validate(this.d, value); }
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
  private readonly d: any;
  private readonly table: any;

  constructor(table: DataFrame, d: any) {
    /** @member {DataFrame} */
    this.table = table;
    this.d = d;
  }

  /** Removes specified rows
   * @param {number} idx
   * @param {number} [count=1] - Number of rows to remove.
   * @param notify - Whether a change notification should be fired. */
  removeAt(idx: number, count: number = 1, notify: boolean = true): void {
    api.grok_RowList_RemoveAt(this.d, idx, count, notify);
  }

  /** Inserts empty rows at the specified position
   * @param {number} idx
   * @param {number} [count=1] - Number of rows to insert.
   * @param notify - Whether a change notification should be fired. */
  insertAt(idx: number, count: number = 1, notify: boolean = true): void {
    api.grok_RowList_InsertAt(this.d, idx, count, notify);
  }

  /** Appends a new row with the specified values
   * @param values - List of values (length and types should match columns)
   * @param notify - Whether a change notification should be fired.
   * @returns {Row} */
  addNew(values: any[] | null = null, notify: boolean = true): Row {
    return new Row(this.table, api.grok_RowList_AddNew(this.d, values, notify));
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
   * @param values - List of values (length and types should match columns) */
  setValues(idx: number, values: any[]): void {
    api.grok_RowList_SetValues(this.d, idx, values);
  }

  /**
   * Creates a query matcher.
   * @param {String|Object} query
   * @returns {RowMatcher}
   * */
  match(query: string | object): RowMatcher {
    return toJs(api.grok_RowList_Match(this.d, query));
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

  /**
   * @param {RowPredicate} rowPredicate
   * */
  select(rowPredicate: RowPredicate): void {
    this._applyPredicate(this.table.selection, rowPredicate);
  }

  /**
   * @param {RowPredicate} rowPredicate
   * */
  filter(rowPredicate: RowPredicate): void {
    this._applyPredicate(this.table.filter, rowPredicate);
  }

  /** Viewers that filter rows should subscribe to DataFrame.onRowsFiltering event.
   * When filtering conditions are changed, viewers should call requestFilter(). */
  requestFilter(): void {
    api.grok_RowList_RequestFilter(this.d);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
  }
}

/** Represents a table cell. */
export class Cell {
  public d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Corresponding table.
   * @returns {DataFrame} */
  get dataFrame(): DataFrame {
    return new DataFrame(api.grok_Cell_Get_DataFrame(this.d));
  }

  /** Corresponding row.
   * @returns {Row} */
  get row(): Row {
    return new Row(this.dataFrame, this.rowIndex);
  }

  /** Index of the corresponding row.
   * @returns {number} */
  get rowIndex(): number {
    return api.grok_Cell_Get_RowIndex(this.d);
  }

  /** Corresponding column.
   * @returns {Column} */
  get column(): Column {
    return toJs(api.grok_Cell_Get_Column(this.d));
  }

  /** Cell value.
   * @returns {*} */
  get value(): any {
    return api.grok_Cell_Get_Value(this.d);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
  }
}

/**
 * Efficient bit storage and manipulation.
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 */
export class BitSet {
  public d: any;

  /** Creates a {BitSet} from the specified Dart object. */
  constructor(d: any) {
    this.d = d;
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

  toBinaryString(): string {
    return api.grok_BitSet_ToBinaryString(this.d);
  }

  /** Number of bits in a bitset
   * @type {number} */
  get length(): number {
    return api.grok_BitSet_Get_Length(this.d);
  }

  /** Number of set bits
   * @type {number} */
  get trueCount(): number {
    return api.grok_BitSet_Get_TrueCount(this.d);
  }

  /** Number of unset bits
   * @type {number}*/
  get falseCount(): number {
    return api.grok_BitSet_Get_FalseCount(this.d);
  }

  /** Clones a bitset
   *  @returns {BitSet} */
  clone(): BitSet {
    return new BitSet(api.grok_BitSet_Clone(this.d));
  }

  /** Inverts a bitset.
   * @returns {BitSet} */
  invert(): BitSet {
    api.grok_BitSet_Invert(this.d);
    return this;
  }

  /** Sets all bits to x
   * @param {boolean} x
   * @param {boolean} notify
   * @returns {BitSet} */
  setAll(x: boolean, notify: boolean = true): BitSet {
    api.grok_BitSet_SetAll(this.d, x, notify);
    return this;
  }

  /** Finds the first index of value x, going forward from i-th position.
   * @param {number} i - index
   * @param {boolean} x
   * @returns {number} */
  findNext(i: number, x: boolean): number {
    return api.grok_BitSet_FindNext(this.d, i, x);
  }

  /** Finds the first index of value x, going forward from i-th position, or -1 if not found.
   * @param {number} i - Index to start searching from.
   * @param {boolean} x - Value to search for.
   * @returns {number} */
  findPrev(i: number, x: boolean): number {
    return api.grok_BitSet_FindPrev(this.d, i, x);
  }

  /** Gets i-th bit
   * @param {number} i
   * @returns {boolean} */
  get(i: number): boolean {
    return api.grok_BitSet_GetBit(this.d, i);
  }

  /** Sets i-th bit to x
   * @param {number} i
   * @param {boolean} x
   * @param {boolean} notify */
  set(i: number, x: boolean, notify: boolean = true): void {
    api.grok_BitSet_SetBit(this.d, i, x, notify);
  }

  /** Sets [i]-th bit to [value], does not check bounds */
  setFast(i: number, value: boolean): void {
    let buf = api.grok_BitSet_GetBuffer(this.d);
    let idx = (i | 0) / 0x20;

    if (value)
      buf[idx] |= 1 << (i & 0x1f);
    else
      buf[idx] &= ~(1 << (i & 0x1f));
  }

  /** Sets all bits by setting i-th bit to the results of f(i)
   * @param {Function} f
   * @returns {BitSet} */
  init(f: IndexPredicate): BitSet {
    let buf = api.grok_BitSet_Get_Buffer(this.d);
    let length = this.length;

    for (let i = 0; i < length; i++)
      buf[i] = 0;

    for (let i = 0; i < length; i++) {
      let idx = (i / 0x20) | 0;
      if (f(i))
        buf[idx] |= 1 << (i & 0x1f);
    }

    api.grok_BitSet_Set_Buffer(this.d, buf);
    this.fireChanged();
    return this;
  }

  /** Indexes of all set bits. The result is cached.
   *  @returns {Int32Array} */
  getSelectedIndexes(): Int32Array {
    return api.grok_BitSet_GetSelectedIndexes(this.d);
  }

  /** Copies the content from the other {BitSet}.
   * @param {BitSet} b - BitSet to copy from.
   * @returns {BitSet} */
  copyFrom(b: BitSet): BitSet {
    api.grok_BitSet_CopyFrom(this.d, b.d);
    return this;
  }

  fireChanged(): void {
    api.grok_BitSet_FireChanged(this.d);
  }

  /** @returns {Observable} - fires when the bitset gets changed. */
  get onChanged(): Observable<any> {
    return observeStream(api.grok_BitSet_Changed(this.d));
  }

  /** Finds the value of similarity between two BitSets.
   * @param {BitSet} b - second BitSet.
   * @param {SimilarityMetric} metric - similarity metric to use.
   * @returns {number} */
  similarityTo(b: BitSet, metric: SimilarityMetric = SIMILARITY_METRIC.TANIMOTO): number {
    return api.grok_BitSet_SimilarityTo(this.d, b.d, metric);
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
  }

  handleClick(rowIndexPredicate: IndexPredicate, mouseEvent: any) {
    api.grok_Utils_SelectRowsWhere(this.d, rowIndexPredicate, mouseEvent.ctrlKey, mouseEvent.shiftKey, mouseEvent.metaKey);
  }
}


/** Represents basic descriptive statistics calculated for a {Column}.
 *  See samples: {@link https://public.datagrok.ai/js/samples/data-frame/stats} */
export class Stats {
  private readonly d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Calculates statistics for the specified column.
   * @param {Column} col
   * @param {BitSet} mask
   * @returns {Stats} */
  static fromColumn(col: Column, mask: BitSet | null = null): Stats {
    return new Stats(api.grok_Stats_FromColumn(col.d, toDart(mask)));
  }

  /** Total number of values (including missing values). */
  get totalCount(): number {
    return api.grok_Stats_Get_TotalCount(this.d);
  }

  /** Number of missing (empty) values. */
  get missingValueCount(): number {
    return api.grok_Stats_Get_MissingValueCount(this.d);
  }

  /** Number of non-empty values. */
  get valueCount(): number {
    return api.grok_Stats_Get_ValueCount(this.d);
  }

  /** @returns {number} - minimum */
  get min(): number {
    return api.grok_Stats_Get_Min(this.d);
  }

  /** @returns {number} - maximum */
  get max(): number {
    return api.grok_Stats_Get_Max(this.d);
  }

  /** @returns {number} - sum */
  get sum(): number {
    return api.grok_Stats_Get_Sum(this.d);
  }

  /** @returns {number} - average */
  get avg(): number {
    return api.grok_Stats_Get_Avg(this.d);
  }

  /** @returns {number} - standard deviation */
  get stdev(): number {
    return api.grok_Stats_Get_Stdev(this.d);
  }

  /** @returns {number} - variance */
  get variance(): number {
    return api.grok_Stats_Get_Variance(this.d);
  }

  /** @returns {number} - skewness */
  get skew(): number {
    return api.grok_Stats_Get_Skew(this.d);
  }

  /** @returns {number} - kurtosis */
  get kurt(): number {
    return api.grok_Stats_Get_Kurt(this.d);
  }

  /** @returns {number} - median value */
  get med(): number {
    return api.grok_Stats_Get_Med(this.d);
  }

  /** @returns {number} - first quartile */
  get q1(): number {
    return api.grok_Stats_Get_Q1(this.d);
  }

  /** @returns {number} - second quartile */
  get q2(): number {
    return api.grok_Stats_Get_Q2(this.d);
  }

  /** @returns {number} - third quartile */
  get q3(): number {
    return api.grok_Stats_Get_Q3(this.d);
  }

  /** Pearson correlation
   * @param {Column} otherColumn
   * @returns {number} */
  corr(otherColumn: Column): number { return api.grok_Stats_Corr(this.d, otherColumn.d); }

  /** Spearman correlation
   * @param {Column} otherColumn
   * @returns {number} */
  spearmanCorr(otherColumn: Column): number { return api.grok_Stats_SpearmanCorr(this.d, otherColumn.d); }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
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
  private readonly d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Performs the aggregation
   *  @returns {DataFrame} */
  aggregate(): DataFrame {
    return new DataFrame(api.grok_GroupByBuilder_Aggregate(this.d));
  }

  /**
   * Performs the aggregation
   * @param {AggregationType} agg - Aggregation type.
   * @param {string} colName - Column name.
   * @param {string} resultColName - Name of the resulting column. Default value is agg(colName).
   * @returns {GroupByBuilder}
   * */
  add(agg: AggregationType, colName?: string | null, resultColName?: string | null): GroupByBuilder {
    api.grok_GroupByBuilder_Add(this.d, agg, colName, resultColName);
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
    return api.grok_GroupByBuilder_GetGroups(this.d);
  }

  /**
   * Specifies the filter for the source rows.
   * @input {String|Object} pattern
   * @returns {GroupByBuilder}
   **/
  where(pattern: string | object): GroupByBuilder {
    api.grok_GroupByBuilder_Where(this.d, pattern);
    return this;
  }

  /**
   * @param {BitSet} bitset
   * @returns {GroupByBuilder} */
  whereRowMask(bitset: BitSet): GroupByBuilder {
    if (bitset != null)
      api.grok_GroupByBuilder_WhereBitSet(this.d, bitset.d);
    return this;
  }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.d);
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
}


//static const int None = -2147483648;
//static const double None = 2.6789344063684636e-34;
//
// class JsColumn {
//
//   /// init
//   void init(Int32Array values);
//
//   /// Removes [count] element at position [idx]. Returns the element at idx-th position.
//   TValue removeAt(int idx, [int count = 1]);
//
//   /// Inserts [count] elements, initialized to the empty value, at position [idx].
//   void insertAt(int idx, [int count = 1]);
//
//   /// Compares elements at positions [idx1] and [idx2].
//   int compare(int idx1, int idx2);
//
//   /// Compares values.
//   int compareValues(TValue v1, TValue v2);
//
//   /// Sets [idx]-th element to None.
//   void setNone(int idx);
//
//   /// Is [idx]-th None?
//   bool isNone(int idx);
//
//   /// Returns whether i-th element exists and is not +-Infinity.
//   /// Applicable for numerical columns only.
//   bool isFinite(int idx) => !isNone(idx);
//
//   /// String representation of the [idx]-th element.
//   /// The result depends on the column formatting options, and therefore might lose precision.
//   String toStr(int idx);
//
//   /// Returns a [double] representation of [idx]-th element.
//   double toDouble(int idx);
//
//   TValue getItem(int pos) => isNone(pos) ? null:  this[pos];
//   void   setItem(int pos, TValue value)
//
//   /// Minimum value of a column.
//   double get min
//
//   /// Maximum value of a column.
//   double get max
// }