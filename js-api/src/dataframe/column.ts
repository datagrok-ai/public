/**
 * Column classes for strongly-typed data storage.
 * @module dataframe/column
 */

import {TYPE, FLOAT_NULL, ColumnType, SemType, ColumnAggregationType, TAGS} from "../const";
import {toDart, toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import dayjs from "dayjs";
import {IDartApi} from "../api/grok_api.g";
import type {Comparer} from "./types";
import type {BitSet} from "./bit-set";
import type {Stats} from "./stats";
import type {DataFrame} from "./data-frame";
import {Qnum, QNUM_EXACT} from "./qnum";
import type {ColumnMetaHelper} from "./column-helpers";

declare let DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** Strongly-typed column.
 * Use {@link get} and {@link set} to access elements by index.
 * */
export class Column<T = any, TInit = T> {
  public dart: any;
  public temp: any;
  public tags: any;
  private _meta: ColumnMetaHelper | undefined;

  constructor(dart: any) {
    this.dart = dart;
    this.temp = new MapProxy(api.grok_Column_Get_Temp(this.dart));
    this.tags = new MapProxy(api.grok_Column_Get_Tags(this.dart));
  }

  get meta(): ColumnMetaHelper {
    if (this._meta == undefined) {
      const {ColumnMetaHelper: CMH} = require('./column-helpers');
      this._meta = new CMH(this);
    }
    return this._meta!;
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
   * @see COLUMN_TYPE
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

  /** [array] will be not be copied and will be used as column's storage */
  static fromFloat64Array(name: string, array: Float64Array, length: number | null = null): Column<number> {
    return toJs(api.grok_Column_FromFloat64Array(name, array, length));
  }

  /** Creates BigIntColumn from BigInt64Array / BigUint64Array */
  static fromBigInt64Array(name: string, array: BigInt64Array | BigUint64Array) {
    return this.fromList(TYPE.BIG_INT as ColumnType, name, Array.from(array, (v: any, _) =>
      api.grok_BigIntJs_To_BigInt(v.toString())));
  }

  /**
   * Creates a {@link Column} from the list of values.
   * @param {string} type - @see {@link COLUMN_TYPE}
   * @param {string} name - column name
   * @param {object[]} list - list of values
   * @returns {Column} */
  static fromList(type: ColumnType, name: string, list: any[]): Column {
    if (type === TYPE.DATE_TIME)
      list = list.map((v) => v?.valueOf());
    return toJs(api.grok_Column_FromList(type, name, list));
  }

  /**
   * Crates a {@link Column} of string type from categories and indexes
   * @param name
   * @param categories
   * @param indexes
   */
  static fromIndexes(name: string, categories: string[], indexes: Int32Array): Column {
    return toJs(api.grok_Column_FromIndexes(name, categories, indexes));
  }

  /** Creates a {@link Column} from the bitset.
   * @param {string} name - column name
   * @param {BitSet} bitset - bitset. The resulting boolean column will be of length bitset.length
   * @returns {Column} */
  static fromBitSet(name: string, bitset: BitSet): Column<boolean> {
    return toJs(api.grok_Column_FromBitSet(name, bitset.dart));
  }

  /** Creates an integer column with the specified name and length. */
  static int(name: string, length: number = 0): Column<number> {
    return Column.fromType(TYPE.INT, name, length);
  }

  /** Creates a floating point column with the specified name and length. */
  static float(name: string, length: number = 0): Column<number> {
    return Column.fromType(TYPE.FLOAT, name, length);
  }

  /** Creates a string column with the specified name and length. */
  static string(name: string, length: number = 0): Column<string> {
    return Column.fromType(TYPE.STRING, name, length);
  }

  /** Creates a boolean column with the specified name and length. */
  static bool(name: string, length: number = 0): Column<boolean> {
    return Column.fromType(TYPE.BOOL, name, length);
  }

  /** Creates a datetime column with the specified name and length. */
  static dateTime(name: string, length: number = 0): Column {
    return Column.fromType(TYPE.DATE_TIME, name, length);
  }

  /** Creates a column containing dataframes with the specified name and length. */
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

  /** Is the column numerical (float, int, bigint, qnum)
  * @type {boolean}*/
  get isNumerical(): boolean {
    return api.grok_Column_Get_Is_Numerical(this.dart);
  }

  /** Is the column categorical (string, boolean)
  * @type {boolean}*/
  get isCategorical(): boolean {
    return api.grok_Column_Get_Is_Categorical(this.dart);
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

  // Obsolete. Recommended method is "meta.dialogs".
  get dialogs(): any {
    return this.meta.dialogs;
  }

  // Obsolete. Recommended method is "meta.colors".
  get colors(): any {
    return this.meta.colors;
  }

  // Obsolete. Recommended method is "meta.markers".
  get markers(): any {
    return this.meta.markers;
  }

  /**
   * Initializes all values in the column to [columnInitializer].
   * @param {string | number | boolean | Function} valueInitializer value, or a function that returns value by index
   * */
  init(valueInitializer: string | number | boolean | Date | dayjs.Dayjs | null | ((ind: number) => any)): Column {
    let initType = typeof valueInitializer;
    if (initType === 'function' && this.type === DG.TYPE.DATA_FRAME) {
      // @ts-ignore
      api.grok_Column_Init(this.dart, (i) => toDart(valueInitializer(i)));
    }
    else if (initType === 'function')
      api.grok_Column_Init(this.dart, valueInitializer);

    else if (initType === 'number' || initType === 'string' || initType === 'boolean')
      api.grok_Column_SetAllValues(this.dart, valueInitializer);
    return this;
  }

  /** Performs deep cloning, optionally taking mask of the rows to be included.
   * Note that the cloned colum is not added to this column's dataframe. */
  clone(mask?: BitSet): Column<T> { return new Column(api.grok_Column_Clone(this.dart, toDart(mask))); }

  /** FOR EXPERT USE ONLY!
   *
   * Returns the raw buffer containing data.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/performance/access}
   * Return type depends on the column type:
   * {Int32Array} for ints, {@link INT_NULL} represents null.
   * {Float32Array} or {Float64Array} for floats depending on doublePrecision parameter, {@link FLOAT_NULL} represents null.
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
   *  An empty string is returned if there is no value. */
  getString(i: number): string {
    return api.grok_Column_GetString(this.dart, i);
  }

  /** Returns i-th value as number. None if isNone(i) */
  getNumber(i: number): number {
    return api.grok_Column_GetNumber(this.dart, i);
  }

  /** Attempts to set i-th value by converting a provided string to the corresponding strongly-typed value.
   *  Returns true if text was successfully parsed and set, otherwise false.
   *  Examples: dateColumn.setString('April 1, 2020');
   *            intColumn.setString('42');
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired */
  setString(i: number, str: string, notify: boolean = true): boolean {
    return api.grok_Column_SetString(this.dart, i, str, notify);
  }

  /** Call this method after setting elements without notifications via {@link set} */
  fireValuesChanged() {
    api.grok_Column_FireValuesChanged(this.dart);
  }

  /**
   * Sets [i]-th value to [x], and optionally notifies the dataframe about this change.
   * @param {number} i - Row index.
   * @param {TInit | null} value - Value to set.
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired. Call {@link fireValuesChanged}
   * after you are done modifying the column.
   */
  set(i: number, value: TInit | null, notify: boolean = true): void {
    api.grok_Column_SetValue(this.dart, i, toDart(value), notify);
  }

  /** Returns whether all values in the columns are empty. */
  get isEmpty(): boolean {
    return this.stats.missingValueCount == this.length;
  }

  /** Returns whether i-th value is missing. */
  isNone(i: number): boolean {
    return api.grok_Column_IsNone(this.dart, i);
  }

  /** Returns the value of the specified tag, or null if the tag is not present. */
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
  get stats(): Stats {
    const {Stats} = require('./stats');
    return Stats.fromColumn(this);
  }

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
    return api.grok_Column_ApplyFormula(this.dart, formula, type, treatAsString);
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

  /** @returns {Float32Array} - typed array of float values representing the column.
   * Does not guarantee to perform a copy of the underlying data. */
  asDoubleList(): Float32Array {
    return api.grok_Column_AsDoubleList(this.dart);
  }
}


export class FloatColumn extends Column<number> {
  get doublePrecision(): boolean {
    return api.grok_FloatColumn_GetDoublePrecision(this.dart);
  }
  set doublePrecision(enableDoublePrecision: boolean) {
    api.grok_FloatColumn_SetDoublePrecision(this.dart, enableDoublePrecision);
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

    // @ts-ignore: fallback for the browsers that don't support BigInt, such as Dartium
    return BigInt ? BigInt(v) : parseInt(v);
  }

  /**
   * Sets [i]-th value to [x], and optionally notifies the dataframe about this change.
   */
  set(i: number, value: BigInt | null, notify: boolean = true): void {
    api.grok_BigIntColumn_SetValue(this.dart, i, value?.toString(), notify);
  }
}


type DateTimeInit = dayjs.Dayjs | string | Date | null;

export class DateTimeColumn extends Column<dayjs.Dayjs, DateTimeInit> {

  static getMs(x: any): number | null {
    if (x == '' || x == null)
      return null;
    if (dayjs.isDayjs(x))
      return x.valueOf();
    else {
      const ms = dayjs(x).valueOf();
      if (isNaN(ms))
        throw `"${x}" is not convertable to date time`;
      return ms;
    }
  }

  init(valueInitializer: DateTimeInit | ((ind: number) => DateTimeInit)): Column {
    let initType = typeof valueInitializer;
    const length = this.length;

    if (initType === 'function')
      for (let i = 0; i < length; i++)
        // @ts-ignore
        this.set(i, valueInitializer(i), false);
    else if (initType === 'number' || initType === 'string' || dayjs.isDayjs(valueInitializer)) {
      const ms = DateTimeColumn.getMs(valueInitializer);
      for (let i = 0; i < length; i++)
        api.grok_DateTimeColumn_SetValue(this.dart, i, ms, false);
    }

    this.fireValuesChanged();
    return this;
  }

  get(row: number): dayjs.Dayjs | null {
    let v = api.grok_DateTimeColumn_GetValue(this.dart, row);
    if (v == null)
      return null;
    return dayjs(v);
  }

  set(i: number, value: dayjs.Dayjs | string | Date | null, notify: boolean = true): void {
    // @ts-ignore
    api.grok_DateTimeColumn_SetValue(this.dart, i, DateTimeColumn.getMs(value)?.valueOf(), notify);
  }
}


export class ObjectColumn extends Column<any> {
  /**
   * Gets [i]-th value.
   */
  get(row: number): any | null {
    return toJs(api.grok_Column_GetValue(this.dart, row));
  }
}

/** Column of type [DataFrame]. */
export class DataFrameColumn extends Column<DataFrame> {
  /**  Gets [i]-th value. */
  get(row: number): DataFrame | null {
    return toJs(api.grok_Column_GetValue(this.dart, row));
  }

  /** Returns all values as an array. */
  toList(): Array<DataFrame> {
    return api.grok_Column_ToList(this.dart).map((x: any) => toJs(x));
  }
}
