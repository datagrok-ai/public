/**
 * Row, Cell, and related classes for data access.
 * @module dataframe/row
 */

import {TYPE, IndexPredicate} from "../const";
import {toDart, toJs} from "../wrappers";
import {_toIterable} from "../utils_convert";
import {DartList} from "../proxies";
import wu, {WuIterable} from "wu";
import {IDartApi} from "../api/grok_api.g";
import type {FilterState} from "../viewer";
import type {DataFrame} from "./data-frame";
import type {Column} from "./column";
import type {BitSet} from "./bit-set";
import type {RowPredicate} from "./types";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/**
 * Represents a row in a DataFrame.
 *
 * This class uses a JavaScript Proxy to allow convenient property access using column names:
 *
 * ```typescript
 * const row = df.row(0);
 * console.log(row.age);      // Access column 'age' value
 * row.name = 'Alice';        // Set column 'name' value
 * ```
 *
 * **Why a Proxy?**
 *
 * The Row constructor returns a Proxy that intercepts property access. When you access
 * `row.columnName`, the Proxy checks if it's a built-in property (`table`, `idx`, `cells`,
 * `get`, `toDart`, `toMap`). If not, it delegates to `table.get(columnName, idx)`.
 *
 * This provides a convenient syntax but has implications:
 * - `instanceof Row` checks may not work as expected (the Proxy is returned, not `this`)
 * - TypeScript cannot provide autocomplete for column names at compile time
 * - There's slight overhead from the Proxy indirection
 *
 * **Alternative access patterns:**
 *
 * For performance-critical code or when you need type safety, use these alternatives:
 *
 * ```typescript
 * // Using get() method
 * const age = row.get('age');
 *
 * // Direct column access (most performant for bulk operations)
 * const ageColumn = df.getCol('age');
 * for (let i = 0; i < df.rowCount; i++) {
 *   const age = ageColumn.get(i);
 * }
 *
 * // Using toMap() for a plain object representation
 * const rowData = row.toMap();  // { age: 25, name: 'Alice', ... }
 * ```
 */
export class Row {
  table: DataFrame;
  readonly idx: number;

  /**
   * Creates a {@link Row} from the specified {@link DataFrame} and row index.
   *
   * **Note:** The constructor returns a Proxy, not the Row instance directly.
   * This allows dynamic property access but means `instanceof Row` checks may fail.
   */
  constructor(table: DataFrame, idx: number) {

    /** @member {DataFrame} */
    this.table = table;
    /** @member {number} */
    this.idx = idx;

    // Return a Proxy to enable dynamic column access via property syntax.
    // This intercepts property get/set and routes them through the DataFrame.
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
        // Pass through known Row properties
        if (name == 'cells' || name == 'get' || name == 'toDart' || name == 'toMap' || target.hasOwnProperty(name))
          return target[<any>name];
        // Otherwise, treat as column name access
        return target.table.get(name, target.idx);
      }
    });
  }

  get cells(): Iterable<Cell> { return _toIterable(api.grok_Row_Get_Cells(this.table.dart, this.idx)); }

  /** Returns a JS object with column names as keys and values as values. */
  toMap(): {[index: string]: any} {
    const res: {[index: string]: any} = {};
    for (const column of this.table.columns)
      res[column.name] = column.get(this.idx);
    return res;
  }

  /** Returns this row's value for the specified column */
  [columnName: string]: any

  /** Returns this row's value for the specified column */
  get(columnName: string): any {
    return this.table.getCol(columnName).get(this.idx);
  }

  toDart(): any { return api.grok_Row(this.table.dart, this.idx); }
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
  static supportedTypes: string[] = [TYPE.FLOAT, TYPE.INT, TYPE.BIG_INT, TYPE.STRING, TYPE.BOOL, TYPE.DATE_TIME];

  private readonly dart: any;

  constructor(dart: any) { this.dart = dart; }

  /** Creates a matcher for the specified column. */
  static forColumn(column: Column, pattern: string): ValueMatcher {
    return new ValueMatcher(api.grok_ValueMatcher_ForColumn(column.dart, pattern));
  }

  /** Creates a matcher for the specified data type. */
  static forType(type: TYPE | string, pattern: string): ValueMatcher {
    switch (type) {
      case TYPE.FLOAT:
      case TYPE.INT:
      case TYPE.BIG_INT:
        return ValueMatcher.numerical(pattern);
      case TYPE.STRING:
        return ValueMatcher.string(pattern);
      case TYPE.BOOL:
        return ValueMatcher.bool(pattern);
      case TYPE.DATE_TIME:
        return ValueMatcher.dateTime(pattern);
      default:
        throw `Value matching not supported for type '${type}'`;
    }
  }

  static numerical(pattern: string): ValueMatcher { return new ValueMatcher(api.grok_ValueMatcher_Numerical(pattern)); }
  static string(pattern: string): ValueMatcher { return new ValueMatcher(api.grok_ValueMatcher_String(pattern)); }
  static dateTime(pattern: string): ValueMatcher { return new ValueMatcher(api.grok_ValueMatcher_DateTime(pattern)); }
  static bool(pattern: string): ValueMatcher { return new ValueMatcher(api.grok_ValueMatcher_Bool(pattern)); }

  /** Expression as entered by user (such as '>42') */
  get pattern(): string { return api.grok_ValueMatcher_Get_Pattern(this.dart); }

  /** Operation (such as '<', 'EQUALS', 'BEFORE', etc). */
  get operator(): string { return api.grok_ValueMatcher_Get_Operator(this.dart); }

  /** Whether [x] passes the filter specified by the [expression].
   * See also {@link validate} for the explanation. */
  match(value: any): boolean { return api.grok_ValueMatcher_Match(this.dart, value); }

  /** Validates the specified conditions. Returns null, if valid, error string otherwise */
  validate(value: any): string | null { return api.grok_ValueMatcher_Validate(this.dart, value); }
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
  readonly dart: any;
  readonly table: DataFrame;

  constructor(table: DataFrame, dart: any) {
    /** @member {DataFrame} */
    this.table = table;
    this.dart = dart;
  }

  /** Gets i-th row. DO NOT USE IN PERFORMANCE-CRITICAL CODE! */
  get(i: number): Row { return new Row(this.table, i); }

  /** List of textual descriptions of currently applied filters */
  get filters(): DartList<string> { return DartList.fromDart(api.grok_RowList_Get_Filters(this.dart)); }

  get mouseOverRowFunc(): IndexPredicate {
    return api.grok_RowList_MouseOverRowFunc(this.dart);
  }

  where(indexPredicate: IndexPredicate): WuIterable<number> {
    return wu(_toIterable(api.grok_RowList_Where(this.dart, indexPredicate)));
  }

  indexes(options?: {onlyFiltered?: boolean, onlySelected?: boolean}): WuIterable<number> {
    return wu(_toIterable(api.grok_RowList_Indexes(this.dart, options?.onlyFiltered ?? false, options?.onlySelected ?? false)));
  }

  /** Removes specified rows
   * @param {number} idx
   * @param {number} [count=1] - Number of rows to remove.
   * @param notify - Whether a change notification should be fired. */
  removeAt(idx: number, count: number = 1, notify: boolean = true): void {
    api.grok_RowList_RemoveAt(this.dart, idx, count, notify);
  }

  /** Removes specified rows
   * @param {RowPredicate} rowPredicate */
  removeWhere(rowPredicate: RowPredicate): void {
    api.grok_RowList_RemoveWhereIdx(this.dart, (i: number) => rowPredicate(this.get(i)));
  }

  /** Removes specified rows
   * @param {IndexPredicate} indexPredicate */
  removeWhereIdx(indexPredicate: IndexPredicate): void {
    api.grok_RowList_RemoveWhereIdx(this.dart, indexPredicate);
  }

  /** Inserts empty rows at the specified position
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

  /** Iterates over all rows. */
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
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/row-matching/patterns} */
  match(query: string | object): RowMatcher {
    return toJs(api.grok_RowList_Match(this.dart, query));
  }

  /**
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

  /** Highlights the corresponding rows. */
  highlight(indexPredicate: IndexPredicate | null): void {
    api.grok_RowList_Highlight(this.dart, indexPredicate);
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
    const {DataFrame} = require('./data-frame');
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
  get value(): any { return toJs(api.grok_Cell_Get_Value(this.dart)); }
  set value(x: any) { api.grok_Cell_Set_Value(this.dart, toDart(x)); }

  /** String representation of the value, if both [column] and [row] are defined;
     otherwise, empty string. */
  get valueString(): string { return api.grok_Cell_Get_ValueString(this.dart); }

  /** Whether the cell is empty */
  isNone(): boolean { return this.column.isNone(this.rowIndex); }

  /** @returns {string} */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }
}
