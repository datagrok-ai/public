/**
 * Collection of columns in a DataFrame.
 * @module dataframe/column-list
 */

import {TYPE, ColumnType, SemType} from "../const";
import {toJs} from "../wrappers";
import {_getIterator, _toIterable} from "../utils_convert";
import wu from "wu";
import {IDartApi} from "../api/grok_api.g";
import type {Column, DateTimeColumn} from "./column";
import type {DataFrame} from "./data-frame";
import type {IndexSetter} from "./types";

declare let DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** Columns in a [DataFrame]. */
export class ColumnList {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  get dataFrame(): DataFrame { return toJs(api.grok_ColumnList_Get_DataFrame(this.dart)); }

  /** Number of columns. */
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

  /** Returns the first column that satisfies the specified criteria. */
  firstWhere(predicate: (col: Column) => boolean): Column | undefined {
    for (const col of this)
      if (predicate(col))
        return col;
  }

  /** Returns all columns. */
  get all(): Iterable<Column> {
    return _toIterable(this.dart);
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

  get dateTime(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_DateTime(this.dart));
  }

  get numericalNoDateTime(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_NumericalNoDateTime(this.dart));
  }

  get boolean(): Iterable<Column> {
    return wu(_toIterable(api.grok_ColumnList_Boolean(this.dart)));
  }

  get selected(): Iterable<Column> {
    return _toIterable(api.grok_ColumnList_Selected(this.dart));
  }

  /** Array containing column names. */
  names(): string[] { return api.grok_ColumnList_Names(this.dart); }

  /** Sets column order.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/columns-ordering}
   *  @param {string[]} columnNames - Order of columns. */
  setOrder(columnNames: string[]): void {
    api.grok_ColumnList_SetOrder(this.dart, columnNames);
  }

  /** Creates an array of columns. */
  toList(): Column[] {
    return this.names().map((name: string) => this.byName(name));
  }

  /** Returns a name->column map. Use it when you need to access columns frequently. */
  toMap(): Map<string, Column> {
    const map = new Map();
    for (const col of this)
      map.set(col.name, col);
    return map;
  }

  /** Adds a column, and optionally notifies the parent dataframe.
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired */
  add(column: Column, notify: boolean = true): Column {
    api.grok_ColumnList_Add(this.dart, column.dart, notify);
    return column;
  }

  /** Returns a column with the specified name and type, or creates a new column if it does not exist.
  * @param {string} name - column name
  * @param {string} type - @see {@link COLUMN_TYPE}
  * @returns {Column} */
  getOrCreate(name: string, type: ColumnType): Column {
    return this.contains(name) ?
      this.byName(name) :
      this.add(DG.Column.fromType(type, name, this.dataFrame.rowCount));
  }

  /** Inserts a column, and optionally notifies the parent dataframe.
   * @param {Column} column - column to insert
   * @param {boolean} notify - whether DataFrame's `changed` event should be fired
   * @param {int} index
   * @returns {Column} */
  insert(column: Column, index: number | null = null, notify: boolean = true): Column {
    api.grok_ColumnList_Insert(this.dart, column.dart, index, notify);
    return column;
  }

  /** Adds an empty column of the specified type. */
  addNew(name: string, type: ColumnType): Column {
    return toJs(api.grok_ColumnList_AddNew(this.dart, name, type));
  }

  /** Adds calculated column.
   * @param {string} name
   * @param {string} expression
   * @param {ColumnType} type
   * @param {bool} treatAsString - if true, [expression] is not evaluated as formula and is treated as a regular string value instead
   * @param {bool} subscribeOnChanges - if true, the column will be recalculated when the source columns change
   * @returns {Column} */
  addNewCalculated(name: string, expression: string, type: ColumnType | 'auto' = 'auto', treatAsString: boolean = false, subscribeOnChanges: boolean = true): Promise<Column> {
    return api.grok_ColumnList_AddNewCalculated(this.dart, name, expression, type, treatAsString, subscribeOnChanges);
  }

  _getNewCalculated(name: string, expression: string, type: ColumnType | 'auto' = 'auto', treatAsString: boolean = false): Promise<Column> {
    return api.grok_ColumnList_GetNewCalculated(this.dart, name, expression, type, treatAsString);
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
  addNewDateTime(name: string): DateTimeColumn { return this.addNew(name, TYPE.DATE_TIME); }

  /** Creates and adds a boolean column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewBool(name: string): Column<boolean> { return this.addNew(name, TYPE.BOOL); }

  /** Creates and adds a byte array column
   * {@link https://dev.datagrok.ai/script/samples/javascript/data-frame/modification/add-columns}
   * */
  addNewBytes(name: string): Column<Uint8Array> { return this.addNew(name, TYPE.BYTE_ARRAY); }

  /** Creates and adds a virtual column.
   * @param {string} name - column name
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
  remove(column: string | number | Column, notify: boolean = true): ColumnList {
    const columnName
      = typeof column === 'string' ? column
      : typeof column === 'number' ? this.byIndex(column).name
      : column.name;

    api.grok_ColumnList_Remove(this.dart, columnName, notify);
    return this;
  }

  /** Checks whether this list contains a column with the specified name. The check is case-insensitive. */
  contains(columnName: string): boolean {
    return api.grok_ColumnList_Contains(this.dart, columnName);
  }

  /** Replaces the column with the new column.
   * @param {boolean} notify */
  replace(columnToReplace: Column | string, newColumn: Column, notify: boolean = true): Column {
    return toJs(api.grok_ColumnList_Replace(this.dart, (typeof columnToReplace === 'string') ? columnToReplace:  columnToReplace.dart, newColumn.dart, notify));
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
