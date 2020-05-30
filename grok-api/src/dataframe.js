import * as rxjs from 'rxjs';
import {AGG, TYPE} from "./const";
import {__obs, observeStream} from "./events";
/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 *
 * In the API, the terms "Table" and "DataFrame" are used interchangeably.
 * See usage samples: https://public.datagrok.ai/js/samples/data-frame/manipulate
 */
export class DataFrame {

    constructor(d) {
        this.d = d;
    }

    /** Creates a {@link DataFrame} with the specified number of rows and no columns.
     * @param {number} rowCount
     * @returns {DataFrame} */
    static create(rowCount = 0) { return new DataFrame(grok_DataFrame(rowCount)); }

    /** Creates a {@link DataFrame} from the specified columns. All columns should be of the same length.
     * @param {Column[]} columns
     * @returns {DataFrame} */
    static fromColumns(columns) { return new DataFrame(grok_DataFrame_FromColumns(columns.map((c) => c.d))); }

    /** Constructs {@link DataFrame} from a comma-separated values string
     * @param {string} csv - The content of the comma-separated values file.
     * @returns {DataFrame} */
    static fromCsv(csv) { return grok.data.parseCsv(csv); }

    /** Constructs {@link DataFrame} from the specified JSON string.
     * @param {string} json - JSON document.
     * @returns {DataFrame} */
    static fromJson(json) { return new DataFrame(grok_DataFrame_FromJson(json)); }

    toString() { return `${this.name} (${this.rowCount} rows, ${this.columns.length} columns)` };

    /** Returns number of rows in the table.
     * @returns {number} */
    get rowCount() { return grok_DataFrame_RowCount(this.d); }

    /** Returns a {@link BitSet} with selected rows.
     * @returns {BitSet} */
    get selection() { return new BitSet(grok_DataFrame_Get_Selection(this.d)); }

    /** Returns a {@link BitSet} with rows that pass filter.
     * @returns {BitSet} */
    get filter() { return new BitSet(grok_DataFrame_Get_Filter(this.d)); }

    /** Name of the dataframe.
     * @returns {string}*/
    get name() { return grok_DataFrame_Get_Name(this.d); }
    set name(s) { grok_DataFrame_Set_Name(this.d, s); }

    /** Returns the value of the specified tag, or null if it does not exist.
     * @returns {string} */
    getTag(tag) { return grok_DataFrame_Get_Tag(this.d, tag); }

    /** Sets a tag to the specified value.
     * @param {string} tag - Key.
     * @param {string} value - Value. */
    setTag(tag, value) { grok_DataFrame_Set_Tag(this.d, tag, value); }

    /** List of columns.
     * @returns {ColumnList} */
    get columns() { return new ColumnList(grok_DataFrame_Columns(this.d)); }

    /** List of columns.
     * @returns {RowList} */
    get rows() { return new RowList(this, grok_DataFrame_Rows(this.d)); }

    /** Returns i-th row.
     * @param {number} i - Row index.
     * @returns {Row} */
    row(i) { return new Row(this, i); }

    /** Returns idx-th value of the specified columns.
     * @param {string} name - Column name.
     * @param {number} idx - Row index. */
    get(name, idx) { return this.getCol(name).get(idx); }

    /** Sets idx-th value of the specified columns.
     * @param {string} name - Column name.
     * @param {number} idx - Row index.
     * @param value - Value. */
    set(name, idx, value) { this.getCol(name).set(idx, value); }

    /** Returns a {@link Column} with the specified name.
     * @param {string} name - Column name.
     * @returns {Column} */
    col(name) { return new Column(grok_DataFrame_ColumnByName(this.d, name)); }

    /** Returns a {@link Cell} with the specified name.
     * @param {number} idx - Row index.
     * @param {string} name - Column name.
     * @returns {Cell} */
    cell(idx, name) { return new Cell(grok_DataFrame_Cell(this.d, idx, name)); }

    /** Same as {@link col}, but throws Error if column is not found
     * @param {string} name - Column name.
     * @returns {Column} */
    getCol(name) {
        let c = this.col(name);
        if (c === null)
            throw new Error(`No such column: ${name}`);
        return c;
    }

    /** Exports the content to comma-separated-values format.
     * @returns {string} */
    toCsv() { return grok_DataFrame_ToCsv(this.d); }

    /** Creates a new dataframe from the specified row mask and a list of columns.
     * @param {BitSet} rowMask - Rows to include.
     * @param {string[]} columnIds - Columns to include.
     * @param {boolean} saveSelection - Whether selection should be saved. */
    clone(rowMask = null, columnIds = null, saveSelection = false) {
        return new DataFrame(grok_DataFrame_Clone(this.d, rowMask.d, columnIds, saveSelection));
    }

    /** Current row
     * @returns {Row} */
    get currentRow() { return new Row(this, grok_DataFrame_Get_CurrentRowIdx(this.d)); }
    set currentRow(idx) { grok_DataFrame_Set_CurrentRowIdx(this.d, idx); }

    /** Current column
     * @returns {Column} */
    get currentCol() { return new Column(grok_DataFrame_Get_CurrentCol(this.d)); }
    set currentCol(col) { grok_DataFrame_Set_CurrentCol(this.d, col.d); }

    /** Current cell
     * @returns {Cell} */
    get currentCell() { return new Cell(grok_DataFrame_Get_CurrentCell(this.d)); }
    set currentCell(cell) { grok_DataFrame_Set_CurrentCell(this.d, cell.d); }

    /** Creates a [DataFrame] from list of objects by using object keys as column names,
     * and object values as values.
     * @param {object[]} list - List of objects. **/
    static fromObjects(list) {
        let table = DataFrame.create(list.length);
        if (list.length === 0)
            return;

        let names = Object.keys(list[0]);
        for (let name of names) {
            let strings = list.map((x) => {
                let value = x[name];
                return value === null ? '' : `${value}`;
            });
            table.columns.add(Column.fromStrings(name, strings));
        }

        return table;
    }

    /** Begins building a query, using the specified columns as keys.
     * @param {string[]} columnNames - Names of the columns to be used as keys.
     * @returns {GroupByBuilder}
     *  */
    groupBy(columnNames = []) { return new GroupByBuilder(grok_DataFrame_GroupBy(this.d, columnNames)); }

    /** @returns {Observable} */ _event(event) { return __obs(event, this.d); }

    /** @returns {Observable} */ get onValuesChanged() { return this._event('ddt-values-changed'); }
    /** @returns {Observable} */ get onCurrentRowChanged() { return this._event('ddt-current-row-changed'); }
    /** @returns {Observable} */ get onMouseOverRowChanged() { return this._event('ddt-mouse-over-row-changed'); }
    /** @returns {Observable} */ get onCurrentColChanged() { return this._event('ddt-current-col-changed'); }
    /** @returns {Observable} */ get onMouseOverColChanged() { return this._event('ddt-mouse-over-col-changed'); }
    /** @returns {Observable} */ get onCurrentCellChanged() { return this._event('ddt-current-cell-changed'); }
    /** @returns {Observable} */ get onMouseOverRowGroupChanged() { return this._event('ddt-mouse-over-row-group-changed'); }
    /** @returns {Observable} */ get onNameChanged() { return this._event('ddt-table-name-changed'); }
    /** @returns {Observable} */ get onMetadataChanged() { return this._event('ddt-table-metadata-changed'); }
    /** @returns {Observable} */ get onColumnNameChanged() { return this._event('ddt-table-column-name-changed'); }
    /** @returns {Observable} */ get onColumnSelectionChanged() { return this._event('ddt-column-selection-changed'); }
    /** @returns {Observable} */ get onColumnsChanged() { return this._event('ddt-columns-changed'); }
    /** @returns {Observable} */ get onColumnsAdded() { return this._event('ddt-columns-added'); }
    /** @returns {Observable} */ get onColumnsRemoved() { return this._event('ddt-columns-removed'); }
    /** @returns {Observable} */ get onRowsAdded() { return this._event('ddt-rows-added'); }
    /** @returns {Observable} */ get onRowsRemoved() { return this._event('ddt-rows-removed'); }

    /** @returns {Observable} */ get onSemanticTypeDetecting() { return this._event('ddt-semantic-type-detecting'); }
    /** @returns {Observable} */ get onSemanticTypeDetected() { return this._event('ddt-semantic-type-detected'); }

    /** @returns {Observable} */ get onDataChanged() { return rxjs.concat(this.onValuesChanged, this.onColumnsAdded,
                                                        this.onColumnsRemoved, this.onRowsAdded, this.onRowsRemoved); }
    /** @returns {Observable} */ get onSelectionChanged() { return this.selection.onChanged; }
    /** @returns {Observable} */ get onFilterChanged() { return this.filter.onChanged; }

    fireValuesChanged() { grok_DataFrame_FireValuesChanged(this.d); }
}

/** Represents a row. Allows for quick property access like "row.height". */
export class Row {

    /** Creates a {@link Row} from the specified {@link DataFrame} and row index.
     * @param {DataFrame} table
     * @param {number} idx
     * @returns {Row} */
    constructor(table, idx) {
        this.table = table;
        this.idx = idx;
        let setables = ['table', 'idx'];
        return new Proxy(this, {
            set(target, name, value) {
                if (setables.includes(name)) {
                    target[name] = value;
                    return true;
                }
                target.table.set(name, target.idx, value);
                return true;
            },
            get(target, name) {
                if (setables.includes(name))
                    return target[name];
                return target.table.get(name, target.idx);
            }
        });
    }

    /** Returns th */
    get(name, idx) { return this.table.getCol(name).get(idx); }
}

/** Strongly-typed column. */
export class Column {
    constructor(d) { this.d = d; }
    static fromStrings(name, list) { return new Column(grok_Column_FromStrings(name, list)); }
    static fromType(type, name = null, length = 0) { return new Column(grok_Column_FromType(type, name, length)); }

    /** [array] will be not be copied and will be used as column's storage */
    static fromInt32Array(name, array, length = null) { return new Column(grok_Column_FromInt32Array(name, array, length)); }

    /** [array] will be not be copied and will be used as column's storage */
    static fromFloat32Array(name, array, length = null) { return new Column(grok_Column_FromFloat32Array(name, array, length)); }

    static fromList(type, name, list) { return new Column(grok_Column_FromList(type, name, list)); }

    /** Creates an integer column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static int(name, length = 0) { return Column.fromType(TYPE.INT, name, length); }

    /** Creates a floating point column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static float(name, length = 0) { return Column.fromType(TYPE.FLOAT, name, length); }

    /** Creates a string column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static string(name, length = 0) { return Column.fromType(TYPE.STRING, name, length); }

    /** Creates a boolean column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static bool(name, length = 0) { return Column.fromType(TYPE.BOOL, name, length); }

    /** Creates a datetime column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static dateTime(name, length = 0) { return Column.fromType(TYPE.DATE_TIME, name, length); }

    /** Column data type. */
    get type() { return grok_Column_Get_Type(this.d); }

    /** Number of elements */
    get length() { return grok_Column_Get_Length(this.d); }

    /** Parent table */
    get dataFrame() { return new DataFrame(grok_Column_Get_DataFrame(this.d)); }

    /** Semantic type */
    get semType() { return grok_Column_Get_SemType(this.d); }
    set semType(s) { grok_Column_Set_SemType(this.d, s); }

    /** Name */
    get name() { return grok_Column_Get_Name(this.d); }
    set name(s) { grok_Column_Set_Name(this.d, s); }

    /** Returns the raw buffer containing data. Return type depends on the column type:
     * {Int32Array} for ints, {@link INT_NULL} represents null.
     * {Float64Array} for floats, {@link FLOAT_NULL} represents null.
     * {Float64Array} for datetime, in microseconds since epoch, {@link DATE_TIME_NULL} represents null.
     * {Int32Array} for strings indexes of {@link categories}.
     * {Uint32Array} bit array.
     * @returns {Array} */
    getRawData() { return grok_Column_GetRawData(this.d); }

    /** Gets i-th value */
    get(i) { return grok_Column_GetValue(this.d, i); }

    /**
     * Sets [i]-th value to [x]
     * @param {number} i
     * @param x
     */
    set(i, x) { grok_Column_SetValue(this.d, i, x); }

    /** Returns whether i-th value is missing.
     * @param {number} i - Row index.
     * @returns {boolean} */
    isNone(i) { return grok_Column_IsNone(this.d, i); }

    getTag(tag) { return grok_Column_Get_Tag(this.d, tag); }
    setTag(tag, value) { grok_Column_Set_Tag(this.d, tag, value); }

    compact() { return grok_Column_Compact(this.d); }

    toList() { return grok_Column_ToList(this.d); }

    /** Returns all unique strings in a sorted order. Applicable to string column only.
     * @returns {string[]} */
    get categories() { return grok_Column_Categories(this.d); }

    get min() { return grok_Column_Min(this.d); }

    get max() { return grok_Column_Max(this.d); }

    get stats() { return Stats.fromColumn(this); }

    *values() {
        for (let i = 0; i < this.length; i++) {
            yield this.get(i);
        }
    }
}

/** Columns in a [DataFrame]. */
export class ColumnList {
    constructor(d) { this.d = d; }

    /** Number of elements in the column.
     * @returns {number} */
    get length() { return grok_ColumnList_Length(this.d); }

    /** Column with the corresponding name (case-insensitive).
     * @param {string} name - Column name
     * @returns {Column} */
    byName(name) { return new Column(grok_ColumnList_ByName(this.d, name)); }

    /** Column by index.
     * @param {number} index - Index of the column.
     * @returns {Column} */
    byIndex(index) { return new Column(grok_ColumnList_ByIndex(this.d, index)); }

    /** First column of [semType], or null. */
    bySemType(semType) { return new Column(grok_ColumnList_BySemType(this.d, semType)); }

    /** Array containing column names.
     * @returns {string[]} */
    names() { return grok_ColumnList_Names(this.d); }

    toList() { return this.names().map(name => this.byName(name)); }
    add(column, notify = false) { grok_ColumnList_Add(this.d, column.d, notify); return column; }
    addNew(name, type) { return new Column(grok_ColumnList_AddNew(this.d, name, type)); }
    remove(name) { grok_ColumnList_Remove(this.d, name); }
    contains(name) { return grok_ColumnList_Contains(this.d, name); }
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
    constructor(table, d) { this.table = table; this.d = d; }

    /** Removes specified rows
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to remove.
     * @param notify - Whether a change notification should be fired. */
    removeAt(idx, count = 1, notify = true) { grok_RowList_RemoveAt(this.d, idx, count, notify); }

    /** Inserts empty rows at the specified position
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to insert.
     * @param notify - Whether a change notification should be fired. */
    insertAt(idx, count = 1, notify = true) { grok_RowList_InsertAt(this.d, idx, count, notify); }

    /** Appends a new row with the specified values
     * @param values - List of values (length and types should match columns)
     * @param notify - Whether a change notification should be fired.
     * @returns {Row} */
    addNew(values = null, notify = true) { return new Row(this.table, grok_RowList_AddNew(this.d, values, notify));}

    /** Iterates over all rows. */
    *iterator() {
        for (let i = 0; i < this.table.rowCount; i++)
            yield new Row(this.table, i);
    }

    /** Sets values for the specified row.
     * @param {number} idx - Row index.
     * @param values - List of values (length and types should match columns) */
    setValues(idx, values) { grok_RowList_SetValues(this.d, idx, values); }
}

/** Represents a table cell. */
export class Cell {
    constructor(d) { this.d = d; }

    get dataFrame() { return new DataFrame(grok_Cell_Get_DataFrame(this.d)); }
    get row() { return new Row(this.dataFrame, this.rowIndex); }
    get rowIndex() { return grok_Cell_Get_RowIndex(this.d); }
    get column() { return new Column(grok_Cell_Get_Column(this.d)); }
    get value() { return grok_Cell_Get_Value(this.d); }
}

/** Represents basic descriptive statistics calculated for a {Column}. */
export class Stats {
    constructor(d) { this.d = d; }

    static fromColumn(col) { return new Stats(grok_Stats_FromColumn(col.d)); }

    /** Total number of values (including missing values). */
    get totalCount() { return grok_Stats_Get_TotalCount(this.d); }

    /** Number of missing (empty) values. */
    get missingValueCount() { return grok_Stats_Get_MissingValueCount(this.d); }

    /** Number of non-empty values. */
    get valueCount() { return grok_Stats_Get_ValueCount(this.d); }

    get min() { return grok_Stats_Get_Min(this.d); }
    get max() { return grok_Stats_Get_Max(this.d); }
    get sum() { return grok_Stats_Get_Sum(this.d); }
    get avg() { return grok_Stats_Get_Avg(this.d); }
    get stdev() { return grok_Stats_Get_Stdev(this.d); }
    get variance() { return grok_Stats_Get_Variance(this.d); }
    get skew() { return grok_Stats_Get_Skew(this.d); }
    get kurt() { return grok_Stats_Get_Kurt(this.d); }
    get med() { return grok_Stats_Get_Med(this.d); }
    get q1() { return grok_Stats_Get_Q1(this.d); }
    get q2() { return grok_Stats_Get_Q2(this.d); }
    get q3() { return grok_Stats_Get_Q3(this.d); }
}

/** Fluid API for building an aggregation query against a {@link DataFrame}. */
export class GroupByBuilder {
    constructor(d) { this.d = d; }

    /**
     * Performs the aggregation
     * @returns {DataFrame}
     * */
    aggregate() { return new DataFrame(grok_GroupByBuilder_Aggregate(this.d)); }

    /**
     * Performs the aggregation
     * @param {AggregationType} agg - Aggregation type.
     * @param {string} colName - Column name.
     * @param {string=} resultColName - Name of the resulting column. Default value is agg(colName).
     * @returns {DataFrame}
     * */
    add(agg, colName, resultColName = null) { grok_GroupByBuilder_Add(this.d, agg, colName, resultColName); return this; }

    key(col, resultColName = null) { this.add(AGG.KEY, col, resultColName); return this; }
    pivot(col, resultColName = null) { this.add(AGG.PIVOT, col, resultColName); return this; }
    count(col, resultColName = null) { this.add(AGG.TOTAL_COUNT, col, resultColName); return this; }
    uniqueCount(col, resultColName = null) { this.add(AGG.UNIQUE_COUNT, col, resultColName); return this; }
    missingValueCount(col, resultColName = null) { this.add(AGG.MISSING_VALUE_COUNT, col, resultColName); return this; }
    valueCount(col, resultColName = null) { this.add(AGG.VALUE_COUNT, col, resultColName); return this; }
    min(col, resultColName = null) { this.add(AGG.MIN, col, resultColName); return this; }
    max(col, resultColName = null) { this.add(AGG.MAX, col, resultColName); return this; }
    sum(col, resultColName = null) { this.add(AGG.SUM, col, resultColName); return this; }
    med(col, resultColName = null) { this.add(AGG.MED, col, resultColName); return this; }
    avg(col, resultColName = null) { this.add(AGG.AVG, col, resultColName); return this; }
    stdev(col, resultColName = null) { this.add(AGG.STDEV, col, resultColName); return this; }
    variance(col, resultColName = null) { this.add(AGG.VARIANCE, col, resultColName); return this; }
    q1(col, resultColName = null) { this.add(AGG.Q1, col, resultColName); return this; }
    q2(col, resultColName = null) { this.add(AGG.Q2, col, resultColName); return this; }
    q3(col, resultColName = null) { this.add(AGG.Q3, col, resultColName); return this; }
}

/** Efficient bit storage and manipulation. */
export class BitSet {
    /** Creates a {BitSet} from the specified Dart object. */
    constructor(d) { this.d = d; }

    /** Creates a {BitSet} of the specified length with all bits set to false.
     * @param {number} length - Number of bits.
     * @returns {BitSet} */
    static create(length) { return new BitSet(grok_BitSet(length)); }

    /** Number of bits in a bitset
     * @returns {number} */
    get length() { return grok_BitSet_Get_Length(this.d); }

    /** Number of set bits
     * @returns {number} */
    get trueCount() { return grok_BitSet_Get_TrueCount(this.d); }

    /** Number of unset bits
     * @returns {number}*/
    get falseCount() { return grok_BitSet_Get_FalseCount(this.d); }

    /** Clones a bitset
     *  @returns {BitSet} */
    clone() { return new BitSet(grok_BitSet_Clone(this.d)); }

    /** Inverts a bitset */
    invert() { grok_BitSet_Invert(this.d); }

    /** Sets all bits to x
     * @param {boolean} x */
    setAll(x) { grok_BitSet_SetAll(this.d, x); }

    /** Finds the first index of value x, going forward from i-th position */
    findNext(i, x) { return grok_BitSet_FindNext(this.d, i, x); }

    /** Finds the first index of value x, going forward from i-th position, or -1 if not found.
     * @param {number} i - Index to start searching from.
     * @param {boolean} x - Value to search for.
     * @returns {number}
     * */
    findPrev(i, x) { return grok_BitSet_FindPrev(this.d, i, x); }

    /** Gets i-th bit
     * @param {number} i*/
    get(i) { return grok_BitSet_GetBit(this.d, i); }

    /** Sets i-th bit to x
     * @param {number} i
     * @param {boolean} x */
    set(i, x) { grok_BitSet_SetBit(this.d, i, x); }

    /** Indexes of all set bits. The result is cached.
     *  @returns {Int32Array} */
    getSelectedIndexes() { return grok_BitSet_GetSelectedIndexes(this.d); }

    /** Copies the content from the other {BitSet}.
     * @param {BitSet} b - BitSet to copy from.
     * */
    copyFrom(b) { grok_BitSet_CopyFrom(this.d, b.d); }

    /** @returns {Observable} - fires when the bitset gets changed. */
    get onChanged() { return observeStream(grok_BitSet_Changed(this.d)); }
}
