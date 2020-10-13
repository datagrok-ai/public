import * as rxjs from 'rxjs';
import {AGG, TYPE, COLUMN_TYPE} from "./const";
import {__obs, observeStream} from "./events";
import {toDart, toJs} from "./wrappers";

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
        this.columns = new ColumnList(grok_DataFrame_Columns(this.d));
        this.rows = new RowList(this, grok_DataFrame_Rows(this.d));
        this.filter = new BitSet(grok_DataFrame_Get_Filter(this.d));
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
     * @type {Row} */
    get currentRow() { return new Row(this, grok_DataFrame_Get_CurrentRowIdx(this.d)); }
    set currentRow(idx) { grok_DataFrame_Set_CurrentRowIdx(this.d, idx); }

    /** Current column
     * @type {Column} */
    get currentCol() { return new Column(grok_DataFrame_Get_CurrentCol(this.d)); }
    set currentCol(col) { grok_DataFrame_Set_CurrentCol(this.d, col.d); }

    /** Current cell
     * @type {Cell} */
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

    /** Converts a column with the specified name to [newType],
     * removes the original column from its dataframe and adds the new column to it.
     * @param {string|Column} column
     * @param {string} newType
     * @param {string=} format
     * @returns {Column} */
    changeColumnType(column, newType, format = null) {
        return new Column(grok_DataFrame_ChangeColumnType(this.d, toDart(column), newType, format));
    }

    /** Begins building a query, using the specified columns as keys.
     * @param {string[]} columnNames - Names of the columns to be used as keys.
     * @returns {GroupByBuilder}
     *  */
    groupBy(columnNames = []) { return new GroupByBuilder(grok_DataFrame_GroupBy(this.d, columnNames)); }

    append(t2, inPlace = false) { return new DataFrame(grok_DataFrame_Append(this.d, t2.d, inPlace)); }

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
    /** @returns {Observable} */ get onRowsFiltering() { return this._event('ddt-rows-filtering'); }

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

        /** @member {DataFrame} */
        this.table = table;
        /** @member {number} */
        this.idx = idx;

        const setables = ['table', 'idx'];
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
    static fromType(type, name = null, length = 0) {
        return new Column(grok_Column_FromType(type, name, length));
    }

    /** [array] will be not be copied and will be used as column's storage */
    static fromInt32Array(name, array, length = null) { return new Column(grok_Column_FromInt32Array(name, array, length)); }

    /** [array] will be not be copied and will be used as column's storage */
    static fromFloat32Array(name, array, length = null) { return new Column(grok_Column_FromFloat32Array(name, array, length)); }

    /** Creates a {@link Column} from the list of values.
     * @param {string} type
     * @param {string} name
     * @param {object[]} list
     * @returns {Column} */
    static fromList(type, name, list) { return new Column(grok_Column_FromList(type, name, list)); }

    /** Creates an integer column with the specified name and length.
     * @param {string} name
     * @param {number} length
     * @returns {Column} */
    static int(name, length = 0) { return Column.fromType(TYPE.INT, name, length); }

    /** Creates a floating point column with the specified name and length.
     * @param {string} name
     * @param {number} length
     * @returns {Column} */
    static float(name, length = 0) { return Column.fromType(TYPE.FLOAT, name, length); }

    /** Creates a string column with the specified name and length.
     * @param {string} name
     * @param {number} length
     * @returns {Column} */
    static string(name, length = 0) { return Column.fromType(TYPE.STRING, name, length); }

    /** Creates a boolean column with the specified name and length.
     * @param {string} name
     * @param {number} length
     * @returns {Column} */
    static bool(name, length = 0) { return Column.fromType(TYPE.BOOL, name, length); }

    /** Creates a datetime column with the specified name and length.
     * @param {string} name
     * @param {number} length
     * @returns {Column} */
    static dateTime(name, length = 0) { return Column.fromType(TYPE.DATE_TIME, name, length); }

    /** Creates a qualified number column with the specified name and length.
     *  Initialized values with [values], if it is specified; strips out the qualifier
     *  part if [exact] is true.
     *
     * @param {string} name
     * @param {number} length
     * @param {number[]} values
     * @param {boolean} exact - if true, strips out qualifier from [values].
     * */
    static qnum(name, length = 0, values = null, exact = true) {
        let col = Column.fromType(TYPE.QNUM, name, length);
        if (values !== null) {
            let buffer = col.getRawData();
            for (let i = 0; i < length; i++)
                buffer[i] = exact ? Qnum.exact(values[i]) : values[i];
            col.setRawData(buffer);
        }
        return col;
    }

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
     * {Float32Array} for floats, {@link FLOAT_NULL} represents null.
     * {Float64Array} for qnums, {@link FLOAT_NULL} represents null.
     * {Float64Array} for datetime, in microseconds since epoch, {@link FLOAT_NULL} represents null.
     * {Int32Array} for strings indexes of {@link categories}.
     * {Uint32Array} bit array.
     * @returns {Array} */
    getRawData() { return grok_Column_GetRawData(this.d); }

    setRawData(rawData, notify = true) { grok_Column_SetRawData(this.d, rawData, notify); }

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

    /** Gets the value of the specified tag.
     *  @param {string} tag
     *  @returns {string} */
    getTag(tag) { return grok_Column_Get_Tag(this.d, tag); }

    /** Sets a tag to the specified value.
     * @param {string} tag - Key.
     * @param {string} value - Value. */
    setTag(tag, value) { grok_Column_Set_Tag(this.d, tag, value); }

    /** Compacts the internal column representation.
     *  Currently, it only affects string columns where values were modified. */
    compact() { return grok_Column_Compact(this.d); }

    /** Copies column values to an array.
     *  @ returns {Array} */
    toList() { return grok_Column_ToList(this.d); }

    /** Returns all unique strings in a sorted order. Applicable to string column only.
     * @returns {string[]} */
    get categories() { return grok_Column_Categories(this.d); }

    /** Column's minimum value. The result is cached.
     * @returns {number} */
    get min() { return grok_Column_Min(this.d); }

    /** Column's maximum value. The result is cached.
     * @returns {number} */
    get max() { return grok_Column_Max(this.d); }

    /** Basic descriptive statistics. The result is cached.
     * @returns {Stats} */
    get stats() { return Stats.fromColumn(this); }

    /** An iterator over all values in this column. */
    *values() {
        for (let i = 0; i < this.length; i++) {
            yield this.get(i);
        }
    }

    /** Creates and returns a new column by converting [column] to the specified [newType].
     *  @param {string} newType
     *  @param {string} format
     *  @returns {Column} */
    convertTo(newType, format = null) {
        return new Column(grok_Column_ConvertTo(this.d, newType, format));
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

    /** First column of [semType], or null.
     * @returns {Column} */
    bySemType(semType) {
        var col = grok_ColumnList_BySemType(this.d, semType);
        return col == null ? null : new Column(col);
    }

    /** Finds columns by the corresponding semTypes, or null, if any of the sem types could not be found.
     * @returns {Column[]} */
    bySemTypesExact(semTypes) {
        let columns = [];
        for (let semType of semTypes) {
            let col = this.bySemType(semType);
            if (col == null)
                return null;
            columns.push(col);
        }
        return columns;
    }

    //todo
    //numerical

    /** @returns {Column[]} */
    get categorical() {
        //todo: convert to iterable
        let result = [];
        for (let i = 0; i < this.length; i++)
            if (this.byIndex(i).type === COLUMN_TYPE.STRING)
                result.push(this.byIndex(i));
        return result;
    }

    /** Array containing column names.
     * @returns {string[]} */
    names() { return grok_ColumnList_Names(this.d); }

    /** Creates an array of columns.
     * @returns {Column[]} */
    toList() { return this.names().map(name => this.byName(name)); }

    /** Adds a column, and optionally notifies the parent dataframe.
     * @param {Column} column
     * @param {boolean} notify
     * @returns {Column} */
    add(column, notify = true) {
        grok_ColumnList_Add(this.d, column.d, notify);
        return column;
    }

    /** Adds an empty column of the specified type.
     * @param {string} name
     * @param {ColumnType} type
     * @returns {Column} */
    addNew(name, type) { return new Column(grok_ColumnList_AddNew(this.d, name, type)); }

    /** Adds a virtual column.
     * @param {string} name
     * @param {Function} getValue - value constructor function that accepts int index and returns value
     * @returns {Column} */
    addNewVirtual(name, getValue) { return toJs(grok_ColumnList_AddNewVirtual(this.d, name, getValue)); }

    /** Removes column by name (case-insensitive).
     * @param {string} columnName
     * @returns {ColumnList} */
    remove(columnName) { grok_ColumnList_Remove(this.d, columnName); return this; }

    /** Checks wheter this list contains a column with the specified name. The check is case-insensitive.
     * @returns {boolean} */
    contains(columnName) { return grok_ColumnList_Contains(this.d, columnName); }
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
    constructor(table, d) {
        /** @member {DataFrame} */
        this.table = table; this.d = d;
    }

    /** Removes specified rows
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to remove.
     * @param notify - Whether a change notification should be fired. */
    removeAt(idx, count = 1, notify = true) {
        grok_RowList_RemoveAt(this.d, idx, count, notify);
    }

    /** Inserts empty rows at the specified position
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to insert.
     * @param notify - Whether a change notification should be fired. */
    insertAt(idx, count = 1, notify = true) {
        grok_RowList_InsertAt(this.d, idx, count, notify);
    }

    /** Appends a new row with the specified values
     * @param values - List of values (length and types should match columns)
     * @param notify - Whether a change notification should be fired.
     * @returns {Row} */
    addNew(values = null, notify = true) {
        return new Row(this.table, grok_RowList_AddNew(this.d, values, notify));
    }

    /** Iterates over all rows.
     * @returns {Iterable.<Row>}
     * */
    *[Symbol.iterator]() {
        for (let i = 0; i < this.table.rowCount; i++)
            yield new Row(this.table, i);
    }

    /** Sets values for the specified row.
     * @param {number} idx - Row index.
     * @param values - List of values (length and types should match columns) */
    setValues(idx, values) {
        grok_RowList_SetValues(this.d, idx, values);
    }

    _applyPredicate(bitset, rowPredicate) {
        for (let row of this) {
            bitset.set(row.idx, rowPredicate(row));
        }
    }

    select(rowPredicate) {
        _applyPredicate(this.table.selection, rowPredicate);
    }

    filter(rowPredicate) {
        _applyPredicate(this.table.filter, rowPredicate);
    }

    /** Viewers that filter rows should subscribe to DataFrame.onRowsFiltering event.
     * When filtering conditions are changed, viewers should call requestFilter(). */
    requestFilter() { grok_RowList_RequestFilter(this.d); }
}

/** Represents a table cell. */
export class Cell {
    constructor(d) { this.d = d; }

    /** Corresponding table.
     * @returns {DataFrame} */
    get dataFrame() { return new DataFrame(grok_Cell_Get_DataFrame(this.d)); }

    /** Corresponding row.
     * @returns {Row} */
    get row() { return new Row(this.dataFrame, this.rowIndex); }

    /** Index of the corresponding row.
     * @returns {number} */
    get rowIndex() { return grok_Cell_Get_RowIndex(this.d); }

    /** Corresponding column.
     * @returns {Column} */
    get column() { return new Column(grok_Cell_Get_Column(this.d)); }

    /** Cell value.
     * @returns {*} */
    get value() { return grok_Cell_Get_Value(this.d); }
}

/**
 * Efficient bit storage and manipulation.
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 */
export class BitSet {

    /** Creates a {BitSet} from the specified Dart object. */
    constructor(d) { this.d = d; }

    /** Creates a {BitSet} from the string representing the bitset.
     * @param {string} zerosOnes - A string containing '1' and '0'.
     * @returns {BitSet} */
    static fromString(zerosOnes) { return new BitSet(grok_BitSet_FromString(zerosOnes)); }

    /** Creates a {BitSet} from the ArrayBuffer representing the bitset.
     * @param {ArrayBuffer} buffer - An array containing 1 and 0.
     * @param {Number} bitLength - count of bits.
     * @returns {BitSet} */
    static fromBytes(buffer, bitLength) {
        if (bitLength == null || !Number.isInteger(bitLength) || bitLength < 0)
            bitLength = buffer.byteLength * 8;
        return new BitSet(grok_BitSet_FromBytes(buffer, bitLength));
    }

    /** Creates a {BitSet} of the specified length with all bits set to false.
     * @param {number} length - Number of bits.
     * @returns {BitSet} */
    static create(length) { return new BitSet(grok_BitSet(length)); }

    toBinaryString() { return grok_BitSet_ToBinaryString(this.d); }

    /** Number of bits in a bitset
     * @type {number} */
    get length() { return grok_BitSet_Get_Length(this.d); }

    /** Number of set bits
     * @type {number} */
    get trueCount() { return grok_BitSet_Get_TrueCount(this.d); }

    /** Number of unset bits
     * @type {number}*/
    get falseCount() { return grok_BitSet_Get_FalseCount(this.d); }

    /** Clones a bitset
     *  @returns {BitSet} */
    clone() { return new BitSet(grok_BitSet_Clone(this.d)); }

    /** Inverts a bitset.
     * @returns {BitSet} */
    invert() { grok_BitSet_Invert(this.d); return this; }

    /** Sets all bits to x
     * @param {boolean} x
     * @param {boolean} notify
     * @returns {BitSet} */
    setAll(x, notify = true) { grok_BitSet_SetAll(this.d, x, notify); return this; }

    /** Finds the first index of value x, going forward from i-th position.
     * @param {number} i - index
     * @param {boolean} x
     * @returns {number} */
    findNext(i, x) { return grok_BitSet_FindNext(this.d, i, x); }

    /** Finds the first index of value x, going forward from i-th position, or -1 if not found.
     * @param {number} i - Index to start searching from.
     * @param {boolean} x - Value to search for.
     * @returns {number} */
    findPrev(i, x) { return grok_BitSet_FindPrev(this.d, i, x); }

    /** Gets i-th bit
     * @param {number} i
     * @returns {boolean} */
    get(i) { return grok_BitSet_GetBit(this.d, i); }

    /** Sets i-th bit to x
     * @param {number} i
     * @param {boolean} x
     * @param {boolean} notify */
    set(i, x, notify = true) { grok_BitSet_SetBit(this.d, i, x, notify); }

    /** Sets [i]-th bit to [value], does not check bounds */
    setFast(i, value) {
        let buf = grok_BitSet_GetBuffer(this.d);
        let idx = (i | 0) / 0x20;

        if (value)
            buf[idx] |= 1 << (i & 0x1f);
        else
            buf[idx] &= ~(1 << (i & 0x1f));
    }

    /** Sets all bits by setting i-th bit to the results of f(i)
     * @param {Function} f  */
    init(f) {
        let buf = grok_BitSet_Get_Buffer(this.d);
        let length = this.length;

        for (let i = 0; i < length; i++)
            buf[i] = 0;

        for (let i = 0; i < length; i++) {
            let idx = (i / 0x20) | 0;
            if (f(i))
                buf[idx] |= 1 << (i & 0x1f);
        }

        grok_BitSet_Set_Buffer(this.d, buf);
        this.fireChanged();
    }

    /** Indexes of all set bits. The result is cached.
     *  @returns {Int32Array} */
    getSelectedIndexes() { return grok_BitSet_GetSelectedIndexes(this.d); }

    /** Copies the content from the other {BitSet}.
     * @param {BitSet} b - BitSet to copy from.
     * @returns {BitSet} */
    copyFrom(b) { grok_BitSet_CopyFrom(this.d, b.d); return this; }

    fireChanged() { grok_BitSet_FireChanged(this.d); }

    /** @returns {Observable} - fires when the bitset gets changed. */
    get onChanged() { return observeStream(grok_BitSet_Changed(this.d)); }

    /** Finds the value of similarity between two BitSets.
     * @param {BitSet} b - second BitSet.
     * @param {string} metric - similarity metric.
     * @returns {number} */
    similarityTo(b, metric) { return grok_BitSet_SimilarityTo(this.d, b.d, metric); }
}


/** Represents basic descriptive statistics calculated for a {Column}.
 *  See samples: {@link https://public.datagrok.ai/js/samples/data-frame/stats} */
export class Stats {
    constructor(d) { this.d = d; }

    /** Calculates statistics for the specified column.
     * @param {Column} col
     * @returns {Stats} */
    static fromColumn(col) { return new Stats(grok_Stats_FromColumn(col.d)); }

    /** Total number of values (including missing values). */
    get totalCount() { return grok_Stats_Get_TotalCount(this.d); }

    /** Number of missing (empty) values. */
    get missingValueCount() { return grok_Stats_Get_MissingValueCount(this.d); }

    /** Number of non-empty values. */
    get valueCount() { return grok_Stats_Get_ValueCount(this.d); }

    /** @returns {number} - minimum */
    get min() { return grok_Stats_Get_Min(this.d); }

    /** @returns {number} - maximum */
    get max() { return grok_Stats_Get_Max(this.d); }

    /** @returns {number} - sum */
    get sum() { return grok_Stats_Get_Sum(this.d); }

    /** @returns {number} - average */
    get avg() { return grok_Stats_Get_Avg(this.d); }

    /** @returns {number} - standard deviation */
    get stdev() { return grok_Stats_Get_Stdev(this.d); }

    /** @returns {number} - variance */
    get variance() { return grok_Stats_Get_Variance(this.d); }

    /** @returns {number} - skewness */
    get skew() { return grok_Stats_Get_Skew(this.d); }

    /** @returns {number} - kurtosis */
    get kurt() { return grok_Stats_Get_Kurt(this.d); }

    /** @returns {number} - median value */
    get med() { return grok_Stats_Get_Med(this.d); }

    /** @returns {number} - first quartile */
    get q1() { return grok_Stats_Get_Q1(this.d); }

    /** @returns {number} - second quartile */
    get q2() { return grok_Stats_Get_Q2(this.d); }

    /** @returns {number} - third quartile */
    get q3() { return grok_Stats_Get_Q3(this.d); }
}

/**
 * Fluid API for building an aggregation query against a {@link DataFrame}.
 * Build a query by calling the following methods: {@link key}, {@link pivot}, {@link count},
 * {@link uniqueCount}, {@link missingValueCount}, {@link valueCount}, {@link min}, {@link max}, {@link sum},
 * {@link avg}, {@link stdev}, {@link variance}, {@link q1}, {@link q2}, {@link q3}.
 *
 * When the query is constructured, execute it by calling {@link aggregate}, which will
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
    constructor(d) { this.d = d; }

    /** Performs the aggregation
     *  @returns {DataFrame} */
    aggregate() { return new DataFrame(grok_GroupByBuilder_Aggregate(this.d)); }

    /**
     * Performs the aggregation
     * @param {AggregationType} agg - Aggregation type.
     * @param {string} colName - Column name.
     * @param {string=} resultColName - Name of the resulting column. Default value is agg(colName).
     * @returns {GroupByBuilder}
     * */
    add(agg, colName, resultColName = null) {
        grok_GroupByBuilder_Add(this.d, agg, colName, resultColName);
        return this;
    }

    /** Adds a key column to group values on. Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    key(srcColName, resultColName = null) { return this.add(AGG.KEY, srcColName, resultColName); }

    /** Adds a column to pivot values on. Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    pivot(srcColName, resultColName = null) { return this.add(AGG.PIVOT, srcColName, resultColName); }

    /** Adds an aggregation that counts rows, including these will null values.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    count(resultColName = 'count') { return this.add(AGG.TOTAL_COUNT, null, resultColName); }

    /** Adds an aggregation that counts number of unique values in the specified column.
     * See also {@link count}, {@link valueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    uniqueCount(srcColName, resultColName = null) { return this.add(AGG.UNIQUE_COUNT, srcColName, resultColName); }

    /** Adds an aggregation that counts number of missing values in the speficied column.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    missingValueCount(srcColName, resultColName = null) { return this.add(AGG.MISSING_VALUE_COUNT, srcColName, resultColName); }

    /** Adds an aggregation that counts rows, including these will null values.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    valueCount(srcColName, resultColName = null) { return this.add(AGG.VALUE_COUNT, srcColName, resultColName); }

    /** Adds an aggregation that calculates minimum value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    min(srcColName, resultColName = null) { return this.add(AGG.MIN, srcColName, resultColName); }

    /** Adds an aggregation that calculates maximum value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    max(srcColName, resultColName = null) { return this.add(AGG.MAX, srcColName, resultColName); }

    /** Adds an aggregation that calculates sum of the values for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    sum(srcColName, resultColName = null) { return this.add(AGG.SUM, srcColName, resultColName); }

    /** Adds an aggregation that calculates median value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    med(srcColName, resultColName = null) { return this.add(AGG.MED, srcColName, resultColName); }

    /** Adds an aggregation that calculates average value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    avg(srcColName, resultColName = null) { return this.add(AGG.AVG, srcColName, resultColName); }

    /** Adds an aggregation that calculates standard deviation for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    stdev(srcColName, resultColName = null) { return this.add(AGG.STDEV, srcColName, resultColName); }

    /** Adds an aggregation that calculates varians for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    variance(srcColName, resultColName = null) { return this.add(AGG.VARIANCE, srcColName, resultColName); }

    /** Adds an aggregation that calculates first quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q1(srcColName, resultColName = null) { return this.add(AGG.Q1, srcColName, resultColName); }

    /** Adds an aggregation that calculates second quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q2(srcColName, resultColName = null) { return this.add(AGG.Q2, srcColName, resultColName); }

    /** Adds an aggregation that calculates third quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q3(srcColName, resultColName = null) { return this.add(AGG.Q3, srcColName, resultColName); }

    /**
     * @param {BitSet} bitset
     * @returns {GroupByBuilder} */
    whereRowMask(bitset) { grok_GroupByBuilder_WhereBitSet(this.d, bitset.d); return this; }
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
    static getQ(x) {
        _qnumBuf.setFloat64(0, x);
        return _qnumBuf.getInt8(7) & 0x03;
    }

    /**
     * Extracts the value from x, stripping the qualifier .
     * See also {@link getQ}
     * @param {number} x
     * @returns {number}
     * */
    static getValue(x) {
        _qnumBuf.setFloat64(0, qnum);
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
    static create(value, q = QNUM_EXACT) {
        _qnumBuf.setFloat64(0, value);
        let last = _qnumBuf.getInt8(7);
        _qnumBuf.setInt8(7, (last & 0xFC) | q);
        return _qnumBuf.getFloat64(0);
    }

    static exact(x) { return Qnum.create(x, QNUM_EXACT) };
    static less(x) { return Qnum.create(x, QNUM_LESS) };
    static greater(x) { return Qnum.create(x, QNUM_GREATER); }

    /**
     * Parses a string into a qualified number.
     * @param {string} s
     * @returns {number}
     * */
    static parse(s) { return grok_Qnum_Parse(s); }

    /**
     * Converts a qualified number to a string representation.
     * @param {number} x
     * @returns {string}
     * */
    static toString(x) { return grok_Qnum_ToString(x); }
}
