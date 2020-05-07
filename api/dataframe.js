
/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 */
class DataFrame {
    constructor(d) { this.d = d; }
    static create(rowCount = 0) { return new DataFrame(grok_DataFrame(rowCount)); }
    static fromColumns(columns) { return new DataFrame(grok_DataFrame_FromColumns(columns.map((c) => c.d))); }
    static fromCsv(csv) { return grok.parseCsv(csv); }
    static fromJson(json) { return new DataFrame(grok_DataFrame_FromJson(json)); }

    toString() { return `${this.name} (${this.rowCount} rows, ${this.columns.length} columns)` };

    get rowCount() { return grok_DataFrame_RowCount(this.d); }
    get selection() { return new BitSet(grok_DataFrame_Get_Selection(this.d)); }
    get filter() { return new BitSet(grok_DataFrame_Get_Filter(this.d)); }

    get name() { return grok_DataFrame_Get_Name(this.d); }
    set name(s) { grok_DataFrame_Set_Name(this.d, s); }

    getTag(tag) { return grok_DataFrame_Get_Tag(this.d, tag); }
    setTag(tag, value) { grok_DataFrame_Set_Tag(this.d, tag, value); }

    get columns() { return new ColumnList(grok_DataFrame_Columns(this.d)); }
    get rows() { return new RowList(this, grok_DataFrame_Rows(this.d)); }

    row(i) { return new Row(this, i); }

    get(name, idx) { return this.getCol(name).get(idx); }
    set(name, idx, value) { this.getCol(name).set(idx, value); }

    col(name) { return new Column(grok_DataFrame_ColumnByName(this.d, name)); }

    cell(idx, name) { return new Cell(grok_DataFrame_Cell(this.d, idx, name)); }

    /// Same as [col], but throws Error if column is not found
    getCol(name) {
        let c = this.col(name);
        if (c === null)
            throw new Error(`No such column: ${name}`);
        return c;
    }

    toCsv() { return grok_DataFrame_ToCsv(this.d); }

    clone(rowMask, columnIds, saveSelection) {
        return new DataFrame(grok_DataFrame_Clone(this.d, rowMask.d, columnIds, saveSelection));
    }

    get currentRow() { return new Row(this, grok_DataFrame_Get_CurrentRowIdx(this.d)); }
    set currentRow(idx) { grok_DataFrame_Set_CurrentRowIdx(this.d, idx); }

    get currentCol() { return new Column(grok_DataFrame_Get_CurrentCol(this.d)); }
    set currentCol(col) { grok_DataFrame_Set_CurrentCol(this.d, col.d); }

    get currentCell() { return new Cell(grok_DataFrame_Get_CurrentCell(this.d)); }
    set currentCell(cell) { grok_DataFrame_Set_CurrentCell(this.d, cell.d); }

    /** Creates a [DataFrame] from list of objects **/
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

    groupBy(columnNames = []) { return new GroupByBuilder(grok_DataFrame_GroupBy(this.d, columnNames)); }

    _event(event, callback) { return _sub(grok_DataFrame_OnEvent(this.d, event, callback)); }

    onValuesChanged(callback) { return this._event('ddt-values-changed', callback); }
    onCurrentRowChanged(callback) { return this._event('ddt-current-row-changed', callback); }
    onMouseOverRowChanged(callback) { return this._event('ddt-mouse-over-row-changed', callback); }
    onCurrentColChanged(callback) { return this._event('ddt-current-col-changed', callback); }
    onMouseOverColChanged(callback) { return this._event('ddt-mouse-over-col-changed', callback); }
    onCurrentCellChanged(callback) { return this._event('ddt-current-cell-changed', callback); }
    onMouseOverRowGroupChanged(callback) { return this._event('ddt-mouse-over-row-group-changed', callback); }
    onNameChanged(callback) { return this._event('ddt-table-name-changed', callback); }
    onMetadataChanged(callback) { return this._event('ddt-table-metadata-changed', callback); }
    onColumnNameChanged(callback) { return this._event('ddt-table-column-name-changed', callback); }
    onColumnSelectionChanged(callback) { return this._event('ddt-column-selection-changed', callback); }
    onColumnsChanged(callback) { return this._event('ddt-columns-changed', callback); }
    onColumnsAdded(callback) { return this._event('ddt-columns-added', callback); }
    onColumnsRemoved(callback) { return this._event('ddt-columns-removed', callback); }
    onRowsRemoved(callback) { return this._event('ddt-rows-removed', callback); }

    onDataChanged(callback) { return _sub(grok_DataFrame_OnDataChanged(this.d, callback)); }
    onFilterChanged(callback) { return _sub(grok_DataFrame_OnFilterChanged(this.d, callback)); }
    onSelectionChanged(callback) { return _sub(grok_DataFrame_OnSelectionChanged(this.d, callback)); }

    fireValuesChanged() { grok_DataFrame_FireValuesChanged(this.d); }
}

/** Represents a row. Allows for quick property access like "row.height". */
class Row {
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

    get(name, idx) { return this.table.getCol(name).get(idx); }
}

/** Strongly-typed column. */
class Column {
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
    static int(name, length = 0) { return Column.fromType(TYPE_INT, name, length); }

    /** Creates a floating point column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static float(name, length = 0) { return Column.fromType(TYPE_FLOAT, name, length); }

    /** Creates a string column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static string(name, length = 0) { return Column.fromType(TYPE_STRING, name, length); }

    /** Creates a boolean column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static bool(name, length = 0) { return Column.fromType(TYPE_BOOL, name, length); }

    /** Creates a datetime column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static dateTime(name, length = 0) { return Column.fromType(TYPE_DATE_TIME, name, length); }

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

    /** Gets i-th value */
    get(i) { return grok_Column_GetValue(this.d, i); }

    /**
     * Sets [i]-th value to [x]
     * @param {number} i
     * @param x
     */
    set(i, x) { grok_Column_SetValue(this.d, i, x); }

    isNone(i) { return grok_Column_IsNone(this.d, i); }

    getTag(tag) { return grok_Column_Get_Tag(this.d, tag); }
    setTag(tag, value) { grok_Column_Set_Tag(this.d, tag, value); }

    compact() { return grok_Column_Compact(this.d); }

    toList() { return grok_Column_ToList(this.d); }

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
class ColumnList {
    constructor(d) { this.d = d; }

    get length() { return grok_ColumnList_Length(this.d); }
    byName(name) { return new Column(grok_ColumnList_ByName(this.d, name)); }
    byIndex(index) { return new Column(grok_ColumnList_ByIndex(this.d, index)); }

    /** First column of [semType], or null. */
    bySemType(semType) { return new Column(grok_ColumnList_BySemType(this.d, semType)); }

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
class RowList {
    constructor(table, d) { this.table = table; this.d = d; }

    removeAt(idx, count = 1, notify = true) { grok_RowList_RemoveAt(this.d, idx, count, notify); }
    insertAt(idx, count = 1, notify = true) { grok_RowList_InsertAt(this.d, idx, count, notify); }
    addNew(values = null, notify = true) { return new Row(this.table, grok_RowList_AddNew(this.d, values, notify));}

    *iterator() {
        for (let i = 0; i < this.table.rowCount; i++)
            yield new Row(this.table, i);
    }

    setValues(idx, values) { grok_RowList_SetValues(this.d, idx, values); }
}

class Cell {
    constructor(d) { this.d = d; }

    get dataFrame() { return new DataFrame(grok_Cell_Get_DataFrame(this.d)); }
    get row() { return new Row(this.dataFrame, this.rowIndex); }
    get rowIndex() { return grok_Cell_Get_RowIndex(this.d); }
    get column() { return new Column(grok_Cell_Get_Column(this.d)); }
    get value() { return grok_Cell_Get_Value(this.d); }
}


class Stats {
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

class GroupByBuilder {
    constructor(d) { this.d = d; }

    aggregate() { return new DataFrame(grok_GroupByBuilder_Aggregate(this.d)); }

    add(agg, col, resultColName = null) { grok_GroupByBuilder_Add(this.d, agg, col, resultColName); return this; }

    key(col, resultColName = null) { this.add(AGG_KEY, col, resultColName); return this; }
    pivot(col, resultColName = null) { this.add(AGG_PIVOT, col, resultColName); return this; }
    count(col, resultColName = null) { this.add(AGG_TOTAL_COUNT, col, resultColName); return this; }
    uniqueCount(col, resultColName = null) { this.add(AGG_UNIQUE_COUNT, col, resultColName); return this; }
    missingValueCount(col, resultColName = null) { this.add(AGG_MISSING_VALUE_COUNT, col, resultColName); return this; }
    valueCount(col, resultColName = null) { this.add(AGG_VALUE_COUNT, col, resultColName); return this; }
    min(col, resultColName = null) { this.add(AGG_MIN, col, resultColName); return this; }
    max(col, resultColName = null) { this.add(AGG_MAX, col, resultColName); return this; }
    sum(col, resultColName = null) { this.add(AGG_SUM, col, resultColName); return this; }
    med(col, resultColName = null) { this.add(AGG_MED, col, resultColName); return this; }
    avg(col, resultColName = null) { this.add(AGG_AVG, col, resultColName); return this; }
    stdev(col, resultColName = null) { this.add(AGG_STDEV, col, resultColName); return this; }
    variance(col, resultColName = null) { this.add(AGG_VARIANCE, col, resultColName); return this; }
    q1(col, resultColName = null) { this.add(AGG_Q1, col, resultColName); return this; }
    q2(col, resultColName = null) { this.add(AGG_Q2, col, resultColName); return this; }
    q3(col, resultColName = null) { this.add(AGG_Q3, col, resultColName); return this; }
}

/** Efficient bit storage and manipulation. */
class BitSet {
    constructor(d) { this.d = d; }

    static create(length) { return new BitSet(grok_BitSet(length)); }

    /** Number of bits in a bitset */
    get length() { return grok_BitSet_Get_Length(this.d); }

    /** Number of set bits */
    get trueCount() { return grok_BitSet_Get_TrueCount(this.d); }

    /** Number of unset bits */
    get falseCount() { return grok_BitSet_Get_FalseCount(this.d); }

    /** Clones a bitset */
    clone() { return new BitSet(grok_BitSet_Clone(this.d)); }

    /** Inverts a bitset */
    invert() { return grok_BitSet_Invert(this.d); }

    /** Sets all bits to x */
    setAll(x) { return grok_BitSet_SetAll(this.d, x); }

    /** Finds the first index of value x, going forward from i-th position */
    findNext(i, x) { return grok_BitSet_FindNext(this.d, i, x); }

    /** Finds the first index of value x, going forward from i-th position */
    findPrev(i, x) { return grok_BitSet_FindPrev(this.d, i, x); }

    /** Gets i-th bit */
    get(i) { return grok_BitSet_GetBit(this.d, i); }

    /** Sets i-th bit to x */
    set(i, x) { return grok_BitSet_SetBit(this.d, i, x); }

    copyFrom(b) { return grok_BitSet_CopyFrom(this.d, b.d); }

    /** Registers a callback that gets invoked when the bitset is changed */
    onChanged(callback) { return _sub(grok_BitSet_OnChanged(this.d, callback)); }
}
