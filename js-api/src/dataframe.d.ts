import {Observable} from 'rxjs';
import {AGG, COLUMN_TYPE, SEMTYPE, TYPE} from './const';

/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 * In the API, the terms "Table" and "DataFrame" are used interchangeably.
 * See usage samples: https://public.datagrok.ai/js/samples/data-frame/manipulate
 */
export class DataFrame {


    /** Returns number of rows in the table.
     * @returns {number} */
    get rowCount(): number;

    /** Returns a {@link BitSet} with selected rows.
     * @returns {BitSet} */
    get selection(): BitSet;

    /** Returns a {@link BitSet} with rows that pass filter.
     * @returns {BitSet} */
    get filter(): BitSet

    /** Name of the dataframe.
     * @returns {string}*/
    get name(): string;

    set name(s: string);

    /** List of columns.
     * @returns {ColumnList} */
    get columns(): ColumnList

    /** List of rows.
     * @returns {RowList} */
    get rows(): RowList

    /** Current row
     * @type {Row} */
    get currentRow(): Row

    //TODO: is it true that it accepts Row?
    set currentRow(idx: Row);

    /** Current column
     * @type {Column} */
    get currentCol(): Column;

    set currentCol(col: Column);

    /** Current cell
     * @type {Cell} */
    get currentCell(): Cell

    set currentCell(cell: Cell)

    /** @returns {Observable} */
    get onValuesChanged(): Observable<any>

    /** @returns {Observable} */
    get onCurrentRowChanged(): Observable<any>

    /** @returns {Observable} */
    get onMouseOverRowChanged(): Observable<any>

    /** @returns {Observable} */
    get onCurrentColChanged(): Observable<any>

    /** @returns {Observable} */ get onMouseOverColChanged(): Observable<any>

    /** @returns {Observable} */ get onCurrentCellChanged(): Observable<any>

    /** @returns {Observable} */ get onMouseOverRowGroupChanged(): Observable<any>

    /** @returns {Observable} */ get onNameChanged(): Observable<any>

    /** @returns {Observable} */ get onMetadataChanged(): Observable<any>

    /** @returns {Observable} */ get onColumnNameChanged(): Observable<any>

    /** @returns {Observable} */ get onColumnSelectionChanged(): Observable<any>

    /** @returns {Observable} */ get onColumnsChanged(): Observable<any>

    /** @returns {Observable} */ get onColumnsAdded(): Observable<any>

    /** @returns {Observable} */ get onColumnsRemoved(): Observable<any>

    /** @returns {Observable} */ get onRowsAdded(): Observable<any>

    /** @returns {Observable} */ get onRowsRemoved(): Observable<any>

    /** @returns {Observable} */ get onSemanticTypeDetecting(): Observable<any>

    /** @returns {Observable} */ get onSemanticTypeDetected(): Observable<any>

    /** @returns {Observable} */ get onDataChanged(): Observable<any>

    /** @returns {Observable} */ get onSelectionChanged(): Observable<any>

    /** @returns {Observable} */ get onFilterChanged(): Observable<any>

    /** Creates a {@link DataFrame} with the specified number of rows and no columns.
     * @param {number} rowCount
     * @returns {DataFrame} */
    static create(rowCount?: number): DataFrame

    /** Creates a {@link DataFrame} from the specified columns. All columns should be of the same length.
     * @param {Column[]} columns
     * @returns {DataFrame} */
    static fromColumns(columns: Column[]): DataFrame

    /** Constructs {@link DataFrame} from a comma-separated values string
     * @param {string} csv - The content of the comma-separated values file.
     * @returns {DataFrame} */
    static fromCsv(csv: string): DataFrame

    /** Constructs {@link DataFrame} from the specified JSON string.
     * @param {string} json - JSON document.
     * @returns {DataFrame} */
    static fromJson(json: string): DataFrame

    /** Creates a [DataFrame] from list of objects by using object keys as column names,
     * and object values as values.
     * @param {object[]} list - List of objects. **/
    static fromObjects(list: object[]): DataFrame[]

    toString(): string;

    /** Returns the value of the specified tag, or null if it does not exist.
     * @returns {string} */
    getTag(tag: string): string;

    /** Sets a tag to the specified value.
     * @param {string} tag - Key.
     * @param {string} value - Value. */
    setTag(tag: string, value: string): void;

    /** Returns i-th row.
     * @param {number} i - Row index.
     * @returns {Row} */
    row(i: number): Row

    /** Returns idx-th value of the specified columns.
     * @param {string} name - Column name.
     * @param {number} idx - Row index. */
    //TODO: any
    get(name: string, idx: number): any;

    /** Sets idx-th value of the specified columns.
     * @param {string} name - Column name.
     * @param {number} idx - Row index.
     * @param value - Value. */
    //TODO: any
    set(name: string, idx: number, value: any): any

    /** Returns a {@link Column} with the specified name.
     * @param {string} name - Column name.
     * @returns {Column} */
    col(name: string): Column

    /** Returns a {@link Cell} with the specified name.
     * @param {number} idx - Row index.
     * @param {string} name - Column name.
     * @returns {Cell} */
    cell(idx: number, name: string): Cell

    /** Same as {@link col}, but throws Error if column is not found
     * @param {string} name - Column name.
     * @returns {Column} */
    getCol(name: string): Column

    /** Exports the content to comma-separated-values format.
     * @returns {string} */
    toCsv(): string;

    /** Creates a new dataframe from the specified row mask and a list of columns.
     * @param {BitSet} rowMask - Rows to include.
     * @param {string[]} columnIds - Columns to include.
     * @param {boolean} saveSelection - Whether selection should be saved. */
    clone(rowMask?: BitSet | null, columnIds?: string[] | null, saveSelection?: boolean): DataFrame

    /** Converts a column with the specified name to [newType],
     * removes the original column from its dataframe and adds the new column to it.
     * @param {string|Column} column
     * @param {string} newType
     * @param {string=} format
     * @returns {Column} */
    changeColumnType(column: string | Column, newType: string,
                     format?: string | null): Column

    /** Begins building a query, using the specified columns as keys.
     * @param {string[]} columnNames - Names of the columns to be used as keys.
     * @returns {GroupByBuilder}
     *  */
    groupBy(columnNames?: string[]): GroupByBuilder

    append(t2: DataFrame, inPlace?: boolean): DataFrame;


    /** @returns {Observable} */
    _event(event: any): Observable<any>

    fireValuesChanged(): void;
}

/** Represents a row. Allows for quick property access like "row.height". */
export class Row {

    /** Creates a {@link Row} from the specified {@link DataFrame} and row index.
     * @param {DataFrame} table
     * @param {number} idx
     * @returns {Row} */
    constructor(table: DataFrame, idx: number)

    /** Returns th */
    get(name: string, idx: number): any
}

/** Strongly-typed column. */
export class Column {

    /** Column data type. */
    get type(): COLUMN_TYPE

    /** Number of elements */
    get length(): number

    /** Parent table */
    get dataFrame(): DataFrame

    /** Semantic type */
    get semType(): SEMTYPE

    set semType(s: SEMTYPE)

    /** Name */
    get name(): string

    set name(s: string)

    /** Returns all unique strings in a sorted order. Applicable to string column only.
     * @returns {string[]} */
    get categories(): string[];

    /** Column's minimum value. The result is cached.
     * @returns {number} */
    get min(): number //TODO: or string?

    /** Column's maximum value. The result is cached.
     * @returns {number} */
    get max(): number //TODO: or string?

    /** Basic descriptive statistics. The result is cached.
     * @returns {Stats} */
    get stats(): Stats

    static fromStrings(name: string, list: string[]): Column

    static fromType(type: COLUMN_TYPE, name?: string, length?: number): Column

    /*
    /** [array] will be not be copied and will be used as column's storage */
    static fromInt32Array(name: string, array: Int32Array, length: number): Column

    /** [array] will be not be copied and will be used as column's storage */
    static fromFloat32Array(name: string, array: Float32Array, length: number): Column

    /** Creates a {@link Column} from the list of values.
     * @ */
    static fromList(type: TYPE, name: string, list: any[]): Column;

    /** Creates an integer column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static int(name: string, length?: number): Column;

    /** Creates a floating point column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static float(name: string, length?: number): Column

    /** Creates a string column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static string(name: string, length?: number): Column

    /** Creates a boolean column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static bool(name: string, length?: number): Column

    /** Creates a qualified number column with the specified name and length.
     *  Initialized values with [values], if it is specified; strips out the qualifier
     *  part if [exact] is true.
     *
     * @param {string} name
     * @param {number} length
     * @param {number[]} values
     * @param {boolean} exact - if true, strips out qualifier from [values].
     * */
    static qnum(name: string, length?: number, values?: number[] | null,
                exact?: boolean): Column

    /** Creates a datetime column with the specified name and length.
     * @param {string} name
     * @param {number} length */
    static dateTime(name: string, length?: number): Column

    /** Returns the raw buffer containing data. Return type depends on the column type:
     * {Int32Array} for ints, {@link INT_NULL} represents null.
     * {Float64Array} for floats, {@link FLOAT_NULL} represents null.
     * {Float64Array} for datetime, in microseconds since epoch, {@link FLOAT_NULL} represents null.
     * {Int32Array} for strings indexes of {@link categories}.
     * {Uint32Array} bit array.
     * @returns {Array} */
    getRawData(): Int32Array | Float64Array | Uint32Array


    setRawData(rawData: any, notify?: boolean): void

    /** Gets i-th value */
    get(i: number): any

    /**
     * Sets [i]-th value to [x]
     * @param {number} i
     * @param x
     */
    set(i: number, x: any): void

    /** Returns whether i-th value is missing.
     * @param {number} i - Row index.
     * @returns {boolean} */
    isNone(i: number): boolean

    /** Gets the value of the specified tag.
     *  @param {string} tag
     *  @returns {string} */
    getTag(tag: string): string;

    /** Sets a tag to the specified value.
     * @param {string} tag - Key.
     * @param {string} value - Value. */
    setTag(tag: string, value: string): void;

    /** Compacts the internal column representation.
     *  Currently, it only affects string columns where values were modified. */
    compact(): void

    /** Copies column values to an array.
     *  @ returns {Array} */
    toList(): Array<any>

    //TODO: check if iterates
    /** An iterator over all values in this column. */
    values(): Generator<any, void, unknown>


    /** Creates and returns a new column by converting [column] to the specified [newType].
     *  @param {string} newType
     *  @param {string} format
     *  @returns {Column} */
    convertTo(newType: string, format?: string | null): Column
}

/** Columns in a [DataFrame]. */
export class ColumnList {

    /** Number of elements in the column.
     * @returns {number} */
    get length(): number

    /** Column with the corresponding name (case-insensitive).
     * @param {string} name - Column name
     * @returns {Column} */
    byName(name: string): Column

    /** Column by index.
     * @param {number} index - Index of the column.
     * @returns {Column} */
    byIndex(index: number): Column

    /** First column of [semType], or null. */
    bySemType(semType: SEMTYPE): Column | null;

    /** Finds columns by the corresponding semTypes, or null, if any of the sem types could not be found.
     * @returns {Column[]} */
    bySemTypesExact(semTypes: SEMTYPE[]): Column[] | null;

    get categorical(): Column[];

    /** Array containing column names.
     * @returns {string[]} */
    names(): string[]

    /** Creates an array of columns.
     * @returns {Column[]} */
    toList(): Column[]

    /** Adds a column, and optionally notifies the parent dataframe.
     * @param {Column} column
     * @param {boolean} notify
     * @returns {Column} */
    add(column: Column, notify?: boolean): Column

    /** Adds an empty column of the specified type.
     * @param {string} name
     * @param {ColumnType} type
     * @returns {Column} */
    addNew(name: string, type: COLUMN_TYPE): Column

    /** Removes column by name (case-insensitive).
     * @param {string} columnName
     * @returns {ColumnList} */
    remove(columnName: string): ColumnList

    /** Checks wheter this list contains a column with the specified name. The check is case-insensitive.
     * @returns {boolean} */
    contains(columnName: string): boolean
}

/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
export class RowList {
    constructor(table: any, d: any);

    /** Removes specified rows
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to remove.
     * @param notify - Whether a change notification should be fired. */
    removeAt(idx: number, count?: number, notify?: boolean): void

    /** Inserts empty rows at the specified position
     * @param {number} idx
     * @param {number} [count=1] - Number of rows to insert.
     * @param notify - Whether a change notification should be fired. */
    insertAt(idx: number, count?: number, notify?: boolean): void

    /** Appends a new row with the specified values
     * @param values - List of values (length and types should match columns)
     * @param notify - Whether a change notification should be fired.
     * @returns {Row} */
    addNew(values?: any[] | null, notify?: boolean): Row

    /** Iterates over all rows. */
    [Symbol.iterator](): Generator<Row, void, unknown>

    /** Sets values for the specified row.
     * @param {number} idx - Row index.
     * @param values - List of values (length and types should match columns) */
    setValues(idx: number, values: any[]): void

    select(rowPredicate: any): void;

    filter(rowPredicate: any): void;
}

/** Represents a table cell. */
export class Cell {

    /** Corresponding table.
     * @returns {DataFrame} */
    get dataFrame(): DataFrame

    /** Corresponding row.
     * @returns {Row} */
    get row(): Row

    /** Index of the corresponding row.
     * @returns {number} */
    get rowIndex(): number

    /** Corresponding column.
     * @returns {Column} */
    get column(): Column

    /** Cell value.
     * @returns {object} */
    get value(): object
}

/** Represents basic descriptive statistics calculated for a {Column}.
 *  See samples: {@link https://public.datagrok.ai/js/samples/data-frame/stats} */
export class Stats {

    /** Total number of values (including missing values). */
    get totalCount(): number

    /** Number of missing (empty) values. */
    get missingValueCount(): number

    /** Number of non-empty values. */
    get valueCount(): number

    /** @returns {number} - minimum */
    get min(): number

    /** @returns {number} - maximum */
    get max(): number

    /** @returns {number} - sum */
    get sum(): number

    /** @returns {number} - average */
    get avg(): number

    /** @returns {number} - standard deviation */
    get stdev(): number

    /** @returns {number} - variance */
    get variance(): number

    /** @returns {number} - skewness */
    get skew(): number

    /** @returns {number} - kurtosis */
    get kurt(): number

    /** @returns {number} - median value */
    get med(): number

    /** @returns {number} - first quartile */
    get q1(): number

    /** @returns {number} - second quartile */
    get q2(): number

    /** @returns {number} - third quartile */
    get q3(): number

    /** Calculates statistics for the specified column.
     * @param {Column} col
     * @returns {Stats} */
    static fromColumn(col: Column): Stats
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

    /** Performs the aggregation
     *  @returns {DataFrame} */
    aggregate(): DataFrame

    /**
     * Performs the aggregation
     * @param {AggregationType} agg - Aggregation type.
     * @param {string} colName - Column name.
     * @param {string=} resultColName - Name of the resulting column. Default value is agg(colName).
     * @returns {GroupByBuilder}
     * */
    add(agg: AGG, colName: string, resultColName?: string | null): GroupByBuilder

    /** Adds a key column to group values on. Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    key(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds a column to pivot values on. Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    pivot(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that counts rows, including these will null values.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    count(resultColName?: 'count' | string): GroupByBuilder;

    /** Adds an aggregation that counts number of unique values in the specified column.
     * See also {@link count}, {@link valueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    uniqueCount(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that counts number of missing values in the speficied column.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    missingValueCount(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that counts rows, including these will null values.
     * See also {@link count}, {@link valueCount}, {@link uniqueCount}, {@link missingValueCount}
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    valueCount(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates minimum value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    min(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates maximum value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    max(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates sum of the values for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    sum(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates median value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    med(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates average value for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    avg(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates standard deviation for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    stdev(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates varians for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    variance(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates first quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q1(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates second quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q2(srcColName: string, resultColName?: string | null): GroupByBuilder

    /** Adds an aggregation that calculates third quartile for the specified column.
     * Call {@link aggregate} when the query is constructed.
     * @param {string} srcColName - column name in the source table
     * @param {string} [resultColName] - column name in the resulting DataFrame
     * @returns {GroupByBuilder} */
    q3(srcColName: string, resultColName?: string | null): GroupByBuilder
}

/**
 * Efficient bit storage and manipulation.
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation}
 */
export class BitSet {
    /** Number of bits in a bitset
     * @type {number} */
    get length(): number

    /** Number of set bits
     * @type {number} */
    get trueCount(): number

    /** Number of unset bits
     * @type {number}*/
    get falseCount(): number

    /** @returns {Observable} - fires when the bitset gets changed. */
    get onChanged(): Observable<{}>

    /** Creates a {BitSet} of the specified length with all bits set to false.
     * @param {number} length - Number of bits.
     * @returns {BitSet} */
    static create(length: number): BitSet

    /** Clones a bitset
     *  @returns {BitSet} */
    clone(): BitSet

    /** Inverts a bitset.
     * @returns {BitSet} */
    invert(): BitSet

    /** Sets all bits to x
     * @param {boolean} x
     * @param {boolean} notify
     * @returns {BitSet} */
    setAll(x: boolean, notify?: boolean): BitSet

    /** Finds the first index of value x, going forward from i-th position.
     * @param {number} i - index
     * @param {boolean} x
     * @returns {number} */
    findNext(i: number, x: boolean): number

    /** Finds the first index of value x, going forward from i-th position, or -1 if not found.
     * @param {number} i - Index to start searching from.
     * @param {boolean} x - Value to search for.
     * @returns {number} */
    findPrev(i: number, x: boolean): number

    /** Gets i-th bit
     * @param {number} i
     * @returns {boolean} */
    get(i: number): boolean

    /** Sets i-th bit to x
     * @param {number} i
     * @param {boolean} x
     * @param {boolean} notify
     * */
    set(i: number, x: boolean, notify?: boolean): BitSet

    /** Sets [i]-th bit to [value], does not check bounds */
    setFast(i: number, value: boolean): void;

    /** Sets all bits by setting i-th bit to the results of f(i)
     * @param {Function} f  */
    init(f: (i: number) => boolean): void;

    /** Indexes of all set bits. The result is cached.
     *  @returns {Int32Array} */
    getSelectedIndexes(): Int32Array

    /** Copies the content from the other {BitSet}.
     * @param {BitSet} b - BitSet to copy from.
     * @returns {BitSet} */
    copyFrom(b: BitSet): BitSet

    fireChanged(): void
}
