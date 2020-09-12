import {DataFrame} from "./dataframe";
import {JOIN_TYPE, SYNC_TYPE} from "./const";
import {FuncCall} from "./functions";

export const QNUM_LESS: number;
export const QNUM_EXACT: number;
export const QNUM_GREATER: number;

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
    static getQ(x: number): number

    /**
     * Extracts the value from x, stripping the qualifier .
     * See also {@link getQ}
     * @param {number} x
     * @returns {number}
     * */
    static getValue(x: number): number

    /**
     * Creates a QNum value out of the [value] and qualifier [q].
     * @param {number} value
     * @param {number} q
     * @returns {number}
     * */
    static create(value: number, q?: number): number

    static exact(x: number): number
    static less(x: number): number
    static greater(x: number): number

    /**
     * Parses a string into a qualified number.
     * @param {string} s
     * @returns {number}
     * */
    static parse(s: string): number

    /**
     * Converts a qualified number to a string representation.
     * @param {number} x
     * @returns {string}
     * */
    static toString(x: number): string
}


type DemoDatasetName = 'wells' | 'demog' | 'biosensor' | 'random walk';

/** Provides convenient access to demo datasets. */
export class DemoDatasets {

    /**
     * Creates a generic dataset with the defined number of rows and columns.
     * [dataset] allowed values:
     * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
     * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
     * "biosensor": wearable sensor data: time, x, y, z, temp, eda
     * "random walk": random walk data for the specified number of dimensions
     *
     * @returns {DataFrame} */
    getDemoTable(dataset: DemoDatasetName, rows?: number, columns?: number): DataFrame

    /** Wearable sensor data: time, x, y, z, temp, eda
     * @returns {DataFrame}*/
    biosensor(rows?: number): DataFrame

    /** Random walk
     * @returns {DataFrame}*/
    randomWalk(rows ?: number, columns?: number): DataFrame

    /** Demographics
     * @returns {DataFrame}*/
    demog(rows?: number): DataFrame

    /** Plate well data
     * @returns {DataFrame}*/
    wells(rows?: number): DataFrame

    /** Lat/lng of a walk around San Francisco
     * @returns {DataFrame}*/
    geo(rows?: number): DataFrame


    /** Returns a demo dataset with the specified path (relative to the demo root)
     * @example
     * grok.data.getDemoTable("sensors/eeg.csv").then((t) => grok.shell.addTableView(t));
     * @returns {Promise<DataFrame>}*/
    loadDemoTable(path: string): Promise<DataFrame>
}

/** Creating, loading, querying, manipulating, joining tables. */
export class Data {

    /**
     * Creates a generic dataset with the defined number of rows and columns.
     * [dataset] allowed values:
     * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
     * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
     * "biosensor": wearable sensor data: time, x, y, z, temp, eda
     * "random walk": random walk data for the specified number of dimensions
     *
     * @returns {DataFrame} */
    testData(dataset: DemoDatasetName, rows?: number, columns?: number): DataFrame

    getDemoTable(path: string): Promise<DataFrame>

    /**
     * Parses the CSV string.
     * @returns {DataFrame}
     * */
    parseCsv(csv: string): DataFrame

    /**
     * Loads table from the specified URL.
     * @param {string} csvUrl
     * @returns {Promise<DataFrame>}
     * */
    loadTable(csvUrl: string): Promise<DataFrame>

    /**
     * @param {DataFrame} t1
     * @param {DataFrame} t2
     * @param {string[]} keyColumns1
     * @param {string[]} keyColumns2
     * @param {SyncType[]} linkTypes
     * @returns {Promise<DataFrame>}
     * */
    linkTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], linkTypes: SYNC_TYPE[]): Promise<DataFrame>

    /**
     * Opens the visual table comparison tool.
     * @param {DataFrame} t1
     * @param {DataFrame} t2
     * @param {string[]} keyColumns1
     * @param {string[]} keyColumns2
     * @param {string[]} valueColumns1
     * @param {string[]} valueColumns2
     * */
    compareTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[]): void

    /**
     * Merges two tables by the specified key columns.
     * @param {DataFrame} t1
     * @param {DataFrame} t2
     * @param {string[]} keyColumns1
     * @param {string[]} keyColumns2
     * @param {string[]} valueColumns1
     * @param {string[]} valueColumns2
     * @param {JoinType} joinType
     * @param {boolean} inPlace - merges content in-place into the source table
     * */
    joinTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[], joinType: JOIN_TYPE, inPlace: boolean): DataFrame

    /**
     * Merges two tables by the specified key columns.
     * @param {string} id - table GUID
     * @returns {Promise<DataFrame>}
     */
    openTable(id: string): Promise<DataFrame>

    query(queryName: string, queryParameters?: Object | null, adHoc?: boolean, pollingInterval?: number): Promise<DataFrame>

    callQuery(queryName: string, queryParameters?: Object | null, adHoc?: boolean, pollingInterval?: number): Promise<FuncCall>

    detectSemanticTypes(t: any): Promise<void>

}