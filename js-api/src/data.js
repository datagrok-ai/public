import {DataFrame} from "./dataframe";
import {toJs} from "./wrappers";
import {FuncCall} from "./functions";

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
        _qnumBuf.setFloat64(0, value.toDouble());
        let last = _qnumBuf.getInt8(7);
        _qnumBuf.setInt8(7, (last & 0xFC) | q);
        return _qnumBuf.getFloat64(0);
    }

    static exact(x) { return create(x, QNUM_EXACT) };
    static less(x) { return create(x, QNUM_LESS) };
    static greater(x) { return create(x, QNUM_GREATER); }

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
    getDemoTable(dataset, rows = 10000, columns = 3) {
        return new DataFrame(grok_TestData(dataset, rows, columns));
    }

    /** Wearable sensor data: time, x, y, z, temp, eda
     * @returns {DataFrame}*/
    biosensor(rows = 10000) { return this.getDemoTable('biosensor', rows); }

    /** Random walk
     * @returns {DataFrame}*/
    randomWalk(rows = 10000, columns = 3) { return this.getDemoTable('random walk', rows, columns); }

    /** Demographics
     * @returns {DataFrame}*/
    demog(rows = 10000) { return this.getDemoTable('demog', rows); }

    /** Plate well data
     * @returns {DataFrame}*/
    wells(rows = 10000) { return this.getDemoTable('wells', rows); }

    /** Returns a demo dataset with the specified path (relative to the demo root)
     * @example
     * grok.data.getDemoTable("sensors/eeg.csv").then((t) => grok.shell.addTableView(t));
     * @returns {Promise<DataFrame>}*/
    loadDemoTable(path) {
        return new Promise((resolve, reject) => grok_GetDemoTable(path, (t) => resolve(toJs(t))));
    }
}

/** Creating, loading, querying, manipulating, joining tables. */
export class Data {

    constructor() {
        this.demo = new DemoDatasets();
    }

    /**
     * Creates a generic dataset with the defined number of rows and columns.
     * [dataset] allowed values:
     * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
     * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
     * "biosensor": wearable sensor data: time, x, y, z, temp, eda
     * "random walk": random walk data for the specified number of dimensions
     *
     * @returns {DataFrame} */
    testData(dataset, rows = 10000, columns = 3) {
        return new DataFrame(grok_TestData(dataset, rows, columns));
    }

    getDemoTable(path) {
        return new Promise((resolve, reject) => grok_GetDemoTable(path, (t) => resolve(toJs(t))));
    }

    /**
     * Parses the CSV string.
     * @returns {DataFrame}
     * */
    parseCsv(csv) { return new DataFrame(grok_ParseCsv(csv)); }

    /**
     * Loads table from the specified URL.
     * @param {string} csvUrl
     * @returns {Promise<DataFrame>}
     * */
    loadTable(csvUrl) {
        return new Promise((resolve, reject) => grok_LoadDataFrame(csvUrl, (t) => resolve(toJs(t, false))));
    }

    /**
     * @param {DataFrame} t1
     * @param {DataFrame} t2
     * @param {string[]} keyColumns1
     * @param {string[]} keyColumns2
     * @param {SyncType[]} linkTypes
     * @returns {Promise<DataFrame>}
     * */
    linkTables(t1, t2, keyColumns1, keyColumns2, linkTypes) { grok_LinkTables(t1.d, t2.d, keyColumns1, keyColumns2, linkTypes); };

    /**
     * Opens the visual table comparison tool.
     * @param {DataFrame} t1
     * @param {DataFrame} t2
     * @param {string[]} keyColumns1
     * @param {string[]} keyColumns2
     * @param {string[]} valueColumns1
     * @param {string[]} valueColumns2
     * */
    compareTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2) {
        grok_CompareTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2);
    };

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
    joinTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace) {
        return new DataFrame(grok_JoinTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
    }

    /**
     * Merges two tables by the specified key columns.
     * @param {string} id - table GUID
     * @returns {Promise<DataFrame>}
     */
    openTable(id) { return new Promise((resolve, reject) => grok_OpenTable(id, (t) => resolve(new DataFrame(t)))); }

    query(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t) => resolve(new DataFrame(t))));
    }

    callQuery(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c) => resolve(new FuncCall(c))));
    }

    detectSemanticTypes(t) { return new Promise((resolve, reject) => grok_DetectSemanticTypes(t.d, (_) => resolve())); }

}
