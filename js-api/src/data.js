import {DataFrame} from "./dataframe";
import {toJs} from "./wrappers";
import {FuncCall} from "./functions";

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
    biosensor(rows = 10000) { return getDemoTable('biosensor', rows); }


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

    detectSemanticTypes(t) { return new Promise((resolve, reject) => grok_DetectSematicTypes(t.d, (_) => resolve())); }

}
