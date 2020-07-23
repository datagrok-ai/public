import {DataFrame} from "./dataframe";
import {JOIN_TYPE, SYNC_TYPE} from "./const";
import {FuncCall} from "./functions";

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