import { DataFrame } from "./dataframe";
import { toJs } from "./wrappers";
import { FuncCall } from "./functions";
let api = window;
/** Provides convenient file shares access **/
export class Files {
  /** Reads table from file. If file contains more than one table - reads first
   * @param {string} path
   * @returns {Promise<DataFrame>}*/
  openTable(path) {
    return new Promise((resolve, reject) => api.grok_Files_OpenTable(path, (t) => resolve(toJs(t)), (e) => reject(e)));
  }
  /** Reads all tables from file
   * @param {string} path
   * @returns {Promise<Array<DataFrame>>}*/
  openTables(path) {
    return new Promise((resolve, reject) => api.grok_Files_OpenTables(path, (t) => resolve(t.map(toJs)), (e) => reject(e)));
  }
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
   * "geo": geographic coordinates given as latitude/longitude pairs: lat, lng, value
   *
   * @returns {DataFrame} */
  getDemoTable(dataset, rows = 10000, columns = 3) {
    return new DataFrame(api.grok_TestData(dataset, rows, columns));
  }
  /** Wearable sensor data: time, x, y, z, temp, eda
   * @returns {DataFrame}*/
  biosensor(rows = 10000) {
    return this.getDemoTable('biosensor', rows);
  }
  /** Random walk
   * @returns {DataFrame}*/
  randomWalk(rows = 10000, columns = 3) {
    return this.getDemoTable('random walk', rows, columns);
  }
  /** Demographics
   * @returns {DataFrame}*/
  demog(rows = 10000) {
    return this.getDemoTable('demog', rows);
  }
  /** Plate well data
   * @returns {DataFrame}*/
  wells(rows = 10000) {
    return this.getDemoTable('wells', rows);
  }
  /** Lat/lng of a walk around San Francisco
   * @returns {DataFrame}*/
  geo(rows = 10000) {
    return this.getDemoTable('geo', rows);
  }
  /** Returns a demo dataset with the specified path (relative to the demo root)
   * @example
   * grok.data.getDemoTable("sensors/eeg.csv").then((t) => grok.shell.addTableView(t));
   * @returns {Promise<DataFrame>}*/
  loadDemoTable(path) {
    return new Promise((resolve, reject) => api.grok_GetDemoTable(path, (t) => resolve(toJs(t)), (e) => reject(e)));
  }
}
/** Creating, loading, querying, manipulating, joining tables. */
export class Data {
  constructor() {
    /** @member {DemoDatasets} */
    this.demo = new DemoDatasets();
    this.files = new Files();
  }
  /**
   * Creates a generic dataset with the defined number of rows and columns.
   * [dataset] allowed values:
   * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
   * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
   * "biosensor": wearable sensor data: time, x, y, z, temp, eda
   * "random walk": random walk data for the specified number of dimensions
   * "geo": geographic coordinates given as latitude/longitude pairs: lat, lng, value
   *
   * @returns {DataFrame} */
  testData(dataset, rows = 10000, columns = 3) {
    return new DataFrame(api.grok_TestData(dataset, rows, columns));
  }
  getDemoTable(path) {
    return new Promise((resolve, reject) => api.grok_GetDemoTable(path, (t) => resolve(toJs(t)), (e) => reject(e)));
  }
  /**
   * @typedef {Object} CsvImportOptions
   * @property {string} delimiter
   * @property {string} decimalSeparator
   * @property {Object} thousandSeparator
   **/
  /**
   * Parses the CSV string.
   * @param {string} csv - The content of the comma-separated values file.
   * @param {CsvImportOptions} options
   * @returns {DataFrame}
   * */
  parseCsv(csv, options) {
    return new DataFrame(api.grok_ParseCsv(csv, options));
  }
  /**
   * Loads table from the specified URL.
   * @param {string} csvUrl
   * @returns {Promise<DataFrame>}
   * */
  loadTable(csvUrl) {
    return new Promise((resolve, reject) => api.grok_LoadDataFrame(csvUrl, (t) => resolve(toJs(t, false)), (e) => reject(e)));
  }
  /**
   * @param {DataFrame} t1
   * @param {DataFrame} t2
   * @param {string[]} keyColumns1
   * @param {string[]} keyColumns2
   * @param {SyncType[]} linkTypes
   * */
  linkTables(t1, t2, keyColumns1, keyColumns2, linkTypes) {
    api.grok_LinkTables(t1.d, t2.d, keyColumns1, keyColumns2, linkTypes);
  }
  ;
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
    api.grok_CompareTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2);
  }
  ;
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
   * @returns {DataFrame}
   * */
  joinTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace) {
    return new DataFrame(api.grok_JoinTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }
  /**
   * Opens a table by its id.
   * @param {string} id - table GUID
   * @returns {Promise<DataFrame>}
   */
  openTable(id) {
    return new Promise((resolve, reject) => api.grok_OpenTable(id, (t) => resolve(new DataFrame(t)), (e) => reject(e)));
  }
  query(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
    return new Promise((resolve, reject) => api.grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t) => resolve(toJs(t)), (e) => reject(e)));
  }
  callQuery(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
    return new Promise((resolve, reject) => api.grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c) => resolve(new FuncCall(c)), (e) => reject(e)));
  }
  detectSemanticTypes(t) {
    return new Promise((resolve, reject) => api.grok_DetectSemanticTypes(t.d, (_) => resolve(), (e) => reject(e)));
  }
}
export class Detector {
  /**
   * Calls [check] function against a random subset of the column values, returns true
   * if all checks return true. Useful for the efficient auto-detection of the column semantic type.
   *
   * @param {Column} column
   * @param {StringPredicate} check
   * @param {number} min - minimum number of categories. Returns false if less than that.
   * @param {number} max - number of checks to make
   * @returns {boolean}
   * */
  static sampleCategories(column, check, min = 5, max = 10) {
    let categories = column.categories;
    if (categories.length > min)
      return false;
    for (let i = 0; i < Math.max(max, categories.length); i++) {
      let value = Math.floor(Math.random() * categories.length);
      if (value !== null && !check(categories[value]))
        return false;
    }
    return true;
  }
}
