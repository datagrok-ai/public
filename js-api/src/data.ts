import {DataFrame, Column} from "./dataframe";
import {toJs} from "./wrappers";
import {FuncCall} from "./functions";
import {TYPE, DemoDatasetName, JoinType, SyncType, CsvImportOptions, StringPredicate} from "./const";

let api = <any>window;

/** Provides convenient file shares access **/
export class Files {

  /** Reads a table from file. If file contains more than one table, reads the first one.
   * @param {string} path
   * @returns {Promise<DataFrame>}*/
  openTable(path: string): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_Files_OpenTable(path, (t: any) => resolve(toJs(t)), (e: any) => reject(e)));
  }

  /** Reads all tables from file
   * @param {string} path
   * @returns {Promise<Array<DataFrame>>}*/
  openTables(path: string): Promise<Array<DataFrame>> {
    return new Promise((resolve, reject) => api.grok_Files_OpenTables(path, (t: any) => resolve(t.map(toJs)), (e: any) => reject(e)));
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
  getDemoTable(dataset: DemoDatasetName, rows: number = 10000, columns: number = 3): DataFrame {
    return toJs(api.grok_TestData(dataset, rows, columns));
  }

  /** Wearable sensor data: time, x, y, z, temp, eda
   * @returns {DataFrame}*/
  biosensor(rows: number = 10000): DataFrame {
    return this.getDemoTable('biosensor', rows);
  }

  /** Random walk
   * @returns {DataFrame}*/
  randomWalk(rows: number = 10000, columns: number = 3): DataFrame {
    return this.getDemoTable('random walk', rows, columns);
  }

  /** Demographics
   * @returns {DataFrame}*/
  demog(rows: number = 10000): DataFrame {
    return this.getDemoTable('demog', rows);
  }

  /** Demographics
   * @returns {DataFrame}*/
  molecules(rows: number = 10000): DataFrame {
    return this.getDemoTable('molecules', rows);
  }

  /** Plate well data
   * @returns {DataFrame}*/
  wells(rows: number = 10000): DataFrame {
    return this.getDemoTable('wells', rows);
  }

  /** Lat/lng of a walk around San Francisco
   * @returns {DataFrame}*/
  geo(rows: number = 10000): DataFrame {
    return this.getDemoTable('geo', rows);
  }

  /** Dose-response data */
  doseResponse(rows: number = 10000): DataFrame {
    return this.getDemoTable('dose-response', rows);
  }

  /** Returns a demo dataset with the specified path (relative to the demo root)
   * @example
   * grok.data.getDemoTable("sensors/eeg.csv").then((t) => grok.shell.addTableView(t));
   * @returns {Promise<DataFrame>}*/
  loadDemoTable(path: string): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_GetDemoTable(path, (t: any) => resolve(toJs(t)), (e: any) => reject(e)));
  }
}

/**
 * Creating, loading, querying, manipulating, joining tables.
 * */

export class Data {
  public demo: DemoDatasets;
  public files: Files;

  constructor() {
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
  testData(dataset: DemoDatasetName, rows: number = 10000, columns: number = 3): DataFrame {
    return new DataFrame(api.grok_TestData(dataset, rows, columns));
  }

  getDemoTable(path: string): Promise<DataFrame> {
    return api.grok_GetDemoTable(path);
  }

  /**
   * Parses the CSV string.
   * @param {string} csv - The content of the comma-separated values file.
   * @param {CsvImportOptions} options
   * @returns {DataFrame}
   * */
  parseCsv(csv: string, options?: CsvImportOptions): DataFrame {
    return new DataFrame(api.grok_ParseCsv(csv, options));
  }

  /**
   * Loads table from the specified URL.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/load-csv}
   * @param {string} csvUrl
   * @returns {Promise<DataFrame>}
   * */
  loadTable(csvUrl: string): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_LoadDataFrame(csvUrl, (t: any) => resolve(toJs(t, false)), (e: any) => reject(e)));
  }

  /**
   * Links tables by the specified key columns using the specified link types (such as "current row to filter", see {@link DG.SYNC_TYPE}).
   * Tables are synchronized on the first change, set the {@link initialSync} option to reflect the current table state according to the sync type.
   * */
  linkTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], linkTypes: SyncType[], initialSync: boolean = false): void {
    api.grok_LinkTables(t1.dart, t2.dart, keyColumns1, keyColumns2, linkTypes, initialSync);
  };

  /**
   * Opens the visual table comparison tool.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/compare-tables}
   * */
  compareTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[]): void {
    api.grok_CompareTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2);
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
   * @returns {DataFrame}
   * */
  joinTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[], joinType: JoinType, inPlace: boolean): DataFrame {
    return new DataFrame(api.grok_JoinTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }

  /**
   * Opens a table by its id.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/open-table-by-id}
   * @param {string} id - table GUID
   * @returns {Promise<DataFrame>}
   */
  openTable(id: string): Promise<DataFrame> {
    return new Promise((resolve, reject) => api.grok_OpenTable(id, (t: any) => resolve(new DataFrame(t)), (e: any) => reject(e)));
  }

  /** Executes a parameterized query.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/parameterized-query} 
   */
  query(queryName: string, queryParameters: object | null = null, adHoc: boolean = false, pollingInterval: number = 1000): Promise<DataFrame | any> {
    return new Promise((resolve, reject) => api.grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t: any) => resolve(toJs(t)), (e: any) => reject(e)));
  }

  callQuery(queryName: string, queryParameters: object | null = null, adHoc: boolean = false, pollingInterval: number = 1000): Promise<FuncCall> {
    return new Promise((resolve, reject) => api.grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c: any) => resolve(new FuncCall(c)), (e: any) => reject(e)));
  }

  detectSemanticTypes(t: DataFrame): Promise<void> {
    return new Promise((resolve, reject) => api.grok_DetectSemanticTypes(t.dart, (_: any) => resolve(), (e: any) => reject(e)));
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
   * @param {number} ratio - [0-1] range: minimum allowed number of the success/total checks.
   * @param minStringLength - values shorter than that are not considered checks
   * @returns {boolean}
   * */
  static sampleCategories(column: Column, check: StringPredicate, min: number = 5, max: number = 10, ratio: number = 1, minStringLength: number = 1): boolean {
    if (column.type !== TYPE.STRING)
      return false;

    let categories = column.categories;
    if (categories.length < min)
      return false;
    
    if (!categories.some(cat => cat.length >= minStringLength)) 
      return false;
      
    let checksCount = Math.min(max, categories.length);
    let maxFailed = checksCount * (1 - ratio);
    let failed = 0;

    for (let i = 0; i < checksCount; i++) {
      let categoryIndex = Math.floor(categories.length * i / checksCount);
      let s = categories[categoryIndex];
      if (s.length >= minStringLength && !check(s) && ++failed > maxFailed)
        return false;
    }

    return true;
  }
}
