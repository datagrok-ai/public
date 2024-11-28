import {Column, DataFrame} from "./dataframe";
import {toJs} from "./wrappers";
import {FuncCall, Functions} from "./functions";
import {CsvImportOptions, DemoDatasetName, JOIN_TYPE, JoinType, StringPredicate, SyncType, TYPE} from "./const";
import {DataConnection, TableQueryBuilder} from "./entities";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = <any>window;

/** Provides convenient file shares access **/
export class Files {

  /** Reads a table from file. If file contains more than one table, reads the first one.
   * @param {string} path
   * @returns {Promise<DataFrame>}*/
  openTable(path: string): Promise<DataFrame> {
    return api.grok_Files_OpenTable(path);
  }

  /** Reads all tables from file
   * @param {string} path
   * @returns {Promise<Array<DataFrame>>}*/
  openTables(path: string): Promise<Array<DataFrame>> {
    return api.grok_Files_OpenTables(path);
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

  /** Demographics:
   * - subj: int
   * - study: string
   * - site: string
   * - age: int
   * - sex: string
   * - race: string
   * - disease: string
   * - started: date
   * - height: float
   * - weight: float */
  demog(rows: number = 10000): DataFrame { return this.getDemoTable('demog', rows); }

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
    return api.grok_GetDemoTable(path);
  }
}

export class Db {

  /** Executes a specified {@link sql} against the specified {@link connectionId}.
   * @param {string} connectionId - fully-qualified connection name (see [nqName])
   * @param {string} sql - SQL statement
   */
  async query(connectionId: string, sql: string): Promise<DataFrame> {
    let connection: DataConnection = await new Functions().eval(connectionId);
    let q = connection.query('adhoc', sql);
    return await q.apply();
  }

  /**
   * Creates {@link TableQueryBuilder} that can be used to construct sql queries.
   * @param connectionId - fully-qualified connection name (see [nqName])
   * @param tableName - database table name
   */
  buildQuery(connectionId: string, tableName: string): TableQueryBuilder {
    return TableQueryBuilder.from(tableName, connectionId);
  }
}

/**
 * Creating, loading, querying, manipulating, joining tables.
 * */

export class Data {
  public demo: DemoDatasets = new DemoDatasets();
  public files: Files = new Files();
  public db: Db = new Db();

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
    return toJs(api.grok_TestData(dataset, rows, columns));
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
    return toJs(api.grok_ParseCsv(csv, options));
  }

  /**
   * Loads table from the specified URL.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/load-csv}
   * @param {string} csvUrl
   * @returns {Promise<DataFrame>}
   * */
  loadTable(csvUrl: string): Promise<DataFrame> {
    return api.grok_LoadDataFrame(csvUrl);
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
  compareTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[]):
    { changedColumns: number,
      changedValues: number,
      diffTable: DataFrame,
      missingRows: number,
      t1MissingRows: number,
      t2MissingRows: number } {
    console.log(t1.toCsv());
    console.log(t2.toCsv());
    let r = toJs(api.grok_CompareTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2));
    console.log(r.diffTable);
    return r;
  };

  /**
   * Merges two tables by the specified key columns.
   * @param {DataFrame} t1 - a table to join
   * @param {DataFrame} t2 - a table to join
   * @param {string[]} keyColumns1 - key column names from the first table
   * @param {string[]} keyColumns2 - key column names from the second table
   * @param {string[]} valueColumns1 - column names to copy from the first table.
   * Pass null to add all columns, an empty array [] to not add any columns, or an array with column names to add them specifically.
   * @param {string[]} valueColumns2 - column names to copy from the second table
   * @param {JoinType} joinType - inner, outer, left, or right. See [DG.JOIN_TYPE]
   * @param {boolean} inPlace - merges content in-place into the source table
   * @returns {DataFrame}
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/join-link/join-tables}
   * */
  joinTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[] | null = null, valueColumns2: string[] | null = null, joinType: JoinType = JOIN_TYPE.INNER, inPlace: boolean = false): DataFrame {
    return new DataFrame(api.grok_JoinTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }

  /**
   * Opens a table by its id.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/open-table-by-id}
   * @param {string} id - table GUID
   * @returns {Promise<DataFrame>}
   */
  openTable(id: string): Promise<DataFrame> {
    return api.grok_OpenTable(id);
  }

  /** Executes a parameterized query.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/parameterized-query}
   */
  query(queryName: string, queryParameters: object | null = null, adHoc: boolean = false): Promise<DataFrame | any> {
    return api.grok_Query(queryName, queryParameters, adHoc);
  }

  callQuery(queryName: string, queryParameters: object | null = null, adHoc: boolean = false): Promise<FuncCall> {
    return api.grok_CallQuery(queryName, queryParameters, adHoc);
  }

  detectSemanticTypes(t: DataFrame): Promise<void> {
    return api.grok_DetectSemanticTypes(t.dart);
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
