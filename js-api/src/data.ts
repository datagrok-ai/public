import {Column, DataFrame} from "./dataframe";
import {toJs} from "./wrappers";
import {FuncCall, Functions} from "./functions";
import {CsvImportOptions, DemoDatasetName, JOIN_TYPE, JoinType, StringPredicate, SyncType, TYPE} from "./const";
import {ColumnInfo, DataConnection, TableInfo, TableQueryBuilder, Property} from "./entities";
import {IDartApi} from "./api/grok_api.g";
import { Grid } from "./grid";
import {Tags} from "./api/ddt.api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

declare let grok: any;
declare let DG: any;

/** Provides convenient file shares access **/
export class Files {

  /** Reads a table from file. If file contains more than one table, reads the first one. */
  openTable(path: string): Promise<DataFrame> {
    return api.grok_Files_OpenTable(path);
  }

  /** Reads all tables from file */
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

  async getInfo(connection: DataConnection): Promise<DbInfo> {
    return grok.dapi.connections.getDatabaseInfo(connection);
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
   * */
  parseCsv(csv: string, options?: CsvImportOptions): DataFrame {
    return toJs(api.grok_ParseCsv(csv, options));
  }

  /**
   * Loads table from the specified URL.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/load-csv} */
  loadTable(csvUrl: string): Promise<DataFrame> {
    return api.grok_LoadDataFrame(csvUrl);
  }

  /**
   * Links tables by the specified key columns using the specified link types (such as "current row to filter", see {@link DG.SYNC_TYPE}).
   * Tables are synchronized on the first change, set the {@link initialSync} option to reflect the current table state according to the sync type.
   * */
  linkTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], linkTypes: SyncType[], initialSync: boolean = false, filterAllOnNoRowsSelected = false): void {
    api.grok_LinkTables(t1.dart, t2.dart, keyColumns1, keyColumns2, linkTypes, initialSync, filterAllOnNoRowsSelected);
  };

  /**
   * Opens the visual table comparison tool.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/compare-tables}
   * */
  compareTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[], valueColumns2: string[], addToWorkspace: boolean = false):
    { changedColumns: number,
      changedValues: number,
      diffTable: DataFrame,
      diffGrid: Grid,
      missingRows: number,
      t1MissingRows: number,
      t2MissingRows: number } {
    return toJs(api.grok_CompareTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2, addToWorkspace));
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
   * Sample: {@link https://public.datagrok.ai/js/samples/data-frame/join-link/join-tables}
   * */
  joinTables(t1: DataFrame, t2: DataFrame, keyColumns1: string[], keyColumns2: string[], valueColumns1: string[] | null = null, valueColumns2: string[] | null = null, joinType: JoinType = JOIN_TYPE.INNER, inPlace: boolean = false): DataFrame {
    return new DataFrame(api.grok_JoinTables(t1.dart, t2.dart, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
  }

  /**
   * Opens a table by its id.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/open-table-by-id}
   * @param {string} id - table GUID
   */
  openTable(id: string): Promise<DataFrame> {
    return api.grok_OpenTable(id);
  }

  /**
   * Executes a Datagrok query and returns the result as the specified type.
   *
   * @template T - The expected return type of the query (defaults to {@link DataFrame}).
   *
   * @param queryName - The fully qualified name of the query to execute
   *   (e.g., `"DbTests:PostgresqlScalarInt"`).
   * @param queryParameters - Optional key-value pairs that will be passed as
   *   parameters to the query. Defaults to `null` if no parameters are needed.
   * @param adHoc - **Deprecated.** Indicates whether the query should be
   *   executed in ad-hoc mode. This parameter will be removed soon and should
   *   not be used.
   *
   * @returns A promise that resolves to the query result, converted to type `T`.
   *
   * @example
   * // Returns a number
   * const count = await grok.data.query<number>("DbTests:PostgresqlScalarInt");
   *
   * @example
   * // Returns a string
   * const text = await grok.data.query<string>("DbTests:PostgresqlScalarStr");
   *
   * @example
   * // Returns a DataFrame (default)
   * const df = await grok.data.query("DbTests:PostgresqlTable");
   *
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/parameterized-query}
   */
  async query<T = DataFrame>(queryName: string,
                             queryParameters: object | null = null,
                             /**
                              * @deprecated Parameter adHoc will be removed soon.
                              */
                             adHoc: boolean = false): Promise<T> {
    return toJs(await api.grok_CallFunc(queryName, queryParameters, true, null));
  }

  callQuery(queryName: string,
            queryParameters: object | null = null,
            /**
             * @deprecated Parameter adHoc will be removed soon.
             */
            adHoc: boolean = false): Promise<FuncCall> {
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

/**
 * Represents metadata for a specific database schema in Datagrok.
 *
 * `DbSchemaInfo` provides access to schema-level information such as tables,
 * user and LLM comments, and utilities for annotating tables and columns.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbSchemaInfo {
  public dart: any;

  /**
   * Creates a new `DbSchemaInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a database schema.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * Schema name (e.g., `"public"`, `"dbo"`).
   *
   * @returns The name of this schema.
   */
  get name(): string {
    return api.grok_DbSchemaInfo_Get_Name(this.dart);
  }

  /**
   * @returns A comment string, or an empty string if none exists.
   */
  get comment(): string | undefined {
    return api.grok_DbSchemaInfo_Get_Prop(this.dart, Tags.DbComment);
  }

  /**
   * @returns User-defined comment used by AI query builder.
   */
  get llmComment(): string | undefined {
    return api.grok_DbSchemaInfo_Get_Prop(this.dart, Tags.LlmComment);
  }

  /**
   * The parent `DataConnection` this schema belongs to.
   *
   * @returns A {@link DataConnection} instance.
   */
  get connection(): DataConnection {
    return api.grok_DbSchemaInfo_Get_Connection(this.dart);
  }

  /**
   * Lists all tables belonging to this schema.
   *
   * @returns A promise resolving to an array of {@link TableInfo} objects.
   */
  getTables(): Promise<TableInfo[]> {
    return api.grok_DbSchemaInfo_Get_Tables(this.dart);
  }

  // ─────────────────────────────── SETTERS ───────────────────────────────

  setComment(comment: string): Promise<void> {
    return api.grok_DbSchemaInfo_Set_Prop(this.dart, Tags.DbComment, comment);
  }

  setLlmComment(llmComment: string): Promise<void> {
    return api.grok_DbSchemaInfo_Set_Prop(this.dart, Tags.LlmComment, llmComment);
  }

  // ─────────────────────────────── ANNOTATION METHODS ───────────────────────────────

  /**
   * Adds or updates metadata for a specific table in this schema.
   * It's just only meta annotation which is stored in Datagrok and real database is not affected.
   *
   * @param table - Either a {@link TableInfo} instance or a table name.
   * @param props - Table-level properties such as comments, row counts, and domains.
   *
   * @returns A promise that resolves when the annotation is saved.
   *
   * @example
   * ```ts
   * await schema.annotateTable("customers", {
   *   comment: "Master table containing customer profiles",
   *   rowCount: 120345,
   * });
   * ```
   */
  annotateTable(
      table: TableInfo | string,
      props: DbTableProperties,
  ): Promise<void> {
    return api.grok_DbSchemaInfo_AnnotateTable(
        this.dart,
        table instanceof TableInfo ? table.friendlyName : table,
        props,
    );
  }

  /**
   * Adds or updates metadata for a column within a table.
   * It's just only meta annotation which is stored in Datagrok and real database is not affected.
   *
   * @param table - A {@link TableInfo} instance or table name.
   * @param column - A {@link ColumnInfo} instance or column name.
   * @param props - Column-level annotations such as uniqueness, min/max values, and comments.
   *
   * @returns A promise that resolves when the annotation is saved.
   *
   * @example
   * ```ts
   * await schema.annotateColumn("orders", "order_id", {
   *   isUnique: true,
   *   comment: "Primary key for the orders table",
   * });
   * ```
   */
  annotateColumn(
      table: TableInfo | string,
      column: ColumnInfo | string,
      props: DbColumnProperties,
  ): Promise<void> {
    return api.grok_DbSchemaInfo_AnnotateColumn(
        this.dart,
        table instanceof TableInfo ? table.friendlyName : table, //  table.friendlyName corresponds to real db table name
        column instanceof ColumnInfo ? column.name : column, // column.name corresponds to real db col name
        props,
    );
  }
}


/**
 * Represents metadata for a database relation in Datagrok.
 *
 * `DbRelationInfo` describes how two tables are connected, including the
 * schema/table pairs, participating columns, cardinality, and whether this
 * relation is considered a primary navigation path.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbRelationInfo implements DbRelationProperties {
  public dart: any;

  /**
   * Creates a new `DbRelationInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a relation.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * User-defined comment describing this relation.
   */
  get comment(): string | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbComment);
  }

  /**
   * LLM-generated annotation or AI summary for this relation.
   */
  get llmComment(): string | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.LlmComment);
  }

  /**
   * The relation cardinality (e.g., `"one-to-many"`, `"many-to-one"`).
   */
  get cardinality(): 'one-to-one' | 'one-to-many' | 'many-to-one' | 'many-to-many' | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelCardinality);
  }

  /**
   * Indicates whether this relation is treated as a primary navigation path
   * (e.g., preferred for automated join suggestions).
   */
  get isPrimaryPath(): boolean | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelIsPrimaryPath);
  }

  /**
   * Connection object this relation belongs to.
   */
  get connection(): DataConnection {
    return toJs(api.grok_DbRelationInfo_Get_Connection(this.dart));
  }

  /**
   * Schema name of the source table of this relation.
   */
  get fromSchema(): string | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelFromSchema);
  }

  /**
   * Source table name of this relation.
   */
  get fromTable(): string | undefined {
      return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelFromTable);
  }

  /**
   * List of source table columns participating in the relation.
   */
  get fromColumns(): string[] | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelFromColumns);
  }

  /**
   * Schema name of the target table of this relation.
   */
  get toSchema(): string | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelToSchema);
  }

  /**
   * Target table name of this relation.
   */
  get toTable(): string | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelToTable);
  }

  /**
   * List of target table columns participating in the relation.
   */
  get toColumns(): string[] | undefined {
    return api.grok_DbRelationInfo_Get_Prop(this.dart, Tags.DbRelToColumns);
  }
}


export interface DbRelationProperties {
  comment?: string;
  llmComment?: string;
  cardinality?: 'one-to-one' | 'one-to-many' | 'many-to-one' | 'many-to-many';
  isPrimaryPath?: boolean;
  fromSchema?: string;
  toSchema?: string;
}

export interface DbColumnProperties {
  comment?: string;
  llmComment?: string;
  isUnique?: boolean;
  min?: number;
  max?: number;
  values?: string[];
  sampleValues?: string[];
  uniqueCount?: number;
  quality?: string; // Semantic type of this column
}

export interface DbTableProperties {
  comment?: string;
  llmComment?: string;
  domains?: string;
  rowCount?: number;
}

/**
 * Represents metadata for a database connection in Datagrok.
 *
 * Provides access to high-level database metadata such as schemas,
 * relations, comments, and LLM-generated annotations. It also allows updating
 * connection-level annotations and adding new relation metadata.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbInfo {
  public dart: any;

  /**
   * Creates a new `DbInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a database connection.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * The name of the database.
   *
   * @returns The database name.
   */
  get name(): string {
    return api.grok_DbInfo_Get_Name(this.dart);
  }

  /**
   * @returns Connection comment text, or an empty string if none exists.
   */
  get comment(): string | undefined {
    return api.grok_DbInfo_Get_Prop(this.dart, Tags.DbComment);
  }

  /**
   * @returns User-defined comment used by AI query builder.
   */
  get llmComment(): string | undefined {
    return api.grok_DbInfo_Get_Prop(this.dart, Tags.LlmComment);
  }

  /**
   * Lists all schemas available for this database connection.
   *
   * @returns A promise resolving to an array of {@link DbSchemaInfo} objects.
   */
  getSchemas(): Promise<DbSchemaInfo[]> {
    return api.grok_DbInfo_Get_Schemas(this.dart);
  }

  /**
   * Lists all known relations for this connection.
   *
   * @returns A promise resolving to an array of {@link DbRelationInfo} objects.
   */
  getRelations(): Promise<DbRelationInfo[]> {
    return api.grok_DbInfo_Get_Relations(this.dart);
  }

  /**
   * The underlying `DataConnection` object.
   *
   * @returns A {@link DataConnection} wrapper.
   */
  get connection(): DataConnection {
    return toJs(api.grok_DbInfo_Get_Connection(this.dart));
  }

  // ─────────────────────────────── SETTERS ───────────────────────────────

  setComment(comment: string): Promise<void> {
    return api.grok_DbInfo_Set_Prop(this.dart, Tags.DbComment, comment);
  }

  setLlmComment(comment: string): Promise<void> {
    return api.grok_DbInfo_Set_Prop(this.dart, Tags.LlmComment, comment);
  }


  /**
   * Creates and registers a logical (metadata-only) relation between two tables
   * within this database connection.
   *
   * This method does **not** modify the physical database schema. The relation is
   * stored only as a Datagrok metadata annotation and is used for navigation,
   * visualization, and AI-assisted features.
   *
   * @param fromTable - Source (child) table name.
   * @param fromColumns - Column(s) in the source table that participate in the relation.
   * @param toTable - Target (parent) table name.
   * @param toColumns - Column(s) in the target table that are referenced.
   * @param props - Optional relation metadata (e.g., cardinality, comments, primary path flag).
   *
   * @returns A promise that resolves to the newly created {@link DbRelationInfo}.
   *
   * @example
   * ```ts
   * const rel = await dbInfo.addRelation(
   *   'orders',
   *   ['customer_id'],
   *   'customers',
   *   ['id'],
   *   {
   *     cardinality: 'many-to-one'
   *   }
   * );
   * ```
   */
  addRelation(fromTable: string, fromColumns: string[], toTable: string, toColumns: string[],
              props?: DbRelationProperties): Promise<DbRelationInfo> {
    return api.grok_DbInfo_AddRelation(this.dart, fromTable, fromColumns, toTable, toColumns, props);
  }

  /**
   * Removes all the meta data associated with the connection in Datagrok.
   */
  clearProperties(): Promise<void> {
    return api.grok_DbInfo_ClearProperties(this.dart);
  }
}

/**
 * Represents data source, such as PostgreSQL, Oracle, S3, etc.
 */
export class ConnectionDataSource {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /**
   * Use {@link DataConnection#dataSource} as input.
   */
  static byType(type: string): ConnectionDataSource | undefined {
    return api.grok_ConnectionDataSource_ByType(type);
  }

  /**
   * Data source name (such as `'Oracle'`).
   * {@link DataConnection#dataSource} refers to this value.
   */
  get type(): string {
    return api.grok_ConnectionDataSource_Type(this.dart);
  }

  /**
   * Category (such as 'database').
   */
  get category(): string | undefined {
    return api.grok_ConnectionDataSource_Category(this.dart);
  }

  /**
   * Free-text description.
   */
  get description(): string | undefined {
    return api.grok_ConnectionDataSource_Description(this.dart);
  }

  /**
   * Indicates if this data provider can work right from the browser.
   * When true, a query always executes on a server. Examples: Oracle, Postgres, MS SQL
   * When false, a query runs on a server or in a browser, depending on the context (queries initiated
   * from the browser would run in a browser).
   */
  get requiresServer(): boolean {
    return api.grok_ConnectionDataSource_RequiresServer(this.dart);
  }

  /**
   * Comment start, depends on SQL engine. Applicable only to database data sources.
   */
  get commentStart(): string | undefined {
    return api.grok_ConnectionDataSource_CommentStart(this.dart);
  }

  /**
   * Name brackets, required if name contains spaces, depends on SQL engine. Applicable only to database data sources.
   */
  get nameBrackets(): string | undefined {
    return api.grok_ConnectionDataSource_NameBrackets(this.dart);
  }

  /**
   * Whether database schemas can be browsed from UI. Applicable only to database data sources.
   */
  get canBrowseSchema(): boolean {
    return api.grok_ConnectionDataSource_CanBrowseSchema(this.dart);
  }

  /**
   * Applicable only to database data sources.
   */
  get queryLanguage(): string | undefined {
    return api.grok_ConnectionDataSource_QueryLanguage(this.dart);
  }

  get connectionTemplate(): Property[] {
    return api.grok_ConnectionDataSource_ConnectionTemplate(this.dart);
  }

  get credentialsTemplate(): Property[] {
    return api.grok_ConnectionDataSource_CredentialsTemplate(this.dart);
  }
}
