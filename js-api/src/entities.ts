// noinspection JSUnusedGlobalSymbols

import {
  AGG, AggregationType,
  COLUMN_TYPE,
  ColumnType,
  JOIN_TYPE,
  JoinType,
  ScriptingLanguage,
  SemType,
  Type,
  TYPE,
  USER_STATUS
} from "./const";
import { FuncCall } from "./functions";
import {toDart, toJs} from "./wrappers";
import type {FileSource} from "./dapi";
import {MapProxy} from "./proxies";
import {DataFrame} from "./dataframe";
import {PackageLogger} from "./logger";
import dayjs from "dayjs";
import {IDartApi} from "./api/grok_api.g";
import {DataSourceType} from "./api/grok_shared.api.g";
import { Tags } from "./api/ddt.api.g";
import {View, ViewBase} from "./views/view";
import {InputType} from "./api/d4.api.g";
import {Observable} from "rxjs";
import {observeStream} from "./events";

declare var grok: any;
declare var DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


type PropertyGetter<TSource = any, TProp = any> = (a: TSource) => TProp;
type PropertySetter<TSource = any, TProp = any> = (a: TSource, value: TProp) => void;
type ValueValidator<T> = (value: T) => string;

/**
 * Represents basic properties of database connection.
 *
 * Note: The actual list of supported parameters varies by connector. For a complete list,
 * see {@link https://datagrok.ai/help/access/databases/connectors/}
 */
export interface DatabaseConnectionProperties {
  /** Database server hostname or IP address */
  server?: string;
  /** Database server port number */
  port?: number;
  /** Database name */
  db?: string;
  /** Full connection string (alternative to individual parameters) */
  connString?: string;
}

/**
 * Represents connection cache properties.
 * See also: {@link https://datagrok.ai/help/develop/how-to/function_results_cache}
 */
export interface DataConnectionCacheProperties {
  /** Whether to cache query results */
  cacheResults?: boolean;
  /** Whether to cache database schema information */
  cacheSchema?: boolean;
  /** Cron expression for cache invalidation schedule */
  cacheInvalidateSchedule?: string;
}

/**
 * Represents data connection properties.
 *
 * The available properties depend on the {@link DataSourceType}. Common properties are
 * defined explicitly, but each connector may support additional parameters.
 *
 * For a complete list of supported parameters:
 * - Database connectors: {@link https://datagrok.ai/help/access/databases/connectors/}
 * - File shares: {@link https://datagrok.ai/help/access/files/shares/}
 *
 * @example
 * // PostgreSQL connection
 * const props: DataConnectionProperties = {
 *   dataSource: 'Postgres',
 *   server: 'localhost',
 *   port: 5432,
 *   db: 'mydb',
 *   login: 'user',
 *   password: 'pass'
 * };
 *
 * @example
 * // S3 file share
 * const props: DataConnectionProperties = {
 *   dataSource: 'S3',
 *   accessKey: 'AKIA...',
 *   secretKey: '...',
 *   region: 'us-east-1',  // connector-specific parameter
 *   bucket: 'my-bucket'   // connector-specific parameter
 * };
 */
export interface DataConnectionProperties extends DatabaseConnectionProperties, DataConnectionCacheProperties {
  /** Data source type identifier (e.g., 'Postgres', 'MySQL', 'S3', 'Git') */
  dataSource: string;
  /** Login/username for authentication */
  login?: string;
  /** Password for authentication */
  password?: string;
  /** AWS access key (for AWS-based connectors) */
  accessKey?: string;
  /** AWS secret key (for AWS-based connectors) */
  secretKey?: string;
  /** AWS region (for AWS-based connectors like S3, Athena) */
  region?: string;
  /** SSL/TLS mode (for database connectors) */
  ssl?: boolean | string;
  /**
   * Additional connector-specific parameters.
   * The available parameters vary by {@link dataSource} type.
   */
  [x: string]: string | number | boolean | undefined;
}

type FieldPredicate = {field: string, pattern: string, dataType: string};
type FieldOrder = {field: string, asc?: boolean};
type GroupAggregation = {aggType: string, colName: string, resultColName?: string, function?: string};
export type TableJoin = {leftTableName: string, rightTableName: string, rightTableAlias?: string, joinType: JoinType, leftTableKeys: string[], rightTableKeys: string[]};
type DockerContainerStatus = 'stopped' | 'started' | 'pending change' | 'changing' | 'error' | 'checking';

/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {

  public dart: any;

  /** @constructs Entity*/
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id(): string { return api.grok_Entity_Get_Id(this.dart); }
  set id(x: string) { api.grok_Entity_Set_Id(this.dart, x); }

  /** Generates new {@link id} for this entity. */
  newId(): void { api.grok_Entity_New_Id(this.dart); }

  /** Entity friendly name
   *  @type {string} */
  get friendlyName(): string { return api.grok_Entity_Get_FriendlyName(this.dart); }
  set friendlyName(x: string) { api.grok_Entity_Set_FriendlyName(this.dart, x); }

  /** Entity short name */
  get name(): string { return api.grok_Entity_Get_Name(this.dart); }
  set name(x: string) { api.grok_Entity_Set_Name(this.dart, x); }

  /** Entity full-qualified name */
  get nqName(): string { return api.grok_Entity_Get_nqName(this.dart); }

  /** Entity path */
  get path(): string { return api.grok_Entity_Path(this.dart); }

  /** Time when entity was created **/
  get createdOn(): dayjs.Dayjs { return dayjs(api.grok_Entity_Get_CreatedOn(this.dart)); }

  /** Time when entity was updated **/
  get updatedOn(): dayjs.Dayjs | null {
    const d = api.grok_Entity_Get_UpdatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  /** Who created entity **/
  get author(): User { return toJs(api.grok_Entity_Get_Author(this.dart)); }

  /** Entity type name **/
  get entityType(): string { return api.grok_Entity_Get_EntityType(this.dart); }

  /** Gets entity properties */
  getProperties(): Promise<{[index: string]: any}> {
    return api.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.dart, this.dart);
  }

  /** Sets entity properties */
  setProperties(props: {[index: string]: any}): Promise<any> {
    return api.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.dart, this.dart, props);
  }

  /** Returns a string representing the object */
  toString(): string { return api.grok_Object_ToString(this.dart); }

  hasTag(tag: string): boolean { return api.grok_Entity_Has_Tag(this.dart, tag); }

  /** Adds a specified tag */
  tag(tag: string): boolean { return api.grok_Entity_Tag(this.dart, tag); }

  /** Removes a specified tag */
  unTag(tag: string): boolean { return api.grok_Entity_UnTag(this.dart, tag); }
}

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {
  /** @constructs User*/
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new user.
   * Note that it's just a client object, it won't be saved in the database. */
  static create(): User { return new User(api.grok_User_From_Id(null)); };

  static fromId(id: string): User { return new User(api.grok_User_From_Id(id)); }

  /** Returns current user. */
  static current(): User { return new User(api.grok_User()); }

  /** First name */
  get firstName(): string { return api.grok_User_Get_FirstName(this.dart); }
  set firstName(name: string) {api.grok_User_Set_FirstName(this.dart, name);}

  /** Last name */
  get lastName(): string { return api.grok_User_Get_LastName(this.dart); }
  set lastName(name: string) {api.grok_User_Set_LastName(this.dart, name);}

  get status(): string & USER_STATUS { return api.grok_User_Get_Status(this.dart); }
  set status(name: string & USER_STATUS) {api.grok_User_Set_Status(this.dart, name);}

  /** Email */
  get email(): string | null { return api.grok_User_Get_Email(this.dart); }
  set email(email: string | null) {api.grok_User_Set_Email(this.dart, email);}

  /** Picture URL */
  get picture(): string | object { return api.grok_User_Get_Picture(this.dart); }

  /** User home project */
  get project(): Project { return toJs(api.grok_User_Get_Project(this.dart)); }

  /** User home folder connection */
  get home(): DataConnection { return toJs(api.grok_User_Get_Storage(this.dart)); }

  /** Login */
  get login(): string { return api.grok_User_Get_Login(this.dart); }
  set login(login: string) { api.grok_User_Set_Login(this.dart, login); }

  toMarkup(): string { return api.grok_User_ToMarkup(this.dart); }

  /** Security Group */
  get group(): Group { return toJs(api.grok_User_Get_Group(this.dart)); }

  /** Date when user joined */
  get joined(): dayjs.Dayjs { return dayjs(api.grok_User_Get_Joined(this.dart)); }

  static get defaultUsersIds() {
    return {
      "Test": "ca1e672e-e3be-40e0-b79b-d2c68e68d380",
      "Admin": "878c42b0-9a50-11e6-c537-6bf8e9ab02ee",
      "System": "3e32c5fa-ac9c-4d39-8b4b-4db3e576b3c3",
    } as const;
  }

  static get test(): User { return new User(api.grok_User_Test()); }

  static get admin(): User { return new User(api.grok_User_Admin()); }

  static get system(): User { return new User(api.grok_User_System()); }
}


/**
 * Represents a user session in the Datagrok platform.
 * @extends Entity
 * */
export class UserSession extends Entity {
  /** @constructs UserSession*/
  constructor(dart: any) {
    super(dart);
  }

  /** Entity ID (GUID) */
  get id(): string { return api.grok_Entity_Get_Id(this.dart); }
  set id(x: string) { api.grok_Entity_Set_Id(this.dart, x); }

  /** */
  get type(): 'internal' | 'guest' { return api.grok_UserSession_Get_Type(this.dart); }

  /** External Token */
  get externalToken(): string { return api.grok_UserSession_Get_ExternalToken(this.dart); }

  /** User */
  get user(): User { return toJs(api.grok_UserSession_Get_User(this.dart)); }
}


/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/datagrok/functions/function}
 * */
export class Func extends Entity {
  public aux: any;
  public options: { [key: string]: any; };

  constructor(dart: any) {
    super(dart);
    this.aux = new MapProxy(api.grok_Func_Get_Aux(this.dart));
    // @ts-ignore
    this.options = new MapProxy(api.grok_Func_Get_Options(this.dart));
  }

  get description(): string { return api.grok_Func_Get_Description(this.dart); }

  get type(): string { return api.grok_Func_Get_Type(this.dart); }

  get path(): string { return api.grok_Func_Get_Path(this.dart); }

  /** Help URL. */
  get helpUrl(): string { return api.grok_Func_Get_HelpUrl(this.dart); }

  set helpUrl(url: string) { api.grok_Func_Set_HelpUrl(this.dart, url); }

  /** A package this function belongs to. */
  get package(): Package { return api.grok_Func_Get_Package(this.dart); }

  /** Indicates that the function (or script) is already vector, meaning it
   * accepts vector input (an entire column) and processes it in a single call,
   * rather than being executed separately for each scalar element (row) */
  get isVectorFunc(): boolean { return api.grok_Func_Get_IsVectorFunc(this.dart); }

  /** Returns {@link FuncCall} object in a stand-by state */
  prepare(parameters: {[name: string]: any} = {}): FuncCall {
    return toJs(api.grok_Func_Prepare(this.dart, parameters));
  };

  async prepareAsync(parameters: {[name: string]: any} = {}): Promise<FuncCall> {
    const call = await api.grok_Func_PrepareAsync(this.dart, parameters);
    return toJs(call);
  };

  /** Input parameters */
  get inputs(): Property[] {
    return toJs(api.grok_Func_Get_InputParams(this.dart));
  }

  /** Output parameters */
  get outputs(): Property[] {
    return toJs(api.grok_Func_Get_OutputParams(this.dart));
  }

  /**
   *  Executes the function with the specified {@link parameters}, and returns result.
   *  If necessary, the corresponding package will be loaded as part of the call.
   * */
  async apply(parameters: {[name: string]: any} | any[] = {}): Promise<any> {
    parameters ??= {};
    if (Array.isArray(parameters)) {
      if (parameters.length != this.inputs.length)
        throw `${this.name}: expected ${this.inputs.length} parameters, got ${parameters.length}`;
      const params: any = {};
      for (let i = 0; i < this.inputs.length; i++)
        params[this.inputs[i].name] = parameters[i];
      parameters = params;
    }

    return (await (this.prepare(parameters)).call()).getOutputParamValue();
  }

  /** Executes the function synchronously, and returns the result.
   *  If the function is asynchronous, throws an exception. */
  applySync(parameters: {[name: string]: any} = {}): any {
    return this.prepare(parameters).callSync().getOutputParamValue();
  }

  /** Returns functions with the specified attributes. */
  static find(params?: { package?: string, name?: string, tags?: string[], meta?: any, returnType?: string, returnSemType?: string}): Func[] {
    return api.grok_Func_Find(params?.package, params?.name, params?.tags, params?.meta, params?.returnType, params?.returnSemType);
  }

  /**
   * @deprecated Use find, it's the same now but does not make a server query and synchronous.
   */
  static async findAll(params?: { package?: string, name?: string, tags?: string[], meta?: any, returnType?: string, returnSemType?: string}): Promise<Func[]> {
    let functions = Func.find(params);
    let queries = await grok.dapi.queries.include('params,connection').filter(`name="${params?.name}"`).list();
    let scripts = await grok.dapi.scripts.include('params').filter(`name="${params?.name}"`).list();

    return [...functions, ...queries, ...scripts];
  }

  /** Returns a function with the specified name. */
  static byName(name: string): Func {
    return Func.find({ name: name})[0];
  }
}

export interface ProjectOpenOptions {
  closeAll: boolean;
  openViews: 'all' | 'saved' | 'none';
}

/** Represents a project */
export class Project extends Entity {
  constructor(dart: any) {
    super(dart);

    this.options = new MapProxy(api.grok_Project_Get_Options(this.dart), 'options');
  }

  public options: any;

  static create(): Project {return toJs(api.grok_Project_From_Id(null)); };

  get pictureUrl(): string {
    return api.grok_PictureMixin_Get_PictureUrl(this.dart);
  }

  get path(): string {
    return api.grok_Project_Get_Path(this.dart);
  }

  get isOnServer(): string {
    return api.grok_Project_Get_IsOnServer(this.dart);
  }

  get isLocal(): string {
    return api.grok_Project_Get_IsLocal(this.dart);
  }

  /** Project description */
  get description(): string {
    return api.grok_Project_Description(this.dart);
  }

  /** Project changes flag */
  get isDirty(): boolean {
    return api.grok_Project_IsDirty(this.dart);
  }

  /** Project is empty flag */
  get isEmpty(): boolean {
    return api.grok_Project_IsEmpty(this.dart);
  }

  get isDashboard(): boolean {
    return api.grok_Project_IsDashboard(this.dart);
  }

  get isPackage(): boolean {
    return api.grok_Project_IsPackage(this.dart);
  }

  toMarkup(): string {
    return api.grok_Project_ToMarkup(this.dart);
  }

  /** Opens the project in workspace */
  open(options?: ProjectOpenOptions): Promise<Project> {
    let openViews = options?.openViews ?? 'all';
    return api.grok_Project_Open(this.dart, options?.closeAll ?? false, openViews == 'all' || openViews == 'saved', openViews == 'all');
  }

  get links(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, true));
  }

  get children(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, false));
  }

  addLink(entity: Entity | DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, true);
  }

  addChild(entity: Entity |DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, false);
  }

  removeLink(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

  removeChild(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/access#data-query}
 * */
export class DataQuery extends Func {
  /** @constructs DataQuery*/
  constructor(dart: any) {
    super(dart);
  }

  /** @deprecated Use FuncCall.adHoc instead **/
  get adHoc(): boolean { return api.grok_Query_Get_AdHoc(this.dart); }
  /** @deprecated Use FuncCall.adHoc instead **/
  set adHoc(a: boolean) { api.grok_Query_Set_AdHoc(this.dart, a); }

  /** Query text */
  get query(): string { return api.grok_Query_Query(this.dart); }
  set query(q: string) { api.grok_Query_Set_Query(this.dart, q); }

  get connection(): DataConnection { return toJs(api.grok_Query_Get_Connection(this.dart)); }
  set connection(c: DataConnection) { api.grok_Query_Set_Connection(this.dart, toDart(c)); }

  get postProcessScript(): string {return api.grok_Query_Get_PostProcessScript(this.dart);}
  set postProcessScript(script: string) { api.grok_Query_Set_PostProcessScript(this.dart, script);}
  /** Executes query
   * @returns {Promise<DataFrame>} */
  async executeTable(): Promise<DataFrame> { return toJs(await api.grok_Query_ExecuteTable(this.dart)); }
}

/** Represents a table query
 * @extends DataQuery */
export class TableQuery extends DataQuery {
  /** @constructs TableQuery */
  constructor(dart: any) { super(dart); }

  /** Creates a TableQuery
   * @param {DataConnection} connection - DataConnection to query table from
   * @returns {TableQuery} */
  static create(connection: DataConnection): TableQuery {
    return toJs(api.grok_TableQuery_Create(connection.dart));
  }

  /** Table name
   * @type {string} */
  get table(): string { return api.grok_TableQuery_GetTable(this.dart); }
  set table(tableName: string) { api.grok_TableQuery_SetTable(this.dart, tableName); }

  /** Fields array
   * @type {string[]} */
  get fields(): string[] { return api.grok_TableQuery_GetFields(this.dart); }
  set fields(fields: string[]) { api.grok_TableQuery_SetFields(this.dart, fields); }

  /** Executes query */
  async executeTable(): Promise<DataFrame> { return toJs(await api.grok_TableQuery_ExecuteTable(this.dart)); }

  /** Where clauses */
  get where(): FieldPredicate[] { return toJs(api.grok_TableQuery_Get_WhereClausesDB(this.dart)); }
  set where(wl: FieldPredicate[]) { api.grok_TableQuery_Set_WhereClausesDB(this.dart, wl.map(param => toDart(param))); }

  /** Aggregation clauses {queryPartParams} */
  get aggregations(): GroupAggregation[] { return toJs(api.grok_TableQuery_Get_AggregationsDB(this.dart)); }
  set aggregations(wl: GroupAggregation[]) { api.grok_TableQuery_Set_AggregationsDB(this.dart, wl.map(param => toDart(param))); }

  get joins(): TableJoin[] { return toJs(api.grok_TableQuery_Get_Joins(this.dart)); }
  set joins(j: TableJoin[]) { api.grok_TableQuery_Set_Joins(this.dart, j.map(join => toDart(join))); }

  /** Having clauses */
  get having(): FieldPredicate[] { return toJs(api.grok_TableQuery_Get_HavingDB(this.dart)); }
  set having(wl: FieldPredicate[]) { api.grok_TableQuery_Set_HavingDB(this.dart, wl.map(param => toDart(param))); }

  /** Order By clauses */
  get orderBy(): FieldOrder[] { return toJs(api.grok_TableQuery_Get_OrderByDB(this.dart)); }
  set orderBy(wl: FieldOrder[]) { api.grok_TableQuery_Set_OrderByDB(this.dart, wl.map(param => toDart(param))); }

  set limit(rows: number | undefined) { api.grok_TableQuery_Set_Limit(this.dart, rows); }
  get limit(): number | undefined { return api.grok_TableQuery_Get_Limit(this.dart); }

  /** Creates {@link TableQueryBuilder} from table name
   * @param {string} table - Table name
   * @returns {TableQueryBuilder} */
  static from(table: string): TableQueryBuilder {return toJs(api.grok_TableQuery_From(table)); }

  /** Creates {@link TableQueryBuilder} from {@link TableInfo}
   * @param {TableInfo} table - TableInfo object
   * @returns {TableQueryBuilder} */
  static fromTable(table: TableInfo): TableQueryBuilder {return toJs(api.grok_TableQuery_FromTable(table.dart)); }
}

/** Table query builder that works with database tables */
export class TableQueryBuilder {
  dart: any;

  /** @constructs TableQueryBuilder */
  constructor(dart: any) { this.dart = dart; }

  /** Creates {@link TableQueryBuilder} from table name
   * @param {string} table - Table name
   * @param connection - {@link DataConnection} that {@link TableQuery} will use after the build. Can be passed lately directly to {@link TableQuery}.
   * @returns {TableQueryBuilder} */
  static from(table: string, connection?: DataConnection | string): TableQueryBuilder {
    let conn: any = connection === undefined ? connection : connection instanceof DataConnection ? connection.dart : connection;
    return toJs(api.grok_DbTableQueryBuilder_From(table, conn));
  }

  /** Creates {@link TableQueryBuilder} from {@link TableInfo}
   * @param {TableInfo} table - TableInfo object
   * @returns {TableQueryBuilder} */
  static fromTable(table: TableInfo): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_FromTable(table.dart));
  }

  /** Selects all fields of the table
   * @returns {TableQueryBuilder} */
  selectAll(): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_SelectAll(this.dart)); }

  /** Selects specified fields of the table. All fields will be selected if not provided.
   * @param {string[]} fields - Array of fields to select
   * @returns {TableQueryBuilder} */
  select(fields: string[]): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_Select(this.dart, fields)); }

  /**
   * Performs aggregation
   * @param {AggregationType} type - Aggregation type.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  selectAggr(type: AggregationType, field: string | null, fieldAlias: string | null): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_SelectAggr(this.dart, type, field, fieldAlias));
  }

  /**
   * Adds an aggregation that calculates average value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  avg(field: string, fieldAlias: string | null = 'avg'): TableQueryBuilder {
    return this.selectAggr(AGG.AVG, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates minimum value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  min(field: string, fieldAlias: string | null = 'min'): TableQueryBuilder {
    return this.selectAggr(AGG.MIN, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates maximum value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  max(field: string, fieldAlias: string | null = 'max'): TableQueryBuilder {
    return this.selectAggr(AGG.MAX, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates sum of the values for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  sum(field: string, fieldAlias: string | null = 'sum'): TableQueryBuilder {
    return this.selectAggr(AGG.SUM, field, fieldAlias);
  }

  /**
   * Adds an aggregation that counts rows. Equivalent to count(*).
   * @param {string} fieldAlias - Name of the resulting column. Default value is count.
   */
  count(fieldAlias: string = 'count'): TableQueryBuilder {
    return this.selectAggr(AGG.TOTAL_COUNT, null, fieldAlias);
  }

  /** Adds an aggregation that counts rows for the specified column. Equivalent to count(field).
   * @param field - Column name.
   * @param fieldAlias Name of the resulting column. Default value is count.
   * @returns {GroupByBuilder} */
  valueCount(field: string, fieldAlias: string = 'count'): TableQueryBuilder {
    return this.selectAggr(AGG.VALUE_COUNT, field, fieldAlias);
  }

  /**
   * Adds an aggregation that counts rows with null values
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  nulls(field: string, fieldAlias: string | null = 'nulls'): TableQueryBuilder {
    return this.selectAggr(AGG.MISSING_VALUE_COUNT, field, fieldAlias);
  }

  /**
   * Performs join operation of main table or table specified in {@link leftTable} to table specified in {@link rightTable}.
   * Specify joining fields of the main table (or table specified in {@link leftTable}) in {@link leftTableKeys} and joining fields of {@link rightTable}
   * in {@link rightTableKeys}.
   * @param rightTable {string}
   * @param joinType {JoinType}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param rightTableAlias {string} - use to specify desired alias for the joining table and apply this alias to fields specified in {@link select}.
   * @param leftTable {string}
   */
  join(rightTable: string, joinType: JoinType, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Join(this.dart, rightTable, joinType, leftTableKeys, rightTableKeys, rightTableAlias ?? null, leftTable ?? null));
  }

  /**
   * Performs left join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  leftJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.LEFT, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs right join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  rightJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.RIGHT, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs inner join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  innerJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.INNER, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs outer join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  outerJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.OUTER, leftTableKeys, rightTableKeys, leftTable)
  }

  /** Groups rows that have the same values into summary values
   * @param {string[]} fields - Array of fields to group by
   * @returns {TableQueryBuilder} */
  groupBy(fields: string[]): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_GroupBy(this.dart, fields));
  }

  /** Rotates a table-valued expression by turning the unique values from one column in the expression into multiple
   * columns in the output
   * @param {string[]} fields - Array of fields to pivot on
   * @returns {TableQueryBuilder} */
  pivotOn(fields: string[]): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_PivotOn(this.dart, fields));
  }

  /** Adds a where clause to the query.
   * @param {string} field - Field name. If you join to other tables use correct alias or table name as prefix, e.g. <table alias>.<field>
   * @param {string} pattern - Pattern to test field values against
   * @param columnType Should be provided, if this builder wasn't created from {@link TableInfo} or alias is used. Defaults to {@link TYPE.STRING}
   * @returns {TableQueryBuilder} */
  where(field: string, pattern: string, columnType?: Type): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Where(this.dart, field, pattern, columnType ?? null));
  }

  /** Adds a having clause to the query. Use only with {@link groupBy}
   * @param {string} field - Field name. If you join to other tables use correct alias or table name as prefix, e.g. <table alias>.<field>
   * @param {string} pattern - Pattern to test field values against
   * @param columnType Should be provided, if this builder wasn't created from {@link TableInfo} or alias is used. Defaults to {@link TYPE.STRING}
   * @returns {TableQueryBuilder} */
  having(field: string, pattern: string, columnType?: Type): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Having(this.dart, field, pattern, columnType ?? null));
  }

  /** Sorts results in ascending or descending order
   * @param {string} field - Field to sort based on
   * @param {boolean} asc - Sort in ascending order
   * @returns {TableQueryBuilder} */
  sortBy(field: string, asc: boolean = true): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_SortBy(this.dart, field, asc));
  }

  /** Selects limited number of records
   * @param {TableQueryBuilder} n - Number of records to select
   * @returns {TableQueryBuilder} */
  limit(n: number): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_Limit(this.dart, n)); }

  /** Builds a query
   * @returns {TableQuery} */
  build(): TableQuery { return toJs(api.grok_DbTableQueryBuilder_Build(this.dart)); }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access}
 * */
export class DataJob extends Func {
  /** @constructs DataJob */
  constructor(dart: any) {
    super(dart);
  }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/access#data-connection}
 * */
export class DataConnection extends Entity {
  parameters: any;

  /** @constructs DataConnection */
  constructor(dart: any) {
    super(dart);
    this.parameters = new MapProxy(api.grok_DataConnection_Get_Parameters(this.dart), 'parameters');
  }

  get credentials(): Credentials {
    return toJs(api.grok_DataConnection_Get_Credentials(this.dart));
  }

  get dataSource(): string {
    return api.grok_DataConnection_Get_DataSource(this.dart);
  }

  /** Collection of parameters: server, database, endpoint, etc. */
  // get parameters(): DataConnectionParams { return api.grok_DataConnection_Parameters(this.dart); }

  /** Tests the connection, returns "ok" on success or an error message on error
   * @returns {Promise<string>}*/
  test(): Promise<string> {
    return api.grok_DataConnection_Test(this.dart);
  }

  /**
   * Creates {@link DataQuery} using this connection. Can be used only with database connections.
   * @param name - name of the query
   * @param sql - text of the query
   */
  query(name: string, sql: string): DataQuery {
    return toJs(api.grok_DataConnection_Query(this.dart, name, sql));
  }

  buildQuery(table: string): TableQueryBuilder {
    return TableQueryBuilder.from(table, this);
  }

  tableQuery(): TableQuery {
    return toJs(api.grok_TableQuery_Create(this.dart));
  }

  /** Creates a data connection. Note that in order to be used, it has to be saved first using {@link DataConnectionsDataSource}
   * @param {string} name - Connection name
   * @param {DataConnectionProperties} parameters - Connection properties
   * */
   static create(name: string, parameters: DataConnectionProperties): DataConnection {
    return toJs(api.grok_DataConnection_Create(name, parameters.dataSource, parameters));
  }
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/-info}
 * */
export class Model extends Entity {
  /** @constructs Model */
  constructor(dart: any) {
    super(dart);
  }
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/compute/jupyter-notebook}
 * */
export class Notebook extends Entity {

  /** @constructs Notebook */
  constructor(dart: any) {
    super(dart);
  }

  /** Environment name */
  get environment(): string { return api.grok_Notebook_Get_Environment(this.dart); }
  set environment(e: string) { api.grok_Notebook_Set_Environment(this.dart, e); }

  /** Description */
  get description(): string { return api.grok_Notebook_Get_Description(this.dart); }
  set description(e: string) { api.grok_Notebook_Set_Description(this.dart, e); }

  get notebook(): {[_:string]: any} { return api.grok_Notebook_Get_Notebook_Content(this.dart); }
  set notebook(data: {[_:string]: any}) { api.grok_Notebook_Set_Notebook_Content(this.dart, data); }
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
  /** @constructs TableInfo */
  public tags: {[key: string]: any};

  constructor(dart: any) {
    super(dart);
    this.tags = new MapProxy(api.grok_TableInfo_Get_Tags(this.dart), 'tags');
  }

  static fromDataFrame(t: DataFrame): TableInfo {return toJs(api.grok_DataFrame_Get_TableInfo(t.dart)); }

  get dataFrame(): DataFrame { return toJs(api.grok_TableInfo_Get_DataFrame(this.dart)); }

  get columns(): ColumnInfo[] { return toJs(api.grok_TableInfo_Get_Columns(this.dart)); }
}


/** @extends Entity
 * Represents Column metadata */
export class ColumnInfo extends Entity {
  public tags: {[key: string]: any};
  /** @constructs ColumnInfo */
  constructor(dart: any) {
    super(dart);
    this.tags = new MapProxy(api.grok_ColumnInfo_Get_Tags(this.dart), 'tags');
  }

  get type(): string { return toJs(api.grok_ColumnInfo_Get_Type(this.dart)); }

  /** Returns reference information if the column is referencing another column in another table */
  get referenceInfo(): {table: string, column: string} | null {
    return this.tags[Tags.ReferencesTable] && this.tags[Tags.ReferencesColumn] ? {
      table: this.tags[Tags.ReferencesTable],
      column: this.tags[Tags.ReferencesColumn]
    } : null;
  }
}

/** @extends Entity
 * Allows for files handling in JS-based info panels
 * {@link https://datagrok.ai/help/discover/info-panels}
 * */
export class FileInfo extends Entity {
  /** @constructs FileInfo */
  constructor(dart: any) {
    super(dart);
  }

  get connection(): DataConnection { return api.grok_FileInfo_Get_Connection(toJs(this.dart)); }

  /** Returns path, i.e. `geo/dmv_offices.csv` */
  get path(): string { return api.grok_FileInfo_Get_Path(this.dart); }

  /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv` */
  get fullPath(): string { return api.grok_FileInfo_Get_FullPath(this.dart); }

  /** Returns URL path, i.e. `Demo.TestJobs.Files.DemoFiles/geo/dmv_offices.csv` */
  get viewPath(): string { return api.grok_FileInfo_Get_ViewPath(this.dart); }

  /** Returns file extension, i.e. `csv` */
  get extension(): string { return api.grok_FileInfo_Get_Extension(this.dart); }

  /** Returns file name, i.e. `dmv_offices.csv` */
  get fileName(): string { return api.grok_FileInfo_Get_FileName(this.dart); }

  /** Returns file URL */
  get url(): string { return api.grok_FileInfo_Get_Url(this.dart); }

  /** Checks if file */
  get isFile(): boolean { return api.grok_FileInfo_Get_IsFile(this.dart); }

  /** Checks if directory */
  get isDirectory(): boolean { return api.grok_FileInfo_Get_IsDirectory(this.dart); }

  get updatedOn(): dayjs.Dayjs | null {
    const d = api.grok_FileInfo_Get_UpdatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  /** @returns {Promise<string>} */
  // readAsString(): Promise<string> {
  //   return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsString(this.dart, (x: any) => resolve(x), (x: any) => reject(x)));
  // }

  get data(): Uint8Array {return api.grok_FileInfo_Get_Data(this.dart);}

  readAsString(): Promise<string> {
    return api.grok_FileInfo_ReadAsString(this.dart);
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes(): Promise<Uint8Array> {
    return api.grok_FileInfo_ReadAsBytes(this.dart);
  }

  static fromBytes(path: string, data: Uint8Array): FileInfo {
    if (!path)
      throw new Error('Path can\'t be null or empty');
    return api.grok_FileInfo_FromBytes(path, data);
  }

  static fromString(path: string, data: string): FileInfo {
    if (!path)
      throw new Error('Path can\'t be null or empty');
    return api.grok_FileInfo_FromString(path, data);
  }
}

/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
  /** @constructs Group */
  constructor(dart: any) {
    super(dart);
  }

  static create(name: string): Group { return new Group(api.grok_Group(name)); }

  /** Adds a member to the group
   * @param {Group} m */
  addMember(m: Group): void { api.grok_Group_Add_Member(this.dart, m.dart, false); }

  /** Adds an admin member to the group
   * @param {Group} m */
  addAdminMember(m: Group): void { api.grok_Group_Add_Member(this.dart, m.dart, true); }

  /** Removes a member from the group
   * @param {Group} m */
  removeMember(m: Group): void { api.grok_Group_Remove_Member(this.dart, m.dart); }

  /** Adds the group to another one
   * @param {Group} m */
  includeTo(m: Group): void { api.grok_Group_Add_Membership(this.dart, m.dart, false); }

  /** Adds the group to another one as an admin
   * @param {Group} m */
  includeAdminTo(m: Group): void { api.grok_Group_Add_Membership(this.dart, m.dart, true); }

  /** Removes membership from another group
   * @param {Group} m */
  excludeFrom(m: Group): void { api.grok_Group_Remove_Membership(this.dart, m.dart); }

  /** Returns list of groups that belong to group, with no admin permissions
   * @type {Array<Group>} */
  get members(): Group[] { return toJs(api.grok_Group_Get_Members(this.dart, false)); }

  /** Returns list of groups that belong to group, with admin permissions
   * @type {Array<Group>} */
  get adminMembers(): Group[] { return toJs(api.grok_Group_Get_Members(this.dart, true)); }

  /** Returns list of groups that group belongs to, with no admin permissions
   * @type {Array<Group>} */
  get memberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.dart, false)); }

  /** Returns list of groups that group belongs to, with admin permissions
   * @type {list<Group>} */
  get adminMemberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.dart, true)); }

  /** Personal user group */
  get personal(): boolean { return api.grok_Group_Get_Personal(this.dart); }
  set personal(e: boolean) { api.grok_Group_Set_Personal(this.dart, e); }

  /** Hidden group */
  get hidden(): boolean { return api.grok_Group_Get_Hidden(this.dart); }
  set hidden(e: boolean) { api.grok_Group_Set_Hidden(this.dart, e); }

  /** Returns associated user.
   * Returns `null` if the group is not {@link personal}.
   * See [Groups](https://datagrok.ai/help/govern/access-control/users-and-groups#groups)
   */
  get user(): User { return new User(api.grok_Group_Get_User(this.dart)); }

  static get defaultGroupsIds() {
    return {
      "All users": "a4b45840-9a50-11e6-9cc9-8546b8bf62e6",
      "Developers": "ba9cd191-9a50-11e6-9cc9-910bf827f0ab",
      "Need to create": "00000000-0000-0000-0000-000000000000",
      "Test": "ca1e672e-e3be-40e0-b79b-8546b8bf62e6",
      "Admin": "a4b45840-9a50-11e6-c537-6bf8e9ab02ee",
      "System": "a4b45840-ac9c-4d39-8b4b-4db3e576b3c3",
      "Administrators": "1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5",
    } as const;
  }

  static get allUsers(): Group { return new Group(api.grok_Group_AllUsers()); }

  static get developers(): Group { return new Group(api.grok_Group_Developers()); }

  static get needToCreate(): Group { return new Group(api.grok_Group_NeedToCreate()); }

  static get test(): Group { return new Group(api.grok_Group_Test()); }

  static get admin(): Group { return new Group(api.grok_Group_Admin()); }

  static get system(): Group { return new Group(api.grok_Group_System()); }

  static get administrators(): Group { return new Group(api.grok_Group_Administrators()); }
}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
  public static readonly vecInputTableName = 'in_vec_table';
  public static readonly vecOutputTableName = 'out_vec_table';

  /** @constructs Script */
  constructor(dart: any) {
    super(dart);
  }

  static create(script: string): Script { return new Script(api.grok_Script_Create(script)); }

  static fromParams(inputs: Property[], outputs: Property[], script: string = ''): Script {
    return new Script(api.grok_Script_FromParams(inputs.map((i) => toDart(i)), outputs.map((i) => toDart(i)), script));
  }

  /** Script */
  get script(): string { return api.grok_Script_GetScript(this.dart); }
  set script(s: string) { api.grok_Script_SetScript(this.dart, s); }

  /** Script */
  get clientCode(): string { return api.grok_Script_ClientCode(this.dart); }

  /** Script language. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get language(): ScriptingLanguage { return api.grok_Script_GetLanguage(this.dart); }
  set language(s: ScriptingLanguage) { api.grok_Script_SetLanguage(this.dart, s); }

  /** Indicates that the script is already vector, meaning it accepts vector
   * input (an entire column) and processes it in a single call, rather than
   * being executed separately for each scalar element (row) */
  get isVectorFunc(): boolean { return api.grok_Script_Get_IsVectorFunc(this.dart); }

  /** Environment name. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get environment(): string { return api.grok_Script_Get_Environment(this.dart); }
  set environment(s: string) { api.grok_Script_Set_Environment(this.dart, s); }

  /** Reference header parameter. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get reference(): string { return api.grok_Script_Get_Reference(this.dart); }
  set reference(s: string) { api.grok_Script_Set_Reference(this.dart, s); }

  /** Sample table. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get sample(): string { return api.grok_Script_Get_Sample(this.dart); }
  set sample(s: string) { api.grok_Script_Set_Sample(this.dart, s); }

  /** Script tags. See also: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation */
  get tags(): string[] { return api.grok_Script_Get_Tags(this.dart); }
  set tags(tags: string[]) { api.grok_Script_Set_Tags(this.dart, tags); }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/datagrok/solutions/enterprise/security#credentials}
 *  */
export class Credentials extends Entity {
  /**
   * Represents credentials parameters that are hidden from other users and when they are used real values
   * of them will be replaced by values from {@link openParameters}.
   * Parameters can be filled with key-value entries that are suitable for the entity that owns that Credentials.
   */
  public parameters: any;

  constructor(dart: any) {
    super(dart);
    this.parameters = new MapProxy(api.grok_Credentials_Parameters(this.dart), 'parameters');
  }

  /**
   * Represents opened credential parameters. They will be showed, for example, in the connection edit form.
   * They are updated when {@link parameters} are changed and instance is saved.
   * See also {@link CredentialsDataSource}.
   */
  get openParameters(): Record<string, string> { return api.grok_Credentials_OpenParameters(this.dart); }
}

/** Represents a script environment */
export class ScriptEnvironment extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Create instance of ScriptEnvironment */
  static create(name: string): ScriptEnvironment {
    return new ScriptEnvironment(api.grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content */
  get environment(): string { return api.grok_ScriptEnvironment_Environment(this.dart); }
}

export class LogEventType extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Friendly name of the event type */
  get name(): string { return api.grok_LogEventType_Get_Name(this.dart); }

  get comment(): string { return api.grok_LogEventType_Get_Comment(this.dart); }
  set comment(comment) { api.grok_LogEventType_Set_Comment(this.dart, comment); }

  get isError(): boolean { return api.grok_LogEventType_Get_IsError(this.dart); }
  set isError(isError) { api.grok_LogEventType_Set_IsError(this.dart, isError); }

}

export class LogEvent extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Description of the event */
  get description(): string { return api.grok_LogEvent_Get_Description(this.dart); }

  /** Friendly name of the event */
  get name(): string { return api.grok_LogEvent_Get_Name(this.dart); }

  /** Session id of the event */
  get session(): UserSession | string { return toJs(api.grok_LogEvent_Get_Session(this.dart)); }

  /** Parameters of the event
   * @type {Array<LogEventParameterValue>} */
  get parameters(): LogEventParameterValue[] { return toJs(api.grok_LogEvent_Get_Parameters(this.dart)); }

  /** Type of the event
   * @type {LogEventType} */
  get eventType(): LogEventType { return toJs(api.grok_LogEvent_Get_Type(this.dart)); }

  get eventTime(): dayjs.Dayjs { return dayjs(api.grok_LogEvent_Get_EventTime(this.dart)); }
}

export class LogEventParameter extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Name of the parameter */
  get name(): string { return api.grok_LogEventParameter_Get_Name(this.dart); }

  /** Type of the parameter */
  get type(): string { return api.grok_LogEventParameter_Get_Type(this.dart); }
}

export class LogEventParameterValue extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Event of the parameter value
   * @type {LogEvent} */
  get event(): LogEvent {
    return toJs(api.grok_LogEventParameterValue_Get_Event(this.dart));
  }

  /** Parameter of the parameter value
   * @type {LogEventParameter} */
  get parameter(): LogEventParameter {
    return toJs(api.grok_LogEventParameterValue_Get_Parameter(this.dart));
  }

  /** Parameter value
   * @type {string} */
  get value(): string {
    return api.grok_LogEventParameterValue_Get_Value(this.dart);
  }
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package extends Entity {
  _webRoot: string | undefined;
  public _version: string = '';

  constructor(dart: any | undefined = undefined) {
    super(dart);

    if (typeof dart === 'string') {
      this._webRoot = dart;
      this.dart = null;
    }
  }

  /** Override init() method to provide package-specific initialization.
   * It is guaranteed to get called exactly once before the execution of any function below.
   */
  init(): Promise<null> { return Promise.resolve(null); }

  private _name: string = '';

  get webRoot(): string {
    if (this._webRoot === undefined)
      return api.grok_Package_Get_WebRoot(this.dart);
    else
      return this._webRoot;
  }

  set webRoot(x) {
    this._webRoot = x;
  }

  get packageOwner(): string {
    return api.grok_Package_Get_Package_Author(this.dart);
  }

  get version(): string {
    if (this.dart != null)
      return api.grok_Package_Get_Version(this.dart);
    else
      return this._version;
  }

  set version(x) {
    if (this.dart != null)
      api.grok_Package_Set_Version(this.dart, x);
    else
      this._version = x;
  }

  /** Package short name */
  get name(): string {
    if (this.dart != null)
      return api.grok_Entity_Get_Name(this.dart);
    else
      return this._name;
  }

  set name(x) {
    if (this.dart != null)
      api.grok_Entity_Set_Name(this.dart, x);
    else
      this._name = x;
  }

  getModuleName(file: string): string {
    if (this.dart != null)
      return api.grok_Package_GetModuleName(this.dart, file);
    else
      return '';
  }

  getIconUrl(): string {
    return api.grok_Package_GetIconUrl(this.dart);
  }

  /** Returns a JavaScript module for this package. */
  getModule(file: string): any {
    if (this.dart != null)
      return api.grok_Package_GetModule(this.dart, file)();
    else
      return null;
  }

  /** Returns metadata associated with the package.
   * The metadata gets generated when the package is built.
   * It is a concatenation of JSON files located under the /meta folder.
   * See example: /packages/PowerPack. */
  get meta(): {[key: string]: any} | null {
    return (this.dart == null) ? null : toJs(api.grok_Package_Get_Meta(this.dart));
  }

  /** Loads package. */
  async load(options?: {file: string}): Promise<Package> {
    return api.grok_Dapi_Packages_Load(this.dart, options?.file);
  }

  private _logger?: PackageLogger;

  get logger(): PackageLogger {
    if (this._logger)
      return this._logger;
    return this._logger = new PackageLogger(this);
  }

  /** Returns credentials for package. */
  getCredentials(): Promise<Credentials> {
    return api.grok_Package_Get_Credentials(this.name);
  }

  /**
   * @deprecated The {@link getProperties} should not be used. Use {@link settings} instead
   */
  getProperties(): Promise<any> {
    return this.getSettings();
  }

  /**
   * @deprecated The {@link getSettings} should not be used. Use {@link settings} instead
   */
  getSettings(): Promise<Map<string, any>> {
    return api.grok_Package_Get_Settings(this.name);
  }

  /** Returns settings for a package. */
  get settings(): {[index: string]: any} {
    return api.grok_Package_Get_Settings_Sync(this.name);
  }

  /** Updates settings for a package. */
  setSettings(props: Map<string, any>, group: Group): Promise<void> {
    return api.grok_Package_Set_Settings(this.name, props, group?.dart);
  }

  /** Global application data */
  get files(): FileSource {
    return new DG.FileSource(`System:AppData/${this.name}`);
  }

  public async getTests(core: boolean = false) {
    try {
      await this.load({ file: 'package-test.js' });
      let module = this.getModule('package-test.js');
      if (core && module.initAutoTests)
        await module.initAutoTests();
      return module.tests;
    } catch (e: any) {
      this.logger.error(e?.msg ?? 'get module error')
      return undefined;
    }
  }
}


// export class DockerImage extends Entity {
//   constructor(dart: any) {
//     super(dart);
//   }
// }


export class DockerContainer extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get status(): DockerContainerStatus {
    return api.grok_DockerContainer_Status(this.dart);
  }
}


/** Represents a property.
 * See also {@link Property}. */
export interface IProperty {

  /** Property name */
  name?: string;

  /** Property data type. See {@link TYPE}. */
  type?: string;

  /** Property input type */
  inputType?: string;

  /** Whether an empty value is allowed. This is used by validators. */
  nullable?: boolean;

  /** Property description */
  description?: string;

  /** Semantic type */
  semType?: string;

  /** Units of measurement. See also: [postfix] */
  units?: string;

  /** Minimum value. Applicable to numerical properties only */
  min?: number;

  /** Maximum value. Applicable to numerical properties only */
  max?: number;

  /** Step to be used in a slider. Only applies to numerical properties. */
  step?: number;

  /** Whether a slider appears next to the number input. Applies to numerical columns only. */
  showSlider?: boolean;

  /** Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only. */
  showPlusMinus?: boolean;

  /** List of choices. Applicable to string properties only */
  choices?: string[];

  /** Initial value used when initializing UI. See also {@link defaultValue} */
  initialValue?: any;

  /** Default value used for deserialization and cloning. See also {@link initialValue}. */
  defaultValue?: any;

  /** Custom editor (such as slider or text area) */
  editor?: string;

  /** Corresponding category on the context panel */
  category?: string;

  /** Value format, such as '0.000' */
  format?: string;

  /** Whether the property should be editable via the UI */
  userEditable?: boolean;

  /** List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names.
   * Signature: validator(x: DG.Type): string | null.
   * [null] indicates that the value is valid, [string] describes a validation error. */
  validators?: string[];

  /** List of value validators (functions that take a value and return error message or null) */
  valueValidators?: ValueValidator<any>[];

  /** Custom field caption shown in [PropertyGrid]
   * @deprecated The property will be removed soon. Use {@link friendlyName} instead */
  caption?: string;

  /** Custom field friendly name shown in [PropertyGrid] */
  friendlyName?: string;

  /** Name of the corresponding JavaScript field. No need to specify it if it is the same as name. */
  fieldName?: string;

  tags?: any;

  /** Additional options. */
  options?: any;

  /** Filter for columns, can be numerical, categorical or directly a column type (string, int...)
   * Applicable when type = Column */
  columnTypeFilter?: ColumnType | 'numerical' | 'categorical' | null;


  viewer?: string;
}


/** Properties of properties. See also {@link Property.propertyOptions}. */
interface IPropertyMeta {
  applicableTo?: string;
}


/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property implements IProperty {
  public readonly dart: any;
  public options: any;

  constructor(dart: any) {
    this.dart = dart;
    this.options = new MapProxy(api.grok_Property_Get_Options(this.dart));
  }

  /** Property getter is a function that accepts one parameter (item)
   * and returns the property value.
   *
   * @returns {PropertyGetter} */
  get get(): PropertyGetter { return api.grok_Property_Get_Get(this.dart); }
  set get(x: PropertyGetter) { api.grok_Property_Set_Get(this.dart, x); }

  /** Property setter */
  get set(): PropertySetter { return api.grok_Property_Get_Set(this.dart); }
  set set(x: PropertySetter) { api.grok_Property_Set_Set(this.dart, x); }

  /** Property name */
  get name(): string { return api.grok_Property_Get_Name(this.dart); }
  set name(s: string) { api.grok_Property_Set_Name(this.dart, s); }

  /** Custom field caption shown in the UI
   * @deprecated The property will be removed soon. Use {@link friendlyName} instead */
  get caption(): string { return api.grok_Property_Get_Caption(this.dart); }
  set caption(s: string) { api.grok_Property_Set_Caption(this.dart, s); }

  /** Custom field caption shown in the UI */
  get friendlyName(): string { return api.grok_Property_Get_Caption(this.dart); }
  set friendlyName(s: string) { api.grok_Property_Set_Caption(this.dart, s); }

  /** Property category */
  get category(): string { return api.grok_Property_Get_Category(this.dart); }
  set category(s: string) { api.grok_Property_Set_Category(this.dart, s); }

  /** Property type. Same as {@link propertyType} */
  get type(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set type(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Property type. Same as {@link type} */
  get propertyType(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set propertyType(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Property subtype */
  get propertySubType(): TYPE { return api.grok_Property_Get_PropertySubType(this.dart); }

  /** Applies to viewers properties whether to include the property in the layout or not. */
  get includeInLayout(): boolean { return api.grok_Property_Get_IncludeInLayout(this.dart); }
  set includeInLayout(s: boolean) { api.grok_Property_Set_IncludeInLayout(this.dart, s); }

  /** Semantic type */
  get semType(): SemType | string { return api.grok_Property_Get_SemType(this.dart); }
  set semType(s: SemType | string) { api.grok_Property_Set_SemType(this.dart, s); }

  /** Input type. See also {@link InputType} */
  get inputType(): string { return api.grok_Property_Get_InputType(this.dart); }
  set inputType(s: string) { api.grok_Property_Set_InputType(this.dart, s); }

  /** Description */
  get description(): string { return api.grok_Property_Get_Description(this.dart); }
  set description(s: string) { api.grok_Property_Set_Description(this.dart, s); }

  /** Nullable */
  get nullable(): boolean { return api.grok_Property_Get_Nullable(this.dart); }
  set nullable(s: boolean) { api.grok_Property_Set_Nullable(this.dart, s); }

  /** Initial value used when initializing UI */
  get initialValue(): string { return toJs(api.grok_Property_Get_InitialValue(this.dart)); }
  set initialValue(s: string) { api.grok_Property_Set_InitialValue(this.dart, toDart(s)); }

  /** Default value */
  get defaultValue(): any { return toJs(api.grok_Property_Get_DefaultValue(this.dart)); }
  set defaultValue(s: any) { api.grok_Property_Set_DefaultValue(this.dart, toDart(s)); }

  /** Property editor */
  get editor(): string { return api.grok_Property_Get(this.dart, 'editor'); }
  set editor(s: string) { api.grok_Property_Set(this.dart, 'editor', s); }

  /** Units of measurement. */
  get units(): string { return api.grok_Property_Get_Units(this.dart); }
  set units(s: string) { api.grok_Property_Set_Units(this.dart, s); }

  /** Format to be used for displaying the value. */
  get format(): string { return api.grok_Property_Get_Format(this.dart); }
  set format(s: string) { api.grok_Property_Set_Format(this.dart, s); }

  /** Whether a user can edit this property from the UI. */
  get userEditable(): boolean { return api.grok_Property_Get_UserEditable(this.dart); }
  set userEditable(s: boolean) { api.grok_Property_Set_UserEditable(this.dart, s); }

  /** Minimum value. Used when constructing UI (sliders), validating, etc. */
  get min(): number { return api.grok_Property_Get_Min(this.dart); }
  set min(s: number) { api.grok_Property_Set_Min(this.dart, s); }

  /** Maximum value. Used when constructing UI (sliders), validating, etc. */
  get max(): number { return api.grok_Property_Get_Max(this.dart); }
  set max(s: number) { api.grok_Property_Set_Max(this.dart, s); }

  /** Step to be used in a slider. Only applies to numerical properties. */
  get step(): number { return api.grok_Property_Get_Step(this.dart); }
  set step(s: number) { api.grok_Property_Set_Step(this.dart, s); }

  /** Whether a slider appears next to the number input. Applies to numerical columns only. */
  get showSlider(): boolean { return api.grok_Property_Get_ShowSlider(this.dart); }
  set showSlider(s: boolean) { api.grok_Property_Set_ShowSlider(this.dart, s); }

  /** Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only. */
  get showPlusMinus(): boolean { return api.grok_Property_Get_ShowPlusMinus(this.dart); }
  set showPlusMinus(s: boolean) { api.grok_Property_Set_ShowPlusMinus(this.dart, s); }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns {Array<string>} */
  get choices(): string[] { return api.grok_Property_Get_Choices(this.dart); }
  set choices(x: string[]) { api.grok_Property_Set_Choices(this.dart, x); }

  /** Validation conditions to be checked when editing the property.
   * See also {@link PropertyValidator}. */
  get validators(): string[] { return api.grok_Property_Get_Validators(this.dart); }
  set validators(x: string[]) { api.grok_Property_Set_Validators(this.dart, x); }

  get isVectorizable(): boolean { return api.grok_Property_Get_IsVectorizable(this.dart); }

  get vectorName(): string { return api.grok_Property_Get_VectorName(this.dart); }

  /** Column type filter (previously "columnFilter") */
  get columnTypeFilter(): ColumnType | 'numerical' | 'categorical' | null {
    return api.grok_Property_Get_ColumnTypeFilter(this.dart);
  }

  /** Applies the specified options */
  fromOptions(opt?: IProperty): Property {
    if (opt)
      api.grok_Property_Options(this.dart, opt);
    return this;
  }

  /** Creates a property */
  static create(name: string, type: Type,
                getter: PropertyGetter,
                setter: PropertySetter,
                defaultValue: any = null): Property {
    return new Property(api.grok_Property(name, type, getter, setter, toDart(defaultValue)));
  }

  /** Creates an integer property */
  static int(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.INT, getter, setter, defaultValue);
  }

  /** Creates a float property */
  static float(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.FLOAT, getter, setter, defaultValue);
  }

  /** Creates a string property */
  static string(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.STRING, getter, setter, defaultValue);
  }

  /** Creates a bool property */
  static bool(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.BOOL, getter, setter, defaultValue);
  }

  /** Creates a datetime property */
  static dateTime(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.DATE_TIME, getter, setter, defaultValue);
  }

  /** Creates property for the JavaScript objects with the corresponding property name */
  static js(name: string, type: TYPE, options?: IProperty): Property {
    return Property.create(name, type,
      (x: any) => x[name],
      function (x: any, v: any) { x[name] = v; },
      options?.defaultValue).fromOptions(options);
  }

  static jsInt(name: string, options?: IProperty): Property { return Property.js(name, TYPE.INT, options); }
  static jsBool(name: string, options?: IProperty): Property { return Property.js(name, TYPE.BOOL, options); }
  static jsFloat(name: string, options?: IProperty): Property { return Property.js(name, TYPE.FLOAT, options); }
  static jsString(name: string, options?: IProperty): Property { return Property.js(name, TYPE.STRING, options); }
  static jsDateTime(name: string, options?: IProperty): Property { return Property.js(name, TYPE.DATE_TIME, options); }

  static fromOptions(options: IProperty): Property { return Property.js(options.name!, options.type! as TYPE, options); }

  /** Registers the attached (dynamic) property for the specified type.
   * It is editable via the context panel, and gets saved into the view layout as well.
   * Property getter/setter typically uses Widget's "temp" property for storing the value. */
  static registerAttachedProperty(typeName: string, property: Property) {
    api.grok_Property_RegisterAttachedProperty(typeName, property.dart);
  }

  static propertyOptions:{[name in keyof IProperty]: IProperty & IPropertyMeta } = {
    'name': { name: 'name', type: TYPE.STRING, nullable: false },
    'type': { name: 'type', type: TYPE.STRING, nullable: false, description: 'Property data type, such as "int" or "string".' },
    'inputType': { name: 'inputType', type: TYPE.STRING, friendlyName: 'Input type', description: 'Property input type' },
    'nullable': { name: 'nullable', type: TYPE.BOOL, description: 'Whether an empty value is allowed. This is used by validators.' },
    'description': { name: 'description', type: TYPE.STRING, editor: InputType.TextArea, description: 'Property description' },
    'semType': { name: 'semType', type: TYPE.STRING, friendlyName: 'Semantic type', description: 'Semantic type' },
    'units': { name: 'units', type: TYPE.STRING, description: 'Units of measurement. See also: [postfix]' },
    'min': { name: 'min', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Minimum value. Applicable to numerical properties only' },
    'max': { name: 'max', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Maximum value. Applicable to numerical properties only' },
    'step': { name: 'step', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Step to be used in a slider. Only applies to numerical properties.' },
    'showSlider': { name: 'showSlider', applicableTo: TYPE.NUMERICAL, type: TYPE.BOOL, description: 'Whether a slider appears next to the number input. Applies to numerical columns only.' },
    'showPlusMinus': { name: 'showPlusMinus', applicableTo: TYPE.NUMERICAL, type: TYPE.BOOL, description: 'Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only.' },
    'choices': { name: 'choices', applicableTo: TYPE.STRING, type: TYPE.STRING_LIST, description: 'List of choices. Applicable to string properties only' },
    'initialValue': { name: 'initialValue', type: TYPE.STRING, description: 'Initial value used when initializing UI. See also {@link defaultValue}' },
    'defaultValue': { name: 'defaultValue', type: TYPE.OBJECT, description: 'Default value used for deserialization and cloning. See also {@link initialValue}.' },
    'editor': { name: 'editor', type: TYPE.STRING, description: 'Custom editor (such as slider or text area)' },
    'category': { name: 'category', type: TYPE.STRING, description: 'Corresponding category on the context panel' },
    'format': { name: 'format', type: TYPE.STRING, description: 'Value format, such as "0.00"' },
    'userEditable': { name: 'userEditable', type: TYPE.BOOL, description: 'Whether the property should be editable via the UI' },
    'validators': { name: 'validators', type: TYPE.STRING_LIST, description: 'List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names. Signature: validator(x: DG.Type): string | null. [null] indicates that the value is valid, [string] describes a validation error.' },
    'valueValidators': { name: 'valueValidators', type: TYPE.OBJECT, description: 'List of value validators (functions that take a value and return error message or null)' },
    'friendlyName': { name: 'friendlyName', type: TYPE.STRING, description: 'Custom field friendly name shown in [PropertyGrid]' },
    'fieldName': { name: 'fieldName', type: TYPE.STRING, description: 'Name of the corresponding JavaScript field. No need to specify it if it is the same as name.' },
    'tags': { name: 'tags', type: TYPE.MAP, description: 'Additional tags' },
    'options': { name: 'options', type: TYPE.MAP, description: 'Additional options.' },
  }
}


export class HistoryEntry {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  get object(): object { return toJs(api.grok_HistoryEntry_Get_Object(this.dart)); }
  get time(): object { return toJs(api.grok_HistoryEntry_Get_Time(this.dart)); }
}


export class EntityType {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  static create(name: string, matching: string): EntityType {
    return toJs(api.grok_EntityType_Create(toDart(name), toDart(matching)));
  }

  get name(): string { return toJs(api.grok_EntityType_Get_Name(this.dart)); }
  set name(s: string) { api.grok_EntityType_Set_Name(this.dart, toDart(s)); }

  get matching(): string { return toJs(api.grok_EntityType_Get_Matching(this.dart)); }
  set matching(s: string) { api.grok_EntityType_Set_Matching(this.dart, toDart(s)); }
}


/** A dynamic property associated with the entity. */
export class EntityProperty extends Property {
  constructor(dart: any) {
    super(dart);
  }

  static create(name: string, type: string): EntityProperty {
    return toJs(api.grok_EntityProperty_Create(toDart(name), toDart(type)));
  }
}


/** Represents dynamic property schema, associated with the entity type. */
export class Schema {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(name: string): Schema {
    return toJs(api.grok_Schema_Create(toDart(name)));
  }

  get name(): string { return api.grok_Schema_Get_Name(this.dart); }

  /** Schema properties */
  get properties(): EntityProperty[] { return toJs(api.grok_Schema_Get_Properties(this.dart)); }
  set properties(p: EntityProperty[]) { api.grok_Schema_Set_Properties(this.dart, p); }

  /** Entity types associated with this schema. */
  get entityTypes(): EntityType[] { return toJs(api.grok_Schema_Get_EntityTypes(this.dart)); }
  set entityTypes(et: EntityType[]) { api.grok_Schema_Set_EntityTypes(this.dart, et); }
}


export class UserReport extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get isResolved(): boolean {
    return api.grok_UserReport_IsResolved(this.dart);
  }

  get jiraTicket(): string {
    return api.grok_UserReport_JiraTicket(this.dart);
  }

  get assignee(): User {
    return toJs(api.grok_UserReport_Assignee(this.dart));
  }

  get reporter(): User {
    return toJs(api.grok_UserReport_Reporter(this.dart));
  }

  get description(): string {
    return toJs(api.grok_UserReport_Description(this.dart));
  }

  get createdOn(): dayjs.Dayjs {
    return dayjs(api.grok_UserReport_CreatedOn(this.dart));
  }
}


export class UserReportsRule extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  static async showAddDialog(): Promise<void> {
    await api.grok_ReportsRule_Add_Dialog();
  }
}

export class UserNotification {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  get user(): User {
    return toJs(api.grok_UserNotification_User(this.dart));
  }

  get name(): string {
    return toJs(api.grok_UserNotification_Name(this.dart));
  }

  get friendlyName(): string {
    return toJs(api.grok_UserNotification_FriendlyName(this.dart));
  }

  get text(): string {
    return toJs(api.grok_UserNotification_Text(this.dart));
  }

  get data(): string {
    return toJs(api.grok_UserNotification_Data(this.dart));
  }

  get sender(): string {
    return toJs(api.grok_UserNotification_Sender(this.dart));
  }

  get createdAt(): dayjs.Dayjs {
    return dayjs(api.grok_UserNotification_CreatedAt(this.dart));
  }

  get isRead(): boolean {
    return api.grok_UserNotification_IsRead(this.dart);
  }

  get readAt(): dayjs.Dayjs {
    return dayjs(api.grok_UserNotification_ReadAt(this.dart));
  }
}



export class ViewLayout extends Entity {

  /** @constructs ViewLayout */
  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromJson(json));
  }

  static fromViewState(state: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromViewState(state));
  }

  get viewState(): string {
    return api.grok_ViewLayout_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewLayout_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewLayout_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewLayout_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewLayout_ToJson(this.dart);
  }

  get columns(): ColumnInfo[] {
    return toJs(api.grok_ViewLayout_Get_Columns(this.dart));
  }

}

export class ViewInfo extends Entity {

  /** @constructs ViewInfo */
  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromJson(json));
  }

  static fromViewState(state: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromViewState(state));
  }

  get table() : TableInfo {
    return toJs(api.grok_ViewInfo_Get_Table(this.dart));
  }

  /** Only defined within the context of the OnViewLayoutXXX events */
  get view(): View {
    return toJs(api.grok_ViewInfo_Get_View(this.dart));
  }

  get viewState(): string {
    return api.grok_ViewInfo_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewInfo_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewInfo_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewInfo_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewInfo_ToJson(this.dart);
  }
}


export class ProgressIndicator {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create() {
    return toJs(api.grok_ProgressIndicator_Create());
  }

  get percent(): number {
    return api.grok_ProgressIndicator_Get_Percent(this.dart);
  }

  /** Flag indicating whether the operation was canceled by the user. */
  get canceled(): boolean { return api.grok_ProgressIndicator_Get_Canceled(this.dart); }

  get description(): string { return api.grok_ProgressIndicator_Get_Description(this.dart); }
  set description(s: string) { api.grok_ProgressIndicator_Set_Description(this.dart, s); }

  update(percent: number, description: string): void {
    api.grok_ProgressIndicator_Update(this.dart, percent, description);
  }

  log(line: string): void {
    api.grok_ProgressIndicator_Log(this.dart, line);
  }

  get onProgressUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Updated(this.dart));
  }

  get onLogUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Log_Updated(this.dart));
  }

  get onCanceled(): Observable<any> {
    return observeStream(api.grok_Progress_Canceled(this.dart));
  }
}

export type ViewSpecificSearchProvider = {
  isApplicable?: (s: string) => boolean;
  returnType?: string;
  search: (s: string, view?: ViewBase) => Promise<{priority?: number, results: any} | null>;
  getSuggestions?: (s: string) => {priority?: number, suggestionText: string, suggestionValue?: string}[] | null;
  onValueEnter?: (s: string, view?: ViewBase) => Promise<void>;
  name: string;
  description?: string;
  options?: {relatedViewName?: string, widgetHeight?: number, [key: string]: any}
}

/** the home view should be under the name home and rest should match the views that they are applicable to */
export type SearchProvider = {
  [view: string]: ViewSpecificSearchProvider | ViewSpecificSearchProvider[];
}
