// noinspection JSUnusedGlobalSymbols

import {ColumnType, ScriptingLanguage, SemType, Type, TYPE, USER_STATUS} from "./const";
import { FuncCall } from "./functions";
import {toDart, toJs} from "./wrappers";
import {FileSource} from "./dapi";
import {MapProxy} from "./utils";
import {DataFrame} from "./dataframe";
import {PackageLogger} from "./logger";
import dayjs from "dayjs";
import {IDartApi} from "./api/grok_api.g";
import {DataSourceType} from "./api/grok_shared.api.g";
import { Tags } from "../dg";

declare var grok: any;
const api: IDartApi = <any>window;


type PropertyGetter<TSource = any, TProp = any> = (a: TSource) => TProp;
type PropertySetter<TSource = any, TProp = any> = (a: TSource, value: TProp) => void;
type ValueValidator<T> = (value: T) => string;

/**
 * Represents basic properties of database connection. The list of parameters is not complete. You can find all
 * supported parameters for every connector on {@link https://datagrok.ai/help/access/databases/connectors/}
 */
export interface DatabaseConnectionProperties {
  server?: string;
  port?: number;
  db?: string;
  connString?: string;
}

/** Represents connection cache properties
 *  See also: {@link https://datagrok.ai/help/develop/function_results_cache}
 *  */
export interface DataConnectionCacheProperties {
  cacheResults?: boolean;
  cacheSchema?: boolean;
  cacheInvalidateSchedule?: string;
}

/**
 * Represents data connection properties. It can have variable number of parameters depending on {@link DataSourceType}.
 * For a complete list of supported parameters for the data source type use {@link https://datagrok.ai/help/access/databases/connectors/}
 * or {@link https://datagrok.ai/help/access/files/shares/}
 */
export interface DataConnectionProperties extends DatabaseConnectionProperties, DataConnectionCacheProperties {
  dataSource: string;
  login?: string;
  password?: string;
  accessKey?: string;
  secretKey?: string;
  [x: string]: any;
}

type FieldPredicate = {field: string, pattern: string};
type FieldOrder = {field: string, asc?: boolean};
type GroupAggregation = {aggType: string, colName: string, resultColName?: string, function?: string};
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

  static get defaultUsersIds() {
    return {
      "Test": "ca1e672e-e3be-40e0-b79b-d2c68e68d380",
      "Admin": "878c42b0-9a50-11e6-c537-6bf8e9ab02ee",
      "System": "3e32c5fa-ac9c-4d39-8b4b-4db3e576b3c3",
    } as const;
  }
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

  /** A package this function belongs to. */
  get package(): Package { return api.grok_Func_Get_Package(this.dart); }

  /** Returns {@link FuncCall} object in a stand-by state */
  prepare(parameters: {[name: string]: any} = {}): FuncCall {
    return toJs(api.grok_Func_Prepare(this.dart, parameters));
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
  get where(): FieldPredicate[] { return toJs(api.grok_TableQuery_GetWhereClausesDB(this.dart)); }
  set where(wl: FieldPredicate[]) { api.grok_TableQuery_SetWhereClausesDB(this.dart, wl.map(param => toDart(param))); }

  /** Aggregation clauses {queryPartParams} */
  get aggregations(): GroupAggregation[] { return toJs(api.grok_TableQuery_GetAggregationsDB(this.dart)); }
  set aggregations(wl: GroupAggregation[]) { api.grok_TableQuery_SetAggregationsDB(this.dart, wl.map(param => toDart(param))); }

  /** Having clauses */
  get having(): FieldPredicate[] { return toJs(api.grok_TableQuery_GetHavingDB(this.dart)); }
  set having(wl: FieldPredicate[]) { api.grok_TableQuery_SetHavingDB(this.dart, wl.map(param => toDart(param))); }

  /** Order By clauses */
  get orderBy(): FieldOrder[] { return toJs(api.grok_TableQuery_GetOrderByDB(this.dart)); }
  set orderBy(wl: FieldOrder[]) { api.grok_TableQuery_SetOrderByDB(this.dart, wl.map(param => toDart(param))); }

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
   * @returns {TableQueryBuilder} */
  static from(table: string): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_From(table)); }

  /** Creates {@link TableQueryBuilder} from {@link TableInfo}
   * @param {TableInfo} table - TableInfo object
   * @returns {TableQueryBuilder} */
  static fromTable(table: TableInfo): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_FromTable(table.dart));
  }

  /** Selects all fields of the table
   * @returns {TableQueryBuilder} */
  selectAll(): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_SelectAll(this.dart)); }

  /** Selects specified fields of the table
   * @param {string[]} fields - Array of fields to select
   * @returns {TableQueryBuilder} */
  select(fields: string[]): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_Select(this.dart, fields)); }

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

  /** Adds a where clause to the query
   * @param {string} field - Field name
   * @param {string} pattern - Pattern to test field values against
   * @returns {TableQueryBuilder} */
  where(field: string, pattern: string): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Where(this.dart, field, pattern));
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

  /** Creates a data connection. Note that in order to be used, it has to be saved first using {@link DataConnectionsDataSource}
   * @param {string} name - Connection name
   * @param {DataConnectionProperties} parameters - Connection properties
   * @returns {DataConnection}
   * */
   static create(name: string, parameters: DataConnectionProperties): DataConnection {
    return toJs(api.grok_DataConnection_Create(name, parameters.dataSource, parameters));
  }
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
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

  /** Create Notebook on server for edit.
   * @returns {Promise<string>} Current notebook's name */
  edit(): Promise<string> {
    return api.grok_Notebook_Edit(this.dart);
  }

  /** Environment name */
  get environment(): string { return api.grok_Notebook_Get_Environment(this.dart); }
  set environment(e: string) { api.grok_Notebook_Set_Environment(this.dart, e); }

  /** Description */
  get description(): string { return api.grok_Notebook_Get_Description(this.dart); }
  set description(e: string) { api.grok_Notebook_Set_Description(this.dart, e); }

  /** Converts Notebook to HTML code
   * @returns {Promise<string>} */
  toHtml(): Promise<string> {
    return api.grok_Notebook_ToHtml(this.dart);
  }
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

  /** Script */
  get script(): string { return api.grok_Script_GetScript(this.dart); }
  set script(s: string) { api.grok_Script_SetScript(this.dart, s); }

  /** Script */
  get clientCode(): string { return api.grok_Script_ClientCode(this.dart); }

  /** Script language. See also: https://datagrok.ai/help/datagrok/functions/func-params-annotation */
  get language(): ScriptingLanguage { return api.grok_Script_GetLanguage(this.dart); }
  set language(s: ScriptingLanguage) { api.grok_Script_SetLanguage(this.dart, s); }

  /** Environment name. See also: https://datagrok.ai/help/datagrok/functions/func-params-annotation */
  get environment(): string { return api.grok_Script_Get_Environment(this.dart); }
  set environment(s: string) { api.grok_Script_Set_Environment(this.dart, s); }

  /** Reference header parameter. See also: https://datagrok.ai/help/datagrok/functions/func-params-annotation */
  get reference(): string { return api.grok_Script_Get_Reference(this.dart); }
  set reference(s: string) { api.grok_Script_Set_Reference(this.dart, s); }

  /** Sample table. See also: https://datagrok.ai/help/datagrok/functions/func-params-annotation */
  get sample(): string { return api.grok_Script_Get_Sample(this.dart); }
  set sample(s: string) { api.grok_Script_Set_Sample(this.dart, s); }

  /** Script tags. See also: https://datagrok.ai/help/datagrok/functions/func-params-annotation */
  get tags(): string[] { return api.grok_Script_Get_Tags(this.dart); }
  set tags(tags: string[]) { api.grok_Script_Set_Tags(this.dart, tags); }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
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

  /** Create instance of ScriptEnvironment
   * @param {string} name
   * @returns {ScriptEnvironment}
   * */
  static create(name: string): ScriptEnvironment {
    return new ScriptEnvironment(api.grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content */
  get environment(): string { return api.grok_ScriptEnvironment_Environment(this.dart); }

  /** Setup environment */
  setup(): Promise<void> {
    return api.grok_ScriptEnvironment_Setup(this.dart);
  }
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

  /** Source of the event */
  get source(): string { return api.grok_LogEvent_Get_Source(this.dart); }

  /** Session id of the event */
  get session(): UserSession | string { return toJs(api.grok_LogEvent_Get_Session(this.dart)); }

  /** Parameters of the event
   * @type {Array<LogEventParameterValue>} */
  get parameters(): LogEventParameterValue[] { return api.grok_LogEvent_Get_Parameters(this.dart); }

  /** Type of the event
   * @type {LogEventType} */
  get eventType(): LogEventType { return toJs(api.grok_LogEvent_Get_Type(this.dart)); }
}

export class LogEventParameter extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Name of the parameter */
  get name(): string { return api.grok_LogEventParameter_Get_Name(this.dart); }

  /** Type of the parameter */
  get type(): string { return api.grok_LogEventParameter_Get_Type(this.dart); }

  /** Description of the parameter */
  get description(): string { return api.grok_LogEventParameter_Get_Description(this.dart); }

  /** Is the parameter input */
  get isInput(): boolean { return api.grok_LogEventParameter_Get_IsInput(this.dart); }

  /** Is the parameter optional */
  get isOptional(): boolean { return api.grok_LogEventParameter_Get_IsOptional(this.dart); }
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
  get meta(): {[key: string]: any} {
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
   * Deprecated. Use getSettings instead.
   *  Returns properties for a package.
  */
  getProperties(): Promise<any> {
    return this.getSettings();
  }

  /** Returns settings for a package. */
  getSettings(): Promise<Map<string, any>> {
    return api.grok_Package_Get_Settings(this.name);
  }

  /** Updates settings for a package. */
  setSettings(props: Map<string, any>, group: Group): Promise<void> {
    return api.grok_Package_Set_Settings(this.name, props, group?.dart);
  }

  private _files: FileSource | null = null;

  /** Global application data */
  get files(): FileSource {
    if (this._files == null)
      this._files = new FileSource(`System:AppData/${this.name}`);
    return this._files;
  }

  public async getTests(core: boolean = false) {
    try {
      await this.load({ file: 'package-test.js' });
      let module = this.getModule('package-test.js');
      if (core && module.initAutoTests)
        module.initAutoTests();
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


export interface PropertyOptions {

  /** Property name */
  name?: string;

  /** Property type */
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

  /** List of choices. Applicable to string properties only */
  choices?: string[];

  /** Default value (used for deserialization, cloning, etc) */
  defaultValue?: any;

  /** Custom editor (such as slider or text area) */
  editor?: string;

  /** Corresponding category on the context panel */
  category?: string;

  /** Value format */
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

  /** Field postfix shown in [PropertyGrid]. [units] take precedence over the [postfix] value. */
  postfix?: string;

  /** Name of the corresponding JavaScript field. No need to specify it if it is the same as name. */
  fieldName?: string;

  tags?: any;

  /** Filter for columns, can be numerical, categorical or directly a column type (string, int...) */
  columnTypeFilter?: ColumnType | 'numerical' | 'categorical';
}


/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property {
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

  get caption(): string { return api.grok_Property_Get_Caption(this.dart); }
  set caption(s: string) { api.grok_Property_Set_Caption(this.dart, s); }

  /** Property category */
  get category(): string { return api.grok_Property_Get_Category(this.dart); }
  set category(s: string) { api.grok_Property_Set_Category(this.dart, s); }

  /** Property type */
  get propertyType(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set propertyType(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Applies to viewers properties whether to include the property in the layout or not. */
  get includeInLayout(): boolean { return api.grok_Property_Get_IncludeInLayout(this.dart); }
  set includeInLayout(s: boolean) { api.grok_Property_Set_IncludeInLayout(this.dart, s); }

  /** Semantic type */
  get semType(): SemType | string { return api.grok_Property_Get_SemType(this.dart); }
  set semType(s: SemType | string) { api.grok_Property_Set_SemType(this.dart, s); }

  /** Description */
  get description(): string { return api.grok_Property_Get_Description(this.dart); }
  set description(s: string) { api.grok_Property_Set_Description(this.dart, s); }

  /** Nullable */
  get nullable(): boolean { return api.grok_Property_Get_Nullable(this.dart); }
  set nullable(s: boolean) { api.grok_Property_Set_Nullable(this.dart, s); }

  /** Default value */
  get defaultValue(): any { return api.grok_Property_Get_DefaultValue(this.dart); }
  set defaultValue(s: any) { api.grok_Property_Set_DefaultValue(this.dart, s); }

  /** Property editor */
  get editor(): string { return api.grok_Property_Get(this.dart, 'editor'); }
  set editor(s: string) { api.grok_Property_Set(this.dart, 'editor', s); }

  /** Whether a slider appears next to the number input. Applies to numerical columns only. */
  get showSlider(): string { return api.grok_Property_Get_ShowSlider(this.dart); }
  set showSlider(s: string) { api.grok_Property_Set_ShowSlider(this.dart, s); }

  /** Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only. */
  get showPlusMinus(): string { return api.grok_Property_Get_ShowPlusMinus(this.dart); }
  set showPlusMinus(s: string) { api.grok_Property_Set_ShowPlusMinus(this.dart, s); }

  get format(): string { return api.grok_Property_Get_Format(this.dart); }
  set format(s: string) { api.grok_Property_Set_Format(this.dart, s); }

  get userEditable(): boolean { return api.grok_Property_Get_UserEditable(this.dart); }
  set userEditable(s: boolean) { api.grok_Property_Set_UserEditable(this.dart, s); }

  get min(): number { return api.grok_Property_Get_Min(this.dart); }
  set min(s: number) { api.grok_Property_Set_Min(this.dart, s); }

  get max(): number { return api.grok_Property_Get_Max(this.dart); }
  set max(s: number) { api.grok_Property_Set_Max(this.dart, s); }

  get step(): number { return api.grok_Property_Get_Step(this.dart); }
  set step(s: number) { api.grok_Property_Set_Step(this.dart, s); }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns {Array<string>} */
  get choices(): string[] { return api.grok_Property_Get_Choices(this.dart); }
  set choices(x: string[]) { api.grok_Property_Set_Choices(this.dart, x); }

  get isVectorizable(): boolean { return api.grok_Property_Get_IsVectorizable(this.dart); }

  get vectorName(): string { return api.grok_Property_Get_VectorName(this.dart); }

  /** Column type filter */
  get columnFilter(): ColumnType | 'numerical' | 'categorical' | null {
    return api.grok_Property_Get_ColumnTypeFilter(this.dart);
  }

  /** Applies the specified options */
  fromOptions(opt?: PropertyOptions): Property {
    if (opt)
      api.grok_Property_Options(this.dart, opt);
    return this;
  }

  /** Creates a property */
  static create(name: string, type: Type,
                getter: PropertyGetter,
                setter: PropertySetter,
                defaultValue: any = null): Property {
    return new Property(api.grok_Property(name, type, getter, setter, defaultValue));
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
  static js(name: string, type: TYPE, options?: PropertyOptions): Property {
    return Property.create(name, type,
      (x: any) => x[name],
      function (x: any, v: any) { x[name] = v; },
      options?.defaultValue).fromOptions(options);
  }

  static jsInt(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.INT, options); }
  static jsBool(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.BOOL, options); }
  static jsFloat(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.FLOAT, options); }
  static jsString(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.STRING, options); }
  static jsDateTime(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.DATE_TIME, options); }

  static fromOptions(options: PropertyOptions): Property { return Property.js(options.name!, options.type! as TYPE, options); }

  /** Registers the attached (dynamic) property for the specified type.
   * It is editable via the context panel, and gets saved into the view layout as well.
   * Property getter/setter typically uses Widget's "temp" property for storing the value. */
  static registerAttachedProperty(typeName: string, property: Property) {
    api.grok_Property_RegisterAttachedProperty(typeName, property.dart);
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

export class EntityProperty extends Property {
  constructor(dart: any) {
    super(dart);
  };

  static create(name: string, type: string): EntityProperty {
    return toJs(api.grok_EntityProperty_Create(toDart(name), toDart(type)));
  }
}

export class Schema {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  static create(name: string): Schema {
    return toJs(api.grok_Schema_Create(toDart(name)));
  }

  get name(): string { return api.grok_Schema_Get_Name(this.dart); }
  get properties(): EntityProperty[] { return toJs(api.grok_Schema_Get_Properties(this.dart)); }
  set properties(p: EntityProperty[]) { api.grok_Schema_Set_Properties(this.dart, p); }
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
