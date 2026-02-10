import { DataFrame } from "./dataframe";
import {
  Credentials,
  DataConnection,
  DataJob,
  DataQuery,
  DockerContainer,
  Entity,
  Group,
  Model,
  Notebook,
  Project,
  Script,
  ScriptEnvironment,
  TableInfo,
  User,
  LogEvent,
  LogEventType,
  Package,
  UserSession,
  Property,
  FileInfo, ProjectOpenOptions, Func, UserReport, UserReportsRule, ViewLayout, ViewInfo, UserNotification,
} from './entities';
import { DockerImage } from "./api/grok_shared.api.g";
import {toJs, toDart} from "./wrappers";
import {_propsToDart} from "./utils_convert";
import {FuncCall} from "./functions";
import {IDartApi} from "./api/grok_api.g";
import { StickyMeta } from "./sticky_meta";
import {CsvImportOptions} from "./const";
import dayjs from 'dayjs';
import {DbInfo} from "./data";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export class ComponentBuildInfo {
  branch: string = '';
  commit: string = '';
  date: string = '';
  version: string = '';
}

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
 * */
export class Dapi {
  constructor() {
  }

  /** HTTP root for DAPI */
  get root(): string { return api.grok_Dapi_Root();  }
  set root(root: string) { api.grok_Dapi_Set_Root(root); }

  /** Session token */
  get token(): string { return api.grok_Dapi_Get_Token(); }
  set token(token: string | undefined) { api.grok_Dapi_Set_Token(token); }

  /** Retrieves entities from server by list of IDs */
  getEntities(ids: string[]): Promise<Entity[]> {
    return api.grok_Dapi_Entities_GetEntities(ids);
  }

  /** Entities API */
  get entities(): EntitiesDataSource {
    return new EntitiesDataSource(api.grok_Dapi_Entities());
  }

  /** Data Queries API */
  get queries(): HttpDataSource<DataQuery> {
    return new HttpDataSource(api.grok_Dapi_Queries(), 'Function');
  }

  /** Functions API (finding functions and historical function calls) */
  get functions(): FuncsDataSource {
    return new FuncsDataSource(api.grok_Dapi_Functions());
  }

  /** Data Connections API (finding, saving, sharing data connections 
   * and folders, getting database schemas) */
  get connections(): DataConnectionsDataSource {
    return new DataConnectionsDataSource(api.grok_Dapi_Connections());
  }

  /** Credentials API (search, saving, deleting credentials)  */
  get credentials(): CredentialsDataSource {
    return new CredentialsDataSource(api.grok_Dapi_Credentials());
  }

  /** Data Jobs API */
  get jobs(): HttpDataSource<DataJob> {
    return new HttpDataSource(api.grok_Dapi_Jobs(), 'Function');
  }

  /** Jupyter Notebooks API */
  get notebooks(): HttpDataSource<Notebook> {
    return new HttpDataSource(api.grok_Dapi_Notebooks());
  }

  /** Predictive Models API endpoint */
  get models(): HttpDataSource<Model> {
    return new HttpDataSource(api.grok_Dapi_Models());
  }

  /** Packages API endpoint */
  get packages(): HttpDataSource<Package> {
    return new HttpDataSource(api.grok_Dapi_Packages());
  }

  /** View Layouts API (getting applicable view layouts for a dataframe) */
  get layouts(): LayoutsDataSource {
    return new LayoutsDataSource(api.grok_Dapi_Layouts());
  }

  /** View Views API endpoint */
  get views(): ViewsDataSource {
    return new ViewsDataSource(api.grok_Dapi_Views());
  }

  /** Data Table Infos API (finding, uploadeing, deleting tables) */
  get tables(): TablesDataSource {
    return new TablesDataSource(api.grok_Dapi_Tables());
  }

  /** Users API (finding, saving, deleting users, notifications).
   * Also, current user and current session. */
  get users(): UsersDataSource {
    return new UsersDataSource(api.grok_Dapi_Users());
  }

  /** Groups API (finding, saving, deleting groups, adding members, adding admins, including to parent groups, excluding from parent groups) */
  get groups(): GroupsDataSource {
    return new GroupsDataSource(api.grok_Dapi_Groups(), 'Group');
  }

  /** Permissions API (checking, granting, revoking permissions on entities) */
  get permissions(): PermissionsDataSource {
    return new PermissionsDataSource();
  }

  /** Scripts API (finding, saving, deleting scripts) */
  get scripts(): HttpDataSource<Script> {
    return new HttpDataSource(api.grok_Dapi_Scripts(), 'Function');
  }

  /** Projects API endpoint
   *  @type {HttpDataSource<Project>} */
  get projects(): ProjectsDataSource {
    return new ProjectsDataSource(api.grok_Dapi_Projects(), 'Project');
  }

  /** Environments API endpoint
   *  @type {HttpDataSource<ScriptEnvironment>} */
  get environments(): HttpDataSource<ScriptEnvironment> {
    return new HttpDataSource(api.grok_Dapi_Environments());
  }

  get docker(): DockerDataSource {
    return new DockerDataSource();
  }

  get spaces(): SpacesDataSource {
    return new SpacesDataSource(api.grok_Dapi_Spaces());
  }

  /**
   * @deprecated The UserDataStorage should not be used. Use {@link UserSettingsStorage} instead
   */
  get userDataStorage(): UserDataStorage {
    return new UserDataStorage();
  }


  /** Users Files management API endpoint
   *  @type {FilesDataSource} */
  get files(): FilesDataSource {
    return new FilesDataSource();
  }

  get reports(): UserReportsDataSource {
    return new UserReportsDataSource(api.grok_Dapi_User_Reports());
  }

  get rules(): HttpDataSource<UserReportsRule> {
    return new UserReportsRulesDataSource(api.grok_Dapi_User_Reports_Rules());
  }

  /** Proxies URL request via Datagrok server with same interface as "fetch".
   * Useful for cicrumventing CORS restrictions, and for caching results.
   * @see [sample](../../../../../packages/ApiSamples/scripts/dapi/fetch.js)
   * @param {number} maxAge - forces server to send Cache-Control in response with configured max-age directive */
  async fetchProxy(url: string, params?: RequestInit, maxAge?: number): Promise<Response> {
    params ??= {};
    params.headers ??= {};
    params.method ??= 'GET';
    // @ts-ignore
    params.headers['original-url'] = `${url}`;
    // @ts-ignore
    params.headers['original-method'] = params.method;
    if (params.redirect === 'follow')
      (params.headers as any)['follow-redirects'] = true;
    let proxyUrl = `${this.root}/connectors/proxy`;
    if (params.method == 'GET' || params.method == 'HEAD') {
      if (maxAge) {
        // @ts-ignore
        params.headers['dg-cache-control'] = `max-age=${maxAge}`;
        params.cache = 'default';
      }
      proxyUrl = `${proxyUrl}?url=${encodeURI(url)}`;
    }
    if (params.method !== 'GET')
      params.method = 'POST';
    return fetch(proxyUrl, params);
  }

  /** Administering API endpoint
   *  @type {AdminDataSource} */
  get admin(): AdminDataSource {
    return new AdminDataSource(api.grok_Dapi_Admin());
  }

  /** Logging API endpoint
   *  @type {HttpDataSource<LogEvent>} */
  get log(): LogDataSource {
    return new LogDataSource(api.grok_Dapi_Log());
  }

  /** Logging API endpoint
   *  @type {HttpDataSource<LogEventType>} */
  get logTypes(): HttpDataSource<LogEventType> {
    return new HttpDataSource(api.grok_Dapi_LogTypes());
  }

  stickyMeta = new StickyMeta();
}


/**
 * Common functionality for handling collections of entities stored on the server.
 * Works with Datagrok REST API, allows to get filtered and paginated lists of entities,
 * Can be extended with specific methods. (i.e. {@link UsersDataSource})
 */
export class HttpDataSource<T> {
  dart: any;
  clsName: string;

  /** @constructs HttpDataSource */
  constructor(s: any, clsName?: string | null) {
    this.dart = s;
    this.clsName = clsName ?? '';
  }

  /** Returns all entities that satisfy the filtering criteria (see {@link filter}).
   *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  Smart filter: {@link https://datagrok.ai/help/datagrok/smart-search} */
  list(options: {pageSize?: number, pageNumber?: number, filter?: string, order?: string} = {}): Promise<T[]> {
    if (options.pageSize !== undefined)
      this.by(options.pageSize);
    if (options.pageNumber !== undefined)
      this.page(options.pageNumber);
    if (options.filter !== undefined)
      this.filter(options.filter);
    if (options.order !== undefined)
      this.order(options.order);
    return api.grok_DataSource_List(this.dart);
  }

  /** Counts entities that satisfy the filtering criteria (see {@link filter}).
   *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  Smart filter: {@link https://datagrok.ai/help/datagrok/smart-search} */
  count(): Promise<number> {
    return api.grok_DataSource_Count(this.dart);
  }

  /** Returns fist entity that satisfies the filtering criteria (see {@link filter}). */
  first(): Promise<T> {
    return api.grok_DataSource_First(this.dart);
  }

  /** Returns an entity with the specified id.
   *  Throws an exception if an entity does not exist, or is not accessible in the current context.
   *  Sample: {@link https://public.datagrok.ai/js/samples/data-access/save-and-load-df}
   *  @param {string} id - GUID of the corresponding object
   *  @returns `{Promise<object>}` - entity. */
  find(id: string): Promise<T> {
    return api.grok_DataSource_Find(this.dart, id);
  }

  /** Saves an entity. */
  save(e: Entity): Promise<T> {
    return api.grok_DataSource_Save(this.dart, e.dart);
  }

  /** Deletes an entity. */
  delete(e: Entity): Promise<void> {
    return api.grok_DataSource_Delete(this.dart, e.dart);
  }

  /** Turns off package versions isolation. This DataSource will return all entities in all versions, not only the current one **/
  allPackageVersions(): HttpDataSource<T> {
    this.dart = api.grok_DataSource_AllPackageVersions(this.dart);
    return this;
  }

  by(i: number): HttpDataSource<T> {
    this.dart = api.grok_DataSource_By(this.dart, i);
    return this;
  }

  /** Restricts results to the specified page number. See also {@link nextPage}. */
  page(i: number): HttpDataSource<T> {
    this.dart = api.grok_DataSource_Page(this.dart, i);
    return this;
  }

  /** Returns next page of all entities that satisfy the filtering criteria (see {@link filter}).
   *  Works only if pageSize was set during previous list() call
   *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list} */
  nextPage(): HttpDataSource<T> {
    this.dart = api.grok_DataSource_NextPage(this.dart);
    return this;
  }

  /** Applies filter to current request.
   *  Also can be set with {@link list} method "options" parameter
   *  See example: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  Smart filter: {@link https://datagrok.ai/help/datagrok/navigation/views/browse#entity-search} */
  filter(w: string): HttpDataSource<T> {
    this.dart = api.grok_DataSource_WhereSmart(this.dart, w);
    return this;
  }

  /** Instructs data source to return results in the specified order. */
  order(fieldName: string, desc: boolean = false): HttpDataSource<T> {
    this.dart = api.grok_DataSource_Order(this.dart, fieldName, desc);
    return this;
  }

  /** Includes entity in the result */
  include(include: string): HttpDataSource<T> {
    this.dart = api.grok_DataSource_Include(this.dart, _propsToDart(include, this.clsName));
    return this;
  }
}


/**
 * Functionality for handling Users collection from server and working with Users remote endpoint
 * Allows to load current user and list of all Datagrok users with filtering and pagination
 * See example: {@link https://public.datagrok.ai/js/samples/dapi/who-am-i}
 * @extends HttpDataSource
 * */
export class UsersDataSource extends HttpDataSource<User> {

  constructor(s: any) {
    super(s);
  }

  /** Notifications API endpoint */
  get notifications(): NotificationsDataSource {
    return new NotificationsDataSource(api.grok_Dapi_Notifications());
  }

  /** Returns current user */
  current(): Promise<User> {
    return toJs(api.grok_UsersDataSource_Current(this.dart));
    //return new Promise((resolve, reject) => api.grok_UsersDataSource_Current(this.dart, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

  /** Returns current session */
  currentSession(): Promise<UserSession> {
    return api.grok_UsersDataSource_CurrentSession(this.dart);
  }
}

export class AdminDataSource {
  dart: any;
  /** @constructs AdminDataSource*/
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Returns information about the services.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/admin} */
  getServiceInfos(): Promise<ServiceInfo[]> {
    return api.grok_Dapi_Admin_GetServiceInfos(this.dart);
  }

  /**
   * Sends email
   * @param email - message that will be sent using configured SMTP service
   */
  sendEmail(email: Email): Promise<void> {
    if (email.to.length === 0)
      throw new Error('Recipients list shouldn\'t be empty');
    return api.grok_Dapi_Admin_Send_Email(this.dart, toDart(email));
  }

  /** Sends a message to the specified browser client
   * See also {@link ServerMessageTypes} */
  pushMessage(messageType: string, message: object, sessionIds: string[]): Promise<any> {
    return api.grok_Dapi_Admin_PushMessage(this.dart, messageType, api.grok_JSON_decode(JSON.stringify(message)), toDart(sessionIds));
  }
}

/** Email that can be sent using the configured SMTP service */
export interface Email {
  /** Message subject */
  subject: string,

  /** List of recipients */
  to: string [],

  /** Plaintext body */
  text?: string,

  /** HTML body (takes precedence over plain text) */
  html?: string,

  /** Blind carbon copy */
  bcc?: string [],
}

export interface ServiceInfo {
  name: string,
  description: string,
  enabled: boolean,
  key: string,
  time: string,
  status: 'Running' | 'Failed' | 'Stopped',
  type: 'Service' | 'Plugin',
}

/**
 * Functionality for handling groups collection from server
 * Allows to manage {@link Group}
 * @extends HttpDataSource
 * */
export class GroupsDataSource extends HttpDataSource<Group> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any, clsName: string) {
    super(s, clsName);
    this.include('members');
    this.include('memberships');
  }

  /** Creates a new group
   *  @returns {Promise<Group>} - Group. */
  createNew(name: string): Promise<Group> {
    let g = new Group(api.grok_Group(name));
    return this.save(g);
  }

  /** Returns group user
   *  @returns {Promise<Group>} - Group. */
  getUser(group: Group): Promise<Group> {
    return api.grok_Dapi_Get_GroupUser(group.dart);
  }

  /** Adds a member to the group */
  async addMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.addMember(m);
    await this.saveRelations(g);
  }

  /** Adds an admin member to the group */
  async addAdminMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.addAdminMember(m);
    await this.saveRelations(g);
  }

  /** Removes a member from the group */
  async removeMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.removeMember(m);
    await this.saveRelations(g);
  }

  /** Adds the group to another one */
  async includeTo(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.includeTo(parent);
    await this.saveRelations(g);
  }


  /** Adds the group to another one as admin */
  async includeAdminTo(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.includeAdminTo(parent);
    await this.saveRelations(g);
  }

  /** Removes a membership from the group */
  async excludeFrom(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.excludeFrom(parent);
    await this.saveRelations(g);
  }

  /** Saves a group with relations
   *  @returns {Promise<Group>} - Group. */
  saveRelations(e: Group): Promise<Group> {
    return api.grok_GroupsDataSource_Save(this.dart, e.dart);
  }

  /** Looking for groups with similar name */
  async getGroupsLookup(name: string): Promise<Group[]> {
    return toJs(await api.grok_Dapi_Get_GroupsLookup(name));
  }
}


/**
 * Functionality for handling entities collection from server
 * Allows to manage {@link Entity}
 * @extends HttpDataSource
 * */
export class EntitiesDataSource extends HttpDataSource<Entity> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any) {
    super(s);
  }

  /** Returns recent entities */
  getRecentEntities(): Promise<Entity[]> {
    return toJs(api.grok_EntitiesDataSource_GetRecent(this.dart));
  }

  /** Allows to set properties for entities */
  saveProperties(props: Map<Property, any>): Promise<void> {
    return api.grok_EntitiesDataSource_SaveProperties(this.dart, props);
  }

  /** Returns entity properties
   * @returns {Promise<Map>} props */
  getProperties(entity: Entity): Promise<Map<Property, any>> {
    return api.grok_EntitiesDataSource_GetProperties(this.dart, entity.dart);
  }

  /** Deletes entity properties */
  deleteProperties(props: Map<Property, any>): Promise<void> {
    return api.grok_EntitiesDataSource_DeleteProperties(this.dart, props);
  }
}

/**
 * Functionality for handling connections collection from server and working with credentials remote endpoint
 * Allows to manage {@link DataConnection}
 * See also: {@link https://datagrok.ai/help/datagrok/solutions/enterprise/security}
 * @extends HttpDataSource
 * */
export class DataConnectionsDataSource extends HttpDataSource<DataConnection> {
  /** @constructs DataConnectionsDataSource*/
  constructor(s: any) {
    super(s);
  }

  /** Saves the Connections */
  async save(e: DataConnection, options?: {saveCredentials?: boolean}): Promise<DataConnection> {
    options ??= {};
    options.saveCredentials ??= true;
    return toJs(await api.grok_DataConnectionsDataSource_Save(this.dart, e.dart, options!.saveCredentials));
  }

  /** Creates connection to the subdirectory of connection */
  async shareFolder(e: DataConnection, path: string): Promise<DataConnection> {
    return toJs(await api.grok_DataConnectionsDataSource_SubDir(this.dart, e.dart, path));
  }

  async getSchemas(e: DataConnection): Promise<string[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Schemas(this.dart, e.dart));
  }

  async getSchema(e: DataConnection, schemaName: string | null = null, tableName: string | null = null): Promise<TableInfo[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Schema(this.dart, e.dart, schemaName ?? null, tableName ?? null));
  }

  async getUniqueColumnsNames(c: DataConnection, schema: string, table: string): Promise<string[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Unique_Columns(this.dart, c.dart, schema, table));
  }


  async getDatabaseInfo(c: DataConnection, catalog: string | null = null): Promise<DbInfo[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Db_Info(this.dart, c.dart, catalog))
  }
}

/**
 * Functionality for handling functions collection from server
 * Allows managing {@link Func}
 * @extends HttpDataSource
 * */
export class FuncsDataSource extends HttpDataSource<Func> {
  /** @constructs DataConnectionsDataSource*/
  constructor(s: any) {
    super(s);
  }

  public get calls() {
    return new HttpDataSource<FuncCall>(api.grok_Dapi_Function_Calls());
  }
}

/**
 * Functionality for handling credentials collection from server and working with credentials remote endpoint
 * Allows to manage {@link Credentials}
 * See also: {@link https://datagrok.ai/help/datagrok/solutions/enterprise/security#credentials}
 * @extends HttpDataSource
 * */
export class CredentialsDataSource extends HttpDataSource<Credentials> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any) {
    super(s);
  }

  /** Returns credentials for entity */
  forEntity(e: Entity): Promise<Credentials> {
    return api.grok_CredentialsDataSource_ForEntity(this.dart, e.dart);
  }

  /** Saves a credentials.
   * Note, that in order to work correct, credentials should be connected
   * to other entity that owns them. So the best way to modify Credentials is load by {@link forEntity}, change
   * {@link Credentials.parameters} and after that call this method.
   */
  save(c: Credentials): Promise<Credentials> {
    return api.grok_CredentialsDataSource_Save(this.dart, c.dart);
  }
}

/**
 * Functionality for handling layouts collection from server
 * Allows to manage {@link ViewLayout}
 * @extends HttpDataSource
 * */
export class LayoutsDataSource extends HttpDataSource<ViewLayout> {
  /** @constructs LayoutsDataSource*/
  constructor(s: any) {
    super(s);
  }

  /** Returns layouts that applicable to the table */
  getApplicable(t: DataFrame): Promise<ViewLayout[]> {
    return api.grok_LayoutsDataSource_Applicable(this.dart, t.dart);
  }
}

/**
 * Functionality for handling views information from server
 * Allows to manage {@link ViewInfo}
 * @extends HttpDataSource
 * */
export class ViewsDataSource extends HttpDataSource<ViewInfo> {
  /** @constructs ViewsDataSource*/
  constructor(s: any) {
    super(s);
  }
}

export class PermissionsDataSource {
  constructor() {
  };

  /** Gets all the permissions granted on entity
   * @returns {Promise<Map>} permissions
   * */
  // { [key:string]:number; }
  async get(e: Entity): Promise<Map<string, Group[]>> {
    let data = await api.grok_Dapi_Get_Permissions(e.dart);
    data.view = toJs(data.view);
    data.edit = toJs(data.edit);
    return data;
  }

  /** Checks if current user has permission {permission} for entity {e}
   * @param {Entity} e Entity to check permission for
   * @param {'Edit' | 'View' | 'Share' | 'Delete'} permission Permission type
   * @returns {boolean} Result */
  check(e: Entity, permission: 'Edit' | 'View' | 'Share' | 'Delete'): Promise<boolean> {
    return api.grok_Dapi_Check_Permissions(e.dart, permission);
  }

  /** Grants permission on entity to the group
   * @param {boolean} edit allow to edit entity
   * */
  grant(e: Entity, g: Group, edit: boolean): Promise<any> {
    return api.grok_Dapi_Set_Permission(e.dart, g.dart, edit);
  }

  /** Revokes permission on entity from the group
   * */
  revoke(g: Group, e: Entity): Promise<any> {
    return api.grok_Dapi_Delete_Permission(e.dart, g.dart);
  }
}

/**
 * @deprecated The UserDataStorage should not be used. Use {@link UserSettingsStorage} instead
 * Functionality for working with remote Users Data Storage
 * Remote storage allows to save key-value pairs on the Datagrok server for further use
 * */
export class UserDataStorage {
  constructor() {
  }

  /** Saves a single value to Users Data Storage
   * @param {string} name Storage name
   * @param {boolean} currentUser Value should be available only for current user. If false, shared storage is used. */
  postValue(name: string, key: string, value: string, currentUser: boolean = true): Promise<void> {
    return api.grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser);
  }

  /** Saves a map to Users Data Storage, will be appended to existing data
   * @param {boolean} currentUser Value should be available only for current user. If false, shared storage is used. */
  post(name: string, data: any, currentUser: boolean = true): Promise<void> {
    return api.grok_Dapi_UserDataStorage_Post(name, data, currentUser);
  }

  /** Saves a map to Users Data Storage, will replace existing data
   * @param {boolean} currentUser Value should be available only for current user. If false, shared storage is used. */
  put(name: string, data: any, currentUser: boolean = true): Promise<void> {
    return api.grok_Dapi_UserDataStorage_Put(name, data, currentUser);
  }

  /** Retrieves a map from Users Data Storage
   * @param {boolean} currentUser - get a value from a current user storage. If false, shared storage is used.
   * @returns {Promise<Map>} */
  get(name: string, currentUser: boolean = true): Promise<any> {
    return api.grok_Dapi_UserDataStorage_Get(name, currentUser);
  }

  /** Retrieves a single value from Users Data Storage
   * @param {string} name Storage name
   * @param {string} key Value key
   * @param {boolean} currentUser get a value from a current user storage. If false, shared storage is used.
   * @returns {Promise<string>} */
  getValue(name: string, key: string, currentUser: boolean = true): Promise<string> {
    return api.grok_Dapi_UserDataStorage_GetValue(name, key, currentUser);
  }

  /** Removes a single value from Users Data Storage
   * @param {string} name Storage name
   * @param {string} key Value key
   * @param {boolean} currentUser get a value from a current user storage. If false, shared storage is used.
   * @returns {Promise} */
  remove(name: string, key: string, currentUser: boolean = true): Promise<void> {
    return api.grok_Dapi_UserDataStorage_Delete(name, key, currentUser);
  }
}


/**
 * Functionality for working with remote projects
 * @extends HttpDataSource
 * */
export class ProjectsDataSource extends HttpDataSource<Project> {
  /** @constructs TablesDataSource*/
  constructor(s: any, clsName: string) {
    super(s, clsName);
    this.include('children');
  }

  /** Gets recent projects datasource */
  get recent(): HttpDataSource<Project> {
     return new HttpDataSource<Project>(api.grok_Dapi_RecentProjects());
  }

  /** Opens the specified project. */
  open(name: string, options?: ProjectOpenOptions): Promise<Project> {
    return this
      .filter(name)
      .first()
      .then(p => {
        if (p)
          return p.open(options);
        return Promise.reject(`Project ${name} not found`);
      });
  }

  /** Saves the Project */
  save(e: Entity, options?: {saveRelations?: boolean}): Promise<Project> {
    options ??= {};
    options.saveRelations ??= true;
    return toJs(api.grok_ProjectsDataSource_Save(this.dart, e.dart, options!.saveRelations));
  }
}

/**
 * Data source for working with Spaces - hierarchical containers for organizing entities and files.
 * Spaces support nested subspaces, file storage, and can contain various entities like scripts, queries, etc. */
export class SpacesDataSource extends HttpDataSource<Project> {
  constructor(s: any) {
    super(s, 'Project');
  }

  createRootSpace(name: string): Promise<Project> {
    return toJs(api.grok_Dapi_Spaces_CreateRootSpace(this.dart, name));
  }

  /**
   * Checks if a root space (top-level space) with the given name already exists.
   * Use this before creating a new root space to avoid naming conflicts. */
  rootSpaceExists(name: string): Promise<boolean> {
    return api.grok_Dapi_Spaces_RootSpaceExists(this.dart, name);
  }

  /**
   * Returns a SpaceClient for the space with the specified ID.
   * Use the returned client to manage subspaces, entities, and files within the space.
   * @param spaceId - The unique identifier of the space {@link Project.id} */
  id(spaceId: string): SpaceClient {
    return new SpaceClient(api.grok_Dapi_Spaces_Id(this.dart, spaceId));
  }
}

/**
 * Client for working with a specific space.
 * Provides methods for managing subspaces, entities (scripts, queries, connections), and files. */
export class SpaceClient {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /**
   * Adds a child space (subspace) under this space.
   * The subspace will have its own storage with a namespaced path derived from the parent hierarchy
   * (e.g., 'ParentSpace:ChildSpace:').
   * You can provide existing space, then based on the {@link link} parameter the subspace will be moved
   * or reference will be created.
   * @param childSpace - The subspace name or an existing space.
   * @param link - If true, creates a link reference instead of moving the subspace (default: false). When only the subspace name is provided,
   * this parameter is ignored. */
  addSubspace(childSpace: Project | string, link: boolean = false): Promise<Project> {
    return toJs(api.grok_SpaceClient_AddSubspace(this.dart, childSpace instanceof Project ? childSpace.dart : childSpace, link));
  }

  /** Checks if a subspace with the given name exists in the current space. */
  subspaceExists(name: string): Promise<boolean> {
    return api.grok_SpaceClient_SubspaceExists(this.dart, name);
  }

  /**
   * Adds an entity (script, query, connection, or another space) to this space.
   * When link is false (default), the entity is moved from its current location to this space.
   * When link is true, creates a reference without moving - the entity can appear in multiple spaces.
   * Moving a space preserves all its files, nested subspaces, and contained entities.
   * @param entityId - The unique identifier of the entity to add
   * @param link - If true, creates a link reference instead of moving the entity (default: false) */
  addEntity(entityId: string, link: boolean = false): Promise<void> {
    return api.grok_SpaceClient_AddEntity(this.dart, entityId, link);
  }

  /**
   * Removes an entity from this space.
   * Note: If the entity is linked, only reference will be deleted.
   * @param entityId - {@link Entity.id}
   */
  removeEntity(entityId: string): Promise<void> {
    return api.grok_SpaceClient_RemoveEntity(this.dart, entityId);
  }

  /** Returns a client for accessing children of this space */
  get children(): SpaceChildrenClient {
    return new SpaceChildrenClient(api.grok_SpaceClient_Children(this.dart));
  }

  /** Returns a client for accessing files in this space */
  get files(): SpaceFilesClient {
    return new SpaceFilesClient(api.grok_SpaceClient_Files(this.dart));
  }
}

/**
 * Client for querying and filtering children of a space.
 * Children can include subspaces (Projects), files (FileInfo), and various entities
 * like scripts, queries, and connections. */
export class SpaceChildrenClient extends HttpDataSource<Entity> {
  constructor(s: any) {
    super(s);
  }

  /**
   * Filters the children by entity types and whether to include linked items.
   * By default, only directly owned children are returned (not links).
   * @param types - Comma-separated list of entity types to include (e.g., 'Script,DataQuery')
   * @param includeLinked - If true, includes linked references in addition to owned children (default: false)
   * @returns A new SpaceChildrenClient with the filter applied
   */
  filter(types: string, includeLinked: boolean = false): SpaceChildrenClient {
    return new SpaceChildrenClient(api.grok_SpaceChildrenClient_Filter(this.dart, types, includeLinked));
  }
}

/**
 * Client for file operations within a space's storage.
 * Files in a space are stored in the space's dedicated storage connection.
 * Directory operations are synchronized with the space hierarchy - creating a directory
 * creates a corresponding subspace, and renaming a directory renames the subspace. */
export class SpaceFilesClient {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Checks if a file exists */
  exists(file: FileInfo | string): Promise<boolean> {
    return api.grok_SpaceFilesClient_Exists(this.dart, file instanceof FileInfo ? file.dart : file);
  }

  /** Reads file content as bytes */
  readAsBytes(file: FileInfo | string): Promise<Uint8Array> {
    return api.grok_SpaceFilesClient_ReadAsBytes(this.dart, file instanceof FileInfo ? file.dart : file);
  }

  /** Reads file content as string */
  readAsString(file: FileInfo | string): Promise<string> {
    return api.grok_SpaceFilesClient_ReadAsString(this.dart, file instanceof FileInfo ? file.dart : file);
  }

  /** Uploads file content */
  write(file: FileInfo | string, bytes: number[]): Promise<void> {
    return api.grok_SpaceFilesClient_Upload(this.dart, file instanceof FileInfo ? file.dart : file, bytes);
  }

  writeString(file: FileInfo | string, data: string): Promise<void> {
    return api.grok_SpaceFilesClient_UploadString(this.dart, file instanceof FileInfo ? file.dart : file, data);
  }

  /**
   * Creates a directory within this space's storage.
   * This automatically creates a corresponding subspace with the same name (normalized).
   * The subspace will have its own namespaced storage path. */
  createDirectory(file: FileInfo | string): Promise<void> {
    return api.grok_SpaceFilesClient_CreateDirectory(this.dart, file instanceof FileInfo ? file.dart : file);
  }

  /** Renames a file or directory. If renaming directory, underlying subspace will be renamed. */
  rename(file: FileInfo | string, newName: string): Promise<void> {
    return api.grok_SpaceFilesClient_Rename(this.dart, file instanceof FileInfo ? file.dart : file, newName);
  }

  /** Moves files to a new path */
  move(files: (FileInfo | string)[], newPath: FileInfo | string): Promise<void> {
    const filesArg = files.map(f => f instanceof FileInfo ? f.dart : f);
    const newPathArg = newPath instanceof FileInfo ? newPath.dart : newPath;
    return api.grok_SpaceFilesClient_Move(this.dart, filesArg, newPathArg);
  }

  /** Copies files to a destination path */
  copy(files: (FileInfo | string)[], destinationPath: string): Promise<void> {
    const filesArg = files.map(f => f instanceof FileInfo ? f.dart : f);
    return api.grok_SpaceFilesClient_Copy(this.dart, filesArg, destinationPath);
  }

  /** Deletes a file or directory */
  delete(file: FileInfo | string): Promise<void> {
    return api.grok_SpaceFilesClient_Delete(this.dart, file instanceof FileInfo ? file.dart : file);
  }
}

/**
 * Functionality for working with remote tables
 * @extends HttpDataSource
 * */
export class TablesDataSource extends HttpDataSource<TableInfo> {
  /** @constructs TablesDataSource*/
  constructor(s: any) {
    super(s);
  }

  /** Saves a dataframe remotely.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/save-and-load-df} */
  uploadDataFrame(dataFrame: DataFrame): Promise<string> {
    return api.grok_Dapi_TablesDataSource_UploadDataFrame(this.dart, dataFrame.dart);
  }

  /** Loads a dataframe by id.
   * Sample: {@link https://public.datagrok.ai/js/samples/data-access/save-and-load-df}
   * @param {string} id - dataframe id */
  getTable(id: string): Promise<DataFrame> {
    return api.grok_Dapi_TablesDataSource_GetTable(this.dart, id);
  }
}

export class DockerDataSource {
  /**DockerImages API endpoint
   * @type {HttpDataSource<DockerImage>} */
  get dockerImages(): DockerImagesDataSource {
    return new DockerImagesDataSource(api.grok_Dapi_DockerImages());
  }

  /**Dockerfiles API endpoint
   * @type {HttpDataSource<DockerImage>} */
  get dockerContainers(): DockerContainersDataSource {
    return new DockerContainersDataSource(api.grok_Dapi_DockerContainers());
  }

  getServiceLogs(serviceName: string, limit: number): Promise<string> {
    return api.grok_Dapi_DockersDataSource_GetServiceLogs(serviceName, limit);
  }

  getAvailableServices(): Promise<string[]> {
    return api.grok_Dapi_DockersDataSource_GetAvailableServices();
  }
}

/** Functionality to work with Docker images. See also {@link DockerContainersDataSource}.
 * @extends HttpDataSource */
export class DockerImagesDataSource extends HttpDataSource<DockerImage> {

  constructor(s: any) {
    super(s);
  }

  /**
   * Rebuilds Docker image.
   * @param imageId - ID of the {@link DockerImage} to rebuild.
   * @returns {Promise<void>} - promise that resolves with void or throws Exception if something went wrong.
   */
  rebuild(imageId: string): Promise<void> {
    return api.grok_Dapi_DockerImagesDataSource_Rebuild(this.dart, imageId);
  }
}

/** Functionality to work with Docker containers.
 * See help: {@link https://datagrok.ai/help/develop/how-to/docker_containers}.
 * @extends HttpDataSource */
export class DockerContainersDataSource extends HttpDataSource<DockerContainer> {
  constructor(s: any) {
    super(s);
  }

  /**
   * Runs container.
   * @param containerId - ID of the {@link DockerContainer} to be run.
   * @param awaitStart - if [true] promise will not be resolved until the container is started,
   * otherwise, it doesn't wait for start and resolves immediately after the container is queued for start.
   * @returns {Promise<void>} - promise that resolves with void or throws Exception if something went wrong.
   */
  run(containerId: string, awaitStart: boolean = false): Promise<void> {
    return api.grok_Dapi_DockerContainersDataSource_Run(this.dart, containerId, awaitStart);
  }

  /**
   * Stops container.
   * @param containerId - ID of the {@link DockerContainer} to be stopped.
   * @param awaitStop - if [true] promise will not be resolved until the container is stopped,
   * otherwise, it doesn't wait for a stop and resolves immediately after the container is queued for a stop.
   * @returns {Promise<void>} or throws Exception if something went wrong.
   */
  stop(containerId: string, awaitStop: boolean = false): Promise<void> {
    return api.grok_Dapi_DockerContainersDataSource_Stop(this.dart, containerId, awaitStop);
  }

  /**
   * Proxies URL requests to Docker containers via Datagrok server with the same interface as [fetch](https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API). Returns response
   * from the containers as it is. If an error occurs on the server side returns a response with "application/json"
   * Content-Type and JSON body with field "datagrok-error" that describes the cause. If container status is incorrect
   * for performing requests returns a response with a 400 status code. If something goes wrong in the server workflow,
   * it returns a response with a 500 status code. Any other cases are the result of direct requests to the container itself.
   * @param containerId - ID of the {@link DockerContainer} to which the http request should be sent.
   * @param path - URI without scheme and authority component.
   * @param params - parameters of the request.
   * @returns {Promise<Response>} - promise that resolves with [Response](https://developer.mozilla.org/en-US/docs/Web/API/Response).
   */
  async fetchProxy(containerId: string, path: string, params?: RequestInit): Promise<Response> {
    params ??= {};
    params.method ??= 'GET';
    params.credentials = 'include';
    if (!path.startsWith('/')) path = `/${path}`;
    return fetch(`${api.grok_Dapi_Root()}/docker/containers/proxy/${containerId}${path}`, params);
  }

  /**
   * Proxies WebSocket connection to Docker containers via Datagrok server. Returns ready WebSocket that is connected through the server to the
   * Docker container WebSocket endpoint. If container status is incorrect or there is error while establishing WebSocket connection to Docker container, caller will receive an error.
   * After the WebSocket is returned, caller can do anything with it and should take care of reconnection. After the caller closes the connection, server will close
   * proxied Docker WebSocket connection.
   * @param containerId - ID of the {@link DockerContainer} to which the WebSocket connection will be established.
   * @param path - URI without scheme and authority component that points to endpoint inside  the Docker container
   * @param timeout - Timeout in ms for initial connection establishment. Set it to higher values if you are using container with on_demand configuration set to `true`.
   */
  async webSocketProxy(containerId: string, path: string, timeout: number = 60000): Promise<WebSocket> {
    if (!path.startsWith('/')) path = `/${path}`;

    return new Promise<WebSocket>((resolve, reject) => {
      const socket = new WebSocket(`${api.grok_Dapi_WS_Root()}/docker/containers/proxy-ws/${containerId}${path}`);

      let timeoutTimer = setTimeout(() => {
        cleanup();
        socket.close(4001, "Timeout waiting for Docker container");
        reject(new Error("Timeout waiting for Docker container"));
      }, timeout);

      // Wait for the CONNECTED message from the server.
      // Message indicates that this and Docker container WebSockets are connected to each other.
      const onMessage = (event: MessageEvent) => {
        if (event.data === "CONNECTED") {
          clearTimeout(timeoutTimer);
          cleanup();
          resolve(socket);
        }
      };

      const onClose = (event: CloseEvent) => {
        clearTimeout(timeoutTimer);
        cleanup();
        reject(new Error(`Could not open WebSocket connection: ${event.reason}`));
      };

      // Unfortunately, event doesn't have error reference
      const onError = (_: Event) => {
        clearTimeout(timeoutTimer);
        cleanup();
        reject(new Error("WebSocket encountered an error"));
      };

      function cleanup() {
        socket.removeEventListener("message", onMessage);
        socket.removeEventListener("close", onClose);
        socket.removeEventListener("error", onError);
      }

      socket.addEventListener("message", onMessage);
      socket.addEventListener("close", onClose);
      socket.addEventListener("error", onError);
    });
  }

  /**
   * This is the synchronous version of the function {@link webSocketProxy}. Note that the container won't be
   * ready to accept messages immediately after this function returns. If your application logic requires sending
   * a message right away, consider using the asynchronous version.
   * @param containerId - ID of the {@link DockerContainer} to which the WebSocket connection will be established.
   * @param path - URI without scheme and authority component that points to endpoint inside  the Docker container
   * @param timeout - Timeout in ms for initial connection establishment. Set it to higher values if you are using container with on_demand configuration set to `true`.
   */
  webSocketProxySync(containerId: string, path: string, timeout: number = 60000): WebSocket {
    const socket = new WebSocket(`${api.grok_Dapi_WS_Root()}/docker/containers/proxy-ws/${containerId}${path}`);
    new Promise((resolve, reject) => {
      let timeoutTimer = setTimeout(() => {
        cleanup();
        socket.close(4001, "Timeout waiting for Docker container");
        reject(new Error("Timeout waiting for Docker container"));
      }, timeout);

      const onMessage = (event: MessageEvent) => {
        if (event.data === "CONNECTED") {
          clearTimeout(timeoutTimer);
          cleanup();
          resolve(socket);
        }
      };

      const onClose = (event: CloseEvent) => {
        clearTimeout(timeoutTimer);
        cleanup();
        reject(new Error(`Could not open WebSocket connection: ${event.reason}`));
      };

      const onError = (_: Event) => {
        clearTimeout(timeoutTimer);
        cleanup();
        socket.close(4001, "WebSocket encountered an error");
        reject(new Error("WebSocket encountered an error"));
      };

      function cleanup() {
        socket.removeEventListener("message", onMessage);
        socket.removeEventListener("close", onClose);
        socket.removeEventListener("error", onError);
      }

      socket.addEventListener("message", onMessage);
      socket.addEventListener("close", onClose);
      socket.addEventListener("error", onError);
    }).catch((e) => console.error(e));
    return socket;
  }

  /**
   * Returns container's logs or throws Exception with the cause.
   * @param containerId - ID of the {@link DockerContainer} whose logs is to be obtained.
   * @param limit - maximum line count of logs.
   * @returns string - container logs or null if there are no logs.
   */
  getContainerLogs(containerId: string, limit: number = 10000): Promise<string | null> {
    return api.grok_Dapi_DockerContainersDataSource_GetContainerLogs(this.dart, containerId, limit);
  }
}

export class UserReportsDataSource extends HttpDataSource<UserReport> {
  constructor(s: any) {
    super(s);
  }


  getReports(num?: number): Promise<DataFrame> {
    return api.grok_Reports_Get(num);
  }
}

export class UserReportsRulesDataSource extends HttpDataSource<UserReportsRule> {
  constructor(s: any) {
    super(s);
  }
}

export class NotificationsDataSource extends HttpDataSource<UserNotification> {
  constructor(s: any) {
    super(s);
  }

  forCurrentUser(): NotificationsDataSource {
    return new NotificationsDataSource(api.grok_Dapi_Notifications_ForCurrentUser());
  }

  async countUnread(): Promise<number> {
    return api.grok_Dapi_Notifications_CountUnread();
  }
}

export class LogDataSource extends HttpDataSource<LogEvent> {
  constructor(s: any) {
    super(s);
  }

  /** Activity API endpoint
   *  @type {ActivityDataSource} */
  get activity(): ActivityDataSource {
    return new ActivityDataSource(api.grok_Dapi_Activity());
  }

  where(options?: {entityId?: string, start?: dayjs.Dayjs, end?: dayjs.Dayjs, favoritesOnly?: boolean}): LogDataSource {
    return new LogDataSource(api.grok_Dapi_Log_Where(this.dart, options?.entityId ?? '', toDart(options?.start), toDart(options?.end), options?.favoritesOnly ?? false));
  }
}

export class ActivityDataSource extends HttpDataSource<LogEvent> {
  constructor(s: any) {
    super(s);
  }

  where(options?: {userId?: string, start?: dayjs.Dayjs, end?: dayjs.Dayjs}): ActivityDataSource {
    return new ActivityDataSource(api.grok_Dapi_Activity_Where(this.dart, options?.userId ?? '', toDart(options?.start), toDart(options?.end)));
  }
}

/**
 * Provides access to file operations in the Datagrok file system.
 *
 * Allows reading, writing, listing, and managing files and directories
 * in user file shares and data connections.
 *
 * Access via `grok.dapi.files`.
 *
 * @example
 * // List files
 * const files = await grok.dapi.files.list('System:AppData/MyApp');
 *
 * // Read a file
 * const content = await grok.dapi.files.readAsText('System:AppData/MyApp/config.json');
 *
 * // Write a file
 * await grok.dapi.files.writeAsText('System:AppData/MyApp/output.txt', 'Hello, World!');
 */
export class FilesDataSource {
  private readonly root: string;
  constructor(root: string = '') {
    this.root = root;
  };

  /** Checks if a file exists.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  exists(file: FileInfo | string): Promise<boolean> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_Exists(file);
  }

  private setRoot(file: FileInfo | string): string {
    if (typeof (file) == 'string') {
      file = `${this.root}${this.root != '' ? '/' : ''}${file}`;
      return <string>file;
    } else {
      return file.fullPath;
    }
  }

  /** Deletes a file.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  delete(file: FileInfo | string): Promise<void> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_Delete(file);
  }

  /** Moves a file.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  move(files: FileInfo[] | string[], newPath: string): Promise<void> {
    for (let i = 0; i < files.length; i++)
      files[i] = this.setRoot(files[i]);
    newPath = this.setRoot(newPath);

    return api.grok_Dapi_UserFiles_Move(files, newPath);
  }

  /** Renames a file.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  rename(file: FileInfo | string, newName: string): Promise<void> {
    file = this.setRoot(file);
    newName = this.setRoot(newName);

    return api.grok_Dapi_UserFiles_Rename(file, newName);
  }

  /** Lists files according to a search pattern.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files}
   * @param {FileInfo | string} file - folder
   * @param {boolean} recursive - whether to search in folders recursively
   * @param {string} searchPattern - search pattern, such as part of a filename or extension, e.g., "filename-prefix" and "csv" */
  async list(file: FileInfo | string, recursive: boolean = false, searchPattern: string | null = null): Promise<FileInfo[]> {
    file = this.setRoot(file);
    return toJs(await api.grok_Dapi_UserFiles_List(file, recursive, searchPattern, this.root));
  }

  /**
   * Reads the entire contents of a folder and returns an object.
   * The resulting object's keys are the file names relative to the folder path, and the corresponding values are of the Blob type.
   * @param recursive - whether to read files in folders recursively
   * @param ext - files extension
   */
  async readFilesAsBlobs(folder: FileInfo | string, recursive: boolean = false, ext: string | undefined = undefined): Promise<{[key: string]: Blob}> {
    const folderPath = this.setRoot(folder);
    const conn = folderPath.replace(":", ".").split('/')[0];
    let url = `${api.grok_Dapi_Root()}/connectors/connections/${conn}/folder/${folderPath.substring(folderPath.indexOf('/') + 1)}?recursive=${recursive}`;
    if (ext) {
      if (!ext.startsWith('.'))
        ext = `.${ext}`;
      url += `&ext=${ext}`;
    }
    const response = await fetch(url);
    if (response.status == 204)
      return {};
    const formData: FormData = await response.formData();
    const files: {[key: string]: any} = {};

    formData.forEach((value: Blob | string, filename) => {
      if (value instanceof Blob)
        files[filename] = value;
    });

    return files;
  }

  /**
   * Reads the entire contents of a folder and returns an object.
   * The resulting object's keys are the file names relative to the folder path, and the corresponding values are JSON objects.
   * If conversion to a JSON fails, the file will be skipped.
   * @param recursive - whether to read files in folders recursively
   * @param ext - files extension
   */
  async readFilesAsJson(folder: FileInfo | string, recursive: boolean = false, ext: string | undefined = undefined): Promise<{[key: string]: any}> {
    const filesBlobs: {[key: string]: Blob} = await this.readFilesAsBlobs(folder, recursive, ext);
    const jsons: {[key: string]: any} = {};
    for (const [name, blob] of Object.entries(filesBlobs)) {
      try {
        jsons[name] = JSON.parse(await blob.text());
      } catch (_) {}
    }
    return jsons;
  }

  /**
   * Reads the entire contents of a folder and returns an object.
   * The resulting object's keys are the file names relative to the folder path, and the corresponding values are strings.
   * If conversion to a string fails, the file will be skipped.
   * @param recursive - whether to read files in folders recursively
   * @param ext - files extension
   */
  async readFilesAsString(folder: FileInfo | string, recursive: boolean = false, ext: string | undefined = undefined): Promise<{[key: string]: string}> {
    const filesBlobs: {[key: string]: Blob} = await this.readFilesAsBlobs(folder, recursive, ext);
    const files: {[key: string]: string} = {};
    for (const [name, blob] of Object.entries(filesBlobs)) {
      try {
        files[name] = await blob.text();
      } catch (_) {}
    }
    return files;
  }


  /** Reads a file as string.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  readAsText(file: FileInfo | string): Promise<string> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_ReadAsText(file);
  }

  /** Reads CSV as DataFrame.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  async readCsv(file: FileInfo | string, options?: CsvImportOptions): Promise<DataFrame> {
    return DataFrame.fromCsv(await this.readAsText(file), options);
  }

  /** Reads a file as bytes.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  readAsBytes(file: FileInfo | string): Promise<Uint8Array> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_ReadAsBytes(file);
  }

  /** Reads a d42 file as a list of dataframes. */
  async readBinaryDataFrames(file: FileInfo | string): Promise<DataFrame[]> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_ReadBinaryDataFrames(file);
  }

  /** Writes a list of dataframes as a d42 file. */
  async writeBinaryDataFrames(file: FileInfo | string, dataFrames: DataFrame[]): Promise<void> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_WriteBinaryDataFrames(file, dataFrames.map((df) => df.dart));
  }

  /** Creates directory
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  async createDirectory(file: FileInfo | string): Promise<void> {
    file = this.setRoot(file);
    return api.grok_Dapi_UserFiles_CreateDirectory(file);
  }

  /** Writes a file.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  write(file: FileInfo | string, blob?: number[]): Promise<void> {
    if (!blob && ((file instanceof FileInfo && !file.data) || typeof file === 'string'))
      throw new Error('blob parameter should be presented');
    return api.grok_Dapi_UserFiles_Write(toDart(file), blob ?? (file as FileInfo).data);
  }

  /** Writes a text file.
   * Sample: {@link https://public.datagrok.ai/js/samples/dapi/files} */
  writeAsText(file: FileInfo | string, data: string): Promise<void> {
    file = this.setRoot(file);

    return api.grok_Dapi_UserFiles_WriteAsText(file, data);
  }
}

/**
 * @deprecated Use {@link FilesDataSource} instead. This alias is provided for backward compatibility.
 */
export const FileSource = FilesDataSource;
