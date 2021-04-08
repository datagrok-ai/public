import { DataFrame } from "./dataframe";
import {
  Credentials,
  DataConnection,
  DataJob,
  DataQuery,
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
  FileInfo
} from "./entities";
import {ViewLayout} from "./view";
import {toDart, toJs} from "./wrappers";

let api = <any>window;

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
 * */
export class Dapi {
  constructor() {
  }

  /** Retrieves entities from server by list of IDs
   *  @returns {Promise<List<Entity>>} */
  getEntities(ids: string[]): Promise<Entity[]> {
    return new Promise((resolve, reject) => api.grok_Dapi_Entities_GetEntities(ids, (q: any) => {
      return resolve(q.map(toJs));
    }, (e: any) => reject(e)));
  }

  /** Entities API endpoint
   *  @type {EntitiesDataSource} */
  get entities(): EntitiesDataSource {
    return new EntitiesDataSource(api.grok_Dapi_Entities(), (a: any) => new Entity(a));
  }

  /** Data Queries API endpoint
   *  @type {HttpDataSource<DataQuery>} */
  get queries(): HttpDataSource<DataQuery> {
    return new HttpDataSource(api.grok_Dapi_Queries(), (a: any) => new DataQuery(a));
  }

  /** Data Connections API endpoint
   *  @type {HttpDataSource<DataConnection>} */
  get connections(): HttpDataSource<DataConnection> {
    return new HttpDataSource(api.grok_Dapi_Connections(), (a: any) => new DataConnection(a));
  }

  /** Credentials API endpoint
   *  @type {CredentialsDataSource} */
  get credentials(): CredentialsDataSource {
    return new CredentialsDataSource(api.grok_Dapi_Credentials(), (a: any) => new Credentials(a));
  }

  /** Data Jobs API endpoint
   *  @type {HttpDataSource<DataJob>} */
  get jobs(): HttpDataSource<DataJob> {
    return new HttpDataSource(api.grok_Dapi_Jobs(), (a: any) => new DataJob(a));
  }

  /** Jupyter Notebooks API endpoint
   *  @type {HttpDataSource<Notebook>} */
  get notebooks(): HttpDataSource<Notebook> {
    return new HttpDataSource(api.grok_Dapi_Notebooks(), (a: any) => toJs(a));
  }

  /** Predictive Models API endpoint
   *  @type {HttpDataSource<Model>} */
  get models(): HttpDataSource<Model> {
    return new HttpDataSource(api.grok_Dapi_Models(), (a: any) => new Model(a));
  }

  /** Packages API endpoint
   *  @type {HttpDataSource<Package>} */
  get packages(): HttpDataSource<Package> {
    return new HttpDataSource(api.grok_Dapi_Packages(), (a: any) => new Package(a));
  }

  /** View Layouts API endpoint
   *  @type {LayoutsDataSource} */
  get layouts(): LayoutsDataSource {
    return new LayoutsDataSource(api.grok_Dapi_Layouts(), (a: any) => new ViewLayout(a));
  }

  /** Data Table Infos API endpoint
   *  @type {TablesDataSource} */
  get tables(): TablesDataSource {
    return new TablesDataSource(api.grok_Dapi_Tables(), (a: any) => new TableInfo(a));
  }

  /** Users API endpoint
   *  @type {UsersDataSource} */
  get users(): UsersDataSource {
    return new UsersDataSource(api.grok_Dapi_Users(), (a: any) => new User(a));
  }

  /** Groups API endpoint
   *  @type {GroupsDataSource} */
  get groups(): GroupsDataSource {
    return new GroupsDataSource(api.grok_Dapi_Groups(), (a: any) => toJs(a));
  }

  /** Permissions API endpoint
   * @type {PermissionsDataSource} */
  get permissions(): PermissionsDataSource {
    return new PermissionsDataSource();
  }

  /** Scripts API endpoint
   *  @type {HttpDataSource<Script>} */
  get scripts(): HttpDataSource<Script> {
    return new HttpDataSource(api.grok_Dapi_Scripts(), (a: any) => new Script(a));
  }

  /** Projects API endpoint
   *  @type {HttpDataSource<Project>} */
  get projects(): HttpDataSource<Project> {
    return new HttpDataSource(api.grok_Dapi_Projects(), (a: any) => new Project(a));
  }

  /** Environments API endpoint
   *  @type {HttpDataSource<ScriptEnvironment>} */
  get environments(): HttpDataSource<ScriptEnvironment> {
    return new HttpDataSource(api.grok_Dapi_Environments(), (a: any) => new ScriptEnvironment(a));
  }

  /** Users Data Storage API endpoint
   *  @type {UserDataStorage} */
  get userDataStorage(): UserDataStorage {
    return new UserDataStorage();
  }


  /** Users Files management API endpoint
   *  @type {FileSource} */
  get files(): FileSource {
    return new FileSource();
  }

  /** Proxies URL request via Datagrok server with same interface as "fetch".
   * @deprecated
   * @param {string} method
   * @param {string} url
   * @param {Object} headers
   * @param {Object} body
   * @returns {Promise<Object>} */
  async proxyFetch(method: string, url: string, headers: Record<string, string>, body: object = {}): Promise<object> {
    headers['Accept'] = 'application/json';
    headers['original-url'] = `${url}`;
    headers['original-method'] = method;
    let params = {
      method: 'POST',
      headers: headers,
      body: JSON.stringify(body)
    };
    return await fetch('/api/connectors/proxy', params);
  }

  /** Proxies URL request via Datagrok server with same interface as "fetch".
   * @param {String} url
   * @param {Object} params
   * @returns {Promise<Object>} */
  async fetchProxy(url: string, params: RequestInit): Promise<object> {
    if (params == null)
      params = {};
    if (params.headers == null)
      params.headers = {};
    if (params.method == null)
      params.method = 'GET';
    // @ts-ignore
    params.headers['original-url'] = `${url}`;
    // @ts-ignore
    params.headers['original-method'] = params.method;
    params.method = 'POST';
    return fetch('/api/connectors/proxy', params);
  }

  /** Administering API endpoint
   *  @type {AdminDataSource} */
  get admin(): AdminDataSource {
    return new AdminDataSource(api.grok_Dapi_Admin());
  }

  /** Logging API endpoint
   *  @type {HttpDataSource<LogEvent>} */
  get log(): HttpDataSource<LogEvent> {
    return new HttpDataSource(api.grok_Dapi_Log(), (a: any) => new LogEvent(a));
  }

  /** Logging API endpoint
   *  @type {HttpDataSource<LogEventType>} */
  get logTypes(): HttpDataSource<LogEventType> {
    return new HttpDataSource(api.grok_Dapi_LogTypes(), (a: any) => new LogEventType(a));
  }
}


/**
 * Common functionality for handling collections of entities stored on the server.
 * Works with Datagrok REST API, allows to get filtered and paginated lists of entities,
 * Can be extended with specific methods. (i.e. {@link UsersDataSource})
 */
export class HttpDataSource<T> {
  s: any;
  entityToJs: any;

  /** @constructs HttpDataSource */
  constructor(s: any, instance: any) {
    this.s = s;
    this.entityToJs = instance;
  }

  /** Returns all entities that satisfy the filtering criteria (see {@link filter}).
   *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  Smart filter: {@link https://datagrok.ai/help/overview/smart-search}
   *  @param {Object} options
   *  @param {int} options.pageSize
   *  @param {int} options.pageNumber
   *  @param {string} options.order
   *  @param {string} options.filter
   *  @returns {Promise<object[]>}  */
  list(options: {pageSize?: number, pageNumber?: number, filter?: string, order?: string} = {}): Promise<T[]> {
    if (options.pageSize !== undefined)
      this.by(options.pageSize);
    if (options.pageNumber !== undefined)
      this.page(options.pageNumber);
    if (options.filter !== undefined)
      this.filter(options.filter);
    if (options.order !== undefined)
      this.order(options.order);
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_DataSource_List(this.s, (q: any) => resolve(q.map(s)), (e: any) => reject(e)));
  }

  /** Returns fist entity that satisfy the filtering criteria (see {@link filter}).
   *  @returns Promise<object>  */
  first(): Promise<T> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_DataSource_First(this.s, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

  /** Returns an entity with the specified id.
   *  Throws an exception if an entity does not exist, or is not accessible in the current context.
   *  @param {string} id - GUID of the corresponding object
   *  @returns {Promise<object>} - entity. */
  find(id: string): Promise<T> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_DataSource_Find(this.s, id, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

  /** Saves an entity
   *  @param {Entity} e
   *  @returns {Promise} - entity. */
  save(e: Entity): Promise<T> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_DataSource_Save(this.s, e.d, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

  /** Deletes an entity
   *  @param {Entity} e
   *  @returns {Promise} */
  delete(e: Entity): Promise<void> {
    return new Promise((resolve, reject) => api.grok_DataSource_Delete(this.s, e.d, () => resolve(), (e: any) => reject(e)));
  }

  by(i: number): HttpDataSource<T> {
    this.s = api.grok_DataSource_By(this.s, i);
    return this;
  }

  page(i: number): HttpDataSource<T> {
    this.s = api.grok_DataSource_Page(this.s, i);
    return this;
  }

  /** Returns next page of all entities that satisfy the filtering criteria (see {@link filter}).
   *  Works only if pageSize was set during previous list() call
   *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  @returns {HttpDataSource} */
  nextPage(): HttpDataSource<T> {
    this.s = api.grok_DataSource_NextPage(this.s);
    return this;
  }

  /** Applies filter to current request.
   *  Also can be set with {@link list} method "options" parameter
   *  See example: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
   *  Smart filter: {@link https://datagrok.ai/help/overview/smart-search}
   *  @param {string} w
   *  @returns {HttpDataSource} */
  filter(w: string): HttpDataSource<T> {
    this.s = api.grok_DataSource_WhereSmart(this.s, w);
    return this;
  }

  /** Instructs data source to return results in the specified order.
   * @param {string} fieldName
   * @param {boolean} desc
   * @returns {HttpDataSource} */
  order(fieldName: string, desc: boolean = false): HttpDataSource<T> {
    this.s = api.grok_DataSource_Order(this.s, fieldName, desc);
    return this;
  }

  /** Includes entity in the result
   * @param {string} include
   * @returns {HttpDataSource} */
  include(include: string): HttpDataSource<T> {
    this.s = api.grok_DataSource_Include(this.s, include);
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
  /** @constructs UsersDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** Returns current user
   * @returns {Promise<User>} */
  current(): Promise<User> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_UsersDataSource_Current(this.s, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

  /** Returns current session
   * @returns {Promise<UserSession>} */
  currentSession(): Promise<UserSession> {
    return new Promise((resolve, reject) => api.grok_UsersDataSource_CurrentSession(this.s, (q: any) => resolve(toJs(q)), (e: any) => reject(e)));
  }
}

export class AdminDataSource {
  s: any;
  /** @constructs AdminDataSource*/
  constructor(s: any) {
    this.s = s;
  }

  /** Returns information about the services
   *  @returns {Promise<Map>} */
  getServiceInfos(): Promise<object> {
    return new Promise((resolve, reject) => api.grok_Dapi_Admin_GetServiceInfos(this.s, (q: any) => resolve(toJs(q)), (e: any) => reject(e)));
  }
}

/**
 * Functionality for handling groups collection from server
 * Allows to manage {@link Group}
 * @extends HttpDataSource
 * */
export class GroupsDataSource extends HttpDataSource<Group> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** Creates a new group
   *  @param {String}  name
   *  @returns {Promise<Group>} - Group. */
  createNew(name: string): Promise<Group> {
    let g = new Group(api.grok_Group(name));
    return this.save(g);
  }

  /** Returns group user
   *  @param {Group} group
   *  @returns {Promise<Group>} - Group. */
  getUser(group: Group): Promise<Group> {
    return new Promise((resolve, reject) => api.grok_Dapi_Get_GroupUser(group.d, (q: any) => resolve(toJs(q)), (e: any) => reject(e)));
  }

  /** Adds a member to the group
   * @param {Group} g
   * @param {Group} m
   * @returns {Promise} */
  async addMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.addMember(m);
    await this.saveRelations(g);
  }

  /** Adds an admin member to the group
   * @param {Group} g
   * @param {Group} m
   * @returns {Promise} */
  async addAdminMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.addAdminMember(m);
    await this.saveRelations(g);
  }

  /** Removes a member from the group
   * @param {Group} g
   * @param {Group} m
   * @returns {Promise} */
  async removeMember(g: Group, m: Group): Promise<void> {
    g = await this.find(g.id);
    g.removeMember(m);
    await this.saveRelations(g);
  }

  /** Adds the group to another one
   * @param {Group} g
   * @param {Group} parent
   * @returns {Promise} */
  async includeTo(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.includeTo(parent);
    await this.saveRelations(g);
  }


  /** Adds the group to another one as admin
   * @param {Group} g
   * @param {Group} parent
   * @returns {Promise} */
  async includeAdminTo(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.includeAdminTo(parent);
    await this.saveRelations(g);
  }

  /** Removes a membership from the group
   * @param {Group} g
   * @param {Group} parent
   * @returns {Promise} */
  async excludeFrom(g: Group, parent: Group): Promise<void> {
    g = await this.find(g.id);
    g.excludeFrom(parent);
    await this.saveRelations(g);
  }

  /** Saves a group with relations
   *  @param {Group} e
   *  @returns {Promise<Group>} - Group. */
  saveRelations(e: Group): Promise<Group> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_GroupsDataSource_Save(this.s, e.d, (q: any) => resolve(s(q)), (e: any) => reject(e)));
  }

}


/**
 * Functionality for handling entities collection from server
 * Allows to manage {@link Entity}
 * @extends HttpDataSource
 * */
export class EntitiesDataSource extends HttpDataSource<Entity> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** Allows to set properties for entities
   * @param {List<Map>} props
   * @returns {Promise} */
  saveProperties(props: Map<Property, any>): Promise<void> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_SaveProperties(this.s, props, (_: any) => resolve(), (e: any) => reject(e)));
  }

  /** Returns entity properties
   * @param {Entity} entity
   * @returns {Promise<Map>} props */
  getProperties(entity: Entity): Promise<Map<Property, any>> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_GetProperties(this.s, entity.d, (p: Map<Property, any> | PromiseLike<Map<Property, any>>) => resolve(p), (e: any) => reject(e)));
  }

  /** Deletes entity properties
   * @param {List<Map>} props
   * @returns {Promise} */
  deleteProperties(props: Map<Property, any>): Promise<void> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_DeleteProperties(this.s, props, (_: any) => resolve(), (e: any) => reject(e)));
  }
}

/**
 * Functionality for handling credentials collection from server and working with credentials remote endpoint
 * Allows to manage {@link Credentials}
 * See also: {@link https://datagrok.ai/help/govern/security}
 * @extends HttpDataSource
 * */
export class CredentialsDataSource extends HttpDataSource<Credentials> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** Returns credentials for entity
   * @param {Entity} e
   * @returns {Promise<Credentials>} */
  forEntity(e: Entity): Promise<Credentials> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_CredentialsDataSource_ForEntity(this.s, e.d, (c: any) => resolve(s(c)), (e: any) => reject(e)));
  }
}

/**
 * Functionality for handling layouts collection from server
 * Allows to manage {@link ViewLayout}
 * @extends HttpDataSource
 * */
export class LayoutsDataSource extends HttpDataSource<ViewLayout> {
  /** @constructs CredentialsDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** Returns layouts that applicable to the table
   * @param {DataFrame} t
   * @returns {Promise<List<ViewLayout>>} */
  getApplicable(t: DataFrame): Promise<ViewLayout[]> {
    let s = this.entityToJs;
    return new Promise((resolve, reject) => api.grok_LayoutsDataSource_Applicable(this.s, t.d, (q: any[]) => resolve(q.map(s)), (e: any) => reject(e)));
  }
}

export class PermissionsDataSource {
  constructor() {
  };

  /** Gets all the permissions granted on entity
   * @param {Entity} e
   * @returns {Promise<Map>} permissions
   * */
  // { [key:string]:number; }
  get(e: Entity): Promise<Map<string, Group[]>> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_Get_Permissions(e.d, (data: any) => {
        data.view = toJs(data.view);
        data.edit = toJs(data.edit);
        resolve(data);
      }, (e: any) => reject(e)));
  }

  /** Grants permission on entity to the group
   * @param {Entity} e
   * @param {Group} g
   * @param {boolean} edit allow to edit entity
   * @returns {Promise}
   * */
  grant(e: Entity, g: Group, edit: boolean): Promise<any> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_Set_Permission(e.d, g.d, edit, (data: any) => resolve(data), (e: any) => reject(e)));
  }

  /** Revokes permission on entity from the group
   * @param {Entity} e
   * @param {Group} g
   * @returns {Promise}
   * */
  revoke(g: Group, e: Entity): Promise<any> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_Delete_Permission(e.d, g.d, (data: any) => resolve(data), (e: any) => reject(e)));
  }
}

/**
 * Functionality for working with remote Users Data Storage
 * Remote storage allows to save key-value pairs on the Datagrok server for further use
 * */
export class UserDataStorage {
  constructor() {
  }

  /** Saves a single value to Users Data Storage
   * @param {string} name Storage name
   * @param {string} key
   * @param {string} value
   * @param {boolean} currentUser Value should be available only for current user
   * @returns {Promise}*/
  postValue(name: string, key: string, value: string, currentUser: boolean = true): Promise<void> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser, () => resolve(), (e: any) => reject(e)));
  }

  /** Saves a map to Users Data Storage, will be appended to existing data
   * @param {string} name Storage name
   * @param {Map} data
   * @param {boolean} currentUser Value should be available only for current user
   * @returns {Promise}*/
  post(name: string, data: Map<string, string>, currentUser: boolean = true): Promise<void> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_Post(name, data, currentUser, () => resolve(), (e: any) => reject(e)));
  }

  /** Saves a map to Users Data Storage, will replace existing data
   * @param {string} name Storage name
   * @param {Map} data
   * @param {boolean} currentUser Value should be available only for current user
   * @returns {Promise}*/
  put(name: string, data: Map<string, string>, currentUser: boolean = true): Promise<void> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_Put(name, data, currentUser, () => resolve(), (e: any) => reject(e)));
  }

  /** Retrieves a map from Users Data Storage
   * @param {string} name Storage name
   * @param {boolean} currentUser get a value from a current user storage
   * @returns {Promise<Map>} */
  get(name: string, currentUser: boolean = true): Promise<Map<string, string>> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_Get(name, currentUser, (data: any) => resolve(data), (e: any) => reject(e)));
  }

  /** Retrieves a single value from Users Data Storage
   * @param {string} name Storage name
   * @param {string} key Value key
   * @param {boolean} currentUser get a value from a current user storage
   * @returns {Promise<string>} */
  getValue(name: string, key: string, currentUser: boolean = true): Promise<string> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_GetValue(name, key, currentUser, (value: string | PromiseLike<string>) => resolve(value), (e: any) => reject(e)));
  }

  /** Removes a single value from Users Data Storage
   * @param {string} name Storage name
   * @param {string} key Value key
   * @param {boolean} currentUser get a value from a current user storage
   * @returns {Promise} */
  remove(name: string, key: string, currentUser: boolean = true): Promise<void> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_UserDataStorage_Delete(name, key, currentUser, () => resolve(), (e: any) => reject(e)));
  }
}

/**
 * Functionality for working with remote tables
 * @extends HttpDataSource
 * */
export class TablesDataSource extends HttpDataSource<TableInfo> {
  /** @constructs TablesDataSource*/
  constructor(s: any, instance: any) {
    super(s, instance);
  }

  /** @returns {Promise<string>} */
  uploadDataFrame(dataFrame: DataFrame): Promise<string> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_TablesDataSource_UploadDataFrame(dataFrame.d, (id: string | PromiseLike<string>) => resolve(id), (e: any) => reject(e)));
  }

  /**
   * @param {string} id
   * @returns {Promise<DataFrame>} */
  getTable(id: string): Promise<DataFrame> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_TablesDataSource_GetTable(id, (df: any) => resolve(toJs(df)), (e: any) => reject(e)));
  }

}

export class FileSource {
  constructor() {
  };

  /** Check if file exists
   * @param {FileInfo | string} file
   * @returns {Promise<Boolean>} */
  exists(file: FileInfo | string): Promise<boolean> {
    return new Promise((resolve, reject) =>
        api.grok_Dapi_UserFiles_Exists(file, (data: boolean | PromiseLike<boolean>) => resolve(data), (e: any) => reject(e)));
  }

  /** Delete file
   * @param {FileInfo | string} file
   * @returns {Promise} */
  delete(file: FileInfo | string): Promise<void> {
    return new Promise((resolve, reject) => api.grok_Dapi_UserFiles_Delete(file, () => resolve(), (e: any) => reject(e)));
  }

  /** Move file
   * @param {List<FileInfo | string>} files
   * @param {string} newPath
   * @returns {Promise} */
  move(files: FileInfo[] | string[], newPath: string): Promise<void> {
    return new Promise((resolve, reject) => api.grok_Dapi_UserFiles_Move(files, newPath, () => resolve(), (e: any) => reject(e)));
  }

  /** Rename file
   * @param { FileInfo | string} file
   * @param {string} newName
   * @returns {Promise} */
  rename(file: FileInfo | string, newName: string): Promise<void> {
    return new Promise((resolve, reject) => api.grok_Dapi_UserFiles_Rename(file, newName, () => resolve(), (e: any) => reject(e)));
  }

  /** List file
   * @param {FileInfo | string} file
   * @param {boolean} recursive
   * @param {string} searchPattern
   * @returns {Promise<List<FileInfo>>} */
  list(file: FileInfo | string, recursive: boolean, searchPattern: string): Promise<FileInfo[]> {
    return new Promise((resolve, reject) =>
        api.grok_Dapi_UserFiles_List(file, recursive, searchPattern, (data: any) => resolve(toJs(data)), (e: any) => reject(e)));
  }

  /** Read file as string
   * @param {FileInfo | string} file
   * @returns {Promise<String>} */
  readAsText(file: FileInfo | string): Promise<string> {
    return new Promise((resolve, reject) => 
        api.grok_Dapi_UserFiles_ReadAsText(file, (data: string | PromiseLike<string>) => resolve(data), (e: any) => reject(e)));
  }

  /** Read file as bytes
   * @param {FileInfo | string} file
   * @returns {Promise<Uint8Array>} */
  readAsBytes(file: FileInfo | string): Promise<Uint8Array> {
    return new Promise((resolve, reject) =>
        api.grok_Dapi_UserFiles_ReadAsBytes(file, (data: any) => resolve(toJs(data)), (e: any) => reject(e)));
  }

  /** Write file
   * @param {FileInfo | string} file
   * @param {Array<number>} blob
   * @returns {Promise} */
  write(file: FileInfo | string, blob: number[]): Promise<void> {
    return new Promise((resolve, reject) => api.grok_Dapi_UserFiles_Write(file, blob, () => resolve(), (e: any) => reject(e)));
  }

  /** Write file
   * @param {FileInfo | string} file
   * @param {string} data
   * @returns {Promise} */
  writeAsText(file: FileInfo | string, data: string): Promise<void> {
    return new Promise((resolve, reject) => api.grok_Dapi_UserFiles_WriteAsText(file, data, () => resolve(), (e: any) => reject(e)));
  }
}

export class Logger {
  putCallback: any;

  constructor(putCallback?: any) {
    this.putCallback = putCallback;
  }

  /** Saves audit record to Datagrok back-end
   * @param {string} message
   * @param {object} params
   * @param {string} type = 'log'
   * */
  log(message: string, params: object, type?: string): void {
    if (type == null)
      type = 'log';
    let msg = {message: message, params: params, type: type};
    if (this.putCallback != null)
      this.putCallback(msg);

    api.grok_Audit(msg.type, msg.message, toDart(msg.params));
  }
}
