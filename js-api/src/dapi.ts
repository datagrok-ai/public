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
  DomainSchema,
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
import {
  DOMAIN_SYSTEM_COLUMNS,
  DomainAggregateRow,
  DomainAggregateSpec,
  DomainAuditEntry,
  DomainBatchOptions,
  DomainBatchReport,
  DomainDatetimeColumns,
  DomainDeleteReport,
  DomainError,
  DomainFacetKind,
  DomainFacetResultOf,
  DomainFacetsSpec,
  DomainFilter,
  DomainGrant,
  DomainInsertResult,
  DomainNotFoundError,
  DomainOpResult,
  DomainOpResultFor,
  DomainPermission,
  DomainQueryBuilder,
  DomainQuerySpec,
  DomainRestrictError,
  DomainRowInsert,
  DomainSavedFilterInfo,
  DomainTableCapabilities,
  DomainTableClientOptions,
  DomainTableInfo,
  DomainTransactionOp,
  DomainUpdateResult,
  DomainValidationError,
  DomainVersionConflictError,
  domainCall,
  retryOnVersionConflict,
  splitDomainTable,
} from './domains';

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

  /** Domain tables API: generic row CRUD over entity-mapped domain schemas
   * declared by plugin manifests (databases/&lt;schema&gt;/schema.json). */
  get domains(): DomainsDataSource {
    return new DomainsDataSource(api.grok_Dapi_Domains());
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

  /** Server info API endpoint
   *  @type {InfoDataSource} */
  get info(): InfoDataSource {
    return new InfoDataSource(api.grok_Dapi_Info());
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

  /** Returns the configured report email address from admin settings.
   * Used as the default recipient for error reports and admin notifications. */
  async getReportEmail(): Promise<string> {
    return JSON.parse(await api.grok_Dapi_Admin_GetReportEmail(this.dart));
  }

  /**
   * Sends email
   * @param email - message that will be sent using configured SMTP service
   */
  async sendEmail(email: Email): Promise<void> {
    if (email.to.length === 0)
      throw new Error('Recipients list shouldn\'t be empty');
    const fd = new FormData();
    fd.append('subject', email.subject);
    fd.append('to', email.to.join(','));
    if (email.text)
      fd.append('text', email.text);
    if (email.html)
      fd.append('html', email.html);
    if (email.bcc && email.bcc.length)
      fd.append('bcc', email.bcc.join(','));
    const toBlob = (a: EmailAttachment) => a.data instanceof Blob ? a.data : new Blob([a.data as BlobPart], { type: a.contentType ?? 'application/octet-stream' });
    for (const a of email.attachments ?? [])
      fd.append('attachment', toBlob(a), a.name);
    const r = await fetch(`${api.grok_Dapi_Root()}/admin/email`,
        { method: 'POST', body: fd, credentials: 'include' });
    if (!r.ok)
      throw new Error(await r.text());
  }

  /** Sends a message to the specified browser client
   * See also {@link ServerMessageTypes} */
  pushMessage(messageType: string, message: object, sessionIds: string[]): Promise<any> {
    return api.grok_Dapi_Admin_PushMessage(this.dart, messageType, api.grok_JSON_decode(JSON.stringify(message)), toDart(sessionIds));
  }
}

export class InfoDataSource {
  dart: any;
  /** @constructs InfoDataSource */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Returns the latest storage usage snapshot, refreshed hourly on the server. */
  getStorageStats(): Promise<{[key: string]: any}> {
    return api.grok_Dapi_Info_GetStorageStats(this.dart);
  }
}

/** Attachment for an [Email]. */
export interface EmailAttachment {
  /** Filename shown in the recipient's mail client. */
  name: string,
  /** Attachment bytes. Blob is streamed natively; Uint8Array is wrapped in a Blob. */
  data: Blob | Uint8Array,
  /** Optional MIME type override; defaults to Blob.type or 'application/octet-stream'. */
  contentType?: string,
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

  /** Files attached to the message. */
  attachments?: EmailAttachment[],
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

  /** Returns all groups the current user belongs to, including transitive parent groups. */
  async currentUserGroups(): Promise<Group[]> {
    return toJs(await api.grok_Dapi_Get_CurrentUserGroups());
  }

  /** Requests that `requester` be added as a member of `group`.
   *  An admin of `group` must approve before the membership takes effect. */
  async requestMembership(group: Group, requester: Group): Promise<void> {
    return api.grok_GroupsDataSource_RequestMembership(group.id, requester.id);
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

  /** Returns entities favorited by the specified group (or the current user's group if not specified) */
  getFavorites(group?: Group): Promise<Entity[]> {
    return toJs(api.grok_EntitiesDataSource_GetFavorites(this.dart, toDart(group)));
  }

  /** Returns favorites for multiple groups in a single round trip; object keyed by group id. */
  async getFavoritesForGroups(groups: Group[]): Promise<{[gid: string]: Entity[]}> {
    const obj: {[k: string]: any[]} = await api.grok_EntitiesDataSource_GetFavoritesForGroups(this.dart, groups.map((g) => toDart(g)));
    for (const gid of Object.keys(obj))
      obj[gid] = obj[gid].map((e) => toJs(e));
    return obj;
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

  async getSchemas(e: DataConnection, catalog: string | null = null): Promise<string[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Schemas(this.dart, e.dart, catalog));
  }

  async getSchema(e: DataConnection, schemaName: string | null = null, tableName: string | null = null, catalog: string | null = null): Promise<TableInfo[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Schema(this.dart, e.dart, schemaName, tableName, catalog));
  }

  async getUniqueColumnsNames(c: DataConnection, schema: string, table: string): Promise<string[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Unique_Columns(this.dart, c.dart, schema, table));
  }


  async getDatabaseInfo(c: DataConnection, catalog: string | null = null): Promise<DbInfo[]> {
    return toJs(await api.grok_DataConnectionsDataSource_Get_Db_Info(this.dart, c.dart, catalog))
  }

  /** Initiates the OAuth consent flow for a connection whose auth method is 'OAuth'. */
  async requestOAuthConsent(c: DataConnection): Promise<void> {
    await api.grok_DataConnectionsDataSource_RequestOAuthConsent(c.dart);
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
 * Generic client for domain tables — entity-mapped PostgreSQL schemas that plugins
 * declare via `databases/<schema>/schema.json` manifests. Rows are plain JSON objects;
 * column metadata and security come from the server-side schema registry. */
export class DomainsDataSource {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Registered domain schemas; use `.include('tables')` to populate their tables. */
  get schemas(): HttpDataSource<DomainSchema> {
    return new HttpDataSource(api.grok_Dapi_Domains_Schemas(this.dart));
  }

  /** Returns a client for the domain table addressed as `'<schema>.<table>'`, e.g. `'plates.plate'`.
   * Pass a row interface for a typed client: `grok.dapi.domains.table<PlateRow>('plates.plate')`;
   * pass an insert interface as the second generic to also gate {@link DomainTableClient.insert},
   * a column-name union as the third to compile-check filters/columns/sorts, and an expand map
   * as the fourth to compile-check expand keys (`grok api`-generated clients pass all four). */
  table<TRow = any, TInsert = DomainRowInsert<TRow>, TColumn extends string = string,
      TExpand extends {[key: string]: {}} = {[key: string]: {}}>(
    name: string, options?: DomainTableClientOptions): DomainTableClient<TRow, TInsert, TColumn, TExpand> {
    const [schema, table] = splitDomainTable(name);
    return new DomainTableClient<TRow, TInsert, TColumn, TExpand>(this.dart, schema, table, options);
  }

  /** Executes ordered [ops] (≤1000) atomically within domain schema [schema]; any failure rolls the
   * whole transaction back (the error carries `opIndex`). Resolves to ordered per-op results
   * shaped like the corresponding single-op endpoints — a tuple ops literal gets per-op
   * result types (`[{op: 'update', ...}, {op: 'insert', ...}]` →
   * `[DomainUpdateResult, DomainInsertResult]`, no positional casts). */
  transaction<T extends DomainTransactionOp[]>(schema: string, ops: [...T]):
      Promise<{[K in keyof T]: DomainOpResultFor<T[K]>}> {
    return domainCall(api.grok_Dapi_Domains_Transaction(this.dart, schema, ops));
  }

  /** Creates an empty user-managed domain schema (physical PG schema 'usr_<name>');
   * requires the CreateDomainSchema privilege. Add tables via schema(name).apply(). */
  createSchema(name: string, options?: {friendlyName?: string; description?: string}): Promise<{[key: string]: any}> {
    return domainCall(api.grok_Dapi_Domains_CreateSchema(this.dart, name,
      options?.friendlyName ?? null, options?.description ?? null));
  }

  /** Lifecycle handle for a registered domain schema (mirrors table()). */
  schema(name: string): DomainSchemaClient {
    return new DomainSchemaClient(this.dart, name);
  }

  /** Reflection over the domain registry: the runtime {@link Property} metadata,
   * table info with FK-inverted child tables, and batched display-name resolution —
   * what reflective UI components consume to work on any table without codegen. */
  get registry(): DomainRegistryClient {
    return new DomainRegistryClient();
  }

  /** Drops the client-side domain UI caches (table capabilities, row permissions,
   * writable columns, resolved display names) so the next read re-probes server
   * truth. Grant changes made through this client invalidate automatically; call
   * this after out-of-band grant changes (another session, `grok s`, server-side). */
  invalidateUiCaches(): void {
    api.grok_Domains_InvalidateUiCaches();
  }
}

/**
 * Reflection over the client-side domain registry (see {@link DomainsDataSource.registry}):
 * schemas → tables → columns → {@link Property} metadata, loaded once per session and
 * shared with the built-in domain UI. Tables are addressed as `'<schema>.<table>'`. */
export class DomainRegistryClient {
  /** {@link Property} objects describing the table's declared columns — the same
   * runtime metadata (type, semType, choices, min/max, nullable, defaultValue) that
   * drives the built-in row editor and server-side validation, get/set-bound to
   * {@link DomainRow} values. `friendlyName` carries the label the platform renders
   * ('Project' for a `project_id` reference), so reflective forms label fields like
   * the built-in surfaces do. Treat the objects as read-only — they are shared with
   * the client-side registry. Rejects with a {@link DomainValidationError} for
   * unknown tables. */
  async rowProperties(table: string): Promise<Property[]> {
    splitDomainTable(table);
    return toJs(await domainCall(api.grok_DomainRegistry_RowProperties(table)));
  }

  /** Registry metadata of one table: display identity (nameColumn, businessKey,
   * singular/plural names), the schema.json description, security mode, audit flag,
   * and the FK-inverted {@link DomainChildTableRef} list that drives detail-table
   * links. */
  tableInfo(table: string): Promise<DomainTableInfo> {
    const [schema, t] = splitDomainTable(table);
    return domainCall(api.grok_DomainRegistry_TableInfo(schema, t));
  }

  /** Batched display-name resolution for row [ids] of [table] (the standard display
   * identity: name-column value → dash-joined business key → id). Requests coalesce
   * client-side into one narrow fetch per table, and every id is cached afterwards —
   * cheap to call per cell/render. EVERY requested id is present as a key; the
   * unresolvable ones (invisible, unknown, nameless) map to `null` rather than being
   * omitted, so `ids.map((id) => names[id])` is always aligned. */
  resolveNames(table: string, ids: string[]): Promise<{[id: string]: string | null}> {
    const [schema, t] = splitDomainTable(table);
    return domainCall(api.grok_Domains_ResolveNames(schema, t, ids));
  }
}

/**
 * Lifecycle handle for one registered domain schema (see {@link DomainsDataSource.schema}):
 * manifest export, partial applies with dry-run change plans, whole-schema audit, and full
 * purge. Mutations apply to user-managed schemas only (package schemas deploy on publish). */
export class DomainSchemaClient {
  dart: any;

  constructor(dart: any, public readonly name: string) {
    this.dart = dart;
  }

  /** Full manifest reconstructed from the registry (feeds editors; doubles as export). */
  manifest(): Promise<{[key: string]: any}> {
    return domainCall(api.grok_Dapi_Domains_GetManifest(this.dart, this.name));
  }

  /** Applies a partial manifest; dryRun returns the change plan. Named tables replace
   * their current definition wholesale — everything omitted comes from the registry.
   * Destructive plans require confirmDestructive; stale ifVersion → DomainVersionConflictError;
   * a destructive plan without confirmation → DomainError code 'destructive-confirmation-required'
   * (the plan rides in error.body.plan).
   *
   * On a USER-managed schema this requires Edit and `ifVersion` tracks the schema's apply
   * counter. On a PACKAGE-managed schema it is the user-extension path: it requires
   * 'Extend' on the schema entity, `ifVersion` tracks `ext_version` (the plan echoes it as
   * `extVersion`) while the plugin's own version stays frozen, and the writable surface is
   * limited to objects you own — your own tables through `tables` (needs the schema's
   * `extensible.tables` opt-in) and your own columns on plugin tables through `extend`
   * (needs that table's `extensible` opt-in). `extend` is full state for YOUR columns of
   * that table: omitting one you added proposes its drop. Plugin objects are immutable —
   * touching them is a 400 ('plugin-table-not-editable', 'not-extensible',
   * 'schema-not-extensible', 'invalid-extend', 'plugin-schema-not-editable'), and columns
   * you add to a plugin table must stay nullable, non-unique, non-display, and must not
   * cascade deletes. Physical names of extension columns are prefixed server-side; every
   * API surface — reads, writes, filters, and the audit trail — keeps using the logical name,
   * with the exception of raw PostgreSQL error text that quotes an identifier verbatim
   * (a constraint violation naming `fk_<table>_x_<column>`, say). */
  apply(body: {tables?: object; extend?: {[table: string]: {columns: {[column: string]: object}}};
               propertySchemas?: object; dropTables?: string[];
               ifVersion?: string; confirmDestructive?: boolean},
        options?: {dryRun?: boolean}): Promise<{[key: string]: any}> {
    return domainCall(api.grok_Dapi_Domains_ApplySchema(this.dart, this.name, body, options?.dryRun ?? false));
  }

  /** Whole-schema history: all tables' events + 'ddl' registry events, newest first. */
  audit(options?: {limit?: number}): Promise<DomainAuditEntry[]> {
    return domainCall(api.grok_Dapi_Domains_SchemaAudit(this.dart, this.name, options?.limit ?? null));
  }

  /** Fully purges a user-managed schema: data, audit partition, registry, entities, permissions. */
  delete(): Promise<void> {
    return domainCall(api.grok_Dapi_Domains_DeleteSchema(this.dart, this.name));
  }

  /** Direct permission rows on this schema's registry entity (schema-level operation
   * rights, not row-data access — see {@link grant}). Requires Share. */
  grants(): Promise<DomainGrant[]> {
    return domainCall(api.grok_Dapi_Domains_SchemaGrants(this.dart, this.name));
  }

  /** Idempotently grants [permission] on this schema's registry entity to [group] (a group
   * id). Schema grants gate schema-level operations (apply requires Edit, delete requires
   * Delete, sharing requires Share, extending a package schema requires Extend). They do
   * NOT grant access to row data — use `table('s.t').grant()` per table for that.
   * 'Extend' is meaningful on schema entities only: granting it on a table or column
   * schema is a 400. Requires Share. */
  async grant(group: string, permission: DomainPermission): Promise<void> {
    await domainCall(api.grok_Dapi_Domains_SchemaGrant(this.dart, this.name, group, permission));
    api.grok_Domains_InvalidateUiCaches();
  }

  /** Revokes [permission] (or every schema permission when omitted) from [group].
   * Requires Share. */
  async revoke(group: string, permission?: DomainPermission): Promise<void> {
    await domainCall(api.grok_Dapi_Domains_SchemaRevoke(this.dart, this.name, group, permission ?? null));
    api.grok_Domains_InvalidateUiCaches();
  }
}

/**
 * Row CRUD for one domain table. Reads return only rows and columns the current
 * user can see; writes are validated, permission-checked, and audited server-side.
 * Pass a row interface as `TRow` for typed reads/writes, and an insert interface
 * as `TInsert` so {@link insert} enforces required columns (`grok api`-generated
 * clients pass both; see {@link DomainsDataSource.table}). */
export class DomainTableClient<TRow = any, TInsert = DomainRowInsert<TRow>,
    TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}> {
  dart: any;
  private readonly _datetimeOverride: DomainDatetimeColumns | null;

  /** Registry-resolved datetime columns per `'<schema>.<table>'`, shared across client
   * instances (the generated client map creates a fresh client per property access).
   * Holds promises so concurrent resolutions coalesce; failures are not cached. */
  private static _datetimeCache = new Map<string, Promise<DomainDatetimeColumns>>();

  constructor(dart: any, public readonly schema: string, public readonly table: string,
              options?: DomainTableClientOptions) {
    this.dart = dart;
    this._datetimeOverride = options?.datetimeColumns != null || options?.detailDatetimeColumns != null
      ? {own: options.datetimeColumns ?? [], details: options.detailDatetimeColumns ?? {}} : null;
  }

  /**
   * Datetime columns of this table and of its `'details:'` child tables — what JSON
   * reads materialize as dayjs. Resolved from the domain registry (the same column
   * metadata that types them in the first place) and cached per table for the session;
   * explicit {@link DomainTableClientOptions} act as an override and skip the registry.
   * A failed resolution (unregistered table, no session) leaves that read unconverted
   * and is retried on the next one.
   */
  private _datetimes(): Promise<DomainDatetimeColumns> {
    if (this._datetimeOverride != null)
      return Promise.resolve(this._datetimeOverride);
    const address = `${this.schema}.${this.table}`;
    let cached = DomainTableClient._datetimeCache.get(address);
    if (cached == null) {
      cached = (async () => {
        const registry = new DomainRegistryClient();
        // The registry describes DECLARED columns only; the system datetimes exist
        // on every domain table and are prepended by hand.
        const system = ['created_on', 'updated_on'];
        const datetimeNames = (props: Property[]) =>
          props.filter((p) => p.propertyType === 'datetime').map((p) => p.name);
        const props = await registry.rowProperties(address);
        const own = [...system, ...datetimeNames(props)];
        // Master-expand fields ('<fk_column>.<column>') of declared datetime columns:
        // a ref column's semType is its target's '<schema>.<table>' address.
        for (const p of props)
          if (/^\w+\.\w+$/.test(p.semType ?? '')) {
            try {
              for (const c of datetimeNames(await registry.rowProperties(p.semType)))
                own.push(`${p.name}.${c}`);
            } catch (_) { /* an unreadable target leaves its fields unconverted */ }
          }
        const details: {[detailField: string]: string[]} = {};
        for (const child of (await registry.tableInfo(address)).childTables ?? []) {
          try {
            details[child.table] =
              [...system, ...datetimeNames(await registry.rowProperties(`${this.schema}.${child.table}`))];
          } catch (_) { /* an unreadable child leaves its rows unconverted */ }
        }
        return {own, details};
      })();
      DomainTableClient._datetimeCache.set(address, cached);
      cached.catch(() => DomainTableClient._datetimeCache.delete(address));
    }
    return cached.catch(() => ({own: [], details: {}}));
  }

  /** Materializes [datetimes] as dayjs on one JSON row, recursing into the
   * `'details:'` child arrays. */
  private _fromWire(row: any, datetimes: DomainDatetimeColumns): any {
    if (row == null)
      return row;
    for (const c of datetimes.own)
      if (typeof row[c] === 'string')
        row[c] = dayjs(row[c]);
    for (const field of Object.keys(datetimes.details)) {
      const children = row[field];
      if (Array.isArray(children))
        for (const child of children)
          for (const c of datetimes.details[field])
            if (child != null && typeof child[c] === 'string')
              child[c] = dayjs(child[c]);
    }
    return row;
  }

  /** Bare `query()` returns an awaitable {@link DomainQueryBuilder} (it used to resolve
   * all-defaults rows — `await table.query()` behaves identically, and everything chains:
   * `await table.query().where('sku', '=', key).orderBy('created_on', true).top(5)`).
   * Prefer the builder's condition forms and the `cond`/`and`/`or` helpers over
   * template-built filter strings — condition values are bound server-side, so any string
   * value is safe (apostrophes included). */
  query(): DomainQueryBuilder<TRow, TColumn, TExpand, TRow, DataFrame>;
  /** Runs a filtered, sorted, paginated query; resolves to an array of row objects (10k row cap). */
  query(spec: DomainQuerySpec<TColumn, keyof TExpand & string>): Promise<TRow[]>;
  query(spec?: DomainQuerySpec<TColumn, keyof TExpand & string>): any {
    if (spec === undefined)
      return new DomainQueryBuilder<TRow, TColumn, TExpand, TRow, DataFrame>(this);
    const datetimes = this._datetimes();
    return domainCall(api.grok_Dapi_Domains_Query(this.dart, this.schema, this.table, spec))
      .then(async (rows: any) => {
        const dt = await datetimes;
        return rows.map((r: any) => this._fromWire(r, dt));
      });
  }

  /** Runs the same query as {@link query} but resolves to a typed DataFrame (d42 wire format,
   * 10M row cap). Columns carry the db property tags (`dbPropertySchema`/`dbPropertyName`),
   * `.choices`, and semantic types; system columns are untagged. `'details:'` expand is
   * JSON-only — use {@link query}; master expand yields flat `'<fk_column>.<name>'` columns. */
  queryDf(spec: DomainQuerySpec<TColumn, keyof TExpand & string> = {}): Promise<DataFrame> {
    return domainCall(api.grok_Dapi_Domains_QueryDf(this.dart, this.schema, this.table, spec));
  }

  /** Grouped aggregation over the rows and columns visible to the caller (10k row cap);
   * resolves to result rows named by group column / measure alias. Alias measures with `as`
   * (and pass literal groupBy) to get typed result keys; without them, cast or use `aggregateDf`. */
  aggregate<TGroup extends string = never, TAlias extends string = never>(
    spec: DomainAggregateSpec<TColumn, TGroup, TAlias>): Promise<DomainAggregateRow<TGroup | TAlias>[]> {
    return domainCall(api.grok_Dapi_Domains_Aggregate(this.dart, this.schema, this.table, spec));
  }

  /** Fetches one row by id; resolves to null if the row does not exist or is not visible
   * (typed `TRow` for backward compatibility — guard against null, or use `first`). */
  async get(id: string): Promise<TRow> {
    const datetimes = this._datetimes();
    const row = await domainCall(api.grok_Dapi_Domains_GetRow(this.dart, this.schema, this.table, id));
    return this._fromWire(row, await datetimes);
  }

  /** Inserts a single row or a small array of rows; resolves to per-row reports
   * (`{id, created}`, or `{status: 'duplicate', existingId}` on a business-key match —
   * pass `options.errorOnDuplicate` to reject duplicates with a
   * {@link DomainValidationError} instead, its `isDuplicate` set).
   * For tables that declare `"idempotency": true`, pass an `idempotencyKey` (UUID) row field
   * to make retries safe: a replay returns the existing id with `status: 'idempotent-replay'`. */
  insert(rows: TInsert | TInsert[], options?: {errorOnDuplicate?: boolean}): Promise<DomainInsertResult[]> {
    return domainCall(api.grok_Dapi_Domains_Insert(this.dart, this.schema, this.table, rows,
      options?.errorOnDuplicate ?? false));
  }

  /** Partially updates a row; pass `options.version` (the version the client last read) for
   * optimistic concurrency — the update fails with a {@link DomainVersionConflictError} if the
   * row has changed since. Resolves to `{id, version}` (version increments on every update). */
  update(id: string, values: Partial<TRow>, options?: {version?: number}): Promise<DomainUpdateResult> {
    return domainCall(api.grok_Dapi_Domains_Patch(this.dart, this.schema, this.table, id, values, options?.version));
  }

  /** Bulk upload: a DataFrame (sent as d42), a CSV string, an array of row objects, or raw
   * bytes (`options.format`: `'d42'` default, `'parquet'` converted via the Arrow package).
   * `options.mode: 'upsert'` merges by the table's business key. Resolves to the batch report;
   * a failure that carries the per-row report (e.g. an allOrNothing abort) resolves with
   * `error` set, report-less failures reject. */
  batch(data: DataFrame | string | object[] | Uint8Array, options: DomainBatchOptions = {}): Promise<DomainBatchReport> {
    const format = data instanceof DataFrame ? 'df' : typeof data === 'string' ? 'csv' :
      data instanceof Uint8Array ? (options.format ?? 'd42') : 'json';
    return domainCall(api.grok_Dapi_Domains_Batch(this.dart, this.schema, this.table,
      data instanceof DataFrame ? data.dart : data, format, options));
  }

  /** Soft-deletes a row (engine-enforced cascade/restrict/setnull for declared relations;
   * a restrict reference rejects with a {@link DomainRestrictError}). */
  delete(id: string): Promise<void> {
    return domainCall(api.grok_Dapi_Domains_Delete(this.dart, this.schema, this.table, id));
  }

  /** Soft-deletes up to `options.limit` (≤1000, default 1000) matching rows you may delete,
   * oldest first, in ONE transaction; referential actions apply per row and a restrict
   * reference rejects the whole call ({@link DomainRestrictError} — nothing is deleted).
   * The filter is required — an empty one rejects with a {@link DomainValidationError}.
   * Each row runs through the per-row engine (cascade fan-out, audit), so this is not a
   * free bulk sweep: prefer a narrow filter and a modest limit, and loop while `hasMore`
   * for larger sets. A row already gone by its turn (deleted concurrently, or eaten by an
   * earlier row's cascade) is skipped, not an error. */
  deleteWhere(filter: DomainFilter<TColumn>, options?: {limit?: number}): Promise<DomainDeleteReport> {
    return domainCall(api.grok_Dapi_Domains_DeleteWhere(this.dart, this.schema, this.table,
      filter, options?.limit ?? null));
  }

  /** Creates the entities row for a domain row so it can be individually shared. */
  promote(id: string): Promise<{id: string; promoted: boolean}> {
    return domainCall(api.grok_Dapi_Domains_Promote(this.dart, this.schema, this.table, id));
  }

  /** Returns the row's audit trail (before/after diffs, in-transaction with each write). */
  audit(id: string): Promise<DomainAuditEntry[]> {
    return domainCall(api.grok_Dapi_Domains_RowAudit(this.dart, this.schema, this.table, id));
  }

  /** Batched facet computation for filter panels: category counts, histograms, min/max, row count,
   * and column profiling in one round trip — counts under all other filters, bounds under the row
   * predicate only (the stable-axis exception; see {@link DomainFacetsSpec}). All results respect
   * the row predicate and column security. Resolves to `{facets: {<id>: <result>}}` —
   * `'categories'` results as `{categories: DomainFacetCategory[], hasMore}`, `'histogram'` as
   * `{min, max, buckets, totalBuckets, nulls}` — `buckets` counted under the other filters,
   * `totalBuckets` under the row predicate only (datetime bounds are ISO-8601 strings), `'minMax'` as
   * `{min, max}`, `'count'` as `{count}`, `'plan'` as `{columns: [{name, distinct, min?, max?}]}`. */
  facets<TId extends string = string, TKind extends DomainFacetKind = DomainFacetKind>(
    spec: DomainFacetsSpec<TColumn, TId, TKind>): Promise<{facets: {[K in TId]: DomainFacetResultOf<TKind>}}> {
    return domainCall(api.grok_Dapi_Domains_Facets(this.dart, this.schema, this.table, spec));
  }

  /** Row count under [filter] (condition tree or smart string; omit for the whole table). */
  count(filter?: DomainFilter<TColumn>): Promise<number> {
    return domainCall(api.grok_Dapi_Domains_Count(this.dart, this.schema, this.table, filter ?? null));
  }

  /** True when at least one visible row matches [filter]. */
  async exists(filter?: DomainFilter<TColumn>): Promise<boolean> {
    return (await this.count(filter)) > 0;
  }

  /** First matching row or null; shorthand for query({...spec, limit: 1}). */
  async first(spec?: DomainQuerySpec<TColumn, keyof TExpand & string>): Promise<TRow | null> {
    const rows = await this.query({...spec, limit: 1});
    return rows.length === 0 ? null : rows[0];
  }

  /** Business-key (or any equality-set) lookup; ambiguous or absent → null. */
  async getByKey(keyValues: Partial<TRow>): Promise<TRow | null> {
    const datetimes = this._datetimes();
    const row = await domainCall(api.grok_Dapi_Domains_GetByKey(this.dart, this.schema, this.table, keyValues));
    return this._fromWire(row, await datetimes);
  }

  /** Rows for [ids] as a typed DataFrame: the 'id' column plus [fields] (default: all
   * visible columns). Chunked client-side at 100k ids; row predicate + column security
   * apply (missing/invisible ids are absent rows). An empty [ids] list short-circuits to
   * an empty ZERO-COLUMN frame — even when [fields] are requested. */
  fetchFields(ids: string[], fields?: TColumn[]): Promise<DataFrame> {
    return domainCall(api.grok_Dapi_Domains_FetchFields(this.dart, this.schema, this.table, ids, fields ?? null));
  }

  /** {@link aggregate} returning a typed d42 DataFrame (10k row cap, both formats). */
  aggregateDf(spec: DomainAggregateSpec<TColumn>): Promise<DataFrame> {
    return domainCall(api.grok_Dapi_Domains_AggregateDf(this.dart, this.schema, this.table, spec));
  }

  /** Inserts or merges ONE row by the table's business key (requires businessKey).
   * Rides the batch engine in upsert mode; failures (invalid values, missing
   * business-key column) reject with a {@link DomainValidationError}. */
  upsert(row: TInsert): Promise<{id: string; status: 'inserted' | 'updated'}> {
    return domainCall(api.grok_Dapi_Domains_Upsert(this.dart, this.schema, this.table, row));
  }

  /** Table-wide audit trail, newest first; limit clamps to [1, 1000]. */
  auditLog(options?: {limit?: number}): Promise<DomainAuditEntry[]> {
    return domainCall(api.grok_Dapi_Domains_TableAudit(this.dart, this.schema, this.table, options?.limit ?? null));
  }

  /** Subscribes the current user to change notifications for the table (or one row when
   * [id] is given; row watch requires the table's audit trail). Resolves to whether the
   * server confirmed the subscription. */
  watch(id?: string): Promise<boolean> {
    return domainCall(api.grok_Dapi_Domains_Watch(this.dart, this.schema, this.table, id ?? null));
  }

  /** Removes the table (or row, when [id] is given) subscription; resolves to whether the
   * server confirmed the removal. */
  unwatch(id?: string): Promise<boolean> {
    return domainCall(api.grok_Dapi_Domains_Unwatch(this.dart, this.schema, this.table, id ?? null));
  }

  /** Whether the current user watches the table (or row, when [id] is given). */
  isWatching(id?: string): Promise<boolean> {
    return domainCall(api.grok_Dapi_Domains_IsWatching(this.dart, this.schema, this.table, id ?? null));
  }

  /** Effective {@link DomainTableCapabilities} of the CURRENT user on this table:
   * server-truth permission probes on the final securing entity plus the writable-column
   * mirror of column security. Cached per registry generation + user; grant changes made
   * through this client drop the cache automatically, out-of-band changes require
   * {@link DomainsDataSource.invalidateUiCaches}. Rejects with a
   * {@link DomainValidationError} for unknown tables. */
  capabilities(): Promise<DomainTableCapabilities> {
    return domainCall(api.grok_Domains_TableCapabilities(this.schema, this.table));
  }

  /** Direct permission rows on this table's registry entity. Requires Share. */
  grants(): Promise<DomainGrant[]> {
    return domainCall(api.grok_Dapi_Domains_TableGrants(this.dart, this.schema, this.table));
  }

  /** Idempotently grants [permission] on this table to [group] (a group id). Requires Share. */
  async grant(group: string, permission: DomainPermission): Promise<void> {
    await domainCall(api.grok_Dapi_Domains_TableGrant(this.dart, this.schema, this.table, group, permission));
    api.grok_Domains_InvalidateUiCaches();
  }

  /** Revokes [permission] (or all four when omitted) from [group]. Requires Share. */
  async revoke(group: string, permission?: DomainPermission): Promise<void> {
    await domainCall(api.grok_Dapi_Domains_TableRevoke(this.dart, this.schema, this.table,
      group, permission ?? null));
    api.grok_Domains_InvalidateUiCaches();
  }

  /** Restricts [column] to its own single-column property schema and grants [group] View on
   * it — afterwards only grantees (and admins) see the column. Returns the per-column
   * schema {id, name} (a further grants target). Requires Share on the core schema;
   * jsonb columns are managed via their property schema (DomainError). Out of scope on
   * this surface: listing a per-column schema's own grants (keep the returned id), and
   * property-schema grants in general (server-supported; no JS surface yet). */
  async shareColumn(column: TColumn, group: string, permission?: DomainPermission): Promise<{id: string; name: string}> {
    const res = await domainCall(api.grok_Dapi_Domains_ShareColumn(this.dart, this.schema, this.table,
      column, group, permission ?? 'View'));
    api.grok_Domains_InvalidateUiCaches();
    return res;
  }

  /** Restricts [column] without granting anyone (inverse: {@link restoreColumnVisibility}). */
  async restrictColumn(column: TColumn): Promise<{id: string; name: string}> {
    const res = await domainCall(api.grok_Dapi_Domains_RestrictColumn(this.dart, this.schema, this.table, column));
    api.grok_Domains_InvalidateUiCaches();
    return res;
  }

  /** Deletes the per-column schema with its grants; the column rejoins the
   * everyone-visible core schema. */
  async restoreColumnVisibility(column: TColumn): Promise<void> {
    await domainCall(api.grok_Dapi_Domains_RestoreColumnVisibility(this.dart, this.schema, this.table, column));
    api.grok_Domains_InvalidateUiCaches();
  }

  /** Insert-or-update by row identity: no id → insert; id → version-checked partial update
   * (uses row.version when present). System columns other than the addressing id/version are
   * stripped from the payload, so spreading a read row is safe. Resolves to the row merged
   * with the new id/version — NB: re-saving an OLD object after a successful save carries its
   * stale version and rejects with a {@link DomainVersionConflictError}; keep the resolved
   * row. An id-less save whose business key matches an existing row is a TRUE
   * insert-or-update: the values are applied to the existing row as a versioned update
   * (retried on conflict), so save() never resolves a version-less row. An idempotency-key
   * replay applies nothing — the original insert already did — and resolves the existing
   * row's fresh version (rejects {@link DomainNotFoundError} when that row is no longer
   * visible). The unversioned
   * (last-write-wins) update remains available ONLY by deliberately constructing
   * `{id}` without a version — never as a side effect of a duplicate. */
  async save(row: Partial<TRow> & {id?: string; version?: number}): Promise<TRow> {
    const values: any = {};
    for (const k of Object.keys(row))
      if (!(DOMAIN_SYSTEM_COLUMNS as readonly string[]).includes(k))
        values[k] = (row as any)[k];
    if (row.id == null) {
      const [res] = await this.insert(values);
      if (res.status === 'duplicate') {
        // The duplicate insert wrote NOTHING — land the caller's values on the
        // existing row under the optimistic check (fresh read + retry).
        const upd = await this.updateWithRetry(res.id, () => values);
        return {...(row as any), id: res.id, version: upd!.version};
      }
      if (res.status === 'idempotent-replay') {
        // The ORIGINAL insert already applied these values; report its state.
        const fresh = await this.get(res.id);
        if (fresh == null)
          throw new DomainNotFoundError(`Row "${res.id}" not found in ${this.schema}.${this.table}`,
            404, {error: 'not-found', id: res.id});
        return {...(row as any), id: res.id, version: (fresh as any).version};
      }
      return {...(row as any), id: res.id, version: res.version};
    }
    const res = await this.update(row.id, values, row.version != null ? {version: row.version} : undefined);
    return {...(row as any), id: res.id, version: res.version};
  }

  /** Read-modify-write with optimistic retry: fetches the fresh row, applies [mutate], writes
   * with the fresh version; retries on DomainVersionConflictError (`maxRetries` counts
   * retries after the initial attempt — default 5 retries = up to 6 attempts, no backoff).
   * [mutate] returning null skips the write (resolves null). Rejects
   * {@link DomainNotFoundError} when the row is invisible/absent, and {@link DomainError}
   * (code 'no-version') when the fresh row carries no version — the write never silently
   * degrades to an unversioned (last-write-wins) update. For multi-op flows
   * (e.g. a guarded transaction), use `DG.retryOnVersionConflict` directly with the fresh
   * read inside the action. */
  updateWithRetry(id: string, mutate: (fresh: TRow) => Partial<TRow> | null,
      options?: {maxRetries?: number}): Promise<DomainUpdateResult | null> {
    return retryOnVersionConflict(async () => {
      const fresh = await this.get(id);
      if (fresh == null)
        throw new DomainNotFoundError(`Row "${id}" not found in ${this.schema}.${this.table}`,
          404, {error: 'not-found', id: id});
      const version = (fresh as any).version;
      if (version == null)
        throw new DomainError(`updateWithRetry: row "${id}" carries no version — cannot update optimistically`,
          0, {error: 'no-version'});
      const values = mutate(fresh);
      if (values == null)
        return null;
      return await this.update(id, values, {version});
    }, options);
  }

  /** Saved filter presets of this table — shareable entities carrying filter panel states. */
  get filters(): DomainSavedFiltersClient {
    return new DomainSavedFiltersClient(this.dart, this.schema, this.table);
  }
}

/**
 * Saved filter presets of one domain table (see {@link DomainTableClient.filters}): small
 * shareable entities carrying the filter panel's state maps. Sharing uses the standard entity
 * permission machinery ({@link Dapi.permissions}); the server scopes {@link list} to what the
 * caller can see. */
export class DomainSavedFiltersClient {
  dart: any;

  constructor(dart: any, public readonly schema: string, public readonly table: string) {
    this.dart = dart;
  }

  /** Presets of this table visible to the current user, ordered by name. */
  list(): Promise<DomainSavedFilterInfo[]> {
    return domainCall(api.grok_Dapi_Domains_ListFilters(this.dart, this.schema, this.table));
  }

  /** Saves a preset for the caller; pass `options.id` to update an existing preset in place
   * (the original author is preserved). Resolves to the saved preset. */
  save(name: string, states: {[column: string]: any}, options?: {id?: string}): Promise<DomainSavedFilterInfo> {
    return domainCall(api.grok_Dapi_Domains_SaveFilter(this.dart, this.schema, this.table, name, states, options?.id));
  }

  /** Deletes a preset (requires the entity Delete permission). */
  delete(id: string): Promise<void> {
    return domainCall(api.grok_Dapi_Domains_DeleteFilter(this.dart, id));
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
   * Revalidates Docker image (checks that image exists and is pullable).
   * @param imageId - ID of the {@link DockerImage} to revalidate.
   * @returns {Promise<void>} - promise that resolves with void or throws Exception if something went wrong.
   */
  revalidate(imageId: string): Promise<void> {
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

  async markAllAsRead(): Promise<void> {
    return api.grok_Dapi_Notifications_MarkAllAsRead();
  }

  async markAsRead(notificationId: string): Promise<void> {
    return api.grok_Dapi_Notifications_MarkAsRead(notificationId);
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
