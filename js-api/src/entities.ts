import {ColumnType, ScriptLanguage, SemType, Type, TYPE} from "./const";
import { FuncCall } from "./functions";
import {toJs} from "./wrappers";
import {FileSource} from "./dapi";
import {MapProxy} from "./utils";

declare var grok: any;
let api = <any>window;

type PropertyGetter = (a: object) => any;
type PropertySetter = (a: object, value: any) => void;

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

  newId(): void { api.grok_Entity_New_Id(this.dart);}

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
  get createdOn(): string { return api.grok_Entity_Get_CreatedOn(this.dart); }

  /** Time when entity was updated **/
  get updatedOn(): string { return api.grok_Entity_Get_UpdatedOn(this.dart); }

  /** Who created entity **/
  get author(): User { return toJs(api.grok_Entity_Get_Author(this.dart)); }

  /** Entity properties */
  getProperties(): Promise<Map<string, any>> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.dart, this.dart, (p: any) => resolve(p), (e: any) => reject(e)));
  }

  /** Sets entity properties */
  setProperties(props: Map<string, any>): Promise<any> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.dart, this.dart, props, (_: any) => resolve(_), (e: any) => reject(e)));
  }

  /** Returns a string representing the object */
  toString(): string { return api.grok_Object_ToString(this.dart); }

  hasTag(tag: string): boolean {
    return api.grok_Entity_Has_Tag(this.dart, tag);
  }

  tag(tag: string): boolean {
    return api.grok_Entity_Tag(this.dart, tag);
  }

  unTag(tag: string): boolean {
    return api.grok_Entity_UnTag(this.dart, tag);
  }
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

  static create(): User {return new User(api.grok_User_From_Id(null)); };
  static fromId(id: string): User { return new User(api.grok_User_From_Id(id)); }

  static current(): User { return new User(api.grok_User()); }

  /** First name */
  get firstName(): string { return api.grok_User_Get_FirstName(this.dart); }
  set firstName(name: string) {api.grok_User_Set_FirstName(this.dart, name);}

    /** Last name */
  get lastName(): string { return api.grok_User_Get_LastName(this.dart); }
  set lastName(name: string) {api.grok_User_Set_LastName(this.dart, name);}

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
  set login(login: string) {api.grok_User_Set_Login(this.dart, login);}

  toMarkup(): string { return api.grok_User_ToMarkup(this.dart); }

  /** Security Group */
  get group(): Group { return toJs(api.grok_User_Get_Group(this.dart)); }
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

  get type(): string { return api.grok_UserSession_Get_Type(this.dart); }

  /** External Token */
  get externalToken(): string { return api.grok_UserSession_Get_ExternalToken(this.dart); }

  /** User */
  get user(): User { return toJs(api.grok_UserSession_Get_User(this.dart)); }
}

/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/overview/functions/function}
 * */
export class Func extends Entity {
  public aux: any;
  public options: any;

  constructor(dart: any) {
    super(dart);
    this.aux = new MapProxy(api.grok_Func_Get_Aux(this.dart));
    this.options = new MapProxy(api.grok_Func_Get_Options(this.dart));
  }

  get description(): string { return api.grok_Func_Get_Description(this.dart); }

  get type(): string { return api.grok_Func_Get_Type(this.dart); }

  get path(): string { return api.grok_Func_Get_Path(this.dart); }

  get helpUrl(): string { return api.grok_Func_Get_HelpUrl(this.dart); }

  get package(): Package { return api.grok_Func_Get_Package(this.dart); }

  /** Returns {@link FuncCall} object in a stand-by state
   * @param {object} parameters
   * @returns {FuncCall} */
  prepare(parameters: object = {}): FuncCall {
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
   *  Executes the function with the specified {link parameters}, and returns result.
   *  If necessary, the corresponding package will be loaded as part of the call.
   * */
  async apply(parameters: object = {}): Promise<any> {
    return (await this.prepare(parameters).call()).getOutputParamValue();
  }

  /** Returns functions with the specified attributes. */
  static find(params?: { package?: string, name?: string, tags?: string[], returnType?: string}): Func[] {
    return api.grok_Func_Find(params?.package, params?.name, params?.tags, params?.returnType);
  }

  /** Returns functions (including these) with the specified attributes. */
  static async findAll(params?: { package?: string, name?: string, tags?: string[], returnType?: string}): Promise<Func[]> {
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

/** Represents a project */
export class Project extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get pictureUrl(): string { return api.grok_PictureMixin_Get_PictureUrl(this.dart); }
  get path(): string { return api.grok_Project_Get_Path(this.dart); }
  get isOnServer(): string { return api.grok_Project_Get_IsOnServer(this.dart); }
  get isLocal(): string { return api.grok_Project_Get_IsLocal(this.dart); }

  /** Project description */
  get description(): string { return api.grok_Project_Description(this.dart); }

  /** Project changes flag */
  get isDirty(): boolean { return api.grok_Project_IsDirty(this.dart); }

  /** Project is empty flag */
  get isEmpty(): boolean { return api.grok_Project_IsEmpty(this.dart); }

  toMarkup(): string { return api.grok_Project_ToMarkup(this.dart); }

  /** Opens the project in workspace */
  open(options?: {closeAll: boolean}): Promise<Project> {
    return api.grok_Project_Open(this.dart, options?.closeAll ?? false);
  }
}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-query}
 * */
export class DataQuery extends Func {
  /** @constructs DataQuery*/
  constructor(dart: any) {
    super(dart);
  }

  get adHoc(): boolean { return api.grok_Query_Get_AdHoc(this.dart); }
  set adHoc(a: boolean) { api.grok_Query_Set_AdHoc(this.dart, a); }

  /** Query text */
  get query(): string { return api.grok_Query_Query(this.dart); }
  set query(q: string) { api.grok_Query_Set_Query(this.dart, q); }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-job}
 * */
export class DataJob extends Func {
  /** @constructs DataJob */
  constructor(dart: any) {
    super(dart);
  }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/data-connection}
 * */
export class DataConnection extends Entity {
  /** @constructs DataConnection */
  constructor(dart: any) {
    super(dart);
  }

  /** Collection of parameters: server, database, endpoint, etc. */
  get parameters(): object { return api.grok_DataConnection_Parameters(this.dart); }

  /** Tests the connection, returns "ok" on success or an error message on error
   * @returns {Promise<string>}*/
  test(): Promise<string> {
    return new Promise((resolve, reject) => api.grok_DataConnection_test(this.dart, (s: any) => resolve(toJs(s)), (e: any) => reject(e)));
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
    return new Promise((resolve, reject) => api.grok_Notebook_Edit(this.dart, (f: any) => resolve(f), (e: any) => reject(e)));
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
    return new Promise((resolve, reject) => api.grok_Notebook_ToHtml(this.dart, (html: any) => resolve(html), (e: any) => reject(e)));
  }
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
  /** @constructs TableInfo */
  constructor(dart: any) {
    super(dart);
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

  /** @returns {Promise<string>} */
  // readAsString(): Promise<string> {
  //   return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsString(this.dart, (x: any) => resolve(x), (x: any) => reject(x)));
  // }

  readAsString(): Promise<string> {
    return api.grok_FileInfo_ReadAsString(this.dart);
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes(): Promise<Uint8Array> {
    return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsBytes(this.dart, (x: any) => resolve(x), (x: any) => reject(x)));
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

}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
  /** @constructs Script */
  constructor(dart: any) {
    super(dart);
  }

  static create(script: string): Script { return new Script(api.grok_Script_Create(script)); }

  /** Script */
  get script(): string { return api.grok_Script_GetScript(this.dart); }
  set script(s: string) { api.grok_Script_SetScript(this.dart, s); }

  /** Script language. See also: https://datagrok.ai/help/compute/scripting#header-parameters */
  get language(): ScriptLanguage { return api.grok_Script_GetLanguage(this.dart); }
  set language(s: ScriptLanguage) { api.grok_Script_SetLanguage(this.dart, s); }

  /** Environment name. See also: https://datagrok.ai/help/compute/scripting#header-parameters */
  get environment(): string { return api.grok_Script_Get_Environment(this.dart); }
  set environment(s: string) { api.grok_Script_Set_Environment(this.dart, s); }

  /** Reference header parameter. See also: https://datagrok.ai/help/compute/scripting#header-parameters */
  get reference(): string { return api.grok_Script_Get_Reference(this.dart); }
  set reference(s: string) { api.grok_Script_Set_Reference(this.dart, s); }

  /** Sample table. See also: https://datagrok.ai/help/compute/scripting#header-parameters */
  get sample(): string { return api.grok_Script_Get_Sample(this.dart); }
  set sample(s: string) { api.grok_Script_Set_Sample(this.dart, s); }

  /** Script tags. See also: https://datagrok.ai/help/compute/scripting#header-parameters */
  get tags(): string[] { return api.grok_Script_Get_Tags(this.dart); }
  set tags(tags: string[]) { api.grok_Script_Set_Tags(this.dart, tags); }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
 *  */
export class Credentials extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Collection of parameters: login, password, API key, etc. */
  get parameters(): object { return api.grok_Credentials_Parameters(this.dart); }
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
    return new Promise((resolve, reject) => api.grok_ScriptEnvironment_Setup(this.dart, () => resolve(), (e: any) => reject(e)));
  }
}

export class LogEventType extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Friendly name of the event type */
  get name(): string { return api.grok_LogEventType_Get_Name(this.dart); }
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
  public webRoot: string | undefined;
  public version: string;

  constructor(dart: any | undefined = undefined) {
    super(dart);
    if (typeof dart === 'string') {
      this.webRoot = dart;
      this.dart = null;
    }
    this.version = "";
  }

  /** Override init() method to provide package-specific initialization.
   * It is guaranteed to get called exactly once before the execution of any function below.
   */
  init(): Promise<null> { return Promise.resolve(null); }

  private _name: string | undefined;
  /** Package short name
   *  @type {string} */
  get name() {
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

  get meta() { return (this.dart == null) ? null : toJs(api.grok_Package_Get_Meta(this.dart)); }

  /** Loads package
   * @returns {Promise<Package>} */
  async load(): Promise<Package> {
    return new Promise((resolve, reject) =>
      api.grok_Dapi_Packages_Load(this.name, (data: any) => resolve(data), (e: any) => reject(e)));
  }

  /** Returns credentials for package
   * @returns {Promise<Credentials>} */
  getCredentials(): Promise<Credentials> {
    return new Promise((resolve, reject) => api.grok_Package_Get_Credentials(this.name, (c: any) => {
      let cred = toJs(c);
      resolve(cred);
    }, (e: any) => reject(e)));
  }

  /** Returns properties for a package. */
  getProperties(): Promise<Map<string, any>> {
    return new Promise((resolve, reject) => api.grok_Package_Get_Properties(this.name, (props: any) => {
      resolve(toJs(props));
    }, (e: any) => reject(e)));
  }

  private _files: FileSource | null = null;
  get files(): FileSource {
    if (this._files == null)
      this._files = new FileSource(`System:AppData/${this.name}`);
    return this._files;
  }
}


export interface PropertyOptions {

  /** Property name */
  name?: string;

  /** Property type */
  type?: string;

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

  /** List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names.
   * Signature: validator(x: DG.Type): string | null.
   * [null] indicates that the value is valid, [string] describes a validation error. */
  validators?: string[];

  /** Custom field caption shown in [PropertyGrid] */
  caption?: string;

  /** Field postfix shown in [PropertyGrid]. [units] take precedence over the [postfix] value. */
  postfix?: string;
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

  /** Property category */
  get category(): string { return api.grok_Property_Get_Category(this.dart); }
  set category(s: string) { api.grok_Property_Set_Category(this.dart, s); }

  /** Property type */
  get propertyType(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set propertyType(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Semantic type */
  get semType(): SemType | string { return api.grok_Property_Get_SemType(this.dart); }
  set semType(s: SemType | string) { api.grok_Property_Set_SemType(this.dart, s); }

  /** Description */
  get description(): string { return api.grok_Property_Get_Description(this.dart); }
  set description(s: string) { api.grok_Property_Set_Description(this.dart, s); }

  /** Default value */
  get defaultValue(): any { return api.grok_Property_Get_DefaultValue(this.dart); }
  set defaultValue(s: any) { api.grok_Property_Set_DefaultValue(this.dart, s); }

  /** Property editor */
  get editor(): string { return api.grok_Property_Get(this.dart, 'editor'); }
  set editor(s: string) { api.grok_Property_Set(this.dart, 'editor', s); }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns {Array<string>} */
  get choices(): string[] { return api.grok_Property_Get_Choices(this.dart); }
  set choices(x: string[]) { api.grok_Property_Set_Choices(this.dart, x); }

  /** Column type filter */
  get columnFilter(): ColumnType | 'numerical' | 'categorical' | null {
    return api.grok_Property_Get_ColumnTypeFilter(this.dart);
  }

  /** List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names.
   * Signature: validator(x: DG.Type): string | null.
   * [null] indicates that the value is valid, [string] describes a validation error. */
  get validators(): string[] { return api.grok_Property_Get_Validators(this.dart); }
  set validators(x: string[]) { api.grok_Property_Set_Validators(this.dart, x); }

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
      (x: any, v: any) => x[name] = v,
      options?.defaultValue).fromOptions(options);
  }

  static jsInt(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.INT, options); }
  static jsBool(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.BOOL, options); }
  static jsFloat(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.FLOAT, options); }
  static jsString(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.STRING, options); }
  static jsDateTime(name: string, options?: PropertyOptions): Property { return Property.js(name, TYPE.DATE_TIME, options); }

  static fromOptions(options: PropertyOptions): Property { return Property.js(options.name!, options.type! as TYPE, options); }
}


export class DateTime {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static fromDate(date: Date): DateTime {
    return DateTime.fromMillisecondsSinceEpoch(date.getTime());
  }

  static fromMillisecondsSinceEpoch(millisecondsSinceEpoch: number): DateTime {
    return new DateTime(api.grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
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