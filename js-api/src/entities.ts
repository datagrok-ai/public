import {ColumnType, SemType, Type, TYPE} from "./const";
import { FuncCall } from "./functions";
import {toJs} from "./wrappers";

declare var grok: any;
let api = <any>window;

type PropertyGetter = (a: string) => any;
type PropertySetter = (a: string, object:  any) => void;

/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {

  public d: any;

  /** @constructs Entity*/
  constructor(d: any) {
    this.d = d;
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id(): string { return api.grok_Entity_Get_Id(this.d); }
  set id(x: string) { api.grok_Entity_Set_Id(this.d, x); }

  /** Entity friendly name
   *  @type {string} */
  get friendlyName(): string { return api.grok_Entity_Get_FriendlyName(this.d); }
  set friendlyName(x: string) { api.grok_Entity_Set_FriendlyName(this.d, x); }

  /** Entity short name */
  get name(): string { return api.grok_Entity_Get_Name(this.d); }
  set name(x: string) { api.grok_Entity_Set_Name(this.d, x); }

  /** Entity full-qualified name */
  get nqName(): string { return api.grok_Entity_Get_nqName(this.d); }

  /** Entity path */
  get path(): string { return api.grok_Entity_Path(this.d); }

  /** Entity properties */
  getProperties(): Promise<Map<string, any>> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.s, this.d, (p: any) => resolve(p), (e: any) => reject(e)));
  }

  /** Sets entity properties */
  setProperties(props: Map<string, any>): Promise<any> {
    return new Promise((resolve, reject) => api.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.s, this.d, props, (_: any) => resolve(_), (e: any) => reject(e)));
  }

  /** Returns a string representing the object */
  toString(): string { return api.grok_Object_ToString(this.d); }

}

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {
  /** @constructs User*/
  constructor(d: any) {
    super(d);
  }

  static fromId(id: string): User { return new User(api.grok_User_From_Id(id)); }

  static current(): User { return new User(api.grok_User()); }

  /** First name */
  get firstName(): string { return api.grok_User_Get_FirstName(this.d); }

  /** Last name */
  get lastName(): string { return api.grok_User_Get_LastName(this.d); }

  /** Email */
  get email(): string | null { return api.grok_User_Get_Email(this.d); }

  /** Picture URL */
  get picture(): string | object { return api.grok_User_Get_Picture(this.d); }

  /** User home project */
  get project(): Project { return toJs(api.grok_User_Get_Project(this.d)); }

  /** User home folder connection */
  get home(): DataConnection { return toJs(api.grok_User_Get_Storage(this.d)); }

  /** Login */
  get login(): string { return api.grok_User_Get_Login(this.d); }

  toMarkup(): string { return api.grok_User_ToMarkup(this.d); }

  /** Security Group */
  get group(): Group { return toJs(api.grok_User_Get_Group(this.d)); }
}


/**
 * Represents a user session in the Datagrok platform.
 * @extends Entity
 * */
export class UserSession extends Entity {
  /** @constructs UserSession*/
  constructor(d: any) {
    super(d);
  }

  /** Entity ID (GUID) */
  get id(): string { return api.grok_Entity_Get_Id(this.d); }
  set id(x: string) { api.grok_Entity_Set_Id(this.d, x); }

  get type(): string { return api.grok_UserSession_Get_Type(this.d); }

  /** External Token */
  get externalToken(): string { return api.grok_UserSession_Get_ExternalToken(this.d); }

  /** User */
  get user(): User { return toJs(api.grok_UserSession_Get_User(this.d)); }
}

/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/overview/functions/function}
 * */
export class Func extends Entity {
  /** @constructs Func*/
  constructor(d: any) {
    super(d);
  }

  /** Returns {@link FuncCall} object in a stand-by state
   * @param {object} parameters
   * @returns {FuncCall} */
  prepare(parameters: object = {}): FuncCall {
    return toJs(api.grok_Func_Prepare(this.d, parameters));
  };
}

/** Represents a project */
export class Project extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Project description */
  get description(): string { return api.grok_Project_Description(this.d); }

  /** Project changes flag */
  get isDirty(): boolean { return api.grok_Project_IsDirty(this.d); }

  /** Project is empty flag */
  get isEmpty(): boolean { return api.grok_Project_IsEmpty(this.d); }

  toMarkup(): string { return api.grok_Project_ToMarkup(this.d); }
}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-query}
 * */
export class DataQuery extends Func {
  /** @constructs DataQuery*/
  constructor(d: any) {
    super(d);
  }

  /** Query text */
  get query(): string { return api.grok_Query_Query(this.d); }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-job}
 * */
export class DataJob extends Func {
  /** @constructs DataJob */
  constructor(d: any) {
    super(d);
  }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/data-connection}
 * */
export class DataConnection extends Entity {
  /** @constructs DataConnection */
  constructor(d: any) {
    super(d);
  }

  /** Collection of parameters: server, database, endpoint, etc. */
  get parameters(): object { return api.grok_DataConnection_Parameters(this.d); }
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
 * */
export class Model extends Entity {
  /** @constructs Model */
  constructor(d: any) {
    super(d);
  }
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/develop/jupyter-notebook}
 * */
export class Notebook extends Entity {
  /** @constructs Notebook */
  constructor(d: any) {
    super(d);
  }

  /** Create Notebook on server for edit.
   * @returns {Promise<string>} Current notebook's name */
  edit(): Promise<string> {
    return new Promise((resolve, reject) => api.grok_Notebook_Edit(this.d, (f: any) => resolve(f), (e: any) => reject(e)));
  }

  /** Environment name */
  get environment(): string { return api.grok_Notebook_Get_Environment(this.d); }
  set environment(e: string) { api.grok_Notebook_Set_Environment(this.d, e); }

  /** Description */
  get description(): string { return api.grok_Notebook_Get_Description(this.d); }
  set description(e: string) { api.grok_Notebook_Set_Description(this.d, e); }

  /** Converts Notebook to HTML code
   * @returns {Promise<string>} */
  toHtml(): Promise<string> {
    return new Promise((resolve, reject) => api.grok_Notebook_ToHtml(this.d, (html: any) => resolve(html), (e: any) => reject(e)));
  }
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
  /** @constructs TableInfo */
  constructor(d: any) {
    super(d);
  }
}

/** @extends Entity
 * Allows for files handling in JS-based info panels
 * {@link https://datagrok.ai/help/discover/info-panels}
 * */
export class FileInfo extends Entity {
  /** @constructs FileInfo */
  constructor(d: any) {
    super(d);
  }

  /** Returns path, i.e. `geo/dmv_offices.csv` */
  get path(): string { return api.grok_FileInfo_Get_Path(this.d); }

  /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv` */
  get fullPath(): string { return api.grok_FileInfo_Get_FullPath(this.d); }

  /** Returns file extension, i.e. `csv` */
  get extension(): string { return api.grok_FileInfo_Get_Extension(this.d); }

  /** Returns file name, i.e. `dmv_offices.csv` */
  get fileName(): string { return api.grok_FileInfo_Get_FileName(this.d); }

  /** Returns file URL */
  get url(): string { return api.grok_FileInfo_Get_Url(this.d); }

  /** @returns {Promise<string>} */
  readAsString(): Promise<string> {
    return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsString(this.d, (x: any) => resolve(x), (x: any) => reject(x)));
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes(): Promise<Uint8Array> {
    return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsBytes(this.d, (x: any) => resolve(x), (x: any) => reject(x)));
  }
}

/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
  /** @constructs Group */
  constructor(d: any) {
    super(d);
  }

  static create(name: string): Group { return new Group(api.grok_Group(name)); }

  /** Adds a member to the group
   * @param {Group} m */
  addMember(m: Group): void { api.grok_Group_Add_Member(this.d, m.d, false); }

  /** Adds an admin member to the group
   * @param {Group} m */
  addAdminMember(m: Group): void { api.grok_Group_Add_Member(this.d, m.d, true); }

  /** Removes a member from the group
   * @param {Group} m */
  removeMember(m: Group): void { api.grok_Group_Remove_Member(this.d, m.d); }

  /** Adds the group to another one
   * @param {Group} m */
  includeTo(m: Group): void { api.grok_Group_Add_Membership(this.d, m.d, false); }

  /** Adds the group to another one as an admin
   * @param {Group} m */
  includeAdminTo(m: Group): void { api.grok_Group_Add_Membership(this.d, m.d, true); }

  /** Removes membership from another group
   * @param {Group} m */
  excludeFrom(m: Group): void { api.grok_Group_Remove_Membership(this.d, m.d); }

  /** Returns list of groups that belong to group, with no admin permissions
   * @type {Array<Group>} */
  get members(): Group[] { return toJs(api.grok_Group_Get_Members(this.d, false)); }

  /** Returns list of groups that belong to group, with admin permissions
   * @type {Array<Group>} */
  get adminMembers(): Group[] { return toJs(api.grok_Group_Get_Members(this.d, true)); }

  /** Returns list of groups that group belongs to, with no admin permissions
   * @type {Array<Group>} */
  get memberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.d, false)); }

  /** Returns list of groups that group belongs to, with admin permissions
   * @type {list<Group>} */
  get adminMemberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.d, true)); }

  /** Personal user group */
  get personal(): boolean { return api.grok_Group_Get_Personal(this.d); }
  set personal(e: boolean) { api.grok_Group_Set_Personal(this.d, e); }

  /** Hidden group */
  get hidden(): boolean { return api.grok_Group_Get_Hidden(this.d); }
  set hidden(e: boolean) { api.grok_Group_Set_Hidden(this.d, e); }

}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
  /** @constructs Script */
  constructor(d: any) {
    super(d);
  }

  static create(script: string): Script { return new Script(api.grok_Script_Create(script)); }

  /** Script */
  get script(): string { return api.grok_Script_GetScript(this.d); }
  set script(s: string) { api.grok_Script_SetScript(this.d, s); }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
 *  */
export class Credentials extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Collection of parameters: login, password, API key, etc. */
  get parameters(): object { return api.grok_Credentials_Parameters(this.d); }
}

/** Represents a script environment */
export class ScriptEnvironment extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Create instance of ScriptEnvironment
   * @param {string} name
   * @returns {ScriptEnvironment}
   * */
  static create(name: string): ScriptEnvironment {
    return new ScriptEnvironment(api.grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content */
  get environment(): string { return api.grok_ScriptEnvironment_Environment(this.d); }

  /** Setup environment */
  setup(): Promise<void> {
    return new Promise((resolve, reject) => api.grok_ScriptEnvironment_Setup(this.d, () => resolve(), (e: any) => reject(e)));
  }
}

export class LogEventType extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Friendly name of the event type */
  get name(): string { return api.grok_LogEventType_Get_Name(this.d); }
}

export class LogEvent extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Description of the event */
  get description(): string { return api.grok_LogEvent_Get_Description(this.d); }

  /** Friendly name of the event */
  get name(): string { return api.grok_LogEvent_Get_Name(this.d); }

  /** Source of the event */
  get source(): string { return api.grok_LogEvent_Get_Source(this.d); }

  /** Session id of the event */
  get session(): UserSession | string { return toJs(api.grok_LogEvent_Get_Session(this.d)); }

  /** Parameters of the event
   * @type {Array<LogEventParameterValue>} */
  get parameters(): LogEventParameterValue[] { return api.grok_LogEvent_Get_Parameters(this.d); }

  /** Type of the event
   * @type {LogEventType} */
  get eventType(): LogEventType { return toJs(api.grok_LogEvent_Get_Type(this.d)); }
}

export class LogEventParameter extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Name of the parameter */
  get name(): string { return api.grok_LogEventParameter_Get_Name(this.d); }

  /** Type of the parameter */
  get type(): string { return api.grok_LogEventParameter_Get_Type(this.d); }

  /** Description of the parameter */
  get description(): string { return api.grok_LogEventParameter_Get_Description(this.d); }

  /** Is the parameter input */
  get isInput(): boolean { return api.grok_LogEventParameter_Get_IsInput(this.d); }

  /** Is the parameter optional */
  get isOptional(): boolean { return api.grok_LogEventParameter_Get_IsOptional(this.d); }
}

export class LogEventParameterValue extends Entity {
  constructor(d: any) {
    super(d);
  }

  /** Event of the parameter value
   * @type {LogEvent} */
  get event(): LogEvent {
    return toJs(api.grok_LogEventParameterValue_Get_Event(this.d));
  }

  /** Parameter of the parameter value
   * @type {LogEventParameter} */
  get parameter(): LogEventParameter {
    return toJs(api.grok_LogEventParameterValue_Get_Parameter(this.d));
  }

  /** Parameter value
   * @type {string} */
  get value(): string {
    return api.grok_LogEventParameterValue_Get_Value(this.d);
  }
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package extends Entity {
  public webRoot: string | String | undefined;
  public version: string;

  constructor(d: any | undefined = undefined) {
    super(d);
    if (typeof d === 'string' || d instanceof String) {
      this.webRoot = d;
      this.d = null;
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
    if (this.d != null)
      return api.grok_Entity_Get_Name(this.d);
    else
      return this._name;
  }

  set name(x) {
    if (this.d != null)
      api.grok_Entity_Set_Name(this.d, x);
    else
      this._name = x;
  }

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

  /** Returns properties for package
   * @returns {Promise<object>} */
  getProperties(): Promise<Map<string, any>> {
    return new Promise((resolve, reject) => api.grok_Package_Get_Properties(this.name, (c: any) => {
      resolve(c);
    }, (e: any) => reject(e)));
  }
}

/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property {
  public readonly d: any;

  constructor(d: any) {
    this.d = d;
  }

  /** Property getter is a function that accepts one parameter (item)
   * and returns the property value.
   *
   * @returns {PropertyGetter} */
  get get(): PropertyGetter { return api.grok_Property_Get_Get(this.d); }
  set get(x: PropertyGetter) { api.grok_Property_Set_Get(this.d, x); }

  /** Property setter */
  get set(): PropertySetter { return api.grok_Property_Get_Set(this.d); }
  set set(x: PropertySetter) { api.grok_Property_Set_Set(this.d, x); }

  /** Property name */
  get name(): string { return api.grok_Property_Get_Name(this.d); }
  set name(s: string) { api.grok_Property_Set_Name(this.d, s); }

  /** Property type */
  get propertyType(): string { return api.grok_Property_Get_PropertyType(this.d); }
  set propertyType(s: string) { api.grok_Property_Set_PropertyType(this.d, s); }

  /** Semantic type */
  get semType(): SemType | string { return api.grok_Property_Get_SemType(this.d); }
  set semType(s: SemType | string) { api.grok_Property_Set_SemType(this.d, s); }

  /** Description */
  get description(): string { return api.grok_Property_Get_Description(this.d); }
  set description(s: string) { api.grok_Property_Set_Description(this.d, s); }

  /** Default value */
  get defaultValue(): any { return api.grok_Property_Get_DefaultValue(this.d); }
  set defaultValue(s: any) { api.grok_Property_Set_DefaultValue(this.d, s); }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns {Array<string>} */
  get choices(): string[] { return api.grok_Property_Get_Choices(this.d); }
  set choices(x: string[]) { api.grok_Property_Set_Choices(this.d, x); }

  /** Column type filter */
  get columnFilter(): ColumnType | 'numerical' | 'categorical' | null {
    return api.grok_Property_Get_ColumnTypeFilter(this.d);
  }

  /** Creates a property
   * @param {string} name
   * @param {TYPE} type
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static create(name: string, type: Type, getter: PropertyGetter, setter: PropertySetter, defaultValue: any = null): Property {
    return new Property(api.grok_Property(name, type, getter, setter, defaultValue));
  }

  /** Creates an integer property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static int(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.INT, getter, setter, defaultValue);
  }

  /** Creates a float property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static float(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.FLOAT, getter, setter, defaultValue);
  }

  /** Creates a string property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static string(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.STRING, getter, setter, defaultValue);
  }
}


export class DateTime {
  public d: any;

  constructor(d: any) {
    this.d = d;
  }

  static fromDate(date: Date): DateTime {
    return DateTime.fromMillisecondsSinceEpoch(date.getTime());
  }

  static fromMillisecondsSinceEpoch(millisecondsSinceEpoch: number): DateTime {
    return new DateTime(api.grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
  }
}
