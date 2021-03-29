import {TYPE} from "./const";
import {toJs} from "./wrappers";

/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {
  /** @constructs Entity*/
  constructor(d) {
    this.d = d;
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id() {
    return window.grok_Entity_Get_Id(this.d);
  }

  set id(x) {
    window.grok_Entity_Set_Id(this.d, x);
  }

  /** Entity friendly name
   *  @type {string} */
  get friendlyName() {
    return window.grok_Entity_Get_FriendlyName(this.d);
  }

  set friendlyName(x) {
    window.grok_Entity_Set_FriendlyName(this.d, x);
  }

  /** Entity short name
   *  @type {string} */
  get name() {
    return window.grok_Entity_Get_Name(this.d);
  }

  set name(x) {
    window.grok_Entity_Set_Name(this.d, x);
  }

  /** Entity full-qualified name
   *  @type {string} */
  get nqName() {
    return window.grok_Entity_Get_nqName(this.d);
  }

  /** Entity path
   *  @type {string} */
  get path() {
    return window.grok_Entity_Path(this.d);
  }

  /** Returns entity properties
   * @returns {Promise<Map>} props */
  getProperties() {
    return new Promise((resolve, reject) => window.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.s, this.d, (p) => resolve(p), (e) => reject(e)));
  }

  /** Allows to set properties for entity
   * @param {Map} props
   * @returns Promise */
  setProperties(props) {
    return new Promise((resolve, reject) => window.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.s, this.d, props, (_) => resolve(_), (e) => reject(e)));
  }

  /** Returns a string representing the object
   * @returns {string} */
  toString() {
    return window.grok_Object_ToString(this.d);
  }
}

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {
  /** @constructs User*/
  constructor(d) {
    super(d);
  }

  static fromId(id) {
    return new User(window.grok_User_From_Id(id));
  }

  static current() {
    return new User(window.grok_User());
  }

  /** First name
   * @type {string} */
  get firstName() {
    return window.grok_User_Get_FirstName(this.d);
  }

  /** Last name
   * @type {string} */
  get lastName() {
    return window.grok_User_Get_LastName(this.d);
  }

  /** Email
   * @type {string} */
  get email() {
    return window.grok_User_Get_Email(this.d);
  }

  /** Picture URL
   * @type {string} */
  get picture() {
    return window.grok_User_Get_Picture(this.d);
  }

  /** User home project
   * @type {Project} */
  get project() {
    return toJs(window.grok_User_Get_Project(this.d));
  }

  /** User home folder connection
   * @type {Project} */
  get home() {
    return toJs(window.grok_User_Get_Storage(this.d));
  }

  /** Login
   *  @type {string} */
  get login() {
    return window.grok_User_Get_Login(this.d);
  }

  /** */
  toMarkup() {
    return window.grok_User_ToMarkup(this.d);
  }

  /** Security Group
   * @type Group
   */
  get group() {
    return toJs(window.grok_User_Get_Group(this.d));
  }
}

/**
 * Represents a user session in the Datagrok platform.
 * @extends Entity
 * */
export class UserSession extends Entity {
  /** @constructs UserSession*/
  constructor(d) {
    super(d);
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id() {
    return window.grok_Entity_Get_Id(this.d);
  }

  set id(x) {
    window.grok_Entity_Set_Id(this.d, x);
  }

  /** Login
   *  @type {string} */
  get type() {
    return window.grok_UserSession_Get_Type(this.d);
  }

  /** External Token
   *  @type {string} */
  get externalToken() {
    return window.grok_UserSession_Get_ExternalToken(this.d);
  }

  /** User
   *  @type {User} */
  get user() {
    return toJs(window.grok_UserSession_Get_User(this.d));
  }
}

/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/overview/functions/function}
 * */
export class Func extends Entity {
  /** @constructs Func*/
  constructor(d) {
    super(d);
  }

  /** Returns {@link FuncCall} object in a stand-by state
   * @param {object} parameters
   * @returns {FuncCall} */
  prepare(parameters = {}) {
    return toJs(window.grok_Func_Prepare(this.d, parameters));
  }
  ;
}

/** Represents a project */
export class Project extends Entity {
  constructor(d) {
    super(d);
  }

  /** Project description
   *  @type {string} */
  get description() {
    return window.grok_Project_Description(this.d);
  }

  /** Project changes flag
   *  @type {string} */
  get isDirty() {
    return window.grok_Project_IsDirty(this.d);
  }

  /** Project is empty flag
   *  @type {string} */
  get isEmpty() {
    return window.grok_Project_IsEmpty(this.d);
  }

  toMarkup() {
    return window.grok_Project_ToMarkup(this.d);
  }
}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-query}
 * */
export class DataQuery extends Func {
  /** @constructs DataQuery*/
  constructor(d) {
    super(d);
  }

  /** Query text
   *  @type {string} */
  get query() {
    return window.grok_Query_Query(this.d);
  }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-job}
 * */
export class DataJob extends Func {
  /** @constructs DataJob */
  constructor(d) {
    super(d);
  }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/data-connection}
 * */
export class DataConnection extends Entity {
  /** @constructs DataConnection */
  constructor(d) {
    super(d);
  }

  /** Collection of parameters: server, database, endpoint, etc.
   *  @type {object} */
  get parameters() {
    return window.grok_DataConnection_Parameters(this.d);
  }
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
 * */
export class Model extends Entity {
  /** @constructs Model */
  constructor(d) {
    super(d);
  }
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/develop/jupyter-notebook}
 * */
export class Notebook extends Entity {
  /** @constructs Notebook */
  constructor(d) {
    super(d);
  }

  /** Create Notebook on server for edit.
   * @returns {Promise<string>} Current notebook's name */
  edit() {
    return new Promise((resolve, reject) => window.grok_Notebook_Edit(this.d, (f) => resolve(f), (e) => reject(e)));
  }

  /** Environment name
   * @type {string} */
  get environment() {
    return window.grok_Notebook_Get_Environment(this.d);
  }

  set environment(e) {
    window.grok_Notebook_Set_Environment(this.d, e);
  }

  /** Description
   * @type {string} */
  get description() {
    return window.grok_Notebook_Get_Description(this.d);
  }

  set description(e) {
    window.grok_Notebook_Set_Description(this.d, e);
  }

  /** Converts Notebook to HTML code
   * @returns {Promise<string>} */
  toHtml() {
    return new Promise((resolve, reject) => window.grok_Notebook_ToHtml(this.d, (html) => resolve(html), (e) => reject(e)));
  }
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
  /** @constructs TableInfo */
  constructor(d) {
    super(d);
  }
}

/** @extends Entity
 * Allows for files handling in JS-based info panels
 * {@link https://datagrok.ai/help/discover/info-panels}
 * */
export class FileInfo extends Entity {
  /** @constructs FileInfo */
  constructor(d) {
    super(d);
  }

  /** Returns path, i.e. `geo/dmv_offices.csv`
   * @type {string} */
  get path() {
    return window.grok_FileInfo_Get_Path(this.d);
  }

  /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv`
   * @type {string} */
  get fullPath() {
    return window.grok_FileInfo_Get_FullPath(this.d);
  }

  /** Returns file extension, i.e. `csv`
   * @type {string} */
  get extension() {
    return window.grok_FileInfo_Get_Extension(this.d);
  }

  /** Returns file name, i.e. `dmv_offices.csv`
   * @type {string} */
  get fileName() {
    return window.grok_FileInfo_Get_FileName(this.d);
  }

  /** @returns {string} */
  get url() {
    return window.grok_FileInfo_Get_Url(this.d);
  }

  /** @returns {Promise<string>} */
  readAsString() {
    return new Promise((resolve, reject) => window.grok_FileInfo_ReadAsString(this.d, (x) => resolve(x), (x) => reject(x)));
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes() {
    return new Promise((resolve, reject) => window.grok_FileInfo_ReadAsBytes(this.d, (x) => resolve(x), (x) => reject(x)));
  }
}

/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
  /** @constructs Group */
  constructor(d) {
    super(d);
  }

  static create(name) {
    return new Group(window.grok_Group(name));
  }

  /** Adds a member to the group
   * @param {Group} m */
  addMember(m) {
    window.grok_Group_Add_Member(this.d, m.d, false);
  }

  /** Adds an admin member to the group
   * @param {Group} m */
  addAdminMember(m) {
    window.grok_Group_Add_Member(this.d, m.d, true);
  }

  /** Removes a member from the group
   * @param {Group} m */
  removeMember(m) {
    window.grok_Group_Remove_Member(this.d, m.d);
  }

  /** Adds the group to another one
   * @param {Group} m */
  includeTo(m) {
    window.grok_Group_Add_Membership(this.d, m.d, false);
  }

  /** Adds the group to another one as an admin
   * @param {Group} m */
  includeAdminTo(m) {
    window.grok_Group_Add_Membership(this.d, m.d, true);
  }

  /** Removes membership from another group
   * @param {Group} m */
  excludeFrom(m) {
    window.grok_Group_Remove_Membership(this.d, m.d);
  }

  /** Returns list of groups that belong to group, with no admin permissions
   * @type {Array<Group>} */
  get members() {
    return toJs(window.grok_Group_Get_Members(this.d, false));
  }

  /** Returns list of groups that belong to group, with admin permissions
   * @type {Array<Group>} */
  get adminMembers() {
    return toJs(window.grok_Group_Get_Members(this.d, true));
  }

  /** Returns list of groups that group belongs to, with no admin permissions
   * @type {Array<Group>} */
  get memberships() {
    return toJs(window.grok_Group_Get_Memberships(this.d, false));
  }

  /** Returns list of groups that group belongs to, with admin permissions
   * @type {list<Group>} */
  get adminMemberships() {
    return toJs(window.grok_Group_Get_Memberships(this.d, true));
  }

  /** Personal user group
   * @type {boolean} */
  get personal() {
    return window.grok_Group_Get_Personal(this.d);
  }

  set personal(e) {
    window.grok_Group_Set_Personal(this.d, e);
  }

  /** Hidden group
   * @type {boolean} */
  get hidden() {
    return window.grok_Group_Get_Hidden(this.d);
  }

  set hidden(e) {
    window.grok_Group_Set_Hidden(this.d, e);
  }
}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
  /** @constructs Script */
  constructor(d) {
    super(d);
  }

  static create(script) {
    return new Script(window.grok_Script_Create(script));
  }

  /** Script
   * @type {string} */
  get script() {
    return window.grok_Script_GetScript(this.d);
  }

  set script(s) {
    window.grok_Script_SetScript(this.d, s);
  }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
 *  */
export class Credentials extends Entity {
  constructor(d) {
    super(d);
  }

  /** Collection of parameters: login, password, API key, etc.
   *  @type {object} */
  get parameters() {
    return window.grok_Credentials_Parameters(this.d);
  }
}

/** Represents a script environment */
export class ScriptEnvironment extends Entity {
  constructor(d) {
    super(d);
  }

  /** Create instance of ScriptEnvironment
   * @param {string} name
   * @returns {ScriptEnvironment}
   * */
  static create(name) {
    return new ScriptEnvironment(window.grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content
   * @type {string} */
  get environment() {
    return window.grok_ScriptEnvironment_Environment(this.d);
  }

  /** Setup environment */
  setup() {
    return new Promise((resolve, reject) => window.grok_ScriptEnvironment_Setup(this.d, () => resolve(), (e) => reject(e)));
  }
}

export class LogEventType extends Entity {
  constructor(d) {
    super(d);
  }

  /** Friendly name of the event type
   * @type {string} */
  get name() {
    return window.grok_LogEventType_Get_Name(this.d);
  }
}

export class LogEvent extends Entity {
  constructor(d) {
    super(d);
  }

  /** Description of the event
   * @type {string} */
  get description() {
    return window.grok_LogEvent_Get_Description(this.d);
  }

  /** Friendly name of the event
   * @type {string} */
  get name() {
    return window.grok_LogEvent_Get_Name(this.d);
  }

  /** Source of the event
   * @type {string} */
  get source() {
    return window.grok_LogEvent_Get_Source(this.d);
  }

  /** Session id of the event
   * @type {UserSession} */
  get session() {
    return toJs(window.grok_LogEvent_Get_Session(this.d));
  }

  /** Parameters of the event
   * @type {Array<LogEventParameterValue>} */
  get parameters() {
    return window.grok_LogEvent_Get_Parameters(this.d);
  }

  /** Type of the event
   * @type {LogEventType} */
  get eventType() {
    return toJs(window.grok_LogEvent_Get_Type(this.d));
  }
}

export class LogEventParameter extends Entity {
  constructor(d) {
    super(d);
  }

  /** Name of the parameter
   * @type {string} */
  get name() {
    return window.grok_LogEventParameter_Get_Name(this.d);
  }

  /** Type of the parameter
   * @type {string} */
  get type() {
    return window.grok_LogEventParameter_Get_Type(this.d);
  }

  /** Description of the parameter
   * @type {string} */
  get description() {
    return window.grok_LogEventParameter_Get_Description(this.d);
  }

  /** Is the parameter input
   * @type {Boolean} */
  get isInput() {
    return window.grok_LogEventParameter_Get_IsInput(this.d);
  }

  /** Is the parameter optional
   * @type {Boolean} */
  get isOptional() {
    return window.grok_LogEventParameter_Get_IsOptional(this.d);
  }
}

export class LogEventParameterValue extends Entity {
  constructor(d) {
    super(d);
  }

  /** Event of the parameter value
   * @type {LogEvent} */
  get event() {
    return toJs(window.grok_LogEventParameterValue_Get_Event(this.d));
  }

  /** Parameter of the parameter value
   * @type {LogEventParameter} */
  get parameter() {
    return toJs(window.grok_LogEventParameterValue_Get_Parameter(this.d));
  }

  /** Parameter value
   * @type {string} */
  get value() {
    return window.grok_LogEventParameterValue_Get_Value(this.d);
  }
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package extends Entity {
  constructor(d) {
    super(d);
    if (typeof d === 'string' || d instanceof String) {
      this.webRoot = d;
      this.d = null;
    }
    this.version = "";
  }

  /** Override init() method to provide package-specific initialization.
   * It is guaranteed to get called exactly once before the execution of any function below.
   * */

  /*async*/
  init() {
    return Promise.resolve(null);
  }

  /** Package short name
   *  @type {string} */
  get name() {
    if (this.d != null)
      return window.grok_Entity_Get_Name(this.d);
    else
      return this._name;
  }

  set name(x) {
    if (this.d != null)
      window.grok_Entity_Set_Name(this.d, x);
    else
      this._name = x;
  }

  /** load package
   * @returns {Promise<Package>} */
  async load() {
    return new Promise((resolve, reject) => window.grok_Dapi_Packages_Load(this.name, (data) => resolve(data), (e) => reject(e)));
  }

  /** Returns credentials for package
   * @returns {Promise<Credentials>} */
  getCredentials() {
    return new Promise((resolve, reject) => window.grok_Package_Get_Credentials(this.name, (c) => {
      let cred = toJs(c);
      resolve(cred);
    }, (e) => reject(e)));
  }

  /** Returns properties for package
   * @returns {Promise<object>} */
  getProperties() {
    return new Promise((resolve, reject) => window.grok_Package_Get_Properties(this.name, (c) => {
      resolve(c);
    }, (e) => reject(e)));
  }
}

/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property {
  constructor(d) {
    this.d = d;
  }

  /** Property getter is a function that accepts one parameter (item)
   * and returns the property value.
   *
   * @returns {PropertyGetter} */
  get get() {
    return window.grok_Property_Get_Get(this.d);
  }

  /** @param {PropertyGetter} x */
  set get(x) {
    window.grok_Property_Set_Get(this.d, x);
  }

  /** @returns {PropertySetter} */
  get set() {
    return window.grok_Property_Get_Set(this.d);
  }

  /** @param {PropertySetter} x*/
  set set(x) {
    window.grok_Property_Set_Set(this.d, x);
  }

  /** Property name
   *  @type {string} */
  get name() {
    return window.grok_Property_Get_Name(this.d);
  }

  set name(s) {
    window.grok_Property_Set_Name(this.d, s);
  }

  /** Property type
   *  @type {string} */
  get propertyType() {
    return window.grok_Property_Get_PropertyType(this.d);
  }

  set propertyType(s) {
    window.grok_Property_Set_PropertyType(this.d, s);
  }

  /** Semantic type
   *  @type {string} */
  get semType() {
    return window.grok_Property_Get_SemType(this.d);
  }

  set semType(s) {
    window.grok_Property_Set_SemType(this.d, s);
  }

  /** Description */
  get description() {
    return window.grok_Property_Get_Description(this.d);
  }

  set description(s) {
    window.grok_Property_Set_Description(this.d, s);
  }

  /** Default value */
  get defaultValue() {
    return window.grok_Property_Get_DefaultValue(this.d);
  }

  set defaultValue(s) {
    window.grok_Property_Set_DefaultValue(this.d, s);
  }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns ArrayList<string>*/
  get choices() {
    return window.grok_Property_Get_Choices(this.d);
  }

  set choices(x) {
    window.grok_Property_Set_Choices(this.d, x);
  }

  /** Column type filter
   * @type {string}*/
  get columnFilter() {
    return window.grok_Property_Get_ColumnTypeFilter(this.d);
  }

  /** Creates a property
   * @param {string} name
   * @param {TYPE} type
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static create(name, type, getter, setter, defaultValue = null) {
    return new Property(window.grok_Property(name, type, getter, setter, defaultValue));
  }

  /** Creates an integer property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static int(name, getter, setter, defaultValue) {
    return Property.create(name, TYPE.INT, getter, setter, defaultValue);
  }

  /** Creates a float property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static float(name, getter, setter, defaultValue) {
    return Property.create(name, TYPE.FLOAT, getter, setter, defaultValue);
  }

  /** Creates a string property
   * @param {string} name
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static string(name, getter, setter, defaultValue) {
    return Property.create(name, TYPE.STRING, getter, setter, defaultValue);
  }
}

export class SemanticValue {
  constructor(d) {
    this.d = d;
  }

  static fromValueType(value, semType) {
    return new SemanticValue(window.grok_SemanticValue(value, semType));
  }

  get value() {
    return window.grok_SemanticValue_Get_Value(this.d);
  }

  set value(x) {
    window.grok_SemanticValue_Set_Value(this.d, x);
  }

  get semType() {
    return window.grok_SemanticValue_Get_SemType(this.d);
  }

  set semType(x) {
    window.grok_SemanticValue_Set_SemType(this.d, x);
  }
}

export class DateTime {
  constructor(d) {
    this.d = d;
  }

  static fromDate(date) {
    return DateTime.fromMillisecondsSinceEpoch(date.getTime());
  }

  static fromMillisecondsSinceEpoch(millisecondsSinceEpoch) {
    return new DateTime(window.grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
  }
}
