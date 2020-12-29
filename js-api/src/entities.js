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
    return grok_Entity_Get_Id(this.d);
  }

  set id(x) {
    return grok_Entity_Set_Id(this.d, x);
  }

  /** Entity name
   *  @type {string} */
  get name() {
    return grok_Entity_Get_Name(this.d);
  }

  set name(x) {
    return grok_Entity_Set_Name(this.d, x);
  }

  /** Entity path
   *  @type {string} */
  get path() {
    return grok_Entity_Path(this.d);
  }

  /** Returns entity properties
   * @returns {Promise<Map>} props */
  getProperties() {
    return new Promise((resolve, reject) => grok_EntitiesDataSource_GetProperties(grok.dapi.entities.s, this.d, (p) => resolve(p), (e) => reject(e)));
  }

  /** Allows to set properties for entity
   * @param {Map} props
   * @returns Promise */
  setProperties(props) {
    return new Promise((resolve, reject) => grok_EntitiesDataSource_SetProperties(grok.dapi.entities.s, this.d, props, (_) => resolve(_), (e) => reject(e)));
  }

  /** Returns a string representing the object
   * @returns {string} */
  toString() {
    return grok_Object_ToString(this.d);
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
    return new User(grok_User_From_Id(id));
  }

  static current() {
    return new User(grok_User());
  }

  /** First name
   * @type {string} */
  get firstName() {
    return grok_User_Get_FirstName(this.d);
  }

  /** Last name
   * @type {string} */
  get lastName() {
    return grok_User_Get_LastName(this.d);
  }

  /** Email
   * @type {string} */
  get email() {
    return grok_User_Get_Email(this.d);
  }

  /** Picture URL
   * @type {string} */
  get picture() {
    return grok_User_Get_Picture(this.d);
  }

  /** Login
   *  @type {string} */
  get login() {
    return grok_User_Get_Login(this.d);
  }

  /** */
  toMarkup() {
    return grok_User_ToMarkup(this.d);
  }

  /** Security Group
   * @type Group
   */
  get group() {
    return toJs(grok_User_Get_Group(this.d))
  }
}


/**
 * Represents a user session in the Datagrok platform.
 * @extends Entity
 * */
export class UserSession {
  /** @constructs UserSession*/
  constructor(d) {
    this.d = d;
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id() {
    return grok_Entity_Get_Id(this.d);
  }

  set id(x) {
    return grok_Entity_Set_Id(this.d, x);
  }

  /** Login
   *  @type {string} */
  get type() {
    return grok_UserSession_Get_Type(this.d);
  }

  /** External Token
   *  @type {string} */
  get externalToken() {
    return grok_UserSession_Get_ExternalToken(this.d);
  }

  /** User
   *  @type {User} */
  get user() {
    return toJs(grok_UserSession_Get_User(this.d));
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
    return toJs(grok_Func_Prepare(this.d, parameters));
  };
}

/** Represents a project */
export class Project extends Entity {
  constructor(d) {
    super(d);
  }

  /** Entity name
   *  @type {string} */
  get description() {
    return grok_Project_Description(this.d);
  }

  /** Entity name
   *  @type {string} */
  get isDirty() {
    return grok_Project_IsDirty(this.d);
  }

  /** Entity name
   *  @type {string} */
  get isEmpty() {
    return grok_Project_IsEmpty(this.d);
  }

  toMarkup() {
    return grok_Project_ToMarkup(this.d);
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
    return grok_Query_Query(this.d);
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
    return grok_DataConnection_Parameters(this.d);
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
    return new Promise((resolve, reject) => grok_Notebook_Edit(this.d, (f) => resolve(f), (e) => reject(e)));
  }

  /** Environment name
   * @type {string} */
  get environment() {
    return grok_Notebook_Get_Environment(this.d);
  }

  set environment(e) {
    return grok_Notebook_Set_Environment(this.d, e);
  }

  /** Description
   * @type {string} */
  get description() {
    return grok_Notebook_Get_Description(this.d);
  }

  set description(e) {
    return grok_Notebook_Set_Description(this.d, e);
  }

  /** Converts Notebook to HTML code
   * @returns {Promise<string>} */
  toHtml() {
    return new Promise((resolve, reject) => grok_Notebook_ToHtml(this.d, (html) => resolve(html), (e) => reject(e)));
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
    return grok_FileInfo_Get_Path(this.d);
  }

  /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv`
   * @type {string} */
  get fullPath() {
    return grok_FileInfo_Get_FullPath(this.d);
  }

  /** Returns file extension, i.e. `csv`
   * @type {string} */
  get extension() {
    return grok_FileInfo_Get_Extension(this.d);
  }

  /** Returns file name, i.e. `dmv_offices.csv`
   * @type {string} */
  get fileName() {
    return grok_FileInfo_Get_FileName(this.d);
  }

  /** @returns {string} */
  get url() { grok_FileInfo_Get_Url(this.d); }

  /** @returns {Promise<string>} */
  readAsString() {
    return new Promise((resolve, reject) => grok_FileInfo_ReadAsString(this.d, (x) => resolve(x), (x) => reject(x)));
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes() {
    return new Promise((resolve, reject) => grok_FileInfo_ReadAsBytes(this.d, (x) => resolve(x), (x) => reject(x)));
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
    return new Group(grok_Group(name));
  }

  /** Adds a member to the group
   * @param {Group} m */
  addMember(m) {
    grok_Group_Add_Member(this.d, m.d, false);
  }

  /** Adds an admin member to the group
   * @param {Group} m */
  addAdminMember(m) {
    grok_Group_Add_Member(this.d, m.d, true);
  }

  /** Removes a member from the group
   * @param {Group} m */
  removeMember(m) {
    grok_Group_Remove_Member(this.d, m.d);
  }

  /** Adds the group to another one
   * @param {Group} m */
  includeTo(m) {
    grok_Group_Add_Membership(this.d, m.d, false);
  }

  /** Adds the group to another one as an admin
   * @param {Group} m */
  includeAdminTo(m) {
    grok_Group_Add_Membership(this.d, m.d, true);
  }

  /** Removes membership from another group
   * @param {Group} m */
  excludeFrom(m) {
    grok_Group_Remove_Membership(this.d, m.d);
  }

  /** Returns list of groups that belong to group, with no admin permissions
   * @type {List<Group>} */
  get members() {
    return toJs(grok_Group_Get_Members(this.d, false));
  }

  /** Returns list of groups that belong to group, with admin permissions
   * @type {List<Group>} */
  get adminMembers() {
    return toJs(grok_Group_Get_Members(this.d, true));
  }

  /** Returns list of groups that group belongs to, with no admin permissions
   * @type {List<Group>} */
  get memberships() {
    return toJs(grok_Group_Get_Memberships(this.d, false));
  }

  /** Returns list of groups that group belongs to, with admin permissions
   * @type {List<Group>} */
  get adminMemberships() {
    return toJs(grok_Group_Get_Memberships(this.d, true));
  }

  /** Personal user group
   * @type {boolean} */
  get personal() {
    return grok_Group_Get_Personal(this.d);
  }

  set personal(e) {
    return grok_Group_Set_Personal(this.d, e);
  }

  /** Hidden group
   * @type {boolean} */
  get hidden() {
    return grok_Group_Get_Hidden(this.d);
  }

  set hidden(e) {
    return grok_Group_Set_Hidden(this.d, e);
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
    return new Script(grok_Script_Create(script));
  }

  /** Script
   * @type {string} */
  get script() {
    return grok_Script_GetScript(this.d);
  }

  set script(s) {
    return grok_Script_SetScript(this.d, s);
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
    return grok_Credentials_Parameters(this.d);
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
    return new ScriptEnvironment(grok_ScriptEnvironment_Create(name));
  }

  /** Environment yaml file content
   * @type {string} */
  get environment() {
    return grok_ScriptEnvironment_Environment(this.d);
  }

  /** Setup environment */
  setup() {
    return new Promise((resolve, reject) => grok_ScriptEnvironment_Setup(this.d, () => resolve(), (e) => reject(e)));
  }
}

export class LogEventType extends Entity {
  constructor(d) {
    super(d);
  }

  /** Friendly name of the event type
   * @type {string} */
  get name() {
    return grok_LogEventType_Get_Name(this.d);
  }
}

export class LogEvent extends Entity {
  constructor(d) {
    super(d);
  }

  /** Description of the event
   * @type {string} */
  get description() {
    return grok_LogEvent_Get_Description(this.d);
  }

  /** Friendly name of the event
   * @type {string} */
  get name() {
    return grok_LogEvent_Get_Name(this.d);
  }

  /** Source of the event
   * @type {string} */
  get source() {
    return grok_LogEvent_Get_Source(this.d);
  }

  /** Session id of the event
   * @type {UserSession} */
  get session() {
    return toJs(grok_LogEvent_Get_Session(this.d));
  }

  /** Parameters of the event
   * @type {ArrayList<LogEventParameterValue>} */
  get parameters() {
    return grok_LogEvent_Get_Parameters(this.d);
  }

  /** Type of the event
   * @type {LogEventType} */
  get eventType() {
    return toJs(grok_LogEvent_Get_Type(this.d));
  }
}

export class LogEventParameter extends Entity {
  constructor(d) {
    super(d);
  }

  /** Name of the parameter
   * @type {string} */
  get name() {
    return grok_LogEventParameter_Get_Name(this.d);
  }

  /** Type of the parameter
   * @type {string} */
  get type() {
    return grok_LogEventParameter_Get_Type(this.d);
  }

  /** Description of the parameter
   * @type {string} */
  get description() {
    return grok_LogEventParameter_Get_Description(this.d);
  }

  /** Is the parameter input
   * @type {Boolean} */
  get isInput() {
    return grok_LogEventParameter_Get_IsInput(this.d);
  }

  /** Is the parameter optional
   * @type {Boolean} */
  get isOptional() {
    return grok_LogEventParameter_Get_IsOptional(this.d);
  }
}

export class LogEventParameterValue extends Entity {
  constructor(d) {
    super(d);
  }

  /** Event of the parameter value
   * @type {LogEvent} */
  get event() {
    return toJs(grok_LogEventParameterValue_Get_Event(this.d));
  }

  /** Parameter of the parameter value
   * @type {LogEventParameter} */
  get parameter() {
    return toJs(grok_LogEventParameterValue_Get_Parameter(this.d));
  }

  /** Parameter value
   * @type {string} */
  get value() {
    return grok_LogEventParameterValue_Get_Value(this.d);
  }
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package {
  constructor(webRoot) {
    this.webRoot = webRoot;
    this.name = "";
    this.version = "";
  }

  /** Override init() method to provide package-specific initialization.
   * It is guaranteed to get called exactly once before the execution of any function below.
   * */

  /*async*/
  init() {
    return Promise.resolve(null);
  }

  /** Returns credentials for package
   * @returns {Promise<Credentials>} */
  getCredentials() {
    return new Promise((resolve, reject) => grok_Package_Get_Credentials(this.name, (c) => {
      let cred = toJs(c);
      resolve(cred);
    }, (e) => reject(e)));
  }

  /** Returns properties for package
   * @returns {Promise<object>} */
  getProperties() {
    return new Promise((resolve, reject) => grok_Package_Get_Properties(this.name, (c) => {
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

  /** Property getter is a function that acccepts one parameter (item)
   * and returns the property value. */
  get get() {
    return grok_Property_Get_Get(this.d);
  }

  set get(x) {
    return grok_Property_Set_Get(this.d, x);
  }

  /** Property setter */
  get set() {
    return grok_Property_Get_Set(this.d);
  }

  set set(x) {
    return grok_Property_Set_Set(this.d, x);
  }

  /** Property name
   *  @type {string} */
  get name() {
    return grok_Property_Get_Name(this.d);
  }

  set name(s) {
    grok_Property_Set_Name(this.d, s);
  }

  /** Property type
   *  @type {string} */
  get propertyType() {
    return grok_Property_Get_PropertyType(this.d);
  }

  set propertyType(s) {
    grok_Property_Set_PropertyType(this.d, s);
  }

  /** Semantic type
   *  @type {string} */
  get semType() {
    return grok_Property_Get_SemType(this.d);
  }

  set semType(s) {
    grok_Property_Set_SemType(this.d, s);
  }

  /** Description */
  get description() {
    return grok_Property_Get_Description(this.d);
  }

  set description(s) {
    grok_Property_Set_Description(this.d, s);
  }

  /** Default value */
  get defaultValue() {
    return grok_Property_Get_DefaultValue(this.d);
  }

  set defaultValue(s) {
    grok_Property_Set_DefaultValue(this.d, s);
  }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns ArrayList<string>*/
  get choices() {
    return grok_Property_Get_Choices(this.d);
  }

  set choices(x) {
    return grok_Property_Set_Choices(this.d, x);
  }

  /** Column type filter
   * @type {string}*/
  get columnFilter() {
    return grok_Property_Get_ColumnTypeFilter(this.d);
  }

  /** Creates a property
   * @param {string} name
   * @param {TYPE|Type} type
   * @param {function} getter
   * @param {function} setter
   * @param {object} defaultValue
   * @returns Property*/
  static create(name, type, getter, setter, defaultValue = null) {
    return new Property(grok_Property(name, type, getter, setter, defaultValue));
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
    return new SemanticValue(grok_SemanticValue(value, semType));
  }

  get value() {
    return grok_SemanticValue_Get_Value(this.d);
  }

  set value(x) {
    grok_SemanticValue_Set_Value(this.d, x);
  }

  get semType() {
    return grok_SemanticValue_Get_SemType(this.d);
  }

  set semType(x) {
    grok_SemanticValue_Set_SemType(this.d, x);
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
    return new DateTime(grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
  }
}
