var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
  function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
  return new (P || (P = Promise))(function (resolve, reject) {
    function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
    function rejected(value) { try { step(generator['throw'](value)); } catch (e) { reject(e); } }
    function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
    step((generator = generator.apply(thisArg, _arguments || [])).next());
  });
};
define(['require', 'exports', './const', './wrappers'], function (require, exports, const_1, wrappers_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.DateTime = exports.SemanticValue = exports.Property = exports.Package = exports.LogEventParameterValue = exports.LogEventParameter = exports.LogEvent = exports.LogEventType = exports.ScriptEnvironment = exports.Credentials = exports.Script = exports.Group = exports.FileInfo = exports.TableInfo = exports.Notebook = exports.Model = exports.DataConnection = exports.DataJob = exports.DataQuery = exports.Project = exports.Func = exports.UserSession = exports.User = exports.Entity = void 0;
  let api = window;
  /** @class
     * Base class for system objects stored in the database in a structured manner.
     * Contains base properties: id, name and path
     * */
  class Entity {
    /** @constructs Entity*/
    constructor(d) {
      this.d = d;
    }
    /** Entity ID (GUID)
         *  @type {string} */
    get id() {
      return api.grok_Entity_Get_Id(this.d);
    }
    set id(x) {
      api.grok_Entity_Set_Id(this.d, x);
    }
    /** Entity friendly name
         *  @type {string} */
    get friendlyName() {
      return api.grok_Entity_Get_FriendlyName(this.d);
    }
    set friendlyName(x) {
      api.grok_Entity_Set_FriendlyName(this.d, x);
    }
    /** Entity short name
         *  @type {string} */
    get name() {
      return api.grok_Entity_Get_Name(this.d);
    }
    set name(x) {
      api.grok_Entity_Set_Name(this.d, x);
    }
    /** Entity full-qualified name
         *  @type {string} */
    get nqName() {
      return api.grok_Entity_Get_nqName(this.d);
    }
    /** Entity path
         *  @type {string} */
    get path() {
      return api.grok_Entity_Path(this.d);
    }
    /** Returns entity properties
         * @returns {Promise<Map>} props */
    getProperties() {
      return new Promise((resolve, reject) => api.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.s, this.d, (p) => resolve(p), (e) => reject(e)));
    }
    /** Allows to set properties for entity
         * @param {Map} props
         * @returns Promise */
    setProperties(props) {
      return new Promise((resolve, reject) => api.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.s, this.d, props, (_) => resolve(_), (e) => reject(e)));
    }
    /** Returns a string representing the object
         * @returns {string} */
    toString() {
      return api.grok_Object_ToString(this.d);
    }
  }
  exports.Entity = Entity;
  /**
     * Represents a user of the Datagrok platform.
     * @extends Entity
     * */
  class User extends Entity {
    /** @constructs User*/
    constructor(d) {
      super(d);
    }
    static fromId(id) {
      return new User(api.grok_User_From_Id(id));
    }
    static current() {
      return new User(api.grok_User());
    }
    /** First name
         * @type {string} */
    get firstName() {
      return api.grok_User_Get_FirstName(this.d);
    }
    /** Last name
         * @type {string} */
    get lastName() {
      return api.grok_User_Get_LastName(this.d);
    }
    /** Email
         * @type {string} */
    get email() {
      return api.grok_User_Get_Email(this.d);
    }
    /** Picture URL
         * @type {string} */
    get picture() {
      return api.grok_User_Get_Picture(this.d);
    }
    /** User home project
         * @type {Project} */
    get project() {
      return wrappers_1.toJs(api.grok_User_Get_Project(this.d));
    }
    /** User home folder connection
         * @type {Project} */
    get home() {
      return wrappers_1.toJs(api.grok_User_Get_Storage(this.d));
    }
    /** Login
         *  @type {string} */
    get login() {
      return api.grok_User_Get_Login(this.d);
    }
    /** */
    toMarkup() {
      return api.grok_User_ToMarkup(this.d);
    }
    /** Security Group
         * @type Group
         */
    get group() {
      return wrappers_1.toJs(api.grok_User_Get_Group(this.d));
    }
  }
  exports.User = User;
  /**
     * Represents a user session in the Datagrok platform.
     * @extends Entity
     * */
  class UserSession extends Entity {
    /** @constructs UserSession*/
    constructor(d) {
      super(d);
    }
    /** Entity ID (GUID)
         *  @type {string} */
    get id() {
      return api.grok_Entity_Get_Id(this.d);
    }
    set id(x) {
      api.grok_Entity_Set_Id(this.d, x);
    }
    /** Login
         *  @type {string} */
    get type() {
      return api.grok_UserSession_Get_Type(this.d);
    }
    /** External Token
         *  @type {string} */
    get externalToken() {
      return api.grok_UserSession_Get_ExternalToken(this.d);
    }
    /** User
         *  @type {User} */
    get user() {
      return wrappers_1.toJs(api.grok_UserSession_Get_User(this.d));
    }
  }
  exports.UserSession = UserSession;
  /** Represents a function
     * @extends Entity
     * {@link https://datagrok.ai/help/overview/functions/function}
     * */
  class Func extends Entity {
    /** @constructs Func*/
    constructor(d) {
      super(d);
    }
    /** Returns {@link FuncCall} object in a stand-by state
         * @param {object} parameters
         * @returns {FuncCall} */
    prepare(parameters = {}) {
      return wrappers_1.toJs(api.grok_Func_Prepare(this.d, parameters));
    }
    ;
  }
  exports.Func = Func;
  /** Represents a project */
  class Project extends Entity {
    constructor(d) {
      super(d);
    }
    /** Project description
         *  @type {string} */
    get description() {
      return api.grok_Project_Description(this.d);
    }
    /** Project changes flag
         *  @type {string} */
    get isDirty() {
      return api.grok_Project_IsDirty(this.d);
    }
    /** Project is empty flag
         *  @type {string} */
    get isEmpty() {
      return api.grok_Project_IsEmpty(this.d);
    }
    toMarkup() {
      return api.grok_Project_ToMarkup(this.d);
    }
  }
  exports.Project = Project;
  /** Represents a data query
     * @extends Func
     * {@link https://datagrok.ai/help/access/data-query}
     * */
  class DataQuery extends Func {
    /** @constructs DataQuery*/
    constructor(d) {
      super(d);
    }
    /** Query text
         *  @type {string} */
    get query() {
      return api.grok_Query_Query(this.d);
    }
  }
  exports.DataQuery = DataQuery;
  /** Represents a data job
     * @extends Func
     * {@link https://datagrok.ai/help/access/data-job}
     * */
  class DataJob extends Func {
    /** @constructs DataJob */
    constructor(d) {
      super(d);
    }
  }
  exports.DataJob = DataJob;
  /** Represents a data connection
     * @extends Entity
     * {@link https://datagrok.ai/help/access/data-connection}
     * */
  class DataConnection extends Entity {
    /** @constructs DataConnection */
    constructor(d) {
      super(d);
    }
    /** Collection of parameters: server, database, endpoint, etc.
         *  @type {object} */
    get parameters() {
      return api.grok_DataConnection_Parameters(this.d);
    }
  }
  exports.DataConnection = DataConnection;
  /** Represents a predictive model
     * @extends Entity
     * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
     * */
  class Model extends Entity {
    /** @constructs Model */
    constructor(d) {
      super(d);
    }
  }
  exports.Model = Model;
  /** @extends Entity
     * Represents a Jupyter notebook
     * {@link https://datagrok.ai/help/develop/jupyter-notebook}
     * */
  class Notebook extends Entity {
    /** @constructs Notebook */
    constructor(d) {
      super(d);
    }
    /** Create Notebook on server for edit.
         * @returns {Promise<string>} Current notebook's name */
    edit() {
      return new Promise((resolve, reject) => api.grok_Notebook_Edit(this.d, (f) => resolve(f), (e) => reject(e)));
    }
    /** Environment name
         * @type {string} */
    get environment() {
      return api.grok_Notebook_Get_Environment(this.d);
    }
    set environment(e) {
      api.grok_Notebook_Set_Environment(this.d, e);
    }
    /** Description
         * @type {string} */
    get description() {
      return api.grok_Notebook_Get_Description(this.d);
    }
    set description(e) {
      api.grok_Notebook_Set_Description(this.d, e);
    }
    /** Converts Notebook to HTML code
         * @returns {Promise<string>} */
    toHtml() {
      return new Promise((resolve, reject) => api.grok_Notebook_ToHtml(this.d, (html) => resolve(html), (e) => reject(e)));
    }
  }
  exports.Notebook = Notebook;
  /** @extends Entity
     * Represents a Table metadata
     * */
  class TableInfo extends Entity {
    /** @constructs TableInfo */
    constructor(d) {
      super(d);
    }
  }
  exports.TableInfo = TableInfo;
  /** @extends Entity
     * Allows for files handling in JS-based info panels
     * {@link https://datagrok.ai/help/discover/info-panels}
     * */
  class FileInfo extends Entity {
    /** @constructs FileInfo */
    constructor(d) {
      super(d);
    }
    /** Returns path, i.e. `geo/dmv_offices.csv`
         * @type {string} */
    get path() {
      return api.grok_FileInfo_Get_Path(this.d);
    }
    /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv`
         * @type {string} */
    get fullPath() {
      return api.grok_FileInfo_Get_FullPath(this.d);
    }
    /** Returns file extension, i.e. `csv`
         * @type {string} */
    get extension() {
      return api.grok_FileInfo_Get_Extension(this.d);
    }
    /** Returns file name, i.e. `dmv_offices.csv`
         * @type {string} */
    get fileName() {
      return api.grok_FileInfo_Get_FileName(this.d);
    }
    /** @returns {string} */
    get url() {
      return api.grok_FileInfo_Get_Url(this.d);
    }
    /** @returns {Promise<string>} */
    readAsString() {
      return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsString(this.d, (x) => resolve(x), (x) => reject(x)));
    }
    /** @returns {Promise<Uint8Array>} */
    readAsBytes() {
      return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsBytes(this.d, (x) => resolve(x), (x) => reject(x)));
    }
  }
  exports.FileInfo = FileInfo;
  /** @extends Entity
     * Represents a User Group
     * */
  class Group extends Entity {
    /** @constructs Group */
    constructor(d) {
      super(d);
    }
    static create(name) {
      return new Group(api.grok_Group(name));
    }
    /** Adds a member to the group
         * @param {Group} m */
    addMember(m) {
      api.grok_Group_Add_Member(this.d, m.d, false);
    }
    /** Adds an admin member to the group
         * @param {Group} m */
    addAdminMember(m) {
      api.grok_Group_Add_Member(this.d, m.d, true);
    }
    /** Removes a member from the group
         * @param {Group} m */
    removeMember(m) {
      api.grok_Group_Remove_Member(this.d, m.d);
    }
    /** Adds the group to another one
         * @param {Group} m */
    includeTo(m) {
      api.grok_Group_Add_Membership(this.d, m.d, false);
    }
    /** Adds the group to another one as an admin
         * @param {Group} m */
    includeAdminTo(m) {
      api.grok_Group_Add_Membership(this.d, m.d, true);
    }
    /** Removes membership from another group
         * @param {Group} m */
    excludeFrom(m) {
      api.grok_Group_Remove_Membership(this.d, m.d);
    }
    /** Returns list of groups that belong to group, with no admin permissions
         * @type {Array<Group>} */
    get members() {
      return wrappers_1.toJs(api.grok_Group_Get_Members(this.d, false));
    }
    /** Returns list of groups that belong to group, with admin permissions
         * @type {Array<Group>} */
    get adminMembers() {
      return wrappers_1.toJs(api.grok_Group_Get_Members(this.d, true));
    }
    /** Returns list of groups that group belongs to, with no admin permissions
         * @type {Array<Group>} */
    get memberships() {
      return wrappers_1.toJs(api.grok_Group_Get_Memberships(this.d, false));
    }
    /** Returns list of groups that group belongs to, with admin permissions
         * @type {list<Group>} */
    get adminMemberships() {
      return wrappers_1.toJs(api.grok_Group_Get_Memberships(this.d, true));
    }
    /** Personal user group
         * @type {boolean} */
    get personal() {
      return api.grok_Group_Get_Personal(this.d);
    }
    set personal(e) {
      api.grok_Group_Set_Personal(this.d, e);
    }
    /** Hidden group
         * @type {boolean} */
    get hidden() {
      return api.grok_Group_Get_Hidden(this.d);
    }
    set hidden(e) {
      api.grok_Group_Set_Hidden(this.d, e);
    }
  }
  exports.Group = Group;
  /** @extends Func
     * Represents a Script
     * */
  class Script extends Func {
    /** @constructs Script */
    constructor(d) {
      super(d);
    }
    static create(script) {
      return new Script(api.grok_Script_Create(script));
    }
    /** Script
         * @type {string} */
    get script() {
      return api.grok_Script_GetScript(this.d);
    }
    set script(s) {
      api.grok_Script_SetScript(this.d, s);
    }
  }
  exports.Script = Script;
  /** Represents connection credentials
     *  Usually it is a login and a password pair
     *  Passwords are stored in the secured credentials storage
     *  See also: {@link https://datagrok.ai/help/govern/security}
     *  */
  class Credentials extends Entity {
    constructor(d) {
      super(d);
    }
    /** Collection of parameters: login, password, API key, etc.
         *  @type {object} */
    get parameters() {
      return api.grok_Credentials_Parameters(this.d);
    }
  }
  exports.Credentials = Credentials;
  /** Represents a script environment */
  class ScriptEnvironment extends Entity {
    constructor(d) {
      super(d);
    }
    /** Create instance of ScriptEnvironment
         * @param {string} name
         * @returns {ScriptEnvironment}
         * */
    static create(name) {
      return new ScriptEnvironment(api.grok_ScriptEnvironment_Create(name));
    }
    /** Environment yaml file content
         * @type {string} */
    get environment() {
      return api.grok_ScriptEnvironment_Environment(this.d);
    }
    /** Setup environment */
    setup() {
      return new Promise((resolve, reject) => api.grok_ScriptEnvironment_Setup(this.d, () => resolve(), (e) => reject(e)));
    }
  }
  exports.ScriptEnvironment = ScriptEnvironment;
  class LogEventType extends Entity {
    constructor(d) {
      super(d);
    }
    /** Friendly name of the event type
         * @type {string} */
    get name() {
      return api.grok_LogEventType_Get_Name(this.d);
    }
  }
  exports.LogEventType = LogEventType;
  class LogEvent extends Entity {
    constructor(d) {
      super(d);
    }
    /** Description of the event
         * @type {string} */
    get description() {
      return api.grok_LogEvent_Get_Description(this.d);
    }
    /** Friendly name of the event
         * @type {string} */
    get name() {
      return api.grok_LogEvent_Get_Name(this.d);
    }
    /** Source of the event
         * @type {string} */
    get source() {
      return api.grok_LogEvent_Get_Source(this.d);
    }
    /** Session id of the event
         * @type {UserSession} */
    get session() {
      return wrappers_1.toJs(api.grok_LogEvent_Get_Session(this.d));
    }
    /** Parameters of the event
         * @type {Array<LogEventParameterValue>} */
    get parameters() {
      return api.grok_LogEvent_Get_Parameters(this.d);
    }
    /** Type of the event
         * @type {LogEventType} */
    get eventType() {
      return wrappers_1.toJs(api.grok_LogEvent_Get_Type(this.d));
    }
  }
  exports.LogEvent = LogEvent;
  class LogEventParameter extends Entity {
    constructor(d) {
      super(d);
    }
    /** Name of the parameter
         * @type {string} */
    get name() {
      return api.grok_LogEventParameter_Get_Name(this.d);
    }
    /** Type of the parameter
         * @type {string} */
    get type() {
      return api.grok_LogEventParameter_Get_Type(this.d);
    }
    /** Description of the parameter
         * @type {string} */
    get description() {
      return api.grok_LogEventParameter_Get_Description(this.d);
    }
    /** Is the parameter input
         * @type {Boolean} */
    get isInput() {
      return api.grok_LogEventParameter_Get_IsInput(this.d);
    }
    /** Is the parameter optional
         * @type {Boolean} */
    get isOptional() {
      return api.grok_LogEventParameter_Get_IsOptional(this.d);
    }
  }
  exports.LogEventParameter = LogEventParameter;
  class LogEventParameterValue extends Entity {
    constructor(d) {
      super(d);
    }
    /** Event of the parameter value
         * @type {LogEvent} */
    get event() {
      return wrappers_1.toJs(api.grok_LogEventParameterValue_Get_Event(this.d));
    }
    /** Parameter of the parameter value
         * @type {LogEventParameter} */
    get parameter() {
      return wrappers_1.toJs(api.grok_LogEventParameterValue_Get_Parameter(this.d));
    }
    /** Parameter value
         * @type {string} */
    get value() {
      return api.grok_LogEventParameterValue_Get_Value(this.d);
    }
  }
  exports.LogEventParameterValue = LogEventParameterValue;
  /**
     * Represents a package, which is a unit of distribution of content in the Datagrok platform.
     */
  class Package extends Entity {
    constructor(d = undefined) {
      super(d);
      if (typeof d === 'string' || d instanceof String) {
        this.webRoot = d;
        this.d = null;
      }
      this.version = '';
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
    /** load package
         * @returns {Promise<Package>} */
    load() {
      return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => api.grok_Dapi_Packages_Load(this.name, (data) => resolve(data), (e) => reject(e)));
      });
    }
    /** Returns credentials for package
         * @returns {Promise<Credentials>} */
    getCredentials() {
      return new Promise((resolve, reject) => api.grok_Package_Get_Credentials(this.name, (c) => {
        let cred = wrappers_1.toJs(c);
        resolve(cred);
      }, (e) => reject(e)));
    }
    /** Returns properties for package
         * @returns {Promise<object>} */
    getProperties() {
      return new Promise((resolve, reject) => api.grok_Package_Get_Properties(this.name, (c) => {
        resolve(c);
      }, (e) => reject(e)));
    }
  }
  exports.Package = Package;
  /**
     * Strongly-typed property associated with an object.
     * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
     *
     * Samples:
     */
  class Property {
    constructor(d) {
      this.d = d;
    }
    /** Property getter is a function that accepts one parameter (item)
         * and returns the property value.
         *
         * @returns {PropertyGetter} */
    get get() {
      return api.grok_Property_Get_Get(this.d);
    }
    /** @param {PropertyGetter} x */
    set get(x) {
      api.grok_Property_Set_Get(this.d, x);
    }
    get set() {
      return api.grok_Property_Get_Set(this.d);
    }
    set set(x) {
      api.grok_Property_Set_Set(this.d, x);
    }
    /** Property name
         *  @type {string} */
    get name() {
      return api.grok_Property_Get_Name(this.d);
    }
    set name(s) {
      api.grok_Property_Set_Name(this.d, s);
    }
    /** Property type
         *  @type {string} */
    get propertyType() {
      return api.grok_Property_Get_PropertyType(this.d);
    }
    set propertyType(s) {
      api.grok_Property_Set_PropertyType(this.d, s);
    }
    /** Semantic type
         *  @type {string} */
    get semType() {
      return api.grok_Property_Get_SemType(this.d);
    }
    set semType(s) {
      api.grok_Property_Set_SemType(this.d, s);
    }
    /** Description */
    get description() {
      return api.grok_Property_Get_Description(this.d);
    }
    set description(s) {
      api.grok_Property_Set_Description(this.d, s);
    }
    /** Default value */
    get defaultValue() {
      return api.grok_Property_Get_DefaultValue(this.d);
    }
    set defaultValue(s) {
      api.grok_Property_Set_DefaultValue(this.d, s);
    }
    /** List of possible values of that property.
         *  PropertyGrid will use it to populate combo boxes.
         *  @returns ArrayList<string>*/
    get choices() {
      return api.grok_Property_Get_Choices(this.d);
    }
    set choices(x) {
      api.grok_Property_Set_Choices(this.d, x);
    }
    /** Column type filter
         * @type {string}*/
    get columnFilter() {
      return api.grok_Property_Get_ColumnTypeFilter(this.d);
    }
    /** Creates a property
         * @param {string} name
         * @param {TYPE} type
         * @param {function} getter
         * @param {function} setter
         * @param {object} defaultValue
         * @returns Property*/
    static create(name, type, getter, setter, defaultValue = null) {
      return new Property(api.grok_Property(name, type, getter, setter, defaultValue));
    }
    /** Creates an integer property
         * @param {string} name
         * @param {function} getter
         * @param {function} setter
         * @param {object} defaultValue
         * @returns Property*/
    static int(name, getter, setter, defaultValue) {
      return Property.create(name, const_1.TYPE.INT, getter, setter, defaultValue);
    }
    /** Creates a float property
         * @param {string} name
         * @param {function} getter
         * @param {function} setter
         * @param {object} defaultValue
         * @returns Property*/
    static float(name, getter, setter, defaultValue) {
      return Property.create(name, const_1.TYPE.FLOAT, getter, setter, defaultValue);
    }
    /** Creates a string property
         * @param {string} name
         * @param {function} getter
         * @param {function} setter
         * @param {object} defaultValue
         * @returns Property*/
    static string(name, getter, setter, defaultValue) {
      return Property.create(name, const_1.TYPE.STRING, getter, setter, defaultValue);
    }
  }
  exports.Property = Property;
  class SemanticValue {
    constructor(d) {
      this.d = d;
    }
    static fromValueType(value, semType) {
      return new SemanticValue(api.grok_SemanticValue(value, semType));
    }
    get value() {
      return api.grok_SemanticValue_Get_Value(this.d);
    }
    set value(x) {
      api.grok_SemanticValue_Set_Value(this.d, x);
    }
    get semType() {
      return api.grok_SemanticValue_Get_SemType(this.d);
    }
    set semType(x) {
      api.grok_SemanticValue_Set_SemType(this.d, x);
    }
  }
  exports.SemanticValue = SemanticValue;
  class DateTime {
    constructor(d) {
      this.d = d;
    }
    static fromDate(date) {
      return DateTime.fromMillisecondsSinceEpoch(date.getTime());
    }
    static fromMillisecondsSinceEpoch(millisecondsSinceEpoch) {
      return new DateTime(api.grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
    }
  }
  exports.DateTime = DateTime;
});
