import * as ui from "./../ui.js";
import {Functions} from "./functions";
import {TYPE} from "./const";
import {toJs} from "./wrappers";


/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {
    /** @constructs Entity*/
    constructor(d) { this.d = d;}

    /** Entity ID (GUID)
     *  @type {string} */
    get id() { return grok_Entity_Get_Id([this.d]); }
    set id(x) { return grok_Entity_Set_Id([this.d], x); }

    /** Entity name
     *  @type {string} */
    get name() { return grok_Entity_Get_Name([this.d]); }
    set name(x) { return grok_Entity_Set_Name([this.d], x); }

    /** Entity path
     *  @type {string} */
    get path() { return grok_Entity_Path([this.d]); }
}

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {
    /** @constructs User*/
    constructor(d) { super(d); }

    static fromId(id) { return new User(grok_User_From_Id(id)); }

    static current() { return new User(grok_User()); }

    /** First name
     * @type {string} */
    get firstName() { return grok_User_Get_FirstName(this.d); }

    /** Last name
     * @type {string} */
    get lastName() { return grok_User_Get_LastName(this.d); }

    /** Email
     * @type {string} */
    get email() { return grok_User_Get_Email(this.d); }

    /** Picture URL
     * @type {string} */
    get picture() { return grok_User_Get_Picture(this.d); }

    /** Login
     *  @type {string} */
    get login() { return grok_User_Get_Login(this.d); }

    /** */
    toMarkup() { return grok_User_ToMarkup(this.d); }
}

/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/overview/functions/function}
 * */
export class Func extends Entity {
    /** @constructs Func*/
    constructor(d) {super(d)}

    /** Returns {@link FuncCall} object in a stand-by state
     * @param {object} parameters
     * @returns {FuncCall} */
    prepare(parameters = {}) { return grok_Func_Prepare(this.d, parameters); };
}

/** Represents a project */
export class Project extends Entity {
    constructor(d) { super(d); }

    /** Entity name
     *  @type {string} */
    get description() { return grok_Project_Description(this.d); }
    /** Entity name
     *  @type {string} */
    get isDirty() { return grok_Project_IsDirty(this.d); }
    /** Entity name
     *  @type {string} */
    get isEmpty() { return grok_Project_IsEmpty(this.d); }

    toMarkup() { return grok_Project_ToMarkup(this.d); }
}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-query}
 * */
export class DataQuery extends Func {
    /** @constructs DataQuery*/
    constructor(d) { super(d); }

    /** Query text
     *  @type {string} */
    get query() { return grok_Query_Query([this.d]); }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-job}
 * */
export class DataJob extends Func {
    /** @constructs DataJob */
    constructor(d) { super(d); }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/data-connection}
 * */
export class DataConnection extends Entity {
    /** @constructs DataConnection */
    constructor(d) { super(d); }

    /** Collection of parameters: server, database, endpoint, etc.
     *  @type {object} */
    get parameters() { return grok_DataConnection_Parameters([this.d]); }
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
 * */
export class Model extends Entity {
    /** @constructs Model */
    constructor(d) { super(d); }
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/compute/jupyter-notebook}
 * */
export class Notebook extends Entity {
    /** @constructs Notebook */
    constructor(d) { super(d); }

    /** Converts Notebook to HTML code
     * @returns {Promise<string>} */
    toHtml() { return new Promise((resolve, reject) => grok_Notebook_ToHtml(this.d, (html) => resolve(html))); }
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
    /** @constructs TableInfo */
    constructor(d) { super(d); }
}

/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
    /** @constructs Group */
    constructor(d) { super(d); }
}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
    /** @constructs Script */
    constructor(d) { super(d); }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
 *  */
export class Credentials extends Entity {
    constructor(d) { super(d); }

    /** Collection of parameters: login, password, API key, etc.
     *  @type {object} */
    get parameters() { return grok_Credentials_Parameters([this.d]); }
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package {
    constructor(webRoot) {
        this.webRoot = webRoot;
        this.name = "";
    }

    /** Override init() method to provide package-specific initialization.
     * It is guaranteed to get called exactly once before the execution of any function below.
     * */
    /*async*/ init() { return Promise.resolve(null); }

    /** Returns credentials for package
     * @returns {Credentials} */
    getCredentials() {
        return new Promise((resolve, reject) => grok_Package_Get_Credentials(this.name, (c) => resolve(toJs(c))));
    }
}

/**
 * Override this class, and {@link register} an instance to integrate the platform with custom
 * types and objects.
 *
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/meta/meta}
 * */
export class JsEntityMeta {

    /** Type of the object that this meta handles. */
    get type() { throw 'Not defined.'; }

    /**
     * Override this method to check whether this meta class should handle the specified object.
     * @param x - specified object.
     * @returns {boolean}
     * */
    isApplicable(x) { throw 'Not defined.'; }

    /** String representation of the [item], by default item.toString().
     * @param x - item
     * @returns {string} */
    getCaption(x) { return `${x}`; };

    /** Renders icon for the item.
     * @param x - item
     * @returns {Element} */
    renderIcon(x) { return ui.divText(this.getCaption(x)); }

    /** Renders markup for the item.
     * @param x - item
     * @returns {Element} */
    renderMarkup(x) { return ui.divText(this.getCaption(x)); }

    /** Renders tooltip for the item.
     * @param x - item
     * @returns {Element} */
    renderTooltip(x) { return ui.divText(this.getCaption(x)); }

    /** Renders card div for the item.
     * @param x - item
     * @returns {Element} */
    renderCard(x) { return ui.divText(this.getCaption(x)); }

    /** Renders properties list for the item.
     * @param x - item
     * @returns {Element} */
    renderProperties(x) { return ui.divText(this.getCaption(x)); }

    /** Renders view for the item.
     * @param x - item
     * @returns {Element} */
    renderView(x) {return this.renderProperties(x); }

    /** Gets called once upon the registration of meta export class. */
    init() {}

    static register(meta) { grok_Meta_Register(meta); }

    /**
     * Registers a function that takes applicable objects an the only argument.
     * It will be suggested to run in the context menu for that object, and
     * also in the "Actions" pane on the property panel.
     *
     * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     *
     * @param {string} name - function name
     * @param run - a function that takes exactly one parameter
     * */
    registerParamFunc(name, run) {
        new Functions().registerParamFunc(name, this.type, run, this.isApplicable);
    }
}


/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property {
    constructor(d) { this.d = d; }

    /** Property getter is a function that acccepts one parameter (item)
     * and returns the property value. */
    get get() { return grok_Property_Get_Get(this.d); }
    set get(x) { return grok_Property_Set_Get(this.d, x); }

    /** Property setter */
    get set() { return grok_Property_Get_Set(this.d); }
    set set(x) { return grok_Property_Set_Set(this.d, x); }

    /** Property name
     *  @type {string} */
    get name() { return grok_Property_Get_Name(this.d); }
    set name(s) { grok_Property_Set_Name(this.d, s); }

    /** Property type
     *  @type {string} */
    get propertyType() { return grok_Property_Get_PropertyType(this.d); }
    set propertyType(s) { grok_Property_Set_PropertyType(this.d, s); }

    /** Semantic type
     *  @type {string} */
    get semType() { return grok_Property_Get_SemType(this.d); }
    set semType(s) { grok_Property_Set_SemType(this.d, s); }

    /** Description */
    get description() { return grok_Property_Get_Description(this.d); }
    set description(s) { grok_Property_Set_Description(this.d, s); }

    /** Default value */
    get defaultValue() { return grok_Property_Get_DefaultValue(this.d); }
    set defaultValue(s) { grok_Property_Set_DefaultValue(this.d, s); }

    /** List of possible values of that property.
     *  PropertyGrid will use it to populate combo boxes.
     *  @returns ArrayList<string>*/
    get choices() { return grok_Property_Get_Choices(this.d); }
    set choices(x) { return grok_Property_Set_Choices(this.d, x); }

    /** Creates a property
     * @param {string} name
     * @param {TYPE|Type} type
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static create(name, type, getter, setter, defaultValue = null) { return new Property(grok_Property(name, type, getter, setter, defaultValue)); }

    /** Creates an integer property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static int(name, getter, setter, defaultValue) { return Property.create(name, TYPE.INT, getter, setter, defaultValue); }

    /** Creates a float property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static float(name, getter, setter, defaultValue) { return Property.create(name, TYPE.FLOAT, getter, setter, defaultValue); }

    /** Creates a string property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static string(name, getter, setter, defaultValue) { return Property.create(name, TYPE.STRING, getter, setter, defaultValue); }
}

export class SemanticValue {
    constructor(d) { this.d = d; }

    static fromValueType(value, semType) { return new SemanticValue(grok_SemanticValue(value, semType)); }

    get value() { return grok_SemanticValue_Get_Value(this.d); }
    set value(x) { grok_SemanticValue_Set_Value(this.d, x); }

    get semType() { return grok_SemanticValue_Get_SemType(this.d); }
    set semType(x) { grok_SemanticValue_Set_SemType(this.d, x); }
}

export class DateTime {
    constructor(d) { this.d = d; }

    static fromDate(date) { return DateTime.fromMillisecondsSinceEpoch(date.getTime()); }
    static fromMillisecondsSinceEpoch(millisecondsSinceEpoch) {
        return new DateTime(grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
    }
}
