import {SEMTYPE, TYPE} from "./const";
import {toJs} from "./wrappers";


/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {

    /** Entity ID (GUID)
     *  @type {string} */
    get id(): string;

    set id(x: string)

    /** Entity name
     *  @type {string} */
    get name(): string

    set name(x: string)

    /** Entity path
     *  @type {string} */
    get path(): string
    
    /** Returns entity properties
     * @returns {Promise<Map>} props */
    getProperties(): Promise<Map<any, any>>
    /** Allows to set properties for entity
     * @param {Map} props
     * @returns Promise */
    setProperties(props: Map<any, any>): Promise<void>
}

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {

    /** First name
     * @type {string} */
    get firstName(): string

    /** Last name
     * @type {string} */
    get lastName(): string

    /** Email
     * @type {string} */
    get email(): string;

    /** Picture URL
     * @type {string} */
    get picture(): string

    /** Login
     *  @type {string} */
    get login(): string

    static fromId(id: string): User

    static current(): User

    /** */
    toMarkup(): string

    /** Security Group
     * @type Group
     */
    get group(): Group
}

/** Represents a function
 * @extends Entity
 * {@link https://datagrok.ai/help/overview/functions/function}
 * */
export class Func extends Entity {

    /** Returns {@link FuncCall} object in a stand-by state
     * @param {object} parameters
     * @returns {FuncCall} */
    prepare(parameters?: object | null): Function
}

/** Represents a project */
export class Project extends Entity {

    /** Entity name
     *  @type {string} */
    get description(): string

    get isDirty(): boolean

    get isEmpty(): boolean

    toMarkup(): string
}

/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-query}
 * */
export class DataQuery extends Func {

    /** Query text
     *  @type {string} */
    get query(): string
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access/data-job}
 * */
export class DataJob extends Func {
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/data-connection}
 * */
export class DataConnection extends Entity {

    /** Collection of parameters: server, database, endpoint, etc.
     *  @type {object} */
    get parameters(): object
}

/** Represents a predictive model
 * @extends Entity
 * {@link https://datagrok.ai/help/learn/predictive-modeling-info}
 * */
export class Model extends Entity {
}

/** @extends Entity
 * Represents a Jupyter notebook
 * {@link https://datagrok.ai/help/develop/jupyter-notebook}
 * */
export class Notebook extends Entity {

    /** Description
     * @type {string} */
    get description(): string
    set description(d)


    /** Environment name
     * @type {string} */
    get environment(): string

    set environment(e: string)

    /** Create Notebook on server for edit.
     * @returns {Promise<string>} Current notebook's name */
    edit(): Promise<string>

    /** Converts Notebook to HTML code
     * @returns {Promise<string>} */
    toHtml(): Promise<string>
}

/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
    /** @constructs TableInfo */
    constructor(d: any)
}

/** @extends Entity
 * Allows for files handling in JS-based info panels
 * {@link https://datagrok.ai/help/discover/info-panels}
 * */
export class FileInfo extends Entity {
    /** @constructs FileInfo */
    constructor(d: any)
    /** Returns path, i.e. `geo/dmv_offices.csv`
     * @type {string} */
    get path() : string
    /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv`
     * @type {string} */
    get fullPath() : string
    /** Returns file extension, i.e. `csv`
     * @type {string} */
    get extension() : string
    /** Returns file name, i.e. `dmv_offices.csv`
     * @type {string} */
    get fileName() : string
}

/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
    /** @constructs Group */
    constructor(d: any)


    static create(name: string): Group

    /** Adds a member to the group
     * @param {Group} m */
    addMember(m: Group): void

    /** Adds an admin member to the group
     * @param {Group} m */
    addAdminMember(m: Group): void

    /** Removes a member from the group
     * @param {Group} m */
    removeMember(m: Group): void

    /** Adds the group to another one
     * @param {Group} m */
    includeTo(m: Group): void

    /** Adds the group to another one as an admin
     * @param {Group} m */
    includeAdminTo(m: Group): void

    /** Removes membership from another group
     * @param {Group} m */
    excludeFrom(m: Group): void

    /** Returns list of groups that belong to group, with no admin permissions
     * @type {List<Group>} */
    get members(): Group[]

    /** Returns list of groups that belong to group, with admin permissions
     * @type {List<Group>} */
    get adminMembers(): Group[]

    /** Returns list of groups that group belongs to, with no admin permissions
     * @type {List<Group>} */
    get memberships(): Group[]

    /** Returns list of groups that group belongs to, with admin permissions
     * @type {List<Group>} */
    get adminMemberships(): Group[]

    /** Personal user group
     * @type {boolean} */
    get personal(): boolean
    set personal(e: boolean)

    /** Hidden group
     * @type {boolean} */
    get hidden(): boolean
    set hidden(e: boolean);

}

/** @extends Func
 * Represents a Script
 * */
export class Script extends Func {
    /** @constructs Script */
    constructor(d: any)

    static create(script: string): Script
    
    /** Script
     * @type {string} */
    get script() : string
    set script(s)
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/govern/security}
 *  */
export class Credentials extends Entity {
    constructor(d: any)

    /** Collection of parameters: login, password, API key, etc.
     *  @type {object} */
    get parameters(): object
}

/** Represents a script environment */
export class ScriptEnvironment extends Entity {
    constructor(d: any)

    /** Environment yaml file content
     * @type {string} */
    get environment(): string

    /** Create instance of ScriptEnvironment
     * @param {string} name
     * @returns {ScriptEnvironment}
     * */
    static create(name: string): ScriptEnvironment

    /** Setup environment */
    setup(): Promise<void>
}

/**
 * Represents a package, which is a unit of distribution of content in the Datagrok platform.
 */
export class Package {
    webRoot?: string

    constructor(webRoot?: string)

    /** Override init() method to provide package-specific initialization.
     * It is guaranteed to get called exactly once before the execution of any function below.
     * */

    /*async*/
    init(): Promise<null>

    /** Returns credentials for package
     * @returns {Promise<Credentials>} */
    getCredentials(): Promise<Credentials>
}


type PropertyGetter = (item: any) => any
type PropertySetter = (item: any) => void

/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property {
    constructor(d: any)

    /** Property getter is a function that acccepts one parameter (item)
     * and returns the property value. */
    get get(): PropertyGetter

    set get(x: PropertyGetter)

    /** Property setter */
    get set(): PropertySetter

    set set(x: PropertySetter);

    /** Property name
     *  @type {string} */
    get name(): string

    set name(s: string)

    /** Property type
     *  @type {string} */
    get propertyType(): string

    set propertyType(s: string)

    /** Semantic type
     *  @type {string} */
    get semType(): SEMTYPE | string

    set semType(s: SEMTYPE | string)

    /** Description */
    get description(): string

    set description(s: string)

    /** Default value */
    get defaultValue(): any

    set defaultValue(s: any)

    /** List of possible values of that property.
     *  PropertyGrid will use it to populate combo boxes.
     *  @returns ArrayList<string>*/
    get choices(): string[]

    set choices(x: string[])

    /** Creates a property
     * @param {string} name
     * @param {TYPE|Type} type
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static create(name: string, type: TYPE | string,
                  getter: PropertyGetter, setter: PropertySetter,
                  defaultValue?: any | null): Property

    /** Creates an integer property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static int(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue?: number | null): Property

    /** Creates a float property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static float(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: number | null): Property

    /** Creates a string property
     * @param {string} name
     * @param {function} getter
     * @param {function} setter
     * @param {object} defaultValue
     * @returns Property*/
    static string(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: string | null): Property
}

export class SemanticValue {
    constructor(d: any)

    get value(): any

    set value(x: any)

    get semType(): SEMTYPE

    set semType(x: SEMTYPE)

    static fromValueType(value: any, semType: SEMTYPE): SemanticValue
}

export class DateTime {
    constructor(d: any)

    static fromDate(date: Date): DateTime

    static fromMillisecondsSinceEpoch(millisecondsSinceEpoch: number): DateTime
}