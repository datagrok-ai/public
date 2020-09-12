import {
    Credentials,
    DataConnection,
    DataJob,
    DataQuery,
    Entity,
    Group,
    Model,
    Notebook,
    Package,
    Project, Property,
    Script,
    ScriptEnvironment,
    TableInfo,
    User
} from "./entities";
import {ViewLayout} from "./view";
import {DataFrame} from "./dataframe";

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
 * */
export class Dapi {
    constructor()

    get entities(): EntitiesDataSource;
    /** Data Queries API endpoint
     *  @type {HttpDataSource<DataQuery>} */
    get queries(): HttpDataSource<DataQuery>

    /** Data Connections API endpoint
     *  @type {HttpDataSource<DataConnection>} */
    get connections(): HttpDataSource<DataConnection>

    /** Credentials API endpoint
     *  @type {CredentialsDataSource} */
    get credentials(): CredentialsDataSource

    /** Data Jobs API endpoint
     *  @type {HttpDataSource<DataJob>} */
    get jobs(): HttpDataSource<DataJob>

    /** Jupyter Notebooks API endpoint
     *  @type {HttpDataSource<Notebook>} */
    get notebooks(): HttpDataSource<Notebook>

    /** Predictive Models API endpoint
     *  @type {HttpDataSource<Model>} */
    get models(): HttpDataSource<Model>

    /** Packages API endpoint
     *  @type {HttpDataSource<Package>} */
    get packages(): HttpDataSource<Package>

    /** View Layouts API endpoint
     *  @type {HttpDataSource<ViewLayout>} */
    get layouts(): HttpDataSource<ViewLayout>

    /** Data Table Infos API endpoint
     *  @type {HttpDataSource<TableInfo>} */
    get tables(): HttpDataSource<TableInfo>

    /** Users API endpoint
     *  @type {UsersDataSource} */
    get users(): UserDataStorage

    get permissions(): PermissionsDataSource
    
    /** Groups API endpoint
     *  @type {HttpDataSource<Group>} */
    get groups(): GroupsDataSource

    /** Scripts API endpoint
     *  @type {HttpDataSource<Script>} */
    get scripts(): HttpDataSource<Script>

    /** Projects API endpoint
     *  @type {HttpDataSource<Project>} */
    get projects(): HttpDataSource<Project>

    /** Environments API endpoint
     *  @type {HttpDataSource<ScriptEnvironment>} */
    get environments(): HttpDataSource<ScriptEnvironment>

    /** Users Data Storage API endpoint
     *  @type {UserDataStorage} */
    get userDataStorage(): UserDataStorage

    /** Retrieves entities from server by list of IDs
     *  @returns {List<Entity>} */
    getEntities(ids: string): Promise<Entity[]>

    /** Proxies URL request via Datagrok server with same interface as "fetch".
     * @param {string} method
     * @param {string} url
     * @param {Object} headers
     * @param {Object} body
     * @returns {Promise<Object>} */
    proxyFetch(method: string, url: string, headers?: Object, body?: Object): Promise<object>
}


/**
 * Common functionality for handling collections of entities stored on the server.
 * Works with Datagrok REST API, allows to get filtered and paginated lists of entities,
 * Can be extended with specific methods. (i.e. {@link UsersDataSource})
 */
export class HttpDataSource<T> {

    /** @constructs HttpDataSource */
    constructor(s: any, instance: any)

    /** Returns all entities that satisfy the filtering criteria (see {@link filter}).
     *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
     *  Smart filter: {@link https://datagrok.ai/help/overview/smart-search}
     *  @param {Object} options
     *  @param {int} options.pageSize
     *  @param {int} options.pageNumber
     *  @param {string} options.order
     *  @param {string} options.filter
     *  @returns Promise<object[]>  */
    list(options?: {
        pageSize?: number,
        pageNumber?: number,
        order?: keyof T,
        filter?: string
    }): Promise<T[]>

    first(): Promise<T>
    
    /** Returns an entity with the specified id.
     *  Throws an exception if an entity does not exist, or is not accessible in the current context.
     *  @param {string} id - GUID of the corresponding object
     *  @returns {Promise<object>} - entity. */
    find(id: string): Promise<T>

    save(e: T): Promise<void>
    
    delete(e: T): Promise<void>
    
    by(i: number): this

    page(i: number): this

    /** Returns next page of all entities that satisfy the filtering criteria (see {@link filter}).
     *  Works only if pageSize was set during previous list() call
     *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
     *  @returns Promise<object[]>  */
    nextPage(): Promise<T[]>

    /** Applies filter to current request.
     *  Also can be set with {@link list} method "options" parameter
     *  See example: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
     *  Smart filter: {@link https://datagrok.ai/help/overview/smart-search}
     *  @param {string} w
     *  @returns HttpDataSource  */
    filter(w: string): this

    /** Instructs data source to return results in the specified order.
     * @param {string} fieldName
     * @param {boolean} desc
     * @returns {HttpDataSource} */
    order(fieldName: keyof T, desc?: boolean): this
}


/**
 * Functionality for handling Users collection from server and working with Users remote endpoint
 * Allows to load current user and list of all Datagrok users with filtering and pagination
 * See example: {@link https://public.datagrok.ai/js/samples/dapi/who-am-i}
 * @extends HttpDataSource
 * */
export class UsersDataSource extends HttpDataSource<User> {
    /** @constructs UsersDataSource*/
    constructor(s: any, instance: any)

    /** Returns current user
     * @returns {Promise<User>} */
    current(): Promise<User>
}

/**
 * Functionality for handling groups collection from server
 * Allows to manage {@link Group}
 * @extends HttpDataSource
 * */
export class GroupsDataSource extends HttpDataSource<Group> {
    /** @constructs CredentialsDataSource*/
    constructor(s: any, instance: any)

    /** Creates a new group
     *  @param {String}  name
     *  @returns {Promise<Group>} - Group. */
    createNew(name: string): Promise<Group>

    /** Returns group user
     *  @param {Group} group
     *  @returns {Promise<Group>} - Group. */
    getUser(group: Group) : Promise<Group>

    /** Adds a member to the group
     * @param {Group} g
     * @param {Group} m
     * @returns Promise */
    addMember(g: Group, m: Group): Promise<Group>

    /** Adds an admin member to the group
     * @param {Group} g
     * @param {Group} m
     * @returns Promise */
    addAdminMember(g: Group, m: Group): Promise<Group>

    /** Removes a member from the group
     * @param {Group} g
     * @param {Group} m
     * returns Promise */
    removeMember(g: Group, m: Group): Promise<Group>

    /** Adds the group to another one
     * @param {Group} g
     * @param {Group} parent
     * returns Promise */
     includeTo(g: Group, parent: Group): Promise<Group>


    /** Adds the group to another one as admin
     * @param {Group} g
     * @param {Group} parent
     * returns Promise */
    includeAdminTo(g: Group, parent: Group): Promise<Group>

    /** Removes a membership from the group
     * @param {Group} g
     * @param {Group} parent
     * returns Promise */
    excludeFrom(g: Group, parent: Group): Promise<Group>

    /** Saves a group with relations
     *  @param {Group}  e
     *  @returns {Promise<Group>} - Group. */
    saveRelations(e: Group): Promise<Group>

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
    saveProperties(props: Map<Property, any>[]): Promise<void>

    /** Returns entity properties
     * @param {Entity} entity
     * @returns {Promise<Map>} props */
    getProperties(entity: Entity): Promise<Map<Property, any>>
}


/**
 * Functionality for handling credentials collection from server and working with credentials remote endpoint
 * Allows to manage {@link Credentials}
 * See also: {@link https://datagrok.ai/help/govern/security}
 * @extends HttpDataSource
 * */
export class CredentialsDataSource extends HttpDataSource<Credentials> {
    /** @constructs CredentialsDataSource*/
    constructor(s: any, instance: any)

    /** Returns credentials for entity
     * @param {Entity} e
     * @returns {Credentials} */
    forEntity(e: Entity): Promise<Credentials>
}

export class PermissionsDataSource {
    constructor();

    /** Gets all the permission granted on entity
     * @param {Entity} e
     * @returns {Promise<Map>} permissions
     * */
    get(e: Entity): Promise<Map<Group, boolean>>

    /** Grants permission on entity to the group
     * @param {Entity} e
     * @param {Group} g
     * @param {boolean} edit allow to edit entity
     * @returns {Promise}
     * */
    grant(e: Entity, g: Group, edit: boolean): Promise<void>
    /** Revokes permission on entity from the group
     * @param {Entity} e
     * @param {Group} g
     * @returns {Promise}
     * */
    revoke(g: Entity, e: Group): Promise<void>
}

/**
 * Functionality for handling layouts collection from server
 * Allows to manage {@link ViewLayout}
 * @extends HttpDataSource
 * */
export class LayoutsDataSource extends HttpDataSource<ViewLayout> {
    /** @constructs CredentialsDataSource*/
    constructor(s: any, instance: any)

    /** Returns layouts that applicable to the table
     * @param {DataFrame} t
     * @returns {Promise<List<ViewLayout>>} */
    getApplicable(t: DataFrame): Promise<ViewLayout[]>
}

/**
 * Functionality for working with remote Users Data Storage
 * Remote storage allows to save key-value pairs on the Datagrok server for further use
 * */
export class UserDataStorage {

    /** Saves a single value to Users Data Storage
     * @param {string} name Storage name
     * @param {string} key
     * @param {string} value
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    postValue(name: string, key: string, value: string, currentUser?: boolean): Promise<void>

    /** Saves a map to Users Data Storage, will be appended to existing data
     * @param {string} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    post(name: string, data: Map<string, string>, currentUser?: boolean): Promise<void>

    /** Saves a map to Users Data Storage, will replace existing data
     * @param {string} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    put(name: string, data: Map<string, string>, currentUser?: boolean): Promise<void>

    /** Retrieves a map from Users Data Storage
     * @param {string} name Storage name
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise<Map>} */
    get(name: string, currentUser?: boolean): Promise<Map<string, string>>

    /** Retrieves a single value from Users Data Storage
     * @param {string} name Storage name
     * @param {string} key Value key
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise<string>} */
    getValue(name: string, key: string, currentUser?: boolean): Promise<string>

    /** Removes a single value from Users Data Storage
     * @param {string} name Storage name
     * @param {string} key Value key
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise} */
    remove(name: string, key: string, currentUser?: boolean): Promise<void>
}