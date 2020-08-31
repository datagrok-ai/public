import {
    Credentials,
    DataConnection,
    DataJob,
    DataQuery,
    Entity, Group,
    Model,
    Notebook,
    Project, Script, ScriptEnvironment,
    TableInfo,
    User
} from "./entities";
import {ViewLayout} from "./view";
import {toJs} from "./wrappers";

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
 * */
export class Dapi {
    constructor() {}

    /** Retrieves entities from server by list of IDs
     *  @returns {List<Entity>} */
    getEntities(ids) {
        return new Promise((resolve, reject) => grok_Dapi_Entities_GetEntities(ids, (q) => {
            return resolve(q.map(toJs));
        }));
    }

    /** Entities API endpoint
     *  @type {EntitiesDataSource} */
    get entities() { return new EntitiesDataSource(grok_Dapi_Entities(), (a) => new Entity(a)); }

    /** Data Queries API endpoint
     *  @type {HttpDataSource<DataQuery>} */
    get queries() { return new HttpDataSource(grok_Dapi_Queries(), (a) => new DataQuery(a)); }

    /** Data Connections API endpoint
     *  @type {HttpDataSource<DataConnection>} */
    get connections() { return new HttpDataSource(grok_Dapi_Connections(), (a) => new DataConnection(a)); }

    /** Credentials API endpoint
     *  @type {CredentialsDataSource} */
    get credentials() { return new CredentialsDataSource(grok_Dapi_Credentials(), (a) => new Credentials(a)); }

    /** Data Jobs API endpoint
     *  @type {HttpDataSource<DataJob>} */
    get jobs() { return new HttpDataSource(grok_Dapi_Jobs(), (a) => new DataJob(a)); }

    /** Jupyter Notebooks API endpoint
     *  @type {HttpDataSource<Notebook>} */
    get notebooks() { return new HttpDataSource(grok_Dapi_Notebooks(), (a) => toJs(a)); }

    /** Predictive Models API endpoint
     *  @type {HttpDataSource<Model>} */
    get models() { return new HttpDataSource(grok_Dapi_Models(), (a) => new Model(a)); }

    /** Packages API endpoint
     *  @type {HttpDataSource<Package>} */
    get packages() { return new HttpDataSource(grok_Dapi_Packages(), (a) => new Package(a)); }

    /** View Layouts API endpoint
     *  @type {HttpDataSource<ViewLayout>} */
    get layouts() { return new LayoutsDataSource(grok_Dapi_Layouts(), (a) => new ViewLayout(a)); }

    /** Data Table Infos API endpoint
     *  @type {HttpDataSource<TableInfo>} */
    get tables() { return new HttpDataSource(grok_Dapi_Tables(), (a) => new TableInfo(a)); }

    /** Users API endpoint
     *  @type {UsersDataSource} */
    get users() { return new UsersDataSource(grok_Dapi_Users(), (a) => new User(a)); }

    /** Groups API endpoint
     *  @type {HttpDataSource<Group>} */
    get groups() { return new GroupsDataSource(grok_Dapi_Groups(), (a) => new Group(a)); }

    /** Scripts API endpoint
     *  @type {HttpDataSource<Script>} */
    get scripts() { return new HttpDataSource(grok_Dapi_Scripts(), (a) => new Script(a)); }

    /** Projects API endpoint
     *  @type {HttpDataSource<Project>} */
    get projects() { return new HttpDataSource(grok_Dapi_Projects(), (a) => new Project(a)); }

    /** Environments API endpoint
     *  @type {HttpDataSource<ScriptEnvironment>} */
    get environments() { return new HttpDataSource(grok_Dapi_Environments(), (a) => new ScriptEnvironment(a)); }

    /** Users Data Storage API endpoint
     *  @type {UserDataStorage} */
    get userDataStorage() { return new UserDataStorage(); }

    /** Proxies URL request via Datagrok server with same interface as "fetch".
     * @param {string} method
     * @param {string} url
     * @param {Object} headers
     * @param {Object} body
     * @returns {Promise<Object>} */
    async proxyFetch(method, url, headers = {}, body = {}) {
        headers['Accept'] = 'application/json';
        headers['original-url'] = `${url}`;
        headers['original-method'] = method;
        let params = {
            method: 'POST',
            headers: headers,
            body: JSON.stringify(body)
        };
        return await fetch('/api/connectors/proxy', params)
    }
}


/**
 * Common functionality for handling collections of entities stored on the server.
 * Works with Datagrok REST API, allows to get filtered and paginated lists of entities,
 * Can be extended with specific methods. (i.e. {@link UsersDataSource})
 */
export class HttpDataSource {

    /** @constructs HttpDataSource */
    constructor(s, instance) {
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
     *  @returns Promise<object[]>  */
    list(options= {}) {
        if (options.pageSize !== undefined)
            this.by(options.pageSize);
        if (options.pageNumber !== undefined)
            this.page(options.pageNumber);
        if (options.filter !== undefined)
            this.filter(options.filter);
        if (options.order !== undefined)
            this.order(options.order);
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_List(this.s, (q) => resolve(q.map(s))));
    }

    /** Returns fist entity that satisfy the filtering criteria (see {@link filter}).
     *  @returns Promise<object>  */
    first() {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_First(this.s, (q) => resolve(q.map(s))));
    }

    /** Returns an entity with the specified id.
     *  Throws an exception if an entity does not exist, or is not accessible in the current context.
     *  @param {string} id - GUID of the corresponding object
     *  @returns {Promise<object>} - entity. */
    find(id) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_Find(this.s, id, (q) => resolve(s(q)), (e) => reject(e)));
    }

    /** Saves an entity
     *  @param {Entity}  e
     *  @returns {Promise<Entity>} - entity. */
    save(e) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_Save(this.s, e.d, (q) => resolve(s(q)), (e) => reject(e)));
    }

    by(i) {
        this.s = grok_DataSource_By(this.s, i);
        return this;
    }

    page(i) {
        this.s = grok_DataSource_Page(this.s, i);
        return this;
    }

    /** Returns next page of all entities that satisfy the filtering criteria (see {@link filter}).
     *  Works only if pageSize was set during previous list() call
     *  See examples: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
     *  @returns Promise<object[]>  */
    nextPage() {
        this.s = grok_DataSource_NextPage(this.s);
        return this;
    }

    /** Applies filter to current request.
     *  Also can be set with {@link list} method "options" parameter
     *  See example: {@link https://public.datagrok.ai/js/samples/dapi/projects-list}
     *  Smart filter: {@link https://datagrok.ai/help/overview/smart-search}
     *  @param {string} w
     *  @returns HttpDataSource  */
    filter(w) {
        this.s = grok_DataSource_WhereSmart(this.s, w);
        return this;
    }

    /** Instructs data source to return results in the specified order.
     * @param {string} fieldName
     * @param {boolean} desc
     * @returns {HttpDataSource} */
    order(fieldName, desc = false) {
        this.s = grok_DataSource_Order(this.s, fieldName, desc);
        return this;
    }
}


/**
 * Functionality for handling Users collection from server and working with Users remote endpoint
 * Allows to load current user and list of all Datagrok users with filtering and pagination
 * See example: {@link https://public.datagrok.ai/js/samples/dapi/who-am-i}
 * @extends HttpDataSource
 * */
export class UsersDataSource extends HttpDataSource {
    /** @constructs UsersDataSource*/
    constructor(s, instance) {
        super(s, instance);
    }

    /** Returns current user
     * @returns {Promise<User>} */
    current() {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_UsersDataSource_Current(this.s, (q) => resolve(s(q))));
    }
}

/**
 * Functionality for handling groups collection from server
 * Allows to manage {@link Group}
 * @extends HttpDataSource
 * */
export class GroupsDataSource extends HttpDataSource {
    /** @constructs CredentialsDataSource*/
    constructor(s, instance) {
        super(s, instance);
    }

    /** Saves a group with relations
     *  @param {Group}  e
     *  @returns {Promise<Group>} - Group. */
    saveRelations(e) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_GroupsDataSource_Save(this.s, e.d, (q) => resolve(s(q)), (e) => reject(e)));
    }

}


/**
 * Functionality for handling entities collection from server
 * Allows to manage {@link Entity}
 * @extends HttpDataSource
 * */
export class EntitiesDataSource extends HttpDataSource {
    /** @constructs CredentialsDataSource*/
    constructor(s, instance) {
        super(s, instance);
    }

    /** Allows to set properties for entities
     * @param {List<Map>} props
     * @returns {Promise} */
    saveProperties(props) {
        return new Promise((resolve, reject) => grok_EntitiesDataSource_SaveProperties(this.s, props, (_) => resolve()));
    }

    /** Returns entity properties
     * @param {Entity} entity
     * @returns {Promise<Map>} props */
    getProperties(entity) {
        return new Promise((resolve, reject) => grok_EntitiesDataSource_GetProperties(this.s, entity.d, (p) => resolve(p)));
    }
}

/**
 * Functionality for handling credentials collection from server and working with credentials remote endpoint
 * Allows to manage {@link Credentials}
 * See also: {@link https://datagrok.ai/help/govern/security}
 * @extends HttpDataSource
 * */
export class CredentialsDataSource extends HttpDataSource {
    /** @constructs CredentialsDataSource*/
    constructor(s, instance) {
        super(s, instance);
    }

    /** Returns credentials for entity
     * @param {Entity} e
     * @returns {Promise<Credentials>} */
    forEntity(e) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_CredentialsDataSource_ForEntity(this.s, e.d, (c) => resolve(s(c))));
    }
}

/**
 * Functionality for handling layouts collection from server
 * Allows to manage {@link ViewLayout}
 * @extends HttpDataSource
 * */
export class LayoutsDataSource extends HttpDataSource {
    /** @constructs CredentialsDataSource*/
    constructor(s, instance) {
        super(s, instance);
    }

    /** Returns layouts that applicable to the table
     * @param {DataFrame} t
     * @returns {Promise<List<ViewLayout>>} */
    getApplicable(t) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_LayoutsDataSource_Applicable(this.s, t.d, (q) => resolve(q.map(s))));
    }
}

/**
 * Functionality for working with remote Users Data Storage
 * Remote storage allows to save key-value pairs on the Datagrok server for further use
 * */
export class UserDataStorage {
    constructor() {}

    /** Saves a single value to Users Data Storage
     * @param {string} name Storage name
     * @param {string} key
     * @param {string} value
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    postValue(name, key, value, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser, () => resolve()));
    }

    /** Saves a map to Users Data Storage, will be appended to existing data
     * @param {string} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    post(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Post(name, data, currentUser, () => resolve()));
    }

    /** Saves a map to Users Data Storage, will replace existing data
     * @param {string} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user
     * @returns {Promise}*/
    put(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Put(name, data, currentUser, () => resolve()));
    }

    /** Retrieves a map from Users Data Storage
     * @param {string} name Storage name
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise<Map>} */
    get(name, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Get(name, currentUser, (data) => resolve(data)));
    }

    /** Retrieves a single value from Users Data Storage
     * @param {string} name Storage name
     * @param {string} key Value key
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise<string>} */
    getValue(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_GetValue(name, key, currentUser, (value) => resolve(value)));
    }

    /** Removes a single value from Users Data Storage
     * @param {string} name Storage name
     * @param {string} key Value key
     * @param {boolean} currentUser get a value from a current user storage
     * @returns {Promise} */
    remove(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Delete(name, key, currentUser, () => resolve()));
    }
}
