import {
    Credentials,
    DataConnection,
    DataJob,
    DataQuery,
    Entity, Group,
    Model,
    Notebook,
    Project, Script,
    TableInfo,
    User
} from "./entities";
import {ViewLayout} from "./view";
import {toJs} from "./wrappers";

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/js-api/projects-list}
 * */
export class Dapi {
    constructor() {}

    /** Retrieves entities from server by list of IDs
     *  @returns {List<Entity>} */
    getEntities(ids) {
        return new Promise((resolve, reject) => grok_Dapi_Entities_GetEntities(ids, (q) => {
            return resolve(q[0].map(toJs));
        }));
    }

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
    get notebooks() { return new HttpDataSource(grok_Dapi_Notebooks(), (a) => new Notebook(a)); }
    /** Predictive Models API endpoint
     *  @type {HttpDataSource<Model>} */
    get models() { return new HttpDataSource(grok_Dapi_Models(), (a) => new Model(a)); }
    /** Packages API endpoint
     *  @type {HttpDataSource<Package>} */
    get packages() { return new HttpDataSource(grok_Dapi_Packages(), (a) => new Package(a)); }
    /** View Layouts API endpoint
     *  @type {HttpDataSource<ViewLayout>} */
    get layouts() { return new HttpDataSource(grok_Dapi_Layouts(), (a) => new ViewLayout(a)); }
    /** Data Table Infos API endpoint
     *  @type {HttpDataSource<TableInfo>} */
    get tables() { return new HttpDataSource(grok_Dapi_Tables(), (a) => new TableInfo(a)); }
    /** Users API endpoint
     *  @type {UsersDataSource} */
    get users() { return new UsersDataSource(grok_Dapi_Users(), (a) => new User(a)); }
    /** Groups API endpoint
     *  @type {HttpDataSource<Group>} */
    get groups() { return new HttpDataSource(grok_Dapi_Groups(), (a) => new Group(a)); }
    /** Scripts API endpoint
     *  @type {HttpDataSource<Script>} */
    get scripts() { return new HttpDataSource(grok_Dapi_Scripts(), (a) => new Script(a)); }
    /** Projects API endpoint
     *  @type {HttpDataSource<Project>} */
    get projects() { return new HttpDataSource(grok_Dapi_Projects(), (a) => new Project(a)); }
    /** Users Data Storage API endpoint
     *  @type {UserDataStorage} */
    get userDataStorage() { return new UserDataStorage(); }
}


/**
 * Common functionality for handling collections of entities stored on the server.
 */
export class HttpDataSource {
    constructor(s, instance) {
        this.s = s;
        this.entityToJs = instance;
    }

    /** Returns all entities that satisfy the filtering criteria (see {@link filter}).
     *  See examples: {@link https://public.datagrok.ai/js/samples/js-api/projects-list}
     *  @param {Object} options
     *  @param {int} options.pageSize
     *  @param {int} options.pageNumber
     *  @param {String} options.filter
     *  @returns Promise<object[]>  */
    list(options) {
        if (options.pageSize !== undefined)
            this.by(options.pageSize);
        if (options.pageNumber !== undefined)
            this.page(options.pageNumber);
        if (options.filter !== undefined)
            this.filter(options.filter);
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_List(this.s, (q) => resolve(q.map(s))));
    }

    /** Returns an entity with the specified id.
     *  Throws an exception if an entity does not exist, or is not accessible in the current context.
     *  @param {string} id - GUID of the corresponding object
     *  @returns {object} - entity. */
    find(id) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_DataSource_Find(this.s, id, (q) => resolve(s(q[0]))));
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
     *  See examples: {@link https://public.datagrok.ai/js/samples/js-api/projects-list}
     *  @returns Promise<object[]>  */
    nextPage() {
        this.s = grok_DataSource_NextPage(this.s);
        return this;
    }

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
 * @extends HttpDataSource
 * */
export class UsersDataSource extends HttpDataSource {
    constructor(s, instance) {
        super(s, instance);
    }

    /** Returns current user
     * @returns {User} */
    current() {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_UsersDataSource_Current(this.s, (q) => resolve(s(q[0]))));
    }
}


/**
 * Functionality for handling Credentials collection from server and working with Credentials remote endpoint
 * @extends HttpDataSource
 * */
export class CredentialsDataSource extends HttpDataSource {
    constructor(s, instance) {
        super(s, instance);
    }

    /** Returns credentials for entity
     * @param {Entity} e
     * @returns {Credentials} */
    forEntity(e) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_CredentialsDataSource_ForEnrity(this.s, e.d, (c) => resolve(s(c[0]))));
    }
}

/**
 * Functionality for working with remote Users Data Storage
 * */
export class UserDataStorage {
    constructor() {}

    /** Saves a single value to Users Data Storage
     * @param {String} name Storage name
     * @param {String} key
     * @param {String} value
     * @param {boolean} currentUser Value should be available only for current user  */
    postValue(name, key, value, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser, () => resolve()));
    }

    /** Saves a map to Users Data Storage, will be appended to existing data
     * @param {String} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user  */
    post(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Post(name, data, currentUser, () => resolve()));
    }

    /** Saves a map to Users Data Storage, will replace existing data
     * @param {String} name Storage name
     * @param {Map} data
     * @param {boolean} currentUser Value should be available only for current user  */
    put(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Put(name, data, currentUser, () => resolve()));
    }

    /** Retrieves a map from Users Data Storage
     * @param {String} name Storage name
     * @param {boolean} currentUser get a value from a current user sotrage
     * @returns Map*/
    get(name, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Get(name, currentUser, (data) => resolve(data)));
    }

    /** Retrieves a single value from Users Data Storage
     * @param {String} name Storage name
     * @param {String} key Value key
     * @param {boolean} currentUser get a value from a current user storage
     * @returns Map*/
    getValue(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_GetValue(name, key, currentUser, (value) => resolve(value)));
    }

    /** Removes a single value from Users Data Storage
     * @param {String} name Storage name
     * @param {String} key Value key
     * @param {boolean} currentUser get a value from a current user storage*/
    remove(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Delete(name, key, currentUser, () => resolve()));
    }
}
