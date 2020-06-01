import {Credentials, DataConnection, DataQuery, Entity, Project, User} from "./entities";
import {ViewLayout} from "./view";
import {toJs} from "./wrappers";

/**
 * Exposes Datagrok's server-side functionality.
 *
 * See examples: {@link https://public.datagrok.ai/js/samples/grok-api/projects-list}
 * */
export class Dapi {
    constructor() {}

    getEntities(ids) {
        return new Promise((resolve, reject) => grok_Dapi_Entities_GetEntities(ids, (q) => {
            return resolve(q[0].map(toJs));
        }));
    }

    get queries() { return new HttpDataSource(grok_Dapi_Queries(), (a) => new DataQuery(a)); }
    get connections() { return new HttpDataSource(grok_Dapi_Connections(), (a) => new DataConnection(a)); }
    get credentials() { return new CredentialsDataSource(grok_Dapi_Credentials(), (a) => new Credentials(a)); }
    get jobs() { return new HttpDataSource(grok_Dapi_Jobs(), (a) => new Entity(a)); }
    get notebooks() { return new HttpDataSource(grok_Dapi_Notebooks(), (a) => new Entity(a)); }
    get models() { return new HttpDataSource(grok_Dapi_Models(), (a) => new Entity(a)); }
    get packages() { return new HttpDataSource(grok_Dapi_Packages(), (a) => new Entity(a)); }
    get layouts() { return new HttpDataSource(grok_Dapi_Layouts(), (a) => new ViewLayout(a)); }
    get tables() { return new HttpDataSource(grok_Dapi_Tables(), (a) => new Entity(a)); }
    get users() { return new UsersDataSource(grok_Dapi_Users(), (a) => new User(a)); }
    get groups() { return new HttpDataSource(grok_Dapi_Groups(), (a) => new Entity(a)); }
    get scripts() { return new HttpDataSource(grok_Dapi_Scripts(), (a) => new Entity(a)); }
    get projects() { return new HttpDataSource(grok_Dapi_Projects(), (a) => new Project(a)); }
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
     *  See examples: {@link https://public.datagrok.ai/js/samples/grok-api/projects-list}
     *  @returns Promise<object[]>  */
    list() {
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


export class UsersDataSource extends HttpDataSource {
    constructor(s, instance) {
        super(s, instance);
    }

    current() {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_UsersDataSource_Current(this.s, (q) => resolve(s(q[0]))));
    }
}


export class CredentialsDataSource extends HttpDataSource {
    constructor(s, instance) {
        super(s, instance);
    }

    forEntity(e) {
        let s = this.entityToJs;
        return new Promise((resolve, reject) => grok_CredentialsDataSource_ForEnrity(this.s, e.d, (c) => resolve(s(c[0]))));
    }
}


export class UserDataStorage {
    constructor() {}

    postValue(name, key, value, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser, () => resolve()));
    }

    post(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Post(name, data, currentUser, () => resolve()));
    }

    put(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Put(name, data, currentUser, () => resolve()));
    }

    get(name, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Get(name, currentUser, (data) => resolve(data)));
    }

    getValue(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_GetValue(name, key, currentUser, (value) => resolve(value)));
    }

    remove(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Delete(name, key, currentUser, () => resolve()));
    }
}
