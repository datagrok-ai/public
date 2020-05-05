//import { fromEventPattern } from 'rxjs';

/** Obsolete and will be removed soon; use "grok". */
gr = new Grok();

grok = new Grok();

/** Grok entry point, use it to get access to top-level views, tables, methods, etc. */
class Grok {

    constructor() {
        this.events = new GrokEvents();
    }

    /** Current table */
    get t() { return new DataFrame(grok_CurrentTable()); }

    /** Current view */
    set v(view) { grok_Set_CurrentView(view.d); }
    get v() { return View.fromDart(grok_Get_CurrentView()); }

    /** Current project */
    get project() { return new Project(grok_Project()); }

    /** List of table names that are currently open */
    get tableNames() { return grok_TableNames(); }

    /** List of currently open tables table */
    get tables() { return this.tableNames.map(this.tableByName); }

    get balloon() { return new Balloon(); }

    get dapi() { return new Dapi(); }

    /** Current user */
    get user() { return new User(grok_User()); }

    get sidebar() { return new TabControl(grok_Get_Sidebar()); }

    /** Presentaion mode **/
    get presentationMode() { return grok_Get_PresentationMode(); }
    set presentationMode(x) { return grok_Set_PresentationMode(x); }

    /** Hide dock tabs in presentaion mode **/
    get hideTabsInPresentationMode() { return grok_Get_HideTabsInPresentationMode(); }
    set hideTabsInPresentationMode(x) { return grok_Set_HideTabsInPresentationMode(x); }

    get functions() { return new Functions(); }

    dockElement(e, title = null, dockStyle = 'fill', ratio = 0.5) { grok_DockElement(e, title, dockStyle, ratio); }

    tableByName(s) { return new DataFrame(grok_TableByName(s)); }

    /**
     * Creates a generic dataset with the defined number of rows and columns.
     * [dataset] allowed values:
     * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
     * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
     * "biosensor": wearable sensor data: time, x, y, z, temp, eda
     * "random walk": random walk data for the specified number of dimensions
     * **/
    testData(dataset, rows = 10000, columns = 3) { return new DataFrame(grok_TestData(dataset, rows, columns)); }
    getDemoTable(path) {
        return new Promise((resolve, reject) => grok_GetDemoTable(path, (t) => resolve(_wrap(t))));
    }

    newView(name = 'view', children = []) {
        let view = View.create();
        view.name = name;
        ui.appendAll(view.root, children);
        this.addView(view);
        return view;
    }

    getVar(name) {return _toJs(grok_GetVar(name)); };
    setVar(name, value) {return grok_SetVar(name, _toDart(value)); };

    addView(v) { grok_AddView(v.d); }
    addTableView(t, dockStyle = 'fill', width = null) { return new TableView(grok_AddTableView(t.d, dockStyle, width)); }
    getTableView(name) { return new TableView(grok_GetTableView(name)); }
    closeAll() { grok_CloseAll(); }
    route(url) { return View.fromDart(grok_Route(url)); }

    get topMenu() { return new Menu(grok_Get_TopMenu()); }

    parseCsv(csv) { return new DataFrame(grok_ParseCsv(csv)); }

    script(s) {
        return new Promise((resolve, reject) => grok_Script(s, (t) => resolve(_wrap(t, false))));
    }

    loadDataFrame(csvUrl) {
        return new Promise((resolve, reject) => grok_LoadDataFrame(csvUrl, (t) => resolve(_wrap(t, false))));
    }

    linkTables(t1, t2, keyColumns1, keyColumns2, linkTypes) { grok_LinkTables(t1.d, t2.d, keyColumns1, keyColumns2, linkTypes); };
    compareTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2) { grok_CompareTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2); };
    joinTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace) {
        return new DataFrame(grok_JoinTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
    };

    scriptSync(s) { return _wrap(grok_ScriptSync(s), false); }

    openTable(id) { return new Promise((resolve, reject) => grok_OpenTable(id, (t) => resolve(new DataFrame(t)))); }
    query(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t) => resolve(new DataFrame(t))));
    }
    callQuery(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c) => resolve(new FuncCall(c))));
    }
    callFunc(name, parameters = {}, showProgress = false) { return new Promise((resolve, reject) => grok_CallFunc(name, parameters, (out) => resolve(out), showProgress)); }
    evalFunc(name) { return new Promise((resolve, reject) => grok_EvalFunc(name, (out) => resolve(out))); }

    registerFunc(func) { grok_RegisterFunc(func); }

    registerViewer(name, description, createViewer) { grok_RegisterViewer(name, description, createViewer); }

    get currentObject() { return _wrap(grok_Get_CurrentObject(), false); }
    set currentObject(x) { grok_Set_CurrentObject(_toDart(x)); }

    detectSemanticTypes(t) { return new Promise((resolve, reject) => grok_DetectSematicTypes(t.d, (_) => resolve())); }

    onEvent(event, f) { return grok_OnEvent(event, f); };
    onCurrentViewChanged(f) { return grok_OnEvent('d4-current-view-changed', f); }
    onCurrentCellChanged(f) { return grok_OnEvent('d4-current-cell-changed', f); }
    onTableAdded(f) { return grok_OnEvent('d4-table-added', f); }
    onTableRemoved(f) { return grok_OnEvent('d4-table-removed', f); }
    onQueryStarted(f) { return grok_OnEvent('d4-query-started', f); }
    onQueryFinished(f) { return grok_OnEvent('d4-query-finished', f); }

    onViewChanged(f) { return grok_OnEvent('grok-view-changed', (v) => f(View.fromDart(v))); }
    onViewAdded(f) { return grok_OnEvent('grok-view-added', (v) => f(View.fromDart(v))); }
    onViewRemoved(f) { return grok_OnEvent('grok-view-removed', (v) => f(View.fromDart(v))); }
    onViewRenamed(f) { return grok_OnEvent('grok-view-renamed', (v) => f(View.fromDart(v))); }

    onCurrentProjectChanged(f) { return grok_OnEvent('grok-current-project-changed', (a) => f(new Project(a))); }
    onProjectUploaded(f) { return grok_OnEvent('grok-project-uploaded', (a) => f(new Project(a))); }
    onProjectSaved(f) { return grok_OnEvent('grok-project-saved', (a) => f(new Project(a))); }
    onProjectOpened(f) { return grok_OnEvent('grok-project-opened', (a) => f(new Project(a))); }
    onProjectClosed(f) { return grok_OnEvent('grok-project-closed', (a) => f(new Project(a))); }
    onProjectModified(f) { return grok_OnEvent('grok-project-modified', (a) => f(new Project(a))); }
}

class Utils {
    static *range(length) {
        for (let i = 0; i < length; i++)
            yield i;
    }

    /** Returns an 'identity' array where the element in idx-th position is equals to idx. */
    static identity(length) {
        let res = new Array(length);
        for (let i = 0; i < length; i++)
            res[i] = i;
        return res;
    }
}


