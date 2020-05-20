import {GrokPackage, DataConnection, Project} from "./entities.js";
import {ui, Balloon} from "./ui.js";
import {Dapi} from "./dapi.js";
import {DataFrame} from "./dataframe.js";
import {TableView, View} from "./view";
import {User} from "./entities";
import {TabControl} from "./ui";
import {Functions} from "./functions";
import * as rxjs from "rxjs";
import * as rxjsOperators from "rxjs/operators";
import {Events} from "./events";

//import { fromEventPattern } from 'rxjs';
/** Grok entry point, use it to get access to top-level views, tables, methods, etc. */

export class Xplorer {

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

    /** Current user */
    get user() { return new User(grok_User()); }

    get sidebar() { return new TabControl(grok_Get_Sidebar()); }

    /** Presentaion mode **/
    get presentationMode() { return grok_Get_PresentationMode(); }
    set presentationMode(x) { return grok_Set_PresentationMode(x); }

    /** Hide dock tabs in presentaion mode **/
    get hideTabsInPresentationMode() { return grok_Get_HideTabsInPresentationMode(); }
    set hideTabsInPresentationMode(x) { return grok_Set_HideTabsInPresentationMode(x); }

    get topMenu() { return new Menu(grok_Get_TopMenu()); }
    get currentObject() { return _wrap(grok_Get_CurrentObject(), false); }
    set currentObject(x) { grok_Set_CurrentObject(_toDart(x)); }
}

export class Utils {
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

export * from './chem.js';
export * from './const.js';
export * from './dapi.js';
export * from './dataframe.js';
export * from './entities.js';
export * from './events.js';
export * from './functions.js';
export * from './grid.js';
export * from './ml.js';
export * from './ui.js';
export * from './utils.js';
export * from './view.js';
export * from './viewer.js';

export * from "rxjs";
export * from "rxjs/operators";

//export let grok = new Grok();
/*
window.grok = grok;
window.ui = ui;
window.GrokPackage = GrokPackage;*/

export let functions = new Functions();

export let events = new Events();

export let dapi = new Dapi();

export let xp = new Xplorer();

export function info(s) { Balloon.info(s); }

export function dockElement(e, title = null, dockStyle = 'fill', ratio = 0.5) { grok_DockElement(e, title, dockStyle, ratio); }

export function tableByName(s) { return new DataFrame(grok_TableByName(s)); }

export * from './utils.js';

/**
 * Creates a generic dataset with the defined number of rows and columns.
 * [dataset] allowed values:
 * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
 * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
 * "biosensor": wearable sensor data: time, x, y, z, temp, eda
 * "random walk": random walk data for the specified number of dimensions
 *
 * @returns {DataFrame} */
export function testData(dataset, rows = 10000, columns = 3) { return new DataFrame(grok_TestData(dataset, rows, columns)); }
export function getDemoTable(path) {
    return new Promise((resolve, reject) => grok_GetDemoTable(path, (t) => resolve(_wrap(t))));
}

export function newView(name = 'view', children = []) {
    let view = View.create();
    view.name = name;
    ui.appendAll(view.root, children);
    this.addView(view);
    return view;
}

export function getVar(name) {return _toJs(grok_GetVar(name)); }
export function setVar(name, value) {return grok_SetVar(name, _toDart(value)); }

export function addView(v) { grok_AddView(v.d); }
export function addTableView(t, dockStyle = 'fill', width = null) { return new TableView(grok_AddTableView(t.d, dockStyle, width)); }
export function getTableView(name) { return new TableView(grok_GetTableView(name)); }
export function closeAll() { grok_CloseAll(); }
export function route(url) { return View.fromDart(grok_Route(url)); }

/** @returns {DataFrame} */
export function parseCsv(csv) { return new DataFrame(grok_ParseCsv(csv)); }

export function script(s) {
    return new Promise((resolve, reject) => grok_Script(s, (t) => resolve(_wrap(t, false))));
}

export function loadDataFrame(csvUrl) {
    return new Promise((resolve, reject) => grok_LoadDataFrame(csvUrl, (t) => resolve(_wrap(t, false))));
}

export function linkTables(t1, t2, keyColumns1, keyColumns2, linkTypes) { grok_LinkTables(t1.d, t2.d, keyColumns1, keyColumns2, linkTypes); };
export function compareTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2) { grok_CompareTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2); };
export function joinTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace) {
    return new DataFrame(grok_JoinTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
}

export function scriptSync(s) { return _wrap(grok_ScriptSync(s), false); }

export function openTable(id) { return new Promise((resolve, reject) => grok_OpenTable(id, (t) => resolve(new DataFrame(t)))); }
export function query(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
    return new Promise((resolve, reject) => grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t) => resolve(new DataFrame(t))));
}
export function callQuery(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
    return new Promise((resolve, reject) => grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c) => resolve(new FuncCall(c))));
}
export function callFunc(name, parameters = {}, showProgress = false) { return new Promise((resolve, reject) => grok_CallFunc(name, parameters, (out) => resolve(out), showProgress)); }
export function evalFunc(name) { return new Promise((resolve, reject) => grok_EvalFunc(name, (out) => resolve(out))); }

export function registerFunc(func) { grok_RegisterFunc(func); }

export function registerViewer(name, description, createViewer) { grok_RegisterViewer(name, description, createViewer); }

export function detectSemanticTypes(t) { return new Promise((resolve, reject) => grok_DetectSematicTypes(t.d, (_) => resolve())); }

export function onEvent(event, f) { return grok_OnEvent(event, f); }
export function onCurrentViewChanged(f) { return grok_OnEvent('d4-current-view-changed', f); }
export function onCurrentCellChanged(f) { return grok_OnEvent('d4-current-cell-changed', f); }
export function onTableAdded(f) { return grok_OnEvent('d4-table-added', f); }
export function onTableRemoved(f) { return grok_OnEvent('d4-table-removed', f); }
export function onQueryStarted(f) { return grok_OnEvent('d4-query-started', f); }
export function onQueryFinished(f) { return grok_OnEvent('d4-query-finished', f); }

export function onViewChanged(f) { return grok_OnEvent('grok-view-changed', (v) => f(View.fromDart(v))); }
export function onViewAdded(f) { return grok_OnEvent('grok-view-added', (v) => f(View.fromDart(v))); }
export function onViewRemoved(f) { return grok_OnEvent('grok-view-removed', (v) => f(View.fromDart(v))); }
export function onViewRenamed(f) { return grok_OnEvent('grok-view-renamed', (v) => f(View.fromDart(v))); }

export function onCurrentProjectChanged(f) { return grok_OnEvent('grok-current-project-changed', (a) => f(new Project(a))); }
export function onProjectUploaded(f) { return grok_OnEvent('grok-project-uploaded', (a) => f(new Project(a))); }
export function onProjectSaved(f) { return grok_OnEvent('grok-project-saved', (a) => f(new Project(a))); }
export function onProjectOpened(f) { return grok_OnEvent('grok-project-opened', (a) => f(new Project(a))); }
export function onProjectClosed(f) { return grok_OnEvent('grok-project-closed', (a) => f(new Project(a))); }
export function onProjectModified(f) { return grok_OnEvent('grok-project-modified', (a) => f(new Project(a))); }
