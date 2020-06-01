import {DataFrame} from "./dataframe";
import {TableView, View} from "./view";
import {Project, User} from "./entities";
import {toDart, toJs} from "./wrappers";
import {Menu, TabControl} from "./ui_classes";
import {DockManager} from "./docking";

/** Grok entry point, use it to get access to top-level views, tables, methods, etc. */
export class Shell {

    constructor() {
        this.windows = new Windows();
        this.settings = new Settings();
    }

    /** Current table, or null.
     *  @returns {DataFrame} */
    get t() { return new DataFrame(grok_CurrentTable()); }

    /** Current view
     *  @returns {View} */
    get v() { return View.fromDart(grok_Get_CurrentView()); }
    set v(view) { grok_Set_CurrentView(view.d); }

    /** Current project
     *  @returns {Project} */
    get project() { return new Project(grok_Project()); }

    /** List of table names that are currently open
     *  @returns {string[]} */
    get tableNames() { return grok_TableNames(); }

    /** List of currently open tables table
     *  @returns {DataFrame[]}*/
    get tables() { return this.tableNames.map(this.tableByName); }

    /** Current user
     *  @returns {User} */
    get user() { return new User(grok_User()); }

    /** Current object (rendered in the property panel) */
    get o() { return toJs(grok_Get_CurrentObject(), false); }
    set o(x) { grok_Set_CurrentObject(toJs(x)); }

    /** Sidebar.
     *  @returns {TabControl} */
    get sidebar() { return new TabControl(grok_Get_Sidebar()); }

    /** Top menu (not visible by default)
     *  @returns {Menu} */
    get topMenu() { return new Menu(grok_Get_TopMenu()); }

    /** Shows information message (green background)
     * @param {string} s - message */
    info(s) { grok_Balloon(s, 'info'); }

    /** Shows information message (red background)
     * @param {string} s - message */
    error(s) { grok_Balloon(s, 'error'); }

    dockElement(e, title = null, dockStyle = 'fill', ratio = 0.5) {
        grok_DockElement(e, title, dockStyle, ratio);
    }

    /** Opens the view that handles the specified url.
     *  @param {string} url
     *  @returns {View} */
    route(url) { return View.fromDart(grok_Route(url)); }

    /**
     * Adds a view.
     * @param {View} v
     * @returns {View}
     */
    addView(v) {
        grok_AddView(v.d);
        return v;
    }

    /**
     * Adds a new view with the specified name.
     * @param {string } name - view name
     * @param {object[]} children - content to be added by calling {@link ui.appendAll}
     * @returns {View}
     */
    newView(name = 'view', children = []) {
        let view = View.create();
        view.name = name;
        ui.appendAll(view.root, children);
        this.addView(view);
        return view;
    }

    addTableView(t, dockStyle = 'fill', width = null) { return new TableView(grok_AddTableView(t.d, dockStyle, width)); }
    getTableView(name) { return new TableView(grok_GetTableView(name)); }
    closeAll() { grok_CloseAll(); }
    registerViewer(name, description, createViewer) { grok_RegisterViewer(name, description, createViewer); }
    tableByName(s) { return new DataFrame(grok_TableByName(s)); }
    getVar(name) {return _toJs(grok_GetVar(name)); }
    setVar(name, value) {return grok_SetVar(name, toDart(value)); }

    /** @returns {DockManager} */
    get dockManager() {
        return new DockManager(grok_Get_DockManager());
    }
}


/** Controls tool windows */
export class Windows {
    get showSidebar() { return grok_Windows_Get_ShowSidebar(); }
    set showSidebar(x) { return grok_Windows_Set_ShowSidebar(x); }

    get showHelp() { return grok_Windows_Get_ShowSidebar(); }
    set showHelp(x) { return grok_Windows_Set_ShowSidebar(x); }

    get showProperties() { return grok_Windows_Get_ShowProperties(); }
    set showProperties(x) { return grok_Windows_Set_ShowProperties(x); }

    get showVariables() { return grok_Windows_Get_ShowVariables(); }
    set showVariables(x) { return grok_Windows_Set_ShowVariables(x); }

    get showTables() { return grok_Windows_Get_ShowTables(); }
    set showTables(x) { return grok_Windows_Set_ShowTables(x); }

    get showColumns() { return grok_Windows_Get_ShowColumns(); }
    set showColumns(x) { return grok_Windows_Set_ShowColumns(x); }
}


export class Settings {
    /** Hide dock tabs in presentation mode **/
    get hideTabsInPresentationMode() { return grok_Get_HideTabsInPresentationMode(); }
    set hideTabsInPresentationMode(x) { return grok_Set_HideTabsInPresentationMode(x); }

    /** Presentaion mode **/
    get presentationMode() { return grok_Get_PresentationMode(); }
    set presentationMode(x) { return grok_Set_PresentationMode(x); }
}