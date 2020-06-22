import {DataFrame} from "./dataframe";
import {TableView, View} from "./view";
import {Project, User} from "./entities";
import {toDart, toJs} from "./wrappers";
import {Menu, TabControl} from "./widgets";
import {DockManager} from "./docking";

/**
 * Grok visual shell, use it to get access to top-level views, tables, methods, etc.
 * @example
 * grok.shell.addTableView(grok.data.demo.biosensor(1000));
 * */
export class Shell {

    constructor() {
        /** Tool windows
         * @type Windows */
        this.windows = new Windows();

        /** Settings
         * @type {Settings} */
        this.settings = new Settings();
    }

    /** Current table, or null.
     *  @type {DataFrame} */
    get t() { return new DataFrame(grok_CurrentTable()); }

    /** Current view
     *  @type {View} */
    get v() { return View.fromDart(grok_Get_CurrentView()); }
    set v(view) { grok_Set_CurrentView(view.d); }

    /** Current project
     *  @type {Project} */
    get project() { return new Project(grok_Project()); }

    /** Names of the currently open tables
     *  @type {string[]} */
    get tableNames() { return grok_TableNames(); }

    /** List of currently open tables
     *  @type {DataFrame[]}*/
    get tables() { return this.tableNames.map(this.tableByName); }

    /** Current user
     *  @type {User} */
    get user() { return new User(grok_User()); }

    /** Current object (rendered in the property panel) */
    get o() { return toJs(grok_Get_CurrentObject(), false); }
    set o(x) { grok_Set_CurrentObject(toJs(x)); }

    /** @type {TabControl} */
    get sidebar() { return new TabControl(grok_Get_Sidebar()); }

    /** @type {Menu} */
    get topMenu() { return new Menu(grok_Get_TopMenu()); }

    /** Shows information message (green background)
     * @param {string} s - message */
    info(s) { grok_Balloon(s, 'info'); }

    /** Shows information message (red background)
     * @param {string} s - message */
    error(s) { grok_Balloon(s, 'error'); }

    /** Docks element in a separate window.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {HTMLElement} element
     * @param {string} title
     * @param {DockType} dockStyle
     * @param {number} ratio - area to take (relative to parent) */
    dockElement(element, title = null, dockStyle = DG.DOCK_TYPE.FILL, ratio = 0.5) {
        grok_DockElement(element, title, dockStyle, ratio);
    }

    /** Opens the view that handles the specified url.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {string} url
     * @returns {View} */
    route(url) { return View.fromDart(grok_Route(url)); }

    /** Adds an item to the workspace.
     * It could be a DataFrame, a View, a Widget, or an HtmlElement.
     * Throws an error
     * @param {DataFrame | View | List<DataFrame>} item
     * @returns {Shell} */
    add(item) {
        if (item instanceof DataFrame)
            this.addTableView(item);
        else if (item instanceof View)
            this.addView(item);
        else if (item instanceof HtmlElement)
            this.dockElement(item);
        else
            throw 'Unknown type';
        return this;
    }

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

    /** Adds a view for the specified table.
     * @param {DataFrame} table
     * @param {DockType} dockType
     * @param {number} width
     * @returns {TableView} */
    addTableView(table, dockType = DG.DOCK_TYPE.FILL, width = null) {
        return new TableView(grok_AddTableView(table.d, dockType, width));
    }

    /** Returns {@link TableView} for the specified table if it exists, opens a new view if necessary.
     * Search is case-insensitive.
     * @param {string} tableName
     * @returns {TableView} */
    getTableView(tableName) { return new TableView(grok_GetTableView(tableName)); }

    /** Closes everything (views, tables, projects) and returns the platform to the initial state. */
    closeAll() { grok_CloseAll(); }

    /** Registers a view.
     * @param {string} viewTypeName
     * @param {Function} createView - a function that returns {@link ViewBase}
     * @param {string} viewUrlPath */
    registerView(viewTypeName, createView, viewUrlPath = '') { grok_RegisterView(viewTypeName, createView, viewUrlPath); }

    /** Registers a viewer.
     * Sample: {@link https://public.datagrok.ai/js/samples/scripts/functions/custom-viewers}
     * @param {string} viewerTypeName
     * @param {string} description
     * @param {Function} createViewer - a function that returns {@link JsViewer} */
    registerViewer(viewerTypeName, description, createViewer) { grok_RegisterViewer(viewerTypeName, description, createViewer); }

    /** Returns table by its name. Search is case-insensitive.
     * @param {string} tableName
     * @returns {DataFrame} */
    tableByName(tableName) { return new DataFrame(grok_TableByName(tableName)); }

    /** Returns the value of the specified variable. Search is case-insensitive.
     * @param {string} variableName
     * @returns {object} */
    getVar(variableName) {return _toJs(grok_GetVar(variableName)); }

    /** Sets the value of the specified variable, and returns this value.
     * @param {string} variableName
     * @param {object} variableValue */
    setVar(variableName, variableValue) {
        grok_SetVar(variableName, toDart(variableValue));
        return variableValue;
    }

    /** @returns {DockManager} */
    get dockManager() {
        return new DockManager(grok_Get_DockManager());
    }
}


/** Controls tool windows */
export class Windows {

    /** Controls the visibility of the sidebar.
     * @type {boolean} */
    get showSidebar() { return grok_Windows_Get_ShowSidebar(); }
    set showSidebar(x) { return grok_Windows_Set_ShowSidebar(x); }

    /** Controls the visibility of the toolbox.
     * @type {boolean} */
    get showToolbox() { return grok_Windows_Get_ShowToolbox(); }
    set showToolbox(x) { return grok_Windows_Set_ShowToolbox(x); }

    /** Controls the visibility of the console.
     * @type {boolean} */
    get showConsole() { return grok_Windows_Get_ShowConsole(); }
    set showConsole(x) { return grok_Windows_Set_ShowConsole(x); }

    /** Controls the visibility of the help window.
     * @type {boolean} */
    get showHelp() { return grok_Windows_Get_ShowSidebar(); }
    set showHelp(x) { return grok_Windows_Set_ShowSidebar(x); }

    /** Controls the visibility of the properties window.
     * @type {boolean} */
    get showProperties() { return grok_Windows_Get_ShowProperties(); }
    set showProperties(x) { return grok_Windows_Set_ShowProperties(x); }

    /** Controls the visibility of the variables window.
     * @type {boolean} */
    get showVariables() { return grok_Windows_Get_ShowVariables(); }
    set showVariables(x) { return grok_Windows_Set_ShowVariables(x); }

    /** Controls the visibility of the tables window.
     * @type {boolean} */
    get showTables() { return grok_Windows_Get_ShowTables(); }
    set showTables(x) { return grok_Windows_Set_ShowTables(x); }

    /** Controls the visibility of the columns window.
     * @type {boolean} */
    get showColumns() { return grok_Windows_Get_ShowColumns(); }
    set showColumns(x) { return grok_Windows_Set_ShowColumns(x); }
}


/** User-specific platform settings. */
export class Settings {
    /** Jupyter Notebook URL */
    get jupyterNotebook() { return grok_Settings_Get_JupyterNotebook(); }

    /** Jupyter Notebook Token */
    get jupyterNotebookToken() { return grok_Settings_Get_JupyterNotebookToken(); }

    /** Hide dock tabs in presentation mode **/
    get hideTabsInPresentationMode() { return grok_Get_HideTabsInPresentationMode(); }
    set hideTabsInPresentationMode(x) { return grok_Set_HideTabsInPresentationMode(x); }

    /** Presentation mode **/
    get presentationMode() { return grok_Get_PresentationMode(); }
    set presentationMode(x) { return grok_Set_PresentationMode(x); }
}
