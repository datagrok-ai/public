import { DataFrame } from './dataframe';
import { TableView, View } from './view';
import { Project } from './entities';
import { toDart, toJs } from './wrappers';
import { Menu, TabControl } from './widgets';
import { DockManager } from './docking';
import { DOCK_TYPE } from './const';
let api = window;
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
  get t() {
    return new DataFrame(api.grok_CurrentTable());
  }
  /** Current view
     *  @type {View} */
  get v() {
    return View.fromDart(api.grok_Get_CurrentView());
  }
  set v(view) {
    // @ts-ignore
    api.grok_Set_CurrentView(view.d);
  }
  /** Current project
     *  @type {Project} */
  get project() {
    return new Project(api.grok_Project());
  }
  /** Names of the currently open tables
     *  @type {string[]} */
  get tableNames() {
    return api.grok_TableNames();
  }
  /** List of currently open tables
     *  @type {DataFrame[]}*/
  get tables() {
    return this.tableNames.map(this.tableByName);
  }
  /** Adds a table to the workspace.
     * @param {DataFrame} table*/
  addTable(table) {
    api.grok_AddTable(table.d);
  }
  /** Closes a table and removes from the workspace.
     * @param {DataFrame} table */
  closeTable(table) {
    api.grok_CloseTable(table.d);
  }
  /** Current user
     *  @type {User} */
  get user() {
    return toJs(api.grok_User());
  }
  /** Current object (rendered in the property panel) */
  get o() {
    return toJs(api.grok_Get_CurrentObject(), false);
  }
  set o(x) {
    api.grok_Set_CurrentObject(toDart(x));
  }
  /** @type {TabControl} */
  get sidebar() {
    return new TabControl(api.grok_Get_Sidebar());
  }
  /** @type {Menu} */
  get topMenu() {
    return new Menu(api.grok_Get_TopMenu());
  }
  /** @type {HTMLDivElement} */
  get topPanel() {
    return api.grok_Get_TopPanel();
  }
  /** @type {HTMLDivElement} */
  get bottomPanel() {
    return api.grok_Get_BottomPanel();
  }
  /** Shows information message (green background)
     * @param {string} x - message */
  info(x) {
    api.grok_Balloon(typeof x == 'string' ? x : JSON.stringify(x), 'info');
  }
  /** Shows information message (red background)
     * @param {string} s - message */
  error(s) {
    api.grok_Balloon(s, 'error');
  }
  /** Shows warning message (yellow background)
     * @param {string} s - message */
  warning(s) {
    api.grok_Balloon(s, 'warning');
  }
  /** Docks element in a separate window.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {HTMLElement} element
     * @param {string} title
     * @param {DockType} dockStyle
     * @param {number} ratio - area to take (relative to parent) */
  dockElement(element, title = null, dockStyle = DOCK_TYPE.FILL, ratio = 0.5) {
    api.grok_DockElement(element, title, dockStyle, ratio);
  }
  /** Opens the view that handles the specified url.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {string} url
     * @returns {View} */
  route(url) {
    return View.fromDart(api.grok_Route(url));
  }
  /** Adds an item to the workspace.
     * It could be a DataFrame, a View, or an HtmlElement.
     * Throws an error, if the item has a different type.
     * @param {DataFrame | View | HTMLElement} item
     * @returns {Shell} */
  add(item) {
    if (item instanceof DataFrame)
      this.addTableView(item);
    else if (item instanceof View)
      this.addView(item);
    else if (item instanceof HTMLElement)
      this.dockElement(item);
    else
      throw 'Unknown type';
    return this;
  }
  /**
     * Adds a view.
     * @param {View} v
     * @param {DockType} dockType
     * @param {number} width
     * @returns {View}
     */
  addView(v, dockType = DOCK_TYPE.FILL, width = null) {
    // @ts-ignore
    if (window.api.grok_AddView == null) {
      document.body.append(v.root);
      v.root.style.width = '100vw';
      v.root.style.height = '100vh';
    }
    else
    // @ts-ignore
      api.grok_AddView(v.d, dockType, width);
    return v;
  }
  /**
     * Adds a new view with the specified name.
     * @param {string} name - view name
     * @param {Object[]} children - content to be added by calling {@link ui.appendAll}
     * @param {string | ElementOptions | null} options
     * @returns {View}
     */
  newView(name = 'view', children, options) {
    if (children == null)
      children = [];
    let view = View.create(options);
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
  addTableView(table, dockType = DOCK_TYPE.FILL, width = null) {
    return toJs(api.grok_AddTableView(table.d, dockType, width));
  }
  /**
     * Returns {@link TableView} for the specified table if it exists, opens a new view if necessary.
     * Search is case-insensitive.
     * @param {string} tableName
     * @returns {TableView} */
  getTableView(tableName) {
    return new TableView(api.grok_GetTableView(tableName));
  }
  /** Closes everything (views, tables, projects) and returns the platform to the initial state. */
  closeAll() {
    api.grok_CloseAll();
  }
  /** Registers a view.
     * @param {string} viewTypeName
     * @param {Function} createView - a function that returns {@link ViewBase}
     * @param {string} viewUrlPath */
  registerView(viewTypeName, createView, viewUrlPath = '') {
    api.grok_RegisterView(viewTypeName, createView, viewUrlPath);
  }
  /** Registers a viewer.
     * Sample: {@link https://public.datagrok.ai/js/samples/scripts/functions/custom-viewers}
     * @param {string} viewerTypeName
     * @param {string} description
     * @param {Function} createViewer - a function that returns {@link JsViewer} */
  registerViewer(viewerTypeName, description, createViewer) {
    api.grok_RegisterViewer(viewerTypeName, description, createViewer);
  }
  /** Returns a table by its name. Search is case-insensitive.
     * @param {string} tableName
     * @returns {DataFrame} */
  tableByName(tableName) {
    return toJs(api.grok_TableByName(tableName));
  }
  /** Returns the value of the specified variable. Search is case-insensitive.
     * @param {string} variableName
     * @returns {Object} */
  getVar(variableName) {
    return toJs(api.grok_GetVar(variableName));
  }
  /** Sets the value of the specified variable, and returns this value.
     * @param {string} variableName
     * @param {Object} variableValue */
  setVar(variableName, variableValue) {
    api.grok_SetVar(variableName, toDart(variableValue));
    return variableValue;
  }
  /** @returns {DockManager} */
  get dockManager() {
    // @ts-ignore
    return new DockManager(api.grok_Get_DockManager());
  }
  /** Clears dirty flag in scratchpad and open projects. */
  clearDirtyFlag() {
    api.grok_ClearDirtyFlag();
  }
}
/** Controls tool windows */
export class Windows {
  /** Controls the visibility of the sidebar.
     * @type {boolean} */
  get showSidebar() {
    return api.grok_Windows_Get_ShowSidebar();
  }
  set showSidebar(x) {
    api.grok_Windows_Set_ShowSidebar(x);
  }
  /** Controls the visibility of the toolbox.
     * @type {boolean} */
  get showToolbox() {
    return api.grok_Windows_Get_ShowToolbox();
  }
  set showToolbox(x) {
    api.grok_Windows_Set_ShowToolbox(x);
  }
  /** Controls the visibility of the console.
     * @type {boolean} */
  get showConsole() {
    return api.grok_Windows_Get_ShowConsole();
  }
  set showConsole(x) {
    api.grok_Windows_Set_ShowConsole(x);
  }
  /** Controls the visibility of the help window.
     * @type {boolean} */
  get showHelp() {
    return api.grok_Windows_Get_ShowHelp();
  }
  set showHelp(x) {
    api.grok_Windows_Set_ShowHelp(x);
  }
  /** Controls the visibility of the properties window.
     * @type {boolean} */
  get showProperties() {
    return api.grok_Windows_Get_ShowProperties();
  }
  set showProperties(x) {
    api.grok_Windows_Set_ShowProperties(x);
  }
  /** Controls the visibility of the variables window.
     * @type {boolean} */
  get showVariables() {
    return api.grok_Windows_Get_ShowVariables();
  }
  set showVariables(x) {
    api.grok_Windows_Set_ShowVariables(x);
  }
  /** Controls the visibility of the tables window.
     * @type {boolean} */
  get showTables() {
    return api.grok_Windows_Get_ShowTables();
  }
  set showTables(x) {
    api.grok_Windows_Set_ShowTables(x);
  }
  /** Controls the visibility of the columns window.
     * @type {boolean} */
  get showColumns() {
    return api.grok_Windows_Get_ShowColumns();
  }
  set showColumns(x) {
    api.grok_Windows_Set_ShowColumns(x);
  }
}
/** User-specific platform settings. */
export class Settings {
  constructor() {
    return new Proxy({}, {
      get: function (target, prop) {
        return toJs(api.grok_PropMixin_GetPropertyValue(api.grok_Get_Settings(), prop));
      },
      set: function (target, prop, value) {
        if (target.hasOwnProperty(prop))
          return target[prop];
        if (target.hasOwnProperty(prop)) {
          target[prop] = value;
          return true;
        }
        api.grok_PropMixin_SetPropertyValue(api.grok_Get_Settings(), prop, toDart(value));
        return true;
      }
    });
  }
  /** Jupyter Notebook URL */
  get jupyterNotebook() {
    return api.grok_Settings_Get_JupyterNotebook();
  }
  /** Jupyter Notebook Token */
  get jupyterNotebookToken() {
    return api.grok_Settings_Get_JupyterNotebookToken();
  }
  /** Hide dock tabs in presentation mode **/
  get hideTabsInPresentationMode() {
    return api.grok_Get_HideTabsInPresentationMode();
  }
  set hideTabsInPresentationMode(x) {
    api.grok_Set_HideTabsInPresentationMode(x);
  }
  /** Presentation mode **/
  get presentationMode() {
    return api.grok_Get_PresentationMode();
  }
  set presentationMode(x) {
    api.grok_Set_PresentationMode(x);
  }
}
