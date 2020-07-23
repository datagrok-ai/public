import {DataFrame} from "./dataframe";
import {TableView, View, ViewBase} from "./view";
import {Project, User} from "./entities";
import {Menu, TabControl} from "./widgets";
import {DOCK_TYPE, DockManager} from "./docking";
import {JsViewer} from "./viewer";

/**
 * Grok visual shell, use it to get access to top-level views, tables, methods, etc.
 * @example
 * grok.shell.addTableView(grok.data.demo.biosensor(1000));
 * */
export class Shell {

    windows: Windows;
    settings: Settings;

    /** Current table, or null.
     *  @type {DataFrame} */
    get t(): DataFrame | null

    /** Current view
     *  @type {View} */
    get v(): View | null

    set v(view)

    /** Current project
     *  @type {Project} */
    get project(): Project

    /** Names of the currently open tables
     *  @type {string[]} */
    get tableNames(): string[]

    /** List of currently open tables
     *  @type {DataFrame[]}*/
    get tables(): DataFrame[]

    /** Current user
     *  @type {User} */
    get user(): User

    /** Current object (rendered in the property panel) */
    get o(): any

    set o(x)

    /** @type {TabControl} */
    get sidebar(): TabControl

    /** @type {Menu} */
    get topMenu(): Menu

    /** @type {HTMLDivElement} */
    get topPanel(): HTMLDivElement

    /** @type {HTMLDivElement} */
    get bottomPanel(): HTMLDivElement

    /** @returns {DockManager} */
    get dockManager(): DockManager

    /** Shows information message (green background)
     * @param {string} s - message */
    info(s: string): void

    /** Shows information message (red background)
     * @param {string} s - message */
    error(s: string): void

    /** Docks element in a separate window.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {HTMLElement} element
     * @param {string} title
     * @param {DockType} dockStyle
     * @param {number} ratio - area to take (relative to parent) */
    dockElement(element: HTMLElement, title?: string | null, dockStyle?: DOCK_TYPE, ratio?: number): void

    /** Opens the view that handles the specified url.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     * @param {string} url
     * @returns {View} */
    route(url: string): View

    /** Adds an item to the workspace.
     * It could be a DataFrame, a View, a Widget, or an HtmlElement.
     * Throws an error
     * @param {DataFrame | View | List<DataFrame>} item
     * @returns {Shell} */
    add(item: DataFrame | View | HTMLElement): this

    /**
     * Adds a view.
     * @param {View} v
     * @returns {View}
     */
    addView(v: View): View

    /**
     * Adds a new view with the specified name.
     * @param {string } name - view name
     * @param {object[]} children - content to be added by calling {@link ui.appendAll}
     * @returns {View}
     */
    newView(name?: string, children?: object[]): View

    /** Adds a view for the specified table.
     * @param {DataFrame} table
     * @param {DockType} dockType
     * @param {number} width
     * @returns {TableView} */
    addTableView(table: DataFrame, dockType: DOCK_TYPE, width?: number): TableView

    /** Returns {@link TableView} for the specified table if it exists, opens a new view if necessary.
     * Search is case-insensitive.
     * @param {string} tableName
     * @returns {TableView} */
    getTableView(tableName: string): TableView

    /** Closes everything (views, tables, projects) and returns the platform to the initial state. */
    closeAll(): void

    /** Registers a view.
     * @param {string} viewTypeName
     * @param {Function} createView - a function that returns {@link ViewBase}
     * @param {string} viewUrlPath */
    registerView(viewTypeName: string, createView: (...params: any[]) => ViewBase, viewUrlPath?: string): void

    /** Registers a viewer.
     * Sample: {@link https://public.datagrok.ai/js/samples/scripts/functions/custom-viewers}
     * @param {string} viewerTypeName
     * @param {string} description
     * @param {Function} createViewer - a function that returns {@link JsViewer} */
    registerViewer(viewerTypeName: string, description: string, createViewer: (...params: any[]) => JsViewer): void

    /** Returns table by its name. Search is case-insensitive.
     * @param {string} tableName
     * @returns {DataFrame} */
    tableByName(tableName: string): DataFrame

    /** Returns the value of the specified variable. Search is case-insensitive.
     * @param {string} variableName
     * @returns {object} */
    getVar(variableName: string): any

    /** Sets the value of the specified variable, and returns this value.
     * @param {string} variableName
     * @param {object} variableValue */
    setVar<T>(variableName: string, variableValue: T): T
}


/** Controls tool windows */
export class Windows {

    /** Controls the visibility of the sidebar.
     * @type {boolean} */
    get showSidebar(): boolean

    set showSidebar(x)

    /** Controls the visibility of the toolbox.
     * @type {boolean} */
    get showToolbox(): boolean

    set showToolbox(x)

    /** Controls the visibility of the console.
     * @type {boolean} */
    get showConsole(): boolean

    set showConsole(x)

    /** Controls the visibility of the help window.
     * @type {boolean} */
    get showHelp(): boolean

    set showHelp(x)

    /** Controls the visibility of the properties window.
     * @type {boolean} */
    get showProperties(): boolean

    set showProperties(x)

    /** Controls the visibility of the variables window.
     * @type {boolean} */
    get showVariables(): boolean

    set showVariables(x)

    /** Controls the visibility of the tables window.
     * @type {boolean} */
    get showTables(): boolean

    set showTables(x)

    /** Controls the visibility of the columns window.
     * @type {boolean} */
    get showColumns(): boolean

    set showColumns(x)
}


/** User-specific platform settings. */
export class Settings {
    /** Jupyter Notebook URL */
    get jupyterNotebook(): string

    /** Jupyter Notebook Token */
    get jupyterNotebookToken(): string

    /** Hide dock tabs in presentation mode **/
    get hideTabsInPresentationMode(): boolean

    set hideTabsInPresentationMode(x)

    /** Presentation mode **/
    get presentationMode(): boolean

    set presentationMode(x)
}