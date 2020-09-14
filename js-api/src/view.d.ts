import {Menu, ToolboxPage, Widget} from "./widgets";
import {DataFrame} from "./dataframe";
import {Grid} from "./grid";
import {Viewer} from "./viewer";
import {DockManager, DockNode} from "./docking";
import {VIEWER} from "./const";
import {Entity, Script} from "./entities";

/**
 * Subclass ViewBase to implement a Datagrok view in JavaScript.
 * */
export class ViewBase {
    /** @constructs ViewBase
     * @param {Object} params - URL parameters.
     * @param {string} path - URL path.
     * @param {boolean} createHost - Create JS host wrapper. */
    constructor(params: Object, path?: string, createHost?: boolean)

    /** @type {HTMLElement} */
    get root(): HTMLElement

    /** View type
     * @type {string} */
    get type(): string

    /** @returns {string|null} View help URL. */
    get helpUrl(): string | null

    /** View name. It gets shown in the tab handle.
     * @type {string} */
    get name(): string

    set name(s)

    /** @type {string} */
    get description(): string

    set description(s)

    /** @type {Object} */
    get entity(): Object

    set entity(e)

    /** View toolbox.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/toolbox}
     * @type {HTMLElement} */
    get toolbox(): HTMLElement

    set toolbox(x)

    /** View menu.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon}
     *  @type {Menu} */
    get ribbonMenu(): Menu

    set ribbonMenu(menu)

    get closing(): boolean
    set closing(c)

    /** View URI, relative to the platform root. See also {@link basePath}
     * @type {string} */
    get path(): string

    set path(s)

    /** Sets custom view panels on the ribbon.
     * @param {Array<Array<HTMLElement>>} panels
     * @param {boolean} clear Clear all previous before setup
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon} */
    setRibbonPanels(panels: HTMLElement[][], clear?: boolean): void

    /** @returns {HTMLElement} View icon. */
    getIcon(): HTMLElement

    /** @returns {Object} Viewer state map. */
    saveStateMap(): object

    /** Load view state map.
     * @param {Object} stateMap - State map. */
    loadStateMap(stateMap: object): void

    /** Handles URL path.
     * @param  {string} path - URL path. */
    handlePath(path: string): void

    /** Checks if URL path is acceptable.
     * @returns {boolean} "true" if path is acceptable, "false" otherwise.
     * @param {string} path - URL path. */
    acceptsPath(path: string): boolean

    /** Appends an item to this view. Use {@link appendAll} for appending multiple elements.
     * @param {HTMLElement | Widget} item */
    append(item: HTMLElement | Widget): void

    /** Appends multiple elements this view. Use {@link appendAll} for appending multiple elements.
     * @param {object[]} items */
    appendAll(items: (HTMLElement | Widget)[]): void

    /** Detach this view. */
    detach(): void

    /** Closes this view. */
    close(): void
}


/**
 * A view is typically docked in the main document area of the Grok platform.
 * See [TableView], [SketchView], etc
 */
export class View extends ViewBase {

    /** @constructs View */
    constructor(d: any)

    get root(): HTMLElement

    get type(): string

    get path(): string

    set path(s)

    /** View type URI. Note that {@path} is specific to the instance of the view.
     *  @type {string} */
    get basePath(): string

    set basePath(s)

    get description(): string

    set description(s)

    static fromDart(d: any): TableView | View

    /** Creates a new empty view.
     * @returns {View} */
    static create(param?: any): View

    /** Loads previously saved view layout. Only applicable to certain views, such as TableView.
     *  See also {@link saveLayout}
     *  @param {ViewLayout} layout */
    loadLayout(layout: ViewLayout): void

    /** Saves view layout as a string. Only applicable to certain views, such as TableView.
     *  See also {@link loadLayout}
     *  @returns {ViewLayout} */
    saveLayout(): ViewLayout

    close(): void
}


/**
 * A {@link View} that is associated with a {@link DataFrame} and exposes
 * exploratory data analysis functionality. This view gets opened whenever
 * a new table is added to the workspace when a user drag-and-drops a CSV file,
 * or opens a table in any other way.
 * @extends View
 */
export class TableView extends View {
    /** @constructs TableView */
    constructor(d: any)

    /** Associated table, if it exists (for TableView), or null.
     *  @type {DataFrame} */
    get table(): DataFrame | null

    /** @type {Grid} */
    get grid(): Grid

    /** @type {DataFrame} */
    get dataFrame(): DataFrame

    set dataFrame(x)

    /** View toolbox that gets shown on the left, in the sidebar
     *  @type {ToolboxPage} */
    get toolboxPage(): ToolboxPage

    /** A dock node for this view.
     *  Use `grok.shell.dockManager` to manipulate it; {@link dockManager} is for controlling
     *  widows that reside inside this view.
     *  @type {DockNode} */
    get dockNode(): DockNode

    /** View's dock manager. Only defined for DockView descendants such as {@link TableView}, UsersView, etc.
     * @type {DockManager} */
    get dockManager(): DockManager | null

    /** Adds a viewer of the specified type.
     * @param {string} viewerType
     * @param options
     * @returns {Viewer} */
    addViewer(viewerType: VIEWER | string, options?: Object): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/histogram | histogram}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/histogram}
     *  @param options
     *  @returns {Viewer} */
    histogram(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/bar-chart | bar chart}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/bar-chart}
     *  @param options
     *  @returns {Viewer} */
    barChart(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/box-plot | box plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/box-plot}
     *  @param options
     *  @returns {Viewer} */
    boxPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/calendar | calendar}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/calendar}
     *  @param options
     *  @returns {Viewer} */
    calendar(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/correlation-plot | correlation plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/corr-plot}
     *  @param options
     *  @returns {Viewer} */
    corrPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/density-plot | density plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/density-plot}
     *  @param options
     *  @returns {Viewer} */
    densityPlot(options?: Object | null): Viewer

    /** Adds {@link https://datagrok.ai/help/visualize/viewers/filters | filters}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/filters}
     *  @param options
     *  @returns {Viewer} */
    filters(options?: Object | null): Viewer

    /** Adds default {@link https://datagrok.ai/help/visualize/viewers/form | form}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/form}
     *  @param options
     *  @returns {Viewer} */
    form(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/google-map | geographical map}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/google-map}
     *  @param options
     *  @returns {Viewer} */
    googleMap(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/heat-map | heat map}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/heat-map}
     *  @param options
     *  @returns {Viewer} */
    heatMap(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/line-chart | line chart}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/line-chart}
     *  @param options
     *  @returns {Viewer} */
    lineChart(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/shape-map | shape map}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/shape-map}
     *  @param options
     *  @returns {Viewer} */
    shapeMap(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/markup | markup viewer}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/markup}
     *  @param options
     *  @returns {Viewer} */
    markup(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/matrix-plot | matrix plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/matrix-plot}
     *  @param options
     *  @returns {Viewer} */
    matrixPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/network-diagram | network diagram}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/network-diagram}
     *  @param options
     *  @returns {Viewer} */
    networkDiagram(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/pc-plot | parallel coordinates plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/pc-plot}
     *  @param options
     *  @returns {Viewer} */
    pcPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/pie-chart | pie chart}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/pie-chart}
     *  @param options
     *  @returns {Viewer} */
    pieChart(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/scatter-plot | scatter plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot}
     *  @param options
     *  @returns {Viewer} */
    scatterPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/3d-scatter-plot | 3D scatter plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot-3d}
     *  @param options
     *  @returns {Viewer} */
    scatterPlot3d(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/statistics | statistics}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/statistics}
     *  @param options
     *  @returns {Viewer} */
    statistics(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/tile-viewer | tile viewer}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/tile-viewer}
     *  @param options
     *  @returns {Viewer} */
    tileViewer(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/tree-map | tree map}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/tree-map}
     *  @param options
     *  @returns {Viewer} */
    treeMap(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/trellis-plot | trellis plot}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/trellis-plot}
     *  @param options
     *  @returns {Viewer} */
    trellisPlot(options?: Object | null): Viewer

    /** Adds a {@link https://datagrok.ai/help/visualize/viewers/word-cloud | word cloud}.
     *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/word-cloud}
     *  @param options
     *  @returns {Viewer} */
    wordCloud(options?: Object | null): Viewer

    /** Resets view layout, leaving only grid visible. */
    resetLayout(): void

    /** Detaches and closes the view. */
    detach(): void

    /** Detaches all viewers. */
    detachViewers(): void
}


/** Base view for working with collection of objects that reside on the server.
 *  Typically, results are filtered by applying AND operation between two
 *  filters: {@link permanentFilter} (which is programmatically and is not visible)
 *  and {@link searchValue} entered by user.
 *
 *  More details on the smart search syntax: @{link https://datagrok.ai/help/overview/smart-search}
 *
 * @extends View */
export class DataSourceCardView extends View {

    /** @constructs DataSourceCardView */
    constructor(d: any)

    /** User-specified {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
     * @type {string} */
    get searchValue(): string

    set searchValue(s)

    /** Programmatically defined invisible
     * {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
     *  @type {string} */
    get permanentFilter(): string

    set permanentFilter(s)
}


/** Projects view */
export class ProjectsView extends DataSourceCardView {
    /** @constructs ProjectsView */
    constructor(d: any)

    static create(params?: Object): ProjectsView
}

/** Script view */
export class ScriptView extends View {
    /** @constructs ScriptView */
    constructor(d: any)

    static create(script: Script): ScriptView
}

export class ViewLayout extends Entity {
    constructor(d: any)

    static fromJson(json: string): ViewLayout

    toJson(): string
}


export class VirtualView {
    constructor(d: any)

    get root(): HTMLElement

    static create(verticalScroll?: boolean, maxCols?: number): VirtualView

    setData(length: number, renderer: any): void
}