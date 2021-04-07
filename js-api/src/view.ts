import {VIEW_TYPE, VIEWER, ViewerType} from './const';
import {DataFrame} from './dataframe.js';
import * as ui from '../ui';
import {Viewer} from './viewer';
import {DockNode, DockManager} from './docking';
import {Grid} from './grid';
import {Menu, ToolboxPage} from './widgets';
import {Entity, Script} from './entities';
import {toJs} from './wrappers';
import {_options, _toIterable} from './utils';
import { StreamSubscription } from './events';
import $ from "cash-dom";


let api = <any>window;

/**
 * Subclass ViewBase to implement a Datagrok view in JavaScript.
 * */
export class ViewBase {
  d: any;
  subs: StreamSubscription[];
  protected _root: HTMLDivElement;
  private _closing: boolean;

  /** 
   * @constructs ViewBase
   * @param {Object} params - URL parameters.
   * @param {string} path - URL path.
   * @param {boolean} createHost - Create JS host wrapper. */
  constructor(params: object | null = null, path: string = '', createHost: boolean = true) {
    if (createHost)
      this.d = api.grok_View_CreateJsViewHost(this);

    this._root = ui.div([], 'grok-view');
    this._root.tabIndex = 0;

    /** @type {StreamSubscription[]} */
    this.subs = [];  // stream subscriptions - will be canceled when the view is detached

    this._closing = false;
  }

  /** @type {HTMLElement} */
  get root(): HTMLElement {
    return this._root;
  }

  /** View type
   * @type {string} */
  get type(): string {
    return 'js-view-base';
  }

  /** @returns {string|null} View help URL. */
  get helpUrl(): string | null {
    return null;
  }

  /** View name. It gets shown in the tab handle.
   * @type {string} */
  get name(): string {
    // @ts-ignore
    return api.grok_View_Get_Name == null ? this._name : api.grok_View_Get_Name(this.d);
  }

  set name(s: string) {
    if (api.grok_View_Set_Name == null)
    // @ts-ignore
      this._name = s;
    else
      api.grok_View_Set_Name(this.d, s);
  }

  /** @type {string} */
  get description(): string {
    return '';
  }

  set description(s: string) {
  }

  /** @type {Object} */
  get entity(): object | null {
    return null;
  }

  set entity(e: object | null) {
  }

  /** View toolbox.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/toolbox}
   * @type {HTMLElement} */
  get toolbox(): HTMLElement {
    return api.grok_View_Get_Toolbox(this.d);
  }

  set toolbox(x: HTMLElement) {
    api.grok_View_Set_Toolbox(this.d, x);
  }

  /** View menu.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon}
   *  @type {Menu} */
  get ribbonMenu(): Menu {
    return new Menu(api.grok_View_Get_RibbonMenu(this.d));
  }

  set ribbonMenu(menu: Menu) {
    api.grok_View_Set_RibbonMenu(this.d, menu.d);
  }

  get closing(): boolean {
    return this._closing;
  }

  set closing(c: boolean) {
    this._closing = c;
  }

  /** Sets custom view panels on the ribbon.
   * @param {Array<Array<HTMLElement>>} panels
   * @param {boolean} clear Clear all previous before setup
   * Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon} */
  setRibbonPanels(panels: HTMLElement[][], clear: boolean = false): void {
    api.grok_View_SetRibbonPanels(this.d, panels, clear);
  }

  /** @returns {HTMLElement} View icon. */
  getIcon(): HTMLElement | null {
    return null;
  }

  /** @returns {Object} Viewer state map. */
  saveStateMap(): object | null {
    return null;
  }

  /** Load view state map.
   * @param {Object} stateMap - State map. */
  loadStateMap(stateMap: object): void {
  }

  /** 
   * View URI, relative to the platform root. See also {@link basePath}
   * @type {string} */
  get path(): string {
    return '';
  }

  set path(s: string) {
  }

  /** Handles URL path.
   * @param  {string} path - URL path. */
  handlePath(path: string): void {
  }

  /** Checks if URL path is acceptable.
   * @returns {boolean} "true" if path is acceptable, "false" otherwise.
   * @param {string} path - URL path. */
  acceptsPath(path: string): boolean {
    return false;
  }

  /** 
   * Appends an item to this view. Use {@link appendAll} for appending multiple elements.
   * @param {Object} item */
  append(item: HTMLElement): HTMLElement {
    return this.appendAll([ui.render(item)]);
  }

  /**
   * Appends multiple elements this view. Use {@link append} for appending a single element.
   * @param {object[]} items */
  appendAll(items: HTMLElement[]): HTMLElement {
    return ui.appendAll(this.root, items.map(ui.render));
  }

  /** Detaches this view. */
  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Closes this view. */
  close(): void {
    this._closing = true;
    api.grok_View_Close(this.d);
  }
}


/**
 * A view is typically docked in the main document area of the Grok platform.
 * See [TableView], [SketchView], etc
 */
export class View extends ViewBase {

  /** @constructs View */
  constructor(d: any) {
    super(null, '', false);
    this.d = d;
  }

  static fromDart(d: any): View | TableView {
    let type = api.grok_View_Get_Type(d);
    if (type === VIEW_TYPE.TABLE_VIEW)
      return new TableView(d);
    else
      return new View(d);
  }

  /** Creates a new empty view.
   * @param {string | ElementOptions | null} options
   * @returns {View} */
  static create(options?: string | {} | null): View {
    let v = api.grok_View == null ? new View(null) : new View(api.grok_View());
    _options(v.root, 'ui-panel');
    _options(v.root, options);
    return v;
  }

  get box(): boolean {
    return $(this.root).hasClass('ui-box');
  }

  set box(b: boolean) {
    let r = $(this.root);
    r.removeClass('ui-panel').removeClass('ui-box');
    r.addClass(b ? 'ui-box' : 'ui-panel');
  }

  get root(): HTMLElement {
    if (api.grok_View_Get_Root == null)
      return this._root;
    return api.grok_View_Get_Root(this.d);
  }

  get type(): string {
    return api.grok_View_Get_Type(this.d);
  }

  get path(): string {
    return api.grok_View_Get_Path(this.d);
  }

  set path(s: string) {
    api.grok_View_Set_Path(this.d, s);
  }

  /**
   *  View type URI. Note that {@link path} is specific to the instance of the view.
   *  @type {string} */
  get basePath(): string {
    return api.grok_View_Get_BasePath(this.d);
  }

  set basePath(s: string) {
    api.grok_View_Set_BasePath(this.d, s);
  }

  get description(): string {
    return api.grok_View_Get_Description(this.d);
  }

  set description(s: string) {
    api.grok_View_Set_Description(this.d, s);
  }

  /** 
   * Loads previously saved view layout. Only applicable to certain views, such as {@link TableView}.
   *  See also {@link saveLayout}
   *  @param {ViewLayout} layout */
  loadLayout(layout: ViewLayout): void {
    return api.grok_View_Load_Layout(this.d, layout.d);
  }

  /** 
   *  Saves view layout as a string. Only applicable to certain views, such as {@link TableView}.
   *  See also {@link loadLayout}
   *  @returns {ViewLayout} */
  saveLayout(): ViewLayout {
    return new ViewLayout(api.grok_View_Save_Layout(this.d));
  }

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
  constructor(d: any) {
    super(d);
  }

  /** Creates a new table view.
   * @param {DataFrame} table
   * @returns {TableView} */
  static create(table: DataFrame): TableView {
    return new TableView(api.grok_TableView(table.d));
  }

  /** Associated table, if it exists (for TableView), or null.
   *  @type {DataFrame} */
  get table(): DataFrame | null {
    return toJs(api.grok_View_Get_Table(this.d));
  }

  /** @type {Grid} */
  get grid(): Grid {
    return new Grid(api.grok_View_Get_Grid(this.d));
  }

  /** @type {DataFrame} */
  get dataFrame(): DataFrame {
    return toJs(api.grok_View_Get_DataFrame(this.d));
  }

  set dataFrame(x: DataFrame) {
    api.grok_View_Set_DataFrame(this.d, x.d);
  }

  /** View toolbox that gets shown on the left, in the sidebar.
   *  @type {ToolboxPage} */
  get toolboxPage(): ToolboxPage {
    return new ToolboxPage(api.grok_View_Get_ToolboxPage(this.d));
  }

  /** Adds a viewer of the specified type.
   * @param {string | Viewer} v
   * @param options
   * @returns {Viewer} */
  addViewer(v: ViewerType | Viewer, options: object | null = null): Viewer {
    if (typeof v === 'string')
      v = new Viewer(api.grok_View_AddViewerByName(this.d, v));
    else
      api.grok_View_AddViewer(this.d, v.d);
    if (options !== null)
      v.setOptions(options);
    return v;
  }

  /** A dock node for this view.
   *  Use `grok.shell.dockManager` to manipulate it; {@link dockManager} is for controlling
   *  windows that reside inside this view.
   *  @type {DockNode} */
  get dockNode(): DockNode {
    return new DockNode(api.grok_View_Get_DockNode(this.d));
  }

  /** 
   * View's dock manager. Only defined for DockView descendants such as {@link TableView}, UsersView, etc.
   * @type {DockManager} */
  get dockManager(): DockManager {
    return new DockManager(api.grok_View_Get_DockManager(this.d));
  }

  /** 
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/histogram | histogram}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/histogram}
   *  @param options
   *  @returns {Viewer} */
  histogram(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.HISTOGRAM, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/bar-chart | bar chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart}
   *  @param options
   *  @returns {Viewer} */
  barChart(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.BAR_CHART, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/box-plot | box plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/box-plot}
   *  @param options
   *  @returns {Viewer} */
  boxPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.BOX_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/calendar | calendar}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/calendar}
   *  @param options
   *  @returns {Viewer} */
  calendar(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.CALENDAR, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/correlation-plot | correlation plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/corr-plot}
   *  @param options
   *  @returns {Viewer} */
  corrPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.CORR_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/density-plot | density plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/density-plot}
   *  @param options
   *  @returns {Viewer} */
  densityPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.DENSITY_PLOT, options);
  }

  /**
   *  Adds {@link https://datagrok.ai/help/visualize/viewers/filters | filters}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/filters}
   *  @param options
   *  @returns {Viewer} */
  filters(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.FILTERS, options);
  }

  /**
   *  Adds default {@link https://datagrok.ai/help/visualize/viewers/form | form}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/form}
   *  @param options
   *  @returns {Viewer} */
  form(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.FORM, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/google-map | geographical map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/google-map}
   *  @param options
   *  @returns {Viewer} */
  googleMap(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.GOOGLE_MAP, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/heat-map | heat map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/heat-map}
   *  @param options
   *  @returns {Viewer} */
  heatMap(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.HEAT_MAP, options);
  }

  /** 
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/line-chart | line chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart}
   *  @param options
   *  @returns {Viewer} */
  lineChart(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.LINE_CHART, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/shape-map | shape map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/shape-map}
   *  @param options
   *  @returns {Viewer} */
  shapeMap(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.SHAPE_MAP, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/markup | markup viewer}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/markup}
   *  @param options
   *  @returns {Viewer} */
  markup(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.MARKUP, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/matrix-plot | matrix plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot}
   *  @param options
   *  @returns {Viewer} */
  matrixPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.MATRIX_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/network-diagram | network diagram}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/network-diagram}
   *  @param options
   *  @returns {Viewer} */
  networkDiagram(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.NETWORK_DIAGRAM, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pc-plot | parallel coordinates plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot}
   *  @param options
   *  @returns {Viewer} */
  pcPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.PC_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pie-chart | pie chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pie-chart}
   *  @param options
   *  @returns {Viewer} */
  pieChart(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.PIE_CHART, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/scatter-plot | scatter plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot}
   *  @param options
   *  @returns {Viewer} */
  scatterPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.SCATTER_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/3d-scatter-plot | 3D scatter plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot-3d}
   *  @param options
   *  @returns {Viewer} */
  scatterPlot3d(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.SCATTER_PLOT_3D, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/statistics | statistics}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/statistics}
   *  @param options
   *  @returns {Viewer} */
  statistics(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.STATISTICS, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tile-viewer | tile viewer}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer}
   *  @param options
   *  @returns {Viewer} */
  tileViewer(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.TILE_VIEWER, options);
  }

  /** 
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tree-map | tree map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tree-map}
   *  @param options
   *  @returns {Viewer} */
  treeMap(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.TREE_MAP, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/trellis-plot | trellis plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot}
   *  @param options
   *  @returns {Viewer} */
  trellisPlot(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.TRELLIS_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/word-cloud | word cloud}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/word-cloud}
   *  @param options
   *  @returns {Viewer} */
  wordCloud(options: object | null = null): Viewer {
    return this.addViewer(VIEWER.WORD_CLOUD, options);
  }

  /** Resets view layout, leaving only grid visible. */
  resetLayout(): void {
    api.grok_View_ResetLayout(this.d);
  }

  /** Detaches and closes the view. */
  detach(): void {
    api.grok_View_Detach(this.d);
  }

  /** Detaches all viewers. */
  detachViewers(): void {
    api.grok_View_DetachViewers(this.d);
  }

  /** Returns all viewers.
   * @type {Iterable.<Viewer>} */
  get viewers(): Iterable<Viewer> {
    return _toIterable(api.grok_View_Get_Viewers(this.d));
  }
}


/** Base view for working with a collection of objects that reside on the server.
 *  Typically, results are filtered by applying AND operation between two
 *  filters: {@link permanentFilter} (which is set programmatically and is not visible)
 *  and {@link searchValue} entered by the user.
 *
 *  More details on the smart search syntax: {@link https://datagrok.ai/help/overview/smart-search}
 *
 * @extends View */
export class DataSourceCardView extends View {

  /** @constructs DataSourceCardView */
  constructor(d: any) {
    super(d);
  }

  /**
   * User-specified {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
   * @type {string} */
  get searchValue(): string {
    return api.grok_DataSourceCardView_Get_SearchValue(this.d);
  }

  set searchValue(s: string) {
    api.grok_DataSourceCardView_Set_SearchValue(this.d, s);
  }

  /** Programmatically defined invisible
   * {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
   *  @type {string} */
  get permanentFilter(): string {
    return api.grok_DataSourceCardView_Get_PermanentFilter(this.d);
  }

  set permanentFilter(s: string) {
    api.grok_DataSourceCardView_Set_PermanentFilter(this.d, s);
  }
}


/** Projects view */
export class ProjectsView extends DataSourceCardView {
  /** @constructs ProjectsView */
  constructor(d: any) {
    super(d);
  }

  static create(params: object): ProjectsView {
    return new ProjectsView(api.grok_ProjectsView(params));
  }
}

/** Script view */
export class ScriptView extends View {
  /** @constructs ScriptView */
  constructor(d: any) {
    super(d);
  }

  static create(script: Script): ScriptView {
    return new ScriptView(api.grok_ScriptView(script.d));
  }
}


export class ViewLayout extends Entity {

  /** @constructs ViewLayout */
  constructor(d: any) {
    super(d);
  }

  static fromJson(json: string): ViewLayout {
    return new ViewLayout(api.grok_ViewLayout_FromJson(json));
  }

  static fromViewState(state: string): ViewLayout {
    return new ViewLayout(api.grok_ViewLayout_FromViewState(state));
  }

  /** Only defined within the context of the OnViewLayoutXXX events */
  get view(): View {
    return api.grok_ViewLayout_Get_View(this.d);
  }

  get viewState(): string {
    return api.grok_ViewLayout_Get_ViewState(this.d);
  }

  set viewState(state: string) {
    api.grok_ViewLayout_Set_ViewState(this.d, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewLayout_Get_UserDataValue(this.d, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewLayout_Set_UserDataValue(this.d, key, value);
  }

  toJson(): string {
    return api.grok_ViewLayout_ToJson(this.d);
  }
}

export class VirtualView {
  d: any;

  constructor(d: any) {
    this.d = d;
  }

  static create(verticalScroll: boolean = true, maxCols: number = 100): VirtualView {
    return new VirtualView(api.grok_VirtualItemView(verticalScroll, maxCols));
  }

  get root(): HTMLElement {
    return api.grok_VirtualItemView_Get_Root(this.d);
  }

  setData(length: number, renderer: any): void {
    api.grok_VirtualItemView_SetData(this.d, length, renderer);
  }
}
