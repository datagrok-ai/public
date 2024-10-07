import {VIEW_TYPE, VIEWER, ViewerType, ViewType} from '../const';
import {DataFrame} from '../dataframe.js';
import * as ui from '../../ui';
import {FilterGroup, ScatterPlotViewer, Viewer} from '../viewer';
import {DockManager, DockNode} from '../docking';
import {Grid} from '../grid';
import {Menu, ToolboxPage, TreeViewGroup} from '../widgets';
import {ColumnInfo, Entity, Script, TableInfo} from '../entities';
import {toDart, toJs} from '../wrappers';
import {_options, _toIterable, MapProxy} from '../utils';
import $ from "cash-dom";
import {Subscription} from "rxjs";
import {FuncCall} from "../functions";
import {IDartApi} from "../api/grok_api.g";
import {
  IBarChartSettings,
  IBoxPlotSettings,
  ICalendarSettings,
  ICorrelationPlotSettings,
  IDensityPlotSettings,
  IFiltersSettings,
  IFormSettings,
  IGridSettings,
  IHistogramSettings,
  ILineChartSettings,
  IMapViewerSettings,
  IMarkupViewerSettings,
  IMatrixPlotSettings,
  INetworkDiagramSettings,
  IPcPlotSettings,
  IPieChartSettings,
  IScatterPlot3dSettings,
  IScatterPlotSettings,
  IStatsViewerSettings, ITileViewerSettings, ITreeMapSettings, ITrellisPlotSettings
} from "../interfaces/d4";

const api: IDartApi = <any>window;

/**
 * Subclass ViewBase to implement a Datagrok view in JavaScript.
 * */
export class ViewBase {
  dart: any;
  subs: Subscription[];
  private _helpUrl: string | null = null;
  protected _root: HTMLElement;
  private _closing: boolean;

  /** 
   * @constructs ViewBase
   * @param {Object} params - URL parameters.
   * @param {string} path - URL path.
   * @param {boolean} createHost - Create JS host wrapper. */
  constructor(params: object | null = null, path: string = '', createHost: boolean = true) {
    if (createHost)
      this.dart = api.grok_View_CreateJsViewHost(this);

    this.name = 'New view';
    this._root = ui.panel([], 'grok-view');
    this._root.tabIndex = 0;

    /** @type {StreamSubscription[]} */
    this.subs = [];  // stream subscriptions - will be canceled when the view is detached

    this._closing = false;
  }

  /** @type {HTMLElement} */
  get root(): HTMLElement {
    return this._root;
  }

  set root(newRoot: HTMLElement){
    this._root = newRoot
  } 

  get box(): boolean {
    return $(this.root).hasClass('ui-box');
  }

  set box(b: boolean) {
    let r = $(this.root);
    r.removeClass('ui-panel').removeClass('ui-box').removeClass('ui-div');
    r.addClass(b ? 'ui-box' : 'ui-panel');
  }

  /** View type
   * @type {string} */
  get type(): string {
    return 'js-view-base';
  }

  /** @returns {string|null} View help URL. */
  get helpUrl(): string | null {
    return this._helpUrl;
  }

  set helpUrl(url: string | null) {
    this._helpUrl = url;
  }

  // /** View name */
  // get name(): string { return api.grok_View_Get_Name(this.dart); }
  // set name(s: string) { api.grok_View_Set_Name(this.dart, s); }

  protected _name: string = 'New View';

  /** @type {string} */
  get name(): string {
    return this._name;
  }

  set name(s: string) {
    this._name = s;
    api.grok_View_Set_Name(this.dart, s);
  }

  get parentCall(): FuncCall | undefined { return toJs(api.grok_View_Get_ParentCall(this.dart)); }
  set parentCall(s: FuncCall | undefined) { api.grok_View_Set_ParentCall(this.dart, toDart(s)); }

  get parentView(): ViewBase { return toJs(api.grok_View_Get_ParentView(this.dart)); }
  set parentView(s: ViewBase) { api.grok_View_Set_ParentView(this.dart, toDart(s)); }

  get description(): string { return ''; }
  set description(s: string) { }

  get entity(): object | null { return null; }
  set entity(_e: object | null) { }

  /** @deprecated use path instead */
  get basePath(): string { return api.grok_View_Get_BasePath(this.dart); }
  set basePath(s: string) { api.grok_View_Set_BasePath(this.dart, s); }

  /** View toolbox.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/toolbox} */
  get toolbox(): HTMLElement { return api.grok_View_Get_Toolbox(this.dart); }
  set toolbox(x: HTMLElement) { api.grok_View_Set_Toolbox(this.dart, x); }

  /** View menu.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon} */
  get ribbonMenu(): Menu { return new Menu(api.grok_View_Get_RibbonMenu(this.dart)); }
  set ribbonMenu(menu: Menu) { api.grok_View_Set_RibbonMenu(this.dart, menu.dart); }

  /** Status bar panels to be shown on the bottom */
  get statusBarPanels(): HTMLDivElement[] { return api.grok_View_Get_StatusBarPanels(this.dart); }
  set statusBarPanels(panels: HTMLDivElement[]) { api.grok_View_Set_StatusBarPanels(this.dart, panels); }

  /** Whether the view is currently closing. */
  get closing(): boolean { return this._closing; }
  set closing(c: boolean) { this._closing = c; }

  /** Sets custom view panels on the ribbon.
   * @param {Array<Array<HTMLElement>>} panels
   * @param {boolean} clear Clear all previous before setup
   * Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon} */
  setRibbonPanels(panels: HTMLElement[][], clear: boolean = true): void {
    api.grok_View_SetRibbonPanels(this.dart, panels, clear);
  }

  getRibbonPanels(): HTMLElement[][] {
    return api.grok_View_GetRibbonPanels(this.dart);
  }

  /** @returns {HTMLElement} View icon. Override in subclasses. */
  getIcon(): HTMLElement | null { return null; }

  /** @returns {Object} Viewer state map. Override in subclasses. */
  saveStateMap(): object | null { return null; }

  /** Loads view state map. Override in subclasses. */
  loadStateMap(_stateMap: object): void { }

  /** View URI, relative to the view root */
  get path(): string { return api.grok_View_Get_Path(this.dart); }
  set path(s: string) { api.grok_View_Set_Path(this.dart, s); }

  /** Handles URL path. Override in subclasses. */
  handlePath(_urlPath: string): void { }

  /** Checks if URL path is acceptable. Override in subclasses.
   * @returns {boolean} "true" if path is acceptable, "false" otherwise. */
  acceptsPath(_urlPath: string): boolean { return false; }

  /** 
   * Appends an item to this view. Use {@link appendAll} for appending multiple elements.
   * @param {Object} item */
  append(item: any): HTMLElement {
    return this.appendAll([ui.render(item)]);
  }

  /**
   * Appends multiple elements this view. Use {@link append} for appending a single element.
   * @param {object[]} items */
  appendAll(items: HTMLElement[]): HTMLElement {
    return ui.appendAll(this.root, items.map(x => ui.render(x)));
  }

  /** Detaches this view. */
  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Closes this view. */
  close(): void {
    this._closing = true;
    api.grok_View_Close(this.dart);
  }
}


/**
 * A view is typically docked in the main document area of the Grok platform.
 * See [TableView], [SketchView], etc
 */
export class View extends ViewBase {

  /** Contains auxiliary information */
  public temp: any;

  /** @constructs View */
  constructor(dart: any) {
    super(null, '', false);
    this.dart = dart;
    this.temp = new MapProxy(api.grok_Widget_Get_Temp(this.dart));
  }

  static fromDart(dart: any): View | TableView {
    let type = api.grok_View_Get_Type(dart);
    if (type === VIEW_TYPE.TABLE_VIEW)
      return new TableView(dart);
    else
      return new View(dart);
  }

  /** Creates a view for the specified object, if it is registered, or null otherwise. */
  static forObject(x: any): View | null {
    return api.grok_View_ForObject(toDart(x));
  }

  static fromRoot(root: HTMLElement) {
    let view = View.create();
    view.root.appendChild(root);
    return view;
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

  /** Creates one of the standard views based on the view type (such as 'functions') */
  static createByType(viewType: ViewType | string, options?: any): View {
    return new View(api.grok_View_CreateByType(viewType, options));
  }

  static fromViewAsync(getViewAsync: () => Promise<View>, ribbon: boolean = true) {
    return toJs(api.grok_View_FromViewAsync(getViewAsync, ribbon));
  }

  get root(): HTMLElement {
    return api.grok_View_Get_Root(this.dart);
  }

  get type(): string {
    return api.grok_View_Get_Type(this.dart);
  }

  get id(): string {
    return api.grok_View_Get_Id(this.dart);
  }

  get description(): string {
    return api.grok_View_Get_Description(this.dart);
  }

  set description(s: string) {
    api.grok_View_Set_Description(this.dart, s);
  }

  /** @returns {string|null} View help URL. */
  get helpUrl(): string | null {
    return api.grok_View_Get_HelpUrl(this.dart);
  }

  set helpUrl(url: string | null) {
    api.grok_View_Set_HelpUrl(this.dart, url);
  }

  /**
   * Loads previously saved view layout. Only applicable to certain views, such as {@link TableView}.
   *  See also {@link saveLayout} */
  loadLayout(layout: ViewLayout, pickupColumnTags?: boolean): void {
    return api.grok_View_Load_Layout(this.dart, layout.dart, pickupColumnTags);
  }

  /** 
   *  Saves view layout as a string. Only applicable to certain views, such as {@link TableView}.
   *  See also {@link loadLayout}
   *  @returns {ViewLayout} */
  saveLayout(): ViewLayout {
    return toJs(api.grok_View_Save_Layout(this.dart));
  }

  /**
   *  Saves view as a ViewInfo. Only applicable to certain views, such as {@link TableView}.
   *  @returns {ViewInfo} */
  getInfo(): ViewLayout {
    return toJs(api.grok_View_Get_Info(this.dart));
  }


  /** View name. It gets shown in the tab handle.
   * @type {string} */
  get name(): string {
    // @ts-ignore
    return api.grok_View_Get_Name == null ? this._name : api.grok_View_Get_Name(this.dart);
  }

  set name(s: string) {
    if (api.grok_View_Set_Name == null)
      // @ts-ignore
      this._name = s;
    else
      api.grok_View_Set_Name(this.dart, s);
  }

  _onAdded() {
    api.grok_View_OnAdded(this.dart);
  }

  // to be used in [createByType].
  static readonly APPS = 'apps';
  static readonly SETTINGS = 'settings';
  static readonly WELCOME = 'welcome';
  static readonly SCRIPT = 'script';
  static readonly SKETCH = 'sketch';
  static readonly FORUM = 'forum';
  static readonly PROJECTS = 'projects';
  static readonly NOTEBOOKS = 'notebooks';
  static readonly HELP = 'help';
  static readonly OPEN_TEXT = 'text';
  static readonly DATABASES = 'databases';
  static readonly WEB_SERVICES = 'webservices';
  static readonly VIEW_LAYOUTS = 'layouts';
  static readonly FUNCTIONS = 'functions';
  static readonly DATA_CONNECTIONS = 'connections';
  static readonly DATA_JOB_RUNS = 'jobs';
  static readonly FILES = 'files';
  static readonly DATA_QUERY_RUNS = 'queryruns';
  static readonly EMAILS = 'emails';
  static readonly GROUPS = 'groups';
  static readonly MODELS = 'models';
  static readonly QUERIES = 'queries';
  static readonly SCRIPTS = 'scripts';
  static readonly USERS = 'users';
  static readonly PACKAGES = 'packages';
  static readonly PACKAGE_REPOSITORIES = 'repositories';
  static readonly JS_EDITOR = 'js';
  static readonly BROWSE = 'browse';

  static readonly ALL_VIEW_TYPES = [View.APPS, View.SETTINGS, View.WELCOME, View.SCRIPT, View.SKETCH,
    View.FORUM, View.PROJECTS, View.NOTEBOOKS, View.HELP, View.OPEN_TEXT, View.DATABASES,
    View.WEB_SERVICES, View.VIEW_LAYOUTS, View.FUNCTIONS, View.DATA_CONNECTIONS, View.DATA_JOB_RUNS,
    View.FILES, View.DATA_QUERY_RUNS, View.EMAILS, View.GROUPS, View.MODELS, View.QUERIES,
    View.SCRIPTS, View.USERS, View.PACKAGES, View.PACKAGE_REPOSITORIES, View.JS_EDITOR, View.BROWSE];
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
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new table view, and adds it to the workspace if specified */
  static create(table: DataFrame, addToWorkspace: boolean = true): TableView {
    return toJs(api.grok_TableView(table.dart, addToWorkspace));
  }

  /** Associated table, if it exists (for TableView), or null. */
  get table(): DataFrame | null { return toJs(api.grok_View_Get_Table(this.dart)); }

  /** Associated grid (spreadsheet). */
  get grid(): Grid { return toJs(api.grok_View_Get_Grid(this.dart)); }

  /** Returns existing, or creates a new filter group. */
  getFiltersGroup(options?: { createDefaultFilters?: boolean }): FilterGroup {
    return toJs(api.grok_TableView_GetFilters(this.dart, options?.createDefaultFilters ?? true));
  }

  get dataFrame(): DataFrame { return toJs(api.grok_View_Get_DataFrame(this.dart)); }
  set dataFrame(x: DataFrame) { api.grok_View_Set_DataFrame(this.dart, x.dart); }

  /** View toolbox that gets shown on the left, in the sidebar */
  get toolboxPage(): ToolboxPage {
    return new ToolboxPage(api.grok_View_Get_ToolboxPage(this.dart));
  }

  /** Adds a viewer of the specified type. */
  addViewer(v: ViewerType | string | Viewer, options?: any): Viewer {
    if (typeof v === 'string')
      v = toJs(api.grok_View_AddViewerByName(this.dart, v)) as Viewer;
    else
      api.grok_View_AddViewer(this.dart, v.dart);
    if (options)
      v.setOptions(options);
    api.grok_TableView_ProcessNewViewer(this.dart, v.dart);
    return v;
  }

  /** A dock node for this view.
   *  Use `grok.shell.dockManager` to manipulate it; {@link dockManager} is for controlling
   *  windows that reside inside this view.
   *  @type {DockNode} */
  get dockNode(): DockNode {
    return new DockNode(api.grok_View_Get_DockNode(this.dart));
  }

  /** 
   * View's dock manager. Only defined for DockView descendants such as {@link TableView}, UsersView, etc.
   * @type {DockManager} */
  get dockManager(): DockManager {
    return new DockManager(api.grok_View_Get_DockManager(this.dart));
  }

  /** This and some of the following methods are "softly deprecated" (will likely be deprecated in 1.21):
   *  deprecated: use addViewer(Viewer.histogram(options)).
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/histogram | histogram}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/histogram} */
  histogram(options?: Partial<IHistogramSettings>): Viewer {
    return this.addViewer(VIEWER.HISTOGRAM, options);
  }

  /** deprecated: use addViewer(Viewer.barChart(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/bar-chart | bar chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart} */
  barChart(options?: Partial<IBarChartSettings>): Viewer {
    return this.addViewer(VIEWER.BAR_CHART, options);
  }

  /** deprecated: use addViewer(Viewer.boxPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/box-plot | box plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/box-plot} */
  boxPlot(options?: Partial<IBoxPlotSettings>): Viewer {
    return this.addViewer(VIEWER.BOX_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.calendar(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/calendar | calendar}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/calendar}
   *  @param options
   *  @returns {Viewer} */
  calendar(options?: Partial<ICalendarSettings>): Viewer {
    return this.addViewer(VIEWER.CALENDAR, options);
  }

  /** deprecated: use addViewer(Viewer.corrPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/correlation-plot | correlation plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/corr-plot} */
  corrPlot(options?: Partial<ICorrelationPlotSettings>): Viewer {
    return this.addViewer(VIEWER.CORR_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.densityPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/density-plot | density plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/density-plot} */
  densityPlot(options?: Partial<IDensityPlotSettings>): Viewer {
    return this.addViewer(VIEWER.DENSITY_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.filters(options))
   *  Adds {@link https://datagrok.ai/help/visualize/viewers/filters | filters}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/filters} */
  filters(options?: Partial<IFiltersSettings>): Viewer {
    return this.addViewer(VIEWER.FILTERS, options);
  }

  /** deprecated: use addViewer(Viewer.form(options))
   *  Adds default {@link https://datagrok.ai/help/visualize/viewers/form | form}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/form} */
  form(options?: Partial<IFormSettings>): Viewer {
    return this.addViewer(VIEWER.FORM, options);
  }

  /** deprecated: use addViewer(Viewer.heatMap(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/heat-map | heat map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/heat-map} */
  heatMap(options?: Partial<IGridSettings>): Viewer {
    return this.addViewer(VIEWER.HEAT_MAP, options);
  }

  /** deprecated: use addViewer(Viewer.histogram(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/line-chart | line chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart}  */
  lineChart(options?: Partial<ILineChartSettings>): Viewer {
    return this.addViewer(VIEWER.LINE_CHART, options);
  }

  /** deprecated: use addViewer(Viewer.shapeMap(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/shape-map | shape map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/shape-map}  */
  shapeMap(options?: Partial<IMapViewerSettings>): Viewer {
    return this.addViewer(VIEWER.SHAPE_MAP, options);
  }

  /** deprecated: use addViewer(Viewer.markup(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/markup | markup viewer}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/markup} */
  markup(options?: Partial<IMarkupViewerSettings>): Viewer {
    return this.addViewer(VIEWER.MARKUP, options);
  }

  /** deprecated: use addViewer(Viewer.matrixPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/matrix-plot | matrix plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot} */
  matrixPlot(options?: Partial<IMatrixPlotSettings>): Viewer {
    return this.addViewer(VIEWER.MATRIX_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.networkDiagram(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/network-diagram | network diagram}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/network-diagram} */
  networkDiagram(options?: Partial<INetworkDiagramSettings>): Viewer {
    return this.addViewer(VIEWER.NETWORK_DIAGRAM, options);
  }

  /** deprecated: use addViewer(Viewer.pcPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pc-plot | parallel coordinates plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot} */
  pcPlot(options?: Partial<IPcPlotSettings>): Viewer {
    return this.addViewer(VIEWER.PC_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.pieChart(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pie-chart | pie chart}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pie-chart} */
  pieChart(options?: Partial<IPieChartSettings>): Viewer {
    return this.addViewer(VIEWER.PIE_CHART, options);
  }

  /** deprecated: use addViewer(Viewer.scatterPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/scatter-plot | scatter plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot} */
  scatterPlot(options?: Partial<IScatterPlotSettings>): ScatterPlotViewer {
    return <ScatterPlotViewer>this.addViewer(VIEWER.SCATTER_PLOT, options);
  }

  /** deprecated: use addViewer(Viewer.scatterPlot3d(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/3d-scatter-plot | 3D scatter plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot-3d} */
  scatterPlot3d(options?: Partial<IScatterPlot3dSettings>): Viewer {
    return this.addViewer(VIEWER.SCATTER_PLOT_3D, options);
  }

  /** deprecated: use addViewer(Viewer.statistics(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/statistics | statistics}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/statistics} */
  statistics(options?: Partial<IStatsViewerSettings>): Viewer {
    return this.addViewer(VIEWER.STATISTICS, options);
  }

  /** deprecated: use addViewer(Viewer.histogram(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tile-viewer | tile viewer}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer} */
  tileViewer(options?: Partial<ITileViewerSettings>): Viewer {
    return this.addViewer(VIEWER.TILE_VIEWER, options);
  }

  /** deprecated: use addViewer(Viewer.treeMap(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tree-map | tree map}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tree-map} */
  treeMap(options?: Partial<ITreeMapSettings>): Viewer {
    return this.addViewer(VIEWER.TREE_MAP, options);
  }

  /** deprecated: use addViewer(Viewer.trellisPlot(options))
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/trellis-plot | trellis plot}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot} */
  trellisPlot(options?: Partial<ITrellisPlotSettings>): Viewer {
    return this.addViewer(VIEWER.TRELLIS_PLOT, options);
  }

  /**
   *  Adds a {@link https://datagrok.ai/help/visualize/viewers/word-cloud | word cloud}.
   *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/word-cloud}
   *  @deprecated */
  wordCloud(options?: any): Viewer {
    return this.addViewer(VIEWER.WORD_CLOUD, options);
  }

  /** Resets view layout, leaving only grid visible. */
  resetLayout(): void {
    api.grok_View_ResetLayout(this.dart);
  }

  /** Detaches and closes the view. */
  detach(): void {
    api.grok_View_Detach(this.dart);
  }

  /** Detaches all viewers. */
  detachViewers(): void {
    api.grok_View_DetachViewers(this.dart);
  }

  /** Returns all viewers.
   * @type {Iterable.<Viewer>} */
  get viewers(): Iterable<Viewer> {
    return _toIterable(api.grok_View_Get_Viewers(this.dart));
  }

  get syncCurrentObject(): boolean { return api.grok_TableView_Get_SyncCurrentObject(this.dart); }
  set syncCurrentObject(x: boolean) { api.grok_TableView_Set_SyncCurrentObject(this.dart, x); }

  saveState(): string {
    return api.grok_TableView_SaveState(this.dart);
  }

  loadState(x: string, options?: IViewStateApplicationOptions): void {
    api.grok_TableView_LoadState(this.dart, x, options?.pickupColumnTags);
  }
}


export interface IViewStateApplicationOptions {
  pickupColumnTags?: boolean;
}


/** Script view */
export class ScriptView extends View {
  /** @constructs ScriptView */
  constructor(dart: any) {
    super(dart);
  }

  static create(script: Script): ScriptView {
    return new ScriptView(api.grok_ScriptView(script.dart));
  }
}


export class DockView extends View {
  constructor(dart: any) {
    super(dart);
  }

  initDock(): string {
    return api.grok_DockView_InitDock(this.dart);
  }

  _handleResize(): string {
    return api.grok_DockView_HandleResize(this.dart);
  }
}


export class BrowseView extends View {
  constructor(dart: any) {
    super(dart);
  }
  // TODO: add static method to return browse view
  get localTree(): TreeViewGroup { return api.grok_BrowseView_Get_LocalTree(this.dart); }
  get mainTree(): TreeViewGroup { return api.grok_BrowseView_Get_MainTree(this.dart); }

  get preview(): View | null { return toJs(api.grok_BrowseView_Get_Preview(this.dart)); }
  set preview(preview: View | null) { api.grok_BrowseView_Set_Preview(this.dart, preview?.dart); }

  get dockManager(): DockManager { return new DockManager(api.grok_BrowseView_Get_DockManager(this.dart)); }

  get showTree(): boolean { return api.grok_BrowseView_Get_ShowTree(this.dart); }
  set showTree(x: boolean) { api.grok_BrowseView_Set_ShowTree(this.dart, x); }
}


export class ViewLayout extends Entity {

  /** @constructs ViewLayout */
  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromJson(json));
  }

  static fromViewState(state: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromViewState(state));
  }

  get viewState(): string {
    return api.grok_ViewLayout_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewLayout_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewLayout_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewLayout_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewLayout_ToJson(this.dart);
  }

  get columns(): ColumnInfo[] {
    return toJs(api.grok_ViewLayout_Get_Columns(this.dart));
  }

}

export class ViewInfo extends Entity {

  /** @constructs ViewInfo */
  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromJson(json));
  }

  static fromViewState(state: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromViewState(state));
  }

  get table() : TableInfo {
    return toJs(api.grok_ViewInfo_Get_Table(this.dart));
  }

  /** Only defined within the context of the OnViewLayoutXXX events */
  get view(): View {
    return toJs(api.grok_ViewInfo_Get_View(this.dart));
  }

  get viewState(): string {
    return api.grok_ViewInfo_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewInfo_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewInfo_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewInfo_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewInfo_ToJson(this.dart);
  }
}

/** Represents a virtual view, where visual elements are created only when user
 * scrolls them into view. */
export class VirtualView {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(verticalScroll: boolean = true, maxCols: number = 100): VirtualView {
    return new VirtualView(api.grok_VirtualItemView(verticalScroll, maxCols));
  }

  /** Visual root. */
  get root(): HTMLElement {
    return api.grok_VirtualItemView_Get_Root(this.dart);
  }

  /** Sets the number of elements, and a function that renders i-th element. */
  setData(length: number, renderer: (index: number) => HTMLElement): void {
    api.grok_VirtualItemView_SetData(this.dart, length, renderer);
  }

  /** Refreshes i-th element without refreshing the whole view */
  refreshItem(i: number): void {
    api.grok_VirtualItemView_RefreshItem(this.dart, i);
  }
}

