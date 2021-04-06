var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
  Object.defineProperty(o, 'default', { enumerable: true, value: v });
}) : function(o, v) {
  o['default'] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) if (k !== 'default' && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
  __setModuleDefault(result, mod);
  return result;
};
var __importDefault = (this && this.__importDefault) || function (mod) {
  return (mod && mod.__esModule) ? mod : { 'default': mod };
};
define(['require', 'exports', './const', '../ui', './viewer', './docking', './grid', './widgets', './entities', './wrappers', './utils', 'cash-dom'], function (require, exports, const_1, ui, viewer_1, docking_1, grid_1, widgets_1, entities_1, wrappers_1, utils_1, cash_dom_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.VirtualView = exports.ViewLayout = exports.ScriptView = exports.ProjectsView = exports.DataSourceCardView = exports.TableView = exports.View = exports.ViewBase = void 0;
  ui = __importStar(ui);
  cash_dom_1 = __importDefault(cash_dom_1);
  let api = window;
  /**
     * Subclass ViewBase to implement a Datagrok view in JavaScript.
     * */
  class ViewBase {
    /**
         * @constructs ViewBase
         * @param {Object} params - URL parameters.
         * @param {string} path - URL path.
         * @param {boolean} createHost - Create JS host wrapper. */
    constructor(params = null, path = '', createHost = true) {
      if (createHost)
        this.d = api.grok_View_CreateJsViewHost(this);
      this._root = ui.div([], 'grok-view');
      this._root.tabIndex = 0;
      /** @type {StreamSubscription[]} */
      this.subs = []; // stream subscriptions - will be canceled when the view is detached
      this._closing = false;
    }
    /** @type {HTMLElement} */
    get root() {
      return this._root;
    }
    /** View type
         * @type {string} */
    get type() {
      return 'js-view-base';
    }
    /** @returns {string|null} View help URL. */
    get helpUrl() {
      return null;
    }
    /** View name. It gets shown in the tab handle.
         * @type {string} */
    get name() {
      // @ts-ignore
      return api.grok_View_Get_Name == null ? this._name : api.grok_View_Get_Name(this.d);
    }
    set name(s) {
      // @ts-ignore
      if (api.grok_View_Set_Name == null)
      // @ts-ignore
        this._name = s;
      else
        api.grok_View_Set_Name(this.d, s);
    }
    /** @type {string} */
    get description() {
      return '';
    }
    set description(s) {
    }
    /** @type {Object} */
    get entity() {
      return null;
    }
    set entity(e) {
    }
    /** View toolbox.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/toolbox}
         * @type {HTMLElement} */
    get toolbox() {
      return api.grok_View_Get_Toolbox(this.d);
    }
    set toolbox(x) {
      api.grok_View_Set_Toolbox(this.d, x);
    }
    /** View menu.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon}
         *  @type {Menu} */
    get ribbonMenu() {
      return new widgets_1.Menu(api.grok_View_Get_RibbonMenu(this.d));
    }
    set ribbonMenu(menu) {
      api.grok_View_Set_RibbonMenu(this.d, menu.d);
    }
    get closing() {
      return this._closing;
    }
    set closing(c) {
      this._closing = c;
    }
    /** Sets custom view panels on the ribbon.
         * @param {Array<Array<HTMLElement>>} panels
         * @param {boolean} clear Clear all previous before setup
         * Sample: {@link https://public.datagrok.ai/js/samples/ui/views/ribbon} */
    setRibbonPanels(panels, clear = false) {
      api.grok_View_SetRibbonPanels(this.d, panels, clear);
    }
    /** @returns {HTMLElement} View icon. */
    getIcon() {
      return null;
    }
    /** @returns {Object} Viewer state map. */
    saveStateMap() {
      return null;
    }
    /** Load view state map.
         * @param {Object} stateMap - State map. */
    loadStateMap(stateMap) {
    }
    /**
         * View URI, relative to the platform root. See also {@link basePath}
         * @type {string} */
    get path() {
      return '';
    }
    set path(s) {
    }
    /** Handles URL path.
         * @param  {string} path - URL path. */
    handlePath(path) {
    }
    /** Checks if URL path is acceptable.
         * @returns {boolean} "true" if path is acceptable, "false" otherwise.
         * @param {string} path - URL path. */
    acceptsPath(path) {
      return false;
    }
    /**
         * Appends an item to this view. Use {@link appendAll} for appending multiple elements.
         * @param {Object} item */
    append(item) {
      return this.appendAll([ui.render(item)]);
    }
    /**
         * Appends multiple elements this view. Use {@link append} for appending a single element.
         * @param {object[]} items */
    appendAll(items) {
      return ui.appendAll(this.root, items.map(ui.render));
    }
    /** Detaches this view. */
    detach() {
      this.subs.forEach((sub) => sub.unsubscribe());
    }
    /** Closes this view. */
    close() {
      this._closing = true;
      api.grok_View_Close(this.d);
    }
  }
  exports.ViewBase = ViewBase;
  /**
     * A view is typically docked in the main document area of the Grok platform.
     * See [TableView], [SketchView], etc
     */
  class View extends ViewBase {
    /** @constructs View */
    constructor(d) {
      super(null, '', false);
      this.d = d;
    }
    static fromDart(d) {
      let type = api.grok_View_Get_Type(d);
      if (type === const_1.VIEW_TYPE.TABLE_VIEW)
        return new TableView(d);
      else
        return new View(d);
    }
    /** Creates a new empty view.
         * @param {string | ElementOptions | null} options
         * @returns {View} */
    static create(options) {
      // @ts-ignore
      let v = api.grok_View == null ? new View(null) : new View(api.grok_View());
      utils_1._options(v.root, 'ui-panel');
      utils_1._options(v.root, options);
      return v;
    }
    get box() {
      return cash_dom_1.default(this.root).hasClass('ui-box');
    }
    set box(b) {
      let r = cash_dom_1.default(this.root);
      r.removeClass('ui-panel').removeClass('ui-box');
      r.addClass(b ? 'ui-box' : 'ui-panel');
    }
    get root() {
      // @ts-ignore
      if (api.grok_View_Get_Root == null)
        return this._root;
      return api.grok_View_Get_Root(this.d);
    }
    get type() {
      return api.grok_View_Get_Type(this.d);
    }
    get path() {
      return api.grok_View_Get_Path(this.d);
    }
    set path(s) {
      api.grok_View_Set_Path(this.d, s);
    }
    /**
         *  View type URI. Note that {@link path} is specific to the instance of the view.
         *  @type {string} */
    get basePath() {
      return api.grok_View_Get_BasePath(this.d);
    }
    set basePath(s) {
      api.grok_View_Set_BasePath(this.d, s);
    }
    get description() {
      return api.grok_View_Get_Description(this.d);
    }
    set description(s) {
      api.grok_View_Set_Description(this.d, s);
    }
    /**
         * Loads previously saved view layout. Only applicable to certain views, such as {@link TableView}.
         *  See also {@link saveLayout}
         *  @param {ViewLayout} layout */
    loadLayout(layout) {
      return api.grok_View_Load_Layout(this.d, layout.d);
    }
    /**
         *  Saves view layout as a string. Only applicable to certain views, such as {@link TableView}.
         *  See also {@link loadLayout}
         *  @returns {ViewLayout} */
    saveLayout() {
      return new ViewLayout(api.grok_View_Save_Layout(this.d));
    }
  }
  exports.View = View;
  /**
     * A {@link View} that is associated with a {@link DataFrame} and exposes
     * exploratory data analysis functionality. This view gets opened whenever
     * a new table is added to the workspace when a user drag-and-drops a CSV file,
     * or opens a table in any other way.
     * @extends View
     */
  class TableView extends View {
    /** @constructs TableView */
    constructor(d) {
      super(d);
    }
    /** Creates a new table view.
         * @param {DataFrame} table
         * @returns {TableView} */
    static create(table) {
      return new TableView(api.grok_TableView(table.d));
    }
    /** Associated table, if it exists (for TableView), or null.
         *  @type {DataFrame} */
    get table() {
      return wrappers_1.toJs(api.grok_View_Get_Table(this.d));
    }
    /** @type {Grid} */
    get grid() {
      return new grid_1.Grid(api.grok_View_Get_Grid(this.d));
    }
    /** @type {DataFrame} */
    get dataFrame() {
      return wrappers_1.toJs(api.grok_View_Get_DataFrame(this.d));
    }
    set dataFrame(x) {
      api.grok_View_Set_DataFrame(this.d, x.d);
    }
    /** View toolbox that gets shown on the left, in the sidebar.
         *  @type {ToolboxPage} */
    get toolboxPage() {
      return new widgets_1.ToolboxPage(api.grok_View_Get_ToolboxPage(this.d));
    }
    /** Adds a viewer of the specified type.
         * @param {string | Viewer} v
         * @param options
         * @returns {Viewer} */
    addViewer(v, options = null) {
      if (typeof v === 'string')
        v = new viewer_1.Viewer(api.grok_View_AddViewerByName(this.d, v));
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
    get dockNode() {
      return new docking_1.DockNode(api.grok_View_Get_DockNode(this.d));
    }
    /**
         * View's dock manager. Only defined for DockView descendants such as {@link TableView}, UsersView, etc.
         * @type {DockManager} */
    get dockManager() {
      return new docking_1.DockManager(api.grok_View_Get_DockManager(this.d));
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/histogram | histogram}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/histogram}
         *  @param options
         *  @returns {Viewer} */
    histogram(options = null) {
      return this.addViewer(const_1.VIEWER.HISTOGRAM, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/bar-chart | bar chart}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart}
         *  @param options
         *  @returns {Viewer} */
    barChart(options = null) {
      return this.addViewer(const_1.VIEWER.BAR_CHART, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/box-plot | box plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/box-plot}
         *  @param options
         *  @returns {Viewer} */
    boxPlot(options = null) {
      return this.addViewer(const_1.VIEWER.BOX_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/calendar | calendar}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/calendar}
         *  @param options
         *  @returns {Viewer} */
    calendar(options = null) {
      return this.addViewer(const_1.VIEWER.CALENDAR, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/correlation-plot | correlation plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/corr-plot}
         *  @param options
         *  @returns {Viewer} */
    corrPlot(options = null) {
      return this.addViewer(const_1.VIEWER.CORR_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/density-plot | density plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/density-plot}
         *  @param options
         *  @returns {Viewer} */
    densityPlot(options = null) {
      return this.addViewer(const_1.VIEWER.DENSITY_PLOT, options);
    }
    /**
         *  Adds {@link https://datagrok.ai/help/visualize/viewers/filters | filters}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/filters}
         *  @param options
         *  @returns {Viewer} */
    filters(options = null) {
      return this.addViewer(const_1.VIEWER.FILTERS, options);
    }
    /**
         *  Adds default {@link https://datagrok.ai/help/visualize/viewers/form | form}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/form}
         *  @param options
         *  @returns {Viewer} */
    form(options = null) {
      return this.addViewer(const_1.VIEWER.FORM, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/google-map | geographical map}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/google-map}
         *  @param options
         *  @returns {Viewer} */
    googleMap(options = null) {
      return this.addViewer(const_1.VIEWER.GOOGLE_MAP, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/heat-map | heat map}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/heat-map}
         *  @param options
         *  @returns {Viewer} */
    heatMap(options = null) {
      return this.addViewer(const_1.VIEWER.HEAT_MAP, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/line-chart | line chart}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart}
         *  @param options
         *  @returns {Viewer} */
    lineChart(options = null) {
      return this.addViewer(const_1.VIEWER.LINE_CHART, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/shape-map | shape map}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/shape-map}
         *  @param options
         *  @returns {Viewer} */
    shapeMap(options = null) {
      return this.addViewer(const_1.VIEWER.SHAPE_MAP, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/markup | markup viewer}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/markup}
         *  @param options
         *  @returns {Viewer} */
    markup(options = null) {
      return this.addViewer(const_1.VIEWER.MARKUP, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/matrix-plot | matrix plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot}
         *  @param options
         *  @returns {Viewer} */
    matrixPlot(options = null) {
      return this.addViewer(const_1.VIEWER.MATRIX_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/network-diagram | network diagram}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/network-diagram}
         *  @param options
         *  @returns {Viewer} */
    networkDiagram(options = null) {
      return this.addViewer(const_1.VIEWER.NETWORK_DIAGRAM, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pc-plot | parallel coordinates plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot}
         *  @param options
         *  @returns {Viewer} */
    pcPlot(options = null) {
      return this.addViewer(const_1.VIEWER.PC_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/pie-chart | pie chart}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/pie-chart}
         *  @param options
         *  @returns {Viewer} */
    pieChart(options = null) {
      return this.addViewer(const_1.VIEWER.PIE_CHART, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/scatter-plot | scatter plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot}
         *  @param options
         *  @returns {Viewer} */
    scatterPlot(options = null) {
      return this.addViewer(const_1.VIEWER.SCATTER_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/3d-scatter-plot | 3D scatter plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot-3d}
         *  @param options
         *  @returns {Viewer} */
    scatterPlot3d(options = null) {
      return this.addViewer(const_1.VIEWER.SCATTER_PLOT_3D, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/statistics | statistics}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/statistics}
         *  @param options
         *  @returns {Viewer} */
    statistics(options = null) {
      return this.addViewer(const_1.VIEWER.STATISTICS, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tile-viewer | tile viewer}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer}
         *  @param options
         *  @returns {Viewer} */
    tileViewer(options = null) {
      return this.addViewer(const_1.VIEWER.TILE_VIEWER, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/tree-map | tree map}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/tree-map}
         *  @param options
         *  @returns {Viewer} */
    treeMap(options = null) {
      return this.addViewer(const_1.VIEWER.TREE_MAP, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/trellis-plot | trellis plot}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot}
         *  @param options
         *  @returns {Viewer} */
    trellisPlot(options = null) {
      return this.addViewer(const_1.VIEWER.TRELLIS_PLOT, options);
    }
    /**
         *  Adds a {@link https://datagrok.ai/help/visualize/viewers/word-cloud | word cloud}.
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/word-cloud}
         *  @param options
         *  @returns {Viewer} */
    wordCloud(options = null) {
      return this.addViewer(const_1.VIEWER.WORD_CLOUD, options);
    }
    /** Resets view layout, leaving only grid visible. */
    resetLayout() {
      api.grok_View_ResetLayout(this.d);
    }
    /** Detaches and closes the view. */
    detach() {
      api.grok_View_Detach(this.d);
    }
    /** Detaches all viewers. */
    detachViewers() {
      api.grok_View_DetachViewers(this.d);
    }
    /** Returns all viewers.
         * @type {Iterable.<Viewer>} */
    get viewers() {
      return utils_1._toIterable(api.grok_View_Get_Viewers(this.d));
    }
  }
  exports.TableView = TableView;
  /** Base view for working with a collection of objects that reside on the server.
     *  Typically, results are filtered by applying AND operation between two
     *  filters: {@link permanentFilter} (which is set programmatically and is not visible)
     *  and {@link searchValue} entered by the user.
     *
     *  More details on the smart search syntax: {@link https://datagrok.ai/help/overview/smart-search}
     *
     * @extends View */
  class DataSourceCardView extends View {
    /** @constructs DataSourceCardView */
    constructor(d) {
      super(d);
    }
    /**
         * User-specified {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
         * @type {string} */
    get searchValue() {
      return api.grok_DataSourceCardView_Get_SearchValue(this.d);
    }
    set searchValue(s) {
      api.grok_DataSourceCardView_Set_SearchValue(this.d, s);
    }
    /** Programmatically defined invisible
         * {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
         *  @type {string} */
    get permanentFilter() {
      return api.grok_DataSourceCardView_Get_PermanentFilter(this.d);
    }
    set permanentFilter(s) {
      api.grok_DataSourceCardView_Set_PermanentFilter(this.d, s);
    }
  }
  exports.DataSourceCardView = DataSourceCardView;
  /** Projects view */
  class ProjectsView extends DataSourceCardView {
    /** @constructs ProjectsView */
    constructor(d) {
      super(d);
    }
    static create(params) {
      return new ProjectsView(api.grok_ProjectsView(params));
    }
  }
  exports.ProjectsView = ProjectsView;
  /** Script view */
  class ScriptView extends View {
    /** @constructs ScriptView */
    constructor(d) {
      super(d);
    }
    static create(script) {
      return new ScriptView(api.grok_ScriptView(script.d));
    }
  }
  exports.ScriptView = ScriptView;
  class ViewLayout extends entities_1.Entity {
    /** @constructs ViewLayout */
    constructor(d) {
      super(d);
    }
    static fromJson(json) {
      return new ViewLayout(api.grok_ViewLayout_FromJson(json));
    }
    static fromViewState(state) {
      return new ViewLayout(api.grok_ViewLayout_FromViewState(state));
    }
    /** Only defined within the context of the OnViewLayoutXXX events */
    get view() {
      return api.grok_ViewLayout_Get_View(this.d);
    }
    get viewState() {
      return api.grok_ViewLayout_Get_ViewState(this.d);
    }
    set viewState(state) {
      api.grok_ViewLayout_Set_ViewState(this.d, state);
    }
    getUserDataValue(key) {
      return api.grok_ViewLayout_Get_UserDataValue(this.d, key);
    }
    setUserDataValue(key, value) {
      return api.grok_ViewLayout_Set_UserDataValue(this.d, key, value);
    }
    toJson() {
      return api.grok_ViewLayout_ToJson(this.d);
    }
  }
  exports.ViewLayout = ViewLayout;
  class VirtualView {
    constructor(d) {
      this.d = d;
    }
    static create(verticalScroll = true, maxCols = 100) {
      return new VirtualView(api.grok_VirtualItemView(verticalScroll, maxCols));
    }
    get root() {
      return api.grok_VirtualItemView_Get_Root(this.d);
    }
    setData(length, renderer) {
      api.grok_VirtualItemView_SetData(this.d, length, renderer);
    }
  }
  exports.VirtualView = VirtualView;
});
