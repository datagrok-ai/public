/**
 * A view is typically docked in the main document area of the Grok platform.
 * See [TableView], [SketchView], etc
 */
import {
    VIEW_TYPE, VIEWER,
} from "./const";
import {DataFrame} from "./dataframe.js";
import * as ui from "./../ui";
import {Viewer} from "./viewer";
import {DockNode, DockManager} from "./docking";

export class View {
    constructor(d) { this.d = d; }

    static fromDart(d) {
        let type = grok_View_Get_Type(d);
        if (type === VIEW_TYPE.TABLE_VIEW)
            return new TableView(d);
        else
            return new View(d);

    }
    static create() { let v = new View(grok_View()); ui._class(v.root, 'grok-default-view'); return v; }

    get root() { return grok_View_Get_Root(this.d); }

    append(x) { return ui.appendAll(this.root, [x]); }
    appendAll(x) { return ui.appendAll(this.root, x); }

    get type() { return grok_View_Get_Type(this.d); }

    get name() { return grok_View_Get_Name(this.d); }
    set name(s) { return grok_View_Set_Name(this.d, s); }

    get path() { return grok_View_Get_Path(this.d); }
    set path(s) { return grok_View_Set_Path(this.d, s); }

    get basePath() { return grok_View_Get_BasePath(this.d); }
    set basePath(s) { return grok_View_Set_BasePath(this.d, s); }

    get description() { return grok_View_Get_Description(this.d); }
    set description(s) { return grok_View_Set_Description(this.d, s); }

    get viewType() { return grok_View_Get_Type(this.d); }

    get table() { return new DataFrame(grok_View_Get_Table(this.d)); }

    get toolbox() { return grok_View_Get_Toolbox(this.d); }
    set toolbox(x) { return grok_View_Set_Toolbox(this.d, x); }

    get ribbonMenu() { return new Menu(grok_View_Get_RibbonMenu(this.d)); }
    set ribbonMenu(menu) { grok_View_Set_RibbonMenu(this.d, menu.d); }

    loadLayout(layout) { return grok_View_Load_Layout(this.d, layout.d);  }
    saveLayout() { return new ViewLayout(grok_View_Save_Layout(this.d)); }

    setRibbonPanels(panels) { grok_View_SetRibbonPanels(this.d, panels); }

    /** Adds a viewer of the specified type.
     * @param {string} viewerType
     * @returns {Viewer} */
    addViewer(viewerType, options = null) {
        let v = new Viewer(grok_View_AddViewer(this.d, viewerType));
        if (options !== null)
            v.options(options);
        return v;
    }

    /** A dock node for this view.
     *  Use `grok.shell.dockManager` to manipulate it; {@link dockManager} is for controlling
     *  widows that reside inside this view.
     *  @returns {DockNode} */
    get dockNode() { return new DockNode(grok_View_Get_DockNode(this.d)); }

    /** @returns {DockManager} - only for DockView descendants (TableView, Users, etc) */
    get dockManager() { return new DockManager(grok_View_Get_DockManager(this.d)); }

    close() { grok_View_Close(this.d); }
}


export class TableView extends View {
    constructor(d) { super(d); }

    get grid() { return new Grid(grok_View_Get_Grid(this.d)); }

    get dataFrame() { return new DataFrame(grok_View_Get_DataFrame(this.d)); }
    set dataFrame(x) { grok_View_Set_DataFrame(this.d, x.d); }

    get toolboxPage() { return new ToolboxPage(grok_View_Get_ToolboxPage(this.d)); }

    histogram      (options = null) { return this.addViewer(VIEWER.HISTOGRAM, options); }
    barChart       (options = null) { return this.addViewer(VIEWER.BAR_CHART, options); }
    boxPlot        (options = null) { return this.addViewer(VIEWER.BOX_PLOT, options); }
    calendar       (options = null) { return this.addViewer(VIEWER.CALENDAR, options); }
    corrPlot       (options = null) { return this.addViewer(VIEWER.CORR_PLOT, options); }
    densityPlot    (options = null) { return this.addViewer(VIEWER.DENSITY_PLOT, options); }
    filters        (options = null) { return this.addViewer(VIEWER.FILTERS, options); }
    form           (options = null) { return this.addViewer(VIEWER.FORM, options); }
    globe          (options = null) { return this.addViewer(VIEWER.GLOBE, options); }
    googleMap      (options = null) { return this.addViewer(VIEWER.GOOGLE_MAP, options); }
    heatMap        (options = null) { return this.addViewer(VIEWER.HEAT_MAP, options); }
    lineChart      (options = null) { return this.addViewer(VIEWER.LINE_CHART, options); }
    shapeMap       (options = null) { return this.addViewer(VIEWER.SHAPE_MAP, options); }
    markup         (options = null) { return this.addViewer(VIEWER.MARKUP, options); }
    matrixPlot     (options = null) { return this.addViewer(VIEWER.MATRIX_PLOT, options); }
    networkDiagram (options = null) { return this.addViewer(VIEWER.NETWORK_DIAGRAM, options); }
    pcPlot         (options = null) { return this.addViewer(VIEWER.PC_PLOT, options); }
    pieChart       (options = null) { return this.addViewer(VIEWER.PIE_CHART, options); }
    scatterPlot    (options = null) { return this.addViewer(VIEWER.SCATTER_PLOT, options); }
    scatterPlot3d  (options = null) { return this.addViewer(VIEWER.SCATTER_PLOT_3D, options); }
    statistics     (options = null) { return this.addViewer(VIEWER.STATISTICS, options); }
    tileViewer     (options = null) { return this.addViewer(VIEWER.TILE_VIEWER, options); }
    treeMap        (options = null) { return this.addViewer(VIEWER.TREE_MAP, options); }
    trellisPlot    (options = null) { return this.addViewer(VIEWER.TRELLIS_PLOT, options); }
    wordCloud      (options = null) { return this.addViewer(VIEWER.WORD_CLOUD, options); }

    resetLayout() { grok_View_ResetLayout(this.d); }

    detach() { grok_View_Detach(this.d); }
    detachViewers() { grok_View_DetachViewers(this.d); }
}


export class DataSourceCardView extends View {
    set searchValue(s) { return grok_DataSourceCardView_Set_SearchValue(this.d, s); }
    get permanentFilter() { return grok_DataSourceCardView_Get_PermanentFilter(this.d); }
    set permanentFilter(s) { return grok_DataSourceCardView_Set_PermanentFilter(this.d, s); }
}


export class ProjectsView extends DataSourceCardView {
    static create(params) { return new ProjectsView(grok_ProjectsView(params)); }
}


export class ViewLayout {
    constructor(d) { this.d = d; }

    static fromJson(json) { return new ViewLayout(grok_ViewLayout_FromJson(json)); }

    toJson() { return grok_ViewLayout_ToJson(this.d); }
}


export class VirtualView {
    constructor(d) { this.d = d; }

    static create() { return new VirtualView(grok_VirtualItemView()); }

    get root() { return grok_VirtualItemView_Get_Root(this.d); }
    setData(length, renderer) { grok_VirtualItemView_SetData(this.d, length, renderer); }
}
