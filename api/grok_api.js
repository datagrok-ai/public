class Utils {
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


/** Strongly-typed property associated with an object.
 *  Used for reflection, serialization, UI generation, and other introspection-dependent tasks. */
class Property {
    constructor(d) { this.d = d; }

    get get() { return grok_Property_Get_Get(this.d); }
    set get(x) { return grok_Property_Set_Get(this.d, x); }

    get set() { return grok_Property_Get_Set(this.d); }
    set set(x) { return grok_Property_Set_Set(this.d, x); }

    get name() { return grok_Property_Get_Name(this.d); }
    set name(s) { grok_Property_Set_Name(this.d, s); }
    
    get propertyType() { return grok_Property_Get_PropertyType(this.d); }
    set propertyType(s) { grok_Property_Set_PropertyType(this.d, s); }
    
    get semType() { return grok_Property_Get_SemType(this.d); }
    set semType(s) { grok_Property_Set_SemType(this.d, s); }
    
    get description() { return grok_Property_Get_Description(this.d); }
    set description(s) { grok_Property_Set_Description(this.d, s); }

    get defaultValue() { return grok_Property_Get_DefaultValue(this.d); }
    set defaultValue(s) { grok_Property_Set_DefaultValue(this.d, s); }

    static create(name, type, getter, setter, defaultValue) { return new Property(grok_Property(name, type, getter, setter, defaultValue)); }
    static int(name, getter, setter, defaultValue) { return Property.create(name, TYPE_INT, getter, setter, defaultValue); }
    static float(name, getter, setter, defaultValue) { return Property.create(name, TYPE_FLOAT, getter, setter, defaultValue); }
    static string(name, getter, setter, defaultValue) { return Property.create(name, TYPE_STRING, getter, setter, defaultValue); }
}


/** Efficient bit storage and manipulation. */
class BitSet {
    constructor(d) { this.d = d; }

    static create(length) { return new BitSet(grok_BitSet()); }

    /** Number of bits in a bitset */
    get length() { return grok_BitSet_Get_Length(this.d); }

    /** Number of set bits */
    get trueCount() { return grok_BitSet_Get_TrueCount(this.d); }

    /** Number of unset bits */
    get falseCount() { return grok_BitSet_Get_FalseCount(this.d); }

    /** Clones a bitset */
    clone() { return new BitSet(grok_BitSet_Clone(this.d)); }

    /** Inverts a bitset */
    invert() { return grok_BitSet_Invert(this.d); }

    /** Sets all bits to x */
    setAll(x) { return grok_BitSet_SetAll(this.d, x); }

    /** Finds the first index of value x, going forward from i-th position */
    findNext(i, x) { return grok_BitSet_FindNext(this.d, i, x); }

    /** Finds the first index of value x, going forward from i-th position */
    findPrev(i, x) { return grok_BitSet_FindPrev(this.d, i, x); }

    /** Gets i-th bit */
    get(i) { return grok_BitSet_GetBit(this.d, i); }

    /** Sets i-th bit to x */
    set(i, x) { return grok_BitSet_SetBit(this.d, i, x); }

    copyFrom(b) { return grok_BitSet_CopyFrom(this.d, b.d); }

    /** Registers a callback that gets invoked when the bitset is changed */
    onChanged(callback) { return _sub(grok_BitSet_OnChanged(this.d, callback)); }
}


/** Strongly-typed column. */
class Column {
    constructor(d) { this.d = d; }
    static fromStrings(name, list) { return new Column(grok_Column_FromStrings(name, list)); }
    static fromType(type, name = null, length = 0) { return new Column(grok_Column_FromType(type, name, length)); }

    /** [array] will be not be copied and will be used as column's storage */
    static fromInt32Array(name, array, length = null) { return new Column(grok_Column_FromInt32Array(name, array, length)); }

    /** [array] will be not be copied and will be used as column's storage */
    static fromFloat32Array(name, array, length = null) { return new Column(grok_Column_FromFloat32Array(name, array, length)); }

    static int(name, length = 0) { return Column.fromType(TYPE_INT, name, length); }
    static float(name, length = 0) { return Column.fromType(TYPE_FLOAT, name, length); }
    static string(name, length = 0) { return Column.fromType(TYPE_STRING, name, length); }
    static bool(name, length = 0) { return Column.fromType(TYPE_BOOL, name, length); }
    static dateTime(name, length = 0) { return Column.fromType(TYPE_DATE_TIME, name, length); }

    get type() { return grok_Column_Get_Type(this.d); }

    get length() { return grok_Column_Get_Length(this.d); }

    get dataFrame() { return new DataFrame(grok_Column_Get_DataFrame(this.d)); }

    get semType() { return grok_Column_Get_SemType(this.d); }
    set semType(s) { grok_Column_Set_SemType(this.d, s); }

    get name() { return grok_Column_Get_Name(this.d); }
    set name(s) { grok_Column_Set_Name(this.d, s); }

    get(i) { return grok_Column_GetValue(this.d, i); }
    set(i, x) { grok_Column_SetValue(this.d, i, x); }

    isNone(i) { return grok_Column_IsNone(this.d, i); }

    getTag(tag) { return grok_Column_Get_Tag(this.d, tag); }
    setTag(tag, value) { grok_Column_Set_Tag(this.d, tag, value); }

    compact() { return grok_Column_Compact(this.d); }

    toList() { return grok_Column_ToList(this.d); }

    get categories() { return grok_Column_Categories(this.d); }
    get min() { return grok_Column_Min(this.d); }
    get max() { return grok_Column_Max(this.d); }
    get stats() { return Stats.fromColumn(this); }

    *values() {
        for (let i = 0; i < this.length; i++) {
            yield this.get(i);
        }
    }
}


class Stats {
    constructor(d) { this.d = d; }

    static fromColumn(col) { return new Stats(grok_Stats_FromColumn(col.d)); }

    /** Total number of values (including missing values). */
    get totalCount() { return grok_Stats_Get_TotalCount(this.d); }

    /** Number of missing (empty) values. */
    get missingValueCount() { return grok_Stats_Get_MissingValueCount(this.d); }

    /** Number of non-empty values. */
    get valueCount() { return grok_Stats_Get_ValueCount(this.d); }

    get min() { return grok_Stats_Get_Min(this.d); }
    get max() { return grok_Stats_Get_Max(this.d); }
    get sum() { return grok_Stats_Get_Sum(this.d); }
    get avg() { return grok_Stats_Get_Avg(this.d); }
    get stdev() { return grok_Stats_Get_Stdev(this.d); }
    get variance() { return grok_Stats_Get_Variance(this.d); }
    get skew() { return grok_Stats_Get_Skew(this.d); }
    get kurt() { return grok_Stats_Get_Kurt(this.d); }
    get med() { return grok_Stats_Get_Med(this.d); }
    get q1() { return grok_Stats_Get_Q1(this.d); }
    get q2() { return grok_Stats_Get_Q2(this.d); }
    get q3() { return grok_Stats_Get_Q3(this.d); }
}

/** Columns in a [DataFrame]. */
class ColumnList {
    constructor(d) { this.d = d; }

    get length() { return grok_ColumnList_Length(this.d); }
    byName(name) { return new Column(grok_ColumnList_ByName(this.d, name)); }
    names() { return grok_ColumnList_Names(this.d); }
    toList() { return this.names().map(name => this.byName(name)); }
    add(column, notify = false) { grok_ColumnList_Add(this.d, column.d, notify); return column; }
    addNew(name, type) { return new Column(grok_ColumnList_AddNew(this.d, name, type)); }
    remove(name) { grok_ColumnList_Remove(this.d, name); }
    contains(name) { return grok_ColumnList_Contains(this.d, name); }
}


/** Represents a row. Allows for quick property access like "row.height". */
class Row {
    constructor(table, idx) {
        this.table = table;
        this.idx = idx;
        let setables = ['table', 'idx'];
        return new Proxy(this, {
            set(target, name, value) {
                if (setables.includes(name)) {
                    target[name] = value;
                    return true;
                }
                target.table.set(name, target.idx, value);
                return true;
            },
            get(target, name) {
                if (setables.includes(name))
                    return target[name];
                return target.table.get(name, target.idx);
            }
        });
    }

    get(name, idx) { return this.table.getCol(name).get(idx); }
}


/**
 * Represents rows of the [DataFrame].
 *
 * Refrain from accessing data via [RowList] and [Row] in performance-critical scenarios.
 * To maximize performance, get values via [DataFrame.columns], instead.
 */
class RowList {
    constructor(table, d) { this.table = table; this.d = d; }

    removeAt(idx, count = 1, notify = true) { grok_RowList_RemoveAt(this.d, idx, count, notify); }
    insertAt(idx, count = 1, notify = true) { grok_RowList_InsertAt(this.d, idx, count, notify); }
    addNew(values = null, notify = true) { return new Row(this.table, grok_RowList_AddNew(this.d, values, notify));}

    *iterator() {
        for (let i = 0; i < this.table.rowCount; i++)
            yield new Row(this.table, i);
    }

    setValues(idx, values) { grok_RowList_SetValues(this.d, idx, values); }
}


class Cell {
    constructor(d) { this.d = d; }

    get dataFrame() { return new DataFrame(grok_Cell_Get_DataFrame(this.d)); }
    get row() { return new Row(this.dataFrame, this.rowIndex); }
    get rowIndex() { return grok_Cell_Get_RowIndex(this.d); }
    get column() { return new Column(grok_Cell_Get_Column(this.d)); }
    get value() { return grok_Cell_Get_Value(this.d); }
}


/**
 * DataFrame is a high-performance, easy to use tabular structure with
 * strongly-typed columns of different types.
 */
class DataFrame {
    constructor(d) { this.d = d; }
    static create(rowCount = 0) { return new DataFrame(grok_DataFrame(rowCount)); }
    static fromColumns(columns) { return new DataFrame(grok_DataFrame_FromColumns(columns.map((c) => c.d))); }
    static fromCsv(csv) { return grok.parseCsv(csv); }
    static fromJson(json) { return new DataFrame(grok_DataFrame_FromJson(json)); }

    toString() { return `${this.name} (${this.rowCount} rows, ${this.columns.length} columns)` };

    get rowCount() { return grok_DataFrame_RowCount(this.d); }
    get selection() { return new BitSet(grok_DataFrame_Get_Selection(this.d)); }
    get filter() { return new BitSet(grok_DataFrame_Get_Filter(this.d)); }

    get name() { return grok_DataFrame_Get_Name(this.d); }
    set name(s) { grok_DataFrame_Set_Name(this.d, s); }

    getTag(tag) { return grok_DataFrame_Get_Tag(this.d, tag); }
    setTag(tag, value) { grok_DataFrame_Set_Tag(this.d, tag, value); }

    get columns() { return new ColumnList(grok_DataFrame_Columns(this.d)); }
    get rows() { return new RowList(this, grok_DataFrame_Rows(this.d)); }

    row(i) { return new Row(this, i); }

    get(name, idx) { return this.getCol(name).get(idx); }
    set(name, idx, value) { this.getCol(name).set(idx, value); }

    col(name) { return new Column(grok_DataFrame_ColumnByName(this.d, name)); }

    /// Same as [col], but throws Error if column is not found
    getCol(name) {
        let c = this.col(name);
        if (c === null)
            throw new Error(`No such column: ${name}`);
        return c;
    }

    toCsv() { return grok_DataFrame_ToCsv(this.d); }

    clone(rowMask, columnIds, saveSelection) {
        return new DataFrame(grok_DataFrame_Clone(this.d, rowMask.d, columnIds, saveSelection));
    }

    get currentRow() { return new Row(this, grok_DataFrame_Get_CurrentRowIdx(this.d)); }

    /** Creates a [DataFrame] from list of objects **/
    static fromObjects(list) {
        let table = DataFrame.create(list.length);
        if (list.length === 0)
            return;

        let names = Object.keys(list[0]);
        for (let name of names) {
            let strings = list.map((x) => {
                let value = x[name];
                return value === null ? '' : `${value}`;
            });
            table.columns.add(Column.fromStrings(name, strings));
        }

        return table;
    }

    groupBy(columnNames = []) { return new GroupByBuilder(grok_DataFrame_GroupBy(this.d, columnNames)); }
    
    _event(event, callback) { return _sub(grok_DataFrame_OnEvent(this.d, event, callback)); }

    onValuesChanged(callback) { return this._event('ddt-values-changed', callback); }
    onCurrentRowChanged(callback) { return this._event('ddt-current-row-changed', callback); }
    onMouseOverRowChanged(callback) { return this._event('ddt-mouse-over-row-changed', callback); }
    onCurrentColChanged(callback) { return this._event('ddt-current-col-changed', callback); }
    onMouseOverColChanged(callback) { return this._event('ddt-mouse-over-col-changed', callback); }
    onCurrentCellChanged(callback) { return this._event('ddt-current-cell-changed', callback); }
    onMouseOverRowGroupChanged(callback) { return this._event('ddt-mouse-over-row-group-changed', callback); }
    onNameChanged(callback) { return this._event('ddt-table-name-changed', callback); }
    onMetadataChanged(callback) { return this._event('ddt-table-metadata-changed', callback); }
    onColumnNameChanged(callback) { return this._event('ddt-table-column-name-changed', callback); }
    onColumnSelectionChanged(callback) { return this._event('ddt-column-selection-changed', callback); }
    onColumnsChanged(callback) { return this._event('ddt-columns-changed', callback); }
    onColumnsAdded(callback) { return this._event('ddt-columns-added', callback); }
    onColumnsRemoved(callback) { return this._event('ddt-columns-removed', callback); }
    onRowsRemoved(callback) { return this._event('ddt-rows-removed', callback); }

    onDataChanged(callback) { return _sub(grok_DataFrame_OnDataChanged(this.d, callback)); }
    onFilterChanged(callback) { return _sub(grok_DataFrame_OnFilterChanged(this.d, callback)); }
    onSelectionChanged(callback) { return _sub(grok_DataFrame_OnSelectionChanged(this.d, callback)); }

    fireValuesChanged() { grok_DataFrame_FireValuesChanged(this.d); }
}


const AGG_KEY = "key";      // Special case: to be used in the 'group by' statement.
const AGG_PIVOT = "pivot";  // Special case: to be used as a pivot.
const AGG_FIRST = "first";
const AGG_TOTAL_COUNT = "count";
const AGG_VALUE_COUNT = "values";
const AGG_UNIQUE_COUNT = "unique";
const AGG_MISSING_VALUE_COUNT = "nulls";
const AGG_MIN = "min";
const AGG_MAX = "max";
const AGG_SUM = "sum";
const AGG_MED = "med";
const AGG_AVG = "avg";
const AGG_STDEV = "stdev";
const AGG_VARIANCE = "variance";
const AGG_SKEW = "skew";
const AGG_KURT = "kurt";
const AGG_Q1 = "q1";
const AGG_Q2 = "q2";
const AGG_Q3 = "q3";
const AGG_SELECTED_ROWS_COUNT = "#selected";


const CURRENT_ROW_TO_ROW = 'row to row';
const CURRENT_ROW_TO_SELECTION = 'row to selection';
const CURRENT_ROW_TO_FILTER = 'row to filter';
const MOUSE_OVER_ROW_TO_SELECTION = 'mouse-over to selection';
const MOUSE_OVER_ROW_TO_FILTER = 'mouse-over to filter';
const FILTER_TO_FILTER = 'filter to filter';
const FILTER_TO_SELECTION = 'filter to selection';
const SELECTION_TO_FILTER = 'selection to filter';
const SELECTION_TO_SELECTION = 'selection to selection';


const JOIN_TYPE_INNER = 'inner';
const JOIN_TYPE_OUTER = 'outer';
const JOIN_TYPE_LEFT = 'left';
const JOIN_TYPE_RIGHT = 'right';


class GroupByBuilder {
    constructor(d) { this.d = d; }

    aggregate() { return new DataFrame(grok_GroupByBuilder_Aggregate(this.d)); }

    add(agg, col, resultColName = null) { grok_GroupByBuilder_Add(this.d, agg, col, resultColName); return this; }

    key(col, resultColName = null) { this.add(AGG_KEY, col, resultColName); return this; }
    pivot(col, resultColName = null) { this.add(AGG_PIVOT, col, resultColName); return this; }
    count(col, resultColName = null) { this.add(AGG_TOTAL_COUNT, col, resultColName); return this; }
    uniqueCount(col, resultColName = null) { this.add(AGG_UNIQUE_COUNT, col, resultColName); return this; }
    missingValueCount(col, resultColName = null) { this.add(AGG_MISSING_VALUE_COUNT, col, resultColName); return this; }
    valueCount(col, resultColName = null) { this.add(AGG_VALUE_COUNT, col, resultColName); return this; }
    min(col, resultColName = null) { this.add(AGG_MIN, col, resultColName); return this; }
    max(col, resultColName = null) { this.add(AGG_MAX, col, resultColName); return this; }
    sum(col, resultColName = null) { this.add(AGG_SUM, col, resultColName); return this; }
    med(col, resultColName = null) { this.add(AGG_MED, col, resultColName); return this; }
    avg(col, resultColName = null) { this.add(AGG_AVG, col, resultColName); return this; }
    stdev(col, resultColName = null) { this.add(AGG_STDEV, col, resultColName); return this; }
    variance(col, resultColName = null) { this.add(AGG_VARIANCE, col, resultColName); return this; }
    q1(col, resultColName = null) { this.add(AGG_Q1, col, resultColName); return this; }
    q2(col, resultColName = null) { this.add(AGG_Q2, col, resultColName); return this; }
    q3(col, resultColName = null) { this.add(AGG_Q3, col, resultColName); return this; }
}

/** Central event hub. */
class EventBus {

}


class SemanticValue {
    constructor(d) { this.d = d; }

    static fromValueType(value, semType) { return new SemanticValue(grok_SemanticValue(value, semType)); }

    get value() { return grok_SemanticValue_Get_Value(this.d); }
    set value(x) { grok_SemanticValue_Set_Value(this.d, x); }

    get semType() { return grok_SemanticValue_Get_SemType(this.d); }
    set semType(x) { grok_SemanticValue_Set_SemType(this.d, x); }
}


/** Grok User */
class User {
    constructor(d) { this.d = d; }
    static fromId(id) { return new User(grok_User_From_Id(id)); }
    static current() { return new User(grok_User()); }

    get firstName() { return grok_User_Get_FirstName(this.d); }
    get lastName() { return grok_User_Get_LastName(this.d); }
    get email() { return grok_User_Get_Email(this.d); }
    get picture() { return grok_User_Get_Picture(this.d); }
    get login() { return grok_User_Get_Login(this.d); }

    toMarkup() { return grok_User_ToMarkup(this.d); }
}


class Widget {
    constructor(root) { this.root = root; }

    static react(reactComponent) {
        let widget = new Widget(ui.div());
        ReactDOM.render(reactComponent, widget.root);
        return widget;
    }

    //get root() { return grok_Widget_Get_Root(this.r); }
}


/**
 * A view is typically docked in the main document area of the Grok platform.
 * See [TableView], [SketchView], etc
 */
class View {
    constructor(d) { this.d = d; }

    static fromDart(d) {
        let type = grok_View_Get_Type(d);
        if (type === VIEW_TYPE_TABLE_VIEW)
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

    loadLayout(layout) { return grok_View_Load_Layout(this.d, layout.d);  }
    saveLayout() { return new ViewLayout(grok_View_Save_Layout(this.d)); }

    setRibbonPanels(panels) { grok_View_SetRibbonPanels(this.d, panels); }

    addViewer(viewerType, options = null) {
        let v = new Viewer(grok_View_AddViewer(this.d, viewerType));
        if (options !== null)
            v.options(options);
        return v;
    }

    close() { grok_View_Close(this.d); }
}


class TableView extends View {
    constructor(d) { super(d); }

    get grid() { return new Grid(grok_View_Get_Grid(this.d)); }

    get dataFrame() { return new DataFrame(grok_View_Get_DataFrame(this.d)); }
    set dataFrame(x) { grok_View_Set_DataFrame(this.d, x.d); }

    get toolboxPage() { return new ToolboxPage(grok_View_Get_ToolboxPage(this.d)); }

    histogram      (options = null) { return this.addViewer(VIEWER_HISTOGRAM, options); }
    barChart       (options = null) { return this.addViewer(VIEWER_BAR_CHART, options); }
    boxPlot        (options = null) { return this.addViewer(VIEWER_BOX_PLOT, options); }
    calendar       (options = null) { return this.addViewer(VIEWER_CALENDAR, options); }
    corrPlot       (options = null) { return this.addViewer(VIEWER_CORR_PLOT, options); }
    densityPlot    (options = null) { return this.addViewer(VIEWER_DENSITY_PLOT, options); }
    filters        (options = null) { return this.addViewer(VIEWER_FILTERS, options); }
    form           (options = null) { return this.addViewer(VIEWER_FORM, options); }
    globe          (options = null) { return this.addViewer(VIEWER_GLOBE, options); }
    googleMap      (options = null) { return this.addViewer(VIEWER_GOOGLE_MAP, options); }
    heatMap        (options = null) { return this.addViewer(VIEWER_HEAT_MAP, options); }
    lineChart      (options = null) { return this.addViewer(VIEWER_LINE_CHART, options); }
    shapeMap       (options = null) { return this.addViewer(VIEWER_SHAPE_MAP, options); }
    markup         (options = null) { return this.addViewer(VIEWER_MARKUP, options); }
    matrixPlot     (options = null) { return this.addViewer(VIEWER_MATRIX_PLOT, options); }
    networkDiagram (options = null) { return this.addViewer(VIEWER_NETWORK_DIAGRAM, options); }
    pcPlot         (options = null) { return this.addViewer(VIEWER_PC_PLOT, options); }
    pieChart       (options = null) { return this.addViewer(VIEWER_PIE_CHART, options); }
    scatterPlot    (options = null) { return this.addViewer(VIEWER_SCATTER_PLOT, options); }
    scatterPlot3d  (options = null) { return this.addViewer(VIEWER_SCATTER_PLOT_3D, options); }
    statistics     (options = null) { return this.addViewer(VIEWER_STATISTICS, options); }
    tileViewer     (options = null) { return this.addViewer(VIEWER_TILE_VIEWER, options); }
    treeMap        (options = null) { return this.addViewer(VIEWER_TREE_MAP, options); }
    trellisPlot    (options = null) { return this.addViewer(VIEWER_TRELLIS_PLOT, options); }
    wordCloud      (options = null) { return this.addViewer(VIEWER_WORD_CLOUD, options); }
}


class DataSourceCardView extends View {
    set searchValue(s) { return grok_DataSourceCardView_Set_SearchValue(this.d, s); }
    get permanentFilter() { return grok_DataSourceCardView_Get_PermanentFilter(this.d); }
    set permanentFilter(s) { return grok_DataSourceCardView_Set_PermanentFilter(this.d, s); }
}


class ProjectsView extends DataSourceCardView {
    static create(params) { return new ProjectsView(grok_ProjectsView(params)); }
}


class ViewLayout {
    constructor(d) { this.d = d; }

    static fromJson(json) { return new ViewLayout(grok_ViewLayout_FromJson(json)); }

    toJson() { return grok_ViewLayout_ToJson(this.d); }
}


class Accordion {
    constructor(d) { this.d = d; }
    static create() { return new Accordion(grok_Accordion()); }

    get root() { return grok_Accordion_Get_Root(this.d); }
    get panes() { return grok_Accordion_Get_Panes(this.d).map(p => new AccordionPane(p)); }

    addPane(name, getContent, expanded = false, before = null) {
        return grok_Accordion_AddPane(this.d, name, getContent, expanded, before !== null ? before.d : null);
    }
}


class AccordionPane {
    constructor(d) { this.d = d; }

    get expanded() { return grok_AccordionPane_Get_Expanded(this.d); }
    set expanded(v) { return grok_AccordionPane_Set_Expanded(this.d, v); }

    get name() { return grok_AccordionPane_Get_Name(this.d); }
    set name(name) { return grok_AccordionPane_Set_Name(this.d, name); }
}


class ToolboxPage {
    constructor(d) { this.d = d; }

    get accordion() { return new Accordion(grok_ToolboxPage_Get_Accordion(this.d)); }
}


class VirtualView {
    constructor(d) { this.d = d; }

    static create() { return new VirtualView(grok_VirtualItemView()); }

    get root() { return grok_VirtualItemView_Get_Root(this.d); }
    setData(length, renderer) { grok_VirtualItemView_SetData(this.d, length, renderer); }
}


class Dialog {
    constructor(d) { this.d = d; }
    static create(title = '') { return new Dialog(grok_Dialog(title)); }

    onOK(handler) { grok_Dialog_OnOK(this.d, handler); return this; }

    show() { grok_Dialog_Show(this.d); return this; }

    add(x) { grok_Dialog_Add(this.d, x); return this; }
}


class Menu {
    constructor(d) { this.d  = d; }

    static popup() { return new Menu(grok_Menu()); }

    group(s) { return new Menu(grok_Menu_Group(this.d, s)); }
    item(item, onClick) { return new Menu(grok_Menu_Item(this.d, item, onClick)); }
    items(items, onClick) { return new Menu(grok_Menu_Items(this.d, items, onClick)); }
    separator() { return new Menu(grok_Menu_Separator(this.d)); }
    show() { return new Menu(grok_Menu_Show(this.d)); }
}


/** A viewer that is typically docked inside a [TableView]. */
class Viewer {
    constructor(d) { this.d = d; }

    static fromType(viewerType, table, options = null) {
        return new Viewer(grok_Viewer_FromType(viewerType, table.d, _toJson(options)));
    }

    options(map) { grok_Viewer_Options(this.d, JSON.stringify(map)); }
    close() { grok_Viewer_Close(this.d); }
    serialize()  { return grok_Viewer_Serialize(this.d); }

    get root() { return grok_Viewer_Root(this.d); }

    static grid        (t, options = null) { return new Viewer(grok_Viewer_Grid(t.d, _toJson(options))); }
    static histogram   (t, options = null) { return new Viewer(grok_Viewer_Histogram(t.d, _toJson(options))); }
    static barChart    (t, options = null) { return Viewer.fromType(VIEWER_BAR_CHART, t, options);  }
    static boxPlot     (t, options = null) { return new Viewer(grok_Viewer_BoxPlot(t.d, _toJson(options)));  }
    static filters     (t, options = null) { return new Viewer(grok_Viewer_Filters(t.d, _toJson(options))); }
    static scatterPlot (t, options = null) { return new Viewer(grok_Viewer_ScatterPlot(t.d, _toJson(options))); }
}


class JsLookAndFeel {
    constructor() {
        return new Proxy(this, {
            set(target, name, value) {
                target.table.set(name, target.idx, value);
                return true;
            },
            get(target, name) {
                return target.table.get(name, target.idx);
            }
        });
    }

    get(name) { return null; }
}


/** JavaScript implementation of the grok viewer */
class JsViewer {
    constructor() {
        this.root = ui.div();
        this.properties = [];
        this.dataFrameHandle = null;
        this.dataFrame = null;
    }

    onFrameAttached(dataFrameHandle) {}

    onPropertyChanged(property) {}

    getProperties() { return this.properties; }
    getDartProperties() {
        return this.getProperties().map((p) => p.d);
    }

    /** cleanup() will get called when the viewer is disposed **/
    registerCleanup(cleanup) { grok_Widget_RegisterCleanup(this.root, cleanup); }

    _prop(name, type, value) {
        let obj = this;
        let p = Property.create(name, type, () => obj[name], null, value);
        p.set = function(_, x) {
              //console.log("xo");
              obj[name] = x;
              obj.onPropertyChanged(p);
            };

        this.properties.push(p);
        return p.defaultValue;
    }

    int(name, value) { return this._prop(name, TYPE_INT, value); }
    float(name, value) { return this._prop(name, TYPE_FLOAT, value); }
    string(name, value) { return this._prop(name, TYPE_STRING, value); }
    bool(name, value) { return this._prop(name, TYPE_BOOL, value); }
    dateTime(name, value) { return this._prop(name, TYPE_DATE_TIME, value); }
}


class JsEntityMeta {
    get type() { throw 'Not defined.'; }

    isApplicable(x) { throw 'Not defined.'; }

    /** String representation of the [item], by default item.toString(). */
    getCaption(x) { return `${x}`; };

    renderIcon(x) { return ui.divText(this.getCaption(x)); }

    renderMarkup(x) { return ui.divText(this.getCaption(x)); }
    renderTooltip(x) { return ui.divText(this.getCaption(x)); }
    renderCard(x) { return ui.divText(this.getCaption(x)); }
    renderProperties(x) { return ui.divText(this.getCaption(x)); }
    renderView(x) {return this.renderProperties(x); }

    /** Gets called once upon the registration of meta class. */
    init() {}

    static register(meta) { grok_Meta_Register(meta); }

    registerParamFunc(name, run) {
        grok.functions.registerParamFunc(name, this.type, run, this.isApplicable);
    }
}


/** Balloon-style visual notifications. */
class Balloon {
    /** Shows information message (green background) */
    info(s) { grok_Balloon(s, 'info'); }

    /** Shows information message (red background) */
    error(s) { grok_Balloon(s, 'error'); }
}


/** Subscription to an event stream. Call [cancel] to stop listening. */
class StreamSubscription {
    constructor(d) { this.d = d; }

    cancel() { grok_Subscription_Cancel(d); }
}


/** Grok functions */
class Functions {
    register(func) { grok_RegisterFunc(func); }
    registerParamFunc(name, type, run, check = null, description = null) { grok_RegisterParamFunc(name, type, run, check, description); }
}


class Entity {
    constructor(d) { this.d = d;}

    get id() { return grok_Entity_Get_Id([this.d]); }
    set id(x) { return grok_Entity_Set_Id([this.d], x); }
    get name() { return grok_Entity_Get_Name([this.d]); }
    set name(x) { return grok_Entity_Set_Name([this.d], x); }
    get path() { return grok_Entity_Path([this.d]); }
}

/** Represents a project */
class Project extends Entity {
    constructor(d) { super(d); }

    get description() { return grok_Project_Description(this.d); }
    get isDirty() { return grok_Project_IsDirty(this.d); }
    get isEmpty() { return grok_Project_IsEmpty(this.d); }

    toMarkup() { return grok_Project_ToMarkup(this.d); }
}


/** Represents a data query */
class DataQuery  extends Entity {
    constructor(d) { super(d); }

    get query() { return grok_DataQuery_Query([this.d]); }
}

class GrokPackage {
    constructor(webRoot) {
        this.webRoot = webRoot;
    }

    // Override init() method to provide package-specific initialization.
    // It is guaranteed to get called exactly once before the execution of any function below.
    /*async*/ init() { return Promise.resolve(null); }
}


class GridCell {
    constructor(d) { this.d = d; }

    get cellType() { return grok_GridCell_Get_CellType(this.d); }

    get isTableCell() { return grok_GridCell_Get_IsTableCell(this.d); }
    get isRowHeader() { return grok_GridCell_Get_IsRowHeader(this.d); }
    get isColHeader() { return grok_GridCell_Get_IsColHeader(this.d); }

    get tableColumn() { return this.gridColumn.column; }
    get tableRow() { return this.isTableCell ? this.cell.row : null; }
    get tableRowIndex() { return this.isTableCell ? this.cell.rowIndex : null; }

    get gridRow() { return grok_GridCell_Get_GridRow(this.d); }
    get gridColumn() { return new GridColumn(grok_GridCell_Get_GridColumn(this.d)); }

    get customText() { return grok_GridCell_Get_CustomText(this.d); }
    set customText(x) { return grok_GridCell_Set_CustomText(this.d, x); }

    get grid() { return new Grid(grok_GridCell_Get_Grid(this.d)); }
    get cell() {
        return new Cell(grok_GridCell_Get_Cell(this.d));
    }
    get style() { return new GridCellStyle(grok_GridCell_Get_Style(this.d)); }
}


class GridColumn {
    constructor(d) { this.d = d; }

    get column() {
        let col = grok_GridColumn_Get_Column(this.d);
        return col === null ? null : new Column(col);
    }

    get idx() { return grok_GridColumn_Get_Idx(this.d); }

    get name() { return grok_GridColumn_Get_Name(this.d); }
    set name(x) { return grok_GridColumn_Set_Name(this.d, x); }

    get width() { return grok_GridColumn_Get_Width(this.d); }
    set width(x) { return grok_GridColumn_Set_Width(this.d, x); }

    get backColor() { return grok_GridColumn_Get_BackColor(this.d); }
    set backColor(x) { return grok_GridColumn_Set_BackColor(this.d, x); }

    get format() { return grok_GridColumn_Get_Format(this.d); }
    set format(x) { return grok_GridColumn_Set_Format(this.d, x); }

    get cellType() { return grok_GridColumn_Get_CellType(this.d); }
    set cellType(x) { return grok_GridColumn_Set_CellType(this.d, x); }

    get visible() { return grok_GridColumn_Get_Visible(this.d); }
    set visible(x) { return grok_GridColumn_Set_Visible(this.d, x); }

    get categoryColors() { return grok_GridColumn_Get_CategoryColors(this.d); }
    set categoryColors(x) { return grok_GridColumn_Set_CategoryColors(this.d, x); }

    get editable() { return grok_GridColumn_Get_Editable(this.d); }
    set editable(x) { return grok_GridColumn_Set_Editable(this.d, x); }

    get selected() { return grok_GridColumn_Get_Selected(this.d); }
    set selected(x) { return grok_GridColumn_Set_Selected(this.d, x); }
}


class GridColumnList {
    constructor(d) { this.d = d; }

    byName(columnName) { return new GridColumn(grok_GridColumnList_ByName(this.d, columnName)); }

    setOrder(columnNames) { grok_GridColumnList_SetOrder(this.d, columnNames); }
    setVisible(columnNames) { grok_GridColumnList_SetVisible(this.d, columnNames); }
}


class Grid extends Viewer {
    constructor(d) { super(d); }

    /** Sorts data by [columnIds].
     *  Specify sort directions via [asc] array (true = ascending, false = descending)
     *  If [asc] is not specified, sorts in ascending order. */
    sort(columns, orders = null) { grok_Grid_Sort(this.d, columns.map((c) => c instanceof Column ? c.d : c), orders); }

    get columns() { return new GridColumnList(grok_Grid_Get_Columns(this.d)); }

    col(name) { return this.columns.byName(name); }

    onCellPrepare(callback) {
        return _sub(grok_Grid_OnCellPrepare(this.d, (dcell) => { return callback(new GridCell(dcell)); }));
    }

    onCellTooltip(callback) {
        return _sub(grok_Grid_OnCellTooltip(this.d, (dcell, x, y) => { return callback(new GridCell(dcell), x, y); }));
    }

    static create(table) { return new Grid(grok_Grid_Create(table.d)); }
}


class GridCellStyle {
    constructor(d) { this.d = d; }

    get font() { return grok_GridCellStyle_Get_Font(this.d); }
    set font(x) { return grok_GridCellStyle_Set_Font(this.d, x); }

    get textColor() { return grok_GridCellStyle_Get_TextColor(this.d); }
    set textColor(x) { return grok_GridCellStyle_Set_TextColor(this.d, x); }

    get backColor() { return grok_GridCellStyle_Get_BackColor(this.d); }
    set backColor(x) { return grok_GridCellStyle_Set_BackColor(this.d, x); }

    get element() { return grok_GridCellStyle_Get_Element(this.d); }
    set element(x) { return grok_GridCellStyle_Set_Element(this.d, x); }
}


/// Input control for Property.
class InputBase {
    constructor(d) { this.d = d; }

    get root() { return grok_InputBase_Get_Root(this.d); };
    get caption() { return grok_InputBase_Get_Caption(this.d); };
    get format() { return grok_InputBase_Get_Format(this.d); } ;
    get captionLabel() { return grok_InputBase_Get_CaptionLabel(this.d); };
    get input() { return grok_InputBase_Get_Input(this.d); };

    get nullable() { return grok_InputBase_Get_Nullable(this.d); };
    set nullable(v) { return grok_InputBase_Set_Nullable(this.d, v); };

    get value() { return grok_InputBase_Get_Value(this.d); };
    set value(x) { return grok_InputBase_Set_Value(this.d, x); };

    get stringValue() { return grok_InputBase_Get_StringValue(this.d); };
    set stringValue(s) { return grok_InputBase_Set_StringValue(this.d, s); };

    get readOnly() { return grok_InputBase_Get_ReadOnly(this.d); };
    set readOnly(v) { return grok_InputBase_Set_ReadOnly(this.d, v); };

    get enabled() { return grok_InputBase_Get_Enabled(this.d); };
    set enabled(v) { return grok_InputBase_Set_Enabled(this.d, v); };

    /// Occurs when [value] is changed, either by user or programmatically.
    onChanged(callback) { return _sub(grok_InputBase_OnChanged(this.d, callback)); }

    /// Occurs when [value] is changed by user.
    onInput(callback) { return _sub(grok_InputBase_OnInput(this.d, callback)); }

    save() { return grok_InputBase_Save(this.d); };
    load(s) { return grok_InputBase_Load(this.d, s); };

    init() { return grok_InputBase_Init(this.d); };
    fireChanged() { return grok_InputBase_FireChanged(this.d); };
    addCaption(caption) { grok_InputBase_AddCaption(this.d, caption); };
    addPatternMenu(pattern) { grok_InputBase_AddPatternMenu(this.d, pattern); }
    setTooltip(msg) { grok_InputBase_SetTooltip(this.d, msg); };

    static forProperty(property) { return new InputBase(grok_InputBase_ForProperty(property.d)); }
}


class ProgressIndicator {
    constructor(d) { this.d = d; }

    get description() { return grok_ProgressIndicator_Get_Description(this.d); }
    set description(s) { grok_ProgressIndicator_Set_Description(this.d, s); }

    update(percent) { grok_ProgressIndicator_Update(percent); }
}


class TagEditor {
    constructor(d) { this.d = d; }

    static create() { return new TagEditor(grok_TagEditor()); }

    get root() { return grok_TagEditor_Get_Root(this.d); }

    get tags() { return grok_TagEditor_Get_Tags(this.d); }

    addTag(tag, notify = true) { return grok_TagEditor_AddTag(this.d, tag, notify); }
    removeTag(tag) { grok_TagEditor_RemoveTag(this.d, tag); }
    clearTags() { grok_TagEditor_ClearTags(this.d); }

    set acceptsDragDrop(predicate) { grok_TagEditor_Set_AcceptsDragDrop(this.d, (x) => predicate(_wrap(x, false))); };
    set doDrop(action) { grok_TagEditor_Set_DoDrop(this.d, (x) => action(_wrap(x, false))); }

    onChanged(callback) { return _sub(grok_TagEditor_OnChanged(this.d, callback)); }
}


class TagElement {
    constructor(d) { this.d = d; }

    get tag() { return grok_TagElement_Get_Tag(this.d); };
    set tag(x) { return grok_TagElement_Set_Tag(this.d, x); };
}


class DateTime {
    constructor(d) { this.d = d; }

    static fromDate(date) { return DateTime.fromMillisecondsSinceEpoch(date.getTime()); }
    static fromMillisecondsSinceEpoch(millisecondsSinceEpoch) {
        return new DateTime(grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
    }
}


class FuncCall {
    constructor(d) { this.d = d; }

    getParamValue(name) { return grok_FuncCall_Get_Param_Value(this.d, name); }
}


/** Grok entry point, use it to get access to top-level views, tables, methods, etc. */
class Grok {

    /** Current table */
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

    get balloon() { return new Balloon(); }

    get dapi() { return new Dapi(); }

    /** Current user */
    get user() { return new User(grok_User()); }

    /** Presentaion mode **/
    get presentationMode() { return grok_Get_PresentationMode(); }
    set presentationMode(x) { return grok_Set_PresentationMode(x); }

    /** Hide dock tabs in presentaion mode **/
    get hideTabsInPresentationMode() { return grok_Get_HideTabsInPresentationMode(); }
    set hideTabsInPresentationMode(x) { return grok_Set_HideTabsInPresentationMode(x); }

    get functions() { return new Functions(); }

    dockElement(e, title = null, dockStyle = null, ratio = null) { grok_DockElement(e, title, dockStyle, ratio); }

    tableByName(s) { return new DataFrame(grok_TableByName(s)); }

    /**
     * Creates a generic dataset with the defined number of rows and columns.
     * [dataset] allowed values:
     * "wells" - experimental plate wells: barcode, row, col, pos, volume, concentration, role
     * "demog" - clinical study demographics data: subj, study, site, sex, race, disease, start date
     * "biosensor": wearable sensor data: time, x, y, z, temp, eda
     * "random walk": random walk data for the specified number of dimensions
     * **/
    testData(dataset, rows = 10000, columns = 3) { return new DataFrame(grok_TestData(dataset, rows, columns)); }
    getDemoTable(path) {
        return new Promise((resolve, reject) => grok_GetDemoTable(path, (t) => resolve(_wrap(t))));
    }

    newView(name = 'view', children = []) {
        let view = View.create();
        view.name = name;
        ui.appendAll(view.root, children);
        this.addView(view);
        return view;
    }

    getVar(name) {return _toJs(grok_GetVar(name)); };
    setVar(name, value) {return grok_SetVar(name, _toDart(value)); };

    addView(v) { grok_AddView(v.d); }
    addTableView(t) { return new TableView(grok_AddTableView(t.d)); }
    getTableView(name) { return new TableView(grok_GetTableView(name)); }
    closeAll() { grok_CloseAll(); }
    route(url) { return View.fromDart(grok_Route(url)); }

    get topMenu() { return new Menu(grok_Get_TopMenu()); }

    parseCsv(csv) { return new DataFrame(grok_ParseCsv(csv)); }

    script(s) {
        return new Promise((resolve, reject) => grok_Script(s, (t) => resolve(_wrap(t, false))));
    }

    loadDataFrame(csvUrl) {
        return new Promise((resolve, reject) => {
            fetch(csvUrl)
                .then((r) => r.text().then((csv) => resolve(this.parseCsv(csv))));
        });
        // let csv = await r.text();
        // return this.parseCsv(csv);
    }

    linkTables(t1, t2, keyColumns1, keyColumns2, linkTypes) { grok_LinkTables(t1.d, t2.d, keyColumns1, keyColumns2, linkTypes); };
    compareTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2) { grok_CompareTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2); };
    joinTables(t1, t2, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace) {
        return new DataFrame(grok_JoinTables(t1.d, t2.d, keyColumns1, keyColumns2, valueColumns1, valueColumns2, joinType, inPlace));
    };

    scriptSync(s) { return _wrap(grok_ScriptSync(s), false); }

    openTable(id) { return new Promise((resolve, reject) => grok_OpenTable(id, (t) => resolve(new DataFrame(t)))); }
    query(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_Query(queryName, queryParameters, adHoc, pollingInterval, (t) => resolve(new DataFrame(t))));
    }
    callQuery(queryName, queryParameters = null, adHoc = false, pollingInterval = 1000) {
        return new Promise((resolve, reject) => grok_CallQuery(queryName, queryParameters, adHoc, pollingInterval, (c) => resolve(new FuncCall(c))));
    }
    callFunc(name, parameters = {}, showProgress = false) { return new Promise((resolve, reject) => grok_CallFunc(name, parameters, (out) => resolve(out), showProgress)); }
    evalFunc(name) { return new Promise((resolve, reject) => grok_EvalFunc(name, (out) => resolve(out))); }

    registerFunc(func) { grok_RegisterFunc(func); }

    registerViewer(name, description, createViewer) { grok_RegisterViewer(name, description, createViewer); }

    get currentObject() { return _wrap(grok_Get_CurrentObject(), false); }
    set currentObject(x) { grok_Set_CurrentObject(_toDart(x)); }

    detectSemanticTypes(t) { return new Promise((resolve, reject) => grok_DetectSematicTypes(t.d, (_) => resolve())); }

    onEvent(event, f) { return grok_OnEvent(event, f); };
    onCurrentViewChanged(f) { return grok_OnEvent('d4-current-view-changed', f); }
    onCurrentCellChanged(f) { return grok_OnEvent('d4-current-cell-changed', f); }
    onTableAdded(f) { return grok_OnEvent('d4-table-added', f); }
    onTableRemoved(f) { return grok_OnEvent('d4-table-removed', f); }
    onQueryStarted(f) { return grok_OnEvent('d4-query-started', f); }
    onQueryFinished(f) { return grok_OnEvent('d4-query-finished', f); }

    onViewChanged(f) { return grok_OnEvent('grok-view-changed', (v) => f(View.fromDart(v))); }
    onViewAdded(f) { return grok_OnEvent('grok-view-added', (v) => f(View.fromDart(v))); }
    onViewRemoved(f) { return grok_OnEvent('grok-view-removed', (v) => f(View.fromDart(v))); }
    onViewRenamed(f) { return grok_OnEvent('grok-view-renamed', (v) => f(View.fromDart(v))); }

    onCurrentProjectChanged(f) { return grok_OnEvent('grok-current-project-changed', (a) => f(new Project(a))); }
    onProjectUploaded(f) { return grok_OnEvent('grok-project-uploaded', (a) => f(new Project(a))); }
    onProjectSaved(f) { return grok_OnEvent('grok-project-saved', (a) => f(new Project(a))); }
    onProjectOpened(f) { return grok_OnEvent('grok-project-opened', (a) => f(new Project(a))); }
    onProjectClosed(f) { return grok_OnEvent('grok-project-closed', (a) => f(new Project(a))); }
    onProjectModified(f) { return grok_OnEvent('grok-project-modified', (a) => f(new Project(a))); }
}

const TYPE_INT = 'int';
const TYPE_BIG_INT = 'bigint';
const TYPE_FLOAT = 'double';
const TYPE_NUM = 'num';
const TYPE_BOOL = 'bool';
const TYPE_STRING = 'string';
const TYPE_STRING_LIST = 'string_list';
const TYPE_DATE_TIME = 'datetime';
const TYPE_OBJECT = 'object';
const TYPE_PROJECT = 'project';
const TYPE_DATA_FRAME = 'dataframe';
const TYPE_DATA_FRAME_LIST = 'dataframe_list';
const TYPE_CELL = 'cell';
const TYPE_COLUMN = 'column';
const TYPE_COLUMN_LIST = 'column_list';
const TYPE_GRAPHICS = 'graphics';
const TYPE_ROW_FILTER = 'tablerowfiltercall';
const TYPE_COLUMN_FILTER = 'colfiltercall';
const TYPE_BIT_SET = 'bitset';
const TYPE_MAP = 'map';
const TYPE_DYNAMIC = 'dynamic';
const TYPE_VIEWER = 'viewer';  // [ViewerBase] subclasses
const TYPE_LIST = 'list';
const TYPE_SEM_VALUE = 'semantic_value';
const TYPE_FUNC_CALL = 'funccall';
const TYPE_PROPERTY = 'property';
const TYPE_CATEGORICAL = 'categorical';
const TYPE_NUMERICAL = 'numerical';

const TYPE_USER = 'User';
const TYPE_MENU = 'Menu';
const TYPE_SEMANTIC_VALUE = 'semantic_value';

const TYPES_SCALAR = new Set([TYPE_INT, TYPE_BIG_INT, TYPE_FLOAT, TYPE_NUM, TYPE_BOOL, TYPE_STRING]);

const VIEW_TYPE_TABLE_VIEW = 'TableView';


/////// Stats

const STATS_TOTAL_COUNT = "count";
const STATS_VALUE_COUNT = "values";
const STATS_UNIQUE_COUNT = "unique";
const STATS_MISSING_VALUE_COUNT = "nulls";
const STATS_MIN = "min";
const STATS_MAX = "max";
const STATS_SUM = "sum";
const STATS_MED = "med";
const STATS_AVG = "avg";
const STATS_STDEV = "stdev";
const STATS_VARIANCE = "variance";
const STATS_SKEW = "skew";
const STATS_KURT = "kurt";
const STATS_Q1 = "q1";
const STATS_Q2 = "q2";
const STATS_Q3 = "q3";

/////// Tags

const TAGS_DESCRIPTION = 'description';

const TAGS_TOOLTIP = '.tooltip';

/** JSON-encoded list of strings to be used in a cell editor. Applicable for string columns only. */
const TAGS_CHOICES = '.choices';

/** When set to 'true', switches the cell editor to a combo box that only allows to choose values
 from a list of already existing values in the column.
 Applicable for string columns only.
 See also [TAGS_CHOICES]. */
const TAGS_AUTO_CHOICES = '.auto-choices';

////// Viewers
const VIEWER_HISTOGRAM = 'Histogram';
const VIEWER_BAR_CHART = 'Bar chart';
const VIEWER_BOX_PLOT = 'Box plot';
const VIEWER_CALENDAR = 'Calendar';
const VIEWER_CORR_PLOT = 'Correlation plot';
const VIEWER_DENSITY_PLOT = 'Density plot';
const VIEWER_FILTERS = 'Filters';
const VIEWER_FORM = 'Form';
const VIEWER_GLOBE = 'Globe';
const VIEWER_GRID = 'Grid';
const VIEWER_GOOGLE_MAP = 'Google map';
const VIEWER_HEAT_MAP = 'Heat map';
const VIEWER_LINE_CHART = 'Line chart';
const VIEWER_SHAPE_MAP = 'Shape map';
const VIEWER_MARKUP = 'Markup';
const VIEWER_MATRIX_PLOT = 'Matrix plot';
const VIEWER_NETWORK_DIAGRAM = 'Network diagram';
const VIEWER_PC_PLOT = 'PC Plot';
const VIEWER_PIE_CHART = 'Pie chart';
const VIEWER_SCATTER_PLOT = 'scatter plot';
const VIEWER_SCATTER_PLOT_3D = '3d scatter plot';
const VIEWER_STATISTICS = 'Statistics';
const VIEWER_TILE_VIEWER = 'Tile viewer';
const VIEWER_TREE_MAP = 'Tree map';
const VIEWER_TRELLIS_PLOT = 'Trellis plot';
const VIEWER_WORD_CLOUD = 'Word cloud';

function _toJs(d) { return _wrap(d, false); }

/** Instantiates the corresponding JS handler for the Dart object [d]. */
function _wrap(d, check = true) {
    let type = grok_GetType(d);
    if (type === TYPE_DATA_FRAME) return new DataFrame(d);
    if (type === TYPE_COLUMN) return new Column(d);
    if (type === TYPE_PROPERTY) return new Property(d);
    if (type === TYPE_PROJECT) return new Project(d);
    if (type === TYPE_USER) return new User(d);
    if (type === TYPE_SEMANTIC_VALUE) return new SemanticValue(d);
    if (type === TYPE_MENU) return new Menu(d);

    if (check)
        throw `Not supported type: ${type}`;
    
    return d;
}


function _toDart(x) {
    return (typeof x.d !== 'undefined') ? x.d : x;
}


function _toJson(x) {
    return x === null ? null : JSON.stringify(x);
}

function _jsThen(promise, f) {
    promise.then(f);
}

function paramsToJs(params) {
    let result = [];
    for (let i = 0; i < params.length; i++) {
        let type = grok_GetType(params[i]);
        if (type !== null && !TYPES_SCALAR.has(type))
            result.push(_wrap(params[i]));
        else
            result.push(params[i]);
    }

    return result;
}

function callFuncWithDartParameters(f, params) {
    let jsParams = paramsToJs(params);
    return f.apply(null, jsParams);
}

function _sub(d) { return new StreamSubscription(d); }

/** Obsolete and will be removed soon; use "grok". */
gr = new Grok();

grok = new Grok();

window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return grok_Error(message, url, lineNumber, columnNumber, errorObject);
};

let time = function(s, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${s}: ${stop - start}ms`);
    grok.balloon.info(`${s}: ${stop - start}ms`);
    return result;
};

class ui {

    static e (s, cl = null) {
        let x = document.createElement(s);
        if (cl !== null)
            ui._class(x, cl);
        return x;
    }

    static appendAll(root, elements) {
        let fragment = document.createDocumentFragment();
        for (let i = 0; i < elements.length; i++) {
            var e = elements[i];
            if (e instanceof Viewer)
               e = e.root;
            fragment.appendChild(e);
        }
        root.appendChild(fragment);
        return root;
    }

    static _innerText(x, s) { x.innerText = s; return x; }
    static _class(x, s) { x.classList.add(s); return x; }
    static _color(x, s) { x.style.color = s; return x; }
    static _backColor(x, s) { x.style.backgroundColor = s; return x; }

    static canvas() { return document.createElement("CANVAS"); }

    static h1(s) { return ui._innerText(ui.e('h1'), s); }
    static h2(s) { let x = ui.e('h2'); x.innerText = s; return x; }
    static h3(s) { let x = ui.e('h3'); x.innerText = s; return x; }

    static accordion() { return  Accordion.create(); }

    static divText(s) {
        let e = document.createElement('div');
        e.innerText = s;
        return e;
    }

    static iconFA(name, handler, tooltip = null) {
        let i = document.createElement('i');
        i.classList.add('grok-icon');
        i.classList.add('fal');
        i.classList.add(`fa-${name}`);
        i.addEventListener("click", handler);
        if (tooltip !== null)
            ui.tooltip(i, tooltip);
        return i;
    }

    static render(x) { return grok_UI_Render(x); }
    static renderCard(x) { return grok_UI_RenderCard(x); }
    static span(x) { return grok_UI_Span(x); }

    static div(items = [], className = null) { return grok_UI_Div(items, className); }
    static divV(items, className = null) { return grok_UI_DivV(items, className); }
    static divH(items, className = null) { return grok_UI_DivH(items, className); }

    static card(content) { return ui.div([content], 'd4-item-card'); }

    static loader() { return grok_UI_Loader(); }
    static setUpdateIndicator(element, updating = true) { return grok_UI_SetUpdateIndicator(element, updating)};

    static button(text, handler, tooltip = null) { return grok_UI_Button(text, handler, tooltip); }
    static bigButton(text, handler, tooltip = null) { return grok_UI_BigButton(text, handler, tooltip); }

    /** Creates a visual table based on [map]. */
    static tableFromMap(map) { return grok_UI_TableFromMap(map); }

    /** Creates a visual element representing list of [items]. */
    static list(items) { return grok_UI_List(Array.from(items).map(_toDart)); }

    /** Creates a [Dialog]. */
    static dialog(title = '') { return Dialog.create(title); }

    /** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
     * tooltip, and popup menu.
     * Returns [element]. */
    static bind(item, element) { return grok_UI_Bind(item, element); }

    static virtualView(length, renderer) {
        let view = VirtualView.create();
        view.setData(length, renderer);
        return view.root;
    }

    static svgMol(smiles, width = 300, height = 200) {
        let m = OCL.Molecule.fromSmiles(smiles);
        let root = document.createElement('div');
        root.innerHTML = m.toSVG(width, height);
        return root;
    }

    static popupMenu(items) {
        function populate(menu, item) {
            for (let key of Object.keys(item)) {
                let value = item[key];
                if (value instanceof Function)
                    menu.item(key, value);
                else
                    populate(menu.group(key), value);
            }
        }

        let menu = Menu.popup();
        populate(menu, items);
        menu.show();
    }

    static tooltipHide() { grok_Tooltip_Hide(); }

    static tooltip(e, x) { grok_Tooltip_SetOn(e, x); return e; }

    static tooltipShow(content, x, y) { grok_Tooltip_Show(content, x, y); }

    static inputs(inputs) { return ui.div(inputs.map((x) => x.root), 'pure-form,pure-form-aligned');}

    static intInput(name, value) { return new InputBase(grok_IntInput(name, value)); }
    static choiceInput(name, selected, items) { return new InputBase(grok_ChoiceInput(name, selected, items)); }
    static multiChoiceInput(name, value, items) { return new InputBase(grok_MultiChoiceInput(name, value, items)); }
    static stringInput(name, value) { return new InputBase(grok_StringInput(name, value)); }
    static floatInput(name, value) { return new InputBase(grok_FloatInput(name, value)); }
    static dateInput(name, value) { return new InputBase(grok_DateInput(name, value.d)); }
    static boolInput(name, value) { return new InputBase(grok_BoolInput(name, value)); }
    static moleculeInput(name, value) { return new InputBase(grok_MoleculeInput(name, value)); }
}


class Dapi {
    constructor() {}

    getEntities(ids) {
        return new Promise((resolve, reject) => grok_Dapi_Entities_GetEntities(ids, (q) => {
            return resolve(q[0].map(_wrap));
        }));
    }

    get queries() { return new HttpDataSource(grok_Dapi_Queries(), (a) => new DataQuery(a)); }
    get connections() { return new HttpDataSource(grok_Dapi_Connections(), (a) => new Entity(a)); }
    get jobs() { return new HttpDataSource(grok_Dapi_Jobs(), (a) => new Entity(a)); }
    get notebooks() { return new HttpDataSource(grok_Dapi_Notebooks(), (a) => new Entity(a)); }
    get models() { return new HttpDataSource(grok_Dapi_Models(), (a) => new Entity(a)); }
    get packages() { return new HttpDataSource(grok_Dapi_Packages(), (a) => new Entity(a)); }
    get layouts() { return new HttpDataSource(grok_Dapi_Layouts(), (a) => new ViewLayout(a)); }
    get tables() { return new HttpDataSource(grok_Dapi_Tables(), (a) => new Entity(a)); }
    get users() { return new UsersDataSource(grok_Dapi_Users(), (a) => new User(a)); }
    get groups() { return new HttpDataSource(grok_Dapi_Groups(), (a) => new Entity(a)); }
    get scripts() { return new HttpDataSource(grok_Dapi_Scripts(), (a) => new Entity(a)); }
    get projects() { return new HttpDataSource(grok_Dapi_Projects(), (a) => new Project(a)); }
    get userDataStorage() { return new UserDataStorage(); }
}


class HttpDataSource {
    constructor(s, instance) {
        this.s = s;
        this.instance = instance;
    }

    list() {
        let s = this.instance;
        return new Promise((resolve, reject) => grok_DataSource_List(this.s, (q) => resolve(q.map(s))));
    }

    find(id) {
        let s = this.instance;
        return new Promise((resolve, reject) => grok_DataSource_Find(this.s, id, (q) => resolve(s(q[0]))));
    }

    by(i) {
        this.s = grok_DataSource_By(this.s, i);
        return this;
    }

    page(i) {
        this.s = grok_DataSource_Page(this.s, i);
        return this;
    }

    nextPage() {
        this.s = grok_DataSource_NextPage(this.s);
        return this;
    }

    filter(w) {
        this.s = grok_DataSource_WhereSmart(this.s, w);
        return this;
    }

    order(name, desc = false) {
        this.s = grok_DataSource_Order(this.s, name, desc);
        return this;
    }
}


class UsersDataSource extends HttpDataSource {
    constructor(s, instance) {
        super(s, instance);
    }

    current() {
        let s = this.instance;
        return new Promise((resolve, reject) => grok_UsersDataSource_Current(this.s, (q) => resolve(s(q[0]))));
    }
}


class UserDataStorage {
    constructor() {}

    postValue(name, key, value, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_PostValue(name, key, value, currentUser, () => resolve()));
    }

    post(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Post(name, data, currentUser, () => resolve()));
    }

    put(name, data, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Put(name, data, currentUser, () => resolve()));
    }

    get(name, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Get(name, currentUser, (data) => resolve(data)));
    }

    getValue(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_GetValue(name, key, currentUser, (value) => resolve(value)));
    }

    remove(name, key, currentUser = true) {
        return new Promise((resolve, reject) =>
            grok_Dapi_UserDataStorage_Delete(name, key, currentUser, () => resolve()));
    }
}


/** UI Tools **/
class uit {

    static handleResize(element, onChanged) {
       var width = element.clientWidth;
       var height = element.clientHeight;
       let interval = setInterval(() => {
           console.log('.');
           let newWidth = element.clientWidth;
           let newHeight = element.clientHeight;
           if (newWidth !== width || newHeight !== height) {
               width = newWidth;
               height = newHeight;
               onChanged(width, height);
           }
       }, 100);
       return () => clearInterval(interval);
    }
}


class Color {

    static r(c) { return (c >> 16) & 0xFF; }
    static g(c) { return (c >> 8) & 0xFF; }
    static b(c) { return c & 0xFF; }

    /** Returns i-th categorical color (looping over the palette if needed) */
    static getCategoricalColor(i) { return Color.categoricalPalette[i % Color.categoricalPalette.length]; }

    /** Returns either black or white color, depending on which one would be most contrast to the specified [color] */
    static getContrastColor(color) { return grok_Color_GetContrastColor(color); }

    static toRgb(color) { return color === null ? '': `rgb(${Color.r(color)},${Color.g(color)},${Color.b(color)})`; }

    static get categoricalPalette() { return grok_Color_CategoricalPalette(); }

    static scale(x, min, max) {
        return min === max ? min : (x - min) / (max - min);
    }

    static get gray() { return 0xFF808080; }
    static get lightLightGray() { return 0xFFF0F0F0; }
    static get lightGray() { return 0xFFD3D3D3; }
    static get darkGray() { return 0xFF838383; }
    static get blue() { return 0xFF0000FF; }
    static get green() { return 0xFF00FF00; }
    static get darkGreen() { return 0xFF006400; }
    static get black() { return 0xFF000000; }
    static get yellow() { return 0xFFFFFF00; }
    static get white() { return 0xFFFFFFFF; }
    static get red() { return 0xFFFF0000; }
    static get darkRed() { return 0xFF8b0000; }
    static get maroon() { return 0xFF800000; }
    static get olive() { return 0xFF808000; }
    static get orange() { return 0xFFFFA500; }
    static get darkOrange() { return 0xFFFF8C00; }
    static get lightBlue() { return 0xFFADD8E6; }
    static get darkBlue() { return 0xFF0000A0; }
    static get purple() { return 0xFF800080; }
    static get whitesmoke() { return 0xFFF5F5F5; }
    static get navy() { return 0xFF000080; }
    static get cyan() { return 0xFF00ffff; }

    static get filteredRows() { return 0xff1f77b4; }
    static get filteredOutRows() { return Color.lightLightGray; }
    static get selectedRows() { return Color.darkOrange; }
    static get missingValueRows() { return Color.filteredOutRows; }
    static get mouseOverRows() { return 0xFFAAAAAA; }
    static get currentRow() { return 0xFF38B738; }

    static get histogramBar() { return Color.filteredRows; }
    static get barChart() { return 0xFF24A221; }
    static get scatterPlotMarker() { return 0xFF40699c; }
    static get scatterPlotSelection() { return 0x80323232; }
    static get scatterPlotZoom() { return 0x80626200; }

    static get areaSelection() { return Color.lightBlue; }
    static get rowSelection() { return 0x60dcdca0; } 
    static get colSelection() { return 0x60dcdca0; } 
    static get areaZoom() { return 0x80323232; }

    static get gridWarningBackground() { return 0xFFFFB9A7; }

    static get success() { return 0xFF3cb173; }
    static get failure() { return 0xFFeb6767; }    
}


class ml {
    static applyModel(name, table, columnNamesMap, showProgress) {
        if (columnNamesMap === undefined)
            columnNamesMap = new Map();
        return new Promise((resolve, reject) =>
            grok_ML_ApplyModel(name, table.d, (t) => resolve(new DataFrame(t)), columnNamesMap, showProgress));
    }

    static missingValuesImputation(table, impute, data, nearestNeighbours) {
        return new Promise((resolve, reject) =>
            grok_ML_MissingValuesImputation(table.d, impute, data, nearestNeighbours, () => resolve(table)));
    }

    static cluster(table, features, clusters) {
        return new Promise((resolve, reject) =>
            grok_ML_Cluster(table.d, features, clusters, () => resolve(table)));
    }

    static pca(table, features, components, center, scale) {
        return new Promise((resolve, reject) =>
            grok_ML_PCA(table.d, features, components, center, scale, () => resolve(table)));
    }

    static randomData(table, distribution, params, seed) {
        return new Promise((resolve, reject) =>
            grok_ML_RandomData(table.d, distribution, params, seed, () => resolve(table)));
    }
}


class chem {
    static similaritySearch(column, molecule, metric = METRIC_TANIMOTO, limit = 10, minScore = 0.7) {
        return new Promise((resolve, reject) => grok_Chem_SimilaritySearch(column.d, molecule, metric,
            limit, minScore, (t) => resolve(new DataFrame(t))));
    }

    static diversitySearch(column, metric = METRIC_TANIMOTO, limit = 10) {
        return new Promise((resolve, reject) => grok_Chem_DiversitySearch(column.d, metric, limit, (mols) => resolve(mols)));
    }

    static substructureSearch(column, pattern, isSmarts = true) {
        return new Promise((resolve, reject) => grok_Chem_SubstructureSearch(column.d, pattern, isSmarts, (bs) => resolve(new BitSet(bs))));
    }

    static rGroup(table, column, core) {
        return new Promise((resolve, reject) => grok_Chem_RGroup(table.d, column, core, () => resolve(table)));
    }

    static mcs(column) {
        return new Promise((resolve, reject) => grok_Chem_MCS(column.d, (mcs) => resolve(mcs)));
    }

    static descriptors(table, column, descriptors) {
        return new Promise((resolve, reject) => grok_Chem_Descriptors(table.d, column, descriptors, () => resolve(table)));
    }
}


METRIC_TANIMOTO = 'tanimoto';
METRIC_DICE = 'dice';
METRIC_COSINE = 'cosine';
METRIC_SOKAL = 'sokal';
METRIC_RUSSEL = 'russel';
METRIC_ROGOT_GOLDBERG = 'rogot-goldberg';
METRIC_KULCZYNSKI = 'kulczynski';
METRIC_MC_CONNAUGHEY = 'mc-connaughey';
METRIC_ASYMMETRIC = 'asymmetric';
METRIC_BRAUN_BLANQUET = 'braun-blanquet';
