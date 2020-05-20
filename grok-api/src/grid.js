import {EventData} from "./events.js";
import {Viewer} from "./viewer.js";
import {Column} from './dataframe.js';
import {Cell} from "./dataframe";

/** Represents a grid cell */
export class GridCell {
    constructor(d) { this.d = d; }

    /** @returns {string} Cell type */
    get cellType() { return grok_GridCell_Get_CellType(this.d); }

    /** @returns {boolean} Whether this is a table (data) cell (as opposed to special cells like row headers). */
    get isTableCell() { return grok_GridCell_Get_IsTableCell(this.d); }

    /** @returns {boolean} Whether this is a row header. */
    get isRowHeader() { return grok_GridCell_Get_IsRowHeader(this.d); }

    /** @returns {boolean} Whether this is a column header. */
    get isColHeader() { return grok_GridCell_Get_IsColHeader(this.d); }

    /** @returns {Column} Corresponding table column, or null. */
    get tableColumn() { return this.gridColumn.column; }

    /** @returns {Row} Corresponding table row, or null. */
    get tableRow() { return this.isTableCell ? this.cell.row : null; }

    /** @returns {number|null} Index of the corresponding table row. */
    get tableRowIndex() { return this.isTableCell ? this.cell.rowIndex : null; }

    /** @returns {number} Index of the corresponding grid row. */
    get gridRow() { return grok_GridCell_Get_GridRow(this.d); }

    /** @returns {GridColumn} Corresponding grid column. */
    get gridColumn() { return new GridColumn(grok_GridCell_Get_GridColumn(this.d)); }

    /** Custom text to be shown in a cell . */
    get customText() { return grok_GridCell_Get_CustomText(this.d); }
    set customText(x) { return grok_GridCell_Set_CustomText(this.d, x); }

    /** @returns {Grid} this cell belongs to. */
    get grid() { return new Grid(grok_GridCell_Get_Grid(this.d)); }

    /** @returns {Cell} Corresponding table cell. */
    get cell() { return new Cell(grok_GridCell_Get_Cell(this.d)); }

    /** @returns {GridCellStyle} Style to use for rendering. */
    get style() { return new GridCellStyle(grok_GridCell_Get_Style(this.d)); }
}

/** Grid cell rendering args. */
export class GridCellRenderArgs extends EventData {
    constructor(d) { super(d); }

    /** @returns {CanvasRenderingContext2D} */
    get g() { return grok_GridCellRenderArgs_Get_G(this.d); }

    get cell() { return new GridCell(grok_GridCellRenderArgs_Get_Cell(this.d)); }

    /** @returns {Rect} */
    get bounds() { return grok_GridCellRenderArgs_Get_Bounds(this.d); }
}


/** Represents a grid column */
export class GridColumn {
    constructor(d) { this.d = d; }
    static fromDart(d) { return d == null ? null : new GridColumn(d); }

    /** @returns {Column} Corresponding table column, or null. */
    get column() {
        let col = grok_GridColumn_Get_Column(this.d);
        return col === null ? null : new Column(col);
    }

    /** Index of the column.
     *  @returns {number} */
    get idx() { return grok_GridColumn_Get_Idx(this.d); }

    /** @returns {string} Column name. */
    get name() { return grok_GridColumn_Get_Name(this.d); }
    set name(x) { return grok_GridColumn_Set_Name(this.d, x); }

    /** Column width in pixels.
     *  @returns {number} */
    get width() { return grok_GridColumn_Get_Width(this.d); }
    set width(x) { return grok_GridColumn_Set_Width(this.d, x); }

    /** Background column as a 4-byte ARGB number.
     *  @returns {number} */
    get backColor() { return grok_GridColumn_Get_BackColor(this.d); }
    set backColor(x) { return grok_GridColumn_Set_BackColor(this.d, x); }

    /** Column format.
     *  @returns {string} */
    get format() { return grok_GridColumn_Get_Format(this.d); }
    set format(x) { return grok_GridColumn_Set_Format(this.d, x); }

    /** @returns {string} Cell type. */
    get cellType() { return grok_GridColumn_Get_CellType(this.d); }
    set cellType(x) { return grok_GridColumn_Set_CellType(this.d, x); }

    /** Column visibility.
     *  @returns {string} */
    get visible() { return grok_GridColumn_Get_Visible(this.d); }
    set visible(x) { return grok_GridColumn_Set_Visible(this.d, x); }

    /** Custom colors for categories.
     *  @returns {Object.<string, number>} */
    get categoryColors() { return grok_GridColumn_Get_CategoryColors(this.d); }
    set categoryColors(x) { return grok_GridColumn_Set_CategoryColors(this.d, x); }

    /** Whether the column is editable.
     *  @returns {boolean} */
    get editable() { return grok_GridColumn_Get_Editable(this.d); }
    set editable(x) { return grok_GridColumn_Set_Editable(this.d, x); }

    /** Whether the column is selected.
     *  @returns {boolean}  */
    get selected() { return grok_GridColumn_Get_Selected(this.d); }
    set selected(x) { return grok_GridColumn_Set_Selected(this.d, x); }
}


/** Represents grid columns. */
export class GridColumnList {
    constructor(d) { this.d = d; }

    /** Returns a grid column by name, or null if it does not exist.
     *  @param {string} columnName
     *  @returns {GridColumn}  */
    byName(columnName) { return GridColumn.fromDart(grok_GridColumnList_ByName(this.d, columnName)); }

    /** Sets column order.
     *  @param {string[]} columnNames - Order of columns. */
    setOrder(columnNames) { grok_GridColumnList_SetOrder(this.d, columnNames); }

    /** Shows the specified columns (and hides the rest).
     *  @param {string[]} columnNames - Names of the columns to show. */
    setVisible(columnNames) { grok_GridColumnList_SetVisible(this.d, columnNames); }
}


/** High-performance, flexible spreadsheet control */
export class Grid extends Viewer {
    constructor(d) { super(d); }

    /** Sorts data by [columnIds].
     *  Specify sort directions via [asc] array (true = ascending, false = descending)
     *  If [asc] is not specified, sorts in ascending order. */
    sort(columns, orders = null) { grok_Grid_Sort(this.d, columns.map((c) => c instanceof Column ? c.d : c), orders); }

    /** Grid columns.
     *  @returns {GridColumnList} */
    get columns() { return new GridColumnList(grok_Grid_Get_Columns(this.d)); }

    /** Returns a column with the specified name.
     * @param {string} name
     * @returns {GridColumn} */
    col(name) { return this.columns.byName(name); }

    /** @returns {Observable<GridCellRenderArgs>} */
    get onCellRender() { return __obs('d4-grid-cell-render', this.d); }

    onCellPrepare(callback) {
        return _sub(grok_Grid_OnCellPrepare(this.d, (dcell) => { return callback(new GridCell(dcell)); }));
    }

    onCellTooltip(callback) {
        return _sub(grok_Grid_OnCellTooltip(this.d, (dcell, x, y) => { return callback(new GridCell(dcell), x, y); }));
    }

    static create(table) { return new Grid(grok_Grid_Create(table.d)); }
}


export class GridCellStyle {
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
