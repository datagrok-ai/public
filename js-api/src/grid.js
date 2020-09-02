import {Viewer} from "./viewer.js";
import {Column} from './dataframe.js';
import {Cell} from "./dataframe";
import {__obs, _sub, EventData} from "./events";
import {_identityInt32} from "./utils";

let _bytes = new Float64Array(4);

export class Point {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }
}

export class Rect {
    constructor(x, y, width, height) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
    }

    static fromDart(d) {
        grok_Rect_Pack(d, _bytes);
        return new Rect(_bytes[0], _bytes[1], _bytes[2], _bytes[3]);
    }

    get left() { return x; }
    get top() { return y; }
    get right() { return x + width; }
    get bottom() { return y + height; }
}

/** Represents a grid cell */
export class GridCell {
    constructor(d) { this.d = d; }

    /** @returns {GridCell} */
    static fromColumnRow(grid, columnName, gridRow) {
        return new GridCell(grok_Grid_GetCell(grid.d, columnName, gridRow));
    }

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

    get bounds() { return Rect.fromDart(grok_GridCell_Get_Bounds(this.d)); }
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

    /** Row header column.
     *  @returns {GridColumn}  */
    get rowHeader() { return this.byIndex(0); }

    /** Returns a grid column by name, or null if it does not exist.
     *  @param {number} index
     *  @returns {GridColumn}  */
    byIndex(index) { return GridColumn.fromDart(grok_GridColumnList_ByIndex(this.d, index)); }

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

    /** Grid columns.
     *  @returns {GridColumnList} */
    get columns() { return new GridColumnList(grok_Grid_Get_Columns(this.d)); }

    /** Returns a column with the specified name.
     * @param {string} name
     * @returns {GridColumn} */
    col(name) { return this.columns.byName(name); }

    cell(columnName, gridRow) { return GridCell.fromColumnRow(this, columnName, gridRow); }

    /** @returns {Observable<GridCellRenderArgs>} */
    get onCellRender() { return __obs('d4-grid-cell-render', this.d); }

    /** @returns {HTMLCanvasElement} */
    get canvas() { return grok_Grid_Get_Canvas(this.d); }

    /** @returns {HTMLCanvasElement} */
    get overlay() { return grok_Grid_Get_Overlay(this.d); }

    onCellPrepare(callback) {
        return _sub(grok_Grid_OnCellPrepare(this.d, (dcell) => { return callback(new GridCell(dcell)); }));
    }

    onCellTooltip(callback) {
        return _sub(grok_Grid_OnCellTooltip(this.d, (dcell, x, y) => { return callback(new GridCell(dcell), x, y); }));
    }

    hitTest(x, y) { return new GridCell(grok_Grid_HitTest(this.d, x, y)); }

    /** Sorts rows by the specified [columnIds].
     *  Specify sort directions via [asc] array (true = ascending, false = descending)
     *  If [asc] is not specified, sorts in ascending order. */
    sort(columns, orders = null) {
        grok_Grid_Sort(this.d, columns.map((c) => c instanceof Column ? c.d : c), orders);
        return this;
    }

    sortIndexes(indexComparer) {
        let indexes = _identityInt32(this.table.rowCount);
        indexes.sort(indexComparer);
        this.setRowOrder(indexes);
        return this;
    }

    setRowOrder(indexes) {
       grok_Grid_SetRowOrder(this.d, indexes);
       return this;
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


/** Grid cell rendering args. */
export class GridCellRenderArgs extends EventData {
    constructor(d) { super(d); }

    /** @returns {CanvasRenderingContext2D} */
    get g() { return grok_GridCellRenderArgs_Get_G(this.d); }

    get cell() { return new GridCell(grok_GridCellRenderArgs_Get_Cell(this.d)); }

    /** @returns {Rect} */
    get bounds() { return grok_GridCellRenderArgs_Get_Bounds(this.d); }
}


export class GridCellRenderer {
    get name() { throw 'Not implemented'; }
    get cellType() { throw 'Not implemented'; }

    /**
     * @param {CanvasRenderingContext2D} g
     * @param {number} x
     * @param {number} y
     * @param {number} w
     * @param {number} h
     * @param {GridCell} gridCell
     * @param {GridCellStyle} cellStyle
     **/
    render(g, x, y, w, h, gridCell, cellStyle) {
        throw 'Not implemented';
    }

    renderInternal(g, x, y, w, h, gridCell, cellStyle) {
        this.render(g, x, y, w, h, new GridCell(gridCell), new GridCellStyle(cellStyle));
    }
}