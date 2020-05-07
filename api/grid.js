
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
