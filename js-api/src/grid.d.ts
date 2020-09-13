import {Cell, Column, Row} from "./dataframe";
import {Viewer} from "./viewer";
import {Observable} from 'rxjs';
import {EventData, StreamSubscription} from "./events";
import {_identityInt32} from "./utils";

export class Point {
    constructor(x: number, y: number);
    x: number;
    y: number;
}

export class Rect {
    x: number;
    y: number;
    width: number;
    height: number;
    
    constructor(x: number, y: number, width: number, height: number)

    static fromDart(d: any): Rect

    get midX(): number
    get midY(): number

    get left(): number
    get top(): number
    get right(): number
    get bottom(): number
}

/** Represents a grid cell */
export class GridCell {
    constructor(d: any)

    /** @returns {GridCell} */
    static fromColumnRow(grid: Grid, columnName: string, gridRow: number): GridCell
    
    /** @returns {string} Cell type */
    get cellType(): string

    /** @returns {boolean} Whether this is a table (data) cell (as opposed to special cells like row headers). */
    get isTableCell(): boolean

    /** @returns {boolean} Whether this is a row header. */
    get isRowHeader(): boolean

    /** @returns {boolean} Whether this is a column header. */
    get isColHeader(): boolean

    /** @returns {Column} Corresponding table column, or null. */
    get tableColumn(): Column | null

    /** @returns {Row} Corresponding table row, or null. */
    get tableRow(): Row | null

    /** @returns {number|null} Index of the corresponding table row. */
    get tableRowIndex(): number | null

    /** @returns {number} Index of the corresponding grid row. */
    get gridRow(): number

    /** @returns {GridColumn} Corresponding grid column. */
    get gridColumn(): GridColumn

    /** Custom text to be shown in a cell . */
    get customText(): string

    set customText(x: string)

    /** @returns {Grid} this cell belongs to. */
    get grid(): Grid

    /** @returns {Cell} Corresponding table cell. */
    get cell(): Cell

    /** @returns {GridCellStyle} Style to use for rendering. */
    get style(): GridCellStyle
    
    get bounds(): Rect
}

/** Represents a grid column */
export class GridColumn {
    constructor(d: any)

    /** @returns {Column} Corresponding table column, or null. */
    get column(): Column | null

    /** Index of the column.
     *  @returns {number} */
    get idx(): number

    /** @returns {string} Column name. */
    get name(): string

    set name(x: string)

    /** Column width in pixels.
     *  @returns {number} */
    get width(): number

    set width(x: number)

    /** Background column as a 4-byte ARGB number.
     *  @returns {number} */
    get backColor(): number

    set backColor(x: number)

    /** Column format.
     *  @returns {string} */
    get format(): string

    set format(x)

    /** @returns {string} Cell type. */
    get cellType(): string

    set cellType(x)

    /** Column visibility.
     *  @returns {string} */
    get visible(): string

    set visible(x)

    /** Custom colors for categories.
     *  @returns {Object.<string, number>} */
    get categoryColors(): { [key: string]: number }

    set categoryColors(x)

    /** Whether the column is editable.
     *  @returns {boolean} */
    get editable(): boolean

    set editable(x)

    /** Whether the column is selected.
     *  @returns {boolean}  */
    get selected(): boolean

    set selected(x)

    static fromDart(d: null): null

    static fromDart(d: any): GridColumn
}


/** Represents grid columns. */
export class GridColumnList {
    constructor(d: any)

    /** Row header column.
     *  @returns {GridColumn}  */
    get rowHeader(): GridColumn

    /** Returns a grid column by name, or null if it does not exist.
     *  @param {number} index
     *  @returns {GridColumn}  */
    byIndex(index: number): GridColumn
    
    /** Returns a grid column by name, or null if it does not exist.
     *  @param {string} columnName
     *  @returns {GridColumn}  */
    byName(columnName: string): GridColumn | null

    /** Sets column order.
     *  @param {string[]} columnNames - Order of columns. */
    setOrder(columnNames: string[]): void

    /** Shows the specified columns (and hides the rest).
     *  @param {string[]} columnNames - Names of the columns to show. */
    setVisible(columnNames: string[]): void
}


/** High-performance, flexible spreadsheet control */
export class Grid extends Viewer {
    constructor(d: any)
    
    
    /** Grid columns.
     *  @returns {GridColumnList} */
    get columns(): GridColumnList

    /** @returns {Observable<GridCellRenderArgs>} */
    get onCellRender(): Observable<GridCellRenderArgs>

    static create(table: { d: any }): Grid

    /** Sorts data by [columnIds].
     *  Specify sort directions via [asc] array (true = ascending, false = descending)
     *  If [asc] is not specified, sorts in ascending order. */
    sort(columns: Column[], orders?: boolean[] | null): void

    sortIndexes(indexComparer: ((a: number, b: number) => number | undefined)): this

    setRowOrder(indexes: number[]): this
    
    /** Returns a column with the specified name.
     * @param {string} name
     * @returns {GridColumn} */
    col(name: string): GridColumn
    
    cell(columnName: string, gridRow: number): GridCell
    
    /** @returns {HTMLCanvasElement} */
    get canvas(): HTMLCanvasElement

    /** @returns {HTMLCanvasElement} */
    get overlay(): HTMLCanvasElement


    onCellPrepare(callback: (cell: GridCell) => any): StreamSubscription

    onCellTooltip(callback: (cell: GridCell) => any): StreamSubscription
    
    hitTest(x: number, y: number): GridCell
}


export class GridCellStyle {
    constructor(d: any)

    //TODO: font type
    get font(): string

    set font(x)

    get textColor(): string

    set textColor(x)

    get backColor(): string

    set backColor(x)

    get element(): HTMLHtmlElement

    set element(x)
}

/** Grid cell rendering args. */
export class GridCellRenderArgs extends EventData {
    constructor(d: any)

    /** @returns {CanvasRenderingContext2D} */
    get g(): CanvasRenderingContext2D

    get cell(): GridCell

    //TODO: what rect?
    /** @returns {Rect} */
    get bounds(): DOMRect
}


export class GridCellRenderer {
    get name(): string
    get cellType(): string

    /**
     * @param {CanvasRenderingContext2D} g
     * @param {number} x
     * @param {number} y
     * @param {number} w
     * @param {number} h
     * @param {GridCell} gridCell
     * @param {GridCellStyle} cellStyle
     **/
    render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: GridCell, cellStyle: GridCellStyle): void
}