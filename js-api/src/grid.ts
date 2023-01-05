import {Cell, Column, DataFrame, Row} from './dataframe';
import {Viewer} from './viewer';
import {toDart, toJs} from './wrappers';
import {__obs, _sub, EventData, StreamSubscription} from './events';
import {_identityInt32, _toIterable} from './utils';
import {Observable} from 'rxjs';
import {RangeSlider} from './widgets';
import {SemType} from './const';
import {Property} from './entities';


let api = <any>window;
let _bytes = new Float64Array(4);

/** Represents a point. */
export class Point {
  x: number;
  y: number;

  constructor(x: number, y: number) {
    this.x = x;
    this.y = y;
  }

  /** Distance to the specified point. */
  distanceTo(p: Point): number {
    return Math.sqrt((this.x - p.x) * (this.x - p.x) + (this.y - p.y) * (this.y - p.y));
  }
}


/** Represents a rectangle. */
export class Rect {
  x: number;
  y: number;
  width: number;
  height: number;

  constructor(x: number, y: number, width: number, height: number) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  static fromDart(dart: any): Rect {
    api.grok_Rect_Pack(dart, _bytes);
    return new Rect(_bytes[0], _bytes[1], _bytes[2], _bytes[3]);
  }

  /** Rectangle of the specified size, with the specified center */
  fromCenterSize(cx: number, cy: number, width: number, height: number): Rect {
    return new Rect(cx - width / 2, cy - height / 2, width, height);
  }

  /** The midpoint of the rectangle along the x-axis. */
  get midX(): number {
    return this.x + this.width / 2;
  }

  /** The midpoint of the rectangle along the y-axis. */
  get midY(): number {
    return this.y + this.height / 2;
  }

  /** Left border position of the rectangle along the x-axis. */
  get left(): number {
    return this.x;
  }

  /** Top border position of the rectangle along the y-axis. */
  get top(): number {
    return this.y;
  }

  /** Right border position of the rectangle along the x-axis. */
  get right(): number {
    return this.x + this.width;
  }

  /** Bottom border position of the rectangle along the y-axis. */
  get bottom(): number {
    return this.y + this.height;
  }

  // --

  /**
   * Rectangle's top part of height {@link height}.
   * New rectangle with the same top border as the original, but with a height of {@link height}.
   * The result may be larger than the original rectangle along y-axis.
   * @param height {number} new rectangle height
   */
  getTop(height: number): Rect {
    return new Rect(this.left, this.top, this.width, height);
  }

  /**
   * Rectangle's bottom part of height {@link height}.
   * New rectangle with the same bottom border as the original, but with a height of {@link height}.
   * The result may be larger than the original rectangle along y-axis.
   * @param height {number} new rectangle height
   */
  getBottom(height: number): Rect {
    return new Rect(this.left, this.bottom - height, this.width, height);
  }

  /**
   * Rectangle's left part of width {@link width}.
   * New rectangle with the same left border as the original, but with a width of {@link width}.
   * The result may be larger than the original rectangle along x-axis.
   * @param width {number} new rectangle width
   */
  getLeft(width: number): Rect {
    return new Rect(this.left, this.top, width, this.height);
  }

  /**
   * Rectangle's right part of width {@link width}.
   * New rectangle with the same right border as the original, but with a width of {@link width}.
   * The result may be larger than the original rectangle along x-axis.
   * @param width {number} new rectangle width
   */
  getRight(width: number): Rect {
    return new Rect(this.right - width, this.top, width, this.height);
  }

  /**
   * New rectangle with the same top-left corner as the original, but with new {@link width} and {@link height}.
   * The result may be larger than the original rectangle on both axes.
   * @param width {number} new rectangle width
   * @param height {number} new rectangle height
   */
  getTopLeft(width: number, height: number): Rect {
    return new Rect(this.left, this.top, width, height);
  }

  /**
   * New rectangle with the same top-right corner as the original, but with new {@link width} and {@link height}.
   * The result may be larger than the original rectangle on both axes.
   * @param width {number} new rectangle width
   * @param height {number} new rectangle height
   */
  getTopRight(width: number, height: number): Rect {
    return new Rect(this.right - width, this.top, width, height);
  }

  /**
   * New rectangle with the same bottom-left corner as the original,
   * but with new {@link width} and {@link height}.
   * The result may be larger than the original rectangle on both axes.
   * @param width {number} new rectangle width
   * @param height {number} new rectangle height
   */
  getBottomLeft(width: number, height: number): Rect {
    return new Rect(this.left, this.bottom - height, width, height);
  }

  /**
   * New rectangle with the same bottom-right corner as the original,
   * but with new {@link width} and {@link height}.
   * The result may be larger than the original rectangle on both axes.
   * @param width {number} new rectangle width
   * @param height {number} new rectangle height
   */
  getBottomRight(width: number, height: number): Rect {
    return new Rect(this.right - width, this.bottom - height, width, height);
  }

  //--

  /**
   * New rectangle with the same right border as the original, but with the left border clipped by {@link dw}.
   * @param dw {number} delta width
   */
  cutLeft(dw: number): Rect {
    return new Rect(this.left + dw, this.top, this.width - dw, this.height);
  }

  /**
   * New rectangle with the same bottom border as the original, but with the top border clipped by {@link dh}.
   * @param dh {number} delta height
   */
  cutTop(dh: number): Rect {
    return new Rect(this.left, this.top + dh, this.width, this.height - dh);
  }

  /**
   * New rectangle with the same top border as the original, but with the bottom border clipped by {@link dh}.
   * @param dh {number} delta height
   */
  cutBottom(dh: number): Rect {
    return new Rect(this.left, this.top, this.width, this.height - dh);
  }

  /**
   * New rectangle with the same left border as the original, but with the right border clipped by {@link dw}.
   * @param dw {number} delta width
   */
  cutRight(dw: number): Rect {
    return new Rect(this.left, this.top, this.width - dw, this.height);
  }

  // --

  /** New rectangle below the original with a height of {@link height}.
   * @param height {number} new rectangle height
   */
  below(height: number): Rect {
    return new Rect(this.left, this.bottom, this.width, height);
  }

  /** New rectangle above the original with a height of {@link height}.
   * @param height {number} new rectangle height
   */
  above(height: number): Rect {
    return new Rect(this.left, this.top - height, this.width, height);
  }

  /** New rectangle to the left of the original with a width of {@link width}.
   * @param width {number} new rectangle width
   */
  toTheLeft(width: number): Rect {
    return new Rect(this.left - width, this.top, width, this.height);
  }

  /**
   * New rectangle to the right of the original with a width of {@link width}.
   * @param width {number} new rectangle width
   */
  toTheRight(width: number): Rect {
    return new Rect(this.right, this.top, width, this.height);
  }

  /**
   * New rectangle above the original with a height of {@link height}. See {@link above}.
   * @param height {number} new rectangle height
   */
  toTheTop(height: number): Rect {
    return this.above(height);
  }

  /**
   * New rectangle below the original with a height of {@link height}. See {@link below}.
   * @param height {number} new rectangle height
   */
  toTheBottom(height: number): Rect {
    return this.below(height);
  }

  // --

  /**
   * New rectangle with the same top border as the original, but with a height
   * scaled with {@link ratio} of the original height. See also {@link getTop}.
   * @param ratio {number} height scale factor
   */
  getTopScaled(ratio: number): Rect {
    return this.getTop(this.height * ratio);
  }

  /**
   * New rectangle with the same bottom border as the original, but with a height
   * scaled with {@link ratio} of the original height. See also {@link getBottom}.
   * @param ratio {number} height scale factor
   */
  getBottomScaled(ratio: number): Rect {
    return this.getBottom(this.height * ratio);
  }

  /**
   * New rectangle with the same left border as the original, but with a width
   * scaled with {@link ratio} of original width. See also {@link getLeft}.
   * @param ratio {number} width scale factor
   */
  getLeftScaled(ratio: number): Rect {
    return this.getLeft(this.width * ratio);
  }

  /**
   * New rectangle with the same right border as the original, but with a width
   * scaled with {@link ratio} of original width. See also {@link getRight}.
   * @param ratio {number} width scale factor
   */
  getRightScaled(ratio: number): Rect {
    return this.getRight(this.width * ratio);
  }

  // --

  /**
   * Horizontal part of {@link index} of the original rectangle sliced on {@link count} parts.
   * @param count {number} count of parts to divide the rectangle vertically
   * @param index {number} index of part to return
   */
  getTopPart(count: number, index: number): Rect {
    return new Rect(
      this.left, this.top + (this.height / count) * index,
      this.width, this.height / count);
  }

  /**
   * Vertical part of {@link index} of the original rectangle slices on {@link count} parts.
   * @param count {number} count of parts to divide rectangle horizontally
   * @param index {number} index of part to return
   */
  getLeftPart(count: number, index: number): Rect {
    return new Rect(
      this.left + (this.width / count) * index, this.top,
      this.width / count, this.height);
  }

  /**
   * Returns ({@link x}, {@link y}) part of the original rectangle divided into parts by a grid.
   * @param xCount {number} dividing grid size along x-axis
   * @param yCount {number} dividing grid size along y-axis
   * @param x {number} index of part to return along x-axis
   * @param y {number} index of part to return along y-axis
   */
  getGridPart(xCount: number, yCount: number, x: number, y: number): Rect {
    return new Rect(
      this.left + (this.width / xCount) * x, this.top + (this.height / yCount) * y,
      this.width / xCount, this.height / y);
  }

  // --

  /**
   * Inflated rectangle with left and right borders of the original shifted outside by {@link dx}
   * and top and bottom borders shifted outside by {@link dy}.
   * Overall size increased by 2*{@link dx} along x-axis, and by 2*{@link dy} along y-axis.
   * @param dx {number} vertical borders shift delta along x-axis
   * @param dy {number} horizontal borders shift delta along y-axis
   */
  inflate(dx: number, dy: number): Rect {
    return new Rect(
      this.left - dx, this.top - dy,
      this.width + 2 * dx, this.height + 2 * dy);
  }

  /**
   * Inflated rectangle with right border of the original shifted outside by {@link dw}
   * and bottom border shifted outside by {@link dh}.
   * @param dw {number} delta width
   * @param dh {number} delta height
   */
  inflateSize(dw: number, dh: number): Rect {
    return new Rect(this.left, this.top, this.width + dw, this.height + dh);
  }

  /**
   * Inflated rectangle with new width scaled by {@link dxRatio} of the original
   * and height scaled by {@link dyRatio}.
   * @param dxRatio {number} width scale ratio
   * @param dyRatio {number} height scale ratio
   */
  inflateRel(dxRatio: number, dyRatio: number): Rect {
    return this.inflate(this.width * (dxRatio - 1), this.height * (dyRatio - 1));
  }

  /** Square fitted (inscribed) to the original rectangle */
  fitSquare(): Rect {
    const size = Math.min(this.width, this.height);
    return new Rect(
      this.left + (this.width - size) / 2, this.top + (this.height - size) / 2, size, size);
  }

  /** The biggest rectangle that fits within this rect that keeps the width/height aspect ratio
   * Positioned in the center. Useful for rendering images in cells. */
  fit(width: number, height: number): Rect {
    return width / height > this.width / this.height
      ? this.fromCenterSize(this.midX, this.midY, this.width, height * (this.width / width))
      : this.fromCenterSize(this.midX, this.midY, width * (this.height / height), this.height)
  }

  /** Checks if this Rect contains the point (x; y) inside */
  contains(x: number, y: number): boolean {
    return this.left <= x && x <= this.right && this.top <= y && y <= this.bottom;
  }

  /** Checks if this Rect contains the specified point */
  containsPoint(p: Point): boolean {
    return this.contains(p.x, p.y);
  }
}

/** Represents a grid cell */
export class GridCell {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** @returns {GridCell} */
  static fromColumnRow(grid: Grid, columnName: string, gridRow: number): GridCell {
    return new GridCell(api.grok_Grid_GetCell(grid.dart, columnName, gridRow));
  }

  /** Returns a synthecic GridCell that only contains value but no row/col. Useful for rendering values. */
  static fromValue(value: any): GridCell {
    return new GridCell(api.grok_GridCell_FromValue(value));
  }

  /** @returns {string} Cell type */
  get cellType(): string {
    return api.grok_GridCell_Get_CellType(this.dart);
  }

  /** @returns {boolean} Whether this is a table (data) cell (as opposed to special cells like row headers). */
  get isTableCell(): boolean {
    return api.grok_GridCell_Get_IsTableCell(this.dart);
  }

  /** @returns {boolean} Whether this is a row header. */
  get isRowHeader(): boolean {
    return api.grok_GridCell_Get_IsRowHeader(this.dart);
  }

  /** @returns {boolean} Whether this is a column header. */
  get isColHeader(): boolean {
    return api.grok_GridCell_Get_IsColHeader(this.dart);
  }

  /** @returns {Column} Corresponding table column, or null. */
  get tableColumn(): Column | null {
    return this.gridColumn.column;
  }

  /** @returns {Row} Corresponding table row, or null. */
  get tableRow(): Row | null {
    return this.isTableCell || this.isRowHeader ? this.cell.row : null;
  }

  /** @returns {number|null} Index of the corresponding table row. */
  get tableRowIndex(): number | null {
    return this.isTableCell || this.isRowHeader ? this.cell.rowIndex : null;
  }

  /** @returns {number} Index of the corresponding grid row. */
  get gridRow(): number {
    return api.grok_GridCell_Get_GridRow(this.dart);
  }

  /** @returns {GridColumn} Corresponding grid column. */
  get gridColumn(): GridColumn {
    return new GridColumn(api.grok_GridCell_Get_GridColumn(this.dart));
  }

  /** Custom text to be shown in a cell . */
  get customText(): string { return api.grok_GridCell_Get_CustomText(this.dart); }
  set customText(x: string) { api.grok_GridCell_Set_CustomText(this.dart, x); }

  /** @returns {Grid} this cell belongs to. */
  get grid(): Grid {
    return new Grid(api.grok_GridCell_Get_Grid(this.dart));
  }

  /** @returns {Cell} Corresponding table cell. */
  get cell(): Cell {
    return new Cell(api.grok_GridCell_Get_Cell(this.dart));
  }

  /** @returns {GridCellStyle} Style to use for rendering. */
  get style(): GridCellStyle {
    return new GridCellStyle(api.grok_GridCell_Get_Style(this.dart));
  }

  /** Calculates cell background color, according to color-coding. */
  get color(): number { return api.grok_GridCell_Get_Color(this.dart); }

  /** Grid cell bounds.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/cell-bounds}
   */
  get bounds(): Rect {
    return Rect.fromDart(api.grok_GridCell_Get_Bounds(this.dart));
  }

  /** Returns grid cell renderer. */
  get renderer(): GridCellRenderer {
    return api.grok_GridCell_Get_Renderer(this.dart);
  }

  /** Gets or sets HTML element for this grid cell. */
  get element(): HTMLElement { return api.grok_GridCell_Get_Element(this.dart); }
  set element(e: HTMLElement) { api.grok_GridCell_Set_Element(this.dart, e); }
}

/** Represents a grid column */
export class GridColumn {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static fromDart(dart: any): GridColumn | null {
    return dart == null ? null : new GridColumn(dart);
  }

  /** A {@link Grid} this column is associated with */
  get grid(): Grid {
    return toJs(api.grok_GridColumn_Get_Grid(this.dart));
  }

  /** @returns {Column} Corresponding table column, or null. */
  get column(): Column | null {
    let col = api.grok_GridColumn_Get_Column(this.dart);
    return col === null ? null : new Column(col);
  }

  /** Index of the column. */
  get idx(): number {
    return api.grok_GridColumn_Get_Idx(this.dart);
  }

  /** @returns {string} Column name. */
  get name(): string { return api.grok_GridColumn_Get_Name(this.dart); }
  set name(x: string) { api.grok_GridColumn_Set_Name(this.dart, x); }

  /** Column width in pixels.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/resize-columns} */
  get width(): number { return api.grok_GridColumn_Get_Width(this.dart); }
  set width(x: number) { api.grok_GridColumn_Set_Width(this.dart, x); }

  /** Background column as a 4-byte ARGB number. */
  get backColor(): number { return api.grok_GridColumn_Get_BackColor(this.dart); }
  set backColor(x: number) { api.grok_GridColumn_Set_BackColor(this.dart, x); }

  /** Column format.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/html-markup-cells} */
  get format(): string { return api.grok_GridColumn_Get_Format(this.dart); }
  set format(x: string) { api.grok_GridColumn_Set_Format(this.dart, x); }

  /** @returns {string} Cell type. */
  get cellType(): string { return api.grok_GridColumn_Get_CellType(this.dart); }
  set cellType(x: string) { api.grok_GridColumn_Set_CellType(this.dart, x); }

  /** Grid cell renderer. */
  get renderer(): GridCellRenderer { return toJs(api.grok_GridColumn_Get_Renderer(this.dart)); }
  set renderer(x: GridCellRenderer) { api.grok_GridColumn_Set_Renderer(this.dart, toDart(x)); }

  /** Column visibility.  */
  get visible(): boolean { return api.grok_GridColumn_Get_Visible(this.dart); }
  set visible(x: boolean) { api.grok_GridColumn_Set_Visible(this.dart, x); }

  /** Custom colors for categories.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/category-colors}
   *  @returns {Object.<string, number>} */
  get categoryColors(): { [s: string]: number } { return api.grok_GridColumn_Get_CategoryColors(this.dart); }
  set categoryColors(x: { [s: string]: number }) { api.grok_GridColumn_Set_CategoryColors(this.dart, x); }

  /** Whether the column is editable.
   *  @returns {boolean} */
  get editable(): boolean { return api.grok_GridColumn_Get_Editable(this.dart); }
  set editable(x: boolean) { api.grok_GridColumn_Set_Editable(this.dart, x); }

  /** Whether the column is selected. */
  get selected(): boolean { return api.grok_GridColumn_Get_Selected(this.dart); }
  set selected(x: boolean) { api.grok_GridColumn_Set_Selected(this.dart, x); }

  /** Left border (in pixels in the virtual viewport) */
  get left(): number { return api.grok_GridColumn_Get_Left(this.dart); }

  /** Right border (in pixels in the virtual viewport) */
  get right(): number { return api.grok_GridColumn_Get_Right(this.dart); }

  /** Returns all visible cells */
  getVisibleCells(): Iterable<GridCell> { return this.grid.getVisibleCells(this); }

  /** Grid column settings. */
  get settings(): any | null { return api.grok_GridColumn_Get_Settings(this.dart); }
  set settings(s: any | null) { api.grok_GridColumn_Set_Settings(this.dart, s); }

  /** Moves the specified column to the specified position */
  move(position: number) { api.grok_GridColumnList_Move(this.grid.columns.dart, this.dart, position); }

  /** If this column is not entirely visible, scrolls the grid horizontally to show it. */
  scrollIntoView(): void { api.grok_GridColumn_ScrollIntoView(this.dart); }
}

/** Represents grid columns. */
export class GridColumnList {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Row header column.
   *  @returns {GridColumn}  */
  get rowHeader(): GridColumn | null {
    return this.byIndex(0);
  }

  /** Returns a grid column by index, or null if it does not exist.
   *  @param {number} index
   *  @returns {GridColumn}  */
  byIndex(index: number): GridColumn | null {
    return GridColumn.fromDart(api.grok_GridColumnList_ByIndex(this.dart, index));
  }

  /** Returns a grid column by name, or null if it does not exist.
   *  @param {string} columnName
   *  @returns {GridColumn}  */
  byName(columnName: string): GridColumn | null {
    return GridColumn.fromDart(api.grok_GridColumnList_ByName(this.dart, columnName));
  }

  /** Sets column order.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/order-columns}
   *  @param {string[]} columnNames - Order of columns. */
  setOrder(columnNames: string[]): void {
    api.grok_GridColumnList_SetOrder(this.dart, columnNames);
  }

  /** Shows the specified columns (and hides the rest).
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/hide-columns}
   *  @param {string[]} columnNames - Names of the columns to show. */
  setVisible(columnNames: string[]): void {
    api.grok_GridColumnList_SetVisible(this.dart, columnNames);
  }

  /** GridColumnList length.
   *  @returns {number}  */
  get length(): number {
    return api.grok_GridColumnList_Get_Length(this.dart);
  }

  /** Adds a new column to the grid (but not to the underlying dataframe). */
  add(options: {gridColumnName?: string, cellType: string, index?: number}): GridColumn {
    return api.grok_GridColumnList_Add(this.dart, options.cellType, options.gridColumnName);
  }
}


/** High-performance, flexible spreadsheet control */
export class Grid<TSettings = any> extends Viewer<TSettings> {

  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new grid. */
  static create(table: { dart: any; }): Grid {
    return new Grid(api.grok_Grid_Create(table.dart));
  }


  /** Creates a new grid from a list of items (rows) and their properties (columns) */
  static fromProperties(items: any[], props: Property[]): Grid {
    const t = DataFrame.create(items.length);
    for (let p of props)
      t.columns.addNewVirtual(p.name,
        (i: number) => p.get(items[i]), p.propertyType,
        p.set == null ? null : ((i: number, x: any) => p.set(items[i], x)));
    return Grid.create(t);
  }

  /** Grid columns.
   *  @returns {GridColumnList} */
  get columns(): GridColumnList {
    return new GridColumnList(api.grok_Grid_Get_Columns(this.dart));
  }

  /** Returns a column with the specified name.
   * @param {string} name
   * @returns {GridColumn} */
  col(name: string): GridColumn | null {
    return this.columns.byName(name);
  }

  /** Returns a grid cell at the specified position. */
  cell(columnName: string, gridRow: number): GridCell {
    return GridCell.fromColumnRow(this, columnName, gridRow);
  }

  /**
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/custom-cell-rendering-indexes}
   * @returns {Observable<GridCellRenderArgs>} */
  get onCellRender(): Observable<GridCellRenderArgs> {
    return __obs('d4-grid-cell-render', this.dart);
  }

  /** @returns {HTMLCanvasElement} */
  get canvas(): HTMLCanvasElement {
    return api.grok_Grid_Get_Canvas(this.dart);
  }

  /** @returns {HTMLCanvasElement} */
  get overlay(): HTMLCanvasElement {
    return api.grok_Grid_Get_Overlay(this.dart);
  }

  /** List of columns the grid is sorted by. */
  get sortByColumns(): Column[] { return toJs(api.grok_Grid_Get_SortByColumns(this.dart)); }

  /** Sort directions for [sortByColumns]: true = ascending order, false = descending order. */
  get sortTypes(): boolean[] { return api.grok_Grid_Get_SortTypes(this.dart); }

  /**
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/custom-cell-prepare}
   */
  onCellPrepare(callback: (cell: GridCell) => any): StreamSubscription {
    return _sub(api.grok_Grid_OnCellPrepare(this.dart, (dcell: any) => {
      return callback(new GridCell(dcell));
    }));
  }

  /**
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/custom-cell-tooltip}
   */
  onCellTooltip(callback: (cell: GridCell, x: number, y: number) => any): StreamSubscription {
    return _sub(api.grok_Grid_OnCellTooltip(this.dart, (dcell: any, x: any, y: any) => {
      return callback(new GridCell(dcell), x, y);
    }));
  }

  /** Returns a grid cell at the specified position, or null if there is none.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/hit-test} */
  hitTest(x: number, y: number): GridCell | null {
    return toJs(api.grok_Grid_HitTest(this.dart, x, y));
  }

  /** Sorts rows by the specified [columnIds].
   *  Specify sort directions via [asc] array (true = ascending, false = descending)
   *  If [asc] is not specified, sorts in ascending order. */
  sort(columns: string[] | Column[], orders: boolean[] | null = null): Grid {
    api.grok_Grid_Sort(this.dart, (columns as any[]).map((c: string | Column) => c instanceof Column ? c.dart : c), orders);
    return this;
  }

  /** Sorts the rows, using the specified comparer that accepts indexes.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/order-rows-by-comparer}
   * Also, @see setRowOrder */
  sortIndexes(indexComparer: (a: number, b: number) => number): Grid {
    let indexes = _identityInt32(this.table.rowCount);
    indexes.sort(indexComparer);
    this.setRowOrder(<number[]><unknown>indexes);
    return this;
  }

  /** Sets the order or rows in the table.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/order-rows}
   * Also, @see sortIndexes */
  setRowOrder(indexes: number[]): Grid {
    api.grok_Grid_SetRowOrder(this.dart, indexes);
    return this;
  }

  /** Pinned rows.
   *  @returns {Iterable<number>} */
  get pinnedRows(): Iterable<number> {
    return _toIterable(api.grok_Grid_Get_PinnedRows(this.dart));
  }

  /** Returns a grid cell at the specified position.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/scroll-to-cell}
   */
  scrollToCell(column: string | Column, row: number): void {
    api.grok_Grid_ScrollToCell(this.dart, toDart(column), row);
  }

  /** Scrolls the grid to the specified position.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/scroll-to-pixels}
   */
  scrollToPixels(x: number, y: number): void {
    api.grok_Grid_ScrollToPixels(this.dart, x, y);
  }

  /** Causes the grid to repaint.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/advanced/invalidate}
   **/
  invalidate(): void { api.grok_Grid_Invalidate(this.dart); }

  /** Vertical scroll bar.
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/scroll-bars}
   */
  get vertScroll(): RangeSlider {
    return toJs(api.grok_Grid_Get_VertScroll(this.dart));
  }

  /** Horizontal scroll bar */
  get horzScroll(): RangeSlider {
    return toJs(api.grok_Grid_Get_HorzScroll(this.dart));
  }

  /** Column labels height */
  get colHeaderHeight(): number {
    return toJs(api.grok_Grid_Get_ColHeaderHeight(this.dart));
  }

  /** Column labels box */
  get colHeaderBox(): Rect {
    return toJs(api.grok_Grid_Get_ColHeaderRect(this.dart));
  }

  /** Forces the grid to execute calculations that were postponed before the next rendering (such as recalculating layout).
   * Call it in case the client code changes column widths and needs to access the new recalculated layout. */
  runPostponedComputations(): void {
    api.grok_CanvasViewer_RunPostponedComputations(this.dart);
  }

  /**
   * Sample: {@link https://public.datagrok.ai/js/samples/grid/resize-events}
   */
  get onColumnResized(): Observable<any> {
    return __obs('d4-grid-column-resized', this.dart);
  }

  /** Sample: {@link https://public.datagrok.ai/js/samples/grid/resize-events} */
  get onRowsResized(): Observable<any> { return __obs('d4-grid-rows-resized', this.dart); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/grid/order-rows} */
  get onRowsSorted(): Observable<any> { return __obs('d4-grid-rows-sorted', this.dart); }
  get onPinnedRowsChanged(): Observable<any> { return __obs('d4-grid-pinned_rows-changed', this.dart); }

  get onCellValueEdited(): Observable<GridCell> { return __obs('d4-grid-cell-value-edited', this.dart); }
  get onCurrentCellChanged(): Observable<GridCell> { return __obs('d4-grid-current-cell-changed', this.dart); }
  get onCellClick(): Observable<GridCell> { return __obs('d4-grid-cell-click', this.dart); }
  get onCellDoubleClick(): Observable<GridCell> { return __obs('d4-grid-cell-double-click', this.dart); }
  get onCellMouseDown(): Observable<GridCell> { return __obs('d4-grid-cell-mouse-down', this.dart); }
  get onCellKeyDown(): Observable<GridCell> { return __obs('d4-grid-cell-key-down', this.dart); }
  get onRowEnter(): Observable<GridCell> { return __obs('d4-grid-row-enter', this.dart); }

  get onBeforeDrawOverlay(): Observable<EventData> { return __obs('d4-grid-before-draw-overlay', this.dart); }
  get onAfterDrawOverlay(): Observable<EventData> { return __obs('d4-grid-after-draw-overlay', this.dart); }
  get onBeforeDrawContent(): Observable<EventData> { return __obs('d4-grid-before-draw-content', this.dart); }
  get onAfterDrawContent(): Observable<EventData> { return __obs('d4-grid-after-draw-content', this.dart); }

  /** Returns currently visible cells */
  getVisibleCells(column: GridColumn | null = null): Iterable<GridCell> {
    return _toIterable(api.grok_Grid_GetVisibleCells(this.dart, column?.dart))
  }

  /** Resizes the grid to fit the content */
  autoSize(maxWidth: number, maxHeight: number, minWidth?: number, minHeight?: number, autoSizeOnDataChange?: boolean): void {
    api.grok_Grid_AutoSize(this.dart, maxWidth, maxHeight, minWidth, minHeight, autoSizeOnDataChange);
  }
}


/** Represents grid cell style. */
export class GridCellStyle {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Font. Example: 12px Verdana */
  get font(): string { return api.grok_GridCellStyle_Get_Font(this.dart); }
  set font(x: string) { api.grok_GridCellStyle_Set_Font(this.dart, x); }

  /** Text color (RGBA-encoded) */
  get textColor(): number { return api.grok_GridCellStyle_Get_TextColor(this.dart); }
  set textColor(x: number) { api.grok_GridCellStyle_Set_TextColor(this.dart, x); }

  /** Background color (RGBA-encoded) */
  get backColor(): number { return api.grok_GridCellStyle_Get_BackColor(this.dart); }
  set backColor(x: number) { api.grok_GridCellStyle_Set_BackColor(this.dart, x); }

  /** DOM Element to put in the cell */
  get element(): HTMLElement { return api.grok_GridCellStyle_Get_Element(this.dart); }
  set element(x: HTMLElement) { api.grok_GridCellStyle_Set_Element(this.dart, x); }

  /** Vertical text orientation */
  get textVertical(): boolean { return api.grok_GridCellStyle_Get_TextVertical(this.dart); }
  set textVertical(x: boolean) { api.grok_GridCellStyle_Set_TextVertical(this.dart, x); }
}


/** Grid cell rendering args. */
export class GridCellRenderArgs extends EventData {

  constructor(dart: any) {
    super(dart);
  }

  /** @returns {CanvasRenderingContext2D} */
  get g(): CanvasRenderingContext2D {
    return api.grok_GridCellRenderArgs_Get_G(this.dart);
  }

  /** Cell that is being rendered. */
  get cell(): GridCell {
    return new GridCell(api.grok_GridCellRenderArgs_Get_Cell(this.dart));
  }

  /** Cell bounds to render to. */
  get bounds(): Rect {
    return api.grok_GridCellRenderArgs_Get_Bounds(this.dart);
  }
}


export class CanvasRenderer {
  get defaultWidth(): number | null {
    return null;
  }

  get defaultHeight(): number | null {
    return null;
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, value: any, context: any): void {
    throw 'Not implemented';
  }
}


export class GridCellRenderer extends CanvasRenderer {
  get name(): string {
    throw 'Not implemented';
  }

  get cellType(): string {
    throw 'Not implemented';
  }

  renderSettings(gridColumn: GridColumn): Element | null { return null; }

  renderInternal(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: GridCell, cellStyle: GridCellStyle): void {
    this.render(g, x, y, w, h, new GridCell(gridCell), new GridCellStyle(cellStyle));
  }

  static register(renderer: any): void {
    api.grok_GridCellRenderer_Register(renderer);
  }

  onKeyDown(gridCell: GridCell, e: KeyboardEvent): void {}
  onKeyPress(gridCell: GridCell, e: KeyboardEvent): void {}

  onMouseEnter(gridCell: GridCell, e: MouseEvent): void {}
  onMouseLeave(gridCell: GridCell, e: MouseEvent): void {}
  onMouseDown(gridCell: GridCell, e: MouseEvent): void {}
  onMouseUp(gridCell: GridCell, e: MouseEvent): void {}
  onMouseMove(gridCell: GridCell, e: MouseEvent): void {}
  onClick(gridCell: GridCell, e: MouseEvent): void {}
  onDoubleClick(gridCell: GridCell, e: MouseEvent): void {}
}


/** Proxy class for the Dart-based grid cell renderers. */
export class GridCellRendererProxy extends GridCellRenderer {
  dart: any;

  constructor(dart: any) {
    super();
    this.dart = dart;
  }

  get name(): string { return api.grok_GridCellRenderer_Get_CellType(this.dart); }
  get cellType(): string { return api.grok_GridCellRenderer_Get_CellType(this.dart); }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, value: any, context: any): void {
    this.renderInternal(g, x, y, w, h, value, context);
  }

  renderInternal(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: GridCell, cellStyle: GridCellStyle): void {
    api.grok_GridCellRenderer_Render(this.dart, g, x, y, w, h, gridCell.dart);
  }
}


export class SemanticValue<T = any> {
  private readonly dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static fromValueType(value: any, semType: SemType | null, units?: string) {
    const v = new SemanticValue(api.grok_SemanticValue(value, semType));
    if (units)
      v.units = units;
    return v;
  }

  get value(): T { return api.grok_SemanticValue_Get_Value(this.dart); }
  set value(x: T) { api.grok_SemanticValue_Set_Value(this.dart, x); }

  get units(): any { return api.grok_SemanticValue_Get_Units(this.dart); }
  set units(x: any) { api.grok_SemanticValue_Set_Units(this.dart, x); }

  get dataType(): string { return api.grok_SemanticValue_Get_DataType(this.dart); }

  get semType(): string { return api.grok_SemanticValue_Get_SemType(this.dart); }
  set semType(x: string) { api.grok_SemanticValue_Set_SemType(this.dart, x); }

  getMeta(name: string): any { return api.grok_SemanticValue_Get_Meta(this.dart, name); }
  setMeta(name: string, value: any): void { api.grok_SemanticValue_Set_Meta(this.dart, name, toDart(value)); }

  get cell(): Cell { return api.grok_SemanticValue_Get_Cell(this.dart); }
  get gridCell(): GridCell { return this.getMeta('gridCell'); }
  get viewer(): Viewer { return this.getMeta('viewer'); }
}
