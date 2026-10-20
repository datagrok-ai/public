/* Virtualized data table — `VirtualRows`' scroller and selection model with a row shape of pooled
   cells under a sticky header, plus the editor's per-cell verdicts. Cells are pooled with their
   row, so scrolling a wide table costs the visible window and nothing else. `BasicTable` stays the
   small-data control: every row in the DOM, real table semantics. */

import {Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {VirtualRows} from './list.js';
import type {VirtualRowsOptions} from './list.js';
import type {RowsLike} from '../../sources/rows-like.js';

export interface DataTableColumn<T> {
  /** The item's property this column shows — and the column name the cell state is keyed by. */
  name: string;
  /** The header text; the column name by default. */
  header?: string;
  /** A CSS grid track (`'120px'`, `'2fr'`); an equal share of the width by default. */
  width?: string;
  align?: 'left' | 'right' | 'center';
  /** A string becomes the cell text, an element its content; `cell` is the pooled cell element,
   * reset before each call — set `title` and extra classes on it. The item's own `name` value by
   * default. */
  render?: (item: T, index: number, cell: HTMLElement) => HTMLElement | string;
}

/** What the table asks about a cell, keyed by ROW KEY and column — platform-free, so the memory
 * writer and the js-api `DomainFrameEditor` both answer it (dg adapts `IFrameEditor` onto it). */
export interface CellStateLike {
  /** An unsaved edit: the cell reads `u2-cell-changed`, the grid's `DIRTY_CELL_COLOR`. */
  isChanged(key: string, column: string): boolean;
  /** A cell the caller may not write: `u2-cell-readonly`. */
  canEdit?(key: string, column: string): boolean;
  /** A refusal over the cell: `u2-cell-error` (the grid's `INVALID_CELL_COLOR`) and the message
   * as the cell's `title`. */
  errorOf(key: string, column: string): {message: string, kind: string} | null;
}

export interface DataTableOptions<T> extends Omit<VirtualRowsOptions<T>, 'itemHeight'> {
  columns: DataTableColumn<T>[];
  items?: RowsLike<T> | ReadonlySignal<T[]>;
  /** Row and header height in pixels; 24 by default. */
  rowHeight?: number;
  cellState?: CellStateLike;
  /** Enter on the selected row, or a double-click on it. */
  onActivate?(item: T, index: number): void;
  /** Items are live views (a frame row proxy) whose content changes without a new object: every
   * items change renders the window again. Default: on for a `RowsLike` or a signal (a source's
   * rows), off for a plain array — whose rows are rendered again only when their item is a new
   * object, or on {@link VirtualRows.refresh}. */
  liveItems?: boolean;
  /** The text a row is found by (`revealItem`, the state block): the first column's value by
   * default. */
  labelOf?: (item: T) => string;
}

export class DataTable<T> extends VirtualRows<T> {
  private _columns: DataTableColumn<T>[];
  private readonly _cellState: CellStateLike | undefined;
  private readonly _onActivate: ((item: T, index: number) => void) | undefined;
  private readonly _optionKeyOf: ((item: T) => string) | undefined;
  private readonly _liveOption: boolean | undefined;
  private _template: string;
  private _head: HTMLElement;
  private readonly _cellBaseline = new Map<string, string>();
  /** Bumped by a {@link columns} change: a row whose cells are older is given new ones on its next fill. */
  private _shape = 0;
  private readonly _rowShape = new WeakMap<HTMLElement, number>();
  /** Whether the current items are live — decided per {@link setItems}, so the base keeps every
   * row whose item is unchanged and {@link placeRow} asks for a render where they are live. */
  private _liveNow = false;

  constructor(options: DataTableOptions<T>) {
    const first = options.columns[0];
    super({...options, itemHeight: options.rowHeight ?? 24, liveItems: false,
      labelOf: options.labelOf ?? (first === undefined ? undefined : (item: T) => DataTable._text(item, first.name))},
    {host: 'u2-data-table', u2: 'data-table', role: 'grid'});
    this._columns = options.columns;
    this._cellState = options.cellState;
    this._onActivate = options.onActivate;
    this._optionKeyOf = options.keyOf;
    this._liveOption = options.liveItems;
    this._template = DataTable._templateOf(options.columns);
    this._head = this._header();
    this.root.prepend(this._head);
    if (options.items)
      this.setItems(options.items);
  }

  /** `name`, or `name:Header` — the column a spec's string names. */
  static column<T>(spec: string): DataTableColumn<T> {
    const at = spec.indexOf(':');
    return at < 0 ? {name: spec} : {name: spec.slice(0, at), header: spec.slice(at + 1)};
  }

  /** The columns shown, in order; a string entry is read by {@link column}. A change rebuilds the
   * header and gives each row its new cells on its next render — the rows, their keys and the
   * selection stay. */
  get columns(): DataTableColumn<T>[] {
    return this._columns;
  }

  set columns(x: readonly (DataTableColumn<T> | string)[]) {
    this._columns = x.map((c) => typeof c === 'string' ? DataTable.column<T>(c) : c);
    this._template = DataTable._templateOf(this._columns);
    this._shape++;
    const head = this._header();
    this._head.replaceWith(head);
    this._head = head;
    this.refresh();
  }

  setItems(items: readonly T[] | RowsLike<T> | ReadonlySignal<T[]>): void {
    this._liveNow = this._liveOption ?? !Array.isArray(items);
    if (Array.isArray(items) || items instanceof Signal) {
      this.keyOf = this._optionKeyOf;
      super.setItems(items as readonly T[] | ReadonlySignal<T[]>);
      return;
    }
    // a `RowsLike` brings its own identity — and a later one replaces it, never the first forever
    const rows = items as RowsLike<T>;
    this.keyOf = this._optionKeyOf ?? ((item: T) => rows.keyOf(item));
    super.setItems(rows.items as ReadonlySignal<T[]>);
  }

  /** The rows sit under the sticky header, so the header's height is not theirs. */
  protected get viewport(): number {
    return Math.max(0, this.root.clientHeight - this.itemHeight);
  }

  protected createRow(): HTMLElement {
    const row = document.createElement('div');
    row.className = 'u2-data-table-row';
    row.setAttribute('role', 'row');
    row.style.height = `${this.itemHeight}px`;
    this._cells(row);
    return row;
  }

  private _cells(row: HTMLElement): void {
    row.style.gridTemplateColumns = this._template;
    row.replaceChildren(...this._columns.map((column) => {
      const cell = document.createElement('div');
      cell.setAttribute('role', 'gridcell');
      cell.className = DataTable._cellClass(column);
      VirtualRows.capture(cell, this._cellBaseline);
      return cell;
    }));
    this._rowShape.set(row, this._shape);
  }

  /** A kept row at its new place; live items are rendered again, since their content may have
   * changed under the same object. */
  protected placeRow(row: HTMLElement, index: number): boolean {
    // the header is row 1, so the data rows start at 2
    row.setAttribute('aria-rowindex', String(index + 2));
    return !this._liveNow;
  }

  /** The row's cells rendered in place — each restored to a fresh cell's attributes first, so a
   * `render` that set a title, a dataset key or an aria flag does not follow the pool. */
  protected fillRow(row: HTMLElement, index: number): void {
    if (this._rowShape.get(row) !== this._shape)
      this._cells(row);
    else if (this._shape > 0)
      row.style.gridTemplateColumns = this._template;
    this.placeRow(row, index);
    const item = this.items[index];
    const keyOf = this.keyOf;
    const state = this._cellState;
    const key = state === undefined || keyOf === undefined ? undefined : keyOf(item);
    const cells = row.children;
    for (const [i, column] of this._columns.entries()) {
      const cell = cells[i] as HTMLElement;
      VirtualRows.restore(cell, this._cellBaseline);
      cell.className = DataTable._cellClass(column);
      if (i === 0)
        cell.setAttribute('data-u2-part', 'item-label');
      const content = column.render === undefined ? DataTable._text(item, column.name) :
        column.render(item, index, cell);
      if (typeof content === 'string')
        cell.textContent = content;
      else
        cell.replaceChildren(content);
      if (key === undefined)
        continue;
      if (state!.isChanged(key, column.name))
        cell.classList.add('u2-cell-changed');
      if (state!.canEdit !== undefined && !state!.canEdit(key, column.name))
        cell.classList.add('u2-cell-readonly');
      const error = state!.errorOf(key, column.name);
      if (error === null)
        continue;
      cell.classList.add('u2-cell-error');
      cell.title = error.message;
    }
  }

  protected rowKey(e: KeyboardEvent, index: number): boolean {
    if (e.key !== 'Enter')
      return false;
    if (this._onActivate)
      e.preventDefault();
    this.activate(index);
    return true;
  }

  protected activate(index: number): void {
    this._onActivate?.(this.items[index], index);
  }

  private _header(): HTMLElement {
    const header = document.createElement('div');
    header.className = 'u2-data-table-header';
    header.setAttribute('role', 'row');
    header.style.gridTemplateColumns = this._template;
    header.style.height = `${this.itemHeight}px`;
    for (const column of this._columns) {
      const cell = document.createElement('div');
      cell.className = `${DataTable._cellClass(column)} u2-data-table-head`;
      cell.setAttribute('role', 'columnheader');
      cell.textContent = cell.title = column.header ?? column.name;
      header.append(cell);
    }
    return header;
  }

  private static _templateOf<T>(columns: DataTableColumn<T>[]): string {
    return columns.map((c) => c.width ?? 'minmax(0, 1fr)').join(' ');
  }

  private static _cellClass<T>(column: DataTableColumn<T>): string {
    return column.align === undefined ? 'u2-data-table-cell' :
      `u2-data-table-cell u2-data-table-align-${column.align}`;
  }

  private static _text<T>(item: T, name: string): string {
    const value = (item as Record<string, unknown>)[name];
    return value === null || value === undefined ? '' : String(value);
  }
}
