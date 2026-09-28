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
}

export class DataTable<T> extends VirtualRows<T> {
  private readonly _columns: DataTableColumn<T>[];
  private readonly _cellState: CellStateLike | undefined;
  private readonly _onActivate: ((item: T, index: number) => void) | undefined;
  private readonly _optionKeyOf: ((item: T) => string) | undefined;
  private readonly _template: string;
  private readonly _cellBaseline = new Map<string, string>();

  constructor(options: DataTableOptions<T>) {
    super({...options, itemHeight: options.rowHeight ?? 24},
      {host: 'u2-data-table', u2: 'data-table', role: 'grid'});
    this._columns = options.columns;
    this._cellState = options.cellState;
    this._onActivate = options.onActivate;
    this._optionKeyOf = options.keyOf;
    this._template = options.columns.map((c) => c.width ?? 'minmax(0, 1fr)').join(' ');
    this.root.prepend(this._header());
    if (options.items)
      this.setItems(options.items);
  }

  setItems(items: readonly T[] | RowsLike<T> | ReadonlySignal<T[]>): void {
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
    row.style.gridTemplateColumns = this._template;
    for (const column of this._columns) {
      const cell = document.createElement('div');
      cell.setAttribute('role', 'gridcell');
      cell.className = DataTable._cellClass(column);
      row.append(cell);
      VirtualRows.capture(cell, this._cellBaseline);
    }
    return row;
  }

  /** The row's cells rendered in place — each restored to a fresh cell's attributes first, so a
   * `render` that set a title, a dataset key or an aria flag does not follow the pool. */
  protected fillRow(row: HTMLElement, index: number): void {
    // the header is row 1, so the data rows start at 2
    row.setAttribute('aria-rowindex', String(index + 2));
    const item = this.items[index];
    const keyOf = this.keyOf;
    const state = this._cellState;
    const key = state === undefined || keyOf === undefined ? undefined : keyOf(item);
    const cells = row.children;
    for (const [i, column] of this._columns.entries()) {
      const cell = cells[i] as HTMLElement;
      VirtualRows.restore(cell, this._cellBaseline);
      cell.className = DataTable._cellClass(column);
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

  private static _cellClass<T>(column: DataTableColumn<T>): string {
    return column.align === undefined ? 'u2-data-table-cell' :
      `u2-data-table-cell u2-data-table-align-${column.align}`;
  }

  private static _text<T>(item: T, name: string): string {
    const value = (item as Record<string, unknown>)[name];
    return value === null || value === undefined ? '' : String(value);
  }
}
