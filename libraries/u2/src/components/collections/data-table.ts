/* Virtualized data table — `VirtualRows`' scroller and selection model with a row shape of pooled
   cells under a sticky header, plus the editor's per-cell verdicts. Cells are pooled with their
   row, so scrolling a wide table costs the visible window and nothing else. `BasicTable` stays the
   small-data control: every row in the DOM, real table semantics. */

import {Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {icon} from '../display/icon.js';
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
  /** Row height in pixels; 24 by default. */
  rowHeight?: number;
  /** The order (`'name'`, `'!name'` descending, `''` none) the headers show and set by click;
   * absent, the headers do not sort. */
  sort?: Signal<string>;
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

const MAX_HEADER_HEIGHT = 28;
const MIN_COLUMN_WIDTH = 32;

export class DataTable<T> extends VirtualRows<T> {
  private _columns: DataTableColumn<T>[];
  private readonly _cellState: CellStateLike | undefined;
  private readonly _onActivate: ((item: T, index: number) => void) | undefined;
  private readonly _optionKeyOf: ((item: T) => string) | undefined;
  private readonly _liveOption: boolean | undefined;
  private _template: string;
  private _head: HTMLElement;
  private readonly _headerHeight: number;
  private readonly _sort: Signal<string> | undefined;
  private _endResize: (() => void) | undefined;
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
    this._headerHeight = Math.min(this.itemHeight, MAX_HEADER_HEIGHT);
    this._sort = options.sort;
    this._template = DataTable._templateOf(options.columns);
    this._head = this._header();
    this.root.prepend(this._head);
    this._listen();
    const sort = this._sort;
    if (sort !== undefined)
      this.effect(() => this._showSort(sort.value));
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
    return Math.max(0, this.root.clientHeight - this._headerHeight);
  }

  private _resizeColumn(index: number, px: number): void {
    const width = `${Math.max(MIN_COLUMN_WIDTH, Math.round(px))}px`;
    this._columns = this._columns.map((c, i) => i === index ? {...c, width} : c);
    this._template = DataTable._templateOf(this._columns);
    this._head.style.gridTemplateColumns = this._template;
    for (const row of this.content.querySelectorAll('.u2-data-table-row'))
      (row as HTMLElement).style.gridTemplateColumns = this._template;
  }

  /** Header sorting, column resizing, and a cut cell's full text as its tooltip. */
  private _listen(): void {
    const scope = this.scope;
    scope.listen(this.root, 'click', (e) => {
      const target = e.target as Element;
      const head = target.closest('.u2-data-table-head') as HTMLElement | null;
      if (head !== null && this._sort !== undefined && target.closest('.u2-data-table-resize') === null)
        this._sort.value = DataTable._nextSort(this._sort.peek(), head.dataset.column!);
    });
    scope.listen(this.root, 'pointerdown', (e) => {
      const handle = (e.target as Element).closest('.u2-data-table-resize');
      if (handle !== null)
        this._startResize(e as PointerEvent, handle.parentElement!);
    });
    scope.listen(this.content, 'mouseover', (e) => {
      const cell = (e.target as Element).closest('.u2-data-table-cell') as HTMLElement | null;
      if (cell !== null && cell.title === '' && cell.scrollWidth > cell.clientWidth)
        cell.title = cell.textContent ?? '';
    });
    this.own(() => this._endResize?.());
  }

  private _startResize(e: PointerEvent, head: HTMLElement): void {
    e.preventDefault();
    this._endResize?.();
    const index = Array.prototype.indexOf.call(this._head.children, head);
    const from = e.clientX;
    const width = head.getBoundingClientRect().width;
    const move = (m: Event) => this._resizeColumn(index, width + (m as PointerEvent).clientX - from);
    const end = () => {
      document.removeEventListener('pointermove', move);
      document.removeEventListener('pointerup', end);
      this._endResize = undefined;
    };
    document.addEventListener('pointermove', move);
    document.addEventListener('pointerup', end);
    this._endResize = end;
  }

  /** Ascending, then descending, then none. */
  private static _nextSort(sort: string, column: string): string {
    const first = sort.split(',')[0];
    return first === column ? `!${column}` : first === `!${column}` ? '' : column;
  }

  private _showSort(sort: string, header = this._head): void {
    const first = sort.split(',')[0];
    const descending = first.startsWith('!');
    const column = descending ? first.slice(1) : first;
    for (const head of header.children) {
      const cell = head as HTMLElement;
      const sorted = cell.dataset.column === column;
      cell.setAttribute('aria-sort', !sorted ? 'none' : descending ? 'descending' : 'ascending');
      cell.querySelector('.u2-data-table-sort')?.remove();
      if (sorted) {
        cell.insertBefore(icon(descending ? 'arrow-down' : 'arrow-up', {cls: 'u2-data-table-sort'}),
          cell.querySelector('.u2-data-table-resize'));
      }
    }
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
    else
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
    header.style.height = `${this._headerHeight}px`;
    this.root.style.setProperty('--u2-data-table-header-height', `${this._headerHeight}px`);
    header.classList.toggle('u2-data-table-sortable', this._sort !== undefined);
    for (const column of this._columns) {
      const cell = document.createElement('div');
      cell.className = `${DataTable._cellClass(column)} u2-data-table-head`;
      cell.setAttribute('role', 'columnheader');
      cell.dataset.column = column.name;
      const label = document.createElement('span');
      label.className = 'u2-data-table-head-label';
      label.textContent = cell.title = column.header ?? column.name;
      const resize = document.createElement('div');
      resize.className = 'u2-data-table-resize';
      cell.append(label, resize);
      header.append(cell);
    }
    if (this._sort !== undefined)
      this._showSort(this._sort.peek(), header);
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
