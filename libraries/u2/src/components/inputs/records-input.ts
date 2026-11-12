/* RecordsInput: an editable list of typed records — a row per record, an inline u2 input per field under a
   header naming the fields, a remove button per row and an add link under them; `orderable`, a drag handle per
   row (Alt+↑/↓ from the keyboard). Rows are keyed (`each`), so an edit in one never rebuilds another. */
import {Input, InputOptions} from '../../core/input-base.js';
import {signal, computed, untracked, afterFlush, Signal, ReadonlySignal} from '../../core/signals.js';
import {div, span, link} from '../../core/elements.js';
import {each} from '../../core/region.js';
import {icon} from '../display/icon.js';
import {iconButton} from '../actions/buttons.js';
import {Reorder} from './reorder.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface RecordsColumn<T> {
  /** The record field the column edits. */
  name: keyof T & string;
  /** The column's header; the field name by default. */
  header?: string;
  /** A fresh input for one cell, started on the field's value — inline (`inline: true`): the header is its
   * caption. The record is there for a cell that depends on its row. */
  input: (value: any, record: T) => Input<any>;
  /** The column's grid track: `'1fr'` (the default), `'120px'`. */
  width?: string;
}

export interface RecordsInputOptions<T> extends InputOptions<T[]> {
  columns: RecordsColumn<T>[];
  /** A new row's record; `{}` by default. */
  newItem?: () => T;
  /** The add link's text ('Add field'); 'Add' by default. */
  addText?: string;
  /** Rows have a drag handle and move with Alt+↑/↓. */
  orderable?: boolean;
  /** The add link (default true). */
  allowAdd?: boolean;
  /** A remove button per row (default true). */
  allowRemove?: boolean;
  /** What a row is called where its verdict is told ("Stage 3: Value can't be empty"); `Row <n>` by default. */
  rowName?: (record: T, index: number) => string;
  /** The key a row answers to in the items block; the first column's value by default. */
  keyOf?: (record: T, index: number) => string;
}

interface Row<T> {
  id: number;
  record: T;
  cells: Input<any>[];
  el?: HTMLElement;
}

export class RecordsInput<T extends object> extends Input<T[], RecordsInputOptions<T>> {
  // built by createEditor(), which the base constructor calls before subclass initializers would run
  private declare _rows: Signal<Row<T>[]>;
  /** Bumped when a row's cells are built: the verdict reads the cells, which are no signal. */
  private declare _built: Signal<number>;
  private declare _problem: ReadonlySignal<string | null>;
  private declare _body: HTMLElement;
  private declare _nextId: number;

  constructor(options: RecordsInputOptions<T>) {
    super(options, []);
    this.root.dataset.u2 = 'records-input';
  }

  /** The rows as records, in their order. */
  get records(): T[] {
    return this._rows.peek().map((row) => row.record);
  }

  get readOnly(): boolean {
    return super.readOnly;
  }

  /** Every cell's too, and a row added later's. */
  set readOnly(x: boolean) {
    super.readOnly = x;
    for (const cell of this._cells())
      cell.readOnly = x;
  }

  get disabled(): string | null {
    return super.disabled;
  }

  set disabled(reason: string | null) {
    super.disabled = reason;
    for (const cell of this._cells())
      cell.disabled = reason;
  }

  /** The rows, keyed by {@link RecordsInputOptions.keyOf}; a row with a failing cell is `invalid`. */
  getItems(opts?: {all?: boolean}): IItemsState {
    return Input.itemsState(this._rows.peek().map((row, index) => ({key: this._keyOf(row.record, index),
      label: RecordsInput._text(row), index, rendered: true,
      ...row.cells.some((cell) => cell.validity.peek() !== null) ? {invalid: true} : {}})), opts?.all);
  }

  /** Adds a row — `record`, or a new one — at `index` (the end by default), and focuses its first cell. */
  add(record: T = this.options.newItem?.() ?? {} as T, index?: number): void {
    const rows = this._rows.peek().slice();
    const row: Row<T> = {id: this._nextId++, record, cells: []};
    rows.splice(index ?? rows.length, 0, row);
    this._rows.value = rows;
    this._commit();
    // inside a handler's batch the row is built when the batch ends
    const focus = () => row.cells[0]?.focus();
    if (row.cells.length > 0)
      focus();
    else
      afterFlush(this.scope.bind(focus));
  }

  /** Removes the row at `index`. */
  remove(index: number): void {
    const rows = this._rows.peek().slice();
    rows.splice(index, 1);
    this._rows.value = rows;
    this._commit();
  }

  /** Moves the row at `from` to `to`. */
  move(from: number, to: number): void {
    this._place(from, to);
    this._commit();
  }

  protected displayText(): string {
    return this._rows.peek().map(RecordsInput._text).join(', ');
  }

  protected createEditor(): HTMLElement {
    const o = this.options;
    this._nextId = 0;
    this._rows = signal([]);
    this._built = signal(0);
    const editor = div([], 'u2-records');
    editor.setAttribute('role', 'grid');
    const tracks = [o.orderable ? '16px' : null, ...o.columns.map((c) => c.width ?? '1fr'),
      o.allowRemove === false ? null : '24px'].filter((t) => t !== null).join(' ');
    editor.style.setProperty('--u2-records-tracks', tracks);

    const head = div([], 'u2-records-row u2-records-head');
    head.setAttribute('role', 'row');
    if (o.orderable)
      head.append(span(''));
    for (const c of o.columns) {
      const h = span(c.header ?? c.name, 'u2-records-header');
      h.setAttribute('role', 'columnheader');
      head.append(h);
    }
    if (o.allowRemove !== false)
      head.append(span(''));
    this._body = div([], 'u2-records-body');
    editor.append(head, this._body);
    this.effect(() => {
      head.hidden = this._rows.value.length === 0;
    });
    each(this._body, this._rows, (row) => String(row.id), (row) => this._row(row));

    if (o.allowAdd !== false) {
      const add = link(o.addText ?? 'Add', () => this.add());
      add.prepend(icon('plus', {cls: 'u2-records-add-icon'}));
      add.classList.add('u2-records-add');
      add.dataset.u2Part = 'add';
      editor.append(add);
    }

    this._problem = computed(() => {
      this._built.value;
      const rows = this._rows.value;
      for (let i = 0; i < rows.length; i++) {
        for (let c = 0; c < rows[i].cells.length; c++) {
          const message = rows[i].cells[c].validity.value;
          if (message !== null) {
            const header = o.columns.length > 1 ? `${o.columns[c].header ?? o.columns[c].name}: ` : '';
            return `${o.rowName?.(rows[i].record, i) ?? `Row ${i + 1}`}: ${header}${message}`;
          }
        }
      }
      return null;
    });
    // registered here, before the base builds its validation effect, so the effect tracks the cells
    this.addValidator(() => this._problem.value);
    // our own edits read back as the rows they came from: only another writer reloads them
    this.effect(() => {
      const value = this.value ?? [];
      const rows = untracked(() => this._rows.value);
      if (value.length !== rows.length || value.some((record, i) => record !== rows[i].record))
        untracked(() => this._rows.value = value.map((record) => ({id: this._nextId++, record, cells: []})));
    });
    return editor;
  }

  protected resetEditor(): void {
    for (const row of this._rows.peek()) {
      for (const cell of row.cells)
        cell.reset(cell.peek());
    }
  }

  private _row(row: Row<T>): HTMLElement {
    const o = this.options;
    const el = row.el = div([], 'u2-records-row');
    el.setAttribute('role', 'row');
    // the key the items block gives the row, not the row's own id `each()` keys it by
    el.dataset.u2Item = this._keyOf(row.record, this._rows.peek().indexOf(row));
    let handle: HTMLElement | undefined;
    if (o.orderable) {
      handle = Reorder.handle();
      el.append(handle);
    }
    row.cells = o.columns.map((column) => {
      const cell = column.input(row.record[column.name], row.record);
      // fills its column (`fill: true`): a form's widths for its own choices and numbers stay off the cells
      cell.root.classList.add('u2-records-cell', 'u2-input-fill');
      if (cell.label === '')
        cell.label = column.header ?? column.name;
      // the record stays as given until the cell is edited: an editor's own empty value ('' for undefined) is no edit
      let initial = true;
      cell.effect(() => {
        const value = cell.value;
        if (initial || Object.is(value, row.record[column.name])) {
          initial = false;
          return;
        }
        row.record = {...row.record, [column.name]: value};
        untracked(() => this._commit());
      });
      if (this.readOnly)
        cell.readOnly = true;
      if (this.disabled !== null)
        cell.disabled = this.disabled;
      el.append(cell.root);
      return cell;
    });
    if (o.allowRemove !== false) {
      const remove = iconButton('minus', () => this.remove(this._rows.peek().indexOf(row)), {tooltip: 'Remove'});
      remove.classList.add('u2-records-remove');
      el.append(remove);
    }
    if (handle !== undefined) {
      Reorder.wire(el, handle, () => Array.from(this._body.children) as HTMLElement[],
        (from, to) => this._place(from, to), () => this._commit());
    }
    this._built.value++;
    return el;
  }

  private _cells(): Input<any>[] {
    return this._rows?.peek().flatMap((row) => row.cells) ?? [];
  }

  /** Moves a row without committing — a drag commits once, at its end. */
  private _place(from: number, to: number): void {
    const rows = this._rows.peek().slice();
    const [row] = rows.splice(from, 1);
    rows.splice(to, 0, row);
    this._rows.value = rows;
  }

  /** The value from the rows, written only when it differs from what the value holds. */
  private _commit(): void {
    const rows = this._rows.peek();
    for (let i = 0; i < rows.length; i++) {
      if (rows[i].el !== undefined)
        rows[i].el!.dataset.u2Item = this._keyOf(rows[i].record, i);
    }
    const records = this.records;
    const value = this.peek() ?? [];
    if (records.length !== value.length || records.some((record, i) => record !== value[i]))
      this.value = records;
  }

  private _keyOf(record: T, index: number): string {
    const key = this.options.keyOf?.(record, index) ?? record[this.options.columns[0]?.name];
    return key === undefined || key === null || key === '' ? `#${index}` : String(key);
  }

  private static _text(row: Row<object>): string {
    return row.cells.map((cell) => cell.getState().text ?? '').filter((t) => t !== '').join(' ');
  }
}
