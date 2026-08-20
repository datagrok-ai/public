/* The DataFrame column selectors — Dart's `ColumnsInput`, `ColumnsMapInput` and the aggregation
   editor behind `AggregatedColumnsInput` — rebuilt on u2 machinery (Overlay + VirtualList + the
   column renderer) instead of wrapping the platform's `selectColumns` modal, which the
   convergence retires. Values are column NAMES, as everywhere in the pickers (pickers.ts). */
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../core/input-base.js';
import {Scope} from '../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../core/overlay.js';
import {button, divH, divV, link, span} from '../core/elements.js';
import {iconButton} from '../components/actions/buttons.js';
import {VirtualList} from '../components/collections/list.js';
import {TextInput} from '../components/inputs/text-input.js';
import {ChoiceInput} from '../components/inputs/choice-input.js';
import {columnInput, onColumnRenamed} from './pickers.js';
import {columnRenderer} from './column-renderer.js';

const ROW_HEIGHT = 24;

/* Aggregations per column type, as ddt registers them — `first` for every type
   (`column_type_meta.dart:64`), `Stats.aggregations` for the numeric metas
   (`float_meta.dart:181` and its int/bigint/qnum twins), `Stats.counts` for bool
   (`bool_meta.dart:27`), counts plus the string functions for string (`string_meta.dart:25-45`),
   counts plus min/max/avg/range for datetime (`date_time_meta.dart:35-42`) — minus the per-type
   exclusions of Dart's aggregation editor (`aggregated_columns_input.dart:8-12`), which is why
   string and bool come down to the counts and datetime keeps only `range` on top of them.
   `geomean` and `range` have no DG.AGG member. Dart applies those exclusions only in its
   multi-property mode (`isMultiPropEditor`, `:56`); here they always apply, which is also the
   branch {@link defaultAggregation} matches (`total count` for a string column, not Dart's
   single-property `value or unique count`, which that same mode excludes). */
const COUNTS: string[] =
  [DG.AGG.TOTAL_COUNT, DG.AGG.VALUE_COUNT, DG.AGG.UNIQUE_COUNT, DG.AGG.MISSING_VALUE_COUNT];
const NUMERIC: string[] = [DG.AGG.FIRST, ...COUNTS, DG.AGG.MIN, DG.AGG.MAX, DG.AGG.SUM, DG.AGG.MED,
  DG.AGG.AVG, 'geomean', DG.AGG.STDEV, DG.AGG.VARIANCE, DG.AGG.SKEW, DG.AGG.KURT,
  DG.AGG.Q1, DG.AGG.Q2, DG.AGG.Q3];
const AGGREGATIONS: {[type: string]: string[]} = {
  [DG.COLUMN_TYPE.INT]: NUMERIC,
  [DG.COLUMN_TYPE.FLOAT]: NUMERIC,
  [DG.COLUMN_TYPE.BIG_INT]: NUMERIC,
  [DG.COLUMN_TYPE.QNUM]: NUMERIC,
  [DG.COLUMN_TYPE.STRING]: COUNTS,
  [DG.COLUMN_TYPE.BOOL]: COUNTS,
  [DG.COLUMN_TYPE.DATE_TIME]: [...COUNTS, 'range'],
};

/** The aggregations offered for a column of that type. Anything unregistered — object, byte
 * array, a virtual column's own type — counts rows and nothing else. */
export function aggregationsFor(type: string): string[] {
  return AGGREGATIONS[type] ?? COUNTS;
}

/** Dart seeds a new row with `avg` for a numerical column and `total count` otherwise in its
 * multi-property mode (`aggregated_columns_input.dart:30`) — but the same file excludes `avg` for
 * datetimes, which the platform also counts as numerical. So the seed is `avg` where the column's
 * type actually offers it, and the first aggregation it does offer otherwise. */
export function defaultAggregation(column: DG.Column): string {
  const available = aggregationsFor(column.type);
  return column.isNumerical && available.includes(DG.AGG.AVG) ? DG.AGG.AVG : available[0];
}

export interface ColumnsInputOptions extends InputOptions<string[]> {
  table: DG.DataFrame;
  filter?: (column: DG.Column) => boolean;
}

/** Picks any number of columns of one table. The editor is the summary Dart shows —
 * `(3) age, sex, race`, or `(N) All` once everything is in (`columns_input.dart:145-156`) — and
 * clicking it opens a searchable check-list with the selected columns on top (`:56-61`), anchored
 * under the control as Dart's own selector is (`selectColumns(… showNextTo: input, maxHeight: 300)`,
 * `columns_input.dart:62-79`). Dismissal follows that modal (`dialog.dart:486-494`): an outside
 * click or Esc cancels, OK commits, and clicking the control again closes it. A column dropped
 * from the table leaves the value with it. */
export class ColumnsInput extends Input<string[], ColumnsInputOptions> {
  private _table!: DG.DataFrame;
  private _filter: ((column: DG.Column) => boolean) | undefined;
  private _control!: HTMLElement;
  private _summary!: HTMLElement;
  private _unfollow: (() => void) | undefined;
  private _popupScope: Scope | undefined;

  constructor(options: ColumnsInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'columns-input';
  }

  get table(): DG.DataFrame {
    return this._table;
  }

  /** Swaps the source table, and the filter with it (omit it to select from every column): the
   * selection goes with the old table, as Dart's `changeTable` does
   * (`columns_input.dart:118-129`). */
  changeTable(table: DG.DataFrame, filter?: (column: DG.Column) => boolean): void {
    this._closePopup();
    this._table = table;
    this._filter = filter;
    this._follow();
    this.value.value = [];
    this._renderSummary(this.value.peek());
  }

  protected createEditor(): HTMLElement {
    this._table = this.options.table;
    this._filter = this.options.filter;
    const editor = document.createElement('div');
    this._control = editor;
    editor.className = 'u2-columns';
    editor.tabIndex = 0;
    editor.setAttribute('role', 'button');
    editor.setAttribute('aria-haspopup', 'listbox');
    editor.setAttribute('aria-expanded', 'false');
    this._summary = span('', 'u2-columns-summary');
    editor.append(this._summary, span('', 'u2-columns-chevron'));
    this._listen(editor, 'click', () => this._togglePopup());
    this._listen(editor, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.effect(() => this._renderSummary(this.value.value));
    this._follow();
    this.own(() => {
      this._closePopup();
      if (this._unfollow)
        this._unfollow();
    });
    return editor;
  }

  private _names(): string[] {
    const columns = this._table.columns;
    return this._filter ? columns.toList().filter(this._filter).map((c) => c.name) : columns.names();
  }

  private _renderSummary(names: string[]): void {
    const available = this._names().length;
    this._summary.textContent = names.length === 0 ? '(0)' :
      names.length === available ? `(${names.length}) All` : `(${names.length}) ${names.join(', ')}`;
  }

  private _follow(): void {
    if (this._unfollow)
      this._unfollow();
    const subs = [this._table.onColumnsChanged.subscribe(() => this._prune()),
      onColumnRenamed(this._table, (oldName, newName) => this._rename(oldName, newName))];
    this._unfollow = () => {
      for (const sub of subs)
        sub.unsubscribe();
    };
  }

  private _rename(oldName: string, newName: string): void {
    const value = this.value.peek();
    const index = value.indexOf(oldName);
    if (index < 0)
      return;
    const renamed = value.slice();
    renamed[index] = newName;
    Input.system(() => this.value.value = renamed);
  }

  private _prune(): void {
    const available = this._names();
    const value = this.value.peek();
    const kept = value.filter((name) => available.includes(name));
    if (kept.length !== value.length)
      Input.system(() => this.value.value = kept);
    else
      this._renderSummary(value); // the 'All' threshold moves with the column count
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'Enter' && e.key !== ' ')
      return;
    e.preventDefault();
    this._togglePopup();
  }

  private _togglePopup(): void {
    if (!this.enabled) // the editor is a div, so the base's disabled sweep never reaches it
      return;
    if (this._popupScope)
      this._closePopup();
    else
      this._openPopup();
  }

  private _openPopup(): void {
    const scope = new Scope();
    this._popupScope = scope;
    Scope.runWith(scope, () => {
      const picked = new Set(this.value.peek());
      const order = this._ordered();
      const renderer = columnRenderer(this._table);
      const search = new TextInput({placeholder: 'Search…', search: true, inline: true});
      const filtered = (query: string): string[] => {
        const q = query.trim().toLowerCase();
        return q ? order.filter((name) => name.toLowerCase().includes(q)) : order;
      };
      const list = new VirtualList<string>({
        itemHeight: ROW_HEIGHT,
        render: (name, _index, row) => {
          row.classList.toggle('u2-columns-picked', picked.has(name));
          const box = document.createElement('input');
          box.type = 'checkbox';
          box.className = 'u2-input-checkbox';
          box.checked = picked.has(name);
          box.tabIndex = -1;
          // named for a screen reader here rather than through a wrapping <label>, whose native
          // activation would toggle the row a second time on top of the list's own click handler
          box.setAttribute('aria-label', name);
          return divH([box, renderer.listItem(name)], 'u2-columns-option');
        },
      });
      list.root.classList.add('u2-columns-list');
      const count = span('', 'u2-columns-count');
      const refresh = (query = search.value.peek()) => {
        list.setItems(filtered(query));
        count.textContent = `${picked.size} selected`;
      };
      const toggle = (index: number) => {
        const name = filtered(search.value.peek())[index];
        if (name === undefined)
          return;
        if (!picked.delete(name))
          picked.add(name);
        refresh();
      };
      // All and None act on what the search left visible, so they read as 'all of these'
      const setAll = (on: boolean) => {
        for (const name of filtered(search.value.peek())) {
          if (on)
            picked.add(name);
          else
            picked.delete(name);
        }
        refresh();
      };
      const ok = () => {
        this._commit(order.filter((name) => picked.has(name)));
        this._closePopup();
      };
      const move = (delta: number) => {
        const rows = filtered(search.value.peek()).length;
        if (rows === 0)
          return;
        const next = Math.min(Math.max(list.selectedIndex.peek() + delta, 0), rows - 1);
        list.selectedIndex.value = next;
        list.scrollToIndex(next);
      };
      const onClick = (e: Event) => {
        const row = (e.target as HTMLElement).closest('.u2-list-row') as HTMLElement | null;
        if (row)
          toggle(Number(row.dataset.index));
      };
      // Space, not Enter: Enter is OK, and the list already owns the arrows
      const onKeyDown = (e: Event) => {
        if ((e as KeyboardEvent).key !== ' ')
          return;
        e.preventDefault();
        toggle(list.selectedIndex.peek());
      };
      list.root.addEventListener('click', onClick);
      list.root.addEventListener('keydown', onKeyDown);
      scope.own(() => {
        list.root.removeEventListener('click', onClick);
        list.root.removeEventListener('keydown', onKeyDown);
      });
      scope.effect(() => refresh(search.value.value));

      const popup = divV([search,
        divH([link('All', () => setAll(true)), link('None', () => setAll(false)), count],
          'u2-columns-actions'),
        list,
        divH([button('OK', ok, {primary: true}), button('CANCEL', () => this._closePopup())],
          'u2-columns-buttons')], 'u2-columns-popup');
      // the list keeps its own keys while it has the focus; the popup covers the rest, so Enter,
      // arrows and Space work from the search box too, where the focus opens
      const onPopupKeyDown = (e: Event) => {
        const key = (e as KeyboardEvent).key;
        const target = e.target as HTMLElement;
        if (target.tagName === 'BUTTON' || target.tagName === 'A') // these answer for themselves
          return;
        if (key === 'Enter') {
          e.preventDefault();
          ok();
          return;
        }
        if (list.root.contains(target))
          return;
        if (key === 'ArrowDown' || key === 'ArrowUp') {
          e.preventDefault();
          move(key === 'ArrowDown' ? 1 : -1);
          return;
        }
        if (key === ' ' && list.selectedIndex.peek() >= 0) {
          e.preventDefault();
          toggle(list.selectedIndex.peek());
        }
      };
      // outside mousedown, Esc and anchor removal all reach the popup as one close event, and all
      // three mean cancel — the value is only ever written by OK (`dialog.dart:486-494`)
      const onDismiss = () => this._closePopup();
      popup.addEventListener('keydown', onPopupKeyDown);
      popup.addEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
      scope.own(() => {
        popup.removeEventListener('keydown', onPopupKeyDown);
        popup.removeEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
      });
      Overlay.show(this._control, popup, scope);
      this._control.setAttribute('aria-expanded', 'true');
      search.root.querySelector<HTMLElement>('input')?.focus();
    });
  }

  /** Selected columns first, in the order they were picked, then the rest in table order — what
   * Dart hands the selector as `order` (`columns_input.dart:56-61`). */
  private _ordered(): string[] {
    const available = this._names();
    const selected = this.value.peek().filter((name) => available.includes(name));
    return [...selected, ...available.filter((name) => !selected.includes(name))];
  }

  private _commit(names: string[]): void {
    const value = this.value.peek();
    if (names.length !== value.length || names.some((name, i) => name !== value[i]))
      this.value.value = names;
  }

  private _closePopup(): void {
    const scope = this._popupScope;
    this._popupScope = undefined;
    if (!scope)
      return;
    scope.dispose();
    this._control.setAttribute('aria-expanded', 'false');
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}

export function columnsInput(label: string, table: DG.DataFrame,
  options: Omit<ColumnsInputOptions, 'table'> = {}): ColumnsInput {
  return new ColumnsInput({...options, label, table});
}

export interface ColumnKey {
  name: string;
  /** Restricts the row to columns of that type, joining int and double as Dart's `_joinNumType`
   * does (`columns_map_input.dart:120`). */
  type?: string;
}

export interface ColumnsMapInputOptions extends InputOptions<Record<string, string>> {
  table: DG.DataFrame;
  keys: (string | ColumnKey)[];
}

/** Maps each key to a column of one table — Dart's `ColumnsMapInput` (`columns_map_input.dart`)
 * with its modal unrolled into a row per key, since there is nothing else in that modal. Unmapped
 * keys stay out of the value (Dart keeps them as nulls and labels the button `(mapped/total)`). */
export class ColumnsMapInput extends Input<Record<string, string>, ColumnsMapInputOptions> {
  private _pickers!: {key: string, input: ChoiceInput}[];

  constructor(options: ColumnsMapInputOptions) {
    super(options, {});
    this.root.dataset.u2 = 'columns-map-input';
  }

  input(key: string): ChoiceInput | undefined {
    return this._pickers.find((p) => p.key === key)?.input;
  }

  protected createEditor(): HTMLElement {
    const table = this.options.table;
    const body = document.createElement('div');
    body.className = 'u2-columns-map';
    this._pickers = this.options.keys.map((key) => {
      const {name, type} = typeof key === 'string' ? {name: key, type: undefined} : key;
      const input = columnInput(name, table,
        {filter: type === undefined ? undefined : (c) => joinNumType(c.type) === joinNumType(type)});
      body.append(input.root);
      return {key: name, input};
    });
    this._load(this.value.peek());
    this.effect(() => this._commit(this._pickers.map((p) => p.input.value.value)));
    this.effect(() => this._load(this.value.value));
    return body;
  }

  private _load(value: Record<string, string>): void {
    for (const picker of this._pickers)
      picker.input.value.value = value[picker.key] ?? null;
  }

  /** Rewrites the value from the pickers, and only when it actually changed — the two effects
   * would otherwise bounce a new object off each other on every load. */
  private _commit(picked: (string | null)[]): void {
    const value: Record<string, string> = {};
    for (let i = 0; i < this._pickers.length; i++) {
      const name = picked[i];
      if (name !== null)
        value[this._pickers[i].key] = name;
    }
    const current = this.value.peek();
    const keys = Object.keys(value);
    if (keys.length === Object.keys(current).length && keys.every((key) => current[key] === value[key]))
      return;
    this.value.value = value;
  }
}

export function columnsMapInput(label: string, keys: (string | ColumnKey)[], table: DG.DataFrame,
  options: Omit<ColumnsMapInputOptions, 'table' | 'keys'> = {}): ColumnsMapInput {
  return new ColumnsMapInput({...options, label, table, keys});
}

export interface ColumnAggregation {
  column: string;
  aggregation: string;
}

export interface AggregatedColumnsInputOptions extends InputOptions<ColumnAggregation[]> {
  table: DG.DataFrame;
}

/** Column + aggregation rows with per-row add and remove — Dart's `editColumnAggregations`
 * (`aggregated_columns_input.dart:35-62`) as an input instead of a modal. The aggregation list
 * follows the column's type ({@link aggregationsFor}); a column change that invalidates the row's
 * aggregation resets it to the type's default, as Dart does (`:49-52`). */
export class AggregatedColumnsInput extends Input<ColumnAggregation[], AggregatedColumnsInputOptions> {
  private _rows!: ColumnAggregation[];
  private _body!: HTMLElement;
  // the row editors are rebuilt whenever rows are added or removed; each generation owns its
  // inputs and tooltips, and is released before the next one is built
  private _rowScope: Scope | undefined;

  constructor(options: AggregatedColumnsInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'aggregated-columns-input';
  }

  protected createEditor(): HTMLElement {
    const body = document.createElement('div');
    this._body = body;
    body.className = 'u2-aggregations';
    this._rows = AggregatedColumnsInput._copy(this.value.peek());
    this._render();
    // echo suppression by comparison, not by a flag: a row editor writes the value from inside its
    // own effect, and this one is flushed only after that write returns
    this.effect(() => {
      const value = this.value.value;
      if (AggregatedColumnsInput._same(value, this._rows))
        return;
      this._rows = AggregatedColumnsInput._copy(value);
      this._render();
    });
    // the rows name their column, and so do the choices they offer: both are rebuilt on a rename
    const renamed = onColumnRenamed(this.options.table, (oldName, newName) => {
      for (const row of this._rows) {
        if (row.column === oldName)
          row.column = newName;
      }
      this._render();
      Input.system(() => this._commit());
    });
    this.own(() => {
      renamed.unsubscribe();
      this._rowScope?.dispose();
    });
    return body;
  }

  private _render(): void {
    this._rowScope?.dispose();
    const scope = new Scope();
    this._rowScope = scope;
    this._body.textContent = '';
    Scope.runWith(scope, () => {
      if (this._rows.length === 0)
        this._body.append(divH([this._addButton(0)], 'u2-aggregations-row'));
      for (let i = 0; i < this._rows.length; i++)
        this._body.append(this._row(i));
    });
    this.refreshEnabled();
  }

  private _row(index: number): HTMLElement {
    const row = this._rows[index];
    const aggregation = new ChoiceInput({inline: true, nullable: false,
      items: aggregationsFor(this._type(row.column)), value: row.aggregation,
      onChanged: (value) => {
        row.aggregation = value ?? row.aggregation;
        this._commit();
      }});
    const column = new ChoiceInput({inline: true, nullable: false,
      items: this.options.table.columns.names(), value: row.column,
      onChanged: (value) => this._onColumn(row, value!, aggregation)});
    return divH([column.root, aggregation.root, this._addButton(index + 1),
      iconButton('minus', () => this._remove(index), {tooltip: 'Remove aggregation'})],
    'u2-aggregations-row');
  }

  private _addButton(at: number): HTMLElement {
    return iconButton('plus', () => this._add(at), {tooltip: 'Add aggregation'});
  }

  /** A column change re-offers the aggregations its type has, and resets the row's own when the
   * new type does not offer it (`aggregated_columns_input.dart:49-52`). */
  private _onColumn(row: ColumnAggregation, name: string, aggregation: ChoiceInput): void {
    row.column = name;
    const available = aggregationsFor(this._type(name));
    if (!available.includes(row.aggregation))
      row.aggregation = defaultAggregation(this.options.table.columns.byName(name));
    aggregation.setItems(available);
    aggregation.value.value = row.aggregation;
    this._commit();
  }

  private _type(name: string): string {
    const columns = this.options.table.columns;
    return columns.contains(name) ? columns.byName(name).type : '';
  }

  private _add(at: number): void {
    const row = this._newRow();
    if (row === null)
      return;
    this._rows.splice(at, 0, row);
    this._render();
    this._commit();
  }

  private _remove(index: number): void {
    this._rows.splice(index, 1);
    this._render();
    this._commit();
  }

  /** The first numerical column, else the first one — Dart's `_newRow`
   * (`aggregated_columns_input.dart:28-30`). Null for a table with no columns. */
  private _newRow(): ColumnAggregation | null {
    const columns = this.options.table.columns.toList();
    const column = columns.find((c) => c.isNumerical) ?? columns[0];
    return column === undefined ? null : {column: column.name, aggregation: defaultAggregation(column)};
  }

  private _commit(): void {
    if (!AggregatedColumnsInput._same(this.value.peek(), this._rows))
      this.value.value = AggregatedColumnsInput._copy(this._rows);
  }

  private static _copy(rows: ColumnAggregation[]): ColumnAggregation[] {
    return rows.map((row) => ({...row}));
  }

  private static _same(a: ColumnAggregation[], b: ColumnAggregation[]): boolean {
    return a.length === b.length &&
      a.every((row, i) => row.column === b[i].column && row.aggregation === b[i].aggregation);
  }
}

export function aggregatedColumnsInput(label: string, table: DG.DataFrame,
  options: Omit<AggregatedColumnsInputOptions, 'table'> = {}): AggregatedColumnsInput {
  return new AggregatedColumnsInput({...options, label, table});
}

function joinNumType(type: string): string {
  return type === DG.COLUMN_TYPE.INT ? DG.COLUMN_TYPE.FLOAT : type;
}
