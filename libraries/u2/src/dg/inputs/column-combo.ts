/* The single-column combo with the platform ColumnsGrid dropdown — Dart's `ColumnComboBox`
   (`column_combo_box.dart`) rebuilt on the u2 input machinery. The value is the column NAME
   (pickers.ts convention); the popup is the REAL Dart grid (`DG.ColumnGrid.popup`, anchored in
   u2's Overlay) where the platform is present, and an Overlay + VirtualList over
   {@link columnRenderer} rows everywhere else — chosen once per open, so headless hosts never
   touch the platform. */
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../../core/input-base.js';
import {Scope} from '../../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {div, divV, span} from '../../core/elements.js';
import {VirtualList} from '../../components/collections/list.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {onColumnRenamed} from './pickers.js';
import {columnRenderer} from '../entities/column-renderer.js';

const ROW_HEIGHT = 24;

/** The platform door, structural: the dg-stub has no `ColumnGrid`, and a vendored datagrok-api
 * may lack the typing. Asked per member, since each backend rides a different constructor. */
export function platformColumnGrid(member: 'popup' | 'columnSelector'):
    {[k: string]: (table: unknown, options: unknown) => any} | null {
  const grid = (DG as unknown as {ColumnGrid?: Record<string, unknown>}).ColumnGrid;
  return grid != null && typeof grid[member] === 'function' ?
    grid as {[k: string]: (table: unknown, options: unknown) => any} : null;
}

export interface ColumnInputOptions2 extends InputOptions<string | null> {
  table: DG.DataFrame | null;
  filter?: (column: DG.Column) => boolean;
  placeholder?: string;
}

/** Picks one column of one table by name. The editor is the ColumnsInput control shape (value +
 * chevron); clicking it drops the platform's ColumnsGrid down under the control when
 * `DG.ColumnGrid.popup` exists, and the u2-native searchable list otherwise. The input follows
 * its table itself: a dropped column clears the value ({@link Input.system}), a renamed one
 * carries it over. Closed, the arrows step through the passing columns as Dart's combo does
 * (`column_combo_box.dart:203-220`) — no mouse-wheel cycling and no hover-preview writes
 * (divergences #17/#14). */
export class ColumnInput extends Input<string | null, ColumnInputOptions2> {
  private _table!: DG.DataFrame | null;
  private _filter: ((column: DG.Column) => boolean) | undefined;
  private _control!: HTMLElement;
  private _summary!: HTMLElement;
  private _unfollow: (() => void) | undefined;
  private _popupScope: Scope | undefined;

  constructor(options: ColumnInputOptions2) {
    super(options, null);
    this.root.dataset.u2 = 'column-combo';
  }

  get table(): DG.DataFrame | null {
    return this._table;
  }

  /** The filtered column names of the bound table. */
  names(): string[] {
    if (this._table === null)
      return [];
    const columns = this._table.columns;
    return this._filter ? columns.toList().filter(this._filter).map((c) => c.name) : columns.names();
  }

  /** Swaps the source table, and the filter with it: closes the popup, retargets the follow
   * subscriptions and CLEARS the value — the ColumnsInput `changeTable` semantics. A null table
   * is the inert empty state. */
  changeTable(table: DG.DataFrame | null, filter?: (column: DG.Column) => boolean): void {
    this._closePopup();
    this._table = table;
    this._filter = filter;
    this._reflectTable();
    this._follow();
    this.value.value = null;
  }

  protected createEditor(): HTMLElement {
    this._table = this.options.table ?? null;
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
    this._reflectTable();
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

  private _reflectTable(): void {
    const off = this._table === null;
    this._control.setAttribute('aria-disabled', String(off));
    if (off)
      this.box.title = 'Select a table first';
    else
      this.box.removeAttribute('title');
  }

  private _renderSummary(value: string | null): void {
    this._summary.textContent = value ?? this.options.placeholder ?? '';
  }

  private _follow(): void {
    if (this._unfollow)
      this._unfollow();
    if (this._table === null) {
      this._unfollow = undefined;
      return;
    }
    const subs = [
      this._table.onColumnsChanged.subscribe(() => this._prune()),
      onColumnRenamed(this._table, (oldName, newName) => {
        if (this.value.peek() === oldName)
          Input.system(() => this.value.value = newName);
      }),
    ];
    this._unfollow = () => {
      for (const sub of subs)
        sub.unsubscribe();
    };
  }

  private _prune(): void {
    const value = this.value.peek();
    if (value !== null && !this.names().includes(value))
      Input.system(() => this.value.value = null);
  }

  /** A pick is a USER edit: the signal write plus the native `change` a select would have
   * dispatched — what the auto-badge's user-edit detection listens for. The `change` fires even
   * for a same-value pick: confirming an auto-guess is still a user touch. */
  private _commit(name: string | null): void {
    if (name !== this.value.peek())
      this.value.value = name;
    this._control.dispatchEvent(new Event('change', {bubbles: true}));
  }

  /** ArrowUp/Left previous passing column, ArrowDown/Right next (`column_combo_box.dart:161-165`):
   * no wrap, a null value steps forward into the first column and never backward. */
  private _step(delta: number): void {
    const names = this.names();
    if (names.length === 0)
      return;
    const value = this.value.peek();
    const idx = value === null ? -1 : names.indexOf(value);
    if (delta < 0 && idx === -1)
      return;
    const next = idx + delta;
    if (next < 0 || next >= names.length)
      return;
    this._commit(names[next]);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      this._togglePopup();
      return;
    }
    if (this._popupScope)
      return;
    if (e.key === 'ArrowUp' || e.key === 'ArrowLeft' || e.key === 'ArrowDown' || e.key === 'ArrowRight') {
      e.preventDefault();
      this._step(e.key === 'ArrowUp' || e.key === 'ArrowLeft' ? -1 : 1);
      return;
    }
    if (e.key.length === 1 && !e.ctrlKey && !e.metaKey && !e.altKey) {
      e.preventDefault();
      this._openPopup(e.key);
    }
  }

  private _togglePopup(): void {
    if (!this.enabled || this._table === null) // the editor is a div, so the base's disabled sweep never reaches it
      return;
    if (this._popupScope)
      this._closePopup();
    else
      this._openPopup();
  }

  private _openPopup(seed?: string): void {
    if (!this.enabled || this._table === null || this._popupScope)
      return;
    const scope = new Scope();
    this._popupScope = scope;
    Scope.runWith(scope, () => {
      const platform = platformColumnGrid('popup');
      if (platform !== null)
        this._openPlatform(scope, platform, seed);
      else
        this._openNative(scope, seed);
    });
    this._control.setAttribute('aria-expanded', 'true');
  }

  /** The real Dart ColumnsGrid in the u2 Overlay: search visible from open (divergence #15),
   * the current column selection-highlighted, autoSize AFTER attach, a fresh popup per open with
   * `cg.close()` on dismissal. */
  private _openPlatform(scope: Scope,
    platform: {[k: string]: (table: unknown, options: unknown) => any}, seed?: string): void {
    const cg = platform.popup(this._table, {addEmpty: this.options.nullable !== false, filter: this._filter});
    const host = div([], 'u2-column-combo-host');
    host.append(cg.root);
    const onDismiss = () => this._closePopup();
    const onKeyDown = (e: Event) => this._onPopupKey(e as KeyboardEvent, cg, host);
    // capture: the Dart search box swallows a bubble-phase Esc even when empty (text_input.dart:148)
    const onEscCapture = (e: Event) => {
      if ((e as KeyboardEvent).key !== 'Escape')
        return;
      e.stopPropagation();
      this._closePopup();
      this._control.focus();
    };
    host.addEventListener('keydown', onEscCapture, true);
    host.addEventListener('keydown', onKeyDown);
    host.addEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
    scope.own(() => {
      host.removeEventListener('keydown', onEscCapture, true);
      host.removeEventListener('keydown', onKeyDown);
      host.removeEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
      cg.close();
    });
    Overlay.show(this.box, host, scope);
    try {
      cg.grid.props.backColor = 0;
      cg.grid.props.rowHeaderBackColor = 0;
    } catch { /* skinning is best-effort */ }
    cg.grid.autoSize(400, document.documentElement.clientHeight - 80);
    const width = this.box.getBoundingClientRect().width;
    if (cg.grid.root.getBoundingClientRect().width < width)
      cg.grid.root.style.width = `${width}px`;
    const current = this.value.peek();
    if (current !== null) {
      try {
        cg.dfColumns.selection.init((i: number) => cg.nameCol.get(i) === current);
      } catch { /* highlight is best-effort */ }
    }
    cg.showSearch = true;
    const sub = cg.dfColumns.onCurrentRowChanged.subscribe(() => {
      this._commit(cg.currentColumn?.name ?? null);
      this._closePopup();
    });
    scope.own(() => sub.unsubscribe());
    const search = host.querySelector('input') as HTMLInputElement | null;
    if (search !== null) {
      if (seed !== undefined) {
        search.value = seed;
        search.dispatchEvent(new Event('input', {bubbles: true}));
      }
      search.focus();
    }
  }

  /** The Dart popup's keyboard (`column_combo_box.dart:311-326`, `:361-365`): arrows walk the
   * filter-visible rows, Enter picks the walked row — or the grid's unique passing column
   * (`passesFilter` owns the verdict: name, friendlyName and the grid's own filters). */
  private _onPopupKey(e: KeyboardEvent, cg: any, host: HTMLElement): void {
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
      e.preventDefault();
      const df = cg.dfColumns;
      df.mouseOverRowIdx = e.key === 'ArrowDown' ?
        df.filter.findNext(df.mouseOverRowIdx, true) : df.filter.findPrev(df.mouseOverRowIdx, true);
      return;
    }
    if (e.key !== 'Enter')
      return;
    e.preventDefault();
    const df = cg.dfColumns;
    if (df.mouseOverRowIdx > -1) {
      df.currentRowIdx = df.mouseOverRowIdx; // fires onCurrentRowChanged — the pick path
      return;
    }
    const text = (host.querySelector('input') as HTMLInputElement | null)?.value.trim() ?? '';
    if (text === '')
      return;
    const matches = this._table!.columns.toList().filter((c) => cg.passesFilter(c));
    if (matches.length === 1) {
      this._commit(matches[0].name);
      this._closePopup();
    }
  }

  /** The headless/platform-less popup: search + a virtualized single-pick list of
   * {@link columnRenderer} rows, an empty row first when nullable (it answers null), the current
   * value highlighted. Click or Enter picks and closes, Esc closes without a write. */
  private _openNative(scope: Scope, seed?: string): void {
    const renderer = columnRenderer(this._table!);
    const current = this.value.peek();
    const search = new TextInput({placeholder: 'Search…', search: true, inline: true});
    const all = (): (string | null)[] => this.nullable ? [null, ...this.names()] : this.names();
    const filtered = (query: string): (string | null)[] => {
      const q = query.trim().toLowerCase();
      // the empty row always passes, as the Dart search's does (column_combo_box.dart:340)
      return q ? all().filter((name) => name === null || name.toLowerCase().includes(q)) : all();
    };
    const list = new VirtualList<string | null>({
      itemHeight: ROW_HEIGHT,
      render: (name, _index, row) => {
        row.classList.toggle('u2-columns-picked', name !== null && name === current);
        const body = name === null ? div([], 'u2-column-combo-none') : renderer.listItem(name);
        return div([body], 'u2-columns-option');
      },
    });
    list.root.classList.add('u2-columns-list');
    const pick = (index: number) => {
      const name = filtered(search.value.peek())[index];
      if (name === undefined)
        return;
      this._commit(name);
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
        pick(Number(row.dataset.index));
    };
    list.root.addEventListener('click', onClick);
    scope.own(() => list.root.removeEventListener('click', onClick));
    scope.effect(() => {
      list.setItems(filtered(search.value.value));
      // a re-filter drops the walked highlight (the rows under it moved), so Enter falls back
      // to the unique match instead of picking a dead index
      list.selectedIndex.value = -1;
    });
    const at = current === null ? -1 : filtered(search.value.peek()).indexOf(current);
    if (at >= 0) {
      list.selectedIndex.value = at;
      list.scrollToIndex(at);
    }

    const popup = divV([search, list], 'u2-columns-popup');
    popup.classList.add('u2-column-combo-popup');
    const onPopupKeyDown = (e: Event) => {
      const key = (e as KeyboardEvent).key;
      if (key === 'Enter') {
        e.preventDefault();
        const index = list.selectedIndex.peek();
        if (index >= 0) {
          pick(index);
          return;
        }
        const rows = filtered(search.value.peek());
        const matches = rows.filter((n) => n !== null);
        if (matches.length === 1)
          pick(rows.indexOf(matches[0]));
        return;
      }
      if (list.root.contains(e.target as HTMLElement))
        return;
      if (key === 'ArrowDown' || key === 'ArrowUp') {
        e.preventDefault();
        move(key === 'ArrowDown' ? 1 : -1);
      }
    };
    const onDismiss = () => this._closePopup();
    popup.addEventListener('keydown', onPopupKeyDown);
    popup.addEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
    scope.own(() => {
      popup.removeEventListener('keydown', onPopupKeyDown);
      popup.removeEventListener(OVERLAY_CLOSE_EVENT, onDismiss);
    });
    Overlay.show(this.box, popup, scope);
    if (seed !== undefined)
      search.value.value = seed;
    search.root.querySelector<HTMLElement>('input')?.focus();
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
