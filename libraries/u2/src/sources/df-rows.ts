/* `RowsLike` over a DataFrame (GOAL step-back: the frame is the collection, rows a derived view):
   one proxy per row key whose reads go to the frame live and whose writes go to the frame — or to
   the edit state a source routes them through — and an item list rebuilt as the frame's values,
   filter or rows change. A row the editor marked deleted stays an item (its `~state` says so), as
   the platform grid keeps it until the save — a list shows it struck through with Restore. */
import {signal, Signal} from '../core/signals.js';
import type {ReadonlySignal} from '../core/signals.js';
import type {Scope} from '../core/scope.js';
import type {DataFrameLike} from './df-bindings.js';
import {Rows} from './rows-like.js';
import type {RowsLike, RowView} from './rows-like.js';

/** What a row proxy reads and writes through. `keys` is what `Object.keys(row)` enumerates. */
export interface RowReader {
  get(column: string): unknown;
  set(column: string, value: unknown): void;
  keys(): string[];
}

export interface FrameRowsOptions {
  /** The column holding the row key (default `id`); a row with none is keyed by its index. */
  idColumn?: string;
  /** Takes every write made through a row instead of `df.set` — how a source hands its edits to
   * the state that tracks them. */
  onWrite?: (id: string, column: string, value: unknown) => void;
}

export class FrameRows implements RowsLike<RowView> {
  readonly items: ReadonlySignal<readonly RowView[]>;

  private readonly _items: Signal<readonly RowView[]>;
  private readonly _proxies = new Map<string, RowView>();
  private readonly _index = new Map<string, number>();
  private readonly _visible = new Set<string>();
  private readonly _idColumn: string;
  private readonly _onWrite: FrameRowsOptions['onWrite'];

  constructor(private readonly _df: ReadonlySignal<DataFrameLike | undefined>, scope: Scope,
    options: FrameRowsOptions = {}) {
    this._idColumn = options.idColumn ?? 'id';
    this._onWrite = options.onWrite;
    this._items = signal<readonly RowView[]>([]);
    this.items = this._items;
    let subs: {unsubscribe(): void}[] = [];
    const drop = () => {
      for (const s of subs)
        s.unsubscribe();
      subs = [];
    };
    scope.effect(() => {
      const d = _df.value;
      drop();
      if (d !== undefined) {
        const rebuild = () => this._rebuild();
        for (const event of [d.onValuesChanged, d.onFilterChanged, d.onRowsAdded, d.onRowsRemoved,
          d.onColumnsChanged])
          subs.push(event.subscribe(rebuild));
      }
      this._rebuild();
    });
    scope.own(drop);
  }

  keyOf(row: RowView): string {
    return row.id;
  }

  /** An item's row — a filtered-out row is not addressable. */
  byKey(key: string): RowView | undefined {
    return this._visible.has(key) ? this._proxy(key) : undefined;
  }

  /** Re-reads the frame: what a source calls when the writer changed rows without a frame event —
   * the editor removes saved deletes and re-keys saved drafts with `notify: false` (H10). */
  rebuild(): void {
    this._rebuild();
  }

  /** The row index behind a key, -1 when the frame no longer holds it. */
  indexOf(key: string): number {
    return this._index.get(key) ?? -1;
  }

  /** The key of the row at `index`, as {@link items} reports it. */
  keyAt(index: number): string | undefined {
    const d = this._df.peek();
    return d === undefined || index < 0 || index >= d.rowCount ? undefined : FrameRows.keyOf(d, index, this._idColumn);
  }

  /** The key of the row at `index` of `df`: its id cell, the draft key while it has none. */
  static keyOf(df: DataFrameLike, index: number, idColumn = 'id'): string {
    const id = df.columns.byName(idColumn) === null ? null : df.get(idColumn, index);
    return id === null || id === undefined || id === '' ? Rows.draftKey(index) : String(id);
  }

  /** A live row over `reader`: `id` is the key, every other name reads and writes through the
   * reader, `Object.keys` enumerates the reader's columns. */
  static proxy(id: () => string, reader: RowReader): RowView {
    return new Proxy({} as RowView, {
      get: (_, prop) => prop === 'id' ? id() : typeof prop === 'string' ? reader.get(prop) : undefined,
      set: (_, prop, value) => {
        if (typeof prop === 'string' && prop !== 'id')
          reader.set(prop, value);
        return true;
      },
      has: (_, prop) => prop === 'id' || (typeof prop === 'string' && reader.keys().includes(prop)),
      ownKeys: () => ['id', ...reader.keys().filter((k) => k !== 'id')],
      getOwnPropertyDescriptor: (_, prop) => typeof prop !== 'string' ? undefined :
        {configurable: true, enumerable: true, writable: prop !== 'id',
          value: prop === 'id' ? id() : reader.get(prop)},
    });
  }

  private _proxy(key: string): RowView {
    let proxy = this._proxies.get(key);
    if (proxy === undefined) {
      const reader: RowReader = {
        get: (column) => this._read(key, column),
        set: (column, value) => {
          if (this._onWrite)
            this._onWrite(key, column, value);
          else {
            const at = this.indexOf(key);
            if (at >= 0)
              this._df.peek()!.set(column, at, value);
          }
        },
        // the service columns are read by name, never enumerated, so a spread of a row carries no
        // editing state into a payload (H7)
        keys: () => (this._df.peek()?.columns.names() ?? []).filter((name) => !Rows.isService(name)),
      };
      proxy = FrameRows.proxy(() => key, reader);
      this._proxies.set(key, proxy);
    }
    return proxy;
  }

  /** A cell with no value is null — never the column type's sentinel `get` answers for it. */
  private _read(key: string, column: string): unknown {
    const d = this._df.peek();
    const at = this.indexOf(key);
    const col = d === undefined || at < 0 ? null : d.columns.byName(column);
    return col === null ? undefined : col.isNone?.(at) ? null : d!.get(column, at);
  }

  private _rebuild(): void {
    const d = this._df.peek();
    this._index.clear();
    this._visible.clear();
    const items: RowView[] = [];
    if (d !== undefined) {
      const filter = d.filter as {get?(i: number): boolean} | null | undefined;
      const state = d.columns.byName(Rows.STATE) === null ? null : Rows.STATE;
      for (let i = 0; i < d.rowCount; i++) {
        const key = this.keyAt(i)!;
        this._index.set(key, i);
        // the editor masks the rows it marked deleted out of the frame's filter (cooperative
        // filtering); they stay items here — struck through, restorable — until the save
        if (typeof filter?.get === 'function' && !filter.get(i) &&
            !(state !== null && d.get(state, i) === 'deleted'))
          continue;
        this._visible.add(key);
        items.push(this._proxy(key));
      }
    }
    for (const key of [...this._proxies.keys()]) {
      if (!this._index.has(key))
        this._proxies.delete(key);
    }
    this._items.value = items;
  }
}

/** Key-addressed live rows over a (possibly repointing) frame signal; subscriptions live on `scope`. */
export function frameRows(df: ReadonlySignal<DataFrameLike | undefined>, scope: Scope,
  options?: FrameRowsOptions): FrameRows {
  return new FrameRows(df, scope, options);
}
