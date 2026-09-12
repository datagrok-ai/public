/* An in-memory `DataFrameLike` over records — what `MemoryDomainBackend.frame()` hands out, so a
   source over the memory backend is a frame host exactly like one over the platform: the same
   row projection, the same edit-state contract, the gallery and the headless tests structurally
   on the platform path. The records are the frame (`rows`); the edit state writes them and
   fires the events a host follows. */
import {Emitter} from '../core/emitter.js';
import type {ColumnLike, DataFrameLike} from './df-bindings.js';

export interface MemoryColumn {
  name: string;
  type: string;
  semType?: string | null;
}

export class MemoryFrame implements DataFrameLike {
  readonly rows: Record<string, unknown>[];
  readonly columns: DataFrameLike['columns'];
  readonly selection: unknown = null;
  readonly filter: unknown = null;
  readonly onCurrentRowChanged = new Emitter<unknown>();
  readonly onValuesChanged = new Emitter<unknown>();
  readonly onSelectionChanged = new Emitter<unknown>();
  readonly onFilterChanged = new Emitter<unknown>();
  readonly onColumnsChanged = new Emitter<unknown>();
  readonly onRowsAdded = new Emitter<unknown>();
  readonly onRowsRemoved = new Emitter<unknown>();

  private _current = -1;

  constructor(columns: MemoryColumn[], rows: Record<string, unknown>[]) {
    this.rows = rows;
    const list: ColumnLike[] = columns.map((c) => ({name: c.name, type: c.type, semType: c.semType ?? null,
      isNone: (row: number) => this.rows[row]?.[c.name] === null || this.rows[row]?.[c.name] === undefined}));
    this.columns = {
      names: () => list.map((c) => c.name),
      byName: (name) => list.find((c) => c.name === name) ?? null,
    };
  }

  get rowCount(): number {
    return this.rows.length;
  }

  get currentRowIdx(): number {
    return this._current;
  }

  set currentRowIdx(index: number) {
    if (index === this._current)
      return;
    this._current = index;
    this.onCurrentRowChanged.fire(index);
  }

  get(column: string, row: number): unknown {
    return this.rows[row]?.[column];
  }

  set(column: string, row: number, value: unknown): void {
    this.rows[row][column] = value;
    this.onValuesChanged.fire({column, row});
  }

  /** How many subscriptions are still held on the frame — what a disposal test asserts on. */
  get subscriberCount(): number {
    return [this.onCurrentRowChanged, this.onValuesChanged, this.onSelectionChanged, this.onFilterChanged,
      this.onColumnsChanged, this.onRowsAdded, this.onRowsRemoved].reduce((n, e) => n + e.count, 0);
  }
}
