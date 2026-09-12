/* The editing-state protocol a `DomainSource` drives (STATE-CONTRACT): what is dirty, what is
   valid, one write path, drafts, deletes, discard and save. The platform implementation is the
   js-api `DomainFrameEditor` (dg); `MemoryEditState` is the same machine over a `MemoryFrame` —
   the memory backend's writer, which the tests and the gallery run. */
import {signal, computed, ReadonlySignal} from '../core/signals.js';
import {Emitter} from '../core/emitter.js';
import type {ObservableLike} from '../core/widget-like.js';
import type {IProperty} from '../core/property-like.js';
import type {Access} from '../core/access.js';
import type {DomainTableLike, DomainTransactionOpLike} from './domain-backend.js';
import type {MemoryFrame} from './memory-frame.js';
import {FrameRows} from './df-rows.js';
import {Rows} from './rows-like.js';

export type RowState = '' | 'new' | 'modified' | 'deleted';

export interface EditState {
  readonly isDirty: ReadonlySignal<boolean>;
  /** Changed cells, as the platform editor counts: an edited row its cells, a new or deleted
   * row once, a pristine draft not at all. */
  readonly changeCount: ReadonlySignal<number>;
  /** The first blocking problem, null when the batch may be saved. */
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  /** The key of the row that changed, null for a change touching several. */
  readonly onChanged: ObservableLike<string | null>;
  isChanged(key: string, column: string): boolean;
  errorOf(key: string, column: string): string | null;
  setValue(key: string, column: string, value: unknown): void;
  /** Adds a draft; answers the key `RowsLike.byKey` finds it under. `pristine` is the insert-form
   * contract: the draft is saved with the batch but arms nothing until its first write. */
  newRow(values?: Record<string, unknown>, options?: {pristine?: boolean}): string;
  markDeleted(key: string): void;
  /** Undoes {@link markDeleted}: the row is back to what it was before — edited or clean. */
  unmarkDeleted(key: string): void;
  discard(): void;
  /** Resolves to whether the batch landed; a backend refusal is thrown as it came. */
  save(): Promise<boolean>;
  /** Detaches from whatever it tracks — the frame, the platform editor. */
  dispose(): void;
}

type Row = Record<string, unknown>;

const EMPTY = 'Value can\'t be empty';
const isEmpty = (v: unknown) => v === null || v === undefined || v === '';

/** Over the records of a `MemoryFrame`, keyed the way `FrameRows` keys them (the id cell, the
 * draft key for a row without one): originals per (row, column), the row's state in its `~state`
 * cell, the batch built in row order as the js-api editor builds it, every write announced
 * through `onChanged` — a save removes deleted rows and re-keys drafts without a frame event, as
 * the editor does (H10). No conflict UI — a refused transaction is thrown to the caller with the
 * batch still pending. */
export class MemoryEditState implements EditState {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  readonly onChanged = new Emitter<string | null>();

  private readonly _originals = new Map<Row, Map<string, unknown>>();
  private readonly _pristine = new WeakSet<Row>();
  private readonly _version = signal(0);
  private readonly _saving = signal(false);

  /** `access` is what a row is written under (`Access.row`): a draft under `insert`, an existing
   * row under `edit` as its own columns narrow it — the columns a save may send. */
  constructor(private readonly _table: DomainTableLike, readonly df: MemoryFrame, private readonly _access: Access) {
    this.changeCount = computed(() => {
      this._version.value;
      return this.df.rows.reduce((n, row) => n + this._contribution(row), 0);
    });
    this.isDirty = computed(() => this.changeCount.value > 0);
    this.validity = computed(() => {
      this._version.value;
      for (const row of this.df.rows) {
        for (const prop of this._table.properties) {
          const error = this._error(row, prop);
          if (error !== null)
            return error;
        }
      }
      return null;
    });
    this.isSaving = this._saving;
  }

  /** The frame row behind a key, -1 when the frame does not hold it. */
  indexOf(key: string): number {
    const df = this.df;
    return Rows.draftRow(key, df.rowCount) ?? df.rows.findIndex((_, i) => FrameRows.keyOf(df, i) === key);
  }

  isChanged(key: string, column: string): boolean {
    const row = this._row(key);
    return row !== undefined && (row[Rows.STATE] === 'new' || this._originals.get(row)?.has(column) === true);
  }

  errorOf(key: string, column: string): string | null {
    const prop = this._table.properties.find((p) => p.name === column);
    const row = this._row(key);
    return prop === undefined || row === undefined ? null : this._error(row, prop);
  }

  setValue(key: string, column: string, value: unknown): void {
    const row = this._row(key);
    if (row === undefined || Object.is(row[column], value))
      return;
    const state = row[Rows.STATE] as RowState | undefined;
    if (state !== 'new') {
      const originals = this._originals.get(row) ?? new Map<string, unknown>();
      this._originals.set(row, originals);
      if (!originals.has(column))
        originals.set(column, row[column]);
      else if (Object.is(originals.get(column), value))
        originals.delete(column);
      if (state !== 'deleted')
        row[Rows.STATE] = originals.size > 0 ? 'modified' : '';
    }
    this._pristine.delete(row);
    row[column] = value;
    this._touch(key);
  }

  newRow(values: Row = {}, options: {pristine?: boolean} = {}): string {
    const row: Row = {...values, id: null, [Rows.STATE]: 'new'};
    if (options.pristine)
      this._pristine.add(row);
    this.df.rows.push(row);
    this.df.onRowsAdded.fire(undefined);
    const key = Rows.draftKey(this.df.rowCount - 1);
    this._touch(key);
    return key;
  }

  markDeleted(key: string): void {
    const row = this._row(key);
    if (row === undefined)
      return;
    if (row[Rows.STATE] === 'new')
      this._remove((r) => r === row);
    else
      row[Rows.STATE] = 'deleted';
    this._touch(key);
  }

  unmarkDeleted(key: string): void {
    const row = this._row(key);
    if (row === undefined || row[Rows.STATE] !== 'deleted')
      return;
    row[Rows.STATE] = (this._originals.get(row)?.size ?? 0) > 0 ? 'modified' : '';
    this._touch(key);
  }

  discard(): void {
    for (const [row, originals] of this._originals) {
      for (const [column, original] of originals)
        row[column] = original;
    }
    this._remove((row) => row[Rows.STATE] === 'new');
    this._settle();
  }

  async save(): Promise<boolean> {
    if (this.validity.peek() !== null || this._saving.peek())
      return false;
    this._saving.value = true;
    try {
      const pending = this._ops();
      const results = pending.length === 0 ? [] : await this._table.transaction(pending.map((p) => p.op));
      const removed = new Set<Row>();
      for (const [i, {op, row}] of pending.entries()) {
        if (op.op === 'delete')
          removed.add(row);
        else
          Object.assign(row, results[i]);
      }
      // as the platform editor's _applyResults: rows leave and drafts are re-keyed with no frame
      // event — the host rebuilds on onChanged (H10)
      this.df.rows.splice(0, this.df.rowCount, ...this.df.rows.filter((row) => !removed.has(row)));
      this._settle();
      return true;
    } finally {
      this._saving.value = false;
    }
  }

  dispose(): void {
    this.onChanged.clear();
  }

  /** The schema's own rules on one value: required, choices, min and max — what the backend
   * refuses a transaction on, and what the state reports per cell before it is sent. */
  static problemOf(prop: IProperty, value: unknown): string | null {
    if (isEmpty(value))
      return prop.nullable === false ? EMPTY : null;
    const choices = prop.choices;
    if (choices !== undefined && choices !== null && choices.length > 0 && !choices.includes(String(value)))
      return `Must be one of: ${choices.join(', ')}`;
    if (typeof value === 'number') {
      if (prop.min !== undefined && prop.min !== null && value < prop.min)
        return `Must be at least ${prop.min}`;
      if (prop.max !== undefined && prop.max !== null && value > prop.max)
        return `Must be at most ${prop.max}`;
    }
    return null;
  }

  private _ops(): {op: DomainTransactionOpLike, row: Row}[] {
    const table = this._table.address.split('.').pop()!;
    const pending: {op: DomainTransactionOpLike, row: Row}[] = [];
    for (const row of this.df.rows) {
      const state = row[Rows.STATE] as RowState | undefined;
      const values: Row = {};
      const writable = this._writable(row);
      if (state === 'deleted')
        pending.push({row, op: {op: 'delete', table, id: String(row.id)}});
      else if (state === 'new') {
        for (const name of writable) {
          if (row[name] != null)
            values[name] = row[name];
        }
        pending.push({row, op: {op: 'insert', table, values}});
      } else if (state === 'modified') {
        for (const name of this._originals.get(row)?.keys() ?? []) {
          if (writable.includes(name))
            values[name] = row[name];
        }
        if (Object.keys(values).length > 0) {
          const version = row.version;
          pending.push({row, op: {op: 'update', table, id: String(row.id), values,
            ...(typeof version === 'number' ? {expectedVersion: version} : {})}});
        }
      }
    }
    return pending;
  }

  private _writable(row: Row): string[] {
    const access = this._access.row(row);
    return this._table.properties.map((p) => p.name!).filter((name) => access.field(name) === 'editable');
  }

  private _remove(drop: (row: Row) => boolean): void {
    const kept = this.df.rows.filter((row) => !drop(row));
    if (kept.length === this.df.rowCount)
      return;
    this.df.rows.splice(0, this.df.rowCount, ...kept);
    this.df.onRowsRemoved.fire(undefined);
  }

  /** Every row clean. */
  private _settle(): void {
    for (const row of this.df.rows) {
      row[Rows.STATE] = '';
      this._pristine.delete(row);
    }
    this._originals.clear();
    this._touch(null);
  }

  private _row(key: string): Row | undefined {
    const at = this.indexOf(key);
    return at < 0 ? undefined : this.df.rows[at];
  }

  private _contribution(row: Row): number {
    switch (row[Rows.STATE]) {
      case 'new': return this._pristine.has(row) ? 0 : 1;
      case 'deleted': return 1;
      case 'modified': return this._originals.get(row)?.size ?? 0;
      default: return 0;
    }
  }

  /** A rule broken on a draft, or on a cell of an existing row the batch will send. */
  private _error(row: Row, prop: IProperty): string | null {
    const state = row[Rows.STATE];
    if (!state || state === 'deleted')
      return null;
    const changed = state === 'new' || this._originals.get(row)?.has(prop.name!) === true;
    return changed ? MemoryEditState.problemOf(prop, row[prop.name!]) : null;
  }

  private _touch(key: string | null): void {
    this._version.value = this._version.peek() + 1;
    this.onChanged.fire(key);
  }
}
