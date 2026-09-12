/* A domain table as data (GOAL ruling 2): the one object every domain control binds to. The
   collection is the frame the backend hands out with its writer attached; `rows` is the
   key-addressed view lists and pickers read, `currentRow` what a form edits, and every edit goes
   through the table's `EditState` — drafts are rows, save is one transaction through the
   session. Platform-free: everything reaches the platform through `backends.domain`. Row
   selection arrives with the grid in phase 2. */
import {signal, computed, batch, Signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {Access} from '../core/access.js';
import {Filters} from '../core/filter/index.js';
import type {FilterGroup, FilterProperty, FilterSchema} from '../core/filter/index.js';
import {backends, requireBackend} from './backends.js';
import {subBind} from './sub-bind.js';
import {Rows} from './rows-like.js';
import type {RowsLike, RowView} from './rows-like.js';
import {FrameRows} from './df-rows.js';
import type {EditState} from './edit-state.js';
import type {DataFrameLike} from './df-bindings.js';
import {DomainBackendError} from './domain-backend.js';
import type {DomainBackend, DomainFrameLike, DomainQueryLike, DomainTableInfoLike, DomainTableLike}
  from './domain-backend.js';
import {SingleSession} from './session.js';
import type {DomainSession} from './session.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';
import type {ComponentEnv, ComponentStart} from '../spec/registry.js';

export type DomainSourceState = 'idle' | 'loading' | 'ready' | 'error';

/** The filter schema of the table plus what the controls read off it: every row property
 * (a ref column's `ref` is its target address) and the table's naming. */
export interface DomainSchema extends FilterSchema {
  info: DomainTableInfoLike;
}

export interface DomainSourceOptions {
  /** `'<schema>.<table>'`. */
  table?: string;
  /** A smart-filter string, or a filter tree built in code. */
  query?: string | FilterGroup;
  pageSize?: number;
  /** Ask for the per-row access columns with every row (default true); the table-level access is
   * always fetched. */
  withAccess?: boolean;
  /** What every draft starts with — a parent's id on a child table. */
  defaults?: Record<string, unknown>;
  /** A source that loads nothing and holds one pristine draft — what a create form binds to;
   * `save()` inserts it. */
  draft?: boolean;
  /** The unit of work `save`/`discard` go through; the source's own session of one by default. */
  session?: DomainSession;
}

type Row = Record<string, unknown>;

const NO_ENV: ComponentEnv = {designTime: false, subBinds: {}, resolve: () => null};
const NO_INFO: DomainTableInfoLike = {nameColumn: null, singularName: '', pluralName: '', businessKey: []};
const REF_ADDRESS = /^\w+\.\w+$/;

export class DomainSource extends Component implements BindSource, ComponentStart {
  readonly table: string;
  readonly pageSize: number;
  readonly withAccess: boolean;
  readonly defaults: Record<string, unknown>;
  readonly isDraft: boolean;
  /** A signal, so the panel edits it live and a bound path drives it. */
  readonly query: Signal<string | FilterGroup>;
  readonly df: ReadonlySignal<DataFrameLike | undefined>;
  readonly rows: RowsLike<RowView>;
  /** Two-way; mirrors the frame's current row. */
  readonly currentRow: Signal<RowView | null>;
  readonly state: ReadonlySignal<DomainSourceState>;
  /** What the last load or save threw — a `DomainBackendError` keeps its `code`. */
  readonly error: ReadonlySignal<unknown>;
  /** How many rows match the query, once counted. */
  readonly total: ReadonlySignal<number | null>;
  readonly access: ReadonlySignal<Access>;
  /** The table's single writer, once loaded — cell state and deletes for the controls. */
  readonly edit: ReadonlySignal<EditState | undefined>;
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  /** "50 of 2,077", "12 issues", "3 unsaved changes" — what a status bar shows. */
  readonly summary: ReadonlySignal<string>;
  /** The unit of work: `save`/`discard` go through it, a paired form follows its events. */
  readonly session: DomainSession;
  /** Bumped by a list on Enter: the form paired through this source takes the focus. */
  readonly activate = signal(0);

  private readonly _env: ComponentEnv;
  private readonly _backend: DomainBackend;
  private readonly _df = signal<DataFrameLike | undefined>(undefined);
  private readonly _edit = signal<EditState | undefined>(undefined);
  private readonly _state = signal<DomainSourceState>('idle');
  private readonly _error = signal<unknown>(undefined);
  private readonly _total = signal<number | null>(null);
  private readonly _access = signal(Access.readOnly);
  private readonly _errorStep: ReadonlySignal<string>;
  /** The `source` step: the resolver walks to a signal, so the source hands itself over in one. */
  private readonly _self = signal<DomainSource>(this);
  /** Bumped on every edit-state change, so everything reading a row live re-reads. */
  private readonly _version = signal(0);
  private readonly _rows: FrameRows;
  private readonly _columns = new Map<string, Signal<unknown>>();
  private readonly _guards = new Set<() => string | null>();
  private _table: DomainTableLike | undefined;
  private _ready: Promise<DomainTableLike> | undefined;
  private _frame: DomainFrameLike | undefined;
  private _schema: DomainSchema = {properties: [], info: NO_INFO};
  private _unwire: (() => void) | undefined;
  private _gen = 0;
  private _loaded = 0;
  private _done = false;

  constructor(options: DomainSourceOptions, env: ComponentEnv = NO_ENV) {
    super();
    this._env = env;
    this.table = options.table ?? '';
    this.pageSize = options.pageSize ?? 50;
    this.withAccess = options.withAccess ?? true;
    this.defaults = {...options.defaults};
    this.isDraft = options.draft ?? false;
    this.query = signal<string | FilterGroup>(options.query ?? '');
    this._backend = requireBackend(this, backends.domain, 'domain tables');
    this.df = this._df;
    this.error = this._error;
    this.access = this._access;
    this.edit = this._edit;
    this.state = this._state;
    this.total = this._total;
    this._errorStep = computed(() => DomainSource._message(this._error.value));
    this.currentRow = signal<RowView | null>(null);
    this._rows = new FrameRows(this._df, this.scope,
      {onWrite: (id, column, value) => this._edit.peek()?.setValue(id, column, value)});
    this.rows = this._rows;
    this.isDirty = computed(() => this._edit.value?.isDirty.value ?? false);
    this.changeCount = computed(() => this._edit.value?.changeCount.value ?? 0);
    this.validity = computed(() => this._edit.value?.validity.value ?? null);
    this.isSaving = computed(() => this._edit.value?.isSaving.value ?? false);
    this.summary = computed(() => this._summary());
    this.session = options.session ?? new SingleSession(this);

    // a row that left the collection — a discarded draft, a saved delete — hands the current row
    // on to whatever sits at its place now (a saved draft under its new key, the next row), or none
    let currentAt = -1;
    this.effect(() => {
      const row = this.currentRow.value;
      const items = this.rows.items.value;
      if (row === null)
        return;
      if (this.rows.byKey(row.id) !== undefined) {
        currentAt = items.findIndex((r) => r.id === row.id);
        return;
      }
      this.currentRow.value = items[Math.min(currentAt, items.length - 1)] ?? null;
    });
    // the frame's current row and `currentRow` are one thing, written in either direction
    this.effect(() => {
      const row = this.currentRow.value;
      const d = this._df.peek();
      if (d === undefined || row === null)
        return;
      const at = this._rows.indexOf(row.id);
      if (at >= 0 && d.currentRowIdx !== at)
        d.currentRowIdx = at;
    });
    this.own(() => {
      this._gen++;
      this._drop();
    });
    this.registerFunction({name: 'refresh', description: 'Load the table again from the first page', inputs: [],
      apply: () => this.refresh()});
    this.registerFunction({name: 'loadMore', description: 'Append the next page', inputs: [],
      apply: () => this.loadMore()});
    this.registerFunction({name: 'save', description: 'Write every pending change as one transaction', inputs: [],
      apply: () => this.save()});
    this.registerFunction({name: 'discard', description: 'Drop every pending change', inputs: [],
      apply: () => this.discard()});
    this.registerFunction({name: 'newRow', description: 'Add a pristine draft row and make it current',
      inputs: [{name: 'values', type: 'object', nullable: true}],
      apply: (params) => this.newRow(params?.values as Row | undefined, {pristine: true})});
  }

  /** Phase two: the query may be bound to an input the form declares after this source; the
   * first page is loaded once it is known. A new query is a new collection — unless edits are
   * pending, which a programmatic re-query never drops silently (STATE-CONTRACT H6; the gate
   * that asks is phase 2). */
  start(): void {
    const bound = subBind(this._env, 'query');
    if (bound !== null)
      this.effect(() => this.query.value = DomainSource._queryOf(bound.value));
    this.effect(() => {
      this.query.value;
      if (!this.isDirty.peek())
        void this.refresh();
    });
  }

  get schema(): DomainSchema {
    return this._schema;
  }

  /** A draft over {@link defaults} and `values`, made current — the row a create form binds to;
   * `pristine` keeps it from arming the dirty gate until its first edit. Needs the table: call
   * after the source is ready. */
  newRow(values: Row = {}, options?: {pristine?: boolean}): RowView {
    const edit = this._edit.peek();
    if (edit === undefined)
      throw new Error(`${this.table}: the table is not loaded yet`);
    const key = edit.newRow({...this.defaults, ...values}, options);
    const row = this.rows.byKey(key);
    if (row === undefined)
      throw new Error(`${this.table}: the draft "${key}" is not among the rows`);
    this.currentRow.value = row;
    return row;
  }

  /** Reloads from the first page; pending changes are dropped. Stale answers — a load a later
   * refresh outran — touch nothing. A draft source loads no rows and starts on a pristine draft. */
  async refresh(): Promise<void> {
    const gen = ++this._gen;
    batch(() => {
      this._error.value = undefined;
      this._state.value = 'loading';
    });
    try {
      const table = await (this._ready ??= this._backend.table(this.table));
      if (gen !== this._gen)
        return;
      this._adopt(table);
      const [access, frame, total] = await Promise.all([
        table.access(), table.frame(this._spec(0)), this.isDraft ? 0 : table.count(this._filter())]);
      if (gen !== this._gen) {
        frame.dispose();
        return;
      }
      this._drop();
      this._frame = frame;
      this._loaded = frame.df.rowCount;
      this._done = this.isDraft || this._loaded < this.pageSize;
      this._wire(frame);
      batch(() => {
        this._access.value = Access.from(access);
        this._df.value = frame.df;
        this._edit.value = frame.edit;
        this._total.value = total;
        this._state.value = 'ready';
      });
      this._syncCurrent();
      if (this.isDraft)
        this.newRow({}, {pristine: true});
    } catch (e) {
      this._ready = undefined;
      if (gen !== this._gen)
        return;
      batch(() => {
        this._error.value = e;
        this._state.value = 'error';
      });
    }
  }

  /** Appends the next page into the same frame; pending changes stay. */
  async loadMore(): Promise<void> {
    const frame = this._frame;
    if (frame === undefined)
      return this.refresh();
    if (this._done || this._state.peek() === 'loading')
      return;
    const gen = this._gen;
    this._state.value = 'loading';
    try {
      const added = await frame.append(this._spec(this._loaded));
      if (gen !== this._gen)
        return;
      this._loaded += added;
      this._done = added < this.pageSize;
      this._state.value = 'ready';
    } catch (e) {
      if (gen !== this._gen)
        return;
      batch(() => {
        this._error.value = e;
        this._state.value = 'error';
      });
    }
  }

  /** Through the session — the one Save every button and shortcut runs. */
  save(): Promise<boolean> {
    return this.session.save();
  }

  discard(): void {
    this.session.discard();
  }

  /** A check Save runs first — a form registers its own, so a refusal names the field ("Cannot
   * save: Title is required"). Answers the unregister. */
  guard(check: () => string | null): () => void {
    this._guards.add(check);
    return () => this._guards.delete(check);
  }

  /** This source's pending changes as one transaction — what its session calls; a refusal is
   * the `error` (the summary and a paired form show it) and answers false. */
  async commit(): Promise<boolean> {
    const edit = this._edit.peek();
    if (edit === undefined || !edit.isDirty.peek())
      return false;
    this._error.value = undefined;
    for (const check of this._guards) {
      const problem = check();
      if (problem !== null)
        return this._refuse(problem);
    }
    const problem = edit.validity.peek();
    if (problem !== null)
      return this._refuse(problem);
    try {
      const saved = await edit.save();
      if (saved && !this.isDraft)
        this._total.value = await this._table!.count(this._filter());
      return saved;
    } catch (e) {
      this._error.value = e;
      return false;
    }
  }

  /** Drops this source's pending changes — what its session calls. */
  revert(): void {
    this._edit.peek()?.discard();
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    switch (name) {
      case '': case 'rows': return this.rows.items as unknown as Signal<unknown>;
      case 'currentRow': return {
        bindStep: (column) => column === '' ? null : this._columnSignal(column),
        bindProps: () => {
          const view = this._view(this.currentRow.peek());
          return this._schema.properties.map((p): BindProp => ({...p, name: p.name,
            writable: view.field(p.name) === 'editable'}));
        },
      };
      case 'source': return this._self as unknown as Signal<unknown>;
      case 'total': return this.total as unknown as Signal<unknown>;
      case 'state': return this.state as unknown as Signal<unknown>;
      case 'error': return this._errorStep as unknown as Signal<unknown>;
      case 'access': return this._access as unknown as Signal<unknown>;
      case 'isDirty': return this.isDirty as unknown as Signal<unknown>;
      default: return super.bindStep(name);
    }
  }

  bindProps(): BindProp[] {
    return [
      {name: 'rows', type: 'object', description: 'The rows loaded so far, drafts included', default: true},
      {name: 'currentRow', type: 'object', walkable: true, description: 'The row a form edits; its columns are steps'},
      {name: 'total', type: 'int', description: 'How many rows match the query'},
      {name: 'state', type: 'string', description: 'idle, loading, ready or error'},
      {name: 'error', type: 'string', description: 'Why the last load or save failed'},
      {name: 'access', type: 'object', description: 'What the caller may do with the table and its fields'},
      {name: 'isDirty', type: 'bool', description: 'Whether there are unsaved changes'},
      {name: 'source', type: 'object', description: 'The source itself — what a domain form, list or picker binds to'},
    ];
  }

  private _refuse(problem: string): false {
    this._error.value = new DomainBackendError('validation', `Cannot save: ${problem}`);
    return false;
  }

  private _adopt(table: DomainTableLike): void {
    if (this._table === table)
      return;
    this._table = table;
    this._schema = {
      properties: table.properties.map((p): FilterProperty => {
        const prop: FilterProperty = {...p, name: p.name!};
        if (REF_ADDRESS.test(p.semType ?? ''))
          prop.ref = p.semType;
        return prop;
      }),
      info: table.info,
    };
  }

  /** Lets go of the current frame and its writer — before a replacement, and at disposal. */
  private _drop(): void {
    this._unwire?.();
    this._unwire = undefined;
    const frame = this._frame;
    this._frame = undefined;
    frame?.dispose();
    batch(() => {
      this._df.value = undefined;
      this._edit.value = undefined;
      this._total.value = null;
      this.currentRow.value = null;
    });
  }

  /** Follows the writer's changes and the frame's current row. */
  private _wire(frame: DomainFrameLike): void {
    const bump = () => this._version.value = this._version.peek() + 1;
    const subs = [
      // the editor removes saved deletes and re-keys saved drafts without a frame event (H10): the
      // row projection follows its onChanged instead
      frame.edit.onChanged.subscribe(() => {
        this._rows.rebuild();
        bump();
      }),
      frame.df.onCurrentRowChanged.subscribe(() => this._syncCurrent()),
      frame.df.onValuesChanged.subscribe(bump),
    ];
    this._unwire = () => {
      for (const s of subs)
        s.unsubscribe();
    };
  }

  private _syncCurrent(): void {
    const d = this._df.peek();
    const key = d === undefined ? undefined : this._rows.keyAt(d.currentRowIdx);
    this.currentRow.value = key === undefined ? null : this._rows.byKey(key) ?? null;
  }

  /** The access a row is written under — its own (`Access.row`: a draft under `insert`, an existing
   * row as its `~can_*` columns say), no row the table's. */
  private _view(row: RowView | null): Access {
    const access = this._access.peek();
    return row === null ? access : access.row(row);
  }

  private _spec(offset: number): DomainQueryLike {
    return {filter: this._filter(), limit: this.isDraft ? 0 : this.pageSize, offset, withAccess: this.withAccess};
  }

  private _filter(): DomainQueryLike['filter'] {
    const q = this.query.peek();
    return typeof q === 'string' ? (q === '' ? undefined : q) : Filters.toDomainTree(q);
  }

  /** A two-way `currentRow.<column>` step: reads the current row live, writes through it. */
  private _columnSignal(column: string): Signal<unknown> {
    let sig = this._columns.get(column);
    if (sig === undefined) {
      const read = () => this.currentRow.peek()?.[column];
      const created = signal(read());
      this.effect(() => {
        this.rows.items.value;
        this._version.value;
        this.currentRow.value;
        created.value = read();
      });
      this.effect(() => {
        const value = created.value;
        const row = this.currentRow.peek();
        if (row !== null && !Component.sameValue(row[column], value))
          row[column] = value;
      });
      this._columns.set(column, sig = created);
    }
    return sig;
  }

  private _summary(): string {
    const error = this._error.value;
    if (error !== undefined)
      return DomainSource._message(error);
    const items = this.rows.items.value;
    const changes = this.changeCount.value;
    const plural = (n: number, one: string, many: string) => `${n.toLocaleString()} ${n === 1 ? one : many}`;
    if (changes > 0) {
      const deleted = items.filter((r) => r[Rows.STATE] === 'deleted').length;
      return deleted === changes ? `${plural(deleted, 'deletion', 'deletions')} pending` :
        plural(changes, 'unsaved change', 'unsaved changes');
    }
    const state = this.state.value;
    const singular = this._schema.info.singularName.toLowerCase() || 'row';
    // a draft is a row, not a row of the table yet; a draft source is a create form until the
    // draft is saved and a row of the table is what it holds
    const loaded = items.filter((r) => !Rows.isDraft(r)).length;
    if (this.isDraft && (loaded === 0 || (this.currentRow.value !== null && Rows.isDraft(this.currentRow.value))))
      return state === 'ready' ? `New ${singular}` : state === 'loading' ? 'Loading…' : '';
    if (state !== 'ready' && loaded === 0)
      return state === 'loading' ? 'Loading…' : '';
    const total = this.total.value;
    return !this.isDraft && total !== null && loaded < total ?
      `${loaded.toLocaleString()} of ${total.toLocaleString()}` :
      plural(loaded, singular, this._schema.info.pluralName.toLowerCase() || 'rows');
  }

  private static _message(error: unknown): string {
    return error === undefined ? '' : error instanceof Error ? error.message : String(error);
  }

  private static _queryOf(value: unknown): string | FilterGroup {
    return typeof value === 'string' ? value : value === null || value === undefined ? '' : value as FilterGroup;
  }
}
