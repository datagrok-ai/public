/* A domain table as data (GOAL ruling 2): the one object every domain control binds to. The
   collection is the frame the backend hands out with its writer attached; `rows` is the
   key-addressed view lists and pickers read, `currentRow` what a form edits, and every edit goes
   through the table's `EditState` — drafts are rows, save is one transaction through the
   session every source of a spec or an app shares. Platform-free: everything reaches the
   platform through `backends.domain`. */
import {signal, computed, batch, Signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {Access} from '../core/access.js';
import {plural} from '../core/text.js';
import {Filters} from '../core/filter/index.js';
import type {FilterGroup, FilterNode, FilterProperty, FilterScalar, FilterSchema} from '../core/filter/index.js';
import {backends, requireBackend} from './backends.js';
import {subBind} from './sub-bind.js';
import {Rows} from './rows-like.js';
import type {DomainRowLike, RowsLike, RowValues, RowView} from './rows-like.js';
import {FrameRows} from './df-rows.js';
import type {EditState} from './edit-state.js';
import type {DataFrameLike} from './df-bindings.js';
import {DomainBackendError} from './domain-backend.js';
import type {DomainBackend, DomainDeletedMode, DomainFrameLike, DomainQueryLike, DomainReadScope,
  DomainTableInfoLike, DomainTableLike} from './domain-backend.js';
import {SharedSession} from './session.js';
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
  /** A case-insensitive substring over the table's searchable columns, ANDed with the query. */
  search?: string;
  /** `'col,!col'` — the row order, `!` for descending; the table's own default when empty. The
   * system columns are orderable (`!updated_on` is what a trash list reads newest-first by),
   * even though they cannot be projected. */
  sort?: string;
  /** Which rows the source answers (default `'exclude'` — the live ones). Anything else is a
   * trash source: the rows carry `~is_deleted` and the access is narrowed to no edit and no
   * insert until they are restored. The backend must declare `restore`. */
  deleted?: DomainDeletedMode;
  /** Whether the collection follows the server: the source probes the table every
   * {@link DomainSourceOptions.liveMs} and refreshes while it is clean, marking itself
   * {@link DomainSource.stale} while it is not. A signal on the source, so a menu toggles it;
   * the timer rides the source's scope, so disposing the source stops it. */
  live?: boolean;
  /** How often a `live` source probes the server (default 30 s). */
  liveMs?: number;
  pageSize?: number;
  /** Ask for the per-row access columns with every row (default true); the table-level access is
   * always fetched. */
  withAccess?: boolean;
  /** Ref columns whose target names ride with the rows (`~caption_<col>`), so a list, a grid and a
   * form show names without a lookup per id. By default every ref column into another domain table
   * that this caller may see; `[]` turns it off. Over a backend that does not project them
   * (`support.captions` false) none are asked for, and a form resolves each per row instead. */
  captions?: string[];
  /** What every draft starts with — a parent's id on a child table. */
  defaults?: Record<string, unknown>;
  /** A source that loads no rows: the frame exists, so drafts can be added — a child collection
   * under a draft parent. */
  empty?: boolean;
  /** `empty` plus one pristine draft — what a create form binds to; `save()` inserts it. */
  draft?: boolean;
  /** The unit of work `save`/`discard` go through: the ambient session of the spec or app being
   * built, else a session of its own. */
  session?: DomainSession;
}

const NO_ENV: ComponentEnv = {designTime: false, subBinds: {}, resolve: () => null};
const NO_INFO: DomainTableInfoLike = {nameColumn: null, singularName: '', pluralName: '', businessKey: [],
  rowAddress: 'businessKey', searchableColumns: [], constraints: [], refFilters: {}, permissions: [], childTables: []};
const REF_ADDRESS = /^\w+\.\w+$/;
/** A draft id as a VALUE of the string query — quoted, which is the only form the grammar takes
 * one in; `~new:` anywhere else is ordinary text. */
const DRAFT_LITERAL = new RegExp(`(['"])${Rows.DRAFT_PREFIX}[^'"]*\\1`);

export class DomainSource<TRow extends DomainRowLike = DomainRowLike> extends Component
  implements BindSource, ComponentStart {
  /** Every source alive — what resolves a draft id to its row across sources ({@link draftOf}). */
  static readonly live = new Set<DomainSource>();
  /** The `code` of a refusal the backend already reported ({@link refuse}). */
  static readonly REFUSED = 'refused';
  /** How often a `live` source probes the server by default. */
  static readonly liveMs = 30000;
  /** How many probes may fail in a row before the poll gives up ({@link DomainSource.live}). */
  static readonly liveFailures = 3;

  readonly table: string;
  readonly pageSize: number;
  readonly withAccess: boolean;
  readonly defaults: Record<string, unknown>;
  readonly isDraft: boolean;
  /** Loads no rows (`draft` implies it). */
  readonly isEmpty: boolean;
  /** A signal, so the panel edits it live and a bound path drives it. */
  readonly query: Signal<string | FilterGroup>;
  /** The search text, re-queried on change like {@link query}. */
  readonly search: Signal<string>;
  /** The row order — see {@link DomainSourceOptions.sort}; a signal, so a mode that wants its
   * own order (the trash, newest first) is a write and the collection re-reads. */
  readonly sort: Signal<string>;
  /** Which rows the collection holds — a signal, so an app's trash mode is a flip of it and the
   * search box, the filters and the list stay bound to the one source. */
  readonly deleted: Signal<DomainDeletedMode>;
  /** Whether the collection follows the server — see {@link DomainSourceOptions.live}. */
  readonly live: Signal<boolean>;
  /** The rows the source holds are no longer the server's: a live probe saw the collection move
   * while this session had changes of its own to protect. Cleared by the next load. */
  readonly stale: ReadonlySignal<boolean>;
  /** A source over deleted rows, or over a table `transaction` does not land on
   * (`support.transaction` false — the whole write path): nothing in it may be edited, inserted
   * or deleted, and `save` refuses. The one place every control's write affordance follows. */
  readonly readOnly: ReadonlySignal<boolean>;
  readonly df: ReadonlySignal<DataFrameLike | undefined>;
  readonly rows: RowsLike<RowView<TRow>>;
  /** Two-way; mirrors the frame's current row. */
  readonly currentRow: Signal<RowView<TRow> | null>;
  /** The frame's selected rows (a grid's); empty over a frame without a selection. */
  readonly selection: ReadonlySignal<readonly RowView<TRow>[]>;
  readonly state: ReadonlySignal<DomainSourceState>;
  /** What the last load or save threw — a `DomainBackendError` keeps its `code`. */
  readonly error: ReadonlySignal<unknown>;
  /** The row the current refusal is about, when one names a row — what a list marks. */
  readonly problemRow: ReadonlySignal<string | null>;
  /** How a refusal the writer raised over one cell is worded; the table's registry sets it (the
   * row's caption and the column's), and without it the writer's own message stands. */
  nameCell: ((row: RowView, column: string, problem: string) => string) | undefined;
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
  private readonly _problemRow = signal<string | null>(null);
  private readonly _total = signal<number | null>(null);
  private readonly _access = signal(Access.readOnly);
  private readonly _selection = signal<readonly RowView<TRow>[]>([]);
  private readonly _isStale = signal(false);
  /** Whether the table's backend writes at all (`support.transaction`), once the table is known. */
  private readonly _writable = signal(true);
  private readonly _errorStep: ReadonlySignal<string>;
  /** The `source` step: the resolver walks to a signal, so the source hands itself over in one. */
  private readonly _self = signal<DomainSource<TRow>>(this);
  /** Bumped on every edit-state change, so everything reading a row live re-reads. */
  private readonly _version = signal(0);
  private readonly _rows: FrameRows;
  private readonly _columns = new Map<string, Signal<unknown>>();
  private readonly _guards = new Set<() => string | null>();
  /** What the caller listed as {@link DomainSourceOptions.captions}; undefined is the default. */
  private readonly _captionOption: string[] | undefined;
  private _table: DomainTableLike | undefined;
  private _ready: Promise<DomainTableLike> | undefined;
  private _frame: DomainFrameLike | undefined;
  private _schema: DomainSchema = {properties: [], info: NO_INFO};
  private _unwire: (() => void) | undefined;
  private _gen = 0;
  private _loaded = 0;
  private _done = false;
  /** Set by {@link markStale}: the loaded window is no longer the server's, so the next load
   * starts over instead of asking for the page after it. */
  private _restart = false;
  /** Set while {@link rebind} rewrites the query: the caller re-reads, and a refresh here would
   * drop the frame the batch just landed in. */
  private _rebinding = false;
  private readonly _liveMs: number;
  private _timer: ReturnType<typeof setInterval> | undefined;
  /** What the last probe answered — the pair the next one is compared against; undefined until
   * the first probe after a load, which only sets the baseline. */
  private _probed: {count: number, last: string | null} | undefined;
  private _probing = false;
  private _failures = 0;

  constructor(options: DomainSourceOptions, env: ComponentEnv = NO_ENV) {
    super();
    this._env = env;
    this.table = options.table ?? '';
    this.pageSize = options.pageSize ?? 50;
    this.withAccess = options.withAccess ?? true;
    this._captionOption = options.captions;
    this.defaults = {...options.defaults};
    this.isDraft = options.draft ?? false;
    this.isEmpty = options.empty ?? this.isDraft;
    this.query = signal<string | FilterGroup>(options.query ?? '');
    this.search = signal(options.search ?? '');
    this.sort = signal(options.sort ?? '');
    this.deleted = signal<DomainDeletedMode>(options.deleted ?? 'exclude');
    this.live = signal(options.live ?? false);
    this.stale = this._isStale;
    this._liveMs = options.liveMs ?? DomainSource.liveMs;
    this.readOnly = computed(() => this.deleted.value !== 'exclude' || !this._writable.value);
    this._backend = requireBackend(this, backends.domain, 'domain tables');
    this.df = this._df;
    this.error = this._error;
    this.problemRow = this._problemRow;
    // the upper bound of a trash source (which keeps the Delete grant: Restore rides it), or of a
    // table nothing can be written to, over whatever the server answered for the table
    this.access = computed(() => this.readOnly.value ?
      this._access.value.narrow({edit: false, insert: false, ...(this._writable.value ? {} : {delete: false})}) :
      this._access.value);
    this.selection = this._selection;
    this.edit = this._edit;
    this.state = this._state;
    this.total = this._total;
    this._errorStep = computed(() => DomainSource._message(this._error.value));
    this.currentRow = signal<RowView<TRow> | null>(null);
    this._rows = new FrameRows(this._df, this.scope,
      {onWrite: (id, column, value) => this._edit.peek()?.setValue(id, column, value)});
    this.rows = this._rows as RowsLike<RowView> as RowsLike<RowView<TRow>>;
    this.isDirty = computed(() => this._edit.value?.isDirty.value ?? false);
    this.changeCount = computed(() => this._edit.value?.changeCount.value ?? 0);
    this.validity = computed(() => this._edit.value?.validity.value ?? null);
    this.isSaving = computed(() => this._edit.value?.isSaving.value ?? false);
    this.summary = computed(() => this._summary());
    this.session = options.session ?? SharedSession.ambient ?? new SharedSession();
    const leave = this.session.add?.(this);
    DomainSource.live.add(this);
    this.own(() => {
      DomainSource.live.delete(this);
      leave?.();
    });

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
    // the registry keeps the name a spec binds a button to; the method is `stageRestoreSelection`
    this.registerFunction({name: 'restoreSelection', description: 'Stage the restore of every selected row',
      inputs: [], apply: () => this.stageRestoreSelection()});
    this.registerFunction({name: 'newRow', description: 'Add a pristine draft row and make it current',
      inputs: [{name: 'values', type: 'object', nullable: true}],
      apply: (params) => this.newRow(params?.values as RowValues<TRow> | undefined, {pristine: true})});
  }

  /** Phase two: the query may be bound to an input the form declares after this source; the
   * first page is loaded once it is known. A new query is a new collection — unless edits are
   * pending, which a programmatic re-query never drops silently (STATE-CONTRACT H6; the gate
   * that asks is phase 2). */
  start(): void {
    const bound = subBind(this._env, 'query');
    if (bound !== null)
      this.effect(() => this.query.value = DomainSource._queryOf(bound.value));
    const boundLive = subBind(this._env, 'live');
    if (boundLive !== null)
      this.effect(() => this.live.value = boundLive.value === true);
    this.effect(() => {
      this.query.value;
      this.search.value;
      this.sort.value;
      this.deleted.value;
      if (!this.isDirty.peek() && !this._rebinding)
        void this.refresh();
    });
    this.effect(() => {
      this._stopPolling();
      this._failures = 0;
      if (this.live.value) {
        this._timer = setInterval(() => void this._probe(), this._liveMs);
        // a Node timer keeps the process alive; a browser timer id ignores this
        (this._timer as {unref?(): void}).unref?.();
      }
    });
    this.own(() => this._stopPolling());
  }

  get schema(): DomainSchema {
    return this._schema;
  }

  /** The draft a `~new:` id names, in whichever live source holds it — how a picker or a readonly
   * reference shows a parent that is not saved yet. */
  static draftOf(id: string): {source: DomainSource, row: RowView} | undefined {
    for (const source of DomainSource.live) {
      const row = source.rows.byKey(id);
      if (row !== undefined)
        return {source, row};
    }
    return undefined;
  }

  /** A draft over {@link defaults} and `values`, made current — the row a create form binds to;
   * `pristine` keeps it from arming the dirty gate until its first edit. Needs the table: call
   * after the source is ready. */
  newRow(values: RowValues<TRow> = {}, options?: {pristine?: boolean}): RowView<TRow> {
    if (this.readOnly.peek())
      throw new DomainBackendError('forbidden', `${this.table}: ${this._readOnlyReason}`);
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
   * refresh outran — touch nothing. An empty source loads no rows; a draft one starts on a
   * pristine draft. */
  async refresh(): Promise<void> {
    const gen = ++this._gen;
    this._restart = false;
    // the rows about to arrive are the server's again, and the next probe baselines on them
    this._probed = undefined;
    batch(() => {
      this._error.value = undefined;
      this._state.value = 'loading';
    });
    try {
      const table = await (this._ready ??= this._backend.table(this.table));
      if (gen !== this._gen)
        return;
      this._adopt(table);
      if (this.deleted.peek() !== 'exclude')
        DomainSource.requireRestore(table);
      // what the caller may do does not depend on the query: a filter the server refuses must not
      // take New away with the rows. It is known BEFORE the spec is built, because the captions
      // default is derived from it and `Access.readOnly` hides nothing — the cost is zero, the
      // platform caches `access()` per table
      const data = await table.access();
      if (gen !== this._gen)
        return;
      this._access.value = Access.from(data);
      const [frame, total] = await Promise.all([
        table.frame(this._spec(0)), this._noRows ? 0 : this._count(table)]);
      if (gen !== this._gen) {
        frame.dispose();
        return;
      }
      // an edit landed on the live frame while this load was in flight: the edited frame stays
      if (this.isDirty.peek()) {
        frame.dispose();
        this._state.value = 'ready';
        return;
      }
      this._drop();
      this._frame = frame;
      this._loaded = frame.df.rowCount;
      this._done = this._noRows || this._loaded < this.pageSize;
      this._wire(frame);
      batch(() => {
        this._df.value = frame.df;
        this._edit.value = frame.edit;
        this._total.value = total;
        this._state.value = 'ready';
        // cleared where the load LANDS, not where it starts: one that bails over pending changes
        // leaves the rows behind the server, and the mark has to say so
        this._isStale.value = false;
      });
      this._syncCurrent();
      this._syncSelection();
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
    if (frame === undefined || this._restart)
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

  /** Stages the restore of soft-deleted rows into the session's unit of work — they are pending
   * changes until Save, shown as such, dropped by Discard, and guarded by the unsaved gate. The
   * write itself rides `save()`; `DomainTable.restore(id)` is still the way to restore a row with
   * no session behind it. */
  stageRestore(ids: readonly string[]): void {
    const table = this._table;
    if (table === undefined)
      throw new DomainBackendError('not-found', `${this.table}: the table is not loaded yet`);
    DomainSource.requireRestore(table);
    const edit = this._edit.peek();
    if (edit === undefined)
      throw new DomainBackendError('not-found', `${this.table}: the table is not loaded yet`);
    for (const id of ids)
      edit.markRestored(id);
  }

  unstageRestore(ids: readonly string[]): void {
    const edit = this._edit.peek();
    if (edit === undefined)
      return;
    for (const id of ids)
      edit.unmarkRestored(id);
  }

  /** {@link stageRestore} over the frame's selected rows — a trash list's bulk Restore. */
  stageRestoreSelection(): void {
    this.stageRestore(this._selection.peek().map((row) => row.id));
  }

  /** Refuses a trash source, or a restore, over a backend that cannot restore: `restore` is the
   * whole soft-delete lifecycle, and a `deleted` list without it shows rows nothing brings back. */
  static requireRestore(table: DomainTableLike): void {
    if (table.restore === undefined) {
      throw new DomainBackendError('unsupported',
        `${table.address}: the backend does not support restoring deleted rows`);
    }
  }

  /** Through the session — the one Save every button and shortcut runs; the trash source's own
   * refusal lives in {@link check}, which the session runs for every source of the batch. */
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

  /** What keeps this source's batch from being sent — a guard's refusal, then the writer's
   * validity — set as the `error` ("Cannot save: Title is required"; the summary and a paired
   * form show it); null when saveable. The session runs it on every source of the batch it is
   * about to send, dirty or not: a pristine parent a child refers to is inserted by that batch,
   * and the app's rules over it hold. */
  check(): string | null {
    if (this.readOnly.peek() && this.pending().some((r) => r[Rows.STATE] !== 'restored'))
      return this._refuse(this._readOnlyReason);
    const edit = this._edit.peek();
    if (edit === undefined)
      return null;
    batch(() => {
      this._error.value = undefined;
      this._problemRow.value = null;
    });
    for (const guard of this._guards) {
      const problem = guard();
      if (problem !== null)
        return this._refuse(problem);
    }
    const problem = edit.validity.peek();
    return problem === null ? null : this._refuse(this._locate(edit, problem) ?? problem);
  }

  /** The cell the writer's `validity` is about, worded by {@link nameCell} — the writer reports
   * the message alone, and "Cannot save: Value can't be empty" names no field. */
  private _locate(edit: EditState, problem: string): string | null {
    const name = this.nameCell;
    if (name === undefined)
      return null;
    for (const row of this.pending()) {
      for (const prop of this._schema.properties) {
        if (edit.errorOf(row.id, prop.name) !== problem)
          continue;
        this.markProblem(row.id);
        return name(row, prop.name, problem);
      }
    }
    return null;
  }

  /** Every row the writer has pending — new, modified or deleted — whatever the frame's filter
   * shows: what the batch is built from, and so what validation and reference discovery are
   * about. */
  pending(): readonly RowView<TRow>[] {
    return this._rows.pending() as RowView<TRow>[];
  }

  /** After a batch this source took part in landed: the loaded window re-read so paging starts
   * from the offset the server agrees with, the total counted again, and every other live source
   * of the same table reloaded — it is holding the rows this batch just rewrote. */
  async afterSave(): Promise<void> {
    // the batch moved the collection itself: the next probe baselines on what we just wrote
    this._probed = undefined;
    await this._rebase();
    // and the re-base has just re-read the window, so the rows are the server's again — a save
    // clears the stale mark as a reload does, or "Data changed — Refresh" would stand for good
    this._isStale.value = false;
    for (const source of DomainSource.live) {
      if (source !== this && source.table === this.table && !source.isDirty.peek())
        source.refresh().catch((e) => source.fail(e));
    }
    if (!this._noRows && this._table !== undefined)
      this._total.value = await this._count(this._table);
  }

  /** After a landed batch: every draft id it assigned taken out of this source's {@link query} and
   * {@link defaults}, so a child collection built under a draft parent (`fk = "~new:…"`, the
   * master–detail shape) names the row the parent became — the re-read, and every draft added
   * later, are about the real id. Answers whether the query changed — what has rows to re-read;
   * re-reading is the caller's ({@link afterSave}), so the rewrite itself loads nothing. */
  rebind(assigned: Record<string, string>): boolean {
    for (const [column, value] of Object.entries(this.defaults)) {
      const real = Rows.real(assigned, value);
      if (real !== undefined)
        this.defaults[column] = real;
    }
    // the rows too: a source outside the batch may hold a pristine child whose FK cell names the
    // draft the batch turned into a row — it keeps its state and re-reads nothing
    this._edit.peek()?.rebind(assigned);
    const query = this.query.peek();
    const rebound = DomainSource._rebound(query, assigned);
    if (rebound === query)
      return false;
    this._rebinding = true;
    try {
      this.query.value = rebound;
    } finally {
      this._rebinding = false;
    }
    return true;
  }

  /** The rows loaded so far, read again from the first one: a batch that inserted, deleted or
   * moved a row leaves the frame holding other rows than the server would answer at those
   * offsets, and the next page would skip or repeat. The current row and the selection are kept
   * by id — a saved draft under the id it was given. */
  private async _rebase(): Promise<void> {
    const table = this._table;
    const frame = this._frame;
    if (table === undefined || frame === undefined || this._noRows)
      return;
    // the window as it was before the batch applied — `_loaded` is what the server answered, and
    // the frame may have just lost every row of it to a delete
    const window = Math.max(this._loaded, frame.df.rowCount) || this.pageSize;
    const gen = this._gen;
    const selected = this._selection.peek().map((row) => row.id);
    // read before the swap: a row the replacement does not hold takes `currentRow` with it
    const current = this.currentRow.peek()?.id ?? null;
    const replacement = await table.frame(this._spec(0, window));
    if (gen !== this._gen) {
      replacement.dispose();
      return;
    }
    // swapped in place, never through `_drop`: a row's proxy is keyed, reads whatever frame the
    // source holds, and everything bound to it stays bound
    this._unwire?.();
    this._unwire = undefined;
    this._frame = replacement;
    this._loaded = replacement.df.rowCount;
    this._done = this._loaded < window;
    this._wire(replacement);
    batch(() => {
      this._df.value = replacement.df;
      this._edit.value = replacement.edit;
    });
    frame.dispose();
    const at = current === null ? -1 : this._rows.indexOf(current);
    if (at >= 0) {
      replacement.df.currentRowIdx = at;
      this._syncCurrent();
    }
    const selection = replacement.df.selection as {set?(i: number, value: boolean): void} | null | undefined;
    if (typeof selection?.set === 'function') {
      for (const id of selected) {
        const index = this._rows.indexOf(id);
        if (index >= 0)
          selection.set(index, true);
      }
    }
    this._syncSelection();
  }

  /** What the last save threw — the session hands every dirty source the backend's refusal. */
  fail(error: unknown): void {
    this._error.value = error;
  }

  /** The batch landed but the re-read after it did not: the changes are safe, the frame holds
   * the rows the transaction wrote, and only the window the next page would follow is lost — so
   * the next load starts from the first page. Said as such, never as a failed save. */
  markStale(error: unknown): void {
    this._restart = true;
    this._error.value = new DomainBackendError('stale',
      `Saved; refresh failed: ${DomainSource._message(error)}`);
  }

  /** A refusal the backend answered with `false` instead of an exception ("Cannot save: …"):
   * carried by the summary like every other, but marked as already reported — the platform
   * editor balloons its own, and nothing says it twice. */
  refuse(problem: string): void {
    this._error.value = new DomainBackendError(DomainSource.REFUSED, `Cannot save: ${problem}`);
  }

  /** A refusal raised over one row before anything was sent — the session's preflight. Shown the
   * way a guard's is ({@link check}), and it names the row. */
  refuseRow(id: string, problem: string): void {
    this.markProblem(id);
    this._refuse(problem);
  }

  /** Names the row the refusal being raised is about — a guard marks it before it answers. */
  markProblem(id: string | null): void {
    this._problemRow.value = id;
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
      case 'search': return this.search as unknown as Signal<unknown>;
      case 'total': return this.total as unknown as Signal<unknown>;
      case 'state': return this.state as unknown as Signal<unknown>;
      case 'error': return this._errorStep as unknown as Signal<unknown>;
      case 'access': return this.access as unknown as Signal<unknown>;
      case 'isDirty': return this.isDirty as unknown as Signal<unknown>;
      default: return super.bindStep(name);
    }
  }

  bindProps(): BindProp[] {
    return [
      {name: 'rows', type: 'object', description: 'The rows loaded so far, drafts included', default: true},
      {name: 'currentRow', type: 'object', walkable: true, description: 'The row a form edits; its columns are steps'},
      {name: 'search', type: 'string', description: 'The search text over the searchable columns', writable: true},
      {name: 'total', type: 'int', description: 'How many rows match the query'},
      {name: 'state', type: 'string', description: 'idle, loading, ready or error'},
      {name: 'error', type: 'string', description: 'Why the last load or save failed'},
      {name: 'access', type: 'object', description: 'What the caller may do with the table and its fields'},
      {name: 'isDirty', type: 'bool', description: 'Whether there are unsaved changes'},
      {name: 'source', type: 'object', description: 'The source itself — what a domain form, list or picker binds to'},
    ];
  }

  private _refuse(problem: string): string {
    this._error.value = new DomainBackendError('validation', `Cannot save: ${problem}`);
    return problem;
  }

  /** Why nothing may be written through this source while {@link readOnly} holds. */
  private get _readOnlyReason(): string {
    return this._writable.peek() ? 'deleted rows are read-only until they are restored' :
      'the table does not accept writes';
  }

  private _adopt(table: DomainTableLike): void {
    if (this._table === table)
      return;
    this._table = table;
    this._writable.value = table.support.transaction;
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
      this._selection.value = [];
      this.currentRow.value = null;
    });
  }

  /** Follows the writer's changes and the frame's current row and selection. */
  private _wire(frame: DomainFrameLike): void {
    const bump = () => this._version.value = this._version.peek() + 1;
    const subs = [
      // the editor removes saved deletes and re-keys saved drafts without a frame event (H10): the
      // row projection follows its onChanged instead
      frame.edit.onChanged.subscribe(() => {
        // a refusal stands in the summary until the next edit or discard takes it back
        batch(() => {
          this._error.value = undefined;
          this._problemRow.value = null;
        });
        this._rows.rebuild();
        bump();
      }),
      // the exact post-save re-point: a current draft is the same row under the id it was given
      frame.edit.onSaved.subscribe(({assigned}) => {
        const row = this.currentRow.peek();
        const real = row === null ? undefined : Rows.real(assigned, row.id);
        if (real !== undefined)
          this.currentRow.value = this.rows.byKey(real) ?? null;
      }),
      frame.df.onCurrentRowChanged.subscribe(() => this._syncCurrent()),
      frame.df.onSelectionChanged.subscribe(() => this._syncSelection()),
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
    this.currentRow.value = key === undefined ? null : this.rows.byKey(key) ?? null;
  }

  private _syncSelection(): void {
    const d = this._df.peek();
    const selected = d?.selection as {get?(i: number): boolean} | null | undefined;
    const rows: RowView<TRow>[] = [];
    if (d !== undefined && typeof selected?.get === 'function') {
      for (let i = 0; i < d.rowCount; i++) {
        const row = selected.get(i) ? this.rows.byKey(this._rows.keyAt(i)!) : undefined;
        if (row !== undefined)
          rows.push(row);
      }
    }
    this._selection.value = rows;
  }

  private _count(table: DomainTableLike): Promise<number> {
    return table.count(this._scope());
  }

  /** One poll of a `live` source: the collection counted and dated in one request, against what
   * the last poll saw. A moved pair refreshes while there is nothing to lose and marks the source
   * {@link stale} while there is. Nothing is asked while the tab is hidden (throttled, and no one
   * is reading it) or mid write-back (the batch moves `updated_on` itself), and a backend with no
   * probe is not polled at all.
   *
   * A probe that fails leaves the last known pair standing, so a one-tick outage cannot baseline
   * a change away; after {@link liveFailures} failures in a row the timer gives up rather than
   * poll forever at something permanently broken (a 403, a table without `updated_on`). `live`
   * stays what the caller asked for — setting it false and true again starts a new timer. */
  private async _probe(): Promise<void> {
    const table = this._table;
    if (table === undefined || table.probe === undefined || this._probing || this._noRows ||
        this.isSaving.peek() || (typeof document !== 'undefined' && document.hidden))
      return;
    this._probing = true;
    const gen = this._gen;
    try {
      const now = await table.probe(this._scope());
      if (gen !== this._gen)
        return;
      this._failures = 0;
      const before = this._probed;
      this._probed = now;
      if (before === undefined || (before.count === now.count && before.last === now.last))
        return;
      if (this.isDirty.peek() || this.isSaving.peek())
        this._isStale.value = true;
      else
        await this.refresh();
    } catch {
      // the baseline stands: the next good probe compares against the last pair that WAS answered,
      // so a change made during an outage is still seen
      if (++this._failures >= DomainSource.liveFailures)
        this._stopPolling();
    } finally {
      this._probing = false;
    }
  }

  private _stopPolling(): void {
    if (this._timer !== undefined)
      clearInterval(this._timer);
    this._timer = undefined;
  }

  /** The access a row is written under — its own (`Access.row`: a draft under `insert`, an existing
   * row as its `~can_*` columns say), no row the table's. */
  private _view(row: RowView | null): Access {
    const access = this.access.peek();
    return row === null ? access : access.row(row);
  }

  /** Loads no rows: an empty (or draft) source, and one whose query still names a draft id. */
  private get _noRows(): boolean {
    return this.isEmpty || DomainSource._namesDraft(this.query.peek());
  }

  private _spec(offset: number, limit = this.pageSize): DomainQueryLike {
    // a query naming a draft id is not sent at all: no saved row can match it, and the server
    // refuses the `~new:` literal on a uuid column — the rows arrive when `rebind` puts the
    // assigned id in
    const deleted = this.deleted.peek();
    if (DomainSource._namesDraft(this.query.peek())) {
      return {limit: 0, offset, withAccess: this.withAccess,
        ...(deleted === 'exclude' ? {} : {deleted})};
    }
    const sort = this.sort.peek();
    const captions = this._captions();
    return {...this._scope(), ...(sort === '' ? {} : {sort}),
      limit: this.isEmpty ? 0 : limit, offset, withAccess: this.withAccess,
      ...(captions.length === 0 ? {} : {captions})};
  }

  /** The ref columns worth a caption: a domain-table reference (`semType` is `<schema>.<table>`) the
   * caller may see. `User`/`Group` refs are NOT domain tables — the server refuses a caption for one
   * — and a column the caller's access hides would fail the whole query with the same refusal. None
   * over a backend that does not project them (`support.captions`): a form resolves each per row. */
  private _captions(): string[] {
    if (this._table?.support.captions !== true)
      return [];
    const listed = this._captionOption;
    const access = this._access.peek();
    const refs = this._schema.properties.filter((p) => REF_ADDRESS.test(p.semType ?? '') &&
      access.field(p.name!) !== 'hidden').map((p) => p.name!);
    return listed === undefined ? refs : listed.filter((name) => refs.includes(name));
  }

  /** What every read of this source is scoped to — the one object `_spec`, `count` and `probe`
   * all take, so none of them can forward the query and forget the search or the trash mode. */
  private _scope(): DomainReadScope {
    const q = this.query.peek();
    const search = this.search.peek();
    const deleted = this.deleted.peek();
    return {filter: typeof q === 'string' ? (q === '' ? undefined : q) : Filters.toDomainTree(q),
      ...(search === '' ? {} : {search}), ...(deleted === 'exclude' ? {} : {deleted})};
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
          (row as RowView)[column] = value;
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
    if (changes > 0) {
      const deleted = items.filter((r) => r[Rows.STATE] === 'deleted').length;
      if (deleted === changes)
        return `${plural(deleted, 'deletion', 'deletions')} pending`;
      const restored = items.filter((r) => r[Rows.STATE] === 'restored').length;
      if (restored === changes)
        return `${plural(restored, 'restore', 'restores')} pending`;
      return plural(changes, 'unsaved change', 'unsaved changes');
    }
    const state = this.state.value;
    const singular = this._schema.info.singularName.toLowerCase() || 'row';
    // a draft is a row, not a row of the table yet; a draft source is a create form until the
    // draft is saved — then it holds the row it made, and says so, until the next New
    const loaded = items.filter((r) => !Rows.isDraft(r)).length;
    if (this.isDraft && (loaded === 0 || (this.currentRow.value !== null && Rows.isDraft(this.currentRow.value))))
      return state === 'ready' ? `New ${singular}` : state === 'loading' ? 'Loading…' : '';
    if (this.isDraft)
      return `${singular.charAt(0).toUpperCase()}${singular.slice(1)} saved`;
    if (state !== 'ready' && loaded === 0)
      return state === 'loading' ? 'Loading…' : '';
    const total = this.total.value;
    const many = this._schema.info.pluralName.toLowerCase() || 'rows';
    if (this.deleted.value === 'only')
      return plural(total ?? loaded, `deleted ${singular}`, `deleted ${many}`);
    return !this.isDraft && total !== null && loaded < total ?
      `${loaded.toLocaleString()} of ${total.toLocaleString()}` :
      plural(loaded, singular, many);
  }

  private static _namesDraft(query: string | FilterGroup): boolean {
    if (typeof query === 'string')
      return DRAFT_LITERAL.test(query);
    let found = false;
    Filters.walk(query, (n) => {
      if (Filters.isGroup(n) || n.value === undefined)
        return;
      for (const v of Array.isArray(n.value) ? n.value : [n.value])
        found = found || (typeof v === 'string' && Rows.isDraft(v));
    });
    return found;
  }

  /** The query with every draft id `assigned` names replaced by the id it was given — the string
   * form by the `~new:` literal, the tree form by the value; the query itself when none is in it. */
  private static _rebound(query: string | FilterGroup, assigned: Record<string, string>): string | FilterGroup {
    if (typeof query === 'string') {
      let text = query;
      // the quoted literal only: `~new:…` inside a longer value is that value's own text
      for (const [draft, id] of Object.entries(assigned)) {
        for (const quote of ['"', '\''])
          text = text.split(`${quote}${draft}${quote}`).join(`${quote}${id}${quote}`);
      }
      return text;
    }
    let hit = false;
    const swap = (v: FilterScalar): FilterScalar => {
      const real = Rows.real(assigned, v);
      if (real === undefined)
        return v;
      hit = true;
      return real;
    };
    const copy = (n: FilterNode): FilterNode => Filters.isGroup(n) ? {...n, nodes: n.nodes.map(copy)} :
      n.value === undefined ? n : {...n, value: Array.isArray(n.value) ? n.value.map(swap) : swap(n.value)};
    const rebound = copy(query) as FilterGroup;
    return hit ? rebound : query;
  }

  private static _message(error: unknown): string {
    return error === undefined ? '' : error instanceof Error ? error.message : String(error);
  }

  private static _queryOf(value: unknown): string | FilterGroup {
    return typeof value === 'string' ? value : value === null || value === undefined ? '' : value as FilterGroup;
  }
}
