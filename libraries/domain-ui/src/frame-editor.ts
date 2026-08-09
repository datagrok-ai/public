/**
 * {@link DomainFrameEditor} — batch editing of domain-table rows on a DataFrame,
 * with the editing state living in the frame itself (the `~state` / `~changes` /
 * `~errors` service columns) rather than in a side store.
 *
 * @module frame-editor
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

/** Row state column: `'' | 'new' | 'modified' | 'deleted'`. */
export const STATE_COLUMN = '~state';
/** JSON column holding the ORIGINAL values of changed cells only (sparse). */
export const CHANGES_COLUMN = '~changes';
/** JSON column holding per-cell {@link DomainCellError}s. */
export const ERRORS_COLUMN = '~errors';

/** The three service columns {@link DomainFrameEditor} attaches — every one of
 * them tagged out of binary AND csv export, so the editing state can never reach
 * a saved project, an export, an upload, or a `batch()` fed from the frame. */
export const SERVICE_COLUMNS: string[] = [STATE_COLUMN, CHANGES_COLUMN, ERRORS_COLUMN];

/** Per-row editing state, as stored in {@link STATE_COLUMN}. */
export type DomainRowState = '' | 'new' | 'modified' | 'deleted';

/** One cell-level problem, as stored in {@link ERRORS_COLUMN}: `'error'` blocks
 * {@link DomainFrameEditor.save}, `'conflict'` (a dismissed version conflict)
 * only marks the cell. */
export interface DomainCellError {
  message: string;
  kind: 'error' | 'conflict';
}

/** What one {@link DomainFrameEditor.save} wrote (all of it in ONE transaction). */
export interface DomainSaveResult {
  inserted: number;
  updated: number;
  deleted: number;
}

/** Options of {@link DomainFrameEditor.attach} / {@link DomainFrameEditor.create}. */
export interface DomainFrameEditorOptions {
  /** The query the frame came from — what {@link DomainFrameEditor.refresh}
   * re-runs, and what {@link DomainFrameEditor.create} runs to build it. */
  query?: DG.DomainQuerySpec;
  /** Pre-probed capabilities, to avoid a second round trip when the caller
   * already has them (a {@link DomainGrid} passes its own). */
  capabilities?: DG.DomainTableCapabilities;
}

/** Marks a first-touch whose original value could not be captured (see
 * {@link DomainFrameEditor.beginEdit}) — the cell still saves, it just cannot be
 * reverted. Object-shaped on purpose: no scalar wire value can collide with it. */
const UNKNOWN_ORIGINAL = {'~unknownOriginal': true};

/** How many pre-edit row snapshots {@link DomainFrameEditor.beginEdit} keeps.
 * One would do for the grid (a cell is edited only while it is current), a
 * handful covers an editor that commits after the focus already moved on. */
const SNAPSHOT_LIMIT = 8;

function _isUnknown(x: any): boolean {
  return x != null && typeof x === 'object' && x['~unknownOriginal'] === true;
}

/** A value as it travels on the wire: dayjs/Date become ISO-8601, non-finite
 * numbers and DG's null sentinels become null. Comparisons and `~changes`
 * entries use this form, so a revert compares like the server would. */
function toWire(v: any): any {
  if (v == null)
    return null;
  if (typeof v === 'object' && typeof v.toISOString === 'function')
    return v.toISOString();
  if (typeof v === 'number' && !isFinite(v))
    return null;
  return v;
}

function wireEquals(a: any, b: any): boolean {
  if (a === b)
    return true;
  if (a == null || b == null)
    return a == null && b == null;
  if (typeof a === 'object' || typeof b === 'object')
    return JSON.stringify(a) === JSON.stringify(b);
  return false;
}

/** The value's platform type name, resolved the dart2js-safe way (`1` is an int
 * here, not the `double` a dart2js `is` check reports) — see
 * {@link validateCellValue}. */
function valueTypeOf(v: any): string {
  if (typeof v === 'string')
    return DG.TYPE.STRING;
  if (typeof v === 'boolean')
    return DG.TYPE.BOOL;
  if (typeof v === 'number')
    return Number.isInteger(v) ? DG.TYPE.INT : DG.TYPE.FLOAT;
  if (Array.isArray(v))
    return DG.TYPE.LIST;
  if (v != null && typeof v === 'object' && typeof v.toISOString === 'function')
    return DG.TYPE.DATE_TIME;
  return DG.TYPE.MAP;
}

/** Numbers formatted the way the server VM prints them, so the messages below
 * match the server's byte-for-byte (the VM prints whole doubles as '1.0'). */
function vmNum(v: number): string {
  return v % 1 === 0 ? v.toFixed(1) : `${v}`;
}

function isNumerical(p: DG.Property): boolean {
  const t = p.propertyType;
  return t === DG.TYPE.INT || t === DG.TYPE.FLOAT || t === DG.TYPE.NUM || t === DG.TYPE.BIG_INT;
}

function validateNumeric(p: DG.Property, v: number): string | null {
  if (p.propertyType === DG.TYPE.INT && v % 1 !== 0)
    return 'Types differ. Expected: int, passed: double';
  const passed = p.propertyType === DG.TYPE.INT ? `${v}` : vmNum(v);
  if (p.min != null && v < p.min)
    return `Value should not be less than ${vmNum(p.min)}, passed: ${passed}`;
  if (p.max != null && v > p.max)
    return `Value should not be more than ${vmNum(p.max)}, passed: ${passed}`;
  return null;
}

/**
 * Validates one cell value against its registry {@link DG.Property} — the same
 * constraints the server re-runs on write, producing the server's exact message
 * texts so an inline marker and a rejected save read identically.
 *
 * Numeric properties take a dart2js-safe path (under dart2js every whole number
 * `is double`, which would reject every integer), and `string_list` values are
 * left to the server, which coerces them. Returns null when the value is fine.
 */
export function validateCellValue(p: DG.Property, value: any): string | null {
  if (!p.nullable && (value == null || value === ''))
    return "Value can't be empty";
  if (value == null || p.propertyType === DG.TYPE.STRING_LIST)
    return null;
  const type = valueTypeOf(value);
  if (isNumerical(p))
    return typeof value === 'number' ? validateNumeric(p, value)
      : `Numerical value expected, passed: ${type}`;
  const expected = p.propertyType === DG.TYPE.NUM ? DG.TYPE.FLOAT : p.propertyType;
  if (type !== expected)
    return `Types differ. Expected: ${p.propertyType}, passed: ${type}`;
  const choices = p.choices;
  if (choices != null && choices.length > 0 && !choices.includes(value))
    return `"${value}" is not one of (${choices.join(', ')})`;
  return null;
}

/** Whether [p] addresses another row rather than carrying a value of its own — a
 * `ref` column (semType `'<schema>.<table>'`) or a `user`/`group` column. Those
 * hold uuids: a picker is their editing path, not a text cell. */
export function isReferenceProperty(p: DG.Property): boolean {
  const semType = `${p.semType ?? ''}`;
  return semType === 'User' || semType === 'Group' || /^[^.]+\.[^.]+$/.test(semType);
}

/** Whether [p] addresses another DOMAIN row — a `ref` column (semType
 * `'<schema>.<table>'`) alone. User / group columns hold uuids too, but the
 * platform's own anchored picker edits them in a grid, so they are not excluded
 * there; every other path ({@link isReferenceProperty}) still treats all three
 * alike. */
export function isDomainRefProperty(p: DG.Property): boolean {
  return /^[^.]+\.[^.]+$/.test(`${p.semType ?? ''}`);
}

/**
 * A widget or page that answers for a set of {@link DomainFrameEditor}s.
 *
 * The dirty rollup is domain-specific — it does not belong on `DG.Widget`, and it
 * cannot be derived from `Widget.children` (which is a DOM view, so it is order-
 * and reparenting-fragile). This is the one structural remainder of the widget
 * model: anything that owns pending changes says so by exposing them here, and a
 * container ({@link AppView}, `domains.view`) rolls its children's up.
 */
export interface IEditorHost {
  /** Every editor whose pending changes this object answers for. Read on demand:
   * an editor that appears after the widget was constructed must show up here. */
  readonly editors: DomainFrameEditor[];
}

/** The editors [widget] answers for, or none — the duck-typed read of
 * {@link IEditorHost} a container uses on children of any class. */
export function editorsOf(widget: any): DomainFrameEditor[] {
  const editors = widget?.editors;
  return Array.isArray(editors) ? editors : [];
}

/** A domain client with its row/insert/column/expand generics erased — what the
 * REFLECTIVE widgets take. A TYPED client fits it, so a generated `<Table>Ui` (or
 * a `DomainTable` handle) keeps its own generics and is still usable here. */
export type AnyDomainTableClient = DG.DomainTableClient<any, any, any, any>;

/**
 * The prefetched table context a SYNCHRONOUS widget factory needs: the typed
 * client plus the registry metadata and capabilities, resolved once by
 * `domains.table(...)` (which is why `form()`, `grid()` and friends need no await).
 *
 * Capabilities are a SNAPSHOT taken when the context was acquired — the contract
 * the grid has always had, moved one level up: a later grant change (or a
 * `grok.dapi.domains.invalidateUiCaches()`) does not reach widgets already built
 * from this context; re-acquire the handle to re-gate.
 */
export interface IDomainTableContext {
  readonly client: AnyDomainTableClient;
  /** Registry {@link DG.Property} metadata of the table's declared columns. */
  readonly properties: DG.Property[];
  readonly info: DG.DomainTableInfo;
  readonly capabilities: DG.DomainTableCapabilities;
  /** `'<schema>.<table>'`. */
  readonly table: string;
}

/**
 * Resolves [client]'s registry metadata and the caller's capabilities in ONE
 * round of requests — the async boundary every synchronous widget factory sits
 * behind (`domains.table()`, `DomainGrid.create()`, `EntityListWidget.create()`).
 */
export async function acquireDomainContext(
  client: AnyDomainTableClient): Promise<IDomainTableContext> {
  const address = `${client.schema}.${client.table}`;
  const [properties, info, capabilities] = await Promise.all([
    grok.dapi.domains.registry.rowProperties(address),
    grok.dapi.domains.registry.tableInfo(address),
    client.capabilities(),
  ]);
  return {client: client, properties: properties, info: info,
    capabilities: capabilities, table: address};
}

/** An empty frame of [properties]' columns — what the local (no round trip)
 * editors of {@link DomainFrameEditor.forContext} / `forRows` start from. */
function _emptyFrame(properties: DG.Property[]): DG.DataFrame {
  const df = DG.DataFrame.create(0);
  for (const p of properties)
    df.columns.add(DG.Column.fromType(columnTypeOf(p), p.name, 0));
  return df;
}

/** The frame column type that holds [p]'s values — what {@link DomainFrameEditor.forRows}
 * builds a local frame from (a `queryDf` result carries the server's own types). */
function columnTypeOf(p: DG.Property): DG.ColumnType {
  switch (p.propertyType) {
    case DG.TYPE.INT: return DG.COLUMN_TYPE.INT;
    case DG.TYPE.BIG_INT: return DG.COLUMN_TYPE.BIG_INT;
    case DG.TYPE.FLOAT:
    case DG.TYPE.NUM: return DG.COLUMN_TYPE.FLOAT;
    case DG.TYPE.BOOL: return DG.COLUMN_TYPE.BOOL;
    case DG.TYPE.DATE_TIME: return DG.COLUMN_TYPE.DATE_TIME;
    default: return DG.COLUMN_TYPE.STRING;
  }
}

/** One pending write of {@link DomainFrameEditor.buildOps}, with the frame row it
 * came from (how a failing `opIndex` finds its row). */
export interface DomainPendingOp {
  op: DG.DomainTransactionOp;
  row: number;
}

/**
 * THE single writer of a domain frame's editing state.
 *
 * It wraps a DataFrame produced by `table.queryDf(...)` and attaches three
 * invisible service columns — {@link STATE_COLUMN}, {@link CHANGES_COLUMN},
 * {@link ERRORS_COLUMN} — that hold everything about the pending batch: which
 * rows are new/modified/deleted, the ORIGINAL value of every changed cell, and
 * the per-cell validation errors. Grids, forms and the save pipeline all read
 * that one state; nothing keeps a parallel store.
 *
 * ```ts
 * const editor = await DomainFrameEditor.create(grok.dapi.domains.table('grit.issue'));
 * editor.setValue(0, 'title', 'New title');     // tracked, validated, highlighted
 * await editor.save();                          // ONE /transaction
 * ```
 *
 * Every service column is tagged out of binary AND csv export, so the state is
 * memory-only: a saved project, `toByteArray()`, `toCsv()`, an export or a
 * `batch()` upload built from the frame never carry it.
 *
 * **Writing.** Go through {@link setValue} (programmatic) or
 * {@link beginEdit} + {@link commitEdit} (an in-grid edit, where the grid has
 * already written the cell). Writing a cell directly on the DataFrame bypasses
 * the tracking and the value is silently NOT saved.
 *
 * **Deleted rows stay in the frame** and are hidden by ANDing them out of the
 * filter bitset on every filter recomputation, so undoing a delete
 * ({@link unmarkDeleted}) is trivial and row order never moves.
 *
 * **Refreshing discards edits — BY DESIGN.** {@link refresh} re-runs the query
 * and rebuilds the frame and its state from scratch; there is no merge and never
 * will be. Deciding whether it is safe to refresh is the CALLER's job: read
 * {@link isDirty} / subscribe to {@link onDirtyChanged} and prompt (save /
 * discard / cancel) before calling it. A component that refreshes on a timer or
 * on a route change without that check WILL eat a user's batch edits.
 */
export class DomainFrameEditor {
  private _df: DG.DataFrame;
  private _subs: rxjs.Subscription[] = [];
  private _properties: DG.Property[];
  private _propByName = new Map<string, DG.Property>();
  private _snapshots = new Map<number, {[column: string]: any}>();
  /** Parsed {@link CHANGES_COLUMN} / {@link ERRORS_COLUMN} objects, per column
   * and row: the JSON columns are read on every cell paint and on every change
   * count, so re-parsing them would cost a full parse per cell per frame. Every
   * writer of those columns goes through {@link _setJson} / {@link _clearRowState},
   * and anything that shifts row indices resets the whole cache. */
  private _parsed = new Map<string, Map<number, any>>();
  /** Running {@link changeCount}: kept per row so a keystroke costs one row
   * recount instead of a scan of the frame. */
  private _contributions = new Map<number, number>();
  /** `'new'` rows nobody has written to yet — see {@link addRow}. Reset with the
   * other per-row caches, which is the safe direction: a pristine row that lost
   * its mark counts as a change, never the other way round. */
  private _pristine = new Set<number>();
  private _changeCount = 0;
  private _suspend = false;
  private _saving = false;
  private _saveDepth = 0;
  private _dirty = false;
  private _query?: DG.DomainQuerySpec;
  private _nameColumn: string | null;

  private readonly _onChanged = new rxjs.Subject<DomainFrameEditor>();
  private readonly _onDirtyChanged = new rxjs.Subject<boolean>();
  private readonly _onSavingChanged = new rxjs.Subject<boolean>();
  private readonly _onSaved = new rxjs.Subject<DomainSaveResult>();
  private readonly _onConflict = new rxjs.Subject<DG.DomainVersionConflictError>();
  private readonly _onRefreshed = new rxjs.Subject<DG.DataFrame>();

  private constructor(
    /** The table the frame's rows belong to. */
    public readonly client: DG.DomainTableClient,
    /** Effective capabilities of the current user, SNAPSHOT when the editor was
     * created — what read-only degradation and the writable-column payload filter
     * derive from. A later `grok.dapi.domains.invalidateUiCaches()` (or a grant
     * change) does NOT reach an existing editor: re-create it to pick the new
     * permissions up. */
    public readonly capabilities: DG.DomainTableCapabilities,
    df: DG.DataFrame, properties: DG.Property[], nameColumn: string | null,
    query?: DG.DomainQuerySpec) {
    this._properties = properties;
    for (const p of properties)
      this._propByName.set(p.name, p);
    this._nameColumn = nameColumn;
    this._query = query;
    this._df = df;
    this._bind(df);
  }

  /** Attaches the editing state to an EXISTING frame of [client]'s rows (a
   * `queryDf` result). Pass `options.query` so {@link refresh} knows what to
   * re-run. */
  static async attach(dataFrame: DG.DataFrame, client: DG.DomainTableClient,
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    const address = `${client.schema}.${client.table}`;
    const [properties, capabilities, info] = await Promise.all([
      grok.dapi.domains.registry.rowProperties(address),
      options?.capabilities != null ? Promise.resolve(options.capabilities) : client.capabilities(),
      grok.dapi.domains.registry.tableInfo(address),
    ]);
    return new DomainFrameEditor(client, capabilities, dataFrame, properties,
      info.nameColumn, options?.query);
  }

  /** Runs `options.query` (everything the caller may see, by default) and
   * attaches to the resulting frame. */
  static async create(client: DG.DomainTableClient,
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    const df = await client.queryDf(options?.query ?? {});
    return DomainFrameEditor.attach(df, client, options);
  }

  /**
   * SYNCHRONOUS: an editor over an EMPTY frame of [context]'s declared columns —
   * everything {@link attach} awaits is already in the prefetched context, so a
   * widget factory can build its editor without a round trip and load the rows
   * afterwards ({@link refresh}, which replaces the frame and fires
   * {@link onRefreshed}).
   */
  static forContext(context: IDomainTableContext,
    options?: DomainFrameEditorOptions): DomainFrameEditor {
    return new DomainFrameEditor(context.client,
      options?.capabilities ?? context.capabilities, _emptyFrame(context.properties),
      context.properties, context.info.nameColumn, options?.query);
  }

  /**
   * Builds an editor over rows that do NOT come from the server: a frame is built
   * locally from the table's declared columns and every entry of [values] is added
   * as a `'new'` row — the INSERT path of a form
   * (`domains.table(...).form({values})`), with no query round trip.
   *
   * The added rows are PRISTINE: an unsaved row the user has not written to yet
   * is not a pending change (see {@link addRow}), so an untouched New form has
   * nothing to prompt about.
   *
   * Everything else is identical to {@link create}: the same single-writer model,
   * the same service columns and export tags, the same validation, and a
   * {@link save} that writes the batch as one transaction. {@link refresh} has
   * nothing to re-run unless `options.query` says otherwise — and re-running it
   * would replace these rows, which is why a form never refreshes.
   */
  static async forRows(client: DG.DomainTableClient, values: {[column: string]: any}[],
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    const properties = await grok.dapi.domains.registry.rowProperties(
      `${client.schema}.${client.table}`);
    const editor = await DomainFrameEditor.attach(_emptyFrame(properties), client, options);
    for (const row of values ?? [])
      editor.addRow(row, {pristine: true});
    return editor;
  }

  /** The frame being edited. It is REPLACED by {@link refresh} — re-read it (or
   * subscribe to {@link onRefreshed}) instead of caching it. */
  get dataFrame(): DG.DataFrame { return this._df; }

  /** `'<schema>.<table>'`. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** Registry {@link DG.Property} metadata of the table's declared columns. */
  get properties(): DG.Property[] { return this._properties; }

  /** The query {@link refresh} re-runs. */
  get query(): DG.DomainQuerySpec | undefined { return this._query; }

  /** Whether anything is pending (a changed cell, a new row, a deleted row). */
  get isDirty(): boolean { return this._dirty; }

  /** Whether a {@link save} is in flight. While it is, the editor refuses every
   * write, {@link discard} and {@link refresh} — see {@link save}. */
  get isSaving(): boolean { return this._saving; }

  /** Number of pending cell changes — what a "N unsaved changes" bar shows. */
  get changeCount(): number { return this._changeCount; }

  /** Fires on every service-state write — the repaint hook for a grid. */
  get onChanged(): rxjs.Observable<DomainFrameEditor> { return this._onChanged; }
  /** Fires when {@link isDirty} flips — what a caller's refresh policy listens to. */
  get onDirtyChanged(): rxjs.Observable<boolean> { return this._onDirtyChanged; }
  /** Fires when {@link isSaving} flips — what a grid locks its editing on. */
  get onSavingChanged(): rxjs.Observable<boolean> { return this._onSavingChanged; }
  /** Fires after a successful {@link save}. */
  get onSaved(): rxjs.Observable<DomainSaveResult> { return this._onSaved; }
  /** Fires when a save hits a version conflict, BEFORE the standard dialog. */
  get onConflict(): rxjs.Observable<DG.DomainVersionConflictError> { return this._onConflict; }
  /** Fires with the NEW frame after {@link refresh} rebuilt it. */
  get onRefreshed(): rxjs.Observable<DG.DataFrame> { return this._onRefreshed; }

  // ─────────────────────── state accessors ─────────────────────────

  /** Editing state of [row]. */
  stateOf(row: number): DomainRowState {
    return (this._col(STATE_COLUMN).get(row) ?? '') as DomainRowState;
  }

  /** ORIGINAL values of [row]'s changed cells, keyed by column (empty when the
   * row is unchanged; always empty for a `'new'` row — all of its values are
   * new). Read-only: the object is the editor's own cached parse, and writing to
   * it changes nothing on the frame. */
  changesOf(row: number): {[column: string]: any} {
    return this._json(CHANGES_COLUMN, row);
  }

  /** Per-cell problems of [row], keyed by column. Read-only, see
   * {@link changesOf}. */
  errorsOf(row: number): {[column: string]: DomainCellError} {
    return this._json(ERRORS_COLUMN, row);
  }

  /** Whether the cell carries a pending change (what highlighting keys on). */
  isChanged(row: number, column: string): boolean {
    return this.stateOf(row) === 'new' || column in this.changesOf(row);
  }

  /** The cell's problem, or null. */
  errorOf(row: number, column: string): DomainCellError | null {
    return this.errorsOf(row)[column] ?? null;
  }

  // ─────────────────────── writing (single writer) ─────────────────────────

  /** Snapshots [row]'s current values so a following {@link commitEdit} knows
   * what the cell held BEFORE the edit. A grid calls this when the cell becomes
   * current — an edit can only start there. Without a snapshot the edit is still
   * tracked and saved, it just cannot be reverted. */
  beginEdit(row: number): void {
    if (row < 0 || row >= this._df.rowCount || this._snapshots.has(row))
      return;
    const values: {[column: string]: any} = {};
    for (const p of this._properties)
      if (this._df.columns.contains(p.name))
        values[p.name] = this._wire(row, p.name);
    if (this._snapshots.size >= SNAPSHOT_LIMIT)
      this._snapshots.delete(this._snapshots.keys().next().value as number);
    this._snapshots.set(row, values);
  }

  /** Writes [value] into the cell AND tracks it — the programmatic write path
   * (a form field, a paste, a fill-down). Refused while a {@link save} is in
   * flight. */
  setValue(row: number, column: string, value: any): void {
    if (this._busy('editing'))
      return;
    const original = column in this.changesOf(row) ? undefined : this._wire(row, column);
    this._write(() => this._df.set(column, row, value));
    this._track(row, column, original);
  }

  /** Tracks a cell the GRID already wrote (its `onCellValueEdited` path); the
   * original comes from the {@link beginEdit} snapshot. */
  commitEdit(row: number, column: string): void {
    if (this._suspend || this._busy('editing'))
      return;
    const snapshot = this._snapshots.get(row);
    const original = snapshot != null && column in snapshot ? snapshot[column] : UNKNOWN_ORIGINAL;
    this._track(row, column, original);
  }

  /**
   * Appends a new, unsaved row (state `'new'`), optionally prefilled; returns
   * its index (-1 when refused because a {@link save} is in flight).
   *
   * `options.pristine` adds it as a row nobody has written to YET: it is part of
   * the batch a {@link save} writes, but it contributes NOTHING to
   * {@link changeCount} / {@link isDirty} until the first {@link setValue} or
   * {@link commitEdit} — the "pristine until touched" contract of an insert form,
   * whose untouched (however prefilled) row must not arm the unsaved-changes gate.
   * A row added by a USER gesture (the grid's Add row) is pending immediately,
   * which is the default.
   */
  addRow(values?: {[column: string]: any}, options?: {pristine?: boolean}): number {
    if (this._busy('adding a row'))
      return -1;
    const row = this._df.rowCount;
    this._write(() => {
      this._df.rows.addNew();
      if (values != null)
        for (const name of Object.keys(values))
          if (this._df.columns.contains(name))
            this._df.set(name, row, values[name]);
    });
    // After the write: adding a row fires onRowsAdded, whose cache reset would
    // otherwise drop the mark this line makes.
    if (options?.pristine === true)
      this._pristine.add(row);
    this._col(STATE_COLUMN).set(row, 'new', false);
    this._recount(row);
    this._validateRow(row);
    // The row was appended while the filter was already computed — recompute it,
    // or the new row is invisible in every filtered view of the frame.
    this._df.rows.requestFilter();
    this._fire();
    return row;
  }

  /** Marks rows deleted: they stay in the frame (order untouched) and are
   * excluded from the filter until {@link save} removes them for real. */
  markDeleted(rows: number | number[]): void {
    this._setDeleted(rows, true);
  }

  /** Undoes {@link markDeleted}, restoring whatever the row was before — a row
   * added in this batch goes back to `'new'`, an edited one back to
   * `'modified'`. */
  unmarkDeleted(rows: number | number[]): void {
    this._setDeleted(rows, false);
  }

  /** Restores one cell to its original value and drops its change entry.
   * Refused while a {@link save} is in flight. */
  revertCell(row: number, column: string): void {
    if (this._busy('reverting'))
      return;
    const changes = this.changesOf(row);
    if (!(column in changes))
      return;
    const original = changes[column];
    if (_isUnknown(original)) {
      grok.log.warning(`${this.table}: cannot revert ${column} — its original value was not captured`);
      return;
    }
    this._write(() => this._df.set(column, row, original));
    delete changes[column];
    this._setJson(CHANGES_COLUMN, row, changes);
    this._setError(row, column, null);
    this._recomputeState(row);
    this._fire();
  }

  /** {@link revertCell} for every changed cell of [row]. */
  revertRow(row: number): void {
    if (this._busy('reverting'))
      return;
    for (const column of Object.keys(this.changesOf(row)))
      this.revertCell(row, column);
  }

  /** Drops the whole pending batch: changed cells go back to their originals,
   * new rows are removed, deleted rows are restored. Refused while a
   * {@link save} is in flight — removing rows under the transaction would make
   * its results land on the wrong ones. */
  discard(): void {
    if (this._busy('discarding'))
      return;
    this._write(() => {
      for (let row = this._df.rowCount - 1; row >= 0; row--) {
        if (this.stateOf(row) === 'new') {
          this._df.rows.removeAt(row, 1, false);
          continue;
        }
        const changes = this.changesOf(row);
        for (const column of Object.keys(changes))
          if (!_isUnknown(changes[column]))
            this._df.set(column, row, changes[column]);
        this._clearRowState(row);
      }
    });
    this._resetCaches();
    this._df.rows.requestFilter();
    this._fire();
  }

  /** Re-runs every cell validator over the pending batch; returns the number of
   * blocking (`kind: 'error'`) cells. Call it before offering Save when values
   * arrived from outside {@link setValue}. Refused (reporting the CURRENT count)
   * while a {@link save} is in flight — it writes the state columns like every
   * other mutator. */
  validate(): number {
    if (this._busy('validating'))
      return this.errorCount;
    for (let row = 0; row < this._df.rowCount; row++)
      if (this.stateOf(row) !== '')
        this._validateRow(row);
    this._fire();
    return this.errorCount;
  }

  /** Number of cells whose problem blocks {@link save}. */
  get errorCount(): number {
    let n = 0;
    for (let row = 0; row < this._df.rowCount; row++) {
      const errors = this.errorsOf(row);
      for (const column of Object.keys(errors))
        if (errors[column].kind === 'error')
          n++;
    }
    return n;
  }

  // ─────────────────────── saving ─────────────────────────

  /** The pending batch as transaction ops, in row order: `'new'` rows insert
   * their writable values, `'modified'` rows update ONLY their changed columns
   * with the row's `expectedVersion`, `'deleted'` rows delete. Exposed so a
   * caller can inspect or extend the payload (a master-detail save appends its
   * own ops to one transaction).
   *
   * An empty cell of a NEW row is LEFT OUT of the insert rather than sent as an
   * explicit null, so the column takes its server-side default; a column with no
   * default and no value is rejected by the server's own nullability check
   * (and by {@link validate} before that). Clearing a cell of a MODIFIED row does
   * send null — that is an edit, not an omission. */
  buildOps(): DomainPendingOp[] {
    const table = this.client.table;
    const writable = this.capabilities.writableColumns;
    const pending: DomainPendingOp[] = [];
    for (let row = 0; row < this._df.rowCount; row++) {
      const state = this.stateOf(row);
      const id = this._wire(row, 'id');
      if (state === 'deleted') {
        // A row that was added and then deleted never reached the server.
        if (id != null)
          pending.push({row: row, op: {op: 'delete', table: table, id: `${id}`}});
      }
      else if (state === 'new') {
        const values: {[column: string]: any} = {};
        for (const name of writable) {
          const v = this._wire(row, name);
          if (v != null)
            values[name] = v;
        }
        pending.push({row: row, op: {op: 'insert', table: table, values: values}});
      }
      else if (state === 'modified') {
        const values: {[column: string]: any} = {};
        for (const name of Object.keys(this.changesOf(row)))
          if (writable.includes(name))
            values[name] = this._wire(row, name);
        if (Object.keys(values).length === 0)
          continue;
        const op: DG.DomainTransactionOp = {op: 'update', table: table, id: `${id}`, values: values};
        const version = this._wire(row, 'version');
        if (typeof version === 'number')
          op.expectedVersion = version;
        pending.push({row: row, op: op});
      }
    }
    return pending;
  }

  /**
   * Writes the whole pending batch as ONE `/transaction`: audit rows share a
   * `tx_id`, and any failure rolls every op back. Resolves to whether the batch
   * landed.
   *
   * Blocking cell errors refuse the save (naming the first one). A version
   * conflict goes through the platform's standard reload/overwrite dialog and
   * the chosen outcome is applied and retried. Server validation errors land on
   * the offending cells as {@link ERRORS_COLUMN} entries and everything stays
   * pending.
   *
   * **The editor is CLOSED while this runs** ({@link isSaving}): every write,
   * {@link discard} and {@link refresh} is refused with a warning instead of
   * being silently lost between the request and its results (and, for the
   * row-removing ones, instead of shifting the rows the results address). A grid
   * bound to the editor locks its own editing off {@link onSavingChanged}.
   */
  async save(): Promise<boolean> {
    if (this._busy('saving'))
      return false;
    const blocking = this._firstBlockingError();
    if (blocking != null) {
      grok.shell.error(`Cannot save: ${blocking}`);
      return false;
    }
    const pending = this.buildOps();
    if (pending.length === 0)
      return true;
    this._setSaving(true);
    this._saveDepth = 0;
    try {
      return await this._runSave(pending);
    } finally {
      this._setSaving(false);
    }
  }

  /**
   * Rebuilds the frame from the server: re-runs [query] (the attached one by
   * default) and re-attaches the service columns from scratch.
   *
   * **Pending edits do NOT survive this — by design.** There is no merge: the
   * control is refresh-agnostic and this method never reconciles old and new
   * state. Deciding WHETHER to refresh while edits are pending is the caller's
   * responsibility — check {@link isDirty} (or subscribe to
   * {@link onDirtyChanged}) and prompt the user to save or discard first.
   *
   * Resolves to the NEW frame, which also arrives on {@link onRefreshed}: a grid
   * bound to the old one must rebind. Refused (resolving to the CURRENT frame,
   * with a warning) while a {@link save} is in flight — the resolved value alone
   * does not tell a refusal from a rebuild, so check {@link isSaving} first (or
   * compare the frame identity) when it matters.
   */
  async refresh(query?: DG.DomainQuerySpec): Promise<DG.DataFrame> {
    if (this._busy('refreshing'))
      return this._df;
    const spec = query ?? this._query;
    const df = await this.client.queryDf(spec ?? {});
    this._query = spec;
    this._unsubscribe();
    this._df = df;
    this._bind(df);
    this._onRefreshed.next(df);
    this._fire();
    return df;
  }

  /** Releases the frame subscriptions. The service columns stay on the frame —
   * drop the frame, or remove them, if it outlives the editor. */
  detach(): void {
    this._unsubscribe();
    this._snapshots.clear();
  }

  // ─────────────────────── internals ─────────────────────────

  private _bind(df: DG.DataFrame): void {
    for (const name of SERVICE_COLUMNS) {
      let col = df.columns.byName(name);
      if (col == null)
        col = df.columns.addNewString(name).init('');
      // Both tags, both writers: the d42 serializer skips one, the CSV writer
      // the other — the state must reach neither.
      col.meta.includeInBinaryExport = false;
      col.meta.includeInCsvExport = false;
    }
    // NB: `_dirty` is deliberately NOT reset here — `_fire()` is what detects the
    // transition, and writing it directly would swallow the true → false edge a
    // refresh() produces (a subscriber that saw `true` would stay stale forever).
    this._resetCaches();
    this._subs.push(df.onRowsFiltering.subscribe(() => this._maskDeleted()));
    // The per-row caches are keyed by ROW INDEX, so anything that adds or removes
    // rows invalidates every index below it. The editor's own removals pass
    // `notify: false` and reset the caches themselves; these two cover the frame's
    // OTHER consumers (a grid's Remove-rows command, a co-owner of the frame).
    this._subs.push(df.onRowsRemoved.subscribe(() => this._resetCaches()));
    this._subs.push(df.onRowsAdded.subscribe(() => this._resetCaches()));
    df.rows.requestFilter();
  }

  private _unsubscribe(): void {
    for (const sub of this._subs)
      sub.unsubscribe();
    this._subs = [];
  }

  private _col(name: string): DG.Column<string> {
    return this._df.columns.byName(name) as DG.Column<string>;
  }

  private _cache(name: string): Map<number, any> {
    let byRow = this._parsed.get(name);
    if (byRow == null)
      this._parsed.set(name, byRow = new Map<number, any>());
    return byRow;
  }

  private _json(name: string, row: number): any {
    const byRow = this._cache(name);
    let value = byRow.get(row);
    if (value !== undefined)
      return value;
    const raw = this._col(name).get(row);
    value = {};
    if (raw != null && raw !== '') {
      try {
        value = JSON.parse(raw);
      } catch (_) { /* corrupt state reads as none */ }
    }
    byRow.set(row, value);
    return value;
  }

  private _setJson(name: string, row: number, value: any): void {
    const empty = value == null || Object.keys(value).length === 0;
    this._col(name).set(row, empty ? '' : JSON.stringify(value), false);
    this._cache(name).set(row, value ?? {});
    this._recount(row);
  }

  /** Drops every per-row cache and re-derives {@link changeCount} — for the
   * paths that add or remove rows, where every index below the change shifts.
   *
   * The pre-edit SNAPSHOTS are keyed by row index too: an external `removeAt`
   * leaves a stale key behind, on which `beginEdit` early-returns and the next
   * `commitEdit` records ANOTHER row's original value. They are dropped here,
   * with everything else that indices invalidate. */
  private _resetCaches(): void {
    this._parsed.clear();
    this._contributions.clear();
    this._snapshots.clear();
    this._pristine.clear();
    this._changeCount = 0;
    for (let row = 0; row < this._df.rowCount; row++)
      this._recount(row);
  }

  /** [row]'s share of {@link changeCount}: a new or deleted row counts once (a
   * PRISTINE new row not at all — see {@link addRow}), an edited one counts its
   * changed cells. */
  private _recount(row: number): void {
    const state = this.stateOf(row);
    const now = state === 'new' ? (this._pristine.has(row) ? 0 : 1)
      : state === 'deleted' ? 1 : Object.keys(this.changesOf(row)).length;
    this._changeCount += now - (this._contributions.get(row) ?? 0);
    if (now === 0)
      this._contributions.delete(row);
    else
      this._contributions.set(row, now);
  }

  /** Refuses a mutation while a save is in flight, saying so out loud: a write
   * landing between the transaction and its results would be wiped by
   * {@link _applyResults}, and a discard/refresh would shift the very rows those
   * results address. */
  private _busy(action: string): boolean {
    if (!this._saving)
      return false;
    grok.shell.warning(`${this.table}: ${action} is not available while the batch is being saved`);
    return true;
  }

  private _setSaving(saving: boolean): void {
    this._saving = saving;
    this._onSavingChanged.next(saving);
  }

  /** Runs [action] with the tracking latch closed, so writes the editor makes
   * itself (reverts, version bumps, reloaded server values) are not tracked as
   * user edits. */
  private _write(action: () => void): void {
    const was = this._suspend;
    this._suspend = true;
    try {
      action();
    } finally {
      this._suspend = was;
    }
  }

  private _wire(row: number, column: string): any {
    const col = this._df.columns.byName(column);
    if (col == null || col.isNone(row))
      return null;
    return toWire(col.get(row));
  }

  private _track(row: number, column: string, original: any): void {
    // The first write into a pristine new row is what makes it a pending change.
    this._pristine.delete(row);
    const invalid = this._invalidValue(row, column);
    if (this.stateOf(row) !== 'new') {
      const changes = this.changesOf(row);
      if (!(column in changes)) {
        if (original !== undefined)
          changes[column] = original;
      }
      // Back to the original AND valid: the cell is clean again (an invalid cell
      // stays pending so its marker survives) — the platform's own semantics.
      if (column in changes && !_isUnknown(changes[column]) &&
          invalid == null && wireEquals(changes[column], this._wire(row, column)))
        delete changes[column];
      this._setJson(CHANGES_COLUMN, row, changes);
    }
    // Read after the bookkeeping: whether the value would be DROPPED depends on
    // the change entry this call has just recorded (or removed).
    const message = this._droppedValue(row, column) ?? invalid;
    this._setError(row, column, message == null ? null : {message: message, kind: 'error'});
    this._recomputeState(row);
    this._fire();
  }

  /** Why the cell cannot be saved as it stands: column security first, then the
   * registry constraints. */
  private _cellProblem(row: number, column: string): string | null {
    return this._droppedValue(row, column) ?? this._invalidValue(row, column);
  }

  /** The cell against its registry {@link DG.Property} alone. An EMPTY cell is
   * empty whatever its column type spells that as — a string column's `''`, an int
   * column's null sentinel — so nullability, and not the sentinel, decides
   * (otherwise an untouched optional choice column reports `""` is not one of ...`,
   * and an untouched optional int reports itself below its minimum). */
  private _invalidValue(row: number, column: string): string | null {
    const property = this._propByName.get(column);
    if (property == null)
      return null;
    const col = this._df.columns.byName(column);
    return validateCellValue(property,
      col == null || col.isNone(row) ? null : this._df.get(column, row));
  }

  /**
   * Whether the cell holds a pending value {@link buildOps} would DROP.
   *
   * It sends writable columns ONLY — for a modified row it drops the change, for
   * a NEW row it drops the value (a prefilled parent FK on a column the caller
   * cannot write included, which would insert a child row with a null FK).
   * Either way the value would vanish without a word, so the cell is marked
   * instead and the save refuses loudly. Untouched cells are not a problem:
   * nothing of theirs would be dropped.
   */
  private _droppedValue(row: number, column: string): string | null {
    if (this.capabilities.writableColumns.includes(column))
      return null;
    const pending = this.stateOf(row) === 'new'
      ? this._wire(row, column) != null : column in this.changesOf(row);
    return pending ? `Column '${column}' is read-only` : null;
  }

  private _setError(row: number, column: string, error: DomainCellError | null): void {
    const errors = this.errorsOf(row);
    if (error == null)
      delete errors[column];
    else
      errors[column] = error;
    this._setJson(ERRORS_COLUMN, row, errors);
  }

  /** Validates every column the FRAME carries — the declared ones against the
   * registry, and the undeclared ones (a column the query returned that the
   * registry does not describe) against column security alone: a prefilled value
   * on one of those would be dropped by {@link buildOps} just as silently. */
  private _validateRow(row: number): void {
    for (const column of this._df.columns.names()) {
      if (SERVICE_COLUMNS.includes(column))
        continue;
      // The same predicate the per-cell path uses, so a prefilled value on a
      // non-writable column is marked here too (and stays marked: a re-validation
      // never clears a problem the cell still has).
      const message = this._cellProblem(row, column);
      const existing = this.errorOf(row, column);
      if (message != null)
        this._setError(row, column, {message: message, kind: 'error'});
      else if (existing != null && existing.kind === 'error')
        this._setError(row, column, null);
    }
  }

  private _recomputeState(row: number): void {
    const state = this.stateOf(row);
    if (state === 'new' || state === 'deleted')
      return;
    const changed = Object.keys(this.changesOf(row)).length > 0;
    this._col(STATE_COLUMN).set(row, changed ? 'modified' : '', false);
    this._recount(row);
  }

  private _clearRowState(row: number): void {
    this._col(STATE_COLUMN).set(row, '', false);
    this._col(CHANGES_COLUMN).set(row, '', false);
    this._col(ERRORS_COLUMN).set(row, '', false);
    this._cache(CHANGES_COLUMN).delete(row);
    this._cache(ERRORS_COLUMN).delete(row);
    this._recount(row);
  }

  private _setDeleted(rows: number | number[], deleted: boolean): void {
    if (this._busy(deleted ? 'deleting' : 'restoring a row'))
      return;
    const list = Array.isArray(rows) ? rows : [rows];
    const state = this._col(STATE_COLUMN);
    for (const row of list) {
      if (row < 0 || row >= this._df.rowCount)
        continue;
      if (deleted)
        state.set(row, 'deleted', false);
      else if (state.get(row) === 'deleted')
        // A row ADDED in this batch has no id — buildOps' own predicate — and its
        // values live nowhere but the frame, so it must go back to 'new'; a 'new'
        // row carries no ~changes by design, which is why the change test below
        // would otherwise land it on '' and make it invisible to save AND discard.
        state.set(row, this._wire(row, 'id') == null ? 'new'
          : Object.keys(this.changesOf(row)).length > 0 ? 'modified' : '', false);
      this._recount(row);
    }
    this._df.rows.requestFilter();
    this._fire();
  }

  /** Cooperative filtering: every participant ANDs its own exclusions into the
   * frame filter while it is being recomputed, so deleted rows stay hidden no
   * matter which other filter runs. */
  private _maskDeleted(): void {
    const state = this._df.columns.byName(STATE_COLUMN);
    if (state == null)
      return;
    for (let row = 0; row < this._df.rowCount; row++)
      if (state.get(row) === 'deleted')
        this._df.filter.set(row, false, false);
  }

  private _firstBlockingError(): string | null {
    for (let row = 0; row < this._df.rowCount; row++) {
      if (this.stateOf(row) === 'deleted')
        continue;
      const errors = this.errorsOf(row);
      for (const column of Object.keys(errors))
        if (errors[column].kind === 'error')
          return errors[column].message;
    }
    return null;
  }

  private _fire(): void {
    this._onChanged.next(this);
    const dirty = this.changeCount > 0;
    if (dirty !== this._dirty) {
      this._dirty = dirty;
      this._onDirtyChanged.next(dirty);
    }
  }

  private async _runSave(pending: DomainPendingOp[]): Promise<boolean> {
    if (pending.length === 0)
      return true;
    // Each conflict outcome resolves exactly one row, so the chain is finite;
    // the cap only guards against a server that keeps reporting the same op.
    if (this._saveDepth++ > 32) {
      grok.shell.error('Cannot save: too many version conflicts in a row');
      return false;
    }
    let results: any[];
    try {
      results = await grok.dapi.domains.transaction(this.client.schema, pending.map((p) => p.op));
    } catch (e: any) {
      return await this._onTransactionError(e, pending);
    }
    this._applyResults(pending, results);
    return true;
  }

  private _applyResults(pending: DomainPendingOp[], results: any[]): void {
    const removed: number[] = [];
    const result: DomainSaveResult = {inserted: 0, updated: 0, deleted: 0};
    this._write(() => {
      for (let i = 0; i < pending.length; i++) {
        const {op, row} = pending[i];
        const r = results[i] ?? {};
        if (op.op === 'delete') {
          removed.push(row);
          result.deleted++;
          continue;
        }
        if (op.op === 'insert') {
          result.inserted++;
          if (r.id != null && this._df.columns.contains('id'))
            this._df.set('id', row, `${r.id}`);
        }
        else
          result.updated++;
        if (r.version != null && this._df.columns.contains('version'))
          this._df.set('version', row, r.version);
        this._clearRowState(row);
      }
      removed.sort((a, b) => b - a);
      for (const row of removed)
        this._df.rows.removeAt(row, 1, false);
    });
    this._resetCaches();
    this._df.rows.requestFilter();
    this._fire();
    this._onSaved.next(result);
    const n = result.inserted + result.updated + result.deleted;
    grok.shell.info(`Saved ${n} row${n === 1 ? '' : 's'}`);
  }

  private async _onTransactionError(e: any, pending: DomainPendingOp[]): Promise<boolean> {
    const index = e?.opIndex;
    const failing = typeof index === 'number' && index >= 0 && index < pending.length
      ? pending[index] : null;
    if (e instanceof DG.DomainVersionConflictError && failing != null)
      return await this._resolveConflict(e, failing);
    if (e instanceof DG.DomainValidationError && failing != null) {
      this._mapValidationError(e, failing);
      grok.shell.error(e.message);
      return false;
    }
    grok.shell.error(e?.message ?? `${e}`);
    return false;
  }

  /** The platform's standard reload/overwrite dialog, with the outcome applied:
   * RELOAD takes the server's values for that ONE row (dropping its edits) and
   * retries the rest; OVERWRITE retries it against the current version. */
  private async _resolveConflict(e: DG.DomainVersionConflictError,
    failing: DomainPendingOp): Promise<boolean> {
    this._onConflict.next(e);
    const row = failing.row;
    const id = `${failing.op.id ?? this._wire(row, 'id')}`;
    const subject = `${this._displayOf(row) ?? id}`;
    const decision = await DG.DomainObjectHandler.showConflictDialog(subject);
    if (decision === 'reload') {
      let fresh: any = null;
      try {
        fresh = await this.client.get(id);
      } catch (_) { /* gone or invisible — reported below */ }
      if (fresh == null)
        grok.shell.error(`${subject} no longer exists.`);
      else
        this._write(() => {
          for (const name of Object.keys(fresh))
            if (this._df.columns.contains(name))
              this._df.set(name, row, fresh[name]);
        });
      this._clearRowState(row);
      // The pre-reload values are gone: a snapshot of them would make the next
      // in-grid edit record an "original" the cell never held.
      this._snapshots.delete(row);
      this._fire();
      return await this._runSave(this.buildOps());
    }
    if (decision === 'overwrite') {
      if (e.currentVersion != null)
        this._write(() => this._df.set('version', row, e.currentVersion));
      return await this._runSave(this.buildOps());
    }
    // Dismissed: nothing is written, and the row's changed cells say why.
    for (const column of Object.keys(this.changesOf(row)))
      this._setError(row, column, {message: e.message, kind: 'conflict'});
    this._fire();
    return false;
  }

  /** `rows[0].errors[{column, message}]` of a rejected op onto that row's cells;
   * an error naming no known column marks every changed cell of the row. */
  private _mapValidationError(e: DG.DomainValidationError, failing: DomainPendingOp): void {
    const row = failing.row;
    const errors = e.rows?.[0]?.errors ?? [];
    let mapped = false;
    for (const columnError of errors)
      if (columnError.column != null && this._df.columns.contains(columnError.column)) {
        this._setError(row, columnError.column, {message: columnError.message, kind: 'error'});
        mapped = true;
      }
    if (!mapped)
      for (const column of Object.keys(this.changesOf(row)))
        this._setError(row, column, {message: e.message, kind: 'error'});
    this._fire();
  }

  /** The row's display value for a dialog caption: the registry's declared name
   * column when the frame carries it — the same identity every other platform
   * surface shows. Null falls the caller back to the id. */
  private _displayOf(row: number): string | null {
    if (this._nameColumn == null || !this._df.columns.contains(this._nameColumn))
      return null;
    const v = this._wire(row, this._nameColumn);
    return v == null || `${v}` === '' ? null : `${v}`;
  }
}
