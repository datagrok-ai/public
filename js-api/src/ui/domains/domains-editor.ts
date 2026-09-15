/**
 * {@link DomainFrameEditor} — batch editing of domain-table rows on a DataFrame,
 * with the editing state living in the frame itself (the `~state` / `~changes` /
 * `~errors` service columns) rather than in a side store — plus the save/discard
 * prompt ({@link confirmDiscardChanges}) that stands between a user gesture and
 * an action which would throw a pending batch away.
 *
 * The editors themselves are refresh-agnostic by design — `refresh()` rebuilds
 * and discards, and deciding whether that is acceptable is the CALLER's job — so
 * the prompt is that decision implemented ONCE: `DomainAppView`,
 * `EntityListWidget` and any app driving an editor of its own route every
 * rebuild through {@link confirmDiscardChanges}.
 *
 * Several editors save as ONE transaction through a {@link DomainSession}
 * (`src/ui/domains/domains-session.ts`); a `DG.Grid` hosts an editor through
 * `Grid.attachEditor` (`src/grid.ts`); the handler, the built-in view and the
 * openers live in `src/domains-ui.ts`.
 *
 * @module domains-editor
 */

import * as rxjs from 'rxjs';

import * as ui from '../../../ui';
import {IDartApi} from '../../api/grok_api.g';
import {COLUMN_TYPE, ColumnType, TYPE} from '../../const';
import {DomainRegistryClient, DomainTableClient} from '../../dapi';
import {Column, DataFrame} from '../../dataframe';
import {DomainAccess, DomainConditionTree, DomainError, DomainQuerySpec, DomainRestrictError,
  DomainTableInfo, DomainTransactionOp, DomainValidationError, DomainVersionConflictError,
  DOMAIN_ACCESS_COLUMNS} from '../../domains';
import {DomainObjectHandler} from '../../domains-ui';
import {Property} from '../../entities/property';
import {IFrameEditor} from '../../grid';
import {Logger} from '../../logger';
import {Utils} from '../../utils';
import {Dialog} from '../../widgets/forms';
import {Balloon} from '../../widgets/menu';
import {DomainSession} from './domains-session';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** The `grok.shell` / `grok.log` / `grok.dapi.domains` surface the code below
 * uses, through the internal classes: no js-api module imports `../grok`. */
const balloon = new Balloon();
const log = Logger.getStatic();
const registry = new DomainRegistryClient();

/** Per-row editing state, as stored in {@link DomainFrameEditor.STATE_COLUMN}. */
export type DomainRowState = '' | 'new' | 'modified' | 'deleted';

/** One cell-level problem, as stored in {@link DomainFrameEditor.ERRORS_COLUMN}: `'error'` blocks
 * {@link DomainFrameEditor.save}, `'conflict'` (a dismissed version conflict)
 * only marks the cell. */
export interface DomainCellError {
  message: string;
  kind: 'error' | 'conflict';
}

/** An edit the host refused because the cell may not be edited
 * ({@link DomainFrameEditor.onRefused}). */
export interface DomainEditRefusal {
  row: number;
  column: string;
  message: string;
}

/** What one {@link DomainFrameEditor.save} (or a {@link DomainSession.save})
 * wrote — all of it in ONE transaction. Ids and versions the server returned are
 * written back into the frame; a host that keeps its own row objects re-queries
 * after a successful save. */
export interface DomainSaveResult {
  inserted: number;
  updated: number;
  deleted: number;
  /** Server id of every inserted row, keyed by the draft id
   * ({@link DomainFrameEditor.draftId}) the row carried before the save. */
  assigned: {[draftId: string]: string};
}

/** Options of {@link DomainFrameEditor.attach} / {@link DomainFrameEditor.create}. */
export interface DomainFrameEditorOptions {
  /** The query the frame came from — what {@link DomainFrameEditor.refresh}
   * re-runs, and what {@link DomainFrameEditor.create} runs to build it. */
  query?: DomainQuerySpec;
  /** Pre-probed access, to avoid a second round trip when the caller
   * already has it (a prefetched {@link IDomainTableContext} passes its own). */
  access?: DomainAccess;
  /** Suppresses the editor's own informational balloons ('Saved N rows') for a
   * host that phrases its own feedback; errors and the conflict dialog always show. */
  quiet?: boolean;
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

/** The draft-id map as a dictionary with NO prototype, taking the OWN keys of every source:
 * a cell (or a caller's object) holding `constructor` / `toString` must read as absent here,
 * not as an inherited Function that would then be written into a frame. */
function toDraftMap(...sources: ({[draftId: string]: string} | undefined)[]): {[draftId: string]: string} {
  const map: {[draftId: string]: string} = Object.create(null);
  for (const source of sources)
    for (const key of Object.keys(source ?? {}))
      map[key] = source![key];
  return map;
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
 * {@link DomainFrameEditor.validateCellValue}. */
function valueTypeOf(v: any): string {
  if (typeof v === 'string')
    return TYPE.STRING;
  if (typeof v === 'boolean')
    return TYPE.BOOL;
  if (typeof v === 'number')
    return Number.isInteger(v) ? TYPE.INT : TYPE.FLOAT;
  if (Array.isArray(v))
    return TYPE.LIST;
  if (v != null && typeof v === 'object' && typeof v.toISOString === 'function')
    return TYPE.DATE_TIME;
  return TYPE.MAP;
}

/** Numbers formatted the way the server VM prints them, so the messages below
 * match the server's byte-for-byte (the VM prints whole doubles as '1.0'). */
function vmNum(v: number): string {
  return v % 1 === 0 ? v.toFixed(1) : `${v}`;
}

function isNumerical(p: Property): boolean {
  const t = p.propertyType;
  return t === TYPE.INT || t === TYPE.FLOAT || t === TYPE.NUM || t === TYPE.BIG_INT;
}

function validateNumeric(p: Property, v: number): string | null {
  if (p.propertyType === TYPE.INT && v % 1 !== 0)
    return 'Types differ. Expected: int, passed: double';
  const passed = p.propertyType === TYPE.INT ? `${v}` : vmNum(v);
  if (p.min != null && v < p.min)
    return `Value should not be less than ${vmNum(p.min)}, passed: ${passed}`;
  if (p.max != null && v > p.max)
    return `Value should not be more than ${vmNum(p.max)}, passed: ${passed}`;
  return null;
}

/**
 * A widget or page that answers for a set of {@link DomainFrameEditor}s.
 *
 * The dirty rollup is domain-specific — it does not belong on `Widget`, and it
 * cannot be derived from `Widget.children` (which is a DOM view, so it is order-
 * and reparenting-fragile). This is the one structural remainder of the widget
 * model: anything that owns pending changes says so by exposing them here, and a
 * container (`DomainAppView`, `domains.view`) rolls its children's up.
 */
export interface IEditorHost {
  /** Every editor whose pending changes this object answers for. Read on demand:
   * an editor that appears after the widget was constructed must show up here. */
  readonly editors: DomainFrameEditor[];
}

/** A domain client with its row/insert/column/expand/update generics erased — what
 * the REFLECTIVE widgets take. A TYPED client fits it, so a generated `<Table>Ui`
 * (or a `DomainTable` handle) keeps its own generics and is still usable here.
 * The update generic must be erased explicitly: its `Partial<TRow>` default is a
 * real object type, and `update()` would then reject a typed client's payload. */
export type AnyDomainTableClient = DomainTableClient<any, any, any, any, any>;

/**
 * The prefetched table context a SYNCHRONOUS widget factory needs: the typed
 * client plus the registry metadata and access, resolved once by
 * `domains.table(...)` (which is why `form()`, `grid()` and friends need no await).
 *
 * Access is a SNAPSHOT taken when the context was acquired — the contract
 * the grid has always had, moved one level up: a later grant change (or a
 * `grok.dapi.domains.invalidateUiCaches()`) does not reach widgets already built
 * from this context; re-acquire the handle to re-gate.
 */
export interface IDomainTableContext {
  readonly client: AnyDomainTableClient;
  /** Registry {@link Property} metadata of the table's declared columns. */
  readonly properties: Property[];
  readonly info: DomainTableInfo;
  readonly access: DomainAccess;
  /** `'<schema>.<table>'`. */
  readonly table: string;
}

/**
 * Resolves [client]'s registry metadata and the caller's access in ONE
 * round of requests — the async boundary every synchronous widget factory sits
 * behind (`domains.table()` and the u2 domain controls).
 */
export async function acquireDomainContext(
  client: AnyDomainTableClient): Promise<IDomainTableContext> {
  const address = `${client.schema}.${client.table}`;
  const [properties, info, access] = await Promise.all([
    registry.rowProperties(address),
    registry.tableInfo(address),
    client.access(),
  ]);
  return {client: client, properties: properties, info: info, access: access, table: address};
}

/** An empty frame of [properties]' columns — what the local (no round trip)
 * editors of {@link DomainFrameEditor.forContext} / `forRows` start from. */
function _emptyFrame(properties: Property[]): DataFrame {
  const df = DataFrame.create(0);
  for (const p of properties)
    df.columns.add(Column.fromType(columnTypeOf(p), p.name, 0));
  return df;
}

/** The frame column type that holds [p]'s values — what {@link DomainFrameEditor.forRows}
 * builds a local frame from (a `queryDf` result carries the server's own types). */
function columnTypeOf(p: Property): ColumnType {
  switch (p.propertyType) {
    case TYPE.INT: return COLUMN_TYPE.INT;
    case TYPE.BIG_INT: return COLUMN_TYPE.BIG_INT;
    case TYPE.FLOAT:
    case TYPE.NUM: return COLUMN_TYPE.FLOAT;
    case TYPE.BOOL: return COLUMN_TYPE.BOOL;
    case TYPE.DATE_TIME: return COLUMN_TYPE.DATE_TIME;
    default: return COLUMN_TYPE.STRING;
  }
}

/** One pending write of {@link DomainFrameEditor.buildOps}, with the frame row it
 * came from (how a failing `opIndex` finds its row). */
export interface DomainPendingOp {
  op: DomainTransactionOp;
  row: number;
}

/**
 * THE single writer of a domain frame's editing state.
 *
 * It wraps a DataFrame produced by `table.queryDf(...)` and attaches three
 * invisible service columns — {@link DomainFrameEditor.STATE_COLUMN},
 * {@link DomainFrameEditor.CHANGES_COLUMN}, {@link DomainFrameEditor.ERRORS_COLUMN} —
 * that hold everything about the pending batch: which
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
export class DomainFrameEditor implements IFrameEditor {
  /** Row state column: `'' | 'new' | 'modified' | 'deleted'`. */
  static readonly STATE_COLUMN = '~state';
  /** JSON column holding the ORIGINAL values of changed cells only (sparse). */
  static readonly CHANGES_COLUMN = '~changes';
  /** JSON column holding per-cell {@link DomainCellError}s. */
  static readonly ERRORS_COLUMN = '~errors';

  /** The three service columns an editor attaches — every one of them tagged out
   * of binary AND csv export, so the editing state can never reach a saved
   * project, an export, an upload, or a `batch()` fed from the frame. */
  static readonly SERVICE_COLUMNS: readonly string[] =
    [DomainFrameEditor.STATE_COLUMN, DomainFrameEditor.CHANGES_COLUMN, DomainFrameEditor.ERRORS_COLUMN];

  /** Prefix of the id {@link addRow} stamps into a row that does not exist on the
   * server yet: `~new:<uuid>`. Another row (of this or of another editor in the
   * same {@link DomainSession}) may hold it in a ref column — {@link buildOps}
   * turns it into the transaction's `$ref` and the server resolves it. */
  static readonly DRAFT_ID_PREFIX = '~new:';

  /** The referential refusal the server sends back, which names the child table and the column
   * pointing here (`repository.dart` ~:2212) — {@link restrictRefusal} says it in the user's words. */
  static readonly LIVE_ROWS = /Live rows in "([^"]+)" reference this row via "([^"]+)"/;

  static isDraftId(v: unknown): boolean {
    return typeof v === 'string' && v.startsWith(DomainFrameEditor.DRAFT_ID_PREFIX);
  }

  static draftId(): string {
    return DomainFrameEditor.DRAFT_ID_PREFIX + Utils.uuid4();
  }

  /** The draft id of every insert of [pending] → the id the server gave it. */
  static assignedOf(pending: DomainPendingOp[], results: any[]): {[draftId: string]: string} {
    const assigned = toDraftMap();
    for (let i = 0; i < pending.length; i++) {
      const {op} = pending[i];
      const id = results[i]?.id;
      if (op.op === 'insert' && op.ref != null && id != null)
        assigned[op.ref] = `${id}`;
    }
    return assigned;
  }

  private _df: DataFrame;
  private _subs: rxjs.Subscription[] = [];
  private _properties: Property[];
  private _propByName = new Map<string, Property>();
  private _snapshots = new Map<number, {[column: string]: any}>();
  /** Parsed {@link DomainFrameEditor.CHANGES_COLUMN} / {@link DomainFrameEditor.ERRORS_COLUMN} objects, per column
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
  private _dirty = false;
  private _query?: DomainQuerySpec;
  /** See {@link DomainFrameEditorOptions.quiet}. */
  readonly quiet: boolean;
  private _info: DomainTableInfo;

  private readonly _onChanged = new rxjs.Subject<DomainFrameEditor>();
  private readonly _onDirtyChanged = new rxjs.Subject<boolean>();
  private readonly _onSavingChanged = new rxjs.Subject<boolean>();
  private readonly _onSaved = new rxjs.Subject<DomainSaveResult>();
  private readonly _onConflict = new rxjs.Subject<DomainVersionConflictError>();
  private readonly _onRefreshed = new rxjs.Subject<DataFrame>();
  private readonly _onRefused = new rxjs.Subject<DomainEditRefusal>();

  private constructor(
    /** The table the frame's rows belong to. */
    public readonly client: DomainTableClient,
    /** Effective access of the current user, SNAPSHOT when the editor was
     * created — what read-only degradation and the writable-column payload filter
     * derive from. A later `grok.dapi.domains.invalidateUiCaches()` (or a grant
     * change) does NOT reach an existing editor: re-create it to pick the new
     * permissions up. */
    public readonly access: DomainAccess,
    df: DataFrame, properties: Property[], info: DomainTableInfo,
    options?: DomainFrameEditorOptions) {
    this._properties = properties;
    for (const p of properties)
      this._propByName.set(p.name, p);
    this._info = info;
    this._query = options?.query;
    this.quiet = options?.quiet === true;
    this._df = df;
    this._bind(df);
  }

  /** The columns [access] lets the caller write, in the server's (declared) order. */
  static writableColumns(access: DomainAccess): string[] {
    return Object.keys(access.fields).filter((c) => access.fields[c] === 'editable');
  }

  /** {@link IFrameEditor}: null when NOTHING in the frame can be edited — no
   * table-level `edit` or `insert`, and no per-row `~can_edit` to override them —
   * otherwise the columns the field rights let anyone write, {@link canEdit}
   * deciding the ROW dimension per cell. A row-mode frame answers the list while
   * the table-level `edit` is false, or the whole grid would lock. */
  get writableColumns(): string[] | null {
    const columns = DomainFrameEditor.writableColumns(this.access);
    if (this.access.can.edit === true || this.access.can.insert === true)
      return columns;
    return this._df.columns.byName(DOMAIN_ACCESS_COLUMNS[0]) != null ? columns : null;
  }

  /** {@link IFrameEditor}: whether [column] of [row] may be edited — a writable
   * column, plus the row's own right: a draft needs `insert`; a persisted row
   * carries `~can_edit` in row mode (a `withAccess` read), and falls back to the
   * table-level `edit` where the frame has no such column. */
  canEdit(row: number, column: string): boolean {
    if (!DomainFrameEditor.writableColumns(this.access).includes(column))
      return false;
    // a row appended here has no server answer in `~can_edit` (a bool column has no null slot,
    // so the cell reads false): the right that governs it is `insert`
    if (this.stateOf(row) === 'new')
      return this.access.can.insert === true;
    const rowRight = this._df.columns.byName(DOMAIN_ACCESS_COLUMNS[0]);
    return rowRight != null ? rowRight.get(row) === true : this.access.can.edit === true;
  }

  /** {@link IFrameEditor}: why [column] of [row] may not be edited, in the words the table
   * itself uses — the column's caption where the field is read-only for everyone, the row's
   * own name where the row is; null when the cell may be edited. */
  refusalOf(row: number, column: string): string | null {
    if (this.canEdit(row, column))
      return null;
    if (!DomainFrameEditor.writableColumns(this.access).includes(column))
      return `${this._captionOf(column)} is read-only`;
    return `${this._displayOf(row) ?? `This ${this._info.singularName.toLowerCase()}`} is read-only for you`;
  }

  /** {@link IFrameEditor}: an edit the HOST refused ({@link refusalOf} said why) — reported on
   * {@link onRefused}, so a host with a status line says it there too and not only in the
   * grid's balloon. Nothing is written to the frame: a refused edit never happened. */
  refuse(row: number, column: string, message: string): void {
    this._onRefused.next({row: row, column: column, message: message});
  }

  /** Attaches the editing state to an EXISTING frame of [client]'s rows (a
   * `queryDf` result). Pass `options.query` so {@link refresh} knows what to
   * re-run. */
  static async attach(dataFrame: DataFrame, client: DomainTableClient,
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    const address = `${client.schema}.${client.table}`;
    const [properties, access, info] = await Promise.all([
      registry.rowProperties(address),
      options?.access != null ? Promise.resolve(options.access) : client.access(),
      registry.tableInfo(address),
    ]);
    return new DomainFrameEditor(client, access, dataFrame, properties, info, options);
  }

  /** {@link attach} for a host that has a frame and a table address but no client
   * of its own — the entry point of the Dart Domain View's grid mode, which owns
   * its frame and hands it over to be edited. */
  static attachTo(dataFrame: DataFrame, schema: string, table: string,
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    return DomainFrameEditor.attach(dataFrame,
      new DomainTableClient(api.grok_Dapi_Domains(), schema, table), options);
  }

  /** Runs `options.query` (everything the caller may see, by default) and
   * attaches to the resulting frame. */
  static async create(client: DomainTableClient,
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
      options?.access ?? context.access, _emptyFrame(context.properties),
      context.properties, context.info, options);
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
  static async forRows(client: DomainTableClient, values: {[column: string]: any}[],
    options?: DomainFrameEditorOptions): Promise<DomainFrameEditor> {
    const properties = await registry.rowProperties(
      `${client.schema}.${client.table}`);
    const editor = await DomainFrameEditor.attach(_emptyFrame(properties), client, options);
    for (const row of values ?? [])
      editor.addRow(row, {pristine: true});
    return editor;
  }

  /**
   * Validates one cell value against its registry {@link Property} — the same
   * constraints the server re-runs on write, producing the server's exact message
   * texts so an inline marker and a rejected save read identically.
   *
   * Numeric properties take a dart2js-safe path (under dart2js every whole number
   * `is double`, which would reject every integer), and `string_list` values are
   * left to the server, which coerces them. Returns null when the value is fine.
   */
  static validateCellValue(p: Property, value: any): string | null {
    if (!p.nullable && (value == null || value === ''))
      return "Value can't be empty";
    if (value == null || p.propertyType === TYPE.STRING_LIST)
      return null;
    const type = valueTypeOf(value);
    if (isNumerical(p))
      return typeof value === 'number' ? validateNumeric(p, value)
        : `Numerical value expected, passed: ${type}`;
    const expected = p.propertyType === TYPE.NUM ? TYPE.FLOAT : p.propertyType;
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
  static isReferenceProperty(p: Property): boolean {
    const semType = `${p.semType ?? ''}`;
    return semType === 'User' || semType === 'Group' || /^[^.]+\.[^.]+$/.test(semType);
  }

  /** The editors [widget] answers for, or none — the duck-typed read of
   * {@link IEditorHost} a container uses on children of any class. */
  static editorsOf(widget: any): DomainFrameEditor[] {
    const editors = widget?.editors;
    return Array.isArray(editors) ? editors : [];
  }

  /** The frame being edited. It is REPLACED by {@link refresh} — re-read it (or
   * subscribe to {@link onRefreshed}) instead of caching it. */
  get dataFrame(): DataFrame { return this._df; }

  /** `'<schema>.<table>'`. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** Registry {@link Property} metadata of the table's declared columns. */
  get properties(): Property[] { return this._properties; }

  /** The query {@link refresh} re-runs. */
  get query(): DomainQuerySpec | undefined { return this._query; }

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
  get onConflict(): rxjs.Observable<DomainVersionConflictError> { return this._onConflict; }
  /** Fires with the NEW frame after {@link refresh} rebuilt it. */
  get onRefreshed(): rxjs.Observable<DataFrame> { return this._onRefreshed; }
  /** Fires when the HOST refused an edit ({@link refuse}) — what a status line shows next to
   * the grid's balloon. */
  get onRefused(): rxjs.Observable<DomainEditRefusal> { return this._onRefused; }

  // ─────────────────────── state accessors ─────────────────────────

  /** Editing state of [row]. */
  stateOf(row: number): DomainRowState {
    return (this._col(DomainFrameEditor.STATE_COLUMN).get(row) ?? '') as DomainRowState;
  }

  /** ORIGINAL values of [row]'s changed cells, keyed by column (empty when the
   * row is unchanged; always empty for a `'new'` row — all of its values are
   * new). Read-only: the object is the editor's own cached parse, and writing to
   * it changes nothing on the frame. */
  changesOf(row: number): {[column: string]: any} {
    return this._json(DomainFrameEditor.CHANGES_COLUMN, row);
  }

  /** Per-cell problems of [row], keyed by column. Read-only, see
   * {@link changesOf}. */
  errorsOf(row: number): {[column: string]: DomainCellError} {
    return this._json(DomainFrameEditor.ERRORS_COLUMN, row);
  }

  /** Whether the cell carries a pending change (what highlighting keys on). */
  isChanged(row: number, column: string): boolean {
    return !column.startsWith('~') && (this.stateOf(row) === 'new' || column in this.changesOf(row));
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

  /** Drops the {@link beginEdit} snapshots. A HOST that rewrites row values IN
   * PLACE (a re-query refreshing the same frame object, with no row-count
   * change and no frame event) MUST call this afterwards: the snapshots are
   * keyed by row index and would otherwise record the OLD rows' values as edit
   * originals (STATE-CONTRACT H8). */
  clearSnapshots(): void {
    this._snapshots.clear();
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
    // Silent throughout: a frame event fired between the append and the stamp would
    // surface a row that is neither keyed nor `'new'` to every other consumer of the
    // frame (a grid, a session, a list rebuilding on each event).
    this._write(() => {
      this._df.rows.addNew(null, false);
      if (values != null)
        for (const name of Object.keys(values))
          if (this._df.columns.contains(name))
            this._col(name).set(row, values[name], false);
      if (this._df.columns.contains('id'))
        this._col('id').set(row, DomainFrameEditor.draftId(), false);
      this._col(DomainFrameEditor.STATE_COLUMN).set(row, 'new', false);
    });
    this._validateRow(row);
    // The row was appended while the filter was already computed — recompute it,
    // or the new row is invisible in every filtered view of the frame.
    this._df.rows.requestFilter();
    // The append's ONE notification, fired now that the row is complete: a zero-row
    // insert raises onRowsAdded (nothing else in the frame API does) without touching
    // the frame.
    this._df.rows.insertAt(row, 0, true);
    // After the notification: its cache reset would otherwise drop the mark this
    // line makes.
    if (options?.pristine === true)
      this._pristine.add(row);
    this._recount(row);
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
      log.warning(`${this.table}: cannot revert ${column} — its original value was not captured`);
      return;
    }
    this._write(() => this._df.set(column, row, original));
    delete changes[column];
    this._setJson(DomainFrameEditor.CHANGES_COLUMN, row, changes);
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
   * their writable values (naming their draft id as the op's `ref`), `'modified'`
   * rows update ONLY their changed columns with the row's `expectedVersion`,
   * `'deleted'` rows delete. Exposed so a caller can inspect or extend the
   * payload; a {@link DomainSession} concatenates several editors' into one
   * transaction.
   *
   * A value equal to a draft id — this editor's or another's — goes out as the
   * `'$<draftId>'` reference the server resolves; a literal leading `$` is
   * escaped as `$$`. An empty cell of a NEW row is LEFT OUT of the insert rather
   * than sent as an explicit null, so the column takes its server-side default; a
   * column with no default and no value is rejected by the server's own
   * nullability check (and by {@link validate} before that). Clearing a cell of a
   * MODIFIED row does send null — that is an edit, not an omission.
   *
   * Every insert carries `onDuplicate: 'error'`: a business-key conflict fails
   * the whole transaction (409) instead of silently merging the new row into an
   * existing one. */
  buildOps(): DomainPendingOp[] {
    const table = this.client.table;
    const writable = DomainFrameEditor.writableColumns(this.access);
    const pending: DomainPendingOp[] = [];
    for (let row = 0; row < this._df.rowCount; row++) {
      const state = this.stateOf(row);
      const id = this._wire(row, 'id');
      if (state === 'deleted') {
        // A row that was added and then deleted never reached the server.
        if (this._isPersisted(row))
          pending.push({row: row, op: {op: 'delete', table: table, id: `${id}`}});
      }
      else if (state === 'new') {
        const values: {[column: string]: any} = {};
        for (const name of writable) {
          const v = this._wire(row, name);
          if (v != null && name !== 'id')
            values[name] = DomainFrameEditor._refValue(v);
        }
        const op: DomainTransactionOp = {op: 'insert', table: table, values: values, onDuplicate: 'error'};
        if (DomainFrameEditor.isDraftId(id))
          op.ref = `${id}`;
        pending.push({row: row, op: op});
      }
      else if (state === 'modified') {
        const values: {[column: string]: any} = {};
        for (const name of Object.keys(this.changesOf(row)))
          if (writable.includes(name))
            values[name] = DomainFrameEditor._refValue(this._wire(row, name));
        if (Object.keys(values).length === 0)
          continue;
        const op: DomainTransactionOp = {op: 'update', table: table, id: `${id}`, values: values};
        const version = this._wire(row, 'version');
        if (typeof version === 'number')
          op.expectedVersion = version;
        pending.push({row: row, op: op});
      }
    }
    return pending;
  }

  /** [v] as the transaction wire wants it: a draft id becomes a `$ref`, a literal
   * leading `$` is doubled (the server's `$$` escape), lists element-wise. */
  private static _refValue(v: any): any {
    if (Array.isArray(v))
      return v.map((x) => DomainFrameEditor._refValue(x));
    return typeof v === 'string' && (DomainFrameEditor.isDraftId(v) || v.startsWith('$')) ? '$' + v : v;
  }

  /**
   * Writes the whole pending batch as ONE `/transaction`: audit rows share a
   * `tx_id`, and any failure rolls every op back. Resolves to whether the batch
   * landed. A session of one — see {@link DomainSession.save} for the flow
   * (blocking errors, the conflict dialog, validation errors on cells) and for
   * saving several editors together.
   *
   * **The editor is CLOSED while this runs** ({@link isSaving}): every write,
   * {@link discard} and {@link refresh} is refused with a warning instead of
   * being silently lost between the request and its results (and, for the
   * row-removing ones, instead of shifting the rows the results address). A grid
   * bound to the editor locks its own editing off {@link onSavingChanged}.
   */
  async save(): Promise<boolean> {
    const session = new DomainSession([this], {quiet: this.quiet});
    try {
      return await session.save();
    } finally {
      session.dispose();
    }
  }

  // ─────────────────────── the save participant protocol (DomainSession) ─────────────────────────

  /** The batch a {@link DomainSession} takes from this editor: null when the
   * editor refuses (a blocking cell error, named in a balloon), else the ops of
   * {@link buildOps} — after a row added and then deleted in the same batch has
   * been resolved LOCALLY (removed from the frame: it never reached the server, so
   * it can neither survive the save as phantom pending state nor hold the batch
   * dirty). */
  prepareSave(): DomainPendingOp[] | null {
    const blocking = this._firstBlockingError();
    if (blocking != null) {
      balloon.error(`Cannot save: ${blocking}`);
      return null;
    }
    let dropped = 0;
    this._write(() => {
      for (let row = this._df.rowCount - 1; row >= 0; row--)
        if (this.stateOf(row) === 'deleted' && !this._isPersisted(row)) {
          this._df.rows.removeAt(row, 1, false);
          dropped++;
        }
    });
    if (dropped > 0) {
      this._resetCaches();
      this._df.rows.requestFilter();
      this._fire();
    }
    return this.buildOps();
  }

  /** Lands this editor's slice of a successful transaction: returned ids and
   * versions into the frame, every cell holding a draft id the transaction
   * resolved (this editor's, or another participant's through [assigned] — the
   * {@link DomainSession} collects the whole batch's map before the first slice
   * is applied), then every OTHER server-assigned column of the rows it wrote
   * ({@link writeBack}), row state cleared, deleted rows removed.
   *
   * {@link onSaved} fires LAST, once the re-read has landed: a host that rebuilds
   * on it (a form's system footer, a list) must see the server's rows — the
   * transaction answers `{id, version, created}` per insert and nothing else, so
   * `number`, `created_on`, `updated_on` and `author_id` exist only after it. */
  async applyResults(pending: DomainPendingOp[], results: any[],
      assigned?: {[draftId: string]: string}): Promise<DomainSaveResult> {
    const removed: number[] = [];
    const result: DomainSaveResult = {inserted: 0, updated: 0, deleted: 0,
      assigned: DomainFrameEditor.assignedOf(pending, results)};
    const ids = toDraftMap(assigned, result.assigned);
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
          // The frame defaults (false) would read the just-created row as locked;
          // the table-level rights stand in until the post-save re-read lands the
          // server's per-row value (or for good, if that re-read fails). Share is
          // per row only in row mode — null elsewhere, as the server answers it.
          for (const column of this._df.columns.names())
            if (column.startsWith('~can_')) {
              const right = column.slice('~can_'.length);
              this._df.set(column, row, right === 'share' && this.access.securityMode !== 'row'
                ? null : this.access.can[right] === true);
            }
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
    this.rebind(ids);
    this._df.rows.requestFilter();
    this._fire();
    await this.writeBack(pending, results);
    this._onSaved.next(result);
    return result;
  }

  /**
   * Rewrites every cell holding a draft id the transaction resolved ([assigned] maps draft id
   * → server id) to the real id, WITHOUT touching the row's editing state: a pristine row
   * stays pristine, a clean row stays clean, and nothing is recorded as a change.
   *
   * {@link applyResults} does this for the editors that took part in the batch; a
   * {@link DomainSession} applies it to the ones that did NOT — a pristine child holding a
   * `~new:` reference to a parent the batch just created keeps a dangling draft id otherwise.
   *
   * Only a cell whose whole value IS a draft id is rewritten — never the `'$~new:…'` wire
   * form {@link buildOps} escapes to, and never a substring.
   */
  rebind(assigned: {[draftId: string]: string}): void {
    const ids = toDraftMap(assigned);
    if (Object.keys(ids).length === 0)
      return;
    let rebound = false;
    this._write(() => {
      for (const name of this._df.columns.names()) {
        if (name === 'id' || DomainFrameEditor.SERVICE_COLUMNS.includes(name)
            || this._df.columns.byName(name).type !== COLUMN_TYPE.STRING)
          continue;
        for (let row = 0; row < this._df.rowCount; row++) {
          const value = this._df.get(name, row);
          const real = DomainFrameEditor.isDraftId(value) ? ids[value] : undefined;
          if (real != null) {
            this._df.set(name, row, real);
            rebound = true;
          }
        }
      }
    });
    if (rebound)
      this._fire();
  }

  /** Re-reads the rows a save inserted or updated and lands EVERY returned column
   * in the frame — server defaults, the autoNumber `number`,
   * `created_on`/`author_id`/`updated_on`, the per-row `~can_*` — so a host
   * reading the frame after a save sees the server's row, not the payload it
   * sent. Run by {@link applyResults} before {@link onSaved}; ENRICHES only — it
   * fails soft, and the ids, versions and resolved draft references
   * {@link applyResults} wrote stand without it. */
  async writeBack(pending: DomainPendingOp[], results: any[]): Promise<void> {
    const ids = pending.map((p, i) => p.op.op === 'insert' ? results[i]?.id
      : p.op.op === 'update' ? p.op.id : null).filter((id) => id != null).map((id) => `${id}`);
    if (ids.length === 0)
      return;
    let rows: any[];
    try {
      rows = await this.client.query({filter: {property: 'id', operator: '=', value: ids},
        withAccess: this._df.columns.contains(DOMAIN_ACCESS_COLUMNS[0])});
    } catch (e) {
      log.warning(`${this.table}: post-save re-read failed — ${e}`);
      return;
    }
    const byId = new Map<string, any>(rows.map((r) => [`${r.id}`, r]));
    this._write(() => {
      for (let row = 0; row < this._df.rowCount; row++) {
        const fresh = byId.get(`${this._wire(row, 'id')}`);
        if (fresh == null)
          continue;
        for (const name of Object.keys(fresh)) {
          if (!this._df.columns.contains(name) || DomainFrameEditor.SERVICE_COLUMNS.includes(name))
            continue;
          const v = fresh[name];
          // jsonb cells travel as objects; the frame holds them as text. dayjs passes through.
          const plain = v != null && typeof v === 'object' && typeof v.toISOString !== 'function';
          this._df.set(name, row, plain ? JSON.stringify(v) : v);
        }
      }
    });
    this._resetCaches();
    this._fire();
  }

  /** The platform's standard reload/overwrite dialog for a version conflict on
   * [failing], with the outcome applied: RELOAD takes the server's values for
   * that ONE row (dropping its edits), OVERWRITE takes the current version.
   * Resolves to whether the batch should be rebuilt and retried; false when the
   * user dismissed it, in which case the row's changed cells say why. */
  async resolveConflict(e: DomainVersionConflictError, failing: DomainPendingOp): Promise<boolean> {
    this._onConflict.next(e);
    const row = failing.row;
    const id = `${failing.op.id ?? this._wire(row, 'id')}`;
    const subject = `${this._displayOf(row) ?? id}`;
    const decision = await DomainObjectHandler.showConflictDialog(subject);
    if (decision === 'reload') {
      let fresh: any = null;
      try {
        fresh = await this.client.get(id);
      } catch (_) { /* gone or invisible — reported below */ }
      if (fresh == null)
        balloon.error(`${subject} no longer exists.`);
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
      return true;
    }
    if (decision === 'overwrite') {
      if (e.currentVersion != null)
        this._write(() => this._df.set('version', row, e.currentVersion));
      return true;
    }
    for (const column of Object.keys(this.changesOf(row)))
      this._setError(row, column, {message: e.message, kind: 'conflict'});
    this._fire();
    return false;
  }

  /** `rows[0].errors[{column, message}]` of a rejected op onto that row's cells;
   * an error naming no known column marks every changed cell of the row. [duplicateText] is what
   * a business-key clash says — {@link refusalFor}'s sentence, so the cell, the status line
   * reading it and the balloon all say the same thing. */
  mapValidationError(e: DomainValidationError, failing: DomainPendingOp, duplicateText?: string): void {
    const row = failing.row;
    const errors = e.rows?.[0]?.errors ?? [];
    const duplicate = duplicateText ?? this._duplicateText();
    let mapped = false;
    for (const columnError of errors)
      if (columnError.column != null && this._df.columns.contains(columnError.column)) {
        const message = columnError.code === 'unique' ? duplicate : columnError.message;
        this._setError(row, columnError.column, {message: message, kind: 'error'});
        mapped = true;
      }
    if (!mapped) {
      const message = e.isDuplicate ? duplicate : e.message;
      for (const column of Object.keys(this.changesOf(row)))
        this._setError(row, column, {message: message, kind: 'error'});
    }
    this._fire();
  }

  /** What a refused save means to the user, in the table's own vocabulary — the business key and
   * the row that already holds it, the children that keep a row alive. The server's own message
   * stands where the refusal is not one this editor can put in the table's words. */
  async refusalFor(e: any, failing: DomainPendingOp): Promise<string> {
    if (e instanceof DomainValidationError && e.isDuplicate)
      return this.duplicateRefusal(failing.row, e);
    if (e instanceof DomainRestrictError)
      return this.restrictRefusal(e);
    return `${e?.message ?? e}`;
  }

  /** "CAS 64-17-5 already belongs to Ethanol": the key's caption and the value at hand, and the
   * row that already holds it — read back through the id the refusal carries, or looked up by
   * the key itself where it carries none: on the UPDATE path the id is the row being SAVED
   * (`repository.dart` `_mapPgError(rt, x, id: id)`), which names the wrong row. */
  async duplicateRefusal(row: number, e?: DomainValidationError): Promise<string> {
    const key = this._info.businessKey ?? [];
    const values = key.map((c) => this._wire(row, c)).filter((v) => v != null && `${v}` !== '');
    if (values.length !== key.length)
      return this._duplicateText();
    const self = `${this._wire(row, 'id') ?? ''}`;
    const id = `${e?.rows?.[0]?.id ?? ''}`;
    const owner = id !== '' && id !== self ? await this._nameOf(id) : await this._keyHolder(key, values, self);
    return owner == null ? this._duplicateText()
      : `${key.map((c) => this._captionOf(c)).join(', ')} ${values.join(', ')} already belongs to ${owner}`;
  }

  /** "Ethanol still has 3 containers; remove or reassign them first": the server names the child
   * table and the column pointing back here, the child table itself the count and the noun. */
  async restrictRefusal(e: DomainError): Promise<string> {
    const named = DomainFrameEditor.LIVE_ROWS.exec(`${e.message}`);
    if (named == null)
      return `${e.message}`;
    const id = `${e.body['id'] ?? ''}`;
    const at = this._rowById(id);
    const subject = (at < 0 ? null : this._displayOf(at)) ?? `This ${this._info.singularName.toLowerCase()}`;
    let count = 0;
    let info: DomainTableInfo | null = null;
    // the count and the plural are the nicety, not the refusal: the sentence stands without them
    try {
      [count, info] = await Promise.all([
        new DomainTableClient(api.grok_Dapi_Domains(), this.client.schema, named[1])
          .count(`${named[2]} = "${id}"`),
        registry.tableInfo(`${this.client.schema}.${named[1]}`)]);
    } catch (x) {
      log.warning(`${this.table}: naming the blocking children failed — ${x}`);
    }
    const noun = ((count === 1 ? info?.singularName : info?.pluralName) ?? named[1])
      .toLowerCase().replace(/_/g, ' ');
    return `${subject} still has ${count === 0 ? '' : `${count} `}${noun};`
      + ` remove or reassign ${count === 1 ? 'it' : 'them'} first`;
  }

  /** Opens/closes the editor around a transaction (see {@link save}); a
   * {@link DomainSession} sets it on every participant. */
  setSaving(saving: boolean): void {
    this._saving = saving;
    this._onSavingChanged.next(saving);
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
  async refresh(query?: DomainQuerySpec): Promise<DataFrame> {
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

  private _bind(df: DataFrame): void {
    for (const name of DomainFrameEditor.SERVICE_COLUMNS) {
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

  private _col(name: string): Column<string> {
    return this._df.columns.byName(name) as Column<string>;
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
   * {@link applyResults}, and a discard/refresh would shift the very rows those
   * results address. */
  private _busy(action: string): boolean {
    if (!this._saving)
      return false;
    balloon.warning(`${this.table}: ${action} is not available while the batch is being saved`);
    return true;
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

  /** The cell as the wire sees it. The grid's DELETE / BACKSPACE path clears a
   * cell by writing `''` — in a column that holds a uuid
   * ({@link isReferenceProperty}) that means "no value", so it goes out as an explicit
   * null and a required reference fails validation as the empty value it is.
   *
   * A DG string column already reports `''` as none, so the coercion below only
   * fires for a reference column backed by something else; it is the DART side
   * (`Cell.value` returns the `''`) that clears refs through this rule today.
   * Kept because the rule is the state contract's, not this column type's. */
  private _wire(row: number, column: string): any {
    const col = this._df.columns.byName(column);
    if (col == null || col.isNone(row))
      return null;
    const value = col.get(row);
    const p = this._propByName.get(column);
    if (value === '' && p != null && DomainFrameEditor.isReferenceProperty(p))
      return null;
    return toWire(value);
  }

  /** Whether [row] exists on the server: it carries an id that is not a draft. */
  private _isPersisted(row: number): boolean {
    const id = this._wire(row, 'id');
    return id != null && !DomainFrameEditor.isDraftId(id);
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
      this._setJson(DomainFrameEditor.CHANGES_COLUMN, row, changes);
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

  /** The cell against its registry {@link Property} alone. An EMPTY cell is
   * empty whatever its column type spells that as — a string column's `''`, an int
   * column's null sentinel — so nullability, and not the sentinel, decides
   * (otherwise an untouched optional choice column reports `""` is not one of ...`,
   * and an untouched optional int reports itself below its minimum). */
  private _invalidValue(row: number, column: string): string | null {
    const property = this._propByName.get(column);
    if (property == null)
      return null;
    const col = this._df.columns.byName(column);
    return DomainFrameEditor.validateCellValue(property,
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
    if (DomainFrameEditor.writableColumns(this.access).includes(column))
      return null;
    const pending = this.stateOf(row) === 'new'
      ? this._wire(row, column) != null && !this._isServiceValue(row, column)
      : column in this.changesOf(row);
    return pending ? `Column '${column}' is read-only` : null;
  }

  /** Whether the cell of a NEW row holds a value the EDITOR put there rather than
   * a user value — the draft id {@link addRow} stamps into `id`, which
   * {@link buildOps} leaves out of the insert and sends as the op's `ref`. Losing
   * it is the design, not a dropped write. */
  private _isServiceValue(row: number, column: string): boolean {
    return column === 'id' && DomainFrameEditor.isDraftId(this._wire(row, column));
  }

  private _setError(row: number, column: string, error: DomainCellError | null): void {
    const errors = this.errorsOf(row);
    if (error == null)
      delete errors[column];
    else
      errors[column] = error;
    this._setJson(DomainFrameEditor.ERRORS_COLUMN, row, errors);
  }

  /** Validates every column the FRAME carries — the declared ones against the
   * registry, and the undeclared ones (a column the query returned that the
   * registry does not describe) against column security alone: a prefilled value
   * on one of those would be dropped by {@link buildOps} just as silently. */
  private _validateRow(row: number): void {
    for (const column of this._df.columns.names()) {
      // Every `~` column is service state (the editor's own, a `withAccess` read's
      // `~can_*`, a relation's id companion) — never data a payload could drop.
      if (column.startsWith('~'))
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
    this._col(DomainFrameEditor.STATE_COLUMN).set(row, changed ? 'modified' : '', false);
    this._recount(row);
  }

  private _clearRowState(row: number): void {
    this._col(DomainFrameEditor.STATE_COLUMN).set(row, '', false);
    this._col(DomainFrameEditor.CHANGES_COLUMN).set(row, '', false);
    this._col(DomainFrameEditor.ERRORS_COLUMN).set(row, '', false);
    this._cache(DomainFrameEditor.CHANGES_COLUMN).delete(row);
    this._cache(DomainFrameEditor.ERRORS_COLUMN).delete(row);
    this._recount(row);
  }

  private _setDeleted(rows: number | number[], deleted: boolean): void {
    if (this._busy(deleted ? 'deleting' : 'restoring a row'))
      return;
    const list = Array.isArray(rows) ? rows : [rows];
    const state = this._col(DomainFrameEditor.STATE_COLUMN);
    for (const row of list) {
      if (row < 0 || row >= this._df.rowCount)
        continue;
      if (deleted)
        state.set(row, 'deleted', false);
      else if (state.get(row) === 'deleted')
        // A row ADDED in this batch never reached the server — buildOps' own
        // predicate — and its values live nowhere but the frame, so it must go back
        // to 'new'; a 'new' row carries no ~changes by design, which is why the
        // change test below would otherwise land it on '' and make it invisible to
        // save AND discard.
        state.set(row, !this._isPersisted(row) ? 'new'
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
    const state = this._df.columns.byName(DomainFrameEditor.STATE_COLUMN);
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

  /** The row's display value for a dialog caption: the registry's declared name
   * column when the frame carries it — the same identity every other platform
   * surface shows. Null falls the caller back to the id. */
  /** The caption a user knows [column] by. */
  private _captionOf(column: string): string {
    const caption = this._propByName.get(column)?.friendlyName;
    return caption == null || caption === '' ? column : caption;
  }

  /** The refusal a business-key clash gets where the row already holding the key is unknown. */
  private _duplicateText(): string {
    const key = (this._info.businessKey ?? []).map((c) => this._captionOf(c)).join(', ');
    return `A ${this._info.singularName.toLowerCase()} with this ${key} already exists`;
  }

  /** The display name of a row of this table, read back by id; null when the table has no name
   * column, or the row is gone or invisible. */
  private async _nameOf(id: string): Promise<string | null> {
    const column = this._info.nameColumn;
    if (column == null)
      return null;
    try {
      const row = await this.client.get(id) as {[column: string]: any};
      const value = row?.[column];
      return value == null || `${value}` === '' ? null : `${value}`;
    } catch (x) {
      log.warning(`${this.table}: naming the row holding the key failed — ${x}`);
      return null;
    }
  }

  /** The display name of the row already holding [values] on the business key, read back by the
   * key itself and skipping [selfId] (the row being saved); null when the table has no name
   * column, or nobody else holds the key where the caller can see it. */
  private async _keyHolder(key: string[], values: any[], selfId: string): Promise<string | null> {
    const column = this._info.nameColumn;
    if (column == null)
      return null;
    const filter: DomainConditionTree = [];
    for (let i = 0; i < key.length; i++) {
      if (filter.length > 0)
        filter.push('and');
      filter.push({property: key[i], operator: '=', value: values[i]});
    }
    try {
      const rows = await this.client.query({filter: filter, columns: [column], limit: 2}) as any[];
      for (const r of rows) {
        if (`${r['id'] ?? ''}` === selfId)
          continue;
        const value = r[column];
        return value == null || `${value}` === '' ? null : `${value}`;
      }
    } catch (x) {
      log.warning(`${this.table}: naming the row holding the key failed — ${x}`);
    }
    return null;
  }

  /** The frame row carrying [id], or -1. */
  private _rowById(id: string): number {
    for (let row = 0; row < this._df.rowCount; row++)
      if (`${this._wire(row, 'id')}` === id)
        return row;
    return -1;
  }

  private _displayOf(row: number): string | null {
    const name = this._info.nameColumn;
    if (name == null || !this._df.columns.contains(name))
      return null;
    const v = this._wire(row, name);
    return v == null || `${v}` === '' ? null : `${v}`;
  }
}

/** What the user chose in {@link promptUnsavedChanges} (dismissing it is `'cancel'`). */
export type UnsavedOutcome = 'save' | 'discard' | 'cancel';

/** Wording of {@link promptUnsavedChanges}. */
export interface UnsavedPromptOptions {
  /** What the user is about to do — 'leave this page', 'change the filter';
   * defaults to 'continue'. */
  action?: string;
  /** Names the pending batch; defaults to the editors' table addresses. */
  subject?: string;
}

/**
 * Asks what to do with the pending changes of [editors] before something that
 * would discard them. Resolves to the user's decision WITHOUT applying it —
 * {@link confirmDiscardChanges} is the one that saves or discards.
 */
export function promptUnsavedChanges(editors: DomainFrameEditor[],
  options?: UnsavedPromptOptions): Promise<UnsavedOutcome> {
  const dirty = (editors ?? []).filter((e) => e != null && e.isDirty);
  const changes = dirty.reduce((n, e) => n + e.changeCount, 0);
  const subject = options?.subject ?? dirty.map((e) => e.table).join(', ');
  const action = options?.action ?? 'continue';
  return new Promise<UnsavedOutcome>((resolve) => {
    // Every path goes through one guard: a double-clicked button (or a close
    // notification arriving after a decision) must repeat the first answer, not
    // overwrite it.
    let decided = false;
    let closed: {unsubscribe(): void} | null = null;
    const decide = (outcome: UnsavedOutcome) => {
      // The onClose subscription outlives the dialog otherwise: its subject is
      // the platform's, not this promise's.
      closed?.unsubscribe();
      closed = null;
      if (!decided) {
        decided = true;
        resolve(outcome);
      }
    };
    const dlg = Dialog.create({title: 'Unsaved changes'});
    dlg.add(ui.divV([
      ui.p(`${changes} unsaved change${changes === 1 ? '' : 's'} in ${subject}.`),
      ui.p(`Save them, discard them, or cancel and do not ${action}.`),
    ], 'ui-hint-block'));
    // Explicit indices: buttons are INSERTED into the command bar (which already
    // carries CANCEL), so the default index 0 would read DISCARD, SAVE, CANCEL.
    dlg.addButton('SAVE', () => {
      decide('save');
      dlg.close();
    }, 0);
    dlg.addButton('DISCARD', () => {
      decide('discard');
      dlg.close();
    }, 1);
    dlg.onCancel(() => decide('cancel'));
    // A dismissal (the X, Esc, a programmatic close) is a cancel — and `onClose`
    // fires for the decided paths too, where the guard keeps the first answer.
    closed = dlg.onClose.subscribe(() => decide('cancel'));
    dlg.show();
  });
}

/**
 * THE gate an app puts in front of a rebuild: resolves to whether the caller may
 * proceed, having saved or discarded whatever was pending.
 *
 * - nothing dirty — proceeds silently, no dialog;
 * - a save in flight — refuses (an editor is closed while its transaction runs,
 *   and prompting would offer a save that cannot be taken); the caller may retry
 *   after `onSavingChanged`;
 * - otherwise it prompts and applies the answer: `save` writes every pending
 *   batch (a FAILED save cancels the navigation — the user keeps their changes
 *   and the error), `discard` drops them, `cancel` stays.
 *
 * A {@link DomainSession} saves as ONE transaction; a plain list of editors
 * saves each on its own (a session of one apiece).
 */
export async function confirmDiscardChanges(editors: DomainFrameEditor[] | DomainSession,
  options?: UnsavedPromptOptions): Promise<boolean> {
  const session = editors instanceof DomainSession ? editors : null;
  const all = (session?.editors ?? editors as DomainFrameEditor[] ?? []).filter((e) => e != null);
  // Mid-save first, and for EVERY editor: one with a transaction in flight
  // refuses writes, discard and refresh anyway, and prompting would offer a
  // save that cannot be taken. Retry after `onSavingChanged`.
  if (all.some((e) => e.isSaving)) {
    balloon.warning('Wait for the batch being saved to finish.');
    return false;
  }
  const dirty = all.filter((e) => e.isDirty);
  if (dirty.length === 0)
    return true;
  const outcome = await promptUnsavedChanges(dirty, options);
  if (outcome === 'cancel')
    return false;
  if (outcome === 'discard') {
    if (session != null)
      session.discard();
    else
      for (const editor of dirty)
        editor.discard();
    return true;
  }
  if (session != null)
    return await session.save();
  for (const editor of dirty)
    if (!(await editor.save()))
      return false;
  return true;
}
