/**
 * {@link DomainForm} — the property form of ONE domain row, as a `DG.Widget`:
 * inputs generated from the registry, reference pickers for FK / user / group
 * columns ({@link RefInput}), sync and async validation, and a save that goes
 * through a one-row {@link DomainFrameEditor} — the same single writer the grid
 * uses, so dirty state, validation, the 409 flow and the transaction are shared.
 *
 * ```ts
 * const issues = await domains.table('grit.issue');
 * const form = issues.form({values: {reporter: me.id}});
 * grok.shell.newView('New issue', [form.root]);
 * ```
 *
 * It is also the feature's machine surface: `getProperties()` / `props` write
 * through the editor exactly as a keystroke does, `getWidgetStatus()` reports every
 * input live, and `getFunctions()` returns real platform Funcs (`Save`, `Discard`,
 * `Reset`) — see {@link DomainForm.getWidgetStatus}.
 *
 * @module form
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {DomainFrameEditor, IDomainTableContext, IEditorHost,
  isReferenceProperty} from './frame-editor';
import {discardFunc, resetFunc, saveFunc} from './actions';
import {applyDomainUiStyles} from './styles';

/** How many suggestions a {@link RefInput} asks the server for. */
export const REF_SUGGESTION_LIMIT = 10;

/** How many users/groups a {@link RefInput} lists to suggest from. */
const PRINCIPAL_PAGE_SIZE = 1000;

const UUID = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

/** Options of `domains.table(...).form()` / `formDialog()`. */
export interface DomainFormOptions {
  /** Values a NEW row starts with — the insert form (the default). */
  values?: {[column: string]: any};
  /** Edit this row instead of creating one. */
  row?: DG.DomainRow;
  /** Edit the row with this id. */
  id?: string;
}


// ─────────────────────── reference inputs ─────────────────────────

/** Where a {@link RefInput}'s suggestions, display names and picker come from —
 * one implementation per reference kind (a domain table, a user, a group). */
interface RefSource {
  suggest(query: string): Promise<{label: string; value: string}[]>;
  /** Display text of a known id; null when it does not resolve. */
  display(id: string): Promise<string | null>;
  /** The full picker, for targets that have one; [anchor] asks for the drop-down
   * flavour, anchored below it. */
  pick?(anchor?: HTMLElement): Promise<{id: string; label: string} | null>;
}

/** Suggestions from the target table's own rows, display names from the registry's
 * batched resolver, and the platform's lookup picker as the full search. */
class DomainRefSource implements RefSource {
  private _nameColumn: string | null | undefined = undefined;

  constructor(private readonly table: string) {}

  private async _column(): Promise<string | null> {
    if (this._nameColumn === undefined)
      this._nameColumn = (await grok.dapi.domains.registry.tableInfo(this.table)).nameColumn;
    return this._nameColumn;
  }

  async suggest(query: string): Promise<{label: string; value: string}[]> {
    const column = await this._column();
    const client = grok.dapi.domains.table(this.table);
    // '~*' is the case-insensitive match; a table with no display-name column can
    // only offer its first rows (the picker is the way in for those).
    const rows: any[] = column == null ? await client.query({limit: REF_SUGGESTION_LIMIT})
      : await client.query().where(column, '~*', query).top(REF_SUGGESTION_LIMIT);
    const ids = rows.map((r) => `${r['id']}`);
    const names = await grok.dapi.domains.registry.resolveNames(this.table, ids);
    return ids.map((id) => ({label: names[id] ?? id, value: id}));
  }

  async display(id: string): Promise<string | null> {
    return (await grok.dapi.domains.registry.resolveNames(this.table, [id]))[id] ?? null;
  }

  async pick(anchor?: HTMLElement): Promise<{id: string; label: string} | null> {
    const row = await DG.DomainObjectHandler.pickRow(this.table, {anchor: anchor});
    return row == null ? null : {id: row.id, label: row.displayName};
  }
}

/** Users and groups: the platform's own lists, fetched once and matched locally
 * (there is no server-side "starts with" endpoint for them). */
class PrincipalRefSource implements RefSource {
  private _all: Promise<{label: string; value: string}[]> | null = null;

  constructor(private readonly kind: 'User' | 'Group') {}

  private _list(): Promise<{label: string; value: string}[]> {
    if (this._all == null)
      this._all = (this.kind === 'User' ? grok.dapi.users.list({pageSize: PRINCIPAL_PAGE_SIZE})
        : grok.dapi.groups.list({pageSize: PRINCIPAL_PAGE_SIZE}))
        .then((items: any[]) => items.map((x) => ({label: `${x.friendlyName ?? x.name}`, value: `${x.id}`})));
    return this._all;
  }

  async suggest(query: string): Promise<{label: string; value: string}[]> {
    const text = `${query}`.toLowerCase();
    return (await this._list())
      .filter((x) => x.label.toLowerCase().includes(text))
      .slice(0, REF_SUGGESTION_LIMIT);
  }

  async display(id: string): Promise<string | null> {
    return (await this._list()).find((x) => x.value === id)?.label ?? null;
  }
}

/** The source for a reference property's semType (`'User'`, `'Group'`, or a
 * `'<schema>.<table>'` address). */
export function refSourceFor(target: string): RefSource {
  return target === 'User' || target === 'Group'
    ? new PrincipalRefSource(target) : new DomainRefSource(target);
}

/** A `User` / `Group` column's visible editor: the platform's own picker — the one
 * the Dart row editor uses (`DomainRowEditor._coreRefInput`) — and which principal
 * its selection points at. */
interface PrincipalInput {
  input: DG.InputBase;
  isUser: boolean;
}

/** The id a principal picker's selection writes: a user column stores the USER's id
 * (the picker hands out the personal GROUP wrapping it — what the Dart editor's
 * `items.first.user?.id` reads), a group column the group's own. */
function principalId(selected: any[], isUser: boolean): string | null {
  if (selected.length === 0)
    return null;
  let id: any = null;
  try {
    const item = isUser ? selected[0]?.user ?? selected[0] : selected[0];
    id = item?.id;
  } catch (_) {
    id = null;      // a group with nothing behind its `user` handle
  }
  return typeof id === 'string' && id !== '' ? id : null;
}

/** What a display text resolved to: an id, or the reason it did not. */
export interface RefResolution {
  id: string | null;
  error: string | null;
}

/** Options of {@link RefInput}. */
export interface RefInputOptions {
  /** What the reference points at: a domain table address, `'User'` or `'Group'`. */
  target: string;
  /** Whether the value may be cleared (offers the clear icon; default true). */
  nullable?: boolean;
  caption?: string;
}

/**
 * The reference picker: a type-ahead over the target's display names, a search
 * icon opening the platform's full lookup picker, and a clear icon for nullable
 * references — the Dart row editor's `_refInput` UX as a reusable JS input.
 *
 * The visible control is {@link input} (a `DG.TypeAhead`, so it composes into any
 * `ui.inputs(...)` form); the VALUE it edits is the referenced row's id, which is
 * what {@link id} reports and {@link setId} writes. {@link resolve} is the
 * display-text → id direction — how a machine write (`props['assignee'] = 'Alex'`)
 * takes the same path a user's pick does.
 */
export class RefInput {
  /** Mount this (or hand it to `ui.inputs`). */
  readonly input: DG.TypeAhead;
  readonly target: string;
  readonly nullable: boolean;

  private readonly _source: RefSource;
  private readonly _onResolved = new rxjs.Subject<RefResolution>();
  private _id: string | null = null;
  /** The text the input shows for {@link _id} — what a blur compares against to
   * tell a typed-over value from the one the input is already pointing at. */
  private _display = '';
  /** Guards a display lookup against a newer write. */
  private _generation = 0;

  constructor(options: RefInputOptions) {
    this.target = options.target;
    this.nullable = options.nullable ?? true;
    this._source = refSourceFor(options.target);
    this.input = ui.typeAhead(options.caption ?? options.target, {
      minLength: 1,
      limit: REF_SUGGESTION_LIMIT,
      source: (query: string) => this._source.suggest(query),
      onSubmit: (_: Event, item?: any) => {
        if (item != null)
          this._apply(`${item['value']}`, `${item['label']}`, true);
      },
    });
    this.input.root.classList.add('domain-ui-ref-input');
    if (this._source.pick != null)
      this.input.addOptions(ui.icons.search(() => this.pick(), `Pick a ${this.target}`));
    if (this.nullable)
      this.input.addOptions(ui.icons.close(() => this._apply(null, '', true), 'Clear'));
    // Typing a name and tabbing away is a pick too — the suggestion list is a
    // convenience, not the only way in. Only a USER gesture reaches this, so it
    // cannot loop with the display lookups {@link _apply} runs.
    this.input.input.addEventListener('blur', () => {
      const text = `${this.input.value ?? ''}`;
      if (text !== this._display)
        this.resolveInto(text);
    });
  }

  /** The referenced row's id, or null. */
  get id(): string | null { return this._id; }

  /** Fires on every outcome of a USER gesture (a picked suggestion, the picker, the
   * clear icon, a typed name resolved on blur): the new id, or the reason it could
   * not be resolved. ONE channel, so a form never has to guess whether the last
   * gesture landed. */
  get onResolved(): rxjs.Observable<RefResolution> { return this._onResolved; }

  /** Points the input at [id], resolving its display text when [display] is not
   * given. [notify] false is the sync-from-the-model direction — it must not write
   * back into the model it came from. */
  setId(id: string | null, display?: string, notify: boolean = true): void {
    this._apply(id, display, notify);
  }

  /** Opens the platform's lookup picker for the target table as a drop-down below
   * this input (right-aligned with its editor box). */
  async pick(): Promise<void> {
    if (this._source.pick == null)
      return;
    const picked = await this._source.pick(this.input.input);
    if (picked != null)
      this._apply(picked.id, picked.label, true);
  }

  /**
   * The display text → id direction, through the SAME suggestion source the
   * type-ahead uses: an id passes through, an unambiguous display name resolves,
   * and anything else comes back as a validation message instead of a throw
   * (`{id: null, error}`) — a machine write is never allowed to fail louder than
   * typing the same thing would.
   */
  async resolve(text: any): Promise<RefResolution> {
    const query = `${text ?? ''}`.trim();
    if (query === '')
      return {id: null, error: this.nullable ? null : 'Value can\'t be empty'};
    if (UUID.test(query) && (await this._displaySafe(query)) != null)
      return {id: query, error: null};
    let items: {label: string; value: string}[];
    try {
      items = await this._source.suggest(query);
    } catch (e: any) {
      return {id: null, error: `Cannot look ${this.target} up: ${e?.message ?? e}`};
    }
    const exact = items.filter((i) => i.label.toLowerCase() === query.toLowerCase());
    const candidates = exact.length > 0 ? exact : items;
    if (candidates.length === 0)
      return {id: null, error: `No ${this.target} matches "${query}"`};
    if (candidates.length > 1)
      return {id: null, error: `"${query}" matches ${candidates.length} ${this.target} rows — ` +
        'type more of the name, or pick one'};
    return {id: candidates[0].value, error: null};
  }

  /** {@link resolve}, applied to the input: an unambiguous match becomes the value,
   * anything else leaves the text alone for the form's validator to explain.
   * Resolves to the outcome, so a caller can wait for it. */
  async resolveInto(text: any): Promise<RefResolution> {
    const resolved = await this.resolve(text);
    if (resolved.error == null)
      this._apply(resolved.id, undefined, true);
    else
      this._onResolved.next(resolved);
    return resolved;
  }

  private async _displaySafe(id: string): Promise<string | null> {
    try {
      return await this._source.display(id);
    } catch (_) {
      return null;
    }
  }

  private _apply(id: string | null, display: string | undefined, notify: boolean): void {
    this._id = id == null || id === '' ? null : id;
    const generation = ++this._generation;
    if (display !== undefined)
      this._setText(display);
    else if (this._id == null)
      this._setText('');
    else
      this._displaySafe(this._id).then((text) => {
        if (generation === this._generation)
          this._setText(text ?? `${this._id}`);
      });
    if (notify)
      this._onResolved.next({id: this._id, error: null});
  }

  private _setText(text: string): void {
    this._display = text;
    this.input.value = text;
  }
}


// ─────────────────────── async validation ─────────────────────────

/**
 * Binds an ASYNC probe to a sync input validator — the standard round trip for a
 * server-side check (a uniqueness probe, a cross-row rule): debounce the input,
 * run the probe, store its answer in a slot the sync validator reads, re-validate.
 * A stale probe never overwrites a newer one.
 *
 * ```ts
 * bindAsyncValidator(form.getInput('number')!,
 *   async (n) => await issues.query().where({number: n}).exists() ? 'Already taken' : null);
 * ```
 *
 * Returns the subscription — unsubscribe it with the input (a detached form does).
 */
export function bindAsyncValidator(input: DG.InputBase, probe: (value: any) => Promise<string | null>,
  options?: {debounce?: number}): {unsubscribe(): void} {
  let error: string | null = null;
  let timer: any = null;
  let generation = 0;
  input.addValidator(() => error);
  const sub = input.onInput.subscribe(() => {
    error = null;
    if (timer != null)
      clearTimeout(timer);
    const id = ++generation;
    timer = setTimeout(async () => {
      let message: string | null;
      try {
        message = await probe(input.value);
      } catch (e: any) {
        message = `${e?.message ?? e}`;
      }
      if (id !== generation)
        return;
      error = message;
      input.validate();
    }, options?.debounce ?? 300);
  });
  return {unsubscribe: () => {
    if (timer != null)
      clearTimeout(timer);
    generation++;
    sub.unsubscribe();
  }};
}


// ─────────────────────── server errors ─────────────────────────

/** What {@link mapServerError} decided: a version conflict resolves to the user's
 * choice (which the CALLER applies), everything else is already reported. */
export type ServerErrorOutcome = 'reload' | 'overwrite' | 'handled';

/** Where {@link mapServerError} puts what it maps. */
export interface ServerErrorContext {
  /** Marks the named field (null clears it) — the per-field markers. */
  setFieldError?: (column: string, message: string | null) => void;
  /** Shows a form-level message (null clears it); a balloon when omitted. */
  banner?: (message: string | null) => void;
  /** What the error is about (a row's display value) — the 409 dialog's caption. */
  subject?: string;
  /** Columns the per-field mapping may address; everything else goes to the banner. */
  columns?: string[];
}

/**
 * THE mapping of a typed `/domains` failure onto the UI, per ARCHITECTURE §9 —
 * so no app hand-rolls error display:
 *
 * | error | UI |
 * |---|---|
 * | `DomainValidationError` with `rows[i].errors[{column}]` | per-field marker; unmapped residue → banner |
 * | cross-field / row-level validation (no column) | banner |
 * | duplicate (`isDuplicate`) | banner naming the conflict |
 * | `DomainVersionConflictError` (409) | the platform's reload/overwrite dialog (its answer is RETURNED) |
 * | `DomainRestrictError` | banner naming the referencing rows |
 * | `DomainForbiddenError` (403) | balloon + capability-cache invalidation, so affordances re-gate |
 * | anything else | balloon |
 */
export async function mapServerError(e: any, ctx?: ServerErrorContext): Promise<ServerErrorOutcome> {
  const banner = ctx?.banner ?? ((message: string | null) => {
    if (message != null)
      grok.shell.error(message);
  });
  if (e instanceof DG.DomainVersionConflictError) {
    const decision = await DG.DomainObjectHandler.showConflictDialog(ctx?.subject ?? 'this row');
    return decision === 'reload' || decision === 'overwrite' ? decision : 'handled';
  }
  if (e instanceof DG.DomainValidationError) {
    const known = ctx?.columns;
    const unmapped: string[] = [];
    let mapped = false;
    for (const row of e.rows ?? [])
      for (const columnError of row.errors ?? []) {
        const column = columnError.column;
        if (column != null && ctx?.setFieldError != null && (known == null || known.includes(column))) {
          ctx.setFieldError(column, columnError.message);
          mapped = true;
        }
        else
          unmapped.push(columnError.message);
      }
    if (e.isDuplicate && unmapped.length === 0 && !mapped)
      unmapped.push(e.message);
    banner(unmapped.length > 0 ? unmapped.join('\n') : mapped ? null : e.message);
    return 'handled';
  }
  if (e instanceof DG.DomainRestrictError || e instanceof DG.DomainNotFoundError) {
    banner(e.message);
    return 'handled';
  }
  if (e instanceof DG.DomainForbiddenError) {
    // The affordances were built from a capability snapshot the server just
    // contradicted — drop it, so the next read re-gates against server truth.
    grok.dapi.domains.invalidateUiCaches();
    grok.shell.error(e.message);
    return 'handled';
  }
  grok.shell.error(`${e?.message ?? e}`);
  return 'handled';
}


// ─────────────────────── the form ─────────────────────────

/** A validator of one form field: sync or async, `null` when the value is fine. */
export type DomainFieldValidator = (value: any) => string | null | Promise<string | null>;

/**
 * The property form of one domain row — see the module doc.
 *
 * **Construction is synchronous**, on a prefetched {@link IDomainTableContext}: the
 * inputs exist immediately (so `getInput`, `replaceInput`, `props` and
 * `getWidgetStatus` all work at once), while the row itself loads behind the
 * platform's update indicator. {@link ready} resolves when the loading is over —
 * INCLUDING when it failed, in which case the form shows the reason instead of
 * rejecting into nobody's catch.
 *
 * **One writer.** Every path into a value — typing, `props['x'] = v`,
 * {@link setValue} — writes through the one-row {@link DomainFrameEditor}, so
 * validation, dirty state and the save transaction cannot disagree with each other.
 */
export class DomainForm extends DG.Widget implements IEditorHost {
  private readonly _context: IDomainTableContext;
  private readonly _options: DomainFormOptions;
  private readonly _banner = ui.div([], 'domain-ui-form-banner');
  private readonly _host = ui.div([], 'domain-ui-form-inputs');
  private readonly _inputs = new Map<string, DG.InputBase>();
  private readonly _refs = new Map<string, RefInput>();
  /** The `User` / `Group` columns among {@link _refs} — the ones whose VISIBLE input
   * is the platform picker rather than the {@link RefInput}'s type-ahead. */
  private readonly _principals = new Map<string, PrincipalInput>();
  /** Guards a picker's seeding against a newer write ({@link _seedPrincipal}). */
  private readonly _seedGeneration = new Map<string, number>();
  private readonly _factories = new Map<string, (p: DG.Property) => DG.InputBase>();
  private readonly _validators = new Map<string, DomainFieldValidator[]>();
  private readonly _validatorGeneration = new Map<string, number>();
  private readonly _customErrors = new Map<string, string | null>();
  private readonly _refErrors = new Map<string, string | null>();
  private readonly _touched = new Set<string>();
  /** Writes made before the editor existed — applied to the row when it does. */
  private readonly _pending: {[column: string]: any} = {};
  /** Everything the form is still waiting for (async validators, ref resolutions):
   * what {@link save} drains before deciding whether the row is valid. */
  private readonly _work = new Set<Promise<any>>();
  private readonly _subs: {unsubscribe(): void}[] = [];
  /** Subscriptions of the CURRENT inputs — dropped and re-made whenever
   * {@link replaceInput} rebuilds them. */
  private _inputSubs: {unsubscribe(): void}[] = [];
  private readonly _ready: Promise<DomainForm>;

  private _editor: DomainFrameEditor | null = null;
  private _formProperties: DG.Property[] | null = null;
  private _loadError: string | null = null;
  private _row = 0;
  /** Whether input value writes are the form syncing itself FROM the model. */
  private _syncing = false;

  constructor(context: IDomainTableContext, options?: DomainFormOptions) {
    super(ui.divV([], 'domain-ui-form'));
    applyDomainUiStyles();
    this._context = context;
    this._options = options ?? {};
    ui.setDisplay(this._banner, false);
    this.root.appendChild(this._banner);
    this.root.appendChild(this._host);
    this._buildInputs();
    // What makes the form discoverable through `DG.Widget.getAll()`: wrapping the
    // widget for Dart is what registers its root in the platform's widget map.
    this.toDart();
    ui.setUpdateIndicator(this.root, true, `Loading ${this._context.info.singularName}...`);
    this._ready = this._load().then(() => this, () => this);
  }

  /** The platform type name — what `DG.Widget.getAll().find((w) => w.type === 'DomainForm')`
   * keys on. */
  get type(): string { return 'DomainForm'; }

  /** The table the row belongs to, `'<schema>.<table>'`. */
  get table(): string { return this._context.table; }

  /** The one-row editor behind the form; null until {@link ready} resolves. */
  get editor(): DomainFrameEditor | null { return this._editor; }

  /** {@link IEditorHost}: the form's pending changes, for a page's dirty rollup. */
  get editors(): DomainFrameEditor[] { return this._editor == null ? [] : [this._editor]; }

  /** Whether the form edits an existing row (as opposed to creating one). */
  get isEditing(): boolean { return this._options.row != null || this._options.id != null; }

  /** Id of the row being edited, or null for an insert form. */
  get rowId(): string | null {
    return this._options.id ?? this._options.row?.id ?? null;
  }

  /** Whether anything has been changed and not yet saved. */
  get isDirty(): boolean { return this._editor?.isDirty === true; }

  /**
   * Resolves when the row has loaded — or when loading FAILED, in which case the
   * form shows the reason and this still resolves (a rejected promise nobody awaits
   * is an unhandled rejection, and a widget's loading is not the caller's business
   * unless it asks). Idempotent: every call returns the same promise.
   */
  get ready(): Promise<this> { return this._ready as Promise<this>; }

  /** The input editing [name], or null when the column is absent (not declared, or
   * not writable for this caller). */
  getInput(name: string): DG.InputBase | null { return this._inputs.get(name) ?? null; }

  /** The {@link RefInput} editing [name], or null when the column is not a reference. */
  getRefInput(name: string): RefInput | null { return this._refs.get(name) ?? null; }

  /**
   * Replaces the generated input for ONE column, keeping the rest of the form —
   * the customization point of §F.1.6c:
   *
   * ```ts
   * form.replaceInput('description', (p) => ui.input.textArea(p.caption));
   * ```
   *
   * The replacement is wired exactly like a generated one (writes go through the
   * editor, validation shows on it), and applies immediately.
   */
  replaceInput(name: string, factory: (p: DG.Property) => DG.InputBase): void {
    this._factories.set(name, factory);
    this._buildInputs();
    this._syncInputs();
  }

  /**
   * Adds a validator to ONE field — sync or async, rendered as the standard
   * per-field marker and reported through {@link getWidgetStatus}:
   *
   * ```ts
   * form.addValidator('title', (s) => s.trim().length < 5 ? 'At least 5 characters' : null);
   * ```
   *
   * An async validator is awaited by {@link save}, so a slow probe cannot let an
   * invalid row through.
   */
  addValidator(name: string, validator: DomainFieldValidator): void {
    const existing = this._validators.get(name);
    if (existing == null)
      this._validators.set(name, [validator]);
    else
      existing.push(validator);
    this._runValidators(name);
  }

  /** Writes [value] into [name] THROUGH the editor — the path a keystroke takes.
   * A reference column accepts a display name as well as an id (see
   * {@link RefInput.resolve}); an unresolvable one lands as a validation error, not
   * as a throw. */
  setValue(name: string, value: any): void {
    if (this._refs.has(name))
      this._writeRef(name, value);
    else {
      this._write(name, value);
      this._syncInput(name);
    }
  }

  /** The current value of [name] — the editor's, or the buffered write when the row
   * is still loading. */
  getValue(name: string): any {
    if (this._editor == null)
      return name in this._pending ? this._pending[name]
        : (this._options.values ?? {})[name] ?? null;
    const column = this._editor.dataFrame.columns.byName(name);
    return column == null || column.isNone(this._row) ? null : column.get(this._row);
  }

  /**
   * Writes the row (one transaction, the editor's own typed-error flows: per-field
   * markers, the duplicate banner, the standard 409 reload/overwrite dialog).
   * Resolves to whether it landed.
   *
   * Waits for the row to load and for every pending async validation / reference
   * resolution first, then refuses on the first blocking error.
   */
  async save(): Promise<boolean> {
    await this._ready;
    if (this._editor == null) {
      grok.shell.error(this._loadError ?? `${this.table}: the form did not load`);
      return false;
    }
    await this.settle();
    const blocking = this._firstError();
    if (blocking != null) {
      grok.shell.error(`Cannot save: ${blocking}`);
      return false;
    }
    return await this._editor.save();
  }

  /** Puts the form back to the values it was opened with. For a one-row form this
   * IS the discard: a new row must survive it (dropping it would leave the form
   * without anything to edit), so it is re-seeded — and pristine again, exactly
   * as it was when the form opened. */
  reset(): void {
    this._clearBanner();
    this._customErrors.clear();
    this._refErrors.clear();
    this._touched.clear();
    if (this._editor == null) {
      for (const name of Object.keys(this._pending))
        delete this._pending[name];
      this._syncInputs();
      return;
    }
    if (this.isEditing)
      this._editor.revertRow(this._row);
    else {
      // Dropping the row and re-adding it (rather than writing every value back)
      // is what restores the PRISTINE state: a reset form must no more arm the
      // unsaved-changes gate than a freshly opened one does.
      this._editor.discard();
      this._editor.addRow(Object.assign({}, this._options.values ?? {}), {pristine: true});
    }
    this._syncInputs();
  }

  /** {@link reset} — a one-row form's discard. */
  discard(): void { this.reset(); }

  /** Waits for every pending async validator and reference resolution. */
  async settle(): Promise<void> {
    // Each round can start new work (a resolution that triggers a validator); the
    // cap keeps a misbehaving validator from spinning here forever.
    for (let i = 0; i < 8 && this._work.size > 0; i++)
      await Promise.all(Array.from(this._work));
  }

  // ─────────────────────── the machine surface ─────────────────────────

  /**
   * The registry-derived properties, WIRED THROUGH THE EDITOR: reading one reads
   * the row, writing one takes the single-writer path a keystroke takes. That makes
   * the platform's own `props` bag the form's named-value surface —
   * `form.props['title'] = 'Login page 500s'` — for the property grid, tests and
   * the AI alike.
   *
   * Available immediately (they come from the prefetched registry metadata): a
   * write made before the row has loaded is buffered and applied to it.
   */
  getProperties(): DG.Property[] {
    if (this._formProperties == null)
      this._formProperties = this.editableProperties().map((p) => this._propertyFor(p));
    return this._formProperties;
  }

  /**
   * The form's runtime structure — the platform's widget protocol, with the
   * additive `inputs` field carrying every field's LIVE value and validation
   * (read on the spot, never cached):
   *
   * - `parts` — every input's root by column name, plus `'form'` and `'banner'`;
   * - `description` — natural language derived from the registry (what the form
   *    creates or edits, its fields, and why it is read-only when it is);
   * - `error` — the first pending validation error;
   * - `inputs` — `IInputStatus` per field, with `type: 'ref'` and `ref` naming the
   *    target for reference columns.
   */
  getWidgetStatus(): DG.IWidgetStatus {
    const parts: {[name: string]: Element} = {form: this._host};
    for (const [name, input] of this._inputs)
      parts[name] = input.root;
    if (this._banner.style.display !== 'none')
      parts['banner'] = this._banner;
    return {
      parts: parts,
      hitAreas: {},
      shortcuts: {},
      events: [],
      description: this._describe(),
      error: this._loadError ?? this._firstError(),
      inputs: this.getInputStatuses(),
    };
  }

  /** The `inputs` part of {@link getWidgetStatus}, on its own. */
  getInputStatuses(): DG.IInputStatus[] {
    const statuses: DG.IInputStatus[] = [];
    for (const p of this.editableProperties()) {
      const error = this._errorOf(p.name);
      const reference = isReferenceProperty(p);
      const choices = p.choices;
      statuses.push({
        name: p.name,
        caption: p.friendlyName ?? p.name,
        type: reference ? 'ref' : `${p.propertyType}`,
        semType: p.semType == null || p.semType === '' ? undefined : `${p.semType}`,
        value: this.getValue(p.name),
        choices: choices != null && choices.length > 0 ? choices : undefined,
        required: p.nullable !== true,
        valid: error == null,
        error: error ?? undefined,
        description: p.description == null || p.description === '' ? undefined : p.description,
        ref: reference ? `${p.semType}` : undefined,
      });
    }
    return statuses;
  }

  /** The actions this caller may perform on this form, as REAL platform Funcs —
   * `Save`, `Discard`, `Reset` from the shared vocabulary, or none at all when the
   * form is read-only (no writable column, or no `canInsert`/`canEdit`). */
  getFunctions(): DG.Func[] {
    return this.canWrite ? [saveFunc(), discardFunc(), resetFunc()] : [];
  }

  /** Columns this form edits: the table's declared columns the caller may WRITE
   * (column security), in declared order. */
  editableProperties(): DG.Property[] {
    const writable = this._context.capabilities.writableColumns;
    return this._context.properties.filter((p) => writable.includes(p.name));
  }

  /** Whether the caller may save this form at all — the capability gate every
   * affordance derives from. */
  get canWrite(): boolean {
    const capabilities = this._context.capabilities;
    return this.editableProperties().length > 0 &&
      (this.isEditing ? capabilities.canEdit : capabilities.canInsert);
  }

  detach(): void {
    for (const sub of this._subs.concat(this._inputSubs))
      sub.unsubscribe();
    this._subs.length = 0;
    this._inputSubs = [];
    this._editor?.detach();
    super.detach();
  }

  // ─────────────────────── internals ─────────────────────────

  private async _load(): Promise<void> {
    try {
      const options = {capabilities: this._context.capabilities};
      if (this.isEditing) {
        const id = this.rowId;
        this._editor = await DomainFrameEditor.create(this._context.client, Object.assign(
          {query: {filter: DG.cond('id', '=', id), limit: 1}}, options));
        if (this._editor.dataFrame.rowCount === 0)
          throw new DG.DomainNotFoundError(`${this.table} ${id} does not exist, or you cannot see it`,
            404, {});
      }
      else {
        // Local, no round trip: everything `forRows` awaits is in the prefetched
        // context. The row is PRISTINE — an untouched New form, however prefilled,
        // is not a pending change (see `DomainFrameEditor.addRow`).
        this._editor = DomainFrameEditor.forContext(this._context, options);
        this._editor.addRow(Object.assign({}, this._options.values ?? {}, this._pending),
          {pristine: Object.keys(this._pending).length === 0});
      }
      if (this.isEditing)
        for (const name of Object.keys(this._pending))
          if (this._editor.dataFrame.columns.contains(name))
            this._editor.setValue(this._row, name, this._pending[name]);
      this._subs.push(this._editor.onChanged.subscribe(() => this._refreshErrors()));
      this._syncInputs();
    } catch (e: any) {
      this._loadError = `${e?.message ?? e}`;
      await mapServerError(e, {banner: (m) => this.setBanner(m), subject: this.table});
      this.setBanner(this._loadError);
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  /** (Re)builds every input from the registry properties, honoring the
   * {@link replaceInput} overrides. */
  private _buildInputs(): void {
    for (const sub of this._inputSubs)
      sub.unsubscribe();
    this._inputSubs = [];
    this._inputs.clear();
    this._refs.clear();
    this._principals.clear();
    this._host.innerHTML = '';
    const inputs: DG.InputBase[] = [];
    for (const p of this.editableProperties()) {
      const input = this._createInput(p);
      this._inputs.set(p.name, input);
      inputs.push(input);
    }
    if (inputs.length === 0) {
      this._host.appendChild(ui.divText(
        `You cannot ${this.isEditing ? 'edit' : 'create'} ${this.table} rows.`));
      return;
    }
    this._host.appendChild(ui.inputs(inputs));
  }

  private _createInput(p: DG.Property): DG.InputBase {
    const name = p.name;
    const factory = this._factories.get(name);
    let input: DG.InputBase;
    if (factory != null)
      input = factory(p);
    else if (isReferenceProperty(p)) {
      const target = `${p.semType}`;
      const ref = new RefInput({target: target, nullable: p.nullable === true,
        caption: p.friendlyName ?? name});
      this._refs.set(name, ref);
      // A user / group column is edited with the platform's own picker, exactly as the
      // Dart row editor does; the RefInput stays behind it as the display-text → id
      // resolver every machine write goes through (its own input is never mounted).
      if (target === 'User' || target === 'Group')
        input = this._principalInput(name, target === 'User', p);
      else {
        this._inputSubs.push(ref.onResolved.subscribe((resolved) => {
          if (this._syncing)
            return;
          this._refErrors.set(name, resolved.error);
          if (resolved.error == null)
            this._write(name, resolved.id);
          else
            this._refreshErrors();
        }));
        input = ref.input;
      }
    }
    else
      input = DG.InputBase.forProperty(p);
    input.addValidator(() => this._touched.has(name) ? this._errorOf(name) : null);
    if (this._refs.get(name) == null)
      this._inputSubs.push(input.onChanged.subscribe(() => {
        if (!this._syncing)
          this._write(name, input.value);
      }));
    return input;
  }

  /** The platform's user / group picker, single-select (`UserGroupSelector` on the
   * Dart side): the server-side lookup, the avatar tags and the clear affordance come
   * with it, and a selection writes the principal's id through the editor. */
  private _principalInput(name: string, isUser: boolean, p: DG.Property): DG.InputBase {
    const caption = p.friendlyName ?? name;
    const input: DG.InputBase = isUser ? ui.input.user(caption, {multiValue: false})
      : ui.input.userGroups(caption, {multiValue: false});
    this._principals.set(name, {input: input, isUser: isUser});
    this._inputSubs.push(input.onChanged.subscribe(() => {
      if (this._syncing)
        return;
      const selected: any[] = input.value ?? [];
      const id = principalId(selected, isUser);
      // A selection whose principal cannot be read is REPORTED, not written as empty —
      // the Dart editor refuses the save with the same message.
      this._refErrors.set(name, selected.length > 0 && id == null
        ? `Cannot resolve the selected ${isUser ? 'user' : 'group'}` : null);
      this._setResolverId(name, id);
      this._write(name, id);
    }));
    return input;
  }

  /** Points the picker at the id the row holds — the Dart editor's `_coreRefInput`
   * seeding. Generation-guarded, so a write made while a resolution is in flight wins;
   * a failure warns and leaves the picker empty instead of throwing. */
  private _seedPrincipal(name: string, id: string | null): void {
    const principal = this._principals.get(name)!;
    const generation = (this._seedGeneration.get(name) ?? 0) + 1;
    this._seedGeneration.set(name, generation);
    this._setResolverId(name, id);
    if (id == null) {
      this._setPrincipal(principal, []);
      return;
    }
    this._trackWork((principal.isUser ? grok.dapi.users.find(id) : grok.dapi.groups.find(id))
      .then((found: any) => {
        if (this._seedGeneration.get(name) !== generation)
          return;
        // The picker's items are personal GROUPs: a user resolves through its own.
        const item = found == null ? null : principal.isUser ? found.group ?? found : found;
        this._setPrincipal(principal, item == null ? [] : [item]);
      }, (e: any) => {
        grok.log.warning(`${this.table}.${name}: cannot resolve "${id}" — ${e?.message ?? e}`);
        if (this._seedGeneration.get(name) === generation)
          this._setPrincipal(principal, []);
      }));
  }

  private _setPrincipal(principal: PrincipalInput, items: any[]): void {
    // Restored rather than cleared: the seed of one field runs inside the sync of all
    // of them ({@link _syncInputs}), which must stay syncing for the fields after it.
    const syncing = this._syncing;
    this._syncing = true;
    try {
      principal.input.value = items;
    } finally {
      this._syncing = syncing;
    }
    principal.input.validate();
  }

  /** Keeps the resolver behind a picker pointing at the same id, so `getRefInput(name).id`
   * stays truthful; there is no display text to keep — its input is not mounted. */
  private _setResolverId(name: string, id: string | null): void {
    this._refs.get(name)?.setId(id, '', false);
  }

  private _write(name: string, value: any): void {
    this._touched.add(name);
    if (this._editor == null)
      this._pending[name] = value;
    else if (this._editor.dataFrame.columns.contains(name))
      this._editor.setValue(this._row, name, value);
    else
      // A column the row does not carry (not declared, or dropped by column
      // security): saying so beats writing into nothing, and beats throwing.
      this._customErrors.set(name, `${this.table} has no writable column '${name}'`);
    this._runValidators(name);
    this._refreshErrors();
  }

  /** A reference write: the display text (or an id) resolves through the input's own
   * suggestion source, and the outcome — an id or a validation message — lands on
   * the form. Never throws. */
  private _writeRef(name: string, value: any): void {
    const ref = this._refs.get(name)!;
    this._touched.add(name);
    const generation = (this._validatorGeneration.get(name) ?? 0) + 1;
    this._validatorGeneration.set(name, generation);
    this._trackWork(ref.resolve(value).then((resolved) => {
      if (this._validatorGeneration.get(name) !== generation)
        return;
      this._refErrors.set(name, resolved.error);
      if (this._principals.has(name))
        this._seedPrincipal(name, resolved.id);
      else {
        this._syncing = true;
        try {
          ref.setId(resolved.id, undefined, false);
        } finally {
          this._syncing = false;
        }
      }
      this._write(name, resolved.id);
      ref.input.validate();
    }));
  }

  private _runValidators(name: string): void {
    const validators = this._validators.get(name);
    if (validators == null || validators.length === 0) {
      this._customErrors.set(name, null);
      return;
    }
    const generation = (this._validatorGeneration.get(name) ?? 0) + 1;
    this._validatorGeneration.set(name, generation);
    const value = this.getValue(name);
    const pending: Promise<string | null>[] = [];
    let message: string | null = null;
    for (const validator of validators) {
      let result: any;
      try {
        result = validator(value);
      } catch (e: any) {
        result = `${e?.message ?? e}`;
      }
      if (result != null && typeof result.then === 'function')
        pending.push(result as Promise<string | null>);
      else if (message == null && result != null)
        message = `${result}`;
    }
    this._customErrors.set(name, message);
    if (pending.length === 0)
      return;
    this._trackWork(Promise.all(pending).then((results) => {
      if (this._validatorGeneration.get(name) !== generation)
        return;
      const async = results.find((r) => r != null);
      if (message == null && async != null) {
        this._customErrors.set(name, `${async}`);
        this._refreshErrors();
      }
    }, () => {}));
  }

  private _trackWork(work: Promise<any>): void {
    this._work.add(work);
    work.then(() => this._work.delete(work), () => this._work.delete(work));
  }

  /** A field's problem, most specific first: a custom validator, a reference that
   * did not resolve, then the editor's own cell error. */
  private _errorOf(name: string): string | null {
    const custom = this._customErrors.get(name);
    if (custom != null)
      return custom;
    const reference = this._refErrors.get(name);
    if (reference != null)
      return reference;
    // Untouched required fields of a NEW row are already marked in the editor; they
    // are a real blocker for save, and reported as such — the INPUTS only show a
    // marker once the field has been touched (see `_createInput`).
    return this._editor?.errorOf(this._row, name)?.message ?? null;
  }

  private _firstError(): string | null {
    for (const p of this.editableProperties()) {
      const error = this._errorOf(p.name);
      if (error != null)
        return `${p.friendlyName ?? p.name}: ${error}`;
    }
    return null;
  }

  private _refreshErrors(): void {
    for (const input of this._inputs.values())
      input.validate();
  }

  private _syncInputs(): void {
    this._syncing = true;
    try {
      for (const [name, input] of this._inputs) {
        const value = this.getValue(name);
        const ref = this._refs.get(name);
        if (this._principals.has(name))
          this._seedPrincipal(name, value == null ? null : `${value}`);
        else if (ref != null)
          ref.setId(value == null ? null : `${value}`, undefined, false);
        else
          input.value = value;
      }
    } finally {
      this._syncing = false;
    }
    this._refreshErrors();
  }

  private _syncInput(name: string): void {
    const input = this._inputs.get(name);
    if (input == null)
      return;
    this._syncing = true;
    try {
      input.value = this.getValue(name);
    } finally {
      this._syncing = false;
    }
    input.validate();
  }

  /** Shows a form-level message (a cross-field validation error, a duplicate, a
   * load failure); null hides the banner. */
  setBanner(message: string | null): void {
    this._banner.textContent = message ?? '';
    ui.setDisplay(this._banner, message != null && message !== '');
  }

  private _clearBanner(): void { this.setBanner(null); }

  private _propertyFor(p: DG.Property): DG.Property {
    const name = p.name;
    // The setter is ASSIGNED, not passed to create(): a setter handed to
    // `grok_Property` comes back from `Property.set` as a Dart closure JS cannot
    // call ('set is not a function' the first time `props[x] = v` runs), which is
    // why the platform's own `Widget.addProperty` assigns it too.
    const property = DG.Property.create(name, p.propertyType,
      () => this.getValue(name), null as any, null);
    property.set = (_: any, value: any) => this.setValue(name, value);
    const options: DG.IProperty = {friendlyName: p.friendlyName, nullable: p.nullable};
    if (p.description != null && p.description !== '')
      options.description = p.description;
    if (p.semType != null && p.semType !== '')
      options.semType = `${p.semType}`;
    if (p.choices != null && p.choices.length > 0)
      options.choices = p.choices;
    if (p.min != null)
      options.min = p.min;
    if (p.max != null)
      options.max = p.max;
    return property.fromOptions(options);
  }

  /** Natural language derived from the registry — what the AI reads to know what
   * this form is for. */
  private _describe(): string {
    const subject = this._context.info.singularName;
    if (!this.canWrite)
      return `A read-only ${subject} form (${this.table}): you cannot ` +
        `${this.isEditing ? 'edit' : 'create'} rows of this table.`;
    const properties = this.editableProperties();
    const fields = properties.map((p) => p.friendlyName ?? p.name).join(', ');
    const required = properties.filter((p) => p.nullable !== true).map((p) => p.friendlyName ?? p.name);
    const head = this.isEditing
      ? `Edits a ${subject} of ${this.table}.` : `Creates a ${subject} in ${this.table}.`;
    return `${head} Fields: ${fields}.` +
      (required.length === 0 ? '' : ` Required: ${required.join(', ')}.`);
  }
}
