/**
 * UI surface for entity-mapped domain tables: {@link DomainObjectHandler} — the
 * reflective, per-table {@link ObjectHandler} every domain table gets for free,
 * whose statics also open the platform's row dialogs, picker, conflict dialog
 * and audit/grants panes — plus {@link DomainView}, the built-in Domain View as
 * an embeddable, programmatically created view.
 *
 * The data-plane counterpart is `src/domains.ts` + `grok.dapi.domains`
 * (see {@link DomainTableCapabilities}, `grok.dapi.domains.registry`).
 *
 * @module domains-ui
 */

import * as ui from '../ui';
import {EntityMetaDartProxy, ObjectHandler} from '../ui';
import {DomainRegistryClient, DomainTableClient} from './dapi';
import {domainCall, DomainConditionTree, DomainQueryParams, DomainTableCapabilities,
  splitDomainTable} from './domains';
import {DomainRow} from './entities/domain';
import {Property} from './entities/property';
import {SemanticValue} from './grid';
import {View} from './views/view';
import {CardView} from './views/card_view';
import {IDartApi} from './api/grok_api.g';
import {toDart, toJs} from './wrappers';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** Whether [x] is a plain values map — the shape ONE row of a `query()` result
 * has. A string, an array or a wrapped platform object would be read key by key
 * into nonsense values (`'abc'` becomes `{0: 'a', 1: 'b', 2: 'c'}`) instead of
 * failing, so {@link DomainObjectHandler.rowFrom} rejects them. */
function _isValuesMap(x: any): boolean {
  if (x == null || typeof x !== 'object' || Array.isArray(x))
    return false;
  const proto = Object.getPrototypeOf(x);
  return proto === null || proto === Object.prototype;
}

/** One FK-inverted child ("detail") table of a domain row — what the built-in
 * Entity View shows as a child tab (see {@link DomainObjectHandler.getDetailTabs}). */
export interface DomainDetailTab {
  /** Tab caption: the child table address, suffixed with the FK label when the
   * child references the parent through more than one column. */
  name: string;
  /** Child table address, `'<schema>.<table>'`. */
  table: string;
  /** The child table's FK column referencing the parent row. */
  fkColumn: string;
  /** Condition tree scoping the child table to the parent row — pass it as the
   * `filter` of `grok.dapi.domains.table(tab.table).query/queryDf`. */
  filter: DomainConditionTree;
  /** Opens the child table's Domain View scoped to the parent row. */
  open(): void;
}

/** A declarative action offered for a domain row (see
 * {@link DomainObjectHandler.getRibbonActions}) — the JS mirror of the Entity
 * View ribbon icons. Only actions the current user may perform are returned, so
 * a consumer renders the list as-is without any permission checks of its own. */
export interface DomainAction {
  name: string;
  /** FontAwesome icon name, matching the platform ribbon ('pencil', 'trash-alt'...). */
  icon?: string;
  /** Whether running the action may have changed the row — what a list or a page
   * reloads on. The platform's own actions ({@link
   * DomainObjectHandler.getRibbonActions}) all declare it; a custom action that
   * omits it is up to its consumer to classify. */
  changesRow?: boolean;
  run(): Promise<any> | any;
}

/**
 * The per-table {@link ObjectHandler} for domain rows: reflective (it takes
 * columns, labels, choices and capabilities from the runtime registry, so it
 * works on ANY table without codegen) and delegating (every render member falls
 * through to the platform's Dart meta for that table).
 *
 * ```ts
 * class IssueHandler extends DG.DomainObjectHandler {
 *   constructor() { super('grit.issue'); }
 *   renderCard(x) { return ui.divText(`#${x.values.number} ${x.values.title}`); }
 * }
 * DG.ObjectHandler.register(new IssueHandler());
 * ```
 *
 * Dispatch precedence (unchanged by this class):
 * 1. **Rendering** — the last-registered JS handler applicable to the object
 *    wins `ObjectHandler.forEntity`; the platform's per-table meta demotes
 *    itself behind it regardless of registration order.
 * 2. **Overriding nothing is never a regression** — an override-nothing subclass
 *    renders exactly like the platform, because the defaults below delegate to
 *    the same Dart meta the platform would have used.
 * 3. **Platform commands stay** — Edit/Delete/Share/Watch/History context-menu
 *    commands remain the platform's, permission-gated; a JS handler ADDS actions
 *    (`registerParamFunc`), it never replaces those.
 * 4. **Grid cells** keyed by semType keep their registered cell renderer.
 */
export class DomainObjectHandler<T = DomainRow> extends ObjectHandler<T> {
  /** Domain schema name of {@link table}. */
  readonly schemaName: string;
  /** Table name of {@link table} within {@link schemaName}. */
  readonly tableName: string;
  private _meta: EntityMetaDartProxy | null = null;

  /** [table] addresses the domain table as `'<schema>.<table>'` (also the row
   * entity type and semType); a malformed address throws immediately. */
  constructor(public readonly table: string) {
    super();
    const [schemaName, tableName] = splitDomainTable(table);
    this.schemaName = schemaName;
    this.tableName = tableName;
  }

  get type(): string { return this.table; }

  /** Claims rows of {@link table}, including ones wrapped in a
   * {@link SemanticValue} (how the platform passes cell values around). */
  isApplicable(x: any): boolean {
    return this.rowOf(x) != null;
  }

  /** The {@link DomainRow} behind [x] (unwrapping a {@link SemanticValue}) when
   * it belongs to this table, null otherwise. Use it in overrides instead of
   * casting — handler members receive both shapes. */
  protected rowOf(x: any): DomainRow | null {
    // SemanticValue.value is a sync interop getter: what it returns for an entity
    // is the Dart handle, so it takes a toJs() to become the row again.
    const row = x instanceof SemanticValue ? toJs(x.value) : x;
    return row instanceof DomainRow && row.typeName === this.table ? row : null;
  }

  /** {@link rowOf}, for the members that cannot act without a row: throws a
   * named error naming the table instead of failing later on a null. */
  protected rowOrThrow(x: any): DomainRow {
    const row = this.rowOf(x);
    if (row == null)
      throw new Error(`not a ${this.table} row: ${DomainObjectHandler._describe(x)}`);
    return row;
  }

  /** A rejected argument, readably: a plain object stringifies to
   * '[object Object]', which names neither the value nor its type. A
   * {@link SemanticValue} is described by what it wraps — 'SemanticValue' names
   * the platform's envelope, not the row the caller passed. */
  private static _describe(x: any): string {
    if (x instanceof SemanticValue)
      return DomainObjectHandler._describe(toJs(x.value));
    if (x == null || typeof x !== 'object')
      return `${x}`;
    if (x instanceof DomainRow)
      return `a ${x.typeName} row`;
    const type = x.constructor?.name;
    if (type != null && type !== 'Object')
      return type;
    try {
      return JSON.stringify(x);
    } catch (_) {
      return `${x}`;
    }
  }

  /** The platform's per-table meta, as an {@link EntityMetaDartProxy} — what the
   * default render members delegate to. Never a JS handler (delegating into one
   * would recurse); null only for a table whose address cannot be parsed. Only a
   * real hit is cached: a miss is re-resolved on the next use, so a handler
   * constructed before the platform metas exist starts delegating as soon as
   * they do. */
  protected get dartMeta(): EntityMetaDartProxy | null {
    if (this._meta == null)
      this._meta = toJs(api.grok_DomainMeta_ForType(this.table)) ?? null;
    return this._meta;
  }

  private _delegate(x: T, dart: (m: EntityMetaDartProxy, row: DomainRow) => HTMLElement,
                    fallback: () => HTMLElement): HTMLElement {
    const m = this.dartMeta;
    const row = this.rowOf(x);
    return m == null || row == null ? fallback() : dart(m, row);
  }

  /** Client for this table (`grok.dapi.domains.table(this.table)`). */
  protected get client(): DomainTableClient {
    return new DomainTableClient(api.grok_Dapi_Domains(), this.schemaName, this.tableName);
  }

  // ─────────────────────── reflective metadata ─────────────────────────

  /** Runtime {@link Property} metadata of the table's declared columns (type,
   * semType, choices, min/max, nullable, friendly label), bound to
   * {@link DomainRow} values. Rejects with a `DomainValidationError` for a table
   * that is not registered. */
  getProperties(): Promise<Property[]> {
    return new DomainRegistryClient().rowProperties(this.table);
  }

  /** Effective {@link DomainTableCapabilities} of the current user on this table.
   * Every affordance below is derived from it — consumers never hand-wire
   * permission checks. */
  capabilities(): Promise<DomainTableCapabilities> {
    return this.client.capabilities();
  }

  /** Loads the row by id (the acquisition path for a {@link DomainRow} in JS);
   * null when it does not exist or is not visible. */
  async getById(id: string): Promise<T | null> {
    const m = this.dartMeta;
    return m == null ? null : await api.grok_Meta_GetById(m.dart, id);
  }

  /** A new, unsaved row of this table — what {@link renderEditor} binds to when
   * called without an object. Insert it with
   * `grok.dapi.domains.table(...).insert(row.values)`. */
  newRow(): DomainRow {
    return this.rowFrom(null);
  }

  /** Wraps [values] — one row of a `query()` result — as a {@link DomainRow} of
   * this table, locally and WITHOUT a round trip: the acquisition path a list, a
   * card or a context panel uses to render, address and act on rows it already
   * fetched ({@link getById} is the one that asks the server). Identity
   * (`displayName`, `semValue`, deep links) resolves from the registry, so a row
   * fetched without its business-key columns addresses itself by id. Anything
   * but a plain values map (or null, for a new row) is REJECTED by name: read
   * key by key, a string or an array produces a row of nonsense values that only
   * shows up much later, as a card with no identity. */
  rowFrom(values: {[key: string]: any} | null): DomainRow {
    if (values != null && !_isValuesMap(values))
      throw new Error(`${this.table}: rowFrom() takes one row of a query() result, not ` +
        DomainObjectHandler._describe(values));
    return toJs(api.grok_DomainRow_Create(this.schemaName, this.tableName, values));
  }

  /** Permalink to the row's Entity View (business-key URL when unambiguous, id
   * URL otherwise — the platform's own rule); null for an unsaved row
   * ({@link newRow}), which has no address yet. */
  deepLink(x: T): string | null {
    return DomainObjectHandler.deepLink(this.rowOrThrow(x));
  }

  /** Opens the row's Entity View — the platform's default (double-click) action
   * for a domain row. */
  openRow(x: T): void {
    DomainObjectHandler.openRow(this.rowOrThrow(x));
  }

  // ─────────────────────── rendering (delegating defaults) ─────────────────────────

  getCaption(x: T): string {
    const m = this.dartMeta;
    const row = this.rowOf(x);
    return m == null || row == null ? super.getCaption(x) : m.getCaption(row);
  }

  renderIcon(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderIcon(r, context), () => super.renderIcon(x, context));
  }

  renderMarkup(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderMarkup(r, context), () => super.renderMarkup(x, context));
  }

  renderTooltip(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderTooltip(r, context), () => super.renderTooltip(x, context));
  }

  renderCard(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderCard(r, context), () => super.renderCard(x, context));
  }

  renderProperties(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderProperties(r, context), () => super.renderProperties(x, context));
  }

  renderView(x: T, context: any = null): HTMLElement {
    return this._delegate(x, (m, r) => m.renderView(r, context), () => super.renderView(x, context));
  }

  async renderPreview(x: T, params?: any, path?: string): Promise<View> {
    const m = this.dartMeta;
    const row = this.rowOf(x);
    return m == null || row == null ? super.renderPreview(x, params, path) : m.renderPreview(row, params, path);
  }

  // renderGrid is deliberately NOT overridden: the base
  // {@link ObjectHandler.renderGrid} already delegates to the platform meta for
  // its type, so this handler decorates a grid exactly like the Domain View —
  // and so does every other non-overriding handler.

  /** Reflective property form over the writable columns of [x] (a new row when
   * omitted) — inputs come from {@link getProperties}, and non-writable columns
   * are excluded from the form AND from any payload built off it, mirroring
   * column security. Callers without the corresponding table capability
   * (`canInsert` for a new row, `canEdit` for an existing one) get a read-only
   * explanation instead: column security alone is NOT table Edit. A richer form
   * (async validation, reference pickers, error mapping) is `DomainForm` in
   * `@datagrok-libraries/domain-ui`. */
  async renderEditor(x?: T): Promise<HTMLElement> {
    const [properties, caps] = await Promise.all([this.getProperties(), this.capabilities()]);
    const writable = caps.writableColumns;
    const inputs = properties.filter((p) => writable.includes(p.name));
    const row = (x == null ? null : this.rowOf(x)) ?? this.newRow();
    const creating = row.id == null;
    if (!(creating ? caps.canInsert : caps.canEdit) || inputs.length === 0)
      return ui.divText(creating ? `You cannot create ${this.table} rows.`
        : `You cannot edit ${this.table} rows.`);
    return ui.input.form(toDart(row), inputs);
  }

  // ─────────────────────── detail tables and actions ─────────────────────────

  /** FK-inverted child tables scoped to [x] — the child tabs the built-in Entity
   * View shows. Empty for an unsaved row or a table nothing references. */
  async getDetailTabs(x: T): Promise<DomainDetailTab[]> {
    const row = this.rowOf(x);
    if (row == null || row.id == null)
      return [];
    const children = (await new DomainRegistryClient().tableInfo(this.table)).childTables;
    return children.map((c) => {
      const table = `${c.schema}.${c.table}`;
      // Two FKs to the same parent yield two tabs for one child table —
      // disambiguate by the FK label (the platform's own rule).
      const ambiguous = children.filter((o) => `${o.schema}.${o.table}` === table).length > 1;
      return {
        name: ambiguous ? `${table} (${c.label})` : table,
        table: table,
        fkColumn: c.fkColumn,
        filter: [{property: c.fkColumn, operator: '=', value: row.id}] as DomainConditionTree,
        // The path (and its `q` smart filter) comes from the platform's own
        // builder, and opens through the platform's own navigation — the
        // /domains URL scheme has ONE encode side and ONE open side.
        open: () => DomainObjectHandler.openPath(
          api.grok_DomainMeta_ChildTablePath(c.schema, c.table, c.fkColumn, row.id!)),
      };
    });
  }

  /** Actions available on [x] for the CURRENT user — the JS mirror of the Entity
   * View ribbon plus its default Open action (Open, Edit, Clone, Delete, Share,
   * Watch, History, Copy link), gated by server-truth row permissions
   * ({@link DomainRow.permissions}) and the table's registry metadata. Actions
   * the user may not perform are absent from the list; an unsaved row
   * ({@link newRow}) has none of them — it has no address, history or
   * permissions until it is inserted.
   *
   * The dialog-backed runs delegate to {@link editRow} / {@link cloneRow} /
   * {@link deleteRow} / {@link shareRow} / {@link showHistory} — override one of
   * those to replace a flow, the action list stays the platform's.
   *
   * Every action declares {@link DomainAction.changesRow}: the platform knows
   * which of its own actions write, so a consumer reloads on it instead of
   * guessing from the caption (which changes with the state: 'Watch' ⇄
   * 'Unwatch') or from the icon. */
  async getRibbonActions(x: T): Promise<DomainAction[]> {
    const row = this.rowOf(x);
    if (row?.id == null)
      return [];
    const info = await new DomainRegistryClient().tableInfo(this.table);
    const perms = await row.permissions();
    const res: DomainAction[] = [
      {name: 'Open', icon: 'folder-open', changesRow: false, run: () => this.openRow(x)}];
    if (perms.edit) {
      res.push({name: 'Edit...', icon: 'pencil', changesRow: true, run: () => this.editRow(x)});
      res.push({name: 'Clone', icon: 'clone', changesRow: true, run: () => this.cloneRow(x)});
    }
    if (perms.delete)
      res.push({name: 'Delete', icon: 'trash-alt', changesRow: true, run: () => this.deleteRow(x)});
    // Row-mode tables are the ones that CAN share a row; Share on the row itself
    // is what makes it offerable (the Dart ribbon shows it and explains the
    // denial on click — JS consumers render the list as-is, so gate it here).
    if (info.securityMode === 'row' && perms.share)
      res.push({name: 'Share...', icon: 'share-alt', changesRow: true, run: () => this.shareRow(x)});
    // Row-level watch needs the audit trail as its change source. NB the state
    // is a server round trip per build (the Dart ribbon reads its warmed cache);
    // callers rebuilding the list often should cache it themselves.
    if (info.audit) {
      const watching = await this.client.isWatching(row.id);
      res.push({name: watching ? 'Unwatch' : 'Watch', icon: 'bell', changesRow: false,
        run: () => this.setWatch(x, !watching)});
    }
    res.push({name: 'History', icon: 'history', changesRow: false, run: () => this.showHistory(x)});
    res.push({name: 'Copy link', icon: 'link', changesRow: false, run: () => this.copyLink(x)});
    return res;
  }

  /** Subscribes/unsubscribes the current user to changes of [x]; resolves to the
   * resulting state (the server refuses row watch on tables without an audit
   * trail, and the state then stays unchanged). */
  async setWatch(x: T, value: boolean): Promise<boolean> {
    const row = this.rowOf(x);
    if (row?.id == null)
      return false;
    return value ? await this.client.watch(row.id)
      : !(await this.client.unwatch(row.id));
  }

  /** Copies the row's {@link deepLink} to the clipboard (with the platform's
   * confirmation balloon). */
  copyLink(x: T): void {
    DomainObjectHandler.copyLink(this.rowOrThrow(x));
  }

  /** Opens the platform's edit dialog for [x] — its create dialog when [x] is
   * omitted or unsaved; resolves to whether a row was saved. */
  editRow(x?: T): Promise<boolean> {
    const row = x == null ? null : this.rowOf(x);
    return row?.id == null ? DomainObjectHandler.createRow(this.table) : DomainObjectHandler.editRow(row);
  }

  /** Opens the platform's create dialog prefilled from [x]. */
  cloneRow(x: T): Promise<boolean> { return DomainObjectHandler.cloneRow(this.rowOrThrow(x)); }

  /** Deletes [x] after the standard confirmation. */
  deleteRow(x: T): Promise<boolean> { return DomainObjectHandler.deleteRow(this.rowOrThrow(x)); }

  /** Opens the standard sharing flow for [x] (row-mode tables); resolves when
   * the sharing dialog is SHOWN, not when sharing is done. */
  shareRow(x: T): Promise<void> { return DomainObjectHandler.shareRow(this.rowOrThrow(x)); }

  /** Shows the row's audit history in the platform's History dialog. */
  showHistory(x: T): void { DomainObjectHandler.showHistory(this.rowOrThrow(x)); }

  /** The row's audit history as an embeddable element (the History pane of the
   * built-in context panel). */
  auditPane(x: T): HTMLElement { return DomainObjectHandler.auditPane(this.rowOrThrow(x)); }

  /** Opens the platform's row picker for this table; resolves to the picked row
   * or null. */
  pickRow(): Promise<DomainRow | null> { return DomainObjectHandler.pickRow(this.table); }

  // ─────────────────────── openers (platform dialogs and panes) ─────────────────────────
  // Canonical as STATICS: each one enters the SAME Dart flow the built-in
  // surfaces use, so there is one implementation of every dialog, and no
  // handler instance is needed to open one (a picker for a table nobody wrote a
  // handler for, a conflict dialog inside a save loop). The instance members
  // above delegate here — a subclass customizing a flow overrides ONE of them.

  /** Opens the platform's create dialog for [table] (`'<schema>.<table>'`);
   * resolves to whether a row was saved. Inputs come from the registry, values
   * are validated with the same code the server re-runs, columns the caller
   * cannot write are absent from the form AND the payload. */
  static createRow(table: string): Promise<boolean> {
    const [schema, name] = splitDomainTable(table);
    return domainCall(api.grok_DomainRowEditor_OpenCreate(schema, name));
  }

  /** Opens the platform's edit dialog for [row]; resolves to whether it was
   * saved. A version conflict goes through the standard reload/overwrite dialog
   * — reloading and saving in the reopened dialog still resolves true. */
  static editRow(row: DomainRow): Promise<boolean> {
    return domainCall(api.grok_DomainRowEditor_OpenEdit(toDart(row)));
  }

  /** Opens the create dialog prefilled from [row] (writable columns only). */
  static cloneRow(row: DomainRow): Promise<boolean> {
    return domainCall(api.grok_DomainRowEditor_OpenClone(toDart(row)));
  }

  /** Deletes [row] after the platform's confirmation; resolves to whether it was
   * deleted — false on cancel, and on a failed delete (a restrict violation, no
   * Delete permission, a server error), which the platform reports in a balloon
   * and then closes the dialog. */
  static deleteRow(row: DomainRow): Promise<boolean> {
    return domainCall(api.grok_DomainMeta_ConfirmDelete(toDart(row)));
  }

  /** Opens the standard sharing flow for [row]: the row is promoted to an entity
   * (idempotent, needs Share on it) and the platform sharing dialog follows.
   * Row-mode tables only — elsewhere access comes from the securing table, whose
   * grants are {@link grantsPane}'s subject. Resolves once the dialog is SHOWN
   * (promotion done, permission denial already ballooned) — NOT when the user
   * finishes sharing; subscribe to the platform's sharing events for that. */
  static shareRow(row: DomainRow): Promise<void> {
    return domainCall(api.grok_DomainMeta_ShareRow(toDart(row)));
  }

  /** Opens the platform's lookup picker for [table] (`'<schema>.<table>'`) — the
   * target table's Domain View in a dialog, with search and single select;
   * resolves to the picked row, or null on cancel. */
  static pickRow(table: string): Promise<DomainRow | null> {
    const [schema, name] = splitDomainTable(table);
    return domainCall(api.grok_DomainRowPicker_Pick(schema, name));
  }

  /** The standard optimistic-concurrency dialog for a 409, naming [subject] (the
   * row's display value). Resolves to the user's decision — `'reload'` (discard
   * my changes, take the server's), `'overwrite'` (retry against the current
   * version), or null when the dialog is dismissed. The caller applies it; the
   * platform's own editors resolve their conflicts through this same dialog. */
  static showConflictDialog(subject: string): Promise<'reload' | 'overwrite' | null> {
    return domainCall(api.grok_DomainConflictDialog_Show(subject));
  }

  /** The row's audit trail as an embeddable element: who / when / what, newest
   * first, with the field diff in each row's tooltip and a link opening the full
   * event list as a table. Explains itself when the table has no audit trail. */
  static auditPane(row: DomainRow): HTMLElement {
    return api.grok_Domains_AuditPane(toDart(row));
  }

  /** {@link auditPane} in the platform's History dialog. */
  static showHistory(row: DomainRow): void {
    api.grok_DomainMeta_ShowHistory(toDart(row));
  }

  /** Read-only grants pane of a REGISTRY entity — a schema, table or property
   * schema id (`DomainTable.id`, `grok.dapi.domains.registry`), NOT a domain row:
   * rows are shared through {@link shareRow}. Lists the direct permission rows
   * plus a 'Manage access...' link ([options.readOnly] drops the link); listing
   * requires Share, and non-Share callers get an explanation instead. */
  static grantsPane(entityId: string, name?: string, options?: {readOnly?: boolean}): HTMLElement {
    return api.grok_Domains_GrantsPane(entityId, name ?? null, options?.readOnly === true);
  }

  /** The editable grants dialog for a registry entity (grant/revoke per group);
   * resolves when it is shown. */
  static showGrantsDialog(entityId: string, name?: string): Promise<void> {
    return api.grok_Domains_GrantsDialog(entityId, name ?? null);
  }

  /** Permalink to [row]'s Entity View — the business-key URL when unambiguous,
   * the id URL otherwise; null for an unsaved row, which has no address yet. */
  static deepLink(row: DomainRow | null): string | null {
    return row?.id == null ? null : api.grok_DomainMeta_DeepLink(toDart(row));
  }

  /** Opens [row]'s Entity View — the platform's default (double-click) action. */
  static openRow(row: DomainRow | null): void {
    const link = DomainObjectHandler.deepLink(row);
    if (link != null)
      DomainObjectHandler.openPath(link);
  }

  /** Opens a `/domains` address through the platform's own navigation: row deep
   * links, table views, and the scoped child-table views {@link getDetailTabs}
   * produces. THE navigation entry point for domain addresses — do not re-spell
   * it with `grok.shell.route`, which pushes history differently. */
  static openPath(path: string): void {
    api.grok_DomainMeta_Open(path);
  }

  /** Copies [row]'s {@link deepLink} to the clipboard, with the platform's
   * confirmation balloon; a no-op for an unsaved row. */
  static copyLink(row: DomainRow | null): void {
    if (DomainObjectHandler.deepLink(row) != null)
      api.grok_DomainMeta_CopyLink(toDart(row));
  }
}


/** Options of {@link DomainView.create}. */
export interface DomainViewOptions {
  /** Domain schema name. */
  schema: string;
  /** Table name within {@link schema}. */
  table: string;
  /** Invisible {@link https://datagrok.ai/help/datagrok/smart-search | smart filter}
   * AND-combined with everything the user does — the way to scope an embedded
   * view to a subset (`'project_id = "..."'`). */
  permanentFilter?: string;
  /** Hosted inside another view: brief list, no breadcrumbs, no export ribbon. */
  embedded?: boolean;
}


/**
 * The platform's Domain View for one domain table — the `/domains/<schema>/<table>`
 * page, creatable programmatically and embeddable: card / brief / grid modes,
 * search, sort, paging, the server-backed filter panel, in-grid editing with
 * one-transaction save, and the per-table {@link ObjectHandler}'s rendering and
 * commands. Show it with `grok.shell.addView(v)`, or mount `v.root` anywhere.
 *
 * ```ts
 * const v = DG.DomainView.create({schema: 'grit', table: 'issue'});
 * grok.shell.addView(v);
 * v.showFilters();
 * ```
 *
 * Rows are capped at 1000 (a banner invites refining the filter), and
 * {@link query} reports that cap as its `limit` — running it reproduces exactly
 * what the view shows; raise or drop `limit` on the reported params to go past
 * the cap. The inherited {@link CardView} members that need a JS-defined data
 * source — `meta`, `objectType`, `searchFields` — do not apply here: the view
 * takes them from the domain registry, and `meta` throws to say so.
 *
 * A view the USER opened by navigating to `/domains/<schema>/<table>` is a plain
 * `DG.View` in JS (route-opened views are, platform-wide) — wrap its handle to
 * get this surface: `new DG.DomainView(grok.shell.v.dart)`.
 */
export class DomainView extends CardView {
  declare dart: any;

  constructor(dart: any) {
    super(dart);
  }

  /** Creates the Domain View for a table; the rows load asynchronously.
   * Deliberately independent of the `enableDomainDatabases` client flag that
   * gates the `/domains` pages: calling this API IS the opt-in. */
  static create(options: DomainViewOptions): DomainView {
    return new DomainView(api.grok_DomainView_Create(options));
  }

  /** @throws always — the Domain View takes its object handler from the domain
   * registry (the per-table `ObjectHandler`), not from a JS-assigned one. */
  get meta(): ObjectHandler { throw new Error(DomainView._noMeta); }
  set meta(_: ObjectHandler) { throw new Error(DomainView._noMeta); }

  private static readonly _noMeta =
    'the Domain View takes its meta from the domain registry: register an ' +
    'ObjectHandler (DG.DomainObjectHandler) for the table instead';

  /** Docks the filter panel (server-backed facets and histograms) — the 'Toggle
   * filters' icon; idempotent. Needs the view to be docked. */
  showFilters(): void {
    api.grok_DomainView_ShowFilters(this.dart);
  }

  /** What the user is looking at right now, as `DomainQuery` parameters: the
   * panel facets AND the search string merged into filter elements, the current
   * ordering, and the view's row cap. Serializable — the shape a deep link, a
   * saved filter or an "open in Table View" carries. */
  get query(): DomainQueryParams {
    return toJs(api.grok_DomainView_Get_Query(this.dart));
  }
}
