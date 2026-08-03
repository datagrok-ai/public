/**
 * UI surface for entity-mapped domain tables: {@link DomainObjectHandler} — the
 * reflective, per-table {@link ObjectHandler} every domain table gets for free —
 * and the declarative types its members return.
 *
 * The data-plane counterpart is `src/domains.ts` + `grok.dapi.domains`
 * (see {@link DomainTableCapabilities}, `grok.dapi.domains.registry`).
 *
 * @module domains-ui
 */

import * as ui from '../ui';
import {EntityMetaDartProxy, ObjectHandler} from '../ui';
import {DomainRegistryClient, DomainTableClient} from './dapi';
import {DomainConditionTree, DomainTableCapabilities, splitDomainTable} from './domains';
import {DomainRow} from './entities/domain';
import {Property} from './entities/property';
import {SemanticValue} from './grid';
import {Grid} from './grid';
import {DataFrame} from './dataframe';
import {View} from './views/view';
import {IDartApi} from './api/grok_api.g';
import {toDart, toJs} from './wrappers';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

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
  private _meta: EntityMetaDartProxy | null | undefined;

  /** [table] addresses the domain table as `'<schema>.<table>'` (also the row
   * entity type and semType); a malformed address throws immediately. */
  constructor(public readonly table: string) {
    super();
    const [schemaName, tableName] = splitDomainTable(table);
    this.schemaName = schemaName;
    this.tableName = tableName;
  }

  get type(): string { return this.table; }

  get name(): string { return `${this.table} handler`; }

  /** Claims rows of {@link table}, including ones wrapped in a
   * {@link SemanticValue} (how the platform passes cell values around). */
  isApplicable(x: any): boolean {
    return this.rowOf(x) != null;
  }

  /** The {@link DomainRow} behind [x] (unwrapping a {@link SemanticValue}) when
   * it belongs to this table, null otherwise. Use it in overrides instead of
   * casting — handler members receive both shapes. */
  protected rowOf(x: any): DomainRow | null {
    const row = x instanceof SemanticValue ? x.value : x;
    return row instanceof DomainRow && row.typeName === this.table ? row : null;
  }

  /** The platform's per-table meta, as an {@link EntityMetaDartProxy} — what the
   * default render members delegate to. Resolved on first use and never a JS
   * handler (delegating into one would recurse); null only for a table whose
   * address cannot be parsed. */
  protected get dartMeta(): EntityMetaDartProxy | null {
    if (this._meta === undefined)
      this._meta = toJs(api.grok_DomainMeta_ForType(this.table)) ?? null;
    return this._meta!;
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
    return toJs(api.grok_DomainRow_Create(this.schemaName, this.tableName, null));
  }

  /** Permalink to the row's Entity View (business-key URL when unambiguous, id
   * URL otherwise — the platform's own rule). */
  deepLink(x: T): string | null {
    const row = this.rowOf(x);
    return row == null ? null : api.grok_DomainMeta_DeepLink(toDart(row));
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

  /** Decorates a grid over this table's rows: registry tags, reference captions,
   * `~`/system-column hiding, name-column-first ordering, semType renderer
   * re-resolution. Delegates to the platform meta, so it is identical to what
   * the built-in Domain View grid does — override it to take over completely
   * (see {@link ObjectHandler.renderGrid}). */
  renderGrid(grid: Grid, options?: {items?: DataFrame}): void {
    const m = this.dartMeta;
    if (m != null)
      m.renderGrid(grid, options);
  }

  /** Reflective property form over the writable columns of [x] (a new row when
   * omitted) — inputs come from {@link getProperties}, and non-writable columns
   * are excluded from the form AND from any payload built off it, mirroring
   * column security. A richer form (async validation, reference pickers, error
   * mapping) is `DomainForm` in `@datagrok-libraries/domain-ui`. */
  async renderEditor(x?: T): Promise<HTMLElement> {
    const [properties, caps] = await Promise.all([this.getProperties(), this.capabilities()]);
    const writable = caps.writableColumns;
    const inputs = properties.filter((p) => writable.includes(p.name));
    if (inputs.length === 0)
      return ui.divText(`You cannot edit ${this.table} rows.`);
    const row = (x == null ? null : this.rowOf(x)) ?? this.newRow();
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
      const smartFilter = `${c.fkColumn} = "${row.id}"`;
      return {
        name: ambiguous ? `${table} (${c.label})` : table,
        table: table,
        fkColumn: c.fkColumn,
        filter: [{property: c.fkColumn, operator: '=', value: row.id}] as DomainConditionTree,
        open: () => api.grok_Route(
          `/domains/${c.schema}/${c.table}?q=${encodeURIComponent(smartFilter)}`),
      };
    });
  }

  /** Actions available on [x] for the CURRENT user — the JS mirror of the Entity
   * View ribbon (Edit, Clone, Delete, Share, Watch, History, Copy link), gated by
   * server-truth row permissions ({@link DomainRow.permissions}) and the table's
   * registry metadata. Actions the user may not perform are absent from the list.
   *
   * The dialog-backed runs delegate to {@link editRow} / {@link cloneRow} /
   * {@link deleteRow} / {@link shareRow} / {@link showHistory}, which land with
   * the core openers (`DG.DomainView`) — override them until then. */
  async getRibbonActions(x: T): Promise<DomainAction[]> {
    const row = this.rowOf(x);
    if (row == null)
      return [];
    const info = await new DomainRegistryClient().tableInfo(this.table);
    const perms = await row.permissions();
    const res: DomainAction[] = [];
    if (perms.edit) {
      res.push({name: 'Edit...', icon: 'pencil', run: () => this.editRow(x)});
      res.push({name: 'Clone', icon: 'clone', run: () => this.cloneRow(x)});
    }
    if (perms.delete)
      res.push({name: 'Delete', icon: 'trash-alt', run: () => this.deleteRow(x)});
    if (info.securityMode === 'row')
      res.push({name: 'Share...', icon: 'share-alt', run: () => this.shareRow(x)});
    // Row-level watch needs the audit trail as its change source.
    if (info.audit && row.id != null) {
      const watching = await this.client.isWatching(row.id);
      res.push({name: watching ? 'Unwatch' : 'Watch', icon: 'bell',
        run: () => this.setWatch(x, !watching)});
    }
    res.push({name: 'History', icon: 'history', run: () => this.showHistory(x)});
    res.push({name: 'Copy link', icon: 'link', run: () => this.copyLink(x)});
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

  /** Copies the row's {@link deepLink} to the clipboard. */
  async copyLink(x: T): Promise<void> {
    const link = this.deepLink(x);
    if (link != null)
      await navigator.clipboard.writeText(`${window.location.origin}${link}`);
  }

  /** Opens the platform's create/edit dialog; resolves to whether a row was saved. */
  editRow(x?: T): Promise<boolean> { return this._opener('editRow'); }

  /** Opens the platform's create dialog prefilled from [x]. */
  cloneRow(x: T): Promise<boolean> { return this._opener('cloneRow'); }

  /** Deletes [x] after the standard confirmation. */
  deleteRow(x: T): Promise<boolean> { return this._opener('deleteRow'); }

  /** Opens the standard sharing flow for [x] (row-mode tables). */
  shareRow(x: T): Promise<boolean> { return this._opener('shareRow'); }

  /** Shows the row's audit history. */
  showHistory(x: T): Promise<void> { return this._opener('showHistory'); }

  /** Opens the platform's row picker for this table; resolves to the picked row
   * or null. */
  pickRow(): Promise<DomainRow | null> { return this._opener('pickRow'); }

  /** The dialog openers land with `DG.DomainView` (ui-js-api WO-5); until then a
   * subclass must provide its own. Rejecting (rather than silently doing nothing)
   * keeps a half-wired action visible in the console instead of dead in the UI. */
  private _opener(member: string): Promise<never> {
    return Promise.reject(new Error(
      `DomainObjectHandler.${member}() is not available yet — the core dialog openers ship with ` +
      `DG.DomainView. Override it, or navigate to /domains/${this.schemaName}/${this.tableName}.`));
  }
}
