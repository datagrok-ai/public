/**
 * {@link DomainGrid} — an editable grid over a domain table: the platform's own
 * grid decoration, batch editing through a {@link DomainFrameEditor}, and a
 * save/cancel bar, all permission-gated — plus the widget-action Func vocabulary
 * every domain widget reports through `getFunctions()`.
 *
 * That vocabulary is the domain widgets' actions as REAL platform {@link Func}s
 * ({@link DomainGrid.saveFunc} and friends), so `widget.getFunctions()` returns
 * something typed, annotated, discoverable and invocable through the machinery
 * the platform already has (and already exposes to the AI integrations,
 * `src/widgets/base.ts`). ONE shared vocabulary, not one function per widget
 * instance: `Save`, `Discard` and `Reset` each take the widget they act on
 * (`//input: widget widget`) and call the member of the same name on it, so a
 * hundred open forms register nothing and every one of them reports the same
 * three names. Registration happens on first use and is idempotent — the Func is
 * cached per qualified name.
 *
 * ```ts
 * const save = form.getFunctions().find((f) => f.name === 'Save')!;
 * await save.apply({widget: form});
 * ```
 *
 * @module domains-grid
 */

import * as ui from '../../../ui';
import {EntityMetaDartProxy, ObjectHandler} from '../../../ui';
import {DomainTableClient} from '../../dapi';
import {DataFrame} from '../../dataframe';
// DomainTableCapabilities is referenced by the gating docs below, nowhere in the code.
import {DomainQuerySpec, DomainTableCapabilities} from '../../domains';
import {DomainAction, DomainObjectHandler} from '../../domains-ui';
import {Func} from '../../entities/func';
import {Functions} from '../../functions';
import {Grid, GridCell} from '../../grid';
import {Logger} from '../../logger';
import {IWidgetStatus, Widget} from '../../widgets/base';
import {Balloon} from '../../widgets/menu';
import {acquireDomainContext, DomainFrameEditor, DomainFrameEditorOptions,
  IDomainTableContext, IEditorHost} from './domains-editor';

/** The `grok.shell` / `grok.log` / `grok.functions` surface the code below uses,
 * through the internal classes: no js-api module imports `../grok`. */
const balloon = new Balloon();
const log = Logger.getStatic();
const functions = new Functions();

/** The widget-action Funcs registered so far, by name — what keeps the
 * vocabulary ONE registration per session. */
const _funcs = new Map<string, Func>();

/** The grid-scoped slice of the `domain-ui-*` stylesheet, injected once at
 * runtime rather than imported as a `.css` file (the node bundle null-loads css).
 * The class names stay `domain-ui-*`: they are the ones plugin overrides target.
 * Colours are design tokens only — the in-grid CELL colours below are canvas
 * ints, which are not CSS-token territory. */
const STYLE_ID = 'domain-grid-styles';

const DOMAIN_GRID_CSS = `
.domain-ui-grid {
  height: 100%;
}

/* The grid host inside it: a canvas grid has no intrinsic size, so it needs the
   space the toolbar leaves — in BOTH directions. */
.domain-ui-grid > .ui-box {
  flex-grow: 1;
  min-height: 0;
  width: 100%;
}

.domain-ui-grid-toolbar {
  flex: 0 0 auto;
  align-items: center;
  gap: 4px;
  padding: 4px 8px;
  background-color: var(--grey-1);
  border-bottom: 1px solid var(--grey-2);
}

.domain-ui-save-bar {
  align-items: center;
  gap: 4px;
  margin-left: auto;
}

.domain-ui-edit-count {
  color: var(--text-color-light);
  font-size: 12px;
  white-space: nowrap;
  padding-right: 4px;
}

.domain-ui-grid-toolbar button[disabled] {
  color: var(--grey-3);
  cursor: default;
}
`;

/** Adds the grid stylesheet to the document once — safe to call per component. */
export function applyDomainGridStyles(): void {
  if (typeof document === 'undefined' || document.getElementById(STYLE_ID) != null)
    return;
  const style = document.createElement('style');
  style.id = STYLE_ID;
  style.textContent = DOMAIN_GRID_CSS;
  document.head.appendChild(style);
}

/** Background of a cell with an unsaved edit. Canvas ints, matching the
 * platform's own in-grid editing (compare-tables convention) — cell backgrounds
 * are painted on the canvas, so they are not CSS-token territory. */
const DIRTY_CELL_COLOR = 0xFFFFF3CD;    // soft amber: unsaved edit
const INVALID_CELL_COLOR = 0xFFFFB3B0;  // soft red: invalid / rejected
const CONFLICT_CELL_COLOR = 0xFFFFE0B2; // soft orange: dismissed version conflict

/** Options of {@link DomainGrid.create}. */
export interface DomainGridOptions extends DomainFrameEditorOptions {
  /** Master switch for editing; ANDed with the table's `canEdit` capability, so
   * a table the user may not edit is read-only whatever this says. */
  editable?: boolean;
  /** Show the toolbar (unsaved-change count, Save/Cancel, Add/Delete row).
   * Defaults to true; turn it off to drive the editor from your own ribbon. */
  toolbar?: boolean;
  /** Values every row added through the toolbar starts with — how a detail grid
   * of a master-detail page prefills the foreign key of its parent, so 'Add row'
   * produces a saveable row instead of one failing on a required column. */
  defaults?: {[column: string]: any};
}

/**
 * An editable grid over one domain table.
 *
 * ```ts
 * const grid = await DG.DomainGrid.create(grok.dapi.domains.table('grit.issue'),
 *   {query: {filter: 'status = "open"'}});
 * grok.shell.newView('Issues', [grid.root]);
 * ```
 *
 * {@link DomainGrid.create} is the one-await form from a bare client; the
 * constructor is SYNCHRONOUS on a prefetched context and the rows load inside the
 * grid ({@link DomainGrid.ready} resolves when they have).
 *
 * What it composes:
 * - the platform's grid decoration — the registered {@link ObjectHandler} for
 *   the table gets to customize the grid (captions, ref cells showing display
 *   names, name column first, system columns hidden), exactly as in the built-in
 *   Domain View;
 * - a {@link DomainFrameEditor} as THE writer of the editing state: in-grid
 *   edits are tracked, validated, highlighted (amber) and marked (red, with the
 *   message in the tooltip), and saved as ONE transaction;
 * - permission gating from {@link DomainTableCapabilities}: no `canEdit` means
 *   a read-only grid, no `canInsert` means no Add row, no `canDelete` means no
 *   Delete row, and columns the caller cannot write stay read-only.
 *
 * Capabilities are SNAPSHOT when the grid is built. A later grant change (or a
 * `grok.dapi.domains.invalidateUiCaches()`) affects only grids created after it —
 * an existing grid keeps the gating it was born with; rebuild it to pick the new
 * permissions up.
 *
 * Refreshing is the CALLER's decision: {@link DomainFrameEditor.refresh} rebuilds
 * and discards pending edits by design — check {@link DomainFrameEditor.isDirty}
 * first (this class rebinds the grid to whatever the editor rebuilds).
 */
export class DomainGrid extends Widget implements IEditorHost {
  readonly grid: Grid;
  readonly editor: DomainFrameEditor;

  private readonly _count: HTMLElement;
  private readonly _saveBar: HTMLElement;
  private readonly _toolbar: HTMLElement;
  private readonly _editable: boolean;
  /** Every toolbar button, disabled together while a save is in flight. */
  private readonly _buttons: HTMLButtonElement[] = [];
  private readonly _defaults?: {[column: string]: any};
  private readonly _ready: Promise<DomainGrid>;
  private _loadError: string | null = null;

  /**
   * SYNCHRONOUS over a prefetched context (the rows load afterwards — see
   * {@link ready}), or over an editor the caller already owns (a detail grid of a
   * master-detail page, sharing one save — {@link forEditor}).
   */
  constructor(source: IDomainTableContext | DomainFrameEditor, options?: DomainGridOptions) {
    super(ui.divV([], 'domain-ui-grid'));
    applyDomainGridStyles();
    const own = source instanceof DomainFrameEditor ? null : source;
    options = options ?? {};
    this.editor = own == null ? source as DomainFrameEditor
      : DomainFrameEditor.forContext(own, options);
    this._defaults = options.defaults;
    this._editable = (options.editable ?? true) && this.editor.capabilities.canEdit;
    this.grid = Grid.create(this.editor.dataFrame);
    // A viewer built outside the DOM carries the default 400x300 inline size it
    // was born with; it has to fill whatever host it is mounted in instead.
    this.grid.root.style.width = '100%';
    this.grid.root.style.height = '100%';

    this._count = ui.divText('', 'domain-ui-edit-count');
    this._saveBar = ui.divH([
      this._count,
      this._button('Save', () => this.editor.save(), 'Save every pending change as one transaction'),
      this._button('Cancel', () => this.editor.discard(), 'Discard every pending change'),
    ], 'domain-ui-save-bar');
    this._toolbar = ui.divH([
      this._editable && this.editor.capabilities.canInsert
        ? this._button('Add row', () => this.addRow(), 'Append a new row') : null,
      this._editable && this.editor.capabilities.canDelete
        ? this._button('Delete row', () => this.deleteRow(), 'Mark the current row for deletion') : null,
      this._saveBar,
    ], 'domain-ui-grid-toolbar');
    ui.setDisplay(this._toolbar, (options.toolbar ?? true) && this._editable);
    this.root.appendChild(this._toolbar);
    this.root.appendChild(ui.box(this.grid.root));

    this._decorate();
    this._wire();
    this._refreshBar();
    // Discoverable through `Widget.getAll()`: wrapping the widget for Dart is
    // what registers its root in the platform's widget map (and `detach()`
    // unregisters it, so a rebuilt grid leaves nothing behind).
    this.toDart();
    // An editor of our own starts empty: the query runs now, and the rebind rides
    // the same onRefreshed path a later refresh() takes.
    this._ready = own == null ? Promise.resolve(this) : this._load(options.query);
  }

  /** Builds the grid: resolves the table's registry metadata and capabilities,
   * runs the query, attaches the editor and decorates the grid through the
   * table's registered handler. */
  static async create(client: DomainTableClient, options?: DomainGridOptions): Promise<DomainGrid> {
    const context = await acquireDomainContext(client);
    return await new DomainGrid(context, options).ready;
  }

  /** Builds the grid over an editor the caller already owns (a detail grid of a
   * master-detail form sharing one save). */
  static forEditor(editor: DomainFrameEditor, options?: DomainGridOptions): DomainGrid {
    return new DomainGrid(editor, options ?? {});
  }

  /** The platform type name — what `Widget` discovery and a page's nested
   * status key on. */
  get type(): string { return 'DomainGrid'; }

  /** The frame the grid shows — replaced whenever the editor refreshes. */
  get dataFrame(): DataFrame { return this.editor.dataFrame; }

  /** Whether in-grid editing is on (the requested mode AND the table's `canEdit`). */
  get editable(): boolean { return this._editable; }

  /** {@link IEditorHost}: the pending batch a hosting page answers for. */
  get editors(): DomainFrameEditor[] { return [this.editor]; }

  /** Whether anything in the grid is unsaved. */
  get isDirty(): boolean { return this.editor.isDirty; }

  /**
   * Resolves when the first query has run — or when it FAILED, in which case the
   * grid shows the reason and this still resolves (a rejected promise nobody
   * awaits is an unhandled rejection). Idempotent: every call returns the same
   * promise; a grid built on a caller's editor is ready at once.
   */
  get ready(): Promise<this> { return this._ready as Promise<this>; }

  /** Releases the subscriptions of both the grid and its editor. */
  detach(): void {
    this.editor.detach();
    super.detach();
  }

  /** Applies the platform decoration, then the editing affordances. Runs again
   * after a refresh, because the frame (and the grid's columns) are new. */
  private _decorate(): void {
    DomainGrid.decorate(this.grid, this.editor.table, this.editor.dataFrame);
    // The handler decorates presentation; the service columns are OURS — hide
    // them here too, so a handler that overrides renderGrid entirely still
    // cannot expose (or let the user edit) the editing state.
    for (const name of DomainFrameEditor.SERVICE_COLUMNS) {
      const gc = this.grid.col(name);
      if (gc != null) {
        gc.visible = false;
        gc.editable = false;
      }
    }
    this._applyEditability();
  }

  /**
   * Customizes [grid] the way the platform customizes every grid over [table]'s
   * rows: column captions, reference cells showing display names, the name
   * column first, system and `~` columns hidden.
   *
   * Decoration goes through the handler that WINS dispatch for the table, which
   * is safe for plain handlers — one that does not override `renderGrid` falls
   * through to the platform meta from the JS side too.
   *
   * Which is not a one-liner, because `forEntity` resolves plenty of handlers
   * that must NOT be decorated through. It applies the collapse rule of the
   * built-in Domain View (`DomainView.refreshGrid`) branch for branch — with the
   * Dart one widened, since the Dart rule names `DomainRowMeta` while every Dart
   * meta reaches JS as one wrapping proxy:
   * - a DART meta collapses — the resolved proxy may be the inert GENERIC
   *   `'DomainRow'` fallback (no per-table meta registered, e.g. a session with
   *   domain databases not enabled), whose `renderGrid` early-returns into a raw
   *   grid; a registered per-table meta decorates identically to `own` anyway;
   * - a JS handler WITHOUT a real `renderGrid` collapses — the base method is a
   *   no-op the platform marks `isPlatformDefault`;
   * - everything else wins, INCLUDING a plugin handler that claims the table
   *   through `isApplicable` under a type of its own.
   *
   * `own`'s inherited `renderGrid` reaches the per-table meta regardless of
   * registration, so the collapse never loses decoration.
   */
  static decorate(grid: Grid, table: string, dataFrame?: DataFrame): void {
    const own = new DomainObjectHandler(table);
    const resolved = ObjectHandler.forEntity(own.newRow());
    const winner = resolved != null && DomainGrid._decorates(resolved, table) ? resolved : own;
    winner.renderGrid(grid, {items: dataFrame ?? grid.dataFrame});
  }

  /** Whether [handler] is worth decorating [table]'s grid through — see
   * {@link decorate}. */
  private static _decorates(handler: ObjectHandler, table: string): boolean {
    if (handler instanceof EntityMetaDartProxy)
      return false;
    return DomainGrid._typeOf(handler) === table ||
      (handler.renderGrid as any)?.isPlatformDefault !== true;
  }

  /** A handler's type, defensively: `type` is abstract on the base class and a
   * JS getter may throw. */
  private static _typeOf(handler: ObjectHandler): string | null {
    try {
      return handler.type;
    } catch (_) {
      return null;
    }
  }

  /** Runs the first query into the grid's own editor. Its outcome is
   * {@link ready}; a failure shows up in the grid instead of rejecting. */
  private async _load(query?: DomainQuerySpec): Promise<this> {
    ui.setUpdateIndicator(this.root, true);
    try {
      await this.editor.refresh(query ?? {});
    } catch (e: any) {
      this._loadError = `${e?.message ?? e}`;
      this.root.appendChild(ui.divText(`Cannot list ${this.editor.table}: ${this._loadError}`));
      log.error(`${this.editor.table}: ${this._loadError}`);
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
    return this;
  }

  private _wire(): void {
    // The single-writer contract: the grid never touches the service columns —
    // it reports the edit and the editor owns what happens to the state. Every
    // subscription rides the widget's own `subs`, which `super.detach()` releases.
    this._sub(this.grid.onCurrentCellChanged.subscribe((gc) => {
      const row = gc?.tableRowIndex;
      if (row != null)
        this.editor.beginEdit(row);
    }));
    this._sub(this.grid.onCellValueEdited.subscribe((gc) => {
      const row = gc?.tableRowIndex;
      const column = gc?.tableColumn?.name;
      if (row != null && column != null && !DomainFrameEditor.SERVICE_COLUMNS.includes(column))
        this.editor.commitEdit(row, column);
    }));
    this._sub(this.editor.onChanged.subscribe(() => {
      this._refreshBar();
      this.grid.invalidate();
    }));
    this._sub(this.editor.onRefreshed.subscribe((df) => {
      this.grid.dataFrame = df;
      this._decorate();
      this._refreshBar();
    }));
    // The editor is the single writer and closes while its transaction is in
    // flight — the grid must not offer an edit that would be refused (or lost).
    this._sub(this.editor.onSavingChanged.subscribe(() => {
      this._applyEditability();
      this._refreshBar();
    }));

    this._sub(this.grid.onCellPrepare((gc) => this._prepareCell(gc)));
    this._sub(this.grid.onCellTooltip((gc, x, y) => this._showTooltip(gc, x, y)));
  }

  /** rxjs subscriptions and grid StreamSubscriptions alike — both cancel through
   * `unsubscribe()`, which is all the inherited `subs` list needs. */
  private _sub(subscription: {unsubscribe(): void}): void {
    this.subs.push(subscription as any);
  }

  /** Amber for a pending edit, red for an invalid cell, orange for a row whose
   * version conflict the user dismissed. */
  private _prepareCell(gc: GridCell): void {
    if (!gc.isTableCell)
      return;
    const row = gc.tableRowIndex;
    const column = gc.tableColumn?.name;
    if (row == null || column == null)
      return;
    const error = this.editor.errorOf(row, column);
    if (error != null)
      gc.style.backColor = error.kind === 'conflict' ? CONFLICT_CELL_COLOR : INVALID_CELL_COLOR;
    else if (this.editor.isChanged(row, column))
      gc.style.backColor = DIRTY_CELL_COLOR;
  }

  private _showTooltip(gc: GridCell, x: number, y: number): boolean {
    if (!gc.isTableCell)
      return false;
    const row = gc.tableRowIndex;
    const column = gc.tableColumn?.name;
    const error = row == null || column == null ? null : this.editor.errorOf(row, column);
    if (error == null)
      return false;
    ui.tooltip.show(ui.divText(error.message), x, y);
    return true;
  }

  /** Read-only degradation and column security, in the platform's own shape:
   * the whole grid when the table is not editable (or while a save is in
   * flight), and per column otherwise — non-writable columns stay read-only.
   * Columns holding uuids (`ref`, user, group) edit in place: the platform's
   * grid opens their own picker anchored on the cell. */
  private _applyEditability(): void {
    const editable = this._editable && !this.editor.isSaving;
    this.grid.props.allowEdit = editable;
    for (let i = 0; i < this.grid.columns.length; i++) {
      const gc = this.grid.columns.byIndex(i);
      if (gc != null)
        gc.editable = false;
    }
    if (!editable)
      return;
    const writable = this.editor.capabilities.writableColumns;
    for (const p of this.editor.properties) {
      const gc = this.grid.col(p.name);
      if (gc != null)
        gc.editable = writable.includes(p.name);
    }
  }

  // ─────────────────────── actions (the shared Func vocabulary) ─────────────────────────

  /** Appends a row prefilled with the grid's `defaults` and puts the cursor on
   * it; -1 when the editor refused (a save is in flight). */
  addRow(): number {
    const row = this.editor.addRow(this._defaults);
    if (row < 0)
      return row;
    this.grid.dataFrame.currentRowIdx = row;
    this.editor.beginEdit(row);
    return row;
  }

  /** Marks the current row for deletion; false when there is no current row. */
  deleteRow(): boolean {
    const row = this.grid.dataFrame.currentRowIdx;
    if (row < 0) {
      balloon.warning('Select a row first');
      return false;
    }
    this.editor.markDeleted(row);
    return true;
  }

  /** Writes the pending batch as one transaction. */
  save(): Promise<boolean> { return this.editor.save(); }

  /** Drops the pending batch. */
  discard(): void { this.editor.discard(); }

  /** Re-runs the query — WITHOUT a prompt: refreshing IS discarding, and the gate
   * belongs to whoever owns the page (see {@link DomainFrameEditor.refresh}). */
  async refresh(): Promise<void> {
    await this.editor.refresh();
  }

  // ─────────────────────── the machine surface ─────────────────────────

  /**
   * What a grid reports about itself: its parts, how much is pending, and the
   * first blocking error. Deliberately NO per-cell `inputs` — a grid's machine
   * surface is its status, its functions and its {@link dataFrame}, which is a
   * first-class platform object already.
   */
  getWidgetStatus(): IWidgetStatus {
    const rows = this.editor.dataFrame.rowCount;
    const pending = this.editor.changeCount;
    const errors = this.editor.errorCount;
    return {
      parts: {grid: this.grid.root, toolbar: this._toolbar, saveBar: this._saveBar},
      hitAreas: {},
      shortcuts: {},
      events: [],
      description: `${this._editable ? 'An editable' : 'A read-only'} grid of ` +
        `${this.editor.table}: ${rows} row${rows === 1 ? '' : 's'}, ` +
        `${pending} unsaved change${pending === 1 ? '' : 's'}` +
        `${errors === 0 ? '' : `, ${errors} invalid cell${errors === 1 ? '' : 's'}`}.`,
      error: this._loadError ?? this._firstError(),
    };
  }

  /** The actions this caller may perform, as REAL platform Funcs from the shared
   * vocabulary — permission-gated exactly like the toolbar buttons. */
  getFunctions(): Func[] {
    if (!this._editable)
      return [DomainGrid.refreshFunc()];
    const funcs = [DomainGrid.saveFunc(), DomainGrid.discardFunc()];
    if (this.editor.capabilities.canInsert)
      funcs.push(DomainGrid.addRowFunc());
    if (this.editor.capabilities.canDelete)
      funcs.push(DomainGrid.deleteRowFunc());
    funcs.push(DomainGrid.refreshFunc());
    return funcs;
  }

  /** The first blocking cell error, named by its column. */
  private _firstError(): string | null {
    for (let row = 0; row < this.editor.dataFrame.rowCount; row++) {
      const errors = this.editor.errorsOf(row);
      for (const column of Object.keys(errors))
        if (errors[column].kind === 'error')
          return `${column}: ${errors[column].message}`;
    }
    return null;
  }

  private _button(text: string, action: () => void, tooltip: string): HTMLButtonElement {
    const button = ui.button(text, action, tooltip);
    this._buttons.push(button);
    return button;
  }

  private _refreshBar(): void {
    const n = this.editor.changeCount;
    this._count.textContent = this.editor.isSaving ? 'Saving…' : `${n} unsaved change${n === 1 ? '' : 's'}`;
    ui.setDisplay(this._saveBar, n > 0 || this.editor.isSaving);
    for (const button of this._buttons)
      button.disabled = this.editor.isSaving;
  }

  // ─────────────────────── the widget-action vocabulary ─────────────────────────

  /** Namespace every widget-action Func is registered under — what keeps `Save`
   * from colliding with another package's `Save`. */
  static readonly WIDGET_ACTION_NAMESPACE = 'DomainUi';

  /** Tags every widget-action Func carries: `domain-ui` marks the origin,
   * `widget-action` the calling convention (one `widget` input). */
  static readonly WIDGET_ACTION_TAGS = 'domain-ui,widget-action';

  /**
   * The widget-action Func called [name], registering it on first use.
   *
   * Idempotent by qualified name (`DomainUi:<name>`): the first registration wins and
   * every later call returns the same Func, so a library bundled into several
   * packages registers one vocabulary, not one per bundle.
   */
  static widgetActionFunc(name: string, run: (widget: any) => Promise<boolean>,
    description: string): Func {
    const existing = _funcs.get(name);
    if (existing != null)
      return existing;
    const func = functions.register({
      signature: `bool ${name}(widget widget)`,
      run: run,
      tags: DomainGrid.WIDGET_ACTION_TAGS,
      isAsync: true,
      namespace: DomainGrid.WIDGET_ACTION_NAMESPACE,
      options: {description: description},
    });
    _funcs.set(name, func);
    return func;
  }

  /** Writes the widget's pending changes (one transaction per editor); resolves to
   * whether everything landed. */
  static saveFunc(): Func {
    return DomainGrid.widgetActionFunc('Save', async (widget: any) =>
      typeof widget?.save === 'function' ? await widget.save() === true : false,
    'Saves the pending changes of the widget');
  }

  /** Drops the widget's pending changes. */
  static discardFunc(): Func {
    return DomainGrid.widgetActionFunc('Discard', async (widget: any) => {
      if (typeof widget?.discard !== 'function')
        return false;
      await widget.discard();
      return true;
    }, 'Discards the pending changes of the widget');
  }

  /** Puts the widget back to the values it was opened with. */
  static resetFunc(): Func {
    return DomainGrid.widgetActionFunc('Reset', async (widget: any) => {
      if (typeof widget?.reset !== 'function')
        return false;
      await widget.reset();
      return true;
    }, 'Restores the values the widget was opened with');
  }

  /** Appends an empty row to a grid's pending batch (prefilled with the grid's
   * `defaults`); false when the widget cannot insert. */
  static addRowFunc(): Func {
    return DomainGrid.widgetActionFunc('AddRow', async (widget: any) => {
      if (typeof widget?.addRow !== 'function')
        return false;
      return await widget.addRow() >= 0;
    }, 'Appends a new row to the widget');
  }

  /** Marks the widget's current row for deletion (the save writes it). */
  static deleteRowFunc(): Func {
    return DomainGrid.widgetActionFunc('DeleteRow', async (widget: any) => {
      if (typeof widget?.deleteRow !== 'function')
        return false;
      return await widget.deleteRow() === true;
    }, 'Marks the current row of the widget for deletion');
  }

  /** Re-reads the widget's rows from the server (through the unsaved-changes gate
   * of whatever owns it). */
  static refreshFunc(): Func {
    return DomainGrid.widgetActionFunc('Refresh', async (widget: any) => {
      if (typeof widget?.refresh !== 'function')
        return false;
      await widget.refresh();
      return true;
    }, 'Reloads the widget from the server');
  }

  /** A function name for [caption]: identifier characters only, so an action called
   * 'Edit...' registers as `Edit`. */
  static actionFuncName(caption: string): string {
    const name = `${caption ?? ''}`.replace(/[^A-Za-z0-9_]/g, '');
    return name === '' ? 'Action' : name;
  }
}

/** What a widget must implement to be acted on by the vocabulary. Every member is
 * optional: a widget offers the actions it can perform (see
 * `DomainForm.getFunctions`), and one that cannot is simply not offered. */
export interface IWidgetActions {
  save?(): Promise<boolean>;
  discard?(): void | Promise<void>;
  reset?(): void | Promise<void>;
  addRow?(): number | Promise<number>;
  deleteRow?(): boolean | Promise<boolean>;
  refresh?(): void | Promise<any>;
}

/** Wraps a custom {@link DomainAction} as a widget-action Func, through the same
 * idempotent registration — a plugin's extra action is as discoverable as `Save`. */
export function domainActionFunc(action: DomainAction): Func {
  return DomainGrid.widgetActionFunc(DomainGrid.actionFuncName(action.name), async () => {
    const result = await action.run();
    return result !== false;
  }, `${action.name}`);
}
