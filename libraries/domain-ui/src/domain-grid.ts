/**
 * {@link DomainGrid} — an editable grid over a domain table: the platform's own
 * grid decoration, batch editing through a {@link DomainFrameEditor}, and a
 * save/cancel bar, all permission-gated.
 *
 * @module domain-grid
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DomainFrameEditor, DomainFrameEditorOptions, isReferenceProperty,
  SERVICE_COLUMNS} from './frame-editor';
import {applyDomainUiStyles} from './styles';

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
 * const grid = await DomainGrid.create(grok.dapi.domains.table('grit.issue'));
 * grok.shell.newView('Issues', [grid.root]);
 * ```
 *
 * What it composes:
 * - the platform's grid decoration — the registered {@link DG.ObjectHandler} for
 *   the table gets to customize the grid (captions, ref cells showing display
 *   names, name column first, system columns hidden), exactly as in the built-in
 *   Domain View;
 * - a {@link DomainFrameEditor} as THE writer of the editing state: in-grid
 *   edits are tracked, validated, highlighted (amber) and marked (red, with the
 *   message in the tooltip), and saved as ONE transaction;
 * - permission gating from {@link DG.DomainTableCapabilities}: no `canEdit` means
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
export class DomainGrid {
  /** Mount this: toolbar + grid. */
  readonly root: HTMLElement;
  readonly grid: DG.Grid;
  readonly editor: DomainFrameEditor;

  private readonly _count: HTMLElement;
  private readonly _saveBar: HTMLElement;
  private readonly _toolbar: HTMLElement;
  /** rxjs subscriptions and grid StreamSubscriptions alike — both cancel
   * through `unsubscribe()`. */
  private readonly _subs: {unsubscribe(): void}[] = [];
  private readonly _editable: boolean;
  /** Every toolbar button, disabled together while a save is in flight. */
  private readonly _buttons: HTMLButtonElement[] = [];
  private readonly _defaults?: {[column: string]: any};

  private constructor(editor: DomainFrameEditor, options: DomainGridOptions) {
    applyDomainUiStyles();
    this.editor = editor;
    this._defaults = options.defaults;
    this._editable = (options.editable ?? true) && editor.capabilities.canEdit;
    this.grid = DG.Grid.create(editor.dataFrame);
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
      this._editable && editor.capabilities.canInsert
        ? this._button('Add row', () => this._addRow(), 'Append a new row') : null,
      this._editable && editor.capabilities.canDelete
        ? this._button('Delete row', () => this._deleteCurrentRow(), 'Mark the current row for deletion') : null,
      this._saveBar,
    ], 'domain-ui-grid-toolbar');
    ui.setDisplay(this._toolbar, (options.toolbar ?? true) && this._editable);
    this.root = ui.divV([this._toolbar, ui.box(this.grid.root)], 'domain-ui-grid');

    this._decorate();
    this._wire();
    this._refreshBar();
  }

  /** Builds the grid: probes capabilities, runs the query, attaches the editor
   * and decorates the grid through the table's registered handler. */
  static async create(client: DG.DomainTableClient, options?: DomainGridOptions): Promise<DomainGrid> {
    const editor = await DomainFrameEditor.create(client, options);
    return new DomainGrid(editor, options ?? {});
  }

  /** Builds the grid over an editor the caller already owns (a detail grid of a
   * master-detail form sharing one save). */
  static forEditor(editor: DomainFrameEditor, options?: DomainGridOptions): DomainGrid {
    return new DomainGrid(editor, options ?? {});
  }

  /** The frame the grid shows — replaced whenever the editor refreshes. */
  get dataFrame(): DG.DataFrame { return this.editor.dataFrame; }

  /** Whether in-grid editing is on (the requested mode AND the table's `canEdit`). */
  get editable(): boolean { return this._editable; }

  /** Releases the subscriptions of both the grid and its editor. */
  detach(): void {
    for (const sub of this._subs)
      sub.unsubscribe();
    this._subs.length = 0;
    this.editor.detach();
  }

  /** Applies the platform decoration, then the editing affordances. Runs again
   * after a refresh, because the frame (and the grid's columns) are new. */
  private _decorate(): void {
    DomainGrid.decorate(this.grid, this.editor.table, this.editor.dataFrame);
    // The handler decorates presentation; the service columns are OURS — hide
    // them here too, so a handler that overrides renderGrid entirely still
    // cannot expose (or let the user edit) the editing state.
    for (const name of SERVICE_COLUMNS) {
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
   * built-in Domain View (`DomainView.refreshGrid`), member for member:
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
  static decorate(grid: DG.Grid, table: string, dataFrame?: DG.DataFrame): void {
    const own = new DG.DomainObjectHandler(table);
    const resolved = DG.ObjectHandler.forEntity(own.newRow());
    const winner = resolved != null && DomainGrid._decorates(resolved, table) ? resolved : own;
    winner.renderGrid(grid, {items: dataFrame ?? grid.dataFrame});
  }

  /** Whether [handler] is worth decorating [table]'s grid through — see
   * {@link decorate}. */
  private static _decorates(handler: DG.ObjectHandler, table: string): boolean {
    if (handler instanceof DG.EntityMetaDartProxy)
      return false;
    return DomainGrid._typeOf(handler) === table ||
      (handler.renderGrid as any)?.isPlatformDefault !== true;
  }

  /** A handler's type, defensively: `type` is abstract on the base class and a
   * JS getter may throw. */
  private static _typeOf(handler: DG.ObjectHandler): string | null {
    try {
      return handler.type;
    } catch (_) {
      return null;
    }
  }

  private _wire(): void {
    // The single-writer contract: the grid never touches the service columns —
    // it reports the edit and the editor owns what happens to the state.
    this._subs.push(this.grid.onCurrentCellChanged.subscribe((gc) => {
      const row = gc?.tableRowIndex;
      if (row != null)
        this.editor.beginEdit(row);
    }));
    this._subs.push(this.grid.onCellValueEdited.subscribe((gc) => {
      const row = gc?.tableRowIndex;
      const column = gc?.tableColumn?.name;
      if (row != null && column != null && !SERVICE_COLUMNS.includes(column))
        this.editor.commitEdit(row, column);
    }));
    this._subs.push(this.editor.onChanged.subscribe(() => {
      this._refreshBar();
      this.grid.invalidate();
    }));
    this._subs.push(this.editor.onRefreshed.subscribe((df) => {
      this.grid.dataFrame = df;
      this._decorate();
      this._refreshBar();
    }));
    // The editor is the single writer and closes while its transaction is in
    // flight — the grid must not offer an edit that would be refused (or lost).
    this._subs.push(this.editor.onSavingChanged.subscribe(() => {
      this._applyEditability();
      this._refreshBar();
    }));

    this._subs.push(this.grid.onCellPrepare((gc) => this._prepareCell(gc)));
    this._subs.push(this.grid.onCellTooltip((gc, x, y) => this._showTooltip(gc, x, y)));
  }

  /** Amber for a pending edit, red for an invalid cell, orange for a row whose
   * version conflict the user dismissed. */
  private _prepareCell(gc: DG.GridCell): void {
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

  private _showTooltip(gc: DG.GridCell, x: number, y: number): boolean {
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
   * flight), and per column otherwise — non-writable columns and reference
   * columns (uuids, whose editing path is a picker) stay read-only. */
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
        gc.editable = writable.includes(p.name) && !isReferenceProperty(p);
    }
  }

  private _addRow(): void {
    const row = this.editor.addRow(this._defaults);
    // Refused while a save is in flight (the editor says so itself).
    if (row < 0)
      return;
    this.grid.dataFrame.currentRowIdx = row;
    this.editor.beginEdit(row);
  }

  private _deleteCurrentRow(): void {
    const row = this.grid.dataFrame.currentRowIdx;
    if (row < 0) {
      grok.shell.warning('Select a row first');
      return;
    }
    this.editor.markDeleted(row);
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
}
