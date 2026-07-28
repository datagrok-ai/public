/**
 * Draft `grokky.*` helpers for the `datagrok-grid-customization` skill.
 *
 * These will eventually move into `src/claude/grokky-api.ts`, replacing
 * the silent no-op `sortRows` with a real `applySort` that mutates the
 * grid. Living here keeps the skill review self-contained: SKILL.md cites
 * helpers by name, examples.ts type-checks against them, eval/prompts.json
 * references them.
 *
 * Design notes:
 * - All helpers operate on a `DG.TableView` (the only view subtype that
 *   exposes `.grid`). Non-TableView callers get a clear error, not a
 *   silent miss.
 * - `colorCode` mirrors the polymorphic shape of the `colorCoding` field
 *   in df-and-columns' `setColumnMeta` so both skills teach one vocabulary.
 *   Linear on non-numeric throws (the only loud failure — silent miscolor
 *   would be worse).
 * - `configureGrid` is intentionally a batch operation with each field
 *   independent. Missing column names emit `console.warn` and skip.
 */
import * as DG from 'datagrok-api/dg';

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

function requireTableView(view: any): DG.TableView {
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grid customization requires a TableView (got: ' + (view?.type ?? typeof view) + ')');
  return view as DG.TableView;
}

function gridCol(view: DG.TableView, name: string): DG.GridColumn | null {
  // `view.grid.col(name)` is shorthand for `view.grid.columns.byName(name)`.
  // Returns null when the name does not match any data column.
  return view.grid.col(name);
}

// ---------------------------------------------------------------------------
// applySort
// ---------------------------------------------------------------------------

export type SortOrder = 'asc' | 'desc' | boolean;

/**
 * Apply a visual sort to the grid. Replaces the silent no-op `sortRows`
 * by actually mutating `view.grid` via `grid.sort(columns, orders)`.
 *
 * `orders[i]` is parallel to `columns[i]`. Missing entries default to
 * `'asc'`. Boolean `true` is `'asc'`, `false` is `'desc'`.
 *
 * **Trap.** `view.grid.sortByColumns = [...]` is a getter-only assignment
 * that silently does nothing. Always call `grid.sort()` (or this helper).
 *
 * @example
 *   applySort(view, ['Activity'], ['desc']);
 *   applySort(view, ['Class', 'Activity'], ['asc', 'desc']);
 *   applySort(view, ['Activity']);  // defaults to asc
 */
export function applySort(
  view: DG.TableView | any,
  columns: string[],
  orders?: SortOrder[],
): void {
  const tv = requireTableView(view);
  if (!Array.isArray(columns))
    throw new Error('applySort: "columns" must be a string[]');
  if (orders && orders.length !== columns.length && orders.length !== 0)
    throw new Error(`applySort: "orders" length (${orders.length}) must match "columns" length (${columns.length})`);
  const ascBools = columns.map((_, i) => {
    const o = orders?.[i] ?? 'asc';
    if (typeof o === 'boolean') return o;
    return o !== 'desc';
  });
  // Skip silently-unknown column names — surface the warning but still
  // sort by whatever does match.
  const known = columns.filter((c) => {
    if (tv.dataFrame.col(c)) return true;
    console.warn(`applySort: column "${c}" not found on dataFrame; skipping`);
    return false;
  });
  const knownBools = columns
    .map((c, i) => [c, ascBools[i]] as const)
    .filter(([c]) => !!tv.dataFrame.col(c))
    .map(([, b]) => b);
  tv.grid.sort(known, knownBools);
}

// ---------------------------------------------------------------------------
// clearSort
// ---------------------------------------------------------------------------

/**
 * Clear any sort applied to the grid.
 *
 * @example
 *   clearSort(view);
 */
export function clearSort(view: DG.TableView | any): void {
  const tv = requireTableView(view);
  tv.grid.sort([], []);
}

// ---------------------------------------------------------------------------
// configureGrid
// ---------------------------------------------------------------------------

export type WidthPolicy = 'Minimal' | 'Compact' | 'Optimal' | 'Maximal';

export type ConfigureGridOpts = {
  /** Set `gridCol.visible = false` for each name. */
  hide?: string[];
  /** Bulk "show only these" — uses `grid.columns.setVisible`. */
  show?: string[];
  /** Reorder via `grid.columns.setOrder`. Unlisted names retain relative position behind the listed set. */
  order?: string[];
  /** Per-column pixel widths. */
  widths?: Record<string, number>;
  /** Global width policy (`DG.ColumnWidthType`). */
  widthPolicy?: WidthPolicy;
  /** Pin (left). Each call appends to the pinned list. */
  pin?: string[];
  /** Unpin. */
  unpin?: string[];
  /** Data-side `col.meta.format` per column. Sugar over `setColumnMeta({format})`. */
  formats?: Record<string, string>;
  /** Grid-level row height (px). */
  rowHeight?: number;
  /** Freeze the first N visible columns (`IGridSettings.frozenColumns`). */
  frozenColumns?: number;
};

/**
 * Batch grid customization. Each field is applied independently; missing
 * column names are skipped with a `console.warn` (never thrown). The order
 * of application is stable:
 *
 *   1. visibility (hide/show)
 *   2. order
 *   3. widths
 *   4. widthPolicy (after explicit widths so a global policy can override)
 *   5. pin / unpin
 *   6. formats (data-side `col.meta.format`)
 *   7. rowHeight
 *   8. frozenColumns
 *
 * @example
 *   configureGrid(view, {hide: ['_id'], pin: ['SMILES'], widths: {SMILES: 250}});
 */
export function configureGrid(view: DG.TableView | any, opts: ConfigureGridOpts): void {
  const tv = requireTableView(view);
  const grid = tv.grid;

  // 1. Visibility — `hide` first, then `show` (show wins because it's the
  //    bulk override). Both tolerate unknown names.
  if (opts.hide) {
    for (const name of opts.hide) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: hide: column "${name}" not found, skipping`); continue; }
      gc.visible = false;
    }
  }
  if (opts.show) {
    const known = opts.show.filter((n) => {
      if (gridCol(tv, n)) return true;
      console.warn(`configureGrid: show: column "${n}" not found, skipping`);
      return false;
    });
    // setVisible expects all desired-visible names.
    grid.columns.setVisible(known);
  }

  // 2. Order — `setOrder` accepts the names; unlisted retain relative position behind.
  if (opts.order) {
    const known = opts.order.filter((n) => {
      if (gridCol(tv, n)) return true;
      console.warn(`configureGrid: order: column "${n}" not found, skipping`);
      return false;
    });
    if (known.length > 0)
      grid.columns.setOrder(known);
  }

  // 3. Widths — per-column pixel widths.
  if (opts.widths) {
    for (const [name, px] of Object.entries(opts.widths)) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: widths: column "${name}" not found, skipping`); continue; }
      if (typeof px !== 'number' || !Number.isFinite(px) || px < 0) {
        console.warn(`configureGrid: widths: "${name}" → ${px} is not a non-negative number, skipping`);
        continue;
      }
      gc.width = px;
    }
  }

  // 4. Width policy — global, applies after per-column widths so callers
  //    can opt into "auto-fit everything".
  if (opts.widthPolicy) {
    const enumVal = (DG as any).ColumnWidthType?.[opts.widthPolicy];
    if (enumVal === undefined)
      console.warn(`configureGrid: widthPolicy "${opts.widthPolicy}" not recognized (expected Minimal | Compact | Optimal | Maximal)`);
    else
      grid.setColumnsWidthType(enumVal);
  }

  // 5. Pin / unpin.
  if (opts.pin) {
    for (const name of opts.pin) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: pin: column "${name}" not found, skipping`); continue; }
      try { gc.pin(); }
      catch (e) { console.warn(`configureGrid: pin "${name}": ${(e as Error).message}`); }
    }
  }
  if (opts.unpin) {
    for (const name of opts.unpin) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: unpin: column "${name}" not found, skipping`); continue; }
      try { gc.unpin(); }
      catch (e) { console.warn(`configureGrid: unpin "${name}": ${(e as Error).message}`); }
    }
  }

  // 6. Formats — DATA-SIDE. col.meta.format persists across viewers.
  if (opts.formats) {
    for (const [name, format] of Object.entries(opts.formats)) {
      const col = tv.dataFrame.col(name);
      if (!col) { console.warn(`configureGrid: formats: column "${name}" not found, skipping`); continue; }
      col.meta.format = format;
    }
  }

  // 7. Row height.
  if (opts.rowHeight !== undefined) {
    if (typeof opts.rowHeight !== 'number' || !Number.isFinite(opts.rowHeight) || opts.rowHeight <= 0)
      console.warn(`configureGrid: rowHeight ${opts.rowHeight} is not a positive number, skipping`);
    else
      grid.setOptions({rowHeight: opts.rowHeight});
  }

  // 8. Frozen columns.
  if (opts.frozenColumns !== undefined) {
    if (typeof opts.frozenColumns !== 'number' || !Number.isInteger(opts.frozenColumns) || opts.frozenColumns < 0)
      console.warn(`configureGrid: frozenColumns ${opts.frozenColumns} is not a non-negative integer, skipping`);
    else
      grid.setOptions({frozenColumns: opts.frozenColumns});
  }
}

// ---------------------------------------------------------------------------
// colorCode
// ---------------------------------------------------------------------------

export type ColorCodeSpec =
  | {
      kind: 'linear';
      range?: (string | number)[];
      min?: number;
      max?: number;
      belowMinColor?: string;
      aboveMaxColor?: string;
    }
  | {
      kind: 'categorical';
      colors: Record<string, string | number>;
      fallbackColor?: string | number;
      matchType?: 'exact' | 'regex';
    }
  | {
      kind: 'conditional';
      rules: Record<string, string | number>;
    }
  | {kind: 'off'};

/**
 * Data-side column color coding. Mirrors the `colorCoding` branch of
 * `setColumnMeta` in df-and-columns — keep shapes identical so both skills
 * teach one vocabulary.
 *
 * Data-side means: writes `col.meta.colors.*`, which propagates to **every**
 * viewer reading the column. Scatter-plot color axis bound to this column
 * picks up the gradient; CSV export preserves the tags.
 *
 * **Throws** on `kind: 'linear'` against a non-numeric column. Silent
 * miscolor is worse than a clear error.
 *
 * @example
 *   colorCode(col, {kind: 'linear', range: ['#00FF00', '#FFFF00', '#FF0000']});
 *   colorCode(col, {kind: 'categorical', colors: {'A': '#FF0000', 'B': '#00FF00'}});
 *   colorCode(col, {kind: 'conditional', rules: {'0-50': '#FF0000', '50-': '#00FF00'}});
 *   colorCode(col, {kind: 'off'});
 */
export function colorCode(col: DG.Column, spec: ColorCodeSpec): void {
  if (!col) throw new Error('colorCode: column is null or undefined');

  if (spec.kind === 'linear') {
    if (!col.isNumerical)
      throw new Error(`colorCode: linear color coding requires a numerical column; "${col.name}" is ${col.type}`);
    const opts: {min?: number; max?: number; belowMinColor?: string; aboveMaxColor?: string} = {};
    if (spec.min !== undefined) opts.min = spec.min;
    if (spec.max !== undefined) opts.max = spec.max;
    if (spec.belowMinColor !== undefined) opts.belowMinColor = spec.belowMinColor;
    if (spec.aboveMaxColor !== undefined) opts.aboveMaxColor = spec.aboveMaxColor;
    (col.meta.colors as any).setLinear(spec.range ?? null, Object.keys(opts).length ? opts : null);
    return;
  }

  if (spec.kind === 'categorical') {
    const opts: {fallbackColor?: string | number; matchType?: string} = {};
    if (spec.fallbackColor !== undefined) opts.fallbackColor = spec.fallbackColor;
    if (spec.matchType !== undefined) opts.matchType = spec.matchType;
    (col.meta.colors as any).setCategorical(spec.colors as any, Object.keys(opts).length ? (opts as any) : null);
    return;
  }

  if (spec.kind === 'conditional') {
    if (!spec.rules || Object.keys(spec.rules).length === 0)
      throw new Error('colorCode: conditional spec requires non-empty `rules`');
    (col.meta.colors as any).setConditional(spec.rules as any);
    return;
  }

  if (spec.kind === 'off') {
    col.meta.colors.setDisabled();
    return;
  }

  throw new Error(`colorCode: unknown kind "${(spec as any).kind}"`);
}

// ---------------------------------------------------------------------------
// resetGrid
// ---------------------------------------------------------------------------

export type ResetGridOpts = {
  /** Re-show every column (`gridCol.visible = true`). Default true. */
  visibility?: boolean;
  /** Restore default widths via `ColumnWidthType.Optimal`. Default true. */
  widths?: boolean;
  /** Clear the sort. Default true. */
  sort?: boolean;
  /** Turn off color coding on every column. Default true. */
  colors?: boolean;
};

/**
 * Selective grid reset. Does **not** close other viewers (use the
 * `datagrok-viewers` skill's `closeAllViewers` for that, or
 * `view.resetLayout()` for the blunt "everything").
 *
 * Default opts: `{visibility: true, widths: true, sort: true, colors: true}`.
 *
 * @example
 *   resetGrid(view);                       // full default reset
 *   resetGrid(view, {sort: true, visibility: false, widths: false, colors: false});
 */
export function resetGrid(view: DG.TableView | any, opts: ResetGridOpts = {}): void {
  const tv = requireTableView(view);
  const grid = tv.grid;
  const o: Required<ResetGridOpts> = {
    visibility: opts.visibility ?? true,
    widths:     opts.widths     ?? true,
    sort:       opts.sort       ?? true,
    colors:     opts.colors     ?? true,
  };

  if (o.visibility) {
    // Skip index 0 (the row header) — only data columns get the toggle.
    const len = grid.columns.length;
    for (let i = 1; i < len; i++) {
      const gc = grid.columns.byIndex(i);
      if (gc) gc.visible = true;
    }
  }

  if (o.widths) {
    const widthType = (DG as any).ColumnWidthType?.Optimal;
    if (widthType !== undefined)
      grid.setColumnsWidthType(widthType);
  }

  if (o.sort)
    clearSort(tv);

  if (o.colors) {
    for (const col of tv.dataFrame.columns.toList()) {
      try { col.meta.colors.setDisabled(); }
      catch (e) { console.warn(`resetGrid: setDisabled on "${col.name}": ${(e as Error).message}`); }
    }
  }
}

// ---------------------------------------------------------------------------
// pinColumn / unpinColumn
// ---------------------------------------------------------------------------

/**
 * Pin one column on the left of the grid. Convenience over
 * `view.grid.col(name).pin()`. Pin is left-only — there's no positional
 * arg.
 *
 * @example
 *   pinColumn(view, 'SMILES');
 */
export function pinColumn(view: DG.TableView | any, name: string): void {
  const tv = requireTableView(view);
  const gc = gridCol(tv, name);
  if (!gc) {
    console.warn(`pinColumn: column "${name}" not found, skipping`);
    return;
  }
  gc.pin();
}

/**
 * Unpin one column.
 *
 * @example
 *   unpinColumn(view, 'SMILES');
 */
export function unpinColumn(view: DG.TableView | any, name: string): void {
  const tv = requireTableView(view);
  const gc = gridCol(tv, name);
  if (!gc) {
    console.warn(`unpinColumn: column "${name}" not found, skipping`);
    return;
  }
  gc.unpin();
}
