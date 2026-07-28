/**
 * Draft `grokky.*` helpers for the `datagrok-selection` skill.
 *
 * These will eventually replace / extend the existing `selectRows` in
 * `src/claude/grokky-api.ts` and add the other selection helpers. Living
 * here keeps the skill review self-contained: SKILL.md cites these by
 * name, examples.ts type-checks against them, eval/prompts.json references
 * them.
 *
 * Polarity convention (this skill's #1 footgun, opposite of filter):
 *   - `df.selection.get(i) === true`    — row is **SELECTED**.
 *   - `df.selection.setAll(false)`      — clears the selection (none selected).
 *   - `df.selection.setAll(true)`       — selects every row (including filtered-out!).
 *   - `selectRows(df, pred, ...)`       — `pred(i) === true` → row goes **IN** the new selection.
 *
 * Contrast with filtering: `df.filter.setAll(true)` clears the filter,
 * `df.filter.setAll(false)` hides every row. Same BitSet class, opposite
 * "clear" idiom — `clearSelection` and `clearFilter` exist to keep them
 * separate.
 */
import * as DG from 'datagrok-api/dg';

// ---------------------------------------------------------------------------
// selectRows — polymorphic across input shape × mode
// ---------------------------------------------------------------------------

/** How a new mask combines with `df.selection`. Default is `'replace'`. */
export type SelectMode = 'replace' | 'add' | 'remove' | 'intersect';

/**
 * Acceptable input shapes for `selectRows`. `ArrayLike<number>` covers
 * `number[]`, `Int32Array`, and `Uint32Array` — so `getSelectedIndexes()`
 * output pipes back in cleanly.
 */
export type SelectInput =
  | DG.BitSet
  | ((i: number) => boolean)
  | ArrayLike<number>;

export type SelectRowsOpts = {
  /** Combine semantics — see `SelectMode`. Default `'replace'`. */
  mode?: SelectMode;
};

/**
 * Set / extend / narrow the selection on `df`. The one polymorphic entry
 * point: accepts a BitSet, a predicate, an index list, or an `Int32Array`,
 * and one of four modes.
 *
 * Polarity: with a predicate, `pred(i) === true` → row goes **IN** the new
 * selection. Same direction as `filterByPredicate`'s "keep" convention. The
 * polarity flip is only at the *clear* level (see `clearSelection`).
 *
 * Mode behaviour:
 * - `'replace'` (default) — predicate → `init` (buffer-direct, single pass);
 *   BitSet → `copyFrom`; indices → `setAll(false)` then batched `set(i, true,
 *   false)` followed by one `fireChanged()`. Empty index list clears.
 * - `'add'`      — union (OR) onto current selection.
 * - `'remove'`   — subtract (AND NOT) from current selection.
 * - `'intersect'`— intersect (AND) with current selection.
 *
 * Empty input shape with mode `'replace'` clears the selection. With the
 * other modes it is a no-op.
 *
 * @example
 *   selectRows(df, (i) => df.getCol('age').get(i) > 40);
 *   selectRows(df, [0, 3, 5, 7]);
 *   selectRows(df, getSelectedIndexes, {mode: 'add'});
 *   selectRows(df, otherBitSet, {mode: 'intersect'});
 */
export function selectRows(
  df: DG.DataFrame,
  input: SelectInput,
  opts: SelectRowsOpts = {},
): void {
  const mode: SelectMode = opts.mode ?? 'replace';
  const sel = df.selection;

  // ---- BitSet input ----
  if (input instanceof DG.BitSet) {
    switch (mode) {
      case 'replace':   sel.copyFrom(input); return;
      case 'add':       sel.or(input); return;
      case 'remove':    sel.andNot(input); return;
      case 'intersect': sel.and(input); return;
    }
  }

  // ---- Predicate input ----
  if (typeof input === 'function') {
    if (mode === 'replace') {
      // Buffer-direct, single notification — the fast path.
      sel.init(input);
      return;
    }
    // For add/remove/intersect we need a materialised mask to combine with.
    const mask = DG.BitSet.create(df.rowCount, input);
    switch (mode) {
      case 'add':       sel.or(mask); return;
      case 'remove':    sel.andNot(mask); return;
      case 'intersect': sel.and(mask); return;
    }
  }

  // ---- ArrayLike<number> input (number[] | Int32Array | Uint32Array) ----
  // `for..of` works on all three; TypeScript narrows via the type predicates above.
  const indices = input as ArrayLike<number>;
  const len = indices.length;

  if (mode === 'replace') {
    sel.setAll(false, false);
    for (let k = 0; k < len; k++) sel.set(indices[k], true, false);
    sel.fireChanged();
    return;
  }
  if (mode === 'add') {
    if (len === 0) return;                   // no-op
    for (let k = 0; k < len; k++) sel.set(indices[k], true, false);
    sel.fireChanged();
    return;
  }
  if (mode === 'remove') {
    if (len === 0) return;                   // no-op
    for (let k = 0; k < len; k++) sel.set(indices[k], false, false);
    sel.fireChanged();
    return;
  }
  // mode === 'intersect' — keep only the listed indices that are currently selected.
  if (len === 0) {
    sel.setAll(false);
    return;
  }
  const keep = new Set<number>();
  for (let k = 0; k < len; k++) keep.add(indices[k]);
  sel.init((i) => keep.has(i) && sel.get(i));
}

// ---------------------------------------------------------------------------
// clearSelection — POLARITY-TRAP KILLER
// ---------------------------------------------------------------------------

/**
 * Deselect every row. Wraps `df.selection.setAll(false)`.
 *
 * **Polarity warning.** This is the **OPPOSITE** of `clearFilter`. The
 * pair to memorise:
 * ```
 *   grokky.clearFilter(view);   // setAll(TRUE)  — every row visible
 *   grokky.clearSelection(df);  // setAll(FALSE) — none selected
 * ```
 * Both helpers exist precisely so Claude doesn't have to track which mask
 * uses which polarity.
 *
 * @example
 *   clearSelection(df);
 */
export function clearSelection(df: DG.DataFrame): void {
  df.selection.setAll(false);
}

// ---------------------------------------------------------------------------
// selectAll
// ---------------------------------------------------------------------------

/**
 * Mark every row in `df` as selected. Wraps `df.selection.setAll(true)`.
 *
 * **This selects every row in the DataFrame, including filtered-out rows.**
 * For "select all currently-visible rows" (only rows passing the filter),
 * use `selectionFromFilter(df)` instead.
 *
 * @example
 *   selectAll(df);                  // every row
 *   selectionFromFilter(df);        // only visible rows
 */
export function selectAll(df: DG.DataFrame): void {
  df.selection.setAll(true);
}

// ---------------------------------------------------------------------------
// invertSelection
// ---------------------------------------------------------------------------

/**
 * Flip the selection mask in place. Selected rows become unselected and
 * vice versa. Wraps `df.selection.invert()`.
 *
 * @example
 *   invertSelection(df);
 */
export function invertSelection(df: DG.DataFrame): void {
  df.selection.invert();
}

// ---------------------------------------------------------------------------
// selectedDf — symmetric to filteredDf
// ---------------------------------------------------------------------------

export type SelectedDfOpts = {
  /** Narrow to a column subset. */
  cols?: string[];
  /** Carry the selection mask onto the clone. */
  withSelection?: boolean;
};

/**
 * Non-destructive — returns a clone of the currently-selected rows. The
 * source DataFrame is unmodified.
 *
 * Distinct from `filteredDf` (uses `df.filter`) and from `filterFromSelection`
 * (mutates `df.filter`, no clone). Same options shape as `filteredDf`.
 *
 * If the selection is empty, returns a zero-row clone (no error, no toast).
 * For user-facing flows, guard with `df.selection.anyTrue` first.
 *
 * @example
 *   const subset = selectedDf(df);
 *   const small  = selectedDf(df, {cols: ['smiles', 'activity']});
 */
export function selectedDf(df: DG.DataFrame, opts: SelectedDfOpts = {}): DG.DataFrame {
  return df.clone(df.selection, opts.cols ?? null, opts.withSelection ?? false);
}

// ---------------------------------------------------------------------------
// Cross-skill bridges
// ---------------------------------------------------------------------------

/**
 * "Show only the selected rows." Copies `df.selection` onto `df.filter`.
 *
 * Fires `onFilterChanged` (and the full filter lifecycle), but does **not**
 * fire `onSelectionChanged` — the selection itself is unchanged.
 *
 * Polarity match: both masks use `true` = "yes" in their respective
 * domains, so no inversion is needed. The names align: visible → selected,
 * selected → visible.
 *
 * If UI filter widgets are active in a `TableView`, their next
 * `onRowsFiltering` pass will overwrite this contribution. To make the
 * filter stick, call `grokky.clearFilter(view)` first (which disables the
 * UI filter group), then `filterFromSelection(df)`.
 *
 * @example
 *   filterFromSelection(df);
 */
export function filterFromSelection(df: DG.DataFrame): void {
  df.filter.copyFrom(df.selection);
}

/**
 * "Select every currently-visible row." Copies `df.filter` onto
 * `df.selection`.
 *
 * Fires `onSelectionChanged`, but does **not** fire a filter event — the
 * filter itself is unchanged.
 *
 * Polarity match: both masks use `true` = "yes" in their respective domains,
 * so no inversion is needed.
 *
 * @example
 *   selectionFromFilter(df);
 */
export function selectionFromFilter(df: DG.DataFrame): void {
  df.selection.copyFrom(df.filter);
}

// ---------------------------------------------------------------------------
// setCurrentRow — distinct from selection
// ---------------------------------------------------------------------------

/**
 * Move the current-row pointer to `idx`. Wraps `df.currentRowIdx = idx`.
 *
 * **Current row is NOT selection.** It is a single-row focus pointer
 * (distinct concept). Setting `currentRowIdx` does not select the row;
 * selecting rows does not move the current row. If the user says
 * "highlight row 5", the intent is ambiguous — pick selection or focus
 * based on surrounding context, or ask.
 *
 * The platform sets `currentRowIdx` on most user interactions (grid click,
 * arrow keys, viewer point click). Programmatic `setCurrentRow` is for
 * scripted navigation.
 *
 * Throws on out-of-range index. Use `df.currentRowIdx = -1` directly for
 * the explicit "no current row" state.
 *
 * @example
 *   setCurrentRow(df, 5);
 */
export function setCurrentRow(df: DG.DataFrame, idx: number): void {
  if (!Number.isInteger(idx))
    throw new Error(`setCurrentRow: idx must be an integer, got ${idx}`);
  if (idx < 0 || idx >= df.rowCount)
    throw new Error(`setCurrentRow: idx ${idx} out of range [0, ${df.rowCount})`);
  df.currentRowIdx = idx;
}

// ---------------------------------------------------------------------------
// describeSelection — JSON summary for debugging / agent reporting
// ---------------------------------------------------------------------------

/** Shape returned by `describeSelection`. */
export type SelectionDescription = {
  /** Number of selected rows. */
  count: number;
  /** Total rows in the DataFrame. */
  total: number;
  /** Indices of selected rows (truncated to first 100). */
  indexes: number[];
  /** Current-row pointer (independent of selection). `-1` means none. */
  currentRowIdx: number;
  /**
   * Column values for the first selected row, for the first ~5 columns
   * (by display order). Absent if nothing is selected.
   */
  sample?: Record<string, unknown>;
};

const MAX_INDEXES_IN_DESCRIBE = 100;
const MAX_SAMPLE_COLUMNS = 5;

/**
 * Read-only snapshot of the selection — useful for "what's selected?"
 * prompts and for agent debug output.
 *
 * The `indexes` array is truncated to the first 100 entries to keep agent
 * output sane; `count` always reflects the real total.
 *
 * @example
 *   describeSelection(df);
 *   // → {count: 3, total: 100, indexes: [5, 7, 11], currentRowIdx: 5,
 *   //    sample: {smiles: 'c1ccccc1', activity: 7.2}}
 */
export function describeSelection(df: DG.DataFrame): SelectionDescription {
  const sel = df.selection;
  const allIndexes = sel.getSelectedIndexes();
  const indexes: number[] = [];
  const cap = Math.min(allIndexes.length, MAX_INDEXES_IN_DESCRIBE);
  for (let k = 0; k < cap; k++) indexes.push(allIndexes[k]);

  const out: SelectionDescription = {
    count: sel.trueCount,
    total: df.rowCount,
    indexes,
    currentRowIdx: df.currentRowIdx,
  };

  if (allIndexes.length > 0) {
    const firstIdx = allIndexes[0];
    const cols = df.columns;
    const sample: Record<string, unknown> = {};
    const colCap = Math.min(cols.length, MAX_SAMPLE_COLUMNS);
    for (let c = 0; c < colCap; c++) {
      const col = cols.byIndex(c);
      sample[col.name] = col.get(firstIdx);
    }
    out.sample = sample;
  }

  return out;
}
