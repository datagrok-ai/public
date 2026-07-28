/**
 * Draft `grokky.*` helpers for the `datagrok-filtering` skill.
 *
 * These will eventually move into `src/claude/grokky-api.ts`. Living here
 * keeps the skill review self-contained: SKILL.md cites these by name,
 * examples.ts type-checks against them, eval/prompts.json references them.
 *
 * Polarity convention (this skill's #1 footgun):
 *   - `filterByPredicate(df, pred)`  — `pred(i) === true` → **KEEP** row i.
 *   - `dropRows(df, pred)`           — `pred(i) === true` → **DROP** row i.
 *   - `df.filter.get(i) === true`    — row passes (visible).
 *   - `df.filter.setAll(true)`       — clears the filter (every row visible).
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

/** Row count above which substructure filtering auto-switches to the BitSet path. */
const LARGE_DF_THRESHOLD = 50_000;

// ---------------------------------------------------------------------------
// filterRows — single-column UI filter (range / values / contains / substructure)
// ---------------------------------------------------------------------------

export type FilterCriteria = {
  /** Numeric range filter. Either bound may be omitted. */
  min?: number;
  max?: number;
  /** Categorical in-set filter (single value = equals). */
  values?: string[];
  /** Free-text substring filter on a string column. */
  contains?: string;
  /** Free-text regex filter on a string column. */
  regex?: string;
  /** Substructure query — SMILES, SMARTS, or molblock. Auto-detected. */
  substructure?: string;
  /** Forwarded to substructure filter widget. */
  searchType?: 'substructure' | 'similar' | 'exact';
};

export type FilterRowsOpts = {
  /**
   * `true` (default) — use the FilterGroup UI path (widget visible in panel,
   * persists in saved layouts).
   * `false` — write the resulting BitSet onto `df.filter` directly (no UI).
   * Substructure on a DF > 50 000 rows defaults to `false` for perf.
   */
  ui?: boolean;
  /** Forwarded to substructure search. */
  molBlockFailover?: string;
};

/**
 * Filter rows on one column. Accepts either a `DG.TableView` (UI path) or a
 * `DG.DataFrame` (BitSet path). Async because the substructure conversion
 * goes through Dart-backed `Chem:convertMolNotation`.
 *
 * Filter polarity: rows **matching** the criteria are kept; the rest are
 * hidden. For multi-column conditions, drop down to `filterByPredicate`.
 *
 * @example
 *   await filterRows(view, 'MW', {min: 200, max: 500});
 *   await filterRows(view, 'category', {values: ['A', 'B']});
 *   await filterRows(view, 'smiles', {substructure: 'c1ccccc1'});
 */
export async function filterRows(
  target: DG.TableView | DG.DataFrame,
  columnName: string,
  criteria: FilterCriteria,
  opts: FilterRowsOpts = {},
): Promise<void> {
  const view = target instanceof DG.TableView ? target : null;
  const df = view ? view.dataFrame : (target as DG.DataFrame);

  if (criteria.substructure !== undefined) {
    await filterSubstructure(target, columnName, criteria.substructure, {
      ui: opts.ui,
      molBlockFailover: opts.molBlockFailover,
      searchType: criteria.searchType,
    });
    return;
  }

  const useUi = (opts.ui ?? true) && !!view;

  if (criteria.min !== undefined || criteria.max !== undefined) {
    if (useUi) {
      const col = df.col(columnName);
      const fg = view!.getFiltersGroup({createDefaultFilters: false});
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.HISTOGRAM,
        column: columnName,
        min: criteria.min ?? col?.min,
        max: criteria.max ?? col?.max,
      } as any);
    } else {
      const col = df.getCol(columnName);
      const lo = criteria.min ?? -Infinity;
      const hi = criteria.max ?? Infinity;
      filterByPredicate(df, (i) => {
        const v = col.get(i) as number;
        return v >= lo && v <= hi;
      });
    }
    return;
  }

  if (criteria.values !== undefined) {
    if (useUi) {
      const fg = view!.getFiltersGroup({createDefaultFilters: false});
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: columnName,
        selected: criteria.values,
      } as any);
    } else {
      const col = df.getCol(columnName);
      const set = new Set<string>(criteria.values);
      filterByPredicate(df, (i) => set.has(String(col.get(i))));
    }
    return;
  }

  if (criteria.contains !== undefined) {
    if (useUi) {
      const fg = view!.getFiltersGroup({createDefaultFilters: false});
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.FREE_TEXT,
        column: columnName,
        value: criteria.contains,
      } as any);
    } else {
      const col = df.getCol(columnName);
      const needle = criteria.contains.toLowerCase();
      filterByPredicate(df, (i) => {
        const v = col.get(i);
        return v != null && String(v).toLowerCase().includes(needle);
      });
    }
    return;
  }

  if (criteria.regex !== undefined) {
    const re = compileRegex(criteria.regex);
    const col = df.getCol(columnName);
    // Regex isn't a first-class FREE_TEXT mode shape — always use BitSet path.
    filterByPredicate(df, (i) => {
      const v = col.get(i);
      return v != null && re.test(String(v));
    });
    return;
  }

  throw new Error('filterRows: criteria must include one of {min,max} / values / contains / regex / substructure');
}

function compileRegex(input: string): RegExp {
  // Accept either '/pattern/flags' or a bare pattern.
  const m = /^\/(.+)\/([a-z]*)$/.exec(input);
  if (m) return new RegExp(m[1], m[2]);
  return new RegExp(input);
}

// ---------------------------------------------------------------------------
// filterByPredicate
// ---------------------------------------------------------------------------

/**
 * Keep rows where `pred(i)` returns true; hide the rest. Buffer-direct via
 * `df.filter.init`, single notification. Replaces (does not intersect) the
 * current filter — `init` zeros the buffer first.
 *
 * Polarity: `pred(i) === true` → row **KEPT**. Opposite of `dropRows`.
 *
 * To intersect with the existing filter, build a fresh BitSet with
 * `DG.BitSet.create(df.rowCount, ...)` and call `df.filter.and(mask)`.
 *
 * @example
 *   filterByPredicate(df, (i) => df.getCol('age').get(i) > 40);
 */
export function filterByPredicate(df: DG.DataFrame, pred: (i: number) => boolean): void {
  df.filter.init(pred);
}

// ---------------------------------------------------------------------------
// filterSubstructure
// ---------------------------------------------------------------------------

export type FilterSubstructureOpts = {
  /**
   * `true` (default) — UI path via FilterGroup widget.
   * `false` — BitSet path via `grok.chem.searchSubstructure`.
   * `'auto'` (or omitted) — UI path under 50 000 rows, BitSet above.
   */
  ui?: boolean | 'auto';
  molBlockFailover?: string;
  searchType?: 'substructure' | 'similar' | 'exact';
};

/**
 * Substructure filter on a molecule column. Accepts SMILES, SMARTS, or
 * molblock — auto-detects via `DG.chem.isMolBlock` / `DG.chem.isSmarts` and
 * converts to molblock via `DG.chem.convert` before handing to the
 * `Chem:substructureFilter` widget (which silently matches zero rows if given
 * raw SMILES/SMARTS — this is the bug the previous `filterRows` had).
 *
 * On large frames (> 50 000 rows) the helper auto-switches to the BitSet
 * path (`grok.chem.searchSubstructure`) which is materially faster but
 * doesn't show a widget. Override with `opts.ui: true | false`.
 *
 * If the Chem package isn't loaded, `console.warn` and treat the query as
 * SMILES — the widget will detect format itself once Chem loads.
 *
 * Empty query → clears the substructure filter (no error).
 *
 * @example
 *   await filterSubstructure(view, 'smiles', 'c1ccccc1');
 *   await filterSubstructure(view, 'smiles', '[#6][!#1]');
 *   await filterSubstructure(view, 'smiles', molblockString);
 */
export async function filterSubstructure(
  target: DG.TableView | DG.DataFrame,
  columnName: string,
  query: string,
  opts: FilterSubstructureOpts = {},
): Promise<void> {
  const view = target instanceof DG.TableView ? target : null;
  const df = view ? view.dataFrame : (target as DG.DataFrame);
  const col = df.getCol(columnName);

  // Empty query — clear the substructure filter rather than erroring.
  if (!query || query.trim() === '') {
    if (view) {
      const fg = view.getFiltersGroup({createDefaultFilters: false});
      // updateOrAdd with empty molBlock effectively disables the filter.
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.SUBSTRUCTURE,
        column: columnName,
        columnName: columnName,
        molBlock: '',
      } as any);
    } else {
      df.filter.setAll(true);
    }
    return;
  }

  const uiMode = opts.ui;
  const wantUi = view != null && (
    uiMode === true ||
    (uiMode === undefined || uiMode === 'auto') && df.rowCount <= LARGE_DF_THRESHOLD
  );

  if (wantUi) {
    let molBlock = query;
    try {
      if (!DG.chem.isMolBlock(query)) {
        const isSmarts = safeIsSmarts(query);
        const src = isSmarts ? DG.chem.Notation.Smarts : DG.chem.Notation.Smiles;
        molBlock = DG.chem.convert(query, src, DG.chem.Notation.MolBlock);
      }
    } catch (e) {
      console.warn(`filterSubstructure: Chem conversion failed (${(e as Error).message}); ` +
        `passing raw query — Chem:substructureFilter will attempt detection.`);
      molBlock = query;
    }
    const fg = view!.getFiltersGroup({createDefaultFilters: false});
    const state: any = {
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: columnName,
      columnName: columnName,
      molBlock,
    };
    if (opts.searchType) state.searchType = opts.searchType;
    if (opts.molBlockFailover) state.molBlockFailover = opts.molBlockFailover;
    fg.updateOrAdd(state);
    return;
  }

  // BitSet path — works for any pattern shape, no widget produced.
  const bs = await grok.chem.searchSubstructure(col, query, opts.molBlockFailover ? {molBlockFailover: opts.molBlockFailover} : {});
  df.filter.copyFrom(bs);
}

function safeIsSmarts(s: string): boolean {
  try {
    return DG.chem.isSmarts(s);
  } catch {
    console.warn('filterSubstructure: DG.chem.isSmarts unavailable (Chem package not loaded?) — treating query as SMILES.');
    return false;
  }
}

// ---------------------------------------------------------------------------
// combineFilters
// ---------------------------------------------------------------------------

/**
 * AND or OR several predicates in a single pass over rows, then assign onto
 * `df.filter`. Faster than building N BitSets and folding them in.
 *
 * Polarity matches `filterByPredicate`: each predicate returns `true` to
 * **KEEP** the row. With `op: 'and'`, all predicates must say keep. With
 * `op: 'or'`, any one is enough.
 *
 * Throws if no predicates are passed — that intent is ambiguous (clear?
 * set-all-false?). Use `clearFilter` / `df.filter.setAll(false)` directly.
 *
 * @example
 *   combineFilters(df, 'and',
 *     (i) => df.getCol('MW').get(i) < 500,
 *     (i) => df.getCol('cLogP').get(i) < 5);
 */
export function combineFilters(
  df: DG.DataFrame,
  op: 'and' | 'or',
  ...predicates: Array<(i: number) => boolean>
): void {
  if (predicates.length === 0)
    throw new Error('combineFilters: need at least one predicate');
  if (op === 'and') {
    df.filter.init((i) => {
      for (const p of predicates) if (!p(i)) return false;
      return true;
    });
  } else {
    df.filter.init((i) => {
      for (const p of predicates) if (p(i)) return true;
      return false;
    });
  }
}

// ---------------------------------------------------------------------------
// clearFilter
// ---------------------------------------------------------------------------

/**
 * Show every row. Polymorphic:
 * - `target: DG.DataFrame` → `df.filter.setAll(true)` only.
 * - `target: DG.TableView` → also disables every UI filter via
 *   `fg.setActive(false)`, so widget contributions don't immediately
 *   re-narrow the now-cleared mask. Widgets stay on screen (just disabled).
 *
 * @example
 *   clearFilter(view);   // disables UI filters + clears mask
 *   clearFilter(df);     // mask only
 */
export function clearFilter(target: DG.DataFrame | DG.TableView): void {
  if (target instanceof DG.TableView) {
    const fg = target.getFiltersGroup({createDefaultFilters: false});
    fg.setActive(false);
    target.dataFrame.filter.setAll(true);
  } else {
    target.filter.setAll(true);
  }
}

// ---------------------------------------------------------------------------
// invertFilter
// ---------------------------------------------------------------------------

/**
 * Flip the current filter mask in place. Same caveat as `filterByPredicate`:
 * if UI filters are active, your inversion is clobbered on the next
 * `onRowsFiltering` cycle — `clearFilter(view)` first if you want it to stick.
 *
 * @example
 *   invertFilter(df);
 */
export function invertFilter(df: DG.DataFrame): void {
  df.filter.invert();
}

// ---------------------------------------------------------------------------
// filteredDf
// ---------------------------------------------------------------------------

export type FilteredDfOpts = {
  /** Narrow to a column subset. */
  cols?: string[];
  /** Carry the selection mask onto the clone. */
  withSelection?: boolean;
};

/**
 * Non-destructive — returns a clone of the currently-visible rows. The
 * source DataFrame is not modified.
 *
 * Distinct from `dropRows` (destructive) and from `filterByPredicate`
 * (mutates `df.filter`, no clone).
 *
 * @example
 *   const subset = filteredDf(df);
 *   const small = filteredDf(df, {cols: ['smiles', 'activity']});
 */
export function filteredDf(df: DG.DataFrame, opts: FilteredDfOpts = {}): DG.DataFrame {
  return df.clone(df.filter, opts.cols ?? null, opts.withSelection ?? false);
}

// ---------------------------------------------------------------------------
// dropRows
// ---------------------------------------------------------------------------

/**
 * **Destructive** — removes rows where `pred(i)` returns true. Rows are gone
 * from the DataFrame; this is **not** the same as filtering.
 *
 * Polarity: `pred(i) === true` → row **DROPPED**. Opposite of
 * `filterByPredicate`, which keeps rows where the predicate is true. Pick
 * carefully or you'll delete the wrong rows.
 *
 * @returns The number of rows actually removed.
 *
 * @example
 *   const removed = dropRows(df, (i) => df.getCol('flag').get(i) === 'bad');
 */
export function dropRows(df: DG.DataFrame, pred: (i: number) => boolean): number {
  const before = df.rowCount;
  df.rows.removeWhereIdx(pred);
  return before - df.rowCount;
}
