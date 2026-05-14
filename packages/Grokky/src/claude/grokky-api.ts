import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type {
  ColumnDescription, CloneDfOpts, ColumnMetaSpec,
  AddColumnSpec, FilterCriteria, FilterRowsOpts, FilteredDfOpts, FilterSubstructureOpts,
  SelectMode, SelectInput, SelectRowsOpts, SelectedDfOpts, SelectionDescription,
  SortOrder, ConfigureGridOpts, ColorCodeSpec, ResetGridOpts, JoinOptions,
} from './grokky-api-types';
import {LARGE_DF_THRESHOLD, MAX_INDEXES_IN_DESCRIBE, MAX_SAMPLE_COLUMNS} from './grokky-api-types';

// ─── DataFrame & columns ──────────────────────────────────────────────────

export function describeColumn(col: DG.Column): ColumnDescription {
  const out: ColumnDescription = {
    name: col.name,
    type: col.type,
    semType: col.semType || null,
    length: col.length,
    missing: col.stats.missingValueCount,
    unique: col.stats.uniqueCount,
  };
  const units = col.meta.units;
  if (units) out.units = units;
  const format = col.meta.format;
  if (format) out.format = format;
  const friendlyName = col.meta.friendlyName;
  if (friendlyName) out.friendlyName = friendlyName;
  const description = col.meta.description;
  if (description) out.description = description;
  if (col.isNumerical) {
    const s = col.stats;
    out.numerical = {min: s.min, max: s.max, avg: s.avg, stdev: s.stdev, med: s.med};
  }
  if (col.type === DG.COLUMN_TYPE.STRING)
    out.categorical = {topCategories: topCategories(col, 5)};
  return out;
}

export function cloneDf(df: DG.DataFrame, opts: CloneDfOpts = {}): DG.DataFrame {
  let mask: DG.BitSet | null = null;
  if (opts.rows === 'filtered')
    mask = df.filter;
  else if (opts.rows === 'selected')
    mask = df.selection;
  else if (opts.rows instanceof DG.BitSet)
    mask = opts.rows;
  return df.clone(mask, opts.cols ?? null, opts.withSelection ?? false, opts.withTags ?? true);
}

export function setColumnMeta(col: DG.Column, meta: ColumnMetaSpec): DG.Column {
  if (meta.semType !== undefined)
    col.semType = meta.semType as string;
  if (meta.format !== undefined)
    col.meta.format = meta.format;
  if (meta.units !== undefined)
    col.meta.units = meta.units;
  if (meta.friendlyName !== undefined)
    col.meta.friendlyName = meta.friendlyName;
  if (meta.description !== undefined)
    col.meta.description = meta.description;
  if (meta.choices !== undefined)
    col.meta.choices = meta.choices;

  if (meta.colorCoding !== undefined)
    colorCode(col, meta.colorCoding);

  if (meta.tags) {
    for (const [k, v] of Object.entries(meta.tags)) {
      if (v === null)
        col.tags.delete(k);
      else col.setTag(k, v);
    }
  }
  return col;
}

export async function addColumn(df: DG.DataFrame, spec: AddColumnSpec): Promise<DG.Column> {
  const sources = [spec.formula, spec.values, spec.init, spec.virtual].filter((s) => s !== undefined);
  if (sources.length > 1)
    throw new Error('addColumn: only one of formula / values / init / virtual may be set');

  const name = (spec.ensureUniqueName ?? true) ? df.columns.getUnusedName(spec.name) : spec.name;
  let col: DG.Column;

  if (spec.formula !== undefined)
    col = await df.columns.addNewCalculated(name, spec.formula, spec.type ?? 'auto');
  else if (spec.values !== undefined) {
    const t = spec.type ?? inferTypeFromValues(spec.values);
    const built = DG.Column.fromList(t, name, spec.values);
    col = df.columns.add(built);
  } else if (spec.virtual !== undefined)
    col = df.columns.addNewVirtual(name, spec.virtual, (spec.type ?? DG.COLUMN_TYPE.OBJECT) as any);
  else if (spec.init !== undefined) {
    col = df.columns.addNew(name, spec.type ?? DG.COLUMN_TYPE.FLOAT);
    col.init(spec.init);
  } else
    col = df.columns.addNew(name, spec.type ?? DG.COLUMN_TYPE.STRING);


  if (spec.insertAt !== undefined && spec.insertAt >= 0) {
    df.columns.remove(col);
    df.columns.insert(col, spec.insertAt);
  }

  if (spec.meta)
    setColumnMeta(col, spec.meta);
  return Promise.resolve(col);
}

function inferTypeFromValues(values: any[]): DG.ColumnType {
  for (const v of values) {
    if (v === null || v === undefined)
      continue;
    if (typeof v === 'number')
      return Number.isInteger(v) ? DG.COLUMN_TYPE.INT : DG.COLUMN_TYPE.FLOAT;
    if (typeof v === 'boolean')
      return DG.COLUMN_TYPE.BOOL;
    if (v instanceof Date)
      return DG.COLUMN_TYPE.DATE_TIME;
    return DG.COLUMN_TYPE.STRING;
  }
  return DG.COLUMN_TYPE.STRING;
}

export function removeColumns(
  df: DG.DataFrame,
  names: string[],
  onMissing: 'skip' | 'throw' | 'warn' = 'skip',
): string[] {
  const removed: string[] = [];
  for (const name of names) {
    const col = df.col(name);
    if (!col) {
      if (onMissing === 'throw')
        throw new Error(`removeColumns: column "${name}" not found`);
      if (onMissing === 'warn')
        console.warn(`removeColumns: column "${name}" not found, skipping`);
      continue;
    }
    df.columns.remove(col);
    removed.push(col.name);
  }
  return removed;
}

export function renameColumn(
  df: DG.DataFrame,
  from: string,
  to: string,
  opts: {ensureUnique?: boolean} = {},
): string {
  const col = df.getCol(from);
  const finalName = opts.ensureUnique ? df.columns.getUnusedName(to) : to;
  if (finalName === col.name)
    return finalName;
  try {
  } catch { /* best-effort */ }
  col.name = finalName;
  return finalName;
}

export function topCategories(col: DG.Column, n: number = 5): Array<{value: string; count: number}> {
  if (col.type !== DG.COLUMN_TYPE.STRING)
    return [];
  const counts = new Map<string, number>();
  const len = col.length;
  for (let i = 0; i < len; i++) {
    const v = col.get(i) as string | null;
    if (v === null || v === undefined)
      continue;
    counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  return Array.from(counts.entries())
    .sort((a, b) => b[1] - a[1])
    .slice(0, n)
    .map(([value, count]) => ({value, count}));
}

// ─── Calculated columns (preserved) ───────────────────────────────────────

export async function addCalculatedColumn(
  df: DG.DataFrame, name: string, formula: string, type: DG.ColumnType | 'auto' = 'auto',
): Promise<DG.Column> {
  const col = await df.columns.addNewCalculated(name, formula, type);
  return col;
}

// ─── Filtering ────────────────────────────────────────────────────────────

function resolveFilterTarget(target: DG.DataFrame | DG.TableView): {df: DG.DataFrame; view: DG.TableView | null} {
  if (target instanceof DG.TableView)
    return {df: target.dataFrame, view: target};
  const tv = grok.shell.tv;
  const view = tv && tv.dataFrame === target ? tv : null;
  return {df: target as DG.DataFrame, view};
}

// onRowsFiltering subscription keeps the predicate alive across UI filter cycles;
// view.subs auto-unsubscribes on view detach.
function installCollaborativeFilter(
  df: DG.DataFrame, view: DG.TableView, apply: (bs: DG.BitSet) => void,
): void {
  const sub = df.onRowsFiltering.subscribe(() => apply(df.filter));
  view.subs.push(sub);
  df.rows.requestFilter();
}

export async function filterRows(
  target: DG.TableView | DG.DataFrame,
  columnName: string,
  criteria: FilterCriteria,
  opts: FilterRowsOpts = {},
): Promise<void> {
  const view = target instanceof DG.TableView ? target : (grok.shell.tv ?? null);
  if (!view)
    throw new Error('filterRows: no TableView available. Pass a TableView, or open a table view first.');
  const df = view.dataFrame;
  const fg = view.getFiltersGroup({createDefaultFilters: false});

  if (criteria.substructure !== undefined) {
    await filterSubstructure(view, columnName, criteria.substructure, {
      molBlockFailover: opts.molBlockFailover,
      searchType: criteria.searchType,
    });
    return;
  }

  if (criteria.min !== undefined || criteria.max !== undefined) {
    const col = df.col(columnName);
    // FILTER_TYPE.HISTOGRAM is the min/max shape; do not swap with CATEGORICAL.
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.HISTOGRAM,
      column: columnName,
      min: criteria.min ?? col?.min,
      max: criteria.max ?? col?.max,
    } as any);
    return;
  }

  if (criteria.values !== undefined) {
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL,
      column: columnName,
      selected: criteria.values,
    } as any);
    return;
  }

  if (criteria.contains !== undefined) {
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.FREE_TEXT,
      column: columnName,
      value: criteria.contains,
    } as any);
    return;
  }

  if (criteria.regex !== undefined) {
    // Regex isn't a first-class FREE_TEXT mode shape — fall back to a collaborative predicate.
    const re = compileRegex(criteria.regex);
    const col = df.getCol(columnName);
    filterByPredicate(view, (i) => {
      const v = col.get(i);
      return v != null && re.test(String(v));
    });
    return;
  }

  throw new Error('filterRows: criteria must include one of {min,max} / values / contains / regex / substructure');
}

function compileRegex(input: string): RegExp {
  const m = /^\/(.+)\/([a-z]*)$/.exec(input);
  if (m)
    return new RegExp(m[1], m[2]);
  return new RegExp(input);
}

/**
 * Keep rows where `pred(i)` returns true; hide the rest. When given a TableView
 * (or when grok.shell.tv binds to the df), subscribes via onRowsFiltering so the
 * predicate persists across UI filter cycles. When given only a DataFrame, applies
 * `df.filter.init(pred)` once — gets wiped on the next UI filter cycle.
 */
export function filterByPredicate(target: DG.DataFrame | DG.TableView, pred: (i: number) => boolean): void {
  const {df, view} = resolveFilterTarget(target);
  if (view) {
    installCollaborativeFilter(df, view, (bs) => bs.init(pred));
    return;
  }
  df.filter.init(pred);
}

export function clearFilter(target: DG.DataFrame | DG.TableView): void {
  if (target instanceof DG.TableView) {
    const fg = target.getFiltersGroup({createDefaultFilters: false});
    fg.setActive(false);
    target.dataFrame.filter.setAll(true);
  } else
    target.filter.setAll(true);
}

/**
 * Flip the current filter mask. With a TableView, installs a collaborative handler
 * that inverts on every filtering pass. With a DataFrame alone, a one-shot invert
 * that the next UI filter cycle will clobber.
 */
export function invertFilter(target: DG.DataFrame | DG.TableView): void {
  const {df, view} = resolveFilterTarget(target);
  if (view) {
    installCollaborativeFilter(df, view, (bs) => bs.invert());
    return;
  }
  df.filter.invert();
}

export function combineFilters(
  target: DG.DataFrame | DG.TableView,
  op: 'and' | 'or',
  ...predicates: Array<(i: number) => boolean>
): void {
  if (predicates.length === 0)
    throw new Error('combineFilters: need at least one predicate');
  const merged = op === 'and' ?
    (i: number) => {
      for (const p of predicates) {
        if (!p(i))
          return false;
      }
      return true;
    } :
    (i: number) => {
      for (const p of predicates) {
        if (p(i))
          return true;
      }
      return false;
    };
  const {df, view} = resolveFilterTarget(target);
  if (view) {
    installCollaborativeFilter(df, view, (bs) => bs.init(merged));
    return;
  }
  df.filter.init(merged);
}

export function filteredDf(df: DG.DataFrame, opts: FilteredDfOpts = {}): DG.DataFrame {
  return df.clone(df.filter, opts.cols ?? null, opts.withSelection ?? false);
}

/**
 * Destructive — removes rows where `pred(i)` returns true. Opposite polarity from
 * `filterByPredicate`. Accepts a view only to keep the call-site shape uniform;
 * row removal does not need the collaborative pattern.
 */
export function dropRows(target: DG.DataFrame | DG.TableView, pred: (i: number) => boolean): number {
  const df = target instanceof DG.TableView ? target.dataFrame : target;
  const before = df.rowCount;
  df.rows.removeWhereIdx(pred);
  return before - df.rowCount;
}

export async function filterSubstructure(
  target: DG.TableView | DG.DataFrame,
  columnName: string,
  query: string,
  opts: FilterSubstructureOpts = {},
): Promise<void> {
  const view = target instanceof DG.TableView ? target : null;
  const df = view ? view.dataFrame : (target as DG.DataFrame);
  const col = df.getCol(columnName);

  if (!query || query.trim() === '') {
    if (view) {
      const fg = view.getFiltersGroup({createDefaultFilters: false});
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.SUBSTRUCTURE,
        column: columnName,
        columnName: columnName,
        molBlock: '',
      } as any);
    } else
      df.filter.setAll(true);

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
    if (opts.searchType)
      state.searchType = opts.searchType;
    if (opts.molBlockFailover)
      state.molBlockFailover = opts.molBlockFailover;
    fg.updateOrAdd(state);
    return;
  }

  const ssOpts = opts.molBlockFailover ? {molBlockFailover: opts.molBlockFailover} : {};
  const bs = await grok.chem.searchSubstructure(col, query, ssOpts);
  df.filter.copyFrom(bs);
}

function safeIsSmarts(s: string): boolean {
  try {
    return DG.chem.isSmarts(s);
  } catch {
    console.warn('filterSubstructure: DG.chem.isSmarts unavailable ' +
      '(Chem package not loaded?) — treating query as SMILES.');
    return false;
  }
}

// ─── Selection ────────────────────────────────────────────────────────────

export function selectRows(
  df: DG.DataFrame,
  input: SelectInput,
  opts: SelectRowsOpts = {},
): void {
  const mode: SelectMode = opts.mode ?? 'replace';
  const sel = df.selection;

  if (input instanceof DG.BitSet) {
    switch (mode) {
    case 'replace': sel.copyFrom(input); return;
    case 'add': sel.or(input); return;
    case 'remove': sel.andNot(input); return;
    case 'intersect': sel.and(input); return;
    }
  }

  if (typeof input === 'function') {
    if (mode === 'replace') {
      sel.init(input);
      return;
    }
    const mask = DG.BitSet.create(df.rowCount, input);
    switch (mode) {
    case 'add': sel.or(mask); return;
    case 'remove': sel.andNot(mask); return;
    case 'intersect': sel.and(mask); return;
    }
  }

  const indices = input as ArrayLike<number>;
  const len = indices.length;

  if (mode === 'replace') {
    sel.setAll(false, false);
    for (let k = 0; k < len; k++) sel.set(indices[k], true, false);
    sel.fireChanged();
    return;
  }
  if (mode === 'add') {
    if (len === 0)
      return;
    for (let k = 0; k < len; k++) sel.set(indices[k], true, false);
    sel.fireChanged();
    return;
  }
  if (mode === 'remove') {
    if (len === 0)
      return;
    for (let k = 0; k < len; k++) sel.set(indices[k], false, false);
    sel.fireChanged();
    return;
  }
  if (len === 0) {
    sel.setAll(false);
    return;
  }
  const keep = new Set<number>();
  for (let k = 0; k < len; k++) keep.add(indices[k]);
  sel.init((i) => keep.has(i) && sel.get(i));
}

// clearSelection / clearFilter polarity is the trap to remember:
//   clearFilter  → setAll(TRUE)   — every row visible
//   clearSelection → setAll(FALSE) — none selected
export function clearSelection(df: DG.DataFrame): void {
  df.selection.setAll(false);
}

export function invertSelection(df: DG.DataFrame): void {
  df.selection.invert();
}

export function selectAll(df: DG.DataFrame): void {
  df.selection.setAll(true);
}

export function selectedDf(df: DG.DataFrame, opts: SelectedDfOpts = {}): DG.DataFrame {
  return df.clone(df.selection, opts.cols ?? null, opts.withSelection ?? false);
}

export function filterFromSelection(df: DG.DataFrame): void {
  df.filter.copyFrom(df.selection);
}

export function selectionFromFilter(df: DG.DataFrame): void {
  df.selection.copyFrom(df.filter);
}

export function setCurrentRow(df: DG.DataFrame, idx: number): void {
  if (!Number.isInteger(idx))
    throw new Error(`setCurrentRow: idx must be an integer, got ${idx}`);
  if (idx < 0 || idx >= df.rowCount)
    throw new Error(`setCurrentRow: idx ${idx} out of range [0, ${df.rowCount})`);
  df.currentRowIdx = idx;
}

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

// ─── Viewers ──────────────────────────────────────────────────────────────

const VIEWER_TYPES: string[] = Object.values(DG.VIEWER);

function canonicalizeViewerType(input: string): string {
  const norm = input.toLowerCase().trim().replace(/\s+/g, ' ');
  const exact = VIEWER_TYPES.find((t) => t.toLowerCase() === norm);
  return exact ?? input;
}

function applyViewerOptions(viewer: DG.Viewer, options: Record<string, any>): void {
  const knownSet = new Set<string>(viewer.getProperties().map((p: any) => p.name));
  const accepted: Record<string, any> = {};
  const unknown: string[] = [];
  for (const [key, val] of Object.entries(options)) {
    if (knownSet.has(key))
      accepted[key] = val;
    else if (knownSet.has(key + 'ColumnName'))
      accepted[key + 'ColumnName'] = val;
    else if (knownSet.has(key + 'ColumnNames'))
      accepted[key + 'ColumnNames'] = val;
    else
      unknown.push(key);
  }
  if (Object.keys(accepted).length > 0)
    viewer.setOptions(accepted);
  if (unknown.length > 0) {
    const quoted = unknown.map((k) => `"${k}"`).join(', ');
    console.warn(`addViewer "${viewer.type}": unknown ${unknown.length === 1 ? 'property' : 'properties'}: ${quoted}`);
  }
}

function resolveTableView(viewArg: DG.TableView | DG.ViewBase | null | undefined): DG.TableView {
  if (viewArg && viewArg instanceof DG.TableView)
    return viewArg;
  if (viewArg) {
    console.warn(
      `grokky viewer helper: \`view\` is a ${viewArg?.constructor?.name ?? 'non-TableView'}; ` +
      `falling back to grok.shell.tv. Pass a TableView explicitly for predictable scoping.`);
  }
  const tv = grok.shell.tv;
  if (!tv)
    throw new Error('grokky viewer helper: no active TableView (grok.shell.tv is null).');
  return tv;
}

export function addViewer(
  viewArg: DG.TableView | DG.ViewBase | null,
  type: string,
  options: Record<string, any> = {},
): DG.Viewer {
  const tv = resolveTableView(viewArg);
  const canonicalType = canonicalizeViewerType(type);
  const v = tv.addViewer(canonicalType);
  applyViewerOptions(v, options);
  return v;
}

export function configureViewer(viewer: DG.Viewer, options: Record<string, any>): void {
  applyViewerOptions(viewer, options);
}

function viewerSlice(view: DG.TableView | null | undefined): DG.Viewer[] {
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grokky viewer helper: findViewer/findViewers requires a TableView.');
  // view.viewers[0] is the grid — skip it.
  return Array.from(view.viewers).slice(1);
}

export function findViewer(
  view: DG.TableView | null,
  pred: string | ((v: DG.Viewer) => boolean),
): DG.Viewer | null {
  const candidates = viewerSlice(view);
  const fn = typeof pred === 'string' ?
    (v: DG.Viewer) => v.type === canonicalizeViewerType(pred) :
    pred;
  return candidates.find(fn) ?? null;
}

export function findViewers(
  view: DG.TableView | null,
  pred?: (v: DG.Viewer) => boolean,
): DG.Viewer[] {
  const candidates = viewerSlice(view);
  return pred ? candidates.filter(pred) : candidates;
}

function safeClose(v: DG.Viewer): boolean {
  try {
    v.close();
    return true;
  } catch (e) {
    console.warn(`grokky.closeViewer: viewer.close() threw (${(e as Error).message}); ` +
      `swallowing — the viewer was likely never attached.`);
    return false;
  }
}

export function closeViewer(
  target: DG.Viewer | string | ((v: DG.Viewer) => boolean),
  view?: DG.TableView | null,
): number {
  if (target instanceof DG.Viewer)
    return safeClose(target) ? 1 : 0;

  if (!view || !(view instanceof DG.TableView))
    throw new Error('grokky.closeViewer: string/predicate target requires a TableView as second arg.');

  const fn = typeof target === 'string' ?
    (v: DG.Viewer) => v.type === canonicalizeViewerType(target) :
    target;
  const candidates = viewerSlice(view).filter(fn);
  let count = 0;
  for (const v of candidates) {
    if (safeClose(v))
      count++;
  }
  return count;
}

export function closeAllViewers(
  view: DG.TableView | null,
  opts: {keepGrid?: boolean} = {},
): number {
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grokky.closeAllViewers: requires a TableView.');
  const keepGrid = opts.keepGrid ?? true;
  const all = Array.from(view.viewers);
  const targets = keepGrid ? all.slice(1) : all;
  let count = 0;
  for (const v of targets) {
    if (safeClose(v))
      count++;
  }
  return count;
}

// ─── Grid customization ───────────────────────────────────────────────────

function requireTableView(view: any): DG.TableView {
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grid customization requires a TableView (got: ' + (view?.type ?? typeof view) + ')');
  return view as DG.TableView;
}

function gridCol(view: DG.TableView, name: string): DG.GridColumn | null {
  return view.grid.col(name);
}

/**
 * Additive by default — appends to the current grid sort, overriding direction
 * for columns that were already in the sort. Pass `{replace: true}` to clear
 * the previous sort first. `view.grid.sortByColumns = [...]` is a getter-only
 * assignment that silently does nothing — always call `grid.sort()`.
 */
export function applySort(
  view: DG.TableView | any,
  columns: string[],
  orders?: SortOrder[],
  opts: {replace?: boolean} = {},
): void {
  const tv = requireTableView(view);
  if (!Array.isArray(columns))
    throw new Error('applySort: "columns" must be a string[]');
  if (orders && orders.length !== columns.length && orders.length !== 0)
    throw new Error(`applySort: "orders" length (${orders.length}) must match "columns" length (${columns.length})`);

  const newCols: string[] = [];
  const newBools: boolean[] = [];
  columns.forEach((c, i) => {
    if (!tv.dataFrame.col(c)) {
      console.warn(`applySort: column "${c}" not found on dataFrame; skipping`);
      return;
    }
    const o = orders?.[i] ?? 'asc';
    newCols.push(c);
    newBools.push(typeof o === 'boolean' ? o : o !== 'desc');
  });

  if (opts.replace) {
    tv.grid.sort(newCols, newBools);
    return;
  }

  const merged = new Map<string, boolean>();
  const prevCols = (tv.grid.sortByColumns ?? []).map((c) => c.name);
  const prevBools = Array.from(tv.grid.sortTypes ?? []);
  prevCols.forEach((c, i) => merged.set(c, prevBools[i] ?? true));
  newCols.forEach((c, i) => merged.set(c, newBools[i]));
  tv.grid.sort(Array.from(merged.keys()), Array.from(merged.values()));
}

export function clearSort(view: DG.TableView | any): void {
  const tv = requireTableView(view);
  tv.grid.sort([], []);
}

export function configureGrid(view: DG.TableView | any, opts: ConfigureGridOpts): void {
  const tv = requireTableView(view);
  const grid = tv.grid;

  if (opts.hide) {
    for (const name of opts.hide) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: hide: column "${name}" not found, skipping`); continue; }
      gc.visible = false;
    }
  }
  if (opts.show) {
    const known = opts.show.filter((n) => {
      if (gridCol(tv, n))
        return true;
      console.warn(`configureGrid: show: column "${n}" not found, skipping`);
      return false;
    });
    grid.columns.setVisible(known);
  }

  if (opts.order) {
    const known = opts.order.filter((n) => {
      if (gridCol(tv, n))
        return true;
      console.warn(`configureGrid: order: column "${n}" not found, skipping`);
      return false;
    });
    if (known.length > 0)
      grid.columns.setOrder(known);
  }

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

  if (opts.widthPolicy) {
    const enumVal = (DG as any).ColumnWidthType?.[opts.widthPolicy];
    if (enumVal === undefined) {
      console.warn(`configureGrid: widthPolicy "${opts.widthPolicy}" not recognized ` +
        '(expected Minimal | Compact | Optimal | Maximal)');
    } else
      grid.setColumnsWidthType(enumVal);
  }

  if (opts.pin) {
    for (const name of opts.pin) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: pin: column "${name}" not found, skipping`); continue; }
      try { gc.pin(); } catch (e) { console.warn(`configureGrid: pin "${name}": ${(e as Error).message}`); }
    }
  }
  if (opts.unpin) {
    for (const name of opts.unpin) {
      const gc = gridCol(tv, name);
      if (!gc) { console.warn(`configureGrid: unpin: column "${name}" not found, skipping`); continue; }
      try { gc.unpin(); } catch (e) { console.warn(`configureGrid: unpin "${name}": ${(e as Error).message}`); }
    }
  }

  // Data-side: col.meta.format persists across viewers.
  if (opts.formats) {
    for (const [name, format] of Object.entries(opts.formats)) {
      const col = tv.dataFrame.col(name);
      if (!col) { console.warn(`configureGrid: formats: column "${name}" not found, skipping`); continue; }
      col.meta.format = format;
    }
  }

  if (opts.rowHeight !== undefined) {
    if (typeof opts.rowHeight !== 'number' || !Number.isFinite(opts.rowHeight) || opts.rowHeight <= 0)
      console.warn(`configureGrid: rowHeight ${opts.rowHeight} is not a positive number, skipping`);
    else
      grid.setOptions({rowHeight: opts.rowHeight});
  }

  if (opts.frozenColumns !== undefined) {
    if (typeof opts.frozenColumns !== 'number' || !Number.isInteger(opts.frozenColumns) || opts.frozenColumns < 0)
      console.warn(`configureGrid: frozenColumns ${opts.frozenColumns} is not a non-negative integer, skipping`);
    else
      grid.setOptions({frozenColumns: opts.frozenColumns});
  }
}

export function colorCode(col: DG.Column, spec: ColorCodeSpec): void {
  if (!col)
    throw new Error('colorCode: column is null or undefined');

  if (spec.kind === 'linear') {
    if (!col.isNumerical)
      throw new Error(`colorCode: linear color coding requires a numerical column; "${col.name}" is ${col.type}`);
    const opts: {min?: number; max?: number; belowMinColor?: string; aboveMaxColor?: string} = {};
    if (spec.min !== undefined)
      opts.min = spec.min;
    if (spec.max !== undefined)
      opts.max = spec.max;
    if (spec.belowMinColor !== undefined)
      opts.belowMinColor = spec.belowMinColor;
    if (spec.aboveMaxColor !== undefined)
      opts.aboveMaxColor = spec.aboveMaxColor;
    (col.meta.colors as any).setLinear(spec.range ?? null, Object.keys(opts).length ? opts : null);
    return;
  }

  if (spec.kind === 'categorical') {
    const opts: {fallbackColor?: string | number; matchType?: string} = {};
    if (spec.fallbackColor !== undefined)
      opts.fallbackColor = spec.fallbackColor;
    if (spec.matchType !== undefined)
      opts.matchType = spec.matchType;
    (col.meta.colors as any).setCategorical(spec.colors as any, Object.keys(opts).length ? (opts as any) : null);
    return;
  }

  if (spec.kind === 'conditional') {
    if (!spec.rules || Object.keys(spec.rules).length === 0)
      throw new Error('colorCode: conditional spec requires non-empty `rules`');
    col.meta.colors.setConditional(spec.rules);
    return;
  }

  if (spec.kind === 'off') {
    col.meta.colors.setDisabled();
    return;
  }

  const exhaustive: never = spec;
  throw new Error(`colorCode: unknown spec ${JSON.stringify(exhaustive)}`);
}

export function resetGrid(view: DG.TableView | any, opts: ResetGridOpts = {}): void {
  const tv = requireTableView(view);
  const grid = tv.grid;
  const o: Required<ResetGridOpts> = {
    visibility: opts.visibility ?? true,
    widths: opts.widths ?? true,
    sort: opts.sort ?? true,
    colors: opts.colors ?? true,
  };

  if (o.visibility) {
    // Skip index 0 (the row header).
    const len = grid.columns.length;
    for (let i = 1; i < len; i++) {
      const gc = grid.columns.byIndex(i);
      if (gc)
        gc.visible = true;
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
      try {
        col.meta.colors.setDisabled();
      } catch (e) {
        console.warn(`resetGrid: setDisabled on "${col.name}": ${(e as Error).message}`);
      }
    }
  }
}

export function pinColumn(view: DG.TableView | any, name: string): void {
  const tv = requireTableView(view);
  const gc = gridCol(tv, name);
  if (!gc) {
    console.warn(`pinColumn: column "${name}" not found, skipping`);
    return;
  }
  gc.pin();
}

export function unpinColumn(view: DG.TableView | any, name: string): void {
  const tv = requireTableView(view);
  const gc = gridCol(tv, name);
  if (!gc) {
    console.warn(`unpinColumn: column "${name}" not found, skipping`);
    return;
  }
  gc.unpin();
}

// ─── Aggregations & joins (preserved) ─────────────────────────────────────

export function aggregateBy(
  df: DG.DataFrame, groupCols: string[], aggMap: Record<string, DG.AggregationType>,
): DG.DataFrame {
  const b = df.groupBy(groupCols);
  for (const [col, agg] of Object.entries(aggMap)) b.add(agg, col);
  return b.aggregate();
}

export function pivot(
  df: DG.DataFrame,
  options: {rows: string[], cols: string, value: string, agg: DG.AggregationType},
): DG.DataFrame {
  return df.groupBy(options.rows).pivot(options.cols).add(options.agg, options.value).aggregate();
}

export function unpivot(
  df: DG.DataFrame, options: {idCols: string[], valueCols: string[]},
): DG.DataFrame {
  return df.unpivot(options.idCols, options.valueCols);
}

export async function joinTables(left: DG.DataFrame, right: DG.DataFrame, opts: JoinOptions): Promise<DG.DataFrame> {
  let lKey = opts.leftKey;
  let rKey = opts.rightKey;
  if (opts.normalizeKey) {
    lKey = `__norm_${opts.leftKey}`;
    rKey = `__norm_${opts.rightKey}`;
    await normalizeInto(left, opts.leftKey, lKey, opts.normalizeKey);
    await normalizeInto(right, opts.rightKey, rKey, opts.normalizeKey);
  }
  return grok.data.joinTables(left, right, [lKey], [rKey], null, null, opts.joinType ?? DG.JOIN_TYPE.INNER, false);
}

async function normalizeInto(
  df: DG.DataFrame, src: string, dst: string,
  kind: NonNullable<JoinOptions['normalizeKey']>,
): Promise<void> {
  const formula =
    kind === 'canonicalSmiles' ? `Chem:canonicalize(\${${src}})` :
      kind === 'lowercase' ? `ToLowerCase(\${${src}})` :
        `Trim(\${${src}})`;
  await df.columns.addNewCalculated(dst, formula);
}

export const grokky = {
  // df & columns
  describeColumn, cloneDf, addColumn, setColumnMeta, removeColumns, renameColumn, topCategories,
  addCalculatedColumn,
  // filter
  filterRows, filterByPredicate, clearFilter, invertFilter, combineFilters, filteredDf, dropRows, filterSubstructure,
  // selection
  selectRows, clearSelection, invertSelection, selectAll, selectedDf,
  filterFromSelection, selectionFromFilter, setCurrentRow, describeSelection,
  // viewers
  addViewer, configureViewer, findViewer, findViewers, closeViewer, closeAllViewers,
  // grid
  applySort, clearSort, configureGrid, colorCode, resetGrid, pinColumn, unpinColumn,
  // aggregations & joins
  aggregateBy, pivot, unpivot, joinTables,
};
