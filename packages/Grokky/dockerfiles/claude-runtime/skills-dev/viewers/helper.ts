/**
 * Draft `grokky.*` helpers for the `datagrok-viewers` skill.
 *
 * These will eventually replace `addViewer` / `configureViewer` in
 * `src/claude/grokky-api.ts` and add the `findViewer` / `findViewers` /
 * `closeViewer` / `closeAllViewers` helpers. Living here keeps the skill
 * review self-contained: SKILL.md cites these by name, examples.ts type-checks
 * against them, eval/prompts.json references them.
 *
 * Conventions:
 *   - `view.viewers[0]` is ALWAYS the grid. `findViewer*` and `closeAllViewers`
 *     skip it. Closing it breaks the view.
 *   - `viewer.close()` on a never-attached viewer **throws**. `closeViewer`
 *     catches and logs — never propagates.
 *   - Type strings fuzzy-match against `DG.VIEWER` via Levenshtein ≤ 3.
 *   - Properties validate against the live `viewer.getProperties()` schema;
 *     unknown keys → console.warn with a did-you-mean hint.
 *   - View resolution: prefer in-scope `view` (TableView); fall back to
 *     `grok.shell.tv` with a warning; throw only if no TableView exists.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

// ---------------------------------------------------------------------------
// Fuzzy-matching primitives (mirrors src/claude/grokky-api.ts; TODO unify)
// ---------------------------------------------------------------------------

const VIEWER_TYPES: string[] = Object.values(DG.VIEWER);

function levenshtein(a: string, b: string): number {
  if (a === b) return 0;
  if (!a.length) return b.length;
  if (!b.length) return a.length;
  const dp = Array.from({length: a.length + 1}, () => new Array<number>(b.length + 1).fill(0));
  for (let i = 0; i <= a.length; i++) dp[i][0] = i;
  for (let j = 0; j <= b.length; j++) dp[0][j] = j;
  for (let i = 1; i <= a.length; i++) {
    for (let j = 1; j <= b.length; j++) {
      const cost = a[i - 1] === b[j - 1] ? 0 : 1;
      dp[i][j] = Math.min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost);
    }
  }
  return dp[a.length][b.length];
}

function closestMatch(candidates: string[], input: string, maxDistance: number): string | null {
  const inputLower = input.toLowerCase();
  let best: string | null = null;
  let bestDist = maxDistance + 1;
  for (const c of candidates) {
    const d = levenshtein(c.toLowerCase(), inputLower);
    if (d < bestDist) { bestDist = d; best = c; }
  }
  return bestDist <= maxDistance ? best : null;
}

function canonicalizeViewerType(input: string): string {
  const norm = input.toLowerCase().trim().replace(/\s+/g, ' ');
  const exact = VIEWER_TYPES.find((t) => t.toLowerCase() === norm);
  if (exact) return exact;
  return closestMatch(VIEWER_TYPES, input, 3) ?? input;
}

// ---------------------------------------------------------------------------
// Property application (preserved from existing grokky-api.ts)
// ---------------------------------------------------------------------------

/**
 * Apply a property bag to a viewer, validating each key against the live
 * `viewer.getProperties()` schema. Three rewrites are attempted before
 * giving up on a key:
 *   1. direct hit (`xColumnName`)
 *   2. column-name suffix (`x` → `xColumnName`, `color` → `colorColumnName`)
 *   3. plural list suffix (`categoryColumns` → `categoryColumnNames`)
 *   4. case-insensitive fallback (`xcolumnname` → `xColumnName`)
 * Unknown keys are logged via `console.warn` with a did-you-mean hint and
 * dropped silently — never throws.
 */
function applyViewerOptions(viewer: DG.Viewer, options: Record<string, any>): void {
  const knownNames = viewer.getProperties().map((p: any) => p.name);
  const knownSet = new Set<string>(knownNames);
  const lowerMap = new Map<string, string>(knownNames.map((n: string) => [n.toLowerCase(), n]));
  const accepted: Record<string, any> = {};
  const unknown: string[] = [];
  for (const [key, val] of Object.entries(options)) {
    if (knownSet.has(key)) { accepted[key] = val; continue; }
    // Column shortcut: x → xColumnName, color → colorColumnName, etc.
    if (knownSet.has(key + 'ColumnName')) { accepted[key + 'ColumnName'] = val; continue; }
    // Plural column list: categoryColumns → categoryColumnNames, etc.
    if (knownSet.has(key + 'ColumnNames')) { accepted[key + 'ColumnNames'] = val; continue; }
    // Case fix
    const caseFix = lowerMap.get(key.toLowerCase());
    if (caseFix) { accepted[caseFix] = val; continue; }
    unknown.push(key);
  }
  if (Object.keys(accepted).length > 0)
    viewer.setOptions(accepted);
  if (unknown.length > 0) {
    const hints = unknown.map((k) => {
      const guess = closestMatch(knownNames, k, 3);
      return guess ? `"${k}" (did you mean "${guess}"?)` : `"${k}"`;
    }).join(', ');
    console.warn(`addViewer "${viewer.type}": unknown ${unknown.length === 1 ? 'property' : 'properties'}: ${hints}`);
  }
}

// ---------------------------------------------------------------------------
// View resolution
// ---------------------------------------------------------------------------

/**
 * Resolve a `DG.TableView` from an in-scope `view` argument, falling back to
 * `grok.shell.tv` with a console.warn. Throws if no TableView exists anywhere.
 */
function resolveTableView(viewArg: DG.TableView | DG.ViewBase | null | undefined): DG.TableView {
  if (viewArg && viewArg instanceof DG.TableView)
    return viewArg;
  if (viewArg) {
    console.warn(
      `grokky viewer helper: \`view\` is a ${viewArg?.constructor?.name ?? 'non-TableView'}; ` +
      `falling back to grok.shell.tv. Pass a TableView explicitly for predictable scoping.`);
  }
  const tv = grok.shell.tv;
  if (!tv) throw new Error('grokky viewer helper: no active TableView (grok.shell.tv is null).');
  return tv;
}

// ---------------------------------------------------------------------------
// addViewer
// ---------------------------------------------------------------------------

/**
 * Attach a built-in viewer to the given TableView. Fuzzy-matches `type`
 * against `DG.VIEWER` (Levenshtein ≤ 3, case-insensitive), validates
 * `options` against the live `viewer.getProperties()` schema, and emits a
 * did-you-mean console warning on unknown keys.
 *
 * View scoping: pass `view` first. If `view` is a `TableView`, the viewer
 * lands there. If `view` is null or a non-TableView, falls back to
 * `grok.shell.tv` with a warning. Throws only if no active TableView exists.
 *
 * Never throws on unknown options keys — they're dropped with a warning.
 *
 * @example
 *   addViewer(view, 'Scatter plot', {x: 'height', y: 'weight', color: 'age'});
 *   addViewer(view, 'Histogram', {valueColumnName: 'activity', splitColumnName: 'class'});
 *   addViewer(view, 'Line chart', {xColumnName: 'date', yColumnNames: ['rev', 'exp']});
 */
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

// ---------------------------------------------------------------------------
// configureViewer
// ---------------------------------------------------------------------------

/**
 * Update properties on an existing viewer. Same schema validation and
 * did-you-mean logic as `addViewer`. Calls `viewer.setOptions(...)` under
 * the hood — a single batched notification.
 *
 * @example
 *   configureViewer(sp, {yAxisType: 'logarithmic', showXHistogram: true});
 */
export function configureViewer(viewer: DG.Viewer, options: Record<string, any>): void {
  applyViewerOptions(viewer, options);
}

// ---------------------------------------------------------------------------
// findViewer / findViewers
// ---------------------------------------------------------------------------

function viewerSlice(view: DG.TableView | null | undefined): DG.Viewer[] {
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grokky viewer helper: findViewer/findViewers requires a TableView.');
  // view.viewers[0] is ALWAYS the grid — skip it.
  return Array.from(view.viewers).slice(1);
}

/**
 * Find the first non-grid viewer matching `pred`. Returns `null` if no match.
 *
 * `pred` may be a function `(v: DG.Viewer) => boolean`, OR a string —
 * in which case it's treated as a type filter (`v.type === pred`,
 * canonicalized).
 *
 * Skips `view.viewers[0]` (the grid). For the grid, use `view.grid`.
 *
 * @example
 *   findViewer(view, (v) => v.type === DG.VIEWER.SCATTER_PLOT);
 *   findViewer(view, 'Scatter plot');  // equivalent shorthand
 */
export function findViewer(
  view: DG.TableView | null,
  pred: string | ((v: DG.Viewer) => boolean),
): DG.Viewer | null {
  const candidates = viewerSlice(view);
  const fn = typeof pred === 'string'
    ? (v: DG.Viewer) => v.type === canonicalizeViewerType(pred)
    : pred;
  return candidates.find(fn) ?? null;
}

/**
 * All non-grid viewers matching `pred`, or all non-grid viewers if `pred`
 * is omitted.
 *
 * @example
 *   findViewers(view);                                         // all non-grid
 *   findViewers(view, (v) => v.type === DG.VIEWER.HISTOGRAM);  // all histograms
 */
export function findViewers(
  view: DG.TableView | null,
  pred?: (v: DG.Viewer) => boolean,
): DG.Viewer[] {
  const candidates = viewerSlice(view);
  return pred ? candidates.filter(pred) : candidates;
}

// ---------------------------------------------------------------------------
// closeViewer
// ---------------------------------------------------------------------------

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

/**
 * Close one viewer, all viewers of a given type, or all viewers matching a
 * predicate. Returns the count actually closed.
 *
 * - `target: DG.Viewer` — close that one viewer. `view` may be omitted.
 * - `target: string` — close every non-grid viewer where `v.type === target`
 *   (fuzzy-canonicalized). Requires `view`.
 * - `target: (v) => boolean` — close every non-grid viewer where pred is true.
 *   Requires `view`.
 *
 * Tolerant: `viewer.close()` throws on never-attached viewers; that's
 * caught, logged, and the count reflects only successful closes.
 *
 * @example
 *   closeViewer(scatterPlot);              // by handle, returns 1
 *   closeViewer('Scatter plot', view);     // every scatter plot, returns N
 *   closeViewer((v) => v.dataFrame !== t, view);  // any viewer not bound to t
 */
export function closeViewer(
  target: DG.Viewer | string | ((v: DG.Viewer) => boolean),
  view?: DG.TableView | null,
): number {
  // Single-handle path — no view required.
  if (target instanceof DG.Viewer)
    return safeClose(target) ? 1 : 0;

  // String / predicate paths — need a view to enumerate.
  if (!view || !(view instanceof DG.TableView))
    throw new Error('grokky.closeViewer: string/predicate target requires a TableView as second arg.');

  const fn = typeof target === 'string'
    ? (v: DG.Viewer) => v.type === canonicalizeViewerType(target)
    : target;
  const candidates = viewerSlice(view).filter(fn);
  let count = 0;
  for (const v of candidates)
    if (safeClose(v)) count++;
  return count;
}

// ---------------------------------------------------------------------------
// closeAllViewers
// ---------------------------------------------------------------------------

/**
 * Close every viewer on the view. `opts.keepGrid` defaults to `true` —
 * `view.viewers[0]` (the grid) is preserved. Returns the count closed.
 *
 * Equivalent to `view.viewers.slice(1).forEach(v => v.close())` plus the
 * never-attached safety net.
 *
 * @example
 *   closeAllViewers(view);                       // reset view, keep grid (count)
 *   closeAllViewers(view, {keepGrid: false});    // nuclear: closes grid too
 */
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
  for (const v of targets)
    if (safeClose(v)) count++;
  return count;
}
