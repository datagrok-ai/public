/**
 * Playwright helpers shared by Viewers/Legend/* spec files.
 *
 * Each helper is a verbatim extraction of a block previously pasted into
 * every Legend spec — same selectors, same sleeps, same JS-API fallbacks.
 * Imported as: `import * as v from '../../helpers/viewers';`.
 *
 * Behavioral contract: a helper must produce the same observable end-state
 * as the inline block it replaces. Do not "improve" sleeps or selectors here
 * — those changes belong in a follow-up that replaces fixed setTimeout with
 * expect.poll across the suite (see review item 5).
 */

import {Page, expect} from '@playwright/test';
import {stepErrors, StepError} from '../spec-login';

// ---------------------------------------------------------------------------
// 1. openTable — canonical open prelude used by Legend + FilterPanel specs.
// ---------------------------------------------------------------------------

export interface OpenTableOptions {
  /** Demo file path. Defaults to SPGI (the standard fixture for Legend specs). */
  path?: string;
  /** When true, prime the Filter Panel (Filters viewer) after open. */
  withFilterPanel?: boolean;
  /** Semantic type hint (unused by the open path; reserved for callers). */
  semType?: 'Molecule' | 'Macromolecule';
  /** Force the OpenFile path (also auto-detected for .sdf/.nwk/.pdb). */
  sdf?: boolean;
  /** Extra settle after the Grid is visible (ms). */
  settleMs?: number;
  /** Cap on the semantic-type-detection race (ms). Defaults to 5000; pass 3000
   * to match the shorter inline prelude that the Viewers specs hand-rolled. */
  semTypeTimeoutMs?: number;
}

/**
 * Open a demo table, attach the default TableView, wait for semantic-type
 * detection + Bio/Chem render settle, and (optionally) prime the Filter Panel.
 *
 * Superset of the prelude block pasted into every Legend / FilterPanel spec.
 * CSV sources go through `readCsv` + `addTableView`; .sdf/.nwk/.pdb (or
 * `sdf: true`) route through the OpenFile function-call recorder (verbatim from
 * openers.ts openTableFromFile). Selenium class + simpleMode +
 * showFiltersIconsConstantly are set the same way as in the inline form so
 * spec behavior is unchanged.
 */
export async function openTable(page: Page, options?: OpenTableOptions): Promise<void> {
  const p = options?.path ?? 'System:DemoFiles/SPGI.csv';
  const useOpenFile = options?.sdf === true || /\.(sdf|nwk|pdb)$/i.test(p);
  const semTypeTimeoutMs = options?.semTypeTimeoutMs ?? 5000;
  await page.evaluate(async ({path, openFile, semTypeTimeoutMs}) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    let df: any;
    if (openFile) {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const fns = DG.Func.find({name: 'OpenFile'});
      if (!fns?.length) throw new Error('OpenFile function not registered');
      const fullPathForOpen = path.replace(/^([^:/]+):/, '$1.');
      await fns[0].prepare({fullPath: fullPathForOpen}).call(undefined, undefined, {processed: false});
      for (let i = 0; i < 24; i++) {
        const tv = grok.shell.tv;
        if (tv?.dataFrame && typeof tv.addViewer === 'function') {
          df = tv.dataFrame;
          break;
        }
        await new Promise((r) => setTimeout(r, 500));
      }
      if (!df) throw new Error(`OpenFile("${path}") did not produce a TableView (12s settle)`);
    } else {
      df = await (window as any).grok.dapi.files.readCsv(path);
      (window as any).grok.shell.addTableView(df);
    }
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, semTypeTimeoutMs);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        const grid = document.querySelector('[name="viewer-Grid"]');
        if (grid?.querySelector('canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, {path: p, openFile: useOpenFile, semTypeTimeoutMs});
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  if (options?.settleMs) await page.waitForTimeout(options.settleMs);
  if (options?.withFilterPanel) await openFilterPanel(page);
}

/**
 * Prime the Filter Panel (Filters viewer) on the active TableView: open the
 * filters group, wait for the first `.d4-filter`, and hover it. Extracted from
 * the withFilterPanel half of the original openTable prelude.
 */
export async function openFilterPanel(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15000});
  // Real DOM gesture on the Filter Panel — exercises the documented
  // hover-to-reveal header-icons interaction and adds an observable
  // Playwright-driven call alongside JS-API operations.
  await page.locator('[name="viewer-Filters"] .d4-filter').first().hover();
}

/**
 * Add a viewer by clicking its ribbon/toolbox icon, then wait for the viewer to
 * attach. Verbatim equivalent of the `querySelector('[name="icon-<icon>"]')`
 * click + `[name="viewer-<name>"]` waitFor block pasted into the Viewers specs.
 * Pass a non-default `timeoutMs` only when a call-site used a different wait.
 */
export async function addViewerByIcon(
  page: Page, iconName: string, viewerName: string, timeoutMs = 5000,
): Promise<void> {
  await page.evaluate((n) => {
    (document.querySelector('[name="icon-' + n + '"]') as HTMLElement).click();
  }, iconName);
  await page.locator('[name="viewer-' + viewerName + '"]').waitFor({timeout: timeoutMs});
  // The DOM node attaches a tick before the viewer is enumerable in
  // grok.shell.tv.viewers (and before its props object exists). A caller that
  // reaches for the viewer via `tv.viewers.find(...)` right after this returns
  // otherwise gets `undefined`. Poll enumerability so downstream prop reads are
  // safe regardless of the render-vs-registration race.
  // `viewerName` is the DOM-attribute form ('Bar-chart'), while viewer.type is
  // the display form ('Bar chart') — compare hyphen/space-insensitively.
  await page.waitForFunction((vn) => {
    const tv = (window as any).grok?.shell?.tv;
    if (!tv) return false;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const v = Array.from(tv.viewers).find((x: any) => norm(x.type) === norm(vn));
    return !!(v && (v as any).props);
  }, viewerName, {timeout: timeoutMs});
}

// ---------------------------------------------------------------------------
// 2. addLegendViewers — uniform viewer-attach + legend-column binding.
// ---------------------------------------------------------------------------

/**
 * Mapping from viewer type to the property that binds the legend's category
 * source. Scalar columns: scatter/histogram/line/bar/pie. Array columns:
 * Trellis + Box plot. Pulled from the per-viewer if/else blocks duplicated
 * across filtering/color-consistency/visibility-and-positioning specs.
 */
export const LEGEND_COLUMN_PROP: Record<string, {prop: string; array: boolean}> = {
  'Scatter plot':  {prop: 'colorColumnName', array: false},
  'Histogram':     {prop: 'splitColumnName', array: false},
  'Line chart':    {prop: 'splitColumnName', array: false},
  'Bar chart':     {prop: 'splitColumnName', array: false},
  'Pie chart':     {prop: 'categoryColumnName', array: false},
  'Trellis plot':  {prop: 'xColumnNames', array: true},
  'Box plot':      {prop: 'categoryColumnNames', array: true},
};

export interface ViewerSpec {
  /** Viewer type (e.g. 'Scatter plot'). */
  type: string;
  /** Override column for this viewer (e.g. 'Series' for scatterplot). Default = options.column. */
  column?: string;
  /** Override the prop (e.g. 'markersColumnName' for scatter markers). */
  prop?: string;
  /** Force value to be an array (overrides LEGEND_COLUMN_PROP.array). */
  array?: boolean;
}

/**
 * Add a list of viewers to the active TableView and bind each one's
 * legend column. Sets `legendVisibility = 'Always'` for every non-Grid
 * viewer (matching the inline pattern in filtering-spec.ts etc).
 *
 * Verbatim equivalent of the `const names = [...]; for (const n of names) tv.addViewer(n); ...`
 * blocks duplicated across 7+ specs.
 */
export async function addLegendViewers(
  page: Page,
  options: {column: string; viewers: (string | ViewerSpec)[]; settleMs?: number},
): Promise<void> {
  const specs: ViewerSpec[] = options.viewers.map((v) => typeof v === 'string' ? {type: v} : v);
  const settleMs = options.settleMs ?? 1500;
  await page.evaluate(async ({s, defaultCol, settle, map}) => {
    const tv = (window as any).grok.shell.tv;
    for (const spec of s) {
      tv.addViewer(spec.type);
      await new Promise((r) => setTimeout(r, 300));
    }
    for (const v of tv.viewers) {
      if (v.type === 'Grid') continue;
      try {
        const spec = s.find((x: ViewerSpec) => x.type === v.type);
        if (spec) {
          const col = spec.column ?? defaultCol;
          const entry = (map as any)[v.type];
          const prop = spec.prop ?? entry?.prop;
          const array = spec.array ?? entry?.array ?? false;
          if (prop) (v.props as any)[prop] = array ? [col] : col;
        }
        try { v.props.legendVisibility = 'Always'; } catch (_) {}
      } catch (_) {}
    }
    await new Promise((r) => setTimeout(r, settle));
  }, {s: specs, defaultCol: options.column, settle: settleMs, map: LEGEND_COLUMN_PROP});
}

// ---------------------------------------------------------------------------
// 2b. pickColumnViaSelector — type-and-search column selector (UI path).
// ---------------------------------------------------------------------------

export interface PickColumnOptions {
  /**
   * Suffix of the column-combobox name attribute. The selector matched is
   * `[name="div-column-combobox-<suffix>"]` — see `references/viewers.md`.
   * Examples (lowercase, multi-word uses double dash):
   *   - `'color'`              Scatter plot color
   *   - `'size'`               Scatter plot size
   *   - `'split'`              Histogram / Bar / Line split
   *   - `'split--by'`          Timelines split-by
   *   - `'category'`           Pie chart category
   *   - `'x'` / `'y'`          Generic XY axis
   *   - `'stack'`              Bar chart stack
   */
  comboboxSuffix: string;
  /** Column name to type into the selector. */
  columnName: string;
  /**
   * Optional viewer type for post-flow verification (e.g. 'Scatter plot').
   * When set together with propName, the helper reads `props[propName]` after
   * the UI flow and reports (via `usedFallback`) whether the JS-API
   * substitution had to be applied.
   */
  viewerType?: string;
  /**
   * Property name to verify the change landed (e.g. 'colorColumnName',
   * 'splitColumnName'). Pass together with viewerType.
   */
  propName?: string;
  /**
   * Opt-in JS-API substitution: when true (and viewerType/propName are set),
   * a UI flow that did not land the column falls back to assigning the prop
   * directly. Default false — a broken UI path is reported, not masked.
   */
  allowFallback?: boolean;
  /**
   * How to wait for the popup to open after the trigger mousedown:
   *   - `'sleep'` (default) — fixed 500ms wait. Matches the original
   *     scatter-plot-spec.ts:25-47 pattern.
   *   - `'backdrop'` — poll for `.d4-column-selector-backdrop` (up to 3s).
   *     Matches density-plot-spec.ts:29-38. More reliable when the popup
   *     init is slow, but the backdrop element doesn't render on every
   *     build/widget — fall back to sleep if it doesn't appear.
   *   - `'either'` — race backdrop (3s) vs 500ms sleep, whichever first.
   *     Robust default for new call sites that don't know which strategy
   *     fits.
   */
  popupWaitStrategy?: 'sleep' | 'backdrop' | 'either';
  /**
   * Optional scope for the column-combobox lookup. When provided, the
   * helper restricts the `[name="div-column-combobox-<suffix>"]` query
   * to descendants of this CSS selector — e.g. `[name="viewer-Density-plot"]`
   * for density-plot, where the same combobox suffix can appear elsewhere
   * on the page (gear panel, second viewer instance, etc).
   */
  scopeSelector?: string;
}

/**
 * Drive the column-selector widget: open the popup, type the column name,
 * press Enter. Verbatim extraction of `setColumnViaSelector` from
 * scatter-plot-spec.ts:25-47 — same mousedown-on-`.d4-column-selector-column`
 * trigger, same first-key + rest-of-name typing rhythm (avoids the timing
 * bug where the popup's async-focused search input drops the first letter),
 * same ArrowDown + Enter commit. The prop-equality JS-API fallback is
 * opt-in via `allowFallback`; the return value reports whether it fired.
 *
 * Assumes the viewer's properties Context Panel (gear) is already open, or
 * the target column-combobox is rendered somewhere on the page. Callers that
 * need to open the gear first should do so before invoking this helper.
 *
 * Reference: `.claude/skills/grok-browser/references/viewers.md` "Column
 * Selectors on Viewers" + `density-plot-run.md` rows 17-19 (UI flow
 * validated against dev).
 */
export async function pickColumnViaSelector(
  page: Page, opts: PickColumnOptions,
): Promise<{usedFallback: boolean}> {
  const selectorName = `div-column-combobox-${opts.comboboxSuffix}`;
  const scope = opts.scopeSelector ?? null;
  const strategy = opts.popupWaitStrategy ?? 'sleep';
  // Open the popup. Mousedown on .d4-column-selector-column is the proven
  // trigger — direct click or focus does NOT reliably open the popup.
  await page.evaluate(({name, sc}) => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const root: Document | Element = sc
      ? (document.querySelector(sc) as Element | null) ?? document
      : document;
    const sel = root.querySelector(`[name="${name}"]`);
    if (!sel) return;
    const colLabel = sel.querySelector('.d4-column-selector-column');
    (colLabel || sel).dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  }, {name: selectorName, sc: scope});
  if (strategy === 'sleep') {
    await page.waitForTimeout(500);
  } else if (strategy === 'backdrop') {
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 3000}).catch(() => {});
  } else { // 'either' — race
    await Promise.race([
      page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
        null, {timeout: 3000}).catch(() => {}),
      page.waitForTimeout(500),
    ]);
  }
  // Type the column name. First key separated by a 100ms wait — the popup
  // focuses its search input asynchronously via Timer.run and the first
  // letter sometimes drops if both keys land in the same tick.
  await page.keyboard.press(opts.columnName[0].toLowerCase());
  await page.waitForTimeout(100);
  if (opts.columnName.length > 1)
    await page.keyboard.type(opts.columnName.slice(1).toLowerCase());
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);

  // Post-flow verify: read the prop back and, only when the caller opted in,
  // substitute via JS API. `usedFallback` lets callers assert the UI path.
  let usedFallback = false;
  if (opts.viewerType && opts.propName) {
    const applied = await page.evaluate(({vt, prop, col}) => {
      const view = (window as any).grok.shell.tv?.viewers?.find((x: any) => x.type === vt);
      return !!view && (view.props as any)[prop] === col;
    }, {vt: opts.viewerType, prop: opts.propName, col: opts.columnName});
    if (!applied && opts.allowFallback === true) {
      await page.evaluate(({vt, prop, col}) => {
        const view = (window as any).grok.shell.tv?.viewers?.find((x: any) => x.type === vt);
        if (view) (view.props as any)[prop] = col;
      }, {vt: opts.viewerType, prop: opts.propName, col: opts.columnName});
      usedFallback = true;
    }
  }
  return {usedFallback};
}

/**
 * Open a viewer's properties Context Panel by clicking the gear icon. The
 * gear is scoped to the viewer's panel-titlebar to disambiguate across
 * multiple viewers on screen.
 *
 * Selectors from `references/viewers.md`:
 *   - viewer container: `[name="viewer-<Type>"]` (spaces preserved or dashed)
 *   - title-bar gear:   `.panel-titlebar [name="icon-font-icon-settings"]`
 */
export async function openViewerGear(page: Page, viewerType: string): Promise<void> {
  await page.evaluate((vt) => {
    const candidates = [
      `[name="viewer-${vt}"]`,
      `[name="viewer-${vt.replace(/\s+/g, '-')}"]`,
    ];
    let vEl: HTMLElement | null = null;
    for (const c of candidates) {
      vEl = document.querySelector(c) as HTMLElement | null;
      if (vEl) break;
    }
    if (!vEl) return;
    const panel = vEl.closest('.panel-base') as HTMLElement | null;
    const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
    gear?.click();
  }, viewerType);
  await page.waitForTimeout(1000);
}

// ---------------------------------------------------------------------------
// 3. readLegend — DOM read of a viewer's legend items.
// ---------------------------------------------------------------------------

export interface LegendInfo {
  /** Count of `.d4-legend-item` elements in this viewer's [name="legend"] host. */
  itemCount: number;
  /** Per-item text labels (from `.d4-legend-value`). */
  labels: string[];
  /** Whether the [name="legend"] host is present at all. */
  legendRendered: boolean;
}

/**
 * Read a viewer's legend (item count + labels) via DOM. Verbatim equivalent
 * of the `sp.root.querySelectorAll('[name="legend"] .d4-legend-item')` snippet
 * duplicated 60+ times across the Legend specs.
 */
export async function readLegend(page: Page, viewerType: string): Promise<LegendInfo> {
  return await page.evaluate((vt) => {
    const tv = (window as any).grok.shell.tv;
    const v = tv?.viewers?.find((x: any) => x.type === vt);
    if (!v) return {itemCount: 0, labels: [], legendRendered: false};
    const legendRoot = v.root.querySelector('[name="legend"]');
    const items = Array.from(v.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
    return {
      itemCount: items.length,
      labels: items.map((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim()),
      legendRendered: !!legendRoot,
    };
  }, viewerType);
}

// ---------------------------------------------------------------------------
// 4. changeLegendItemColor — picker UI flow + JS-API fallback.
// ---------------------------------------------------------------------------

export interface ChangeColorOptions {
  /** Viewer type (e.g. 'Histogram', 'Scatter plot', 'Line chart'). */
  viewerType: string;
  /** Category label (e.g. 'R_ONE'). The dialog name is `dialog-<sanitized>`. */
  category: string;
  /** Target color as rgb tuple — used to locate the swatch via inline style. */
  rgb: [number, number, number];
  /** Target color as hex (e.g. '#1f77b4') — used for the tag-verify assertion. */
  hex: string;
  /** Column carrying `.color-coding-categorical` (e.g. 'Stereo Category'). */
  column: string;
  /**
   * Additive map for fallback. The fallback applies setCategorical with the
   * full map so previously-applied colors are not reset (github-3132 invariant).
   * Default: the single {category: hex} pair.
   */
  additive?: Record<string, string>;
  /** Optional alt viewer-container selectors (e.g. 'Line chart' uses both 'viewer-Line-chart' and 'viewer-Line chart'). */
  altContainerNames?: string[];
}

/**
 * Change a legend item's color via right-click picker. Falls back to
 * `col.meta.colors.setCategorical` if the UI flow does not commit (validated
 * 2026-05-08 — both paths exist in every picker-using spec).
 *
 * Verbatim extraction of the ~50-line block duplicated in legend-github-3132,
 * legend-grok-17278, legend-grok-17438, color-consistency, scatterplot,
 * line-chart, and visibility-and-positioning specs.
 */
export async function changeLegendItemColor(page: Page, opts: ChangeColorOptions): Promise<void> {
  // Dialog name sanitization mirrors the existing specs: `R_ONE` → `dialog-R-ONE`.
  const dlgName = 'dialog-' + opts.category.replace(/[_\s]/g, '-');
  const containers = [
    `viewer-${opts.viewerType}`,
    `viewer-${opts.viewerType.replace(/\s+/g, '-')}`,
    ...(opts.altContainerNames ?? []),
  ];
  const rgbStr = `rgb(${opts.rgb[0]}, ${opts.rgb[1]}, ${opts.rgb[2]})`;
  const lowerHex = opts.hex.toLowerCase();
  let okCommitted = false;
  try {
    let item = page.locator(`[name="${containers[0]}"] [name="legend"] .d4-legend-item`)
      .filter({has: page.locator('.d4-legend-value', {hasText: new RegExp(`^${opts.category}x?$`)})}).first();
    for (let i = 1; i < containers.length && (await item.count()) === 0; i++) {
      item = page.locator(`[name="${containers[i]}"] [name="legend"] .d4-legend-item`)
        .filter({has: page.locator('.d4-legend-value', {hasText: new RegExp(`^${opts.category}x?$`)})}).first();
    }
    await item.scrollIntoViewIfNeeded();
    await item.click({button: 'right', timeout: 5000});
    await page.locator(`.d4-dialog[name="${dlgName}"]`).waitFor({timeout: 5000});
    await page.evaluate(({dn, rgb}) => {
      const dlg = document.querySelector(`.d4-dialog[name="${dn}"]`)!;
      const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
        .find((s) => s.style.backgroundColor === rgb);
      if (sw) {
        const o = {bubbles: true, cancelable: true, view: window, button: 0};
        sw.dispatchEvent(new MouseEvent('mousedown', o));
        sw.dispatchEvent(new MouseEvent('mouseup', o));
        sw.dispatchEvent(new MouseEvent('click', o));
      }
    }, {dn: dlgName, rgb: rgbStr});
    await page.waitForTimeout(200);
    await page.locator(`.d4-dialog[name="${dlgName}"] [name="button-OK"]`).click({timeout: 5000});
    await page.waitForTimeout(700);
    const tag = await page.evaluate(({c, col}) => {
      const dfCol = (window as any).grok.shell.tv.dataFrame.col(col);
      const t = JSON.parse(dfCol.tags['.color-coding-categorical'] ?? '{}');
      return String(t[c] ?? '').toLowerCase();
    }, {c: opts.category, col: opts.column});
    okCommitted = tag === lowerHex || tag.includes(lowerHex.slice(1));
  } catch (_) {
    okCommitted = false;
  }
  if (!okCommitted) {
    const additive = opts.additive ?? {[opts.category]: opts.hex};
    await page.evaluate(({col, add}) => {
      const dfCol = (window as any).grok.shell.tv.dataFrame.col(col);
      const current: Record<string, string> = {};
      try { Object.assign(current, JSON.parse(dfCol.tags['.color-coding-categorical'] ?? '{}')); } catch (_) {}
      Object.assign(current, add);
      dfCol.tags['.color-coding-type'] = 'Categorical';
      dfCol.meta.colors.setCategorical(current, {fallbackColor: '#808080'});
      for (const v of (window as any).grok.shell.tv.viewers)
        if (v.type !== 'Grid') try { v.invalidate?.(); } catch (_) {}
    }, {col: opts.column, add: additive});
    await page.waitForTimeout(800);
  }
  // Final assert: the column tag carries the requested color for the category.
  const final = await page.evaluate(({c, col}) => {
    const dfCol = (window as any).grok.shell.tv.dataFrame.col(col);
    const t = JSON.parse(dfCol.tags['.color-coding-categorical'] ?? '{}');
    return String(t[c] ?? '').toLowerCase();
  }, {c: opts.category, col: opts.column});
  expect(final).toBe(lowerHex);
}

// ---------------------------------------------------------------------------
// 5. clickCanvasFilter — Bar/Pie click-to-filter with multi-position retry.
// ---------------------------------------------------------------------------

/** Geometry positions for Bar / Pie hit-tests. Tuned to headless layouts. */
const BAR_POSITIONS = (w: number, h: number, nCats: number) => [
  {x: w * 0.5, y: h * 0.2},
  {x: w * 0.5, y: h * 0.4},
  {x: w * 0.5, y: h * 0.6},
  {x: 60 + 0.5 * (w - 80) / Math.max(nCats, 1), y: h * 0.6},
  {x: w * 0.3, y: h * 0.3},
  {x: w * 0.7, y: h * 0.3},
];
const PIE_POSITIONS = (w: number, h: number) => [
  {x: w * 0.45, y: h * 0.4},
  {x: w * 0.5, y: h * 0.35},
  {x: w * 0.6, y: h * 0.5},
  {x: w * 0.4, y: h * 0.5},
  {x: w * 0.5, y: h * 0.6},
];

export interface CanvasFilterOptions {
  /** 'Bar chart' or 'Pie chart'. */
  viewerType: 'Bar chart' | 'Pie chart';
  /** Column on which the categorical filter will narrow (e.g. 'Stereo Category'). */
  column: string;
}

export interface CanvasFilterResult {
  /** Filtered-row count after the click (or after JS-API fallback). */
  totalFiltered: number;
  /** For Bar chart: number of categories still represented (1 if narrowed). */
  survivors: number;
  /** Whether the canvas click actually narrowed the filter (false → fallback used). */
  canvasClickWorked: boolean;
}

/**
 * Click-to-filter on Bar / Pie canvas with multi-position retry. If the
 * canvas click does not narrow the filter, falls back to a Filter Panel
 * categorical filter that produces the same user-observable contract.
 *
 * Verbatim extraction of the ~70-line block in filtering-spec.ts and
 * legend-grok-17222-spec.ts.
 */
export async function clickCanvasFilter(
  page: Page, opts: CanvasFilterOptions,
): Promise<CanvasFilterResult> {
  const setup = await page.evaluate(async ({vt, col}) => {
    const tv = (window as any).grok.shell.tv;
    const df = tv.dataFrame;
    df.filter.setAll(true);
    const fg = tv.getFiltersGroup();
    for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
    const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
    if (sp) try { sp.props.filter = ''; } catch (_) {}
    await new Promise((r) => setTimeout(r, 600));
    const v = tv.viewers.find((x: any) => x.type === vt);
    v.props.onClick = 'Filter';
    await new Promise((r) => setTimeout(r, 500));
    const cv = v.root.querySelector('canvas') as HTMLCanvasElement;
    const r = cv.getBoundingClientRect();
    const nCats = df.col(col)?.categories?.length ?? 0;
    return {w: r.width, h: r.height, before: df.filter.trueCount, nCats};
  }, {vt: opts.viewerType, col: opts.column});

  const positions = opts.viewerType === 'Bar chart'
    ? BAR_POSITIONS(setup.w, setup.h, setup.nCats)
    : PIE_POSITIONS(setup.w, setup.h);

  let canvasClickWorked = false;
  let result: CanvasFilterResult = {totalFiltered: setup.before, survivors: 0, canvasClickWorked: false};
  for (const pos of positions) {
    await page.evaluate(async () => {
      (window as any).grok.shell.tv.dataFrame.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 300));
    });
    try {
      await page.locator(`[name="viewer-${opts.viewerType}"] canvas`).first()
        .click({position: {x: pos.x, y: pos.y}, timeout: 3000});
      await page.waitForTimeout(opts.viewerType === 'Bar chart' ? 800 : 900);
      const probe = await page.evaluate((col) => {
        const df = (window as any).grok.shell.tv.dataFrame;
        const c = df.col(col);
        const cats = c.categories;
        const counts: Record<string, number> = {};
        for (const x of cats) counts[x] = 0;
        for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
        const survivors = Object.entries(counts).filter(([_, n]) => n > 0).length;
        return {survivors, totalFiltered: df.filter.trueCount};
      }, opts.column);
      result = {...probe, canvasClickWorked: false};
      const ok = opts.viewerType === 'Bar chart'
        ? probe.survivors === 1 && probe.totalFiltered > 0
        : probe.totalFiltered !== setup.before && probe.totalFiltered > 0;
      if (ok) { canvasClickWorked = true; break; }
    } catch (_) { /* try next position */ }
  }
  if (!canvasClickWorked) {
    // JS-API fallback: post-condition only — narrow to one category via FP.
    const probe = await page.evaluate(async (col) => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.filter.setAll(true);
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      const DG = (window as any).DG;
      const cats = df.col(col).categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: col, selected: [cats[0]]});
      await new Promise((r) => setTimeout(r, 1500));
      const c = df.col(col);
      const counts: Record<string, number> = {};
      for (const x of c.categories) counts[x] = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
      return {
        survivors: Object.entries(counts).filter(([_, n]) => n > 0).length,
        totalFiltered: df.filter.trueCount,
      };
    }, opts.column);
    result = {...probe, canvasClickWorked: false};
  } else {
    result.canvasClickWorked = true;
  }
  return result;
}

// ---------------------------------------------------------------------------
// 6. applyCategoricalFilter — fg.updateOrAdd with CATEGORICAL filter.
// ---------------------------------------------------------------------------

/**
 * Apply a categorical filter via the Filter Panel. Verbatim extraction of
 * the `fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, ...})` block
 * duplicated across filtering / legend-grok-17222 / legend-api specs.
 */
export async function applyCategoricalFilter(
  page: Page, column: string, selected: string[], settleMs = 1500,
): Promise<{filteredCount: number}> {
  return await page.evaluate(async ({col, sel, settle}) => {
    const tv = (window as any).grok.shell.tv;
    const fg = tv.getFiltersGroup();
    const DG = (window as any).DG;
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: col, selected: sel});
    await new Promise((r) => setTimeout(r, settle));
    return {filteredCount: tv.dataFrame.filter.trueCount};
  }, {col: column, sel: selected, settle: settleMs});
}

// ---------------------------------------------------------------------------
// 7. resetFilters — clear all filters and reset df.filter to all-true.
// ---------------------------------------------------------------------------

/**
 * Reset every filter (FP + in-viewer) and set df.filter to all-true.
 * Verbatim extraction of the reset block used at the start of multiple
 * softSteps in filtering / legend-grok-17222 / scatterplot specs.
 */
export async function resetFilters(page: Page, opts: {clearScatterFilter?: boolean} = {}): Promise<void> {
  await page.evaluate(async (clearSp) => {
    const tv = (window as any).grok.shell.tv;
    if (!tv) return;
    const df = tv.dataFrame;
    df.filter.setAll(true);
    const fg = tv.getFiltersGroup();
    for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
    if (clearSp) {
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (sp) try { sp.props.filter = ''; } catch (_) {}
    }
    await new Promise((r) => setTimeout(r, 500));
  }, opts.clearScatterFilter ?? false);
}

// ---------------------------------------------------------------------------
// 8. cleanupShell — closeAll + 500ms settle, the canonical Cleanup softStep.
// ---------------------------------------------------------------------------

/**
 * Reset the shell: closeAll views + 500ms settle. Verbatim equivalent of the
 * Cleanup softStep body duplicated at the end of every Legend spec.
 *
 * Use inside `softStep('Cleanup', ...)` to preserve the soft-failure semantics.
 */
export async function cleanupShell(
  page: Page,
  opts: {clearStereoCategoryColorCoding?: boolean} = {},
): Promise<void> {
  await page.evaluate((clearColors) => {
    if (clearColors) {
      try {
        const col = (window as any).grok.shell.tv?.dataFrame.col('Stereo Category');
        if (col) {
          delete col.tags['.color-coding-categorical'];
          delete col.tags['.color-coding-type'];
        }
      } catch (_) {}
    }
    (window as any).grok.shell.closeAll();
  }, opts.clearStereoCategoryColorCoding ?? false);
  await page.waitForTimeout(500);
}

// ---------------------------------------------------------------------------
// 9. setViewerProps — drive a viewer's props through a set→wait→read-back ladder.
// ---------------------------------------------------------------------------

export interface ViewerPropStep {
  /** Props to assign this step (one or many, in insertion order). */
  set: Record<string, any>;
  /** Settle delay in ms after the assignment (default `delayMs`). */
  wait?: number;
  /** Prop name → its value is collected; array → an object keyed by those names; omitted → nothing. */
  read?: string | string[];
}

/**
 * Collapses the repeated `h.props.x = v; await sleep; r.push(h.props.x)` ladders
 * that fill the viewer specs. Finds the current table view's viewer by `type`,
 * then for each step assigns `set`, waits, and reads back `read`. Returns the
 * collected values in order. Behaviourally identical to the hand-rolled blocks:
 * same set order, same per-step delay, read-back happens after the wait.
 */
export async function setViewerProps(
  page: Page, viewerType: string, steps: ViewerPropStep[], delayMs = 300,
): Promise<any[]> {
  return page.evaluate(async ({viewerType, steps, delayMs}) => {
    const h = Array.from((window as any).grok.shell.tv.viewers)
      .find((x: any) => x.type === viewerType) as any;
    const out: any[] = [];
    for (var step of steps) {
      for (var k of Object.keys(step.set)) h.props[k] = step.set[k];
      await new Promise((res) => setTimeout(res, step.wait ?? delayMs));
      if (step.read === undefined) continue;
      if (Array.isArray(step.read)) {
        const obj: Record<string, any> = {};
        for (var rk of step.read) obj[rk] = h.props[rk];
        out.push(obj);
      } else
        out.push(h.props[step.read]);
    }
    return out;
  }, {viewerType, steps, delayMs});
}

// ---------------------------------------------------------------------------
// 10. finishSpec — trailing soft-step assertion (throws if any softStep failed).
// ---------------------------------------------------------------------------

/**
 * Verbatim equivalent of the `if (stepErrors.length > 0) throw new Error(...)`
 * block at the end of every Legend spec. Reads from the shared `stepErrors`
 * array exported by spec-login. Pass a non-default `prefix` only if a spec
 * wants a different message header.
 */
export function finishSpec(prefix = 'Step failures'): void {
  if (stepErrors.length === 0) return;
  const summary = stepErrors.map((e: StepError) => `- ${e.step}: ${e.error}`).join('\n');
  throw new Error(`${prefix}:\n${summary}`);
}
