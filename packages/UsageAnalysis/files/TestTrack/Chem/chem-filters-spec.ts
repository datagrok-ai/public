/* ---
realizes: [chem.cp.substructure-search-with-filter, filters.cp.chem-and-bio-filters]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import * as v from '../helpers/viewers';
import * as fp from '../helpers/filter-panel';

declare const grok: any;
declare const DG: any;

// Selector recon-notes (class-2: source-verified, not in the grok-browser reference):
//   .chem-mol-box > .mol-host — a rendered scaffold node inside a Scaffold Tree card.
//     Built by ui.divH([ui.div(moleculeHost, 'mol-host')], 'chem-mol-box')
//     (Chem/src/widgets/scaffold-tree.ts:641, second site :2176; styled in Chem/css/chem.css:342).
//     chem.md documents the node canvas as canvas.chem-canvas inside .chem-mol-box but does not
//     name the host div, which is the element that accepts the node click.

test.use(specTestOptions);

const CHEM_PATH = 'System:AppData/Chem/tests/spgi-100.csv';
const CHEM_FULL = 100;
const CHEM_COL = 'Structure';
const CHEM_PROBE = 'CCC(N(C)C)=O';
// Step 9's further edit. Must stay DIFFERENT from CHEM_PROBE, or "the edit reached the grid" and
// "it was withheld" look the same.
const CHEM_PROBE_2 = 'c1ccccc1';
const SCAFFOLD = 'c1ccccc1';
const SUBSTRUCTURE_FILTER = 'Chem:substructureFilter';
const SCAFFOLD_TREE_FILTER = 'Chem:scaffoldTreeFilter';
const SUBSTRUCTURE_SOURCE = 'Chem:Substructure Filter';
const SCAFFOLD_SOURCE = 'Chem:Scaffold Tree Filter';
const SIMILAR_CUTOFF = 0.6;
const AMBIENT_CONSOLE_ERROR = 'Permissions policy violation: compute-pressure';
const CONSOLE_PROBE = 'chem-filters-spec console channel probe';
// Ceiling on the clearing watch, not the length of it: the watch ends as soon as three consecutive
// reads come back clean, and only a reset that never clears spends the whole window.
const SKETCH_CLEAR_OBSERVE_MS = 60_000;
const SKETCH_CLEAR_HOLD_MS = 10_000;

const CHEM_CARD = `[name="viewer-Filters"] .d4-filter:has(.d4-filter-column-name:text-is("${CHEM_COL}"))` +
  `:has([data-source="${SUBSTRUCTURE_SOURCE}"])`;
const SCAFFOLD_CARD = `[name="viewer-Filters"] .d4-filter:has(.d4-filter-column-name:text-is("${CHEM_COL}"))` +
  `:has([data-source="${SCAFFOLD_SOURCE}"])`;
const GRID_CANVAS = 'canvas[name="canvas"]';
const HEADER_COUNTER = '[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator';
const SKETCHER_DIALOG = '.d4-dialog:has(input[placeholder*="SMILES"])';
const ADD_SCAFFOLD_DIALOG = '.d4-dialog:has(input[placeholder*="SMILES"]):has([name="button-Add"])';
const POPUP_HOST = '.d4-popup-host';
const POPUP_CHEM_CARD = `${POPUP_HOST} [data-source="${SUBSTRUCTURE_SOURCE}"]`;

interface ChemReadiness { ok: boolean; detail: string; }

async function chemReadiness(page: Page, column: string): Promise<ChemReadiness> {
  return page.evaluate(async (col) => {
    let detail = '(no probe completed)';
    for (let i = 0; i < 40; i++) {
      const semType = grok.shell.tv?.dataFrame?.col(col)?.semType ?? null;
      const factory = (DG.Func.find({name: 'substructureFilter'}) || [])
        .find((f: any) => f.package?.name === 'Chem');
      const action = (DG.Func.find({meta: {action: 'Use as filter'}}) || [])
        .find((f: any) => f.package?.name === 'Chem');
      let applicable = false;
      let applyErr = '';
      if (factory) {
        try { applicable = !!(await factory.apply()); }
        catch (e: any) { applyErr = e?.message ?? String(e); }
      }
      detail = `semType(${col})=${semType}; Chem:substructureFilter registered=${!!factory}; ` +
        `factory applicable=${applicable}${applyErr ? ` (apply threw: ${applyErr})` : ''}; ` +
        `Chem meta.action "Use as filter" registered=${!!action}`;
      if (semType === 'Molecule' && factory && applicable && action) return {ok: true, detail};
      await new Promise((r) => setTimeout(r, 500));
    }
    return {ok: false, detail};
  }, column);
}

async function openChemFilterPanel(page: Page): Promise<void> {
  try {
    await v.openFilterPanel(page);
  }
  catch (e) {
    const seen = await page.evaluate(() => ({
      panelRoot: !!document.querySelector('[name="viewer-Filters"]'),
      cards: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
    }));
    const molecular = await moleculeColumns(page);
    throw new Error('Chem filters are registered (the readiness gate passed) but the Filter Panel produced no ' +
      `filter card within the shared barrier — panel root present: ${seen.panelRoot}; .d4-filter cards: ` +
      `${seen.cards}; columns reporting semType "Molecule": ${molecular.length ? molecular.join(', ') : '(none)'}. ` +
      `Underlying barrier failure: ${(e as Error).message}`);
  }
}

async function trueCount(page: Page): Promise<number> {
  const c = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
  expect(typeof c, 'trueCount must be a number, not null/undefined').toBe('number');
  return c;
}

async function tableViewCount(page: Page): Promise<number> {
  return page.evaluate(() => Array.from(grok.shell.tableViews).length);
}

async function cardCount(page: Page): Promise<number> {
  return page.evaluate(() => document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length);
}

interface HeaderCounter { present: boolean; text: string; }

async function settledHeaderCounter(page: Page): Promise<HeaderCounter> {
  const read = () => page.evaluate((sel) => {
    const el = document.querySelector(sel);
    return {present: !!el, text: (el?.textContent ?? '').trim()};
  }, HEADER_COUNTER);
  try {
    await page.waitForFunction((sel) => {
      const w = window as any;
      const el = document.querySelector(sel);
      const text = (el?.textContent ?? '').trim();
      if (!el || text.length === 0) { w.__hdrCounter = null; return false; }
      const st = w.__hdrCounter;
      if (st && st.text === text) st.n++;
      else w.__hdrCounter = {text, n: 1};
      return w.__hdrCounter.n >= 3;
    }, HEADER_COUNTER, {timeout: 30_000, polling: 250});
  }
  catch (_) {
    const seen = await read();
    throw new Error('the Filter Panel header active-filter counter never held one value for 750 ms within 30 s, so ' +
      'every reading taken from it below would be a value caught in passing rather than the counter\'s settled ' +
      `state (last read: present=${seen.present}, text="${seen.text}")`);
  }
  return read();
}

async function waitForCardCount(page: Page, expected: number, timeoutMs = 20_000): Promise<void> {
  await page.waitForFunction((n) =>
    document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length === n,
  expected, {timeout: timeoutMs, polling: 200});
}

async function waitForFilterIdle(page: Page, timeoutMs = 90_000): Promise<number> {
  await page.waitForFunction(() => {
    const w = window as any;
    const busy = [...document.querySelectorAll('.chem-filter .grok-loader')]
      .some((l) => (l as HTMLElement).style.display === 'initial');
    if (busy) { w.__filterIdle = null; return false; }
    const tc = w.grok.shell.tv.dataFrame.filter.trueCount;
    if (typeof tc !== 'number') { w.__filterIdle = null; return false; }
    const st = w.__filterIdle;
    if (st && st.tc === tc) st.n++;
    else w.__filterIdle = {tc, n: 1};
    return w.__filterIdle.n >= 3;
  }, null, {timeout: timeoutMs, polling: 250});
  return trueCount(page);
}

async function sustainedCount(page: Page, label: string, attempts = 5): Promise<number> {
  let first = 0;
  let second = -1;
  for (let i = 0; i < attempts && first !== second; i++) {
    first = await waitForFilterIdle(page);
    await page.waitForTimeout(2500);
    second = await trueCount(page);
  }
  expect(second, `${label}: the row count must HOLD across a re-read, not merely be caught in passing ` +
    `(settled on ${first}, re-read ${second} after ${attempts} attempts)`).toBe(first);
  return first;
}

async function installFilterPassCounter(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    w.__filterPassSub?.unsubscribe();
    w.__filterPasses = 0;
    w.__filterPassSub = grok.shell.tv.dataFrame.onRowsFiltered
      .subscribe(() => { w.__filterPasses++; });
  });
}

async function filterPasses(page: Page): Promise<number> {
  const n = await page.evaluate(() => (window as any).__filterPasses);
  expect(typeof n, 'the filtering-pass counter is not installed on the table under test').toBe('number');
  return n;
}

async function countForSearchType(page: Page, type: string, passesBefore: number | null): Promise<number> {
  try {
    await page.waitForFunction(({col, filterType, mode, n}) => {
      const w = window as any;
      const states = grok.shell.tv.getFiltersGroup().getStates(col, filterType) || [];
      return states.some((s: any) => s?.searchType === mode) && (n === null || (w.__filterPasses ?? 0) > n);
    }, {col: CHEM_COL, filterType: SUBSTRUCTURE_FILTER, mode: type, n: passesBefore},
    {timeout: 90_000, polling: 250});
  }
  catch (_) {
    const seen = await page.evaluate(({col, filterType}) => ({
      searchTypes: (grok.shell.tv.getFiltersGroup().getStates(col, filterType) || [])
        .map((s: any) => String(s?.searchType)),
      passes: (window as any).__filterPasses,
    }), {col: CHEM_COL, filterType: SUBSTRUCTURE_FILTER});
    throw new Error(`switching the ${CHEM_COL} card to "${type}" produced no filtering pass on the table within ` +
      `90 s (the card's own states read ${seen.searchTypes.join(', ') || '(none)'}; the frame has run ` +
      `${passesBefore === null ? 'n/a' : seen.passes - passesBefore} filtering passes since the switch). The ` +
      'recompute for this mode is what moves the row count; without an observed pass the count read here is the ' +
      'PREVIOUS mode\'s count carried forward, and every subset and partition claim below would be made about it');
  }
  return sustainedCount(page, `search type "${type}"`);
}

async function menuLeaf(page: Page, group: string, leaf: string): Promise<void> {
  await page.evaluate(async ({g, l}) => {
    const norm = (s: string | null | undefined) => (s ?? '').trim().toLowerCase();
    const labels = () => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])];
    };
    const visible = () => labels().map((i) => (i.textContent ?? '').trim()).join(' | ');
    const find = (text: string) =>
      labels().find((i) => norm(i.textContent) === norm(text))?.closest('.d4-menu-item') ?? null;
    const waitFor = async (text: string) => {
      const deadline = Date.now() + 5000;
      let item = find(text);
      while (!item && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 100));
        item = find(text);
      }
      return item;
    };
    const groupItem = await waitFor(g);
    if (!groupItem)
      throw new Error(`menu: group "${g}" not found in the last popup; visible items: ${visible()}`);
    const b = groupItem.getBoundingClientRect();
    for (const type of ['mouseover', 'mousemove'])
      groupItem.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
    const leafItem = await waitFor(l);
    if (!leafItem)
      throw new Error(`menu: leaf "${l}" not found under "${g}" in the last popup; visible items: ${visible()}`);
    (leafItem as HTMLElement).click();
  }, {g: group, l: leaf});
}

async function installCardFinder(page: Page): Promise<void> {
  await page.evaluate(() => {
    (window as any).__filterCard = (col: string, source: string, viewIdx = 0): HTMLElement | null => {
      const views: any[] = [];
      for (const t of grok.shell.tableViews) views.push(t);
      const view = views[viewIdx];
      if (!view) return null;
      const fv = [...(view.viewers ?? [])].find((x: any) => x.type === 'Filters');
      const root: ParentNode | null = fv?.root ?? view.root ?? null;
      if (!root) return null;
      const card = [...root.querySelectorAll('.d4-filter')].find((c) =>
        (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col &&
        !!c.querySelector(`[data-source="${source}"]`));
      return (card as HTMLElement) ?? null;
    };
  });
}

async function pickAllColumnsAndConfirm(page: Page): Promise<number> {
  const dialog = page.locator('.d4-dialog[name="dialog-Select-columns..."]');
  await dialog.waitFor({timeout: 20_000});
  const readChecked = async () => page.evaluate(() => {
    const dlg = document.querySelector('.d4-dialog[name="dialog-Select-columns..."]');
    const label = [...(dlg?.querySelectorAll('label') ?? [])].find((l) => /checked/.test(l.textContent ?? ''));
    const n = parseInt((label?.textContent ?? '').trim(), 10);
    return Number.isFinite(n) ? n : 0;
  });
  let checked = 0;
  for (let i = 0; i < 15 && checked === 0; i++) {
    await dialog.locator('[name="label-All"]').click();
    await page.waitForTimeout(500);
    checked = await readChecked();
  }
  await dialog.locator('[name="button-OK"]').click();
  await dialog.waitFor({state: 'detached', timeout: 20_000});
  return checked;
}

async function moleculeColumns(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const out: string[] = [];
    for (let i = 0; i < df.columns.length; i++) {
      const c = df.columns.byIndex(i);
      if (c.semType === 'Molecule') out.push(c.name);
    }
    return out;
  });
}

async function gridCellPoint(page: Page, column: string, rowIdx: number):
  Promise<{x: number; y: number; tableRow: number} | null> {
  return page.evaluate(async ({col, row}) => {
    const grid = grok.shell.tv.grid;
    grid.scrollToCell(col, row);
    await new Promise((r) => setTimeout(r, 500));
    const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')].find((g) => !g.closest('.d4-filter'));
    const overlay = mainGrid?.querySelector('[name="overlay"]') as HTMLElement | null;
    const gc = grid.columns.byName(col);
    if (!overlay || !gc) return null;
    const rect = overlay.getBoundingClientRect();
    const xLocal = Math.round((gc.left + gc.right) / 2);
    const band: number[] = [];
    // Grid row != table row whenever a filter or sort is armed; the caller asserts identity against this.
    let tableRow: number | null = null;
    for (let y = Math.round(grid.colHeaderHeight) + 2; y < rect.height; y += 4) {
      const hit = grid.hitTest(xLocal, y);
      if (hit && hit.gridRow === row && hit.gridColumn && hit.gridColumn.name === col && hit.isTableCell) {
        band.push(y);
        if (tableRow === null) tableRow = hit.tableRowIndex;
      }
      else if (band.length > 0)
        break;
    }
    if (band.length === 0 || tableRow === null) return null;
    const yLocal = band[Math.floor(band.length / 2)];
    return {x: Math.round(rect.left + xLocal), y: Math.round(rect.top + yLocal), tableRow};
  }, {col: column, row: rowIdx});
}

interface HeaderPoint {
  x: number;
  y: number;
  // Carried so a miss reports the arithmetic instead of just "nothing happened".
  scrollMin: number;
  virtualMid: number;
  left: number;
  right: number;
}

// `gc.left` / `gc.right` are VIRTUAL grid coordinates, not screen ones: they are measured from
// the start of the whole column run, so they only coincide with screen x while the grid is
// scrolled fully left. This probe calls scrollToCell first, which is exactly what makes them
// diverge — after the scroll `horzScroll.min` is the width scrolled past, and a point computed
// without subtracting it lands that far to the right, possibly outside the window. The press
// then hits nothing, the drag never starts, and the failure surfaces much later as "the drop
// zone never appeared".
async function gridHeaderPoint(page: Page, column: string): Promise<HeaderPoint | null> {
  return page.evaluate(async (col) => {
    const grid = grok.shell.tv.grid;
    grid.scrollToCell(col, 0);
    await new Promise((r) => setTimeout(r, 500));
    const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')].find((g) => !g.closest('.d4-filter'));
    const overlay = mainGrid?.querySelector('[name="overlay"]') as HTMLElement | null;
    const gc = grid.columns.byName(col);
    if (!overlay || !gc) return null;
    const rect = overlay.getBoundingClientRect();
    const scrollMin = grid.horzScroll?.min ?? 0;
    const virtualMid = (gc.left + gc.right) / 2;
    return {
      x: Math.round(rect.left + virtualMid - scrollMin),
      y: Math.round(rect.top + grid.colHeaderHeight / 2),
      scrollMin,
      virtualMid,
      left: rect.left,
      right: rect.right,
    };
  }, column);
}

// Fails at the point of aiming rather than at the consequence. Without it a press outside the
// grid is indistinguishable from a product that ignored the gesture.
function assertHeaderPointOnScreen(pt: HeaderPoint | null, column: string): asserts pt is HeaderPoint {
  expect(pt, `the ${column} grid column header could not be located, so any failure below would be ` +
    'about the probe rather than the product').not.toBeNull();
  expect(pt!.x >= pt!.left && pt!.x <= pt!.right,
    `the computed ${column} header point x=${pt!.x} falls outside the grid overlay ` +
    `[${Math.round(pt!.left)}, ${Math.round(pt!.right)}] — virtual mid ${Math.round(pt!.virtualMid)}, ` +
    `horzScroll.min ${Math.round(pt!.scrollMin)}. Pressing there hits nothing.`).toBe(true);
}

interface ColumnPopup { withinColumn: boolean; title: string; }

async function openColumnOptionsPopup(page: Page, column: string): Promise<ColumnPopup> {
  const pt = await gridHeaderPoint(page, column);
  assertHeaderPointOnScreen(pt, column);
  await page.mouse.move(pt!.x - 60, pt!.y, {steps: 4});
  await page.mouse.move(pt!.x, pt!.y, {steps: 8});
  await page.waitForFunction(() => [...document.querySelectorAll('[column_name]')]
    .some((e) => (e as HTMLElement).getBoundingClientRect().width > 0), null, {timeout: 20_000, polling: 200});
  const withinColumn = await page.evaluate((col) => {
    const icon = [...document.querySelectorAll('[column_name]')]
      .find((e) => (e as HTMLElement).getBoundingClientRect().width > 0) as HTMLElement | null;
    if (!icon) return false;
    const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')].find((g) => !g.closest('.d4-filter'));
    const overlay = mainGrid?.querySelector('[name="overlay"]') as HTMLElement | null;
    const gc = grok.shell.tv.grid.columns.byName(col);
    if (!overlay || !gc) return false;
    const r = overlay.getBoundingClientRect();
    const ir = icon.getBoundingClientRect();
    const cx = ir.left + ir.width / 2;
    const inside = cx >= r.left + gc.left && cx <= r.left + gc.right;
    icon.click();
    return inside;
  }, column);
  await page.locator(POPUP_CHEM_CARD).first().waitFor({timeout: 45_000});
  const title = await page.evaluate(() => {
    const host = document.querySelector('.d4-popup-host');
    return ([...(host?.querySelectorAll('.ui-label') ?? [])][0]?.textContent ?? '').trim();
  });
  return {withinColumn, title};
}

async function chemFilterRoots(page: Page): Promise<{total: number; inPanel: number; outside: number}> {
  return page.evaluate(() => {
    const all = [...document.querySelectorAll('.chem-filter')];
    const inPanel = all.filter((e) => !!e.closest('[name="viewer-Filters"]')).length;
    return {total: all.length, inPanel, outside: all.length - inPanel};
  });
}

interface CardBody {
  column: string;
  canvas: boolean;
  link: boolean;
  w: number;
  h: number;
  ink: number;
  opaque: number;
  molBlock: string | null;
}

async function cardBodies(page: Page): Promise<CardBody[]> {
  return page.evaluate((type) => {
    const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .filter((c) => !!c.querySelector('.chem-filter'));
    const fg = grok.shell.tv.getFiltersGroup();
    return cards.map((c) => {
      const column = (c.querySelector('.d4-filter-column-name')?.textContent ?? '?').trim();
      const cv = c.querySelector('.chem-external-sketcher-canvas') as HTMLCanvasElement | null;
      let ink = -1;
      let opaque = -1;
      if (cv && cv.width > 0 && cv.height > 0) {
        const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        ink = 0;
        opaque = 0;
        for (let i = 0; i < d.length; i += 4) {
          if (d[i + 3] <= 8) continue;
          opaque++;
          if (d[i] < 245 || d[i + 1] < 245 || d[i + 2] < 245) ink++;
        }
      }
      const st = (fg.getStates(column, type) || [])[0];
      return {
        column,
        canvas: !!cv,
        link: !!c.querySelector('.sketch-link'),
        w: cv?.width ?? 0,
        h: cv?.height ?? 0,
        ink,
        opaque,
        molBlock: typeof st?.molBlock === 'string' ? st.molBlock : null,
      };
    });
  }, SUBSTRUCTURE_FILTER);
}

function bodyTable(label: string, bodies: CardBody[]): string {
  const rows = bodies.map((b) => `  ${b.column}: canvas=${b.canvas} link=${b.link} ` +
    `size=${b.w}x${b.h} inkPx=${b.ink} opaquePx=${b.opaque} molAtoms=${molAtomCount(b.molBlock)} ` +
    `molLen=${(b.molBlock ?? '').length}`);
  return `${label}\n${rows.join('\n')}`;
}

interface SketchAreas {
  total: number;
  drawn: number;
  unreadable: number;
  placeholder: number;
  blankCanvas: number;
  drawnColumns: string[];
  unreadableColumns: string[];
  detail: string;
}

// `cardBodies` leaves ink at -1 for a canvas it could not measure (0x0, or no 2d context), and
// `drawn` counts ink > 0 — so an unreadable body would otherwise read as a cleared one. It is
// counted separately and the caller fails on it: the reset resizes this canvas (204x102 -> 100x50),
// so a body caught mid-teardown is a real state, not a hypothetical.
function classifyBodies(bodies: CardBody[]): SketchAreas {
  const drawn = bodies.filter((b) => b.ink > 0);
  const unreadable = bodies.filter((b) => b.canvas && b.ink < 0);
  return {
    total: bodies.length,
    drawn: drawn.length,
    unreadable: unreadable.length,
    placeholder: bodies.filter((b) => b.link).length,
    blankCanvas: bodies.filter((b) => b.canvas && b.ink === 0 && !b.link).length,
    drawnColumns: drawn.map((b) => `${b.column} (${b.ink} px)`),
    unreadableColumns: unreadable.map((b) => `${b.column} (canvas ${b.w}x${b.h}, ink unreadable)`),
    detail: bodies.map((b) => `${b.column}=${b.link ? 'placeholder'
      : b.canvas ? `canvas ${b.w}x${b.h} ink ${b.ink}px` : 'neither'}`).join(', '),
  };
}

async function sketchAreaBodies(page: Page): Promise<SketchAreas> {
  return classifyBodies(await cardBodies(page));
}

interface ClearingObservation {
  clearedAfterMs: number | null;
  observedMs: number;
  samples: number;
  unreadableSamples: number;
  lastDrawn: number;
  lastDrawnColumns: string[];
  lastUnreadableColumns: string[];
}

// A sample only counts toward the cleared run when every body was READABLE as well as unpainted:
// an unreadable canvas is the state a mid-teardown body is in, and counting it as clear is how the
// window would end 1.5 s after a reset that had not finished repainting.
async function observeSketchClearing(page: Page, windowMs: number): Promise<ClearingObservation> {
  const started = Date.now();
  let samples = 0;
  let unreadableSamples = 0;
  let zeroes = 0;
  let clearedAt: number | null = null;
  let last = await sketchAreaBodies(page);
  while (Date.now() - started < windowMs) {
    samples++;
    last = await sketchAreaBodies(page);
    if (last.unreadable > 0) unreadableSamples++;
    if (last.drawn === 0 && last.unreadable === 0) {
      if (clearedAt === null) clearedAt = Date.now() - started;
      if (++zeroes >= 3)
        return {clearedAfterMs: clearedAt, observedMs: Date.now() - started, samples, unreadableSamples,
          lastDrawn: 0, lastDrawnColumns: [], lastUnreadableColumns: []};
    }
    else { zeroes = 0; clearedAt = null; }
    await page.waitForTimeout(500);
  }
  return {clearedAfterMs: null, observedMs: Date.now() - started, samples, unreadableSamples,
    lastDrawn: last.drawn, lastDrawnColumns: last.drawnColumns,
    lastUnreadableColumns: last.unreadableColumns};
}

// Each card goes down through a real click on its own X icon, one card at a time. The icon lives in
// the same hover-revealed `.d4-filter-controls` block as the enable checkbox
// (grok-browser/references/viewers/filters.md:59), so the card is hovered first; without that the
// icon is 0x0 and only an in-page click reaches it, which proves the handler rather than the
// gesture. Step 1's claim is that the cards go down through the X, and this helper is also the
// entry state for Steps 6, 10 and 11.
async function removeSubstructureCards(page: Page, columns: string[]): Promise<number> {
  const target = '[name="viewer-Filters"] .d4-filter:has(.chem-filter)';
  let removed = 0;
  for (let i = 0; i < 50; i++) {
    const before = await page.locator(target).count();
    if (before === 0) break;
    const card = page.locator(target).first();
    const column = (await card.locator('.d4-filter-column-name').first().textContent() ?? '?').trim();
    await card.hover();
    const icon = card.locator('[name="icon-times"]').first();
    expect(await icon.count(),
      `the ${column} substructure card carries no [name="icon-times"] remove icon, so its removal cannot be ` +
      'actuated through the card\'s own control').toBeGreaterThan(0);
    await icon.waitFor({state: 'visible', timeout: 15_000});
    await icon.click();
    await expect.poll(() => page.locator(target).count(),
      {timeout: 15_000, intervals: [200, 300, 500], message:
        `clicking the ${column} substructure card's own X icon must take exactly that card out of the panel — ` +
        `${before} substructure cards were in it before the click`}).toBe(before - 1);
    removed++;
  }
  await page.waitForFunction(({cols, type}) => {
    const fg = grok.shell.tv.getFiltersGroup();
    return cols.every((c: string) => (fg.getStates(c, type) || []).length === 0);
  }, {cols: columns, type: SUBSTRUCTURE_FILTER}, {timeout: 20_000, polling: 200});
  await page.waitForTimeout(2000);
  const zeroStates = await page.evaluate(({cols, type}) => {
    const fg = grok.shell.tv.getFiltersGroup();
    return cols.map((c: string) => (fg.getStates(c, type) || []).length);
  }, {cols: columns, type: SUBSTRUCTURE_FILTER});
  expect(zeroStates, 'no molecular column may still carry a substructure filter state after the cards are taken down')
    .toEqual(columns.map(() => 0));
  expect((await chemFilterRoots(page)).total,
    'no substructure card body may be left anywhere in the document after the cards are taken down').toBe(0);
  return removed;
}

async function removeAllViaHamburger(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(() => cardCount(page),
    {timeout: 20_000, intervals: [300, 600, 1200], message:
      'the panel menu "Remove All" leaf was driven but the panel still carries cards — the reopen below rebuilds ' +
      'one card per eligible column only from an EMPTY panel; a non-empty one is restored from its saved card set ' +
      'instead, and the configuration the step measures would be inherited rather than guaranteed'}).toBe(0);
}

async function closeFilterPanel(page: Page): Promise<void> {
  await page.locator('[name="viewer-Filters"]').first().hover();
  let clicked = false;
  for (const icon of ['icon-times', 'Close']) {
    try { await v.clickViewerTitlebarIcon(page, 'Filters', icon); clicked = true; break; }
    catch (_) { clicked = false; }
  }
  expect(clicked, 'the Filter Panel exposes no title-bar close control, so the panel cannot be closed and the ' +
    'auto-population that the reopen triggers never runs').toBe(true);
  await expect.poll(() => page.locator('[name="viewer-Filters"]').count(),
    {timeout: 15_000, intervals: [300, 600, 1200], message:
      'the Filter Panel title-bar close control was clicked but a panel root is still in the document — a panel ' +
      'that never closed is not reopened either, and the cards below would be the ones already there'}).toBe(0);
}

async function reopenFilterPanel(page: Page): Promise<void> {
  await page.locator('.d4-ribbon-panel [name="icon-filter"]').first().click();
  await page.locator('[name="viewer-Filters"]').first().waitFor({timeout: 20_000});
  await expect.poll(() => cardCount(page),
    {timeout: 30_000, intervals: [400, 800, 1500], message:
      'the ribbon funnel icon was clicked and a panel root appeared, but it carries no card at all — the reopen ' +
      'is the auto-population path this precondition is built on'}).toBeGreaterThan(0);
}

// The enable checkbox is driven with a real click — it carries the central claim of Steps 5, 10, 12
// and 13. It needs the hover: the card's `.d4-filter-controls` block is display:none until the pointer
// is over the card, so the input measures 0x0 and Playwright refuses to click it (measured on dev
// 2026-08-22 — hovered, the click toggles the card and the row count follows in both directions;
// unhovered the click times out and nothing moves).
async function setCardEnabled(page: Page, on: boolean): Promise<boolean> {
  const card = page.locator(CHEM_CARD).first();
  if (await card.count() === 0) return false;
  const read = () => card.evaluate((c) => {
    const cb = c.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null;
    return {checked: cb ? cb.checked : null, disabled: c.classList.contains('d4-filter-disabled')};
  });
  await card.hover();
  const cb = card.locator('input[type="checkbox"].ui-input-editor').first();
  await cb.waitFor({state: 'visible', timeout: 15_000});
  if ((await read()).checked !== on) await cb.click();
  for (let i = 0; i < 30; i++) {
    const s = await read();
    if (s.checked === on && s.disabled === !on) break;
    await page.waitForTimeout(100);
  }
  const st = await read();
  return st.checked === on && st.disabled === !on;
}

interface SubstructureState {
  molBlock: string | null;
  searchType: string | null;
  simCutOff: number | null;
  fp: string | null;
}

async function substructureStates(page: Page, viewIdx = 0): Promise<SubstructureState[]> {
  return page.evaluate(({vi, col, type}) => {
    const views: any[] = [];
    for (const t of grok.shell.tableViews) views.push(t);
    const view = views[vi];
    if (!view) return [];
    const states = view.getFiltersGroup().getStates(col, type) || [];
    return states.map((s: any) => ({
      molBlock: typeof s?.molBlock === 'string' ? s.molBlock : null,
      searchType: typeof s?.searchType === 'string' ? s.searchType : null,
      simCutOff: typeof s?.simCutOff === 'number' ? s.simCutOff : null,
      fp: typeof s?.fp === 'string' ? s.fp : null,
    }));
  }, {vi: viewIdx, col: CHEM_COL, type: SUBSTRUCTURE_FILTER});
}

function molAtomCount(molBlock: string | null): number {
  if (typeof molBlock !== 'string' || molBlock.trim().length === 0) return 0;
  const lines = molBlock.split('\n');
  const v3 = lines.find((l) => l.includes('V30 COUNTS'));
  if (v3) {
    const n = parseInt(v3.trim().split(/\s+/)[3] ?? '0', 10);
    return Number.isFinite(n) ? n : 0;
  }
  if (lines.length < 4) return 0;
  const n = parseInt(lines[3].slice(0, 3).trim(), 10);
  return Number.isFinite(n) ? n : 0;
}

function hasStructure(molBlock: string | null): boolean {
  if (typeof molBlock !== 'string' || molBlock.trim().length === 0) return false;
  const looksLikeMolBlock = molBlock.includes('V2000') || molBlock.includes('V3000')
    || molBlock.split('\n').length >= 4;
  return looksLikeMolBlock ? molAtomCount(molBlock) > 0 : true;
}

async function revealSearchOptions(page: Page): Promise<void> {
  const revealed = await page.evaluate(({col, source}) => {
    const card = (window as any).__filterCard(col, source);
    if (!card) return false;
    if (card.querySelector('.chem-filter-search-type select')) return true;
    const gear = card.querySelector('.chem-search-options-icon') as HTMLElement | null;
    if (!gear) return false;
    gear.click();
    return true;
  }, {col: CHEM_COL, source: SUBSTRUCTURE_SOURCE});
  expect(revealed, `the ${CHEM_COL} substructure card must be present and expose its search options`).toBe(true);
  await page.locator(`${CHEM_CARD} .chem-filter-search-type select`).first()
    .waitFor({state: 'attached', timeout: 15_000});
}

async function setSearchType(page: Page, type: string): Promise<boolean> {
  const select = page.locator(`${CHEM_CARD} .chem-filter-search-type select`).first();
  if (await select.count() === 0) return false;
  const offered = await select.evaluate((s) =>
    [...(s as HTMLSelectElement).options].map((o) => o.value));
  if (!offered.includes(type)) return false;
  await select.selectOption(type);
  return true;
}

async function readSearchType(page: Page, viewIdx = 0): Promise<string | null> {
  return page.evaluate(({vi, col, source}) => {
    const card = (window as any).__filterCard(col, source, vi);
    const sel = card?.querySelector('.chem-filter-search-type select') as HTMLSelectElement | null;
    return sel ? sel.value : null;
  }, {vi: viewIdx, col: CHEM_COL, source: SUBSTRUCTURE_SOURCE});
}

async function setSimilarityCutoff(page: Page, value: number): Promise<boolean> {
  const editor = page.locator(`${CHEM_CARD} .chem-filter-sim-cutoff-editor`).first();
  if (await editor.count() === 0) return false;
  await editor.fill(String(value));
  await editor.press('Enter');
  return true;
}

async function openSketcherDialog(page: Page): Promise<void> {
  const card = page.locator(CHEM_CARD).first();
  const link = card.locator('.sketch-link').first();
  const mini = card.locator('.chem-external-sketcher-canvas').first();
  if (await link.count() > 0 && await link.isVisible()) await link.click();
  else await mini.evaluate((el) => (el.parentElement as HTMLElement).click());
  await page.locator(SKETCHER_DIALOG).first().waitFor({timeout: 15_000});
}

async function enterSmiles(page: Page, smiles: string, dialog = SKETCHER_DIALOG): Promise<void> {
  const input = page.locator(`${dialog} input[placeholder*="SMILES"]`).first();
  await input.fill(smiles);
  await input.press('Enter');
  await page.waitForFunction((probe) => {
    const dlg = [...document.querySelectorAll('.d4-dialog')]
      .find((d) => !!d.querySelector('input[placeholder*="SMILES"]'));
    if (!dlg) return true;
    const inp = dlg.querySelector('input[placeholder*="SMILES"]') as HTMLInputElement | null;
    return !!inp && inp.value.trim() === probe;
  }, smiles, {timeout: 15_000, polling: 200});
}

async function commitSmiles(page: Page, smiles: string): Promise<void> {
  await enterSmiles(page, smiles);
  const ok = page.locator(`${SKETCHER_DIALOG} [name="button-OK"]`).first();
  if (await ok.count() > 0) await ok.click();
}

interface AlignHighlightState { set: boolean; align: boolean; highlight: boolean; }

// The three sketcher checkboxes are driven with real clicks — Step 9's two claims turn on them, and
// "Filter as you draw" is the switch the whole propagation claim measures the consequence of. Unlike
// the card's enable checkbox they need no hover: they are visible as soon as the dialog is up
// (measured on dev 2026-08-22). Returns the read-back state, or null when the control is absent or
// is not a checkbox — the callers assert on that rather than on the click having been attempted.
async function setDialogCheckbox(page: Page, name: string, on: boolean): Promise<boolean | null> {
  const cb = page.locator(`${SKETCHER_DIALOG} .chem-sketcher-filter-options [name="${name}"]`).first();
  if (await cb.count() === 0) return null;
  if (await cb.evaluate((e) => (e as HTMLInputElement).type) !== 'checkbox') return null;
  await cb.waitFor({state: 'visible', timeout: 15_000});
  if (await cb.isChecked() !== on) {
    await cb.click();
    for (let i = 0; i < 20 && await cb.isChecked() !== on; i++)
      await page.waitForTimeout(100);
  }
  return cb.isChecked();
}

async function setAlignHighlight(page: Page, on: boolean): Promise<AlignHighlightState> {
  const align = await setDialogCheckbox(page, 'input-Align', on);
  const highlight = await setDialogCheckbox(page, 'input-Highlight', on);
  if (align === null || highlight === null) return {set: false, align: false, highlight: false};
  return {set: true, align, highlight};
}

async function expectAlignHighlight(page: Page, on: boolean, why: string): Promise<void> {
  const st = await setAlignHighlight(page, on);
  expect(st.set, `${why}: the sketcher must expose its Align and Highlight checkboxes`).toBe(true);
  expect({align: st.align, highlight: st.highlight},
    `${why}: both Align and Highlight must READ ${on ? 'CHECKED' : 'CLEARED'} after being set so — ` +
    'the state they start from is measured, never assumed').toEqual({align: on, highlight: on});
}

async function setFilterAsYouDraw(page: Page, on: boolean): Promise<{set: boolean; value: boolean}> {
  const value = await setDialogCheckbox(page, 'input-Filter-as-you-draw', on);
  return value === null ? {set: false, value: false} : {set: true, value};
}

async function panelCardOrder(page: Page): Promise<Array<{column: string; source: string}>> {
  return page.evaluate(() =>
    [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')].map((c) => ({
      column: (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim(),
      source: (c.querySelector('[data-source]') as HTMLElement | null)?.getAttribute('data-source') ?? '',
    })));
}

// The drag is the gesture under test, so it is a real pointer stream: press on the grid header,
// eight intermediate moves to clear the drag threshold, then a wait on the panel's OWN
// "Add filter" drop zone as the readiness barrier before releasing inside it. Lifted from
// Viewers/FilterPanel/add-remove-entry-points-spec.ts:52 per the consolidation decision; the only
// change is the column. Substituting fg.updateOrAdd here would test the API, not the drop.
async function dragColumnHeaderToPanel(page: Page, column: string): Promise<void> {
  const src = await gridHeaderPoint(page, column);
  assertHeaderPointOnScreen(src, column);
  const panel = await page.evaluate(() => {
    const el = document.querySelector('[name="viewer-Filters"]');
    if (!el) return null;
    const r = el.getBoundingClientRect();
    return {cx: Math.round(r.x + r.width / 2), cy: Math.round(r.y + r.height / 2)};
  });
  expect(panel, 'the Filter Panel root is absent, so the drag has no target').not.toBeNull();
  await page.mouse.move(src!.x, src!.y);
  await page.mouse.down();
  for (let i = 1; i <= 8; i++) {
    await page.mouse.move(
      Math.round(src!.x + (panel!.cx - src!.x) * i / 8),
      Math.round(src!.y + (panel!.cy - src!.y) * i / 8), {steps: 3});
    await page.waitForTimeout(30);
  }
  try {
    await page.waitForFunction(() => [...document.querySelectorAll('.d4-drop-zone')]
      .some((z) => z.parentElement === document.body && (z.textContent ?? '').trim() === 'Add filter'),
    null, {timeout: 5000, polling: 100});
  }
  catch {
    await page.mouse.up();
    throw new Error(`dragColumnHeaderToPanel(${column}): the "Add filter" drop zone never appeared within 5s — ` +
      'the panel did not start accepting the drag');
  }
  const zone = await page.evaluate(() => {
    const dz = [...document.querySelectorAll('.d4-drop-zone')]
      .find((z) => z.parentElement === document.body && (z.textContent ?? '').trim() === 'Add filter')!;
    const r = dz.getBoundingClientRect();
    return {x: Math.round(r.x + r.width / 2), y: Math.round(r.y + r.height / 2)};
  });
  await page.mouse.move(zone.x, zone.y, {steps: 4});
  await page.waitForTimeout(120);
  await page.mouse.up();
  await page.waitForTimeout(700);
}

interface IntersectionPick {
  column: string;
  category: string;
  categoryRows: number;
  intersection: number;
}

// Picks the non-molecular category to cross the substructure filter with FROM THE DATA, so the
// step's inequalities hold by construction instead of by luck: the category must keep some
// substructure hits, drop some others (passOutside > 0 => intersection < substructure-only), and
// reach rows the substructure filter rejects (categoryRows > intersection => intersection <
// category-only). Choosing a category blind is how a cross-filter assertion ends up satisfied by a
// filter that intersects nothing.
async function pickIntersectingCategory(page: Page): Promise<IntersectionPick | null> {
  return page.evaluate((molCol) => {
    const df = grok.shell.tv.dataFrame;
    const passing: boolean[] = [];
    for (let i = 0; i < df.rowCount; i++) passing.push(df.filter.get(i));
    for (const c of df.columns) {
      if (c.name === molCol || c.type !== 'string') continue;
      const cats: string[] = [...c.categories].filter((x: string) => x !== null && x !== undefined && x !== '');
      if (cats.length < 2) continue;
      for (const cat of cats) {
        let intersection = 0;
        let categoryRows = 0;
        let passOutside = 0;
        for (let i = 0; i < df.rowCount; i++) {
          const isCat = String(c.get(i)) === cat;
          if (isCat) categoryRows++;
          if (isCat && passing[i]) intersection++;
          if (!isCat && passing[i]) passOutside++;
        }
        if (intersection > 0 && passOutside > 0 && categoryRows > intersection)
          return {column: c.name, category: cat, categoryRows, intersection};
      }
    }
    return null;
  }, CHEM_COL);
}

// A categorical filter card renders its categories through an INNER d4 grid, not a list
// (CatFilter extends GridFilterBase, core/client/d4/lib/src/viewers/filters/cat_filter.dart).
// Two consequences the earlier probe got wrong, and both made it click empty space: rows sit on
// a fixed pitch under a header band rather than dividing the card height by the category count,
// and the rendered order is the grid's, not `column.categories`. Geometry below matches the
// probe that works in Viewers/FilterPanel/filter-type-selection-modes-spec.ts:153.
const CARD_HEADER_BAND_H = 10;
const CARD_ROW_PITCH = 27;
const CARD_ROW_CENTRE = 13;
const CARD_X_NAME = 60;

interface CardGeometry {
  left: number;
  top: number;
  height: number;
  rows: number;
  viewportH: number;
  onScreen: boolean;
}

async function cardOverlayGeometry(page: Page, column: string): Promise<CardGeometry | null> {
  const geom = await page.evaluate(async (col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === col);
    if (!card) return null;
    // By this point the panel carries ~37 cards and is several viewport-heights tall, so the
    // card under test is usually scrolled out of sight. getBoundingClientRect still returns a
    // plausible rect for an off-screen card, which is how derived row points end up looking
    // fine while every click lands on nothing.
    card.scrollIntoView({block: 'center'});
    await new Promise((r) => setTimeout(r, 400));
    const overlay = card.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null;
    if (!overlay) return null;
    const r = overlay.getBoundingClientRect();
    return {
      left: r.left, top: r.top, height: r.height,
      rows: Math.max(0, Math.floor((r.height - 10) / 27)),
      viewportH: window.innerHeight,
      onScreen: r.top >= 0 && r.bottom <= window.innerHeight,
    };
  }, column);
  return geom;
}

async function selectedCategories(page: Page, column: string): Promise<string[]> {
  return page.evaluate((col) => {
    const st = grok.shell.tv.getFiltersGroup().getStates(col, 'categorical');
    const sel = st?.[0]?.selected;
    return Array.isArray(sel) ? sel.map((s: any) => String(s)) : [];
  }, column);
}

// Dispatched on the overlay element rather than driven with page.mouse, matching the probe that
// works in filter-type-selection-modes-spec.ts:191. A real pointer needs the point to be inside
// the viewport and unobstructed; an element-targeted event does not, which matters here because
// the card sits far down a 37-card panel. The scroll above is belt-and-braces for the same reason.
async function clickCardRow(page: Page, column: string, row: number): Promise<void> {
  await page.evaluate(({col, cx, cy}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === col);
    const overlay = card?.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null;
    if (!overlay) return;
    const r = overlay.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, view: window, button: 0,
      clientX: r.left + cx, clientY: r.top + cy};
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('click', o));
  }, {col: column, cx: CARD_X_NAME, cy: CARD_HEADER_BAND_H + CARD_ROW_PITCH * row + CARD_ROW_CENTRE});
}

// Resolves the row by READ-BACK instead of trusting an index: click a row, ask the filter group
// which category it selected, and keep the one that matches. Leaves the wanted category selected,
// so the caller must NOT click again — a second click would toggle it straight back off.
async function selectCategoryRow(page: Page, column: string, category: string):
  Promise<{row: number; geom: CardGeometry} | null> {
  const geom = await cardOverlayGeometry(page, column);
  expect(geom, `the ${column} card exposes no categorical grid overlay, so no category row can be clicked`)
    .not.toBeNull();
  expect(geom!.rows, `the ${column} card overlay is ${Math.round(geom!.height)} px tall, which holds no ` +
    'category rows at all — the walk below would have nothing to click').toBeGreaterThan(0);
  for (let row = 0; row < geom!.rows; row++) {
    await clickCardRow(page, column, row);
    await page.waitForTimeout(400);
    const sel = await selectedCategories(page, column);
    if (sel.length === 1 && sel[0] === category) return {row, geom: geom!};
    if (sel.length > 0) {
      // Toggle the wrong row back off before trying the next, so selections never accumulate.
      await clickCardRow(page, column, row);
      await page.waitForTimeout(400);
    }
  }
  return null;
}

interface CardEnabledState { found: boolean; disabled: boolean; checked: boolean | null; }

// Reads the enable state WITHOUT touching it — setCardEnabled toggles, so it cannot answer
// "is the card still switched off after the panel came back".
async function readCardEnabled(page: Page): Promise<CardEnabledState> {
  return page.evaluate(({col, source}) => {
    const card = (window as any).__filterCard(col, source) as HTMLElement | null;
    if (!card) return {found: false, disabled: false, checked: null};
    const cb = card.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null;
    return {found: true, disabled: card.classList.contains('d4-filter-disabled'), checked: cb ? cb.checked : null};
  }, {col: CHEM_COL, source: SUBSTRUCTURE_SOURCE});
}

test('Filter panel — Chem package filters', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const consoleErrors: string[] = [];
  page.on('console', (msg) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });
  const errorSet = () => new Set(consoleErrors.filter((t) => !t.includes(AMBIENT_CONSOLE_ERROR)));

  await loginToDatagrok(page);
  await page.evaluate((marker) => console.error(marker), CONSOLE_PROBE);
  await expect.poll(() => consoleErrors.filter((t) => t.includes(CONSOLE_PROBE)).length,
    {timeout: 15_000, intervals: [100, 200, 300, 500], message:
      'positive control: the console listener must receive a deliberately emitted error before any step relies on ' +
      'it — a detached listener would satisfy every "no new console error" assertion below by receiving nothing'})
    .toBeGreaterThan(0);

  await v.openTable(page, {path: CHEM_PATH});
  expect(await page.evaluate(() => grok.shell.tv.dataFrame.rowCount),
    'spgi-100 must open with its full row count').toBe(CHEM_FULL);

  const chem = await chemReadiness(page, CHEM_COL);
  console.log(`Chem readiness: ${chem.detail}`);
  test.skip(!chem.ok, `Chem is not usable on this build — ${chem.detail}`);

  await openChemFilterPanel(page);

  await installCardFinder(page);
  const molColumns = await moleculeColumns(page);
  expect(molColumns, 'spgi-100 must expose molecular columns for the Chem filters to bind to').toContain(CHEM_COL);

  let cardsAfterStep1 = 0;
  let step2Contains = 0;

  await softStep('Scenario 1, Step 1: Add Substructure Filter through the panel menu; record the card count', async () => {
    const removed = await removeSubstructureCards(page, molColumns);
    expect(removed, 'spgi-100 must open with substructure cards, otherwise their removal is not a real gesture')
      .toBe(molColumns.length);
    const cardsBefore = await cardCount(page);

    await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Substructure Filter...');
    const checked = await pickAllColumnsAndConfirm(page);
    expect(checked, 'the Select columns... picker must offer every molecular column').toBe(molColumns.length);
    await page.locator(CHEM_CARD).first().waitFor({state: 'attached', timeout: 20_000});
    const perColumn = await page.evaluate(({cols, type}) => {
      const fg = grok.shell.tv.getFiltersGroup();
      return cols.map((c: string) => (fg.getStates(c, type) || []).length);
    }, {cols: molColumns, type: SUBSTRUCTURE_FILTER});
    expect(perColumn, 'every molecular column must carry exactly one substructure filter card')
      .toEqual(molColumns.map(() => 1));
    cardsAfterStep1 = await cardCount(page);
    expect(cardsAfterStep1, 'confirming the picker must add exactly one card per checked column to the emptied panel')
      .toBe(cardsBefore + checked);
    expect(cardsAfterStep1, 'the panel must hold at least the substructure cards').toBeGreaterThanOrEqual(molColumns.length);
    console.log(`Cards after Add Filter > Substructure Filter...: ${cardsAfterStep1} (picker checked ${checked})`);
  });

  await softStep('Scenario 1, Step 2: commit the probe SMILES through the card sketcher\'s own field; Scenario 1, Step 3: the Filter Panel header counter and trueCount both report the engaged filter, and the card carries the structure', async () => {
    await openSketcherDialog(page);
    expect(await page.locator(`${SKETCHER_DIALOG} input[placeholder*="SMILES"]`).count(),
      'Step 2: the card\'s sketcher dialog must expose its own SMILES field — that field IS the gesture, ' +
      'the programmatic seeding in the Automation notes is only the fallback').toBeGreaterThan(0);
    await commitSmiles(page, CHEM_PROBE);
    await page.waitForFunction(
      (full) => grok.shell.tv.dataFrame.filter.trueCount < full, CHEM_FULL, {timeout: 30_000});
    const active = await waitForFilterIdle(page);
    expect(active).toBeLessThan(CHEM_FULL);
    expect(active, 'the probe must match at least one row, otherwise the partition test below is empty').toBeGreaterThan(0);
    const counter = await settledHeaderCounter(page);
    expect(counter.present,
      'Step 3: the Filter Panel header must expose its active-filter counter — the header channel cannot be read ' +
      'off an absent element').toBe(true);
    expect(counter.text.length,
      'Step 3: the header active-filter counter must render text — an empty counter satisfies any comparison ' +
      'below while reporting nothing').toBeGreaterThan(0);
    expect(counter.text,
      'Step 3: the header active-filter counter reports how many cards are FILTERING, not how many rows pass — ' +
      `exactly one card (the ${CHEM_COL} substructure card that just took the probe) is filtering, and it read ` +
      `"${counter.text}"`).toBe('1');
    const states = await substructureStates(page);
    expect(states.length, 'the sketched probe must be held by a substructure filter card').toBeGreaterThan(0);
    expect(states.some((s) => hasStructure(s.molBlock)),
      'the substructure card must carry the committed structure after the sketch').toBe(true);
    step2Contains = active;
    console.log(`Chem substructure active trueCount: ${active}`);
  });

  let containsCount = 0;

  await softStep('Scenario 1, Step 4: Contains + Not contains partition the dataset exactly (compound invariant)', async () => {
    await revealSearchOptions(page);
    await installFilterPassCounter(page);
    const types = ['Contains', 'Not contains', 'Exact', 'Similar', 'Included in', 'Not included in'];
    const counts: Record<string, number> = {};
    for (const t of types) {
      const shown = await readSearchType(page);
      const passesBefore = await filterPasses(page);
      const applied = await setSearchType(page, t);
      expect(applied, `search type "${t}" must be offered by the Structure card's search-type control`).toBe(true);
      counts[t] = await countForSearchType(page, t, shown === t ? null : passesBefore);
      expect(counts[t], `trueCount for "${t}" must be a real count`).toBeGreaterThanOrEqual(0);
      expect(counts[t], `trueCount for "${t}" must not exceed the full row count`).toBeLessThanOrEqual(CHEM_FULL);
    }
    console.log(`Chem search-type counts: ${JSON.stringify(counts)}`);
    expect(step2Contains > 0 && step2Contains < CHEM_FULL,
      `the narrowing count Step 2 proved for the same probe in the same mode must be a real partial row set, or the ` +
      `Contains pin below anchors on nothing (Step 2 read ${step2Contains})`).toBe(true);
    expect(counts['Contains'],
      'Contains must reproduce exactly the narrowing row count Step 2 already proved for the same probe in the same ' +
      'mode — every partition sum, subset and superset claim below is satisfied by a filter that matches every row ' +
      `(Contains 100 / Not contains 0), so Contains is pinned first (Step 2 read ${step2Contains})`)
      .toBe(step2Contains);
    expect(counts['Contains'] + counts['Not contains'],
      'Contains and Not contains must partition the dataset exactly').toBe(CHEM_FULL);
    expect(counts['Included in'] + counts['Not included in'],
      'Included in and Not included in must partition the dataset exactly').toBe(CHEM_FULL);
    expect(counts['Exact'], 'Exact matches are a subset of Contains').toBeLessThanOrEqual(counts['Contains']);
    expect(counts['Exact'], 'Exact matches are a subset of Included in').toBeLessThanOrEqual(counts['Included in']);
    expect(counts['Similar'], 'Similar (cutoff 0.8) must include every Exact match').toBeGreaterThanOrEqual(counts['Exact']);
    const restoreShown = await readSearchType(page);
    const restorePasses = await filterPasses(page);
    expect(await setSearchType(page, 'Contains')).toBe(true);
    containsCount = await countForSearchType(page, 'Contains', restoreShown === 'Contains' ? null : restorePasses);
    expect(containsCount, 'restoring Contains must reproduce the Contains count measured above').toBe(counts['Contains']);
  });

  await softStep('Scenario 1, Step 5: GROK-18383: add Scaffold Tree Filter — one new card, narrows further', async () => {
    const before = await cardCount(page);
    expect(before, 'no card may appear between Step 1 and Step 5 — the GROK-18383 delta is measured against Step 1')
      .toBe(cardsAfterStep1);

    await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Scaffold Tree Filter...');
    const checked = await pickAllColumnsAndConfirm(page);
    expect(checked, 'the Scaffold Tree picker must offer every molecular column').toBe(molColumns.length);
    await waitForCardCount(page, before + checked, 30_000);
    expect(await cardCount(page)).toBe(before + checked);
    const scaffoldPerColumn = await page.evaluate(({cols, type}) => {
      const fg = grok.shell.tv.getFiltersGroup();
      return cols.map((c: string) => (fg.getStates(c, type) || []).length);
    }, {cols: molColumns, type: SCAFFOLD_TREE_FILTER});
    expect(scaffoldPerColumn, 'GROK-18383: every picked column must carry exactly one scaffold-tree card')
      .toEqual(molColumns.map(() => 1));

    const scaffoldCard = page.locator(SCAFFOLD_CARD).first();
    await scaffoldCard.waitFor({state: 'attached', timeout: 20_000});
    expect(await page.locator(SCAFFOLD_CARD).count(),
      `exactly one scaffold-tree card must sit on ${CHEM_COL}, or the gestures below address an arbitrary card`).toBe(1);
    await scaffoldCard.locator('.chem-scaffold-tree-toolbar .fa-plus').first().click();
    await page.locator(ADD_SCAFFOLD_DIALOG).first().waitFor({timeout: 20_000});
    await enterSmiles(page, SCAFFOLD, ADD_SCAFFOLD_DIALOG);
    await page.locator(`${ADD_SCAFFOLD_DIALOG} [name="button-Add"]`).first().click();
    await page.locator(ADD_SCAFFOLD_DIALOG).first().waitFor({state: 'detached', timeout: 30_000});
    const node = scaffoldCard.locator('.d4-tree-view-group-label .chem-mol-box .mol-host').first();
    await node.waitFor({timeout: 45_000});
    expect(await scaffoldCard.locator('.d4-tree-view-group-label .chem-mol-box').count(),
      'the sketched scaffold must land as exactly one node in the tree card').toBe(1);
    await node.click();
    await expect.poll(() => scaffoldCard.locator('.d4-tree-view-node input[type="checkbox"]:checked').count(),
      {timeout: 30_000, intervals: [200, 300, 500, 1000], message:
        'selecting the scaffold node in the tree card must leave that node CHECKED — an unchecked node is not a ' +
        'criterion and the card would not be filtering at all'}).toBe(1);
    await waitForFilterIdle(page);
    expect(await cardCount(page),
      'selecting a scaffold must update the menu-added card, not add a second one').toBe(before + checked);

    const isolated = await setCardEnabled(page, false);
    expect(isolated, 'unchecking the substructure card must disable it, so the scaffold effect can be read on its own').toBe(true);
    const scaffoldOnly = await waitForFilterIdle(page);
    expect(scaffoldOnly,
      'GROK-18383: the scaffold tree card must actually filter — with the substructure card off it must still narrow below the full row count')
      .toBeLessThan(CHEM_FULL);
    expect(scaffoldOnly, 'the scaffold must keep at least one row, otherwise "narrowing" is indistinguishable from a broken filter')
      .toBeGreaterThan(0);
    expect(await setCardEnabled(page, true),
      're-arming the substructure card must re-enable it before the AND of both cards is read').toBe(true);
    const combined = await waitForFilterIdle(page);
    console.log(`Scaffold-tree: contains=${containsCount} scaffoldOnly=${scaffoldOnly} combined=${combined}`);
    expect(combined, 'the AND of both cards cannot pass more rows than the substructure card alone')
      .toBeLessThanOrEqual(containsCount);
  });

  await softStep('Scenario 1, Step 6: GROK-20001: Current Value > Use as filter adds exactly one card', async () => {
    const removed = await removeSubstructureCards(page, molColumns);
    expect(removed, 'at least one substructure card must exist to be removed before the Use-as-filter gesture').toBeGreaterThan(0);
    const before = await cardCount(page);

    await chemReadiness(page, CHEM_COL);

    const pt = await gridCellPoint(page, CHEM_COL, 0);
    expect(pt, `a ${CHEM_COL} cell in row 0 must be locatable in the grid for the right-click`).not.toBeNull();
    await page.mouse.move(pt!.x, pt!.y, {steps: 4});
    await page.mouse.click(pt!.x, pt!.y, {button: 'right'});
    await page.locator('.d4-menu-popup').first().waitFor({state: 'visible', timeout: 10_000});
    const clicked = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {row: df.currentRowIdx, col: df.currentCol?.name ?? null};
    });
    expect(clicked.col,
      `the right-click must make a ${CHEM_COL} cell current before the menu item is chosen — that menu item reads ` +
      'the CURRENT cell, and nothing here set that cell but the click itself').toBe(CHEM_COL);
    console.log(`Use-as-filter click: aimed at grid row 0 = table row ${pt!.tableRow}, current row after click = ${clicked.row}`);
    expect(clicked.row,
      `the right-click must make current the very cell it was aimed at — the product's own hit test placed grid row 0 ` +
      `of ${CHEM_COL} at table row ${pt!.tableRow}, and a click landing on a neighbour would still yield a valid ` +
      'structure, so only this identity rules the miss out').toBe(pt!.tableRow);
    await menuLeaf(page, 'Current Value', 'Use as filter');

    await page.waitForFunction(({col, type}) =>
      (grok.shell.tv.getFiltersGroup().getStates(col, type) || []).length > 0,
    {col: CHEM_COL, type: SUBSTRUCTURE_FILTER}, {timeout: 20_000, polling: 200});
    const active = await waitForFilterIdle(page);
    expect(await cardCount(page),
      'GROK-20001: Use as filter must add exactly ONE card to the panel, not two').toBe(before + 1);
    const states = await substructureStates(page);
    expect(states.length, 'GROK-20001: exactly one substructure card, not two').toBe(1);
    expect(hasStructure(states[0].molBlock),
      'the card added by Use as filter must carry the clicked cell\'s structure').toBe(true);
    expect(await page.evaluate((row) => grok.shell.tv.dataFrame.filter.get(row), clicked.row),
      `the card must carry the structure of the row the click made current (row ${clicked.row}), so that row must ` +
      'itself still pass the filter the card installed — a card built from a different cell need not keep it')
      .toBe(true);
    expect(active).toBeLessThan(CHEM_FULL);
    expect(active, 'the cell\'s own molecule must match at least itself').toBeGreaterThan(0);
  });

  const rebuildFullSubstructurePanel = async (why: string) => {
    await removeAllViaHamburger(page);
    await closeFilterPanel(page);
    await reopenFilterPanel(page);
    await page.waitForFunction(({cols, type}) => {
      const fg = grok.shell.tv.getFiltersGroup();
      return cols.every((c: string) => (fg.getStates(c, type) || []).length === 1)
        && document.querySelectorAll('[name="viewer-Filters"] .d4-filter .chem-filter').length === cols.length;
    }, {cols: molColumns, type: SUBSTRUCTURE_FILTER}, {timeout: 60_000, polling: 250})
      .catch(() => undefined);
    const shape = await page.evaluate(({cols, type}) => {
      const fg = grok.shell.tv.getFiltersGroup();
      return {
        perColumn: cols.map((c: string) => (fg.getStates(c, type) || []).length),
        chemBodies: document.querySelectorAll('[name="viewer-Filters"] .d4-filter .chem-filter').length,
        cards: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
      };
    }, {cols: molColumns, type: SUBSTRUCTURE_FILTER});
    const targetCards = await page.locator(CHEM_CARD).count();
    expect(shape.perColumn,
      `${why}: emptying the panel with Remove All, closing it and reopening it must bring every one of the ` +
      `${molColumns.length} molecular columns (${molColumns.join(', ')}) back carrying exactly ONE substructure ` +
      'state — the configuration the scenario describes is GUARANTEED here, never inherited from what the earlier ' +
      `steps happened to leave (read ${JSON.stringify(shape.perColumn)})`)
      .toEqual(molColumns.map(() => 1));
    expect(shape.chemBodies,
      `${why}: the rebuilt panel must render one substructure card body per molecular column ` +
      `(${shape.chemBodies} bodies against ${molColumns.length} columns), so a half-built panel cannot be ` +
      'measured as a whole one').toBe(molColumns.length);
    expect(targetCards,
      `${why}: exactly one of those bodies must sit on ${CHEM_COL}, or every card gesture below addresses an ` +
      'arbitrary card').toBe(1);
    expect((await substructureStates(page)).length,
      `${why}: the ${CHEM_COL} column must carry exactly one substructure state after the rebuild`).toBe(1);
    console.log(`${why} precondition: ${shape.cards} cards in the panel, ${shape.chemBodies} substructure bodies ` +
      `on ${molColumns.length} molecular columns (${molColumns.join(', ')}), ${targetCards} on ${CHEM_COL}`);
  };

  const armStructureCard = async (why: string) => {
    await rebuildFullSubstructurePanel(why);
    expect((await substructureStates(page)).map((s) => hasStructure(s.molBlock)),
      `${why}: the rebuilt panel must come up with no structure on ${CHEM_COL} — the probe below is what puts one ` +
      'there, so an inherited structure would make the arming gesture untested').toEqual([false]);
    await openSketcherDialog(page);
    await commitSmiles(page, CHEM_PROBE);
    await waitForFilterIdle(page);
    expect((await substructureStates(page)).map((s) => hasStructure(s.molBlock)),
      `${why}: the ${CHEM_COL} card must carry the probe structure before anything is measured against it`)
      .toEqual([true]);
  };

  await softStep('Scenario 1, Step 7: GROK-18530: cloned view + Similar mode agree with the original', async () => {

    await armStructureCard('Step 7');
    await revealSearchOptions(page);
    expect(await setSearchType(page, 'Similar'), 'the Structure card must offer the Similar search type').toBe(true);
    await waitForFilterIdle(page);
    expect(await setSimilarityCutoff(page, SIMILAR_CUTOFF),
      'the Similar mode must expose its similarity-cutoff editor').toBe(true);
    await waitForFilterIdle(page);
    const orig = (await substructureStates(page))[0];
    expect(orig, 'the original view must hold a substructure filter state before cloning').toBeTruthy();
    expect(orig.searchType, 'the original card must really be in Similar mode before cloning').toBe('Similar');
    expect(orig.simCutOff, 'the original card must really carry the non-default cutoff before cloning').toBe(SIMILAR_CUTOFF);
    expect(hasStructure(orig.molBlock), 'the original card must carry a molecule before cloning').toBe(true);

    expect(await tableViewCount(page),
      'exactly one table view must be open before the clone, or "the second view" is not the one this step made')
      .toBe(1);
    expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Clone View']),
      'GROK-18530: the clone must be driven through the product\'s own View > Layout > Clone View command — that ' +
      'command saves the view layout and rebuilds the view from it (layout.dart:97-114), and the filter-state ' +
      'serialization it runs is the code the bug lived in; opening a second view over the table bypasses it entirely')
      .toBe(true);
    await page.waitForFunction(({col, type}) => {
      const views: any[] = [];
      for (const t of grok.shell.tableViews) views.push(t);
      if (views.length < 2) return false;
      return ((views[1].getFiltersGroup().getStates(col, type)) || []).length > 0;
    }, {col: CHEM_COL, type: SUBSTRUCTURE_FILTER}, {timeout: 60_000, polling: 300});
    await waitForFilterIdle(page);

    const distinct = await page.evaluate(({col, type}) => {
      const views: any[] = [];
      for (const t of grok.shell.tableViews) views.push(t);
      if (views.length < 2) return null;
      const g0 = views[0].getFiltersGroup();
      const g1 = views[1].getFiltersGroup();
      const s0 = (g0.getStates(col, type) || [])[0] ?? null;
      const s1 = (g1.getStates(col, type) || [])[0] ?? null;
      return {
        views: views[0].dart !== views[1].dart,
        groups: !!g0.dart && !!g1.dart && g0.dart !== g1.dart,
        states: !!s0 && !!s1 && s0 !== s1,
      };
    }, {col: CHEM_COL, type: SUBSTRUCTURE_FILTER});
    expect(distinct, 'both table views must be reachable after the clone').not.toBeNull();
    expect(distinct, 'the clone must be a SECOND view with its OWN filter group holding its OWN substructure state — ' +
      'without this every field compared below could be the original\'s own value read twice')
      .toEqual({views: true, groups: true, states: true});

    const clone = (await substructureStates(page, 1))[0];
    expect(clone, 'the clone must hold a substructure filter state').toBeTruthy();
    expect(clone.searchType, 'GROK-18530: the clone must come up in Similar mode').toBe('Similar');
    expect(clone.simCutOff, 'GROK-18530: the clone must carry the original non-default cutoff').toBe(SIMILAR_CUTOFF);
    expect(typeof clone.fp, 'the clone must carry a fingerprint type').toBe('string');
    expect(clone.fp, 'the clone must carry the original fingerprint type').toBe(orig.fp);
    expect(hasStructure(clone.molBlock), 'the clone must carry a molecule').toBe(true);
    expect(clone.molBlock, 'the clone must carry the original molecule').toBe(orig.molBlock);
    expect(await readSearchType(page, 0), 'the original card must display Similar').toBe('Similar');
    expect(await readSearchType(page, 1), 'the cloned card must display Similar').toBe('Similar');
  });

  await softStep('Scenario 1, Step 8: Reset clears filtering and the sketched substructure state', async () => {
    await page.evaluate(() => {
      const tvs: any[] = [];
      for (const t of grok.shell.tableViews) tvs.push(t);
      if (tvs.length > 1) tvs[tvs.length - 1].close();
    });
    await page.waitForFunction(() => {
      const tvs: any[] = [];
      for (const t of grok.shell.tableViews) tvs.push(t);
      return tvs.length === 1;
    }, null, {timeout: 20_000, polling: 200});
    await expect.poll(() => page.locator('[name="viewer-Filters"]').count(),
      {timeout: 30_000, intervals: [200, 300, 500, 1000], message:
        'the clone\'s own Filter Panel must leave the document when the clone view is closed — every card reading ' +
        'below is document-wide, so a lingering second panel would be counted as the original\'s'}).toBe(1);

    const beforeReset = await substructureStates(page);
    expect(beforeReset.length, 'a substructure card must be present before the reset').toBeGreaterThan(0);
    expect(beforeReset.some((s) => hasStructure(s.molBlock)),
      'a drawn structure must be present before the reset — otherwise the clearing check is vacuous').toBe(true);
    expect(beforeReset.some((s) => s.searchType !== 'Contains'),
      'a non-default search mode must be in place before the reset').toBe(true);
    expect(await trueCount(page), 'the panel must be filtering before the reset').toBeLessThan(CHEM_FULL);
    const cardsBeforeReset = await cardCount(page);
    const bodiesBeforeReset = await cardBodies(page);
    const drawnBeforeReset = classifyBodies(bodiesBeforeReset);
    console.log(bodyTable('Step 8 card bodies before the reset:', bodiesBeforeReset));
    expect(drawnBeforeReset.drawn,
      'at least one card body must actually PAINT a structure before the reset — read as non-background pixels on ' +
      'the body\'s own canvas, since a canvas ELEMENT carrying nothing satisfies element-presence while showing ' +
      `nothing; bodies read: ${drawnBeforeReset.detail}`).toBeGreaterThan(0);

    const resetIcon = page.locator(fp.RESET_CRITERIA_ICON).first();
    await resetIcon.click();
    await page.waitForFunction((full) => grok.shell.tv.dataFrame.filter.trueCount === full,
      CHEM_FULL, {timeout: 30_000});
    const after = await waitForFilterIdle(page);
    expect(after, 'the reset must return every row to the filter').toBe(CHEM_FULL);

    const afterStates = await substructureStates(page);
    expect(afterStates.length, 'reset keeps the cards — it clears their criteria').toBe(beforeReset.length);
    expect(afterStates.map((s) => hasStructure(s.molBlock)),
      'no drawn structure may survive the reset').toEqual(afterStates.map(() => false));
    expect(afterStates.map((s) => s.searchType),
      'reset restores the default Contains mode (chem-substructure-filter.ts:295-298)')
      .toEqual(afterStates.map(() => 'Contains'));
    expect(await cardCount(page),
      `reset clears criteria; it must not remove cards — the panel held ${cardsBeforeReset} cards before it`)
      .toBe(cardsBeforeReset);
    const clearing = await observeSketchClearing(page, SKETCH_CLEAR_OBSERVE_MS);
    console.log(`Step 8 sketch-area clearing: ${clearing.clearedAfterMs === null
      ? `NOT cleared — ${clearing.lastDrawn} card body/bodies (${clearing.lastDrawnColumns.join(', ')}) were still ` +
        `painting a structure and ${clearing.lastUnreadableColumns.length} were unreadable ` +
        `(${clearing.lastUnreadableColumns.join(', ') || 'none'}) across the whole ${clearing.observedMs} ms ` +
        `observation window (${clearing.samples} reads)`
      : `cleared after ${clearing.clearedAfterMs} ms and held clear across three consecutive reads ` +
        `(${clearing.unreadableSamples} of ${clearing.samples} reads found an unreadable canvas)`}`);
    // The cleared reading is re-taken after a held interval: the run above ends 1.5 s after the
    // first clean sample, and the claim is that the sketch area STAYS cleared.
    await page.waitForTimeout(SKETCH_CLEAR_HOLD_MS);
    const held = await sketchAreaBodies(page);
    expect({drawn: held.drawn, unreadable: held.unreadable},
      `the sketch areas must STILL be cleared, and still readable, when the bodies are re-read ` +
      `${SKETCH_CLEAR_HOLD_MS} ms after the clearing run ended — a body that repaints or goes unreadable after ` +
      `the run is exactly what a 1.5 s run cannot see (read: ${held.detail})`).toEqual({drawn: 0, unreadable: 0});
    const bodiesAfterReset = await cardBodies(page);
    const cleared = classifyBodies(bodiesAfterReset);
    console.log(bodyTable('Step 8 card bodies after the reset:', bodiesAfterReset));
    console.log(`Step 8 card-body shapes after the reset (GROK-20739, open — recorded, not asserted): ` +
      `${cleared.placeholder} back to the "Sketch" placeholder, ${cleared.blankCanvas} left holding a blank ` +
      `canvas with no placeholder (${cleared.detail})`);
    expect(cleared.total,
      'substructure card bodies must still be in the panel after the reset, or "none of them shows a structure" ' +
      'is true of a panel that has none').toBeGreaterThan(0);
    expect(bodiesAfterReset.filter((b) => b.canvas).length + cleared.placeholder,
      'every substructure card body must still expose its sketch area after the reset — either a canvas or the ' +
      `"Sketch" placeholder — or the pixel channel below is read off nothing (${cleared.detail})`).toBe(cleared.total);
    expect(cleared.unreadableColumns,
      'every card body that still carries a canvas must be READABLE after the reset — an unmeasurable canvas ' +
      '(0x0, or no 2d context) leaves ink at -1, which the "no body still paints a structure" reading below ' +
      'would count as cleared. The reset resizes this canvas, so an unreadable body is a real mid-teardown ' +
      'state and must fail here rather than pass there').toEqual([]);
    expect(cleared.drawn,
      'no substructure card body may still PAINT a drawn structure after the reset — the sketch area must be ' +
      'cleared, not merely detached from filtering (GROK-14028 family: the sketcher input persisted after reset). ' +
      'This is read as non-background pixels on each body\'s own canvas, never as the presence of a canvas ' +
      `element, which survives the reset unpainted; still painting: ${cleared.drawnColumns.join(', ') || '(none)'} ` +
      `— the panel was watched deliberately for ${clearing.observedMs} ms (${clearing.samples} reads, one every ` +
      `500 ms) and ${clearing.clearedAfterMs === null ? 'never cleared within it'
        : `cleared after ${clearing.clearedAfterMs} ms`}`)
      .toBe(0);
    expect(bodiesAfterReset.filter((b) => hasStructure(b.molBlock)).length,
      'no substructure card may still hold a structure in its saved state after the reset — the pixel channel and ' +
      `the saved-state channel must agree that the structure is gone (${bodiesAfterReset
        .map((b) => `${b.column}=${molAtomCount(b.molBlock)} atoms/${(b.molBlock ?? '').length} chars`)
        .join(', ')})`).toBe(0);

    const clearedCounter = await settledHeaderCounter(page);
    expect(clearedCounter.present,
      'the header active-filter counter must still be in the panel header after the reset').toBe(true);
    expect(clearedCounter.text,
      'the header active-filter counter must fall to 0 once the reset leaves no card filtering — Step 3 only ever ' +
      'read it at 1, which a stale or hardcoded indicator satisfies; the reset is the contrast that discriminates')
      .toBe('0');
  });

  await softStep('Scenario 1, Step 9: github-3445: realignment repaints on the spot; "Filter as you draw" gates the sketch', async () => {
    await rebuildFullSubstructurePanel('Step 9');

    const quiesce = async (why: string) => {
      await waitForFilterIdle(page);
      await v.waitForCanvasQuiet(page, 'Grid', {canvasSelector: GRID_CANVAS, stableReads: 3});
      expect(await v.snapshotCanvasColors(page, 'Grid', GRID_CANVAS),
        `the grid canvas must be readable for ${why}`).toBe(true);
      await page.waitForTimeout(1500);
      expect((await v.diffCanvasColors(page, 'Grid', GRID_CANVAS)).deltaPx,
        `negative control: the grid must be pixel-stable before ${why}`).toBe(0);
    };

    const clearDrawnStructure = async (why: string) => {
      const resetIcon = page.locator(fp.RESET_CRITERIA_ICON).first();
      await resetIcon.click();
      await page.waitForFunction((full) => grok.shell.tv.dataFrame.filter.trueCount === full,
        CHEM_FULL, {timeout: 30_000});
      await waitForFilterIdle(page);
      const states = await substructureStates(page);
      expect(states.length, `a substructure card must survive the reset that precedes ${why}`).toBeGreaterThan(0);
      expect(states.map((s) => hasStructure(s.molBlock)),
        `no drawn structure may be left when ${why} enters its own probe`).toEqual(states.map(() => false));
    };

    let appliedAfterOk = 0;

    const countAfterChangeFrom = async (from: number, timeoutMs: number): Promise<number> => {
      await page.waitForFunction((n) => grok.shell.tv.dataFrame.filter.trueCount !== n,
        from, {timeout: timeoutMs, polling: 200})
        .catch(() => undefined);
      return waitForFilterIdle(page);
    };

    const heldNoPropagation = async (label: string, expected: number) => {
      for (const waitMs of [3000, 4000]) {
        await page.waitForTimeout(waitMs);
        expect((await v.diffCanvasColors(page, 'Grid', GRID_CANVAS)).deltaPx,
          `${label}: with "Filter as you draw" CLEARED, the sketch must not repaint the structure ` +
          `column while the dialog stays open (read ${waitMs} ms after the sketch)`).toBe(0);
        expect(await trueCount(page),
          `${label}: with "Filter as you draw" CLEARED, the sketch must not move the row set while ` +
          `the dialog stays open (read ${waitMs} ms after the sketch; expected ${expected})`).toBe(expected);
      }
    };

    await clearDrawnStructure('Half A');
    await openSketcherDialog(page);
    await page.locator(`${SKETCHER_DIALOG} .chem-sketcher-filter-options`).first().waitFor({timeout: 15_000});
    const fydOff = await setFilterAsYouDraw(page, false);
    expect(fydOff.set, 'the sketcher must expose the "Filter as you draw" checkbox (Half A)').toBe(true);
    expect(fydOff.value,
      '"Filter as you draw" must read CLEARED after being set so, BEFORE the probe is entered (Half A) — never inherited').toBe(false);
    await expectAlignHighlight(page, false, 'Half A, opening probe');
    await quiesce('Half A, the opening probe (withheld from the grid until OK)');

    await enterSmiles(page, CHEM_PROBE);
    await heldNoPropagation('Half A, opening probe', CHEM_FULL);
    const okDraw = page.locator(`${SKETCHER_DIALOG} [name="button-OK"]`).first();
    expect(await okDraw.count(), 'the sketcher OK button must be present to apply the Half A probe').toBeGreaterThan(0);
    await okDraw.click();
    await page.locator(SKETCHER_DIALOG).first().waitFor({state: 'detached', timeout: 20_000});
    const filteredA = await countAfterChangeFrom(CHEM_FULL, 60_000);
    expect(filteredA,
      'Half A: confirming the sketcher must engage the substructure filter — an unfiltered grid is not a valid stage for either claim')
      .toBeLessThan(CHEM_FULL);
    expect(filteredA,
      'Half A: the filter must keep rows on screen — an empty grid cannot show a realignment').toBeGreaterThan(0);

    await openSketcherDialog(page);
    await page.locator(`${SKETCHER_DIALOG} .chem-sketcher-filter-options`).first().waitFor({timeout: 15_000});
    const fydStillOff = await setFilterAsYouDraw(page, false);
    expect(fydStillOff.set, 'the reopened sketcher must expose the "Filter as you draw" checkbox (Half A)').toBe(true);
    expect(fydStillOff.value,
      '"Filter as you draw" must still read CLEARED on the reopened dialog (Half A) — the reopen must not change the state under test').toBe(false);
    await expectAlignHighlight(page, false, 'Half A, Claim 1 (the ON flip must be the only cause inside the window)');
    await quiesce('Half A, Claim 1 (realignment with the dialog open)');

    let errorsBefore = errorSet();
    await expectAlignHighlight(page, true, 'Half A, Claim 1 (turning them ON with the dialog open)');
    const deltaAlign = await v.waitForCanvasChange(page, 'Grid',
      {minDelta: 1, timeoutMs: 30_000, canvasSelector: GRID_CANVAS});
    expect(deltaAlign,
      'Claim 1 (github-3445): turning Align + Highlight on over an already-applied structure must repaint ' +
      'the rendered structure column immediately, with the sketcher dialog still open — the alignment path ' +
      'is not gated by "Filter as you draw"').toBeGreaterThan(0);
    expect(await waitForFilterIdle(page),
      `Claim 1: the realignment must repaint only the drawing, not change the row set (expected ${filteredA})`)
      .toBe(filteredA);
    let newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
    expect(newErrors, `Claim 1: no new console error text across the realignment; got: ${newErrors.join(' | ')}`).toEqual([]);
    console.log(`Claim 1: realignment delta ${deltaAlign} px, trueCount held at ${filteredA}`);

    await expectAlignHighlight(page, false, 'Half A, Claim 2 (each claim starts from Align and Highlight OFF)');
    await quiesce('Half A, Claim 2 (further edit with the dialog open)');
    errorsBefore = errorSet();
    await enterSmiles(page, CHEM_PROBE_2);
    expect(await page.locator(SKETCHER_DIALOG).count(),
      'Half A: the sketcher dialog must still be open — the withheld-edit claim only holds while it is up')
      .toBeGreaterThan(0);
    await heldNoPropagation('Claim 2', filteredA);

    const okA = page.locator(`${SKETCHER_DIALOG} [name="button-OK"]`).first();
    expect(await okA.count(), 'the sketcher OK button must be present to confirm the Half A edit').toBeGreaterThan(0);
    await okA.click();
    await page.locator(SKETCHER_DIALOG).first().waitFor({state: 'detached', timeout: 20_000});
    appliedAfterOk = await countAfterChangeFrom(filteredA, 60_000);
    expect(appliedAfterOk,
      `Claim 2: confirming the sketcher with OK must apply the withheld edit — the row set must move off ${filteredA}`)
      .not.toBe(filteredA);
    expect(appliedAfterOk,
      `Claim 2: the confirmed edit must still filter (0 < rows < ${CHEM_FULL}), got ${appliedAfterOk}`)
      .toBeGreaterThan(0);
    expect(appliedAfterOk,
      `Claim 2: the confirmed edit must still filter (0 < rows < ${CHEM_FULL}), got ${appliedAfterOk}`)
      .toBeLessThan(CHEM_FULL);
    newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
    expect(newErrors, `Claim 2 (cleared): no new console error text across the withheld edit and its confirmation; got: ${newErrors.join(' | ')}`).toEqual([]);
    console.log(`Claim 2 (cleared): held at ${filteredA} with the dialog open, ${appliedAfterOk} after OK`);

    await clearDrawnStructure('Half B');
    await openSketcherDialog(page);
    await page.locator(`${SKETCHER_DIALOG} .chem-sketcher-filter-options`).first().waitFor({timeout: 15_000});
    const fydOn = await setFilterAsYouDraw(page, true);
    expect(fydOn.set, 'the sketcher must expose the "Filter as you draw" checkbox (Half B)').toBe(true);
    expect(fydOn.value,
      '"Filter as you draw" must read CHECKED after being set so, BEFORE the probe is entered (Half B) — never inherited').toBe(true);
    await expectAlignHighlight(page, false, 'Half B, opening probe');

    await enterSmiles(page, CHEM_PROBE);
    const filteredB = await countAfterChangeFrom(CHEM_FULL, 60_000);
    expect(filteredB,
      'Half B: with "Filter as you draw" CHECKED the probe must filter as it is entered, with no OK needed')
      .toBeLessThan(CHEM_FULL);
    expect(filteredB,
      'Half B: the filter must keep rows on screen — an empty grid can show neither an edit nor a repaint')
      .toBeGreaterThan(0);
    expect(filteredB,
      `Half B must start from the same row set as Half A (A=${filteredA}, B=${filteredB}) — the same probe in the same mode`).toBe(filteredA);
    expect(await page.locator(SKETCHER_DIALOG).count(),
      'Half B: the sketcher dialog must still be open — a live update is only the claim while the dialog is up').toBeGreaterThan(0);
    await quiesce('Half B, Claim 1 (realignment with the dialog open, "Filter as you draw" CHECKED)');

    errorsBefore = errorSet();
    await expectAlignHighlight(page, true, 'Half B, Claim 1 (turning them ON with the dialog open)');
    const deltaAlignB = await v.waitForCanvasChange(page, 'Grid',
      {minDelta: 1, timeoutMs: 30_000, canvasSelector: GRID_CANVAS});
    expect(deltaAlignB,
      'Claim 1 (github-3445), Half B: the realignment repaint is NOT gated by "Filter as you draw" — with the ' +
      'checkbox CHECKED, turning Align + Highlight on over an already-applied structure must repaint the rendered ' +
      'structure column immediately, exactly as it did with the checkbox cleared').toBeGreaterThan(0);
    expect(await waitForFilterIdle(page),
      `Claim 1, Half B: the realignment must repaint only the drawing, not change the row set (expected ${filteredB})`)
      .toBe(filteredB);
    newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
    expect(newErrors,
      `Claim 1 (Half B): no new console error text across the realignment; got: ${newErrors.join(' | ')}`).toEqual([]);
    console.log(`Claim 1 (Half B): realignment delta ${deltaAlignB} px, trueCount held at ${filteredB}`);

    await expectAlignHighlight(page, false, 'Half B, Claim 2 (each claim starts from Align and Highlight OFF)');
    await quiesce('Half B, Claim 2 (further edit with the dialog open)');

    errorsBefore = errorSet();
    await enterSmiles(page, CHEM_PROBE_2);
    const liveB = await countAfterChangeFrom(filteredB, 45_000);
    expect(await page.locator(SKETCHER_DIALOG).count(),
      'Half B: the sketcher dialog must still be open — a live update is only the claim while it is up')
      .toBeGreaterThan(0);
    expect(liveB,
      `Claim 2 (checked): with "Filter as you draw" CHECKED, a further edit must reach the grid with the ` +
      `dialog still open — the row set must move off ${filteredB}, and it read ${liveB}`).not.toBe(filteredB);
    expect(liveB,
      `Claim 2 (checked): the live edit must reach the SAME row set the cleared half reached only on OK ` +
      `(after OK: ${appliedAfterOk}, live: ${liveB}) — the two halves differ in timing, not in outcome`)
      .toBe(appliedAfterOk);
    const deltaLive = (await v.diffCanvasColors(page, 'Grid', GRID_CANVAS)).deltaPx;
    expect(deltaLive,
      `Claim 2 (checked): the live edit must repaint the grid as well as move the row set, got delta ${deltaLive}`)
      .toBeGreaterThan(0);
    newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
    expect(newErrors, `Claim 2 (checked): no new console error text across the live edit; got: ${newErrors.join(' | ')}`).toEqual([]);
    console.log(`Claim 2 (checked): ${filteredB} -> ${liveB} live with the dialog open, delta ${deltaLive} px`);

    const okB = page.locator(`${SKETCHER_DIALOG} [name="button-OK"]`).first();
    if (await okB.count() > 0) await okB.click();
  });

  await softStep('Scenario 1 Step 10: GROK-14952: the panel filter and the column-popup filter on the same column do not diverge', async () => {
    const staleOk = page.locator(`${SKETCHER_DIALOG} [name="button-OK"]`).first();
    if (await staleOk.count() > 0) {
      await staleOk.click();
      await page.locator(SKETCHER_DIALOG).first().waitFor({state: 'detached', timeout: 20_000}).catch(() => undefined);
    }

    await removeSubstructureCards(page, molColumns);
    const cardsBefore = await cardCount(page);
    expect(await sustainedCount(page, 'the baseline before either surface is armed'),
      'no card may be filtering when the popup surface is armed — otherwise the narrowing read below is not the popup\'s')
      .toBe(CHEM_FULL);

    const popup = await openColumnOptionsPopup(page, CHEM_COL);
    expect(popup.withinColumn,
      `the column-options icon that was clicked must sit inside the ${CHEM_COL} column's own header band, ` +
      'or the popup belongs to a different column').toBe(true);
    expect(popup.title, `the column popup must be the ${CHEM_COL} column's own`).toBe(CHEM_COL);

    const armedRoots = await chemFilterRoots(page);
    expect(armedRoots.outside,
      'the column popup must hold a substructure filter of its own, OUTSIDE the Filter Panel').toBe(1);
    expect(armedRoots.inPanel,
      'the Filter Panel must hold NO substructure card while the popup filter is the only one — otherwise the two ' +
      'readings below could come from a single filter').toBe(0);
    expect((await substructureStates(page)).length,
      `the panel's filter group must carry no ${CHEM_COL} substructure state while the popup filter is armed`).toBe(0);

    const popupInput = page.locator(`${POPUP_HOST} .grok-sketcher-input input`).first();
    await popupInput.fill(CHEM_PROBE);
    await popupInput.press('Enter');
    const popupCount = await sustainedCount(page, 'the popup filter\'s own row count');
    expect(popupCount,
      'the popup substructure filter must actually narrow — a count equal to the full row count would make the ' +
      'agreement below true of two filters that both do nothing').toBeLessThan(CHEM_FULL);
    expect(popupCount, 'the popup substructure filter must keep at least one row').toBeGreaterThan(0);

    const addFilterClicked = await page.evaluate(() => {
      const link = [...document.querySelectorAll('.d4-popup-host label.d4-link-action')]
        .find((e) => (e.textContent ?? '').trim().toLowerCase() === 'add filter'
          && e.getBoundingClientRect().width > 0) as HTMLElement | undefined;
      if (!link) return false;
      link.click();
      return true;
    });
    expect(addFilterClicked, 'the column popup must offer its own "Add filter" link').toBe(true);
    await page.waitForFunction(({col, type}) =>
      (grok.shell.tv.getFiltersGroup().getStates(col, type) || []).length > 0,
    {col: CHEM_COL, type: SUBSTRUCTURE_FILTER}, {timeout: 30_000, polling: 200});

    const transferred = await chemFilterRoots(page);
    expect(transferred.outside,
      'GROK-14952: the popup filter must be torn down when its criterion moves to the panel — a surviving popup ' +
      'filter keeps contributing its own mask and the panel card can no longer switch the filtering off').toBe(0);
    expect(transferred.inPanel, 'the panel must now hold exactly one substructure card').toBe(1);
    expect(await cardCount(page), 'the transfer must add exactly one card to the panel').toBe(cardsBefore + 1);
    const panelStates = await substructureStates(page);
    expect(panelStates.length, `the panel must carry exactly ONE ${CHEM_COL} substructure state after the transfer`).toBe(1);
    expect(hasStructure(panelStates[0].molBlock),
      'the panel card must carry the structure that was configured on the popup surface').toBe(true);

    const panelCount = await sustainedCount(page, 'the panel card\'s own row count');
    expect(panelCount,
      'the panel substructure card must actually narrow, so the agreement is between two real row sets').toBeLessThan(CHEM_FULL);
    expect(panelCount, 'the panel substructure card must keep at least one row').toBeGreaterThan(0);
    expect(panelCount,
      `GROK-14952: the two surfaces must not diverge — the panel card must pass exactly the row set the popup ` +
      `filter passed for the same structure (popup ${popupCount}, panel ${panelCount})`).toBe(popupCount);

    expect(await page.locator(POPUP_HOST).count(),
      'the column popup must still be up before it is dismissed — a dismissal of nothing is inert by default and ' +
      'proves nothing about the row set').toBeGreaterThan(0);
    const away = await page.evaluate((sel) => {
      const r = document.querySelector(sel)!.getBoundingClientRect();
      const x = r.left > window.innerWidth - r.right ? Math.round(r.left / 2)
        : Math.round((r.right + window.innerWidth) / 2);
      return {x, y: Math.round(r.top + r.height / 2)};
    }, POPUP_HOST);
    await page.mouse.click(away.x, away.y);
    await page.locator(POPUP_HOST).first().waitFor({state: 'detached', timeout: 20_000});
    expect(await page.locator(POPUP_HOST).count(),
      `the neutral click at (${away.x}, ${away.y}) must really dismiss the column popup — on a viewport where the ` +
      'click lands on nothing the row-count reading below would hold because nothing happened at all').toBe(0);
    expect(await sustainedCount(page, 'the row set after the popup is dismissed'),
      'dismissing the column popup must be inert — it must not move the row set the panel card now owns').toBe(panelCount);

    expect(await setCardEnabled(page, false),
      'the panel substructure card must switch off through its own enable checkbox and read disabled').toBe(true);
    expect(await sustainedCount(page, 'the row set with the panel card switched off'),
      'GROK-14952: switching the panel card off must lift ALL molecular filtering on that column — a row count still ' +
      'below the full one means a second, popup-side filter is still contributing its mask').toBe(CHEM_FULL);

    expect(await setCardEnabled(page, true),
      'the panel substructure card must switch back on through its own enable checkbox').toBe(true);
    expect(await sustainedCount(page, 'the row set with the panel card switched back on'),
      `re-arming the panel card must restore exactly the row set it held before (${panelCount})`).toBe(panelCount);
    console.log(`GROK-14952: popup-only ${popupCount}, panel-only ${panelCount}, disabled ${CHEM_FULL}`);
  });

  await softStep('Scenario 1, Step 11: dragging the molecular column header onto the panel adds an EMPTY card at the TOP',
    async () => {
      await removeAllViaHamburger(page);
      await closeFilterPanel(page);
      await reopenFilterPanel(page);
      expect(await removeSubstructureCards(page, molColumns),
        'the rebuilt panel must carry substructure cards to take down — with none present the drag below would not ' +
        'be adding a card that was absent').toBeGreaterThan(0);
      expect((await substructureStates(page)).length,
        `${CHEM_COL} must hold no substructure state before the drag, or the drop would UPDATE the existing card ` +
        'and add nothing, leaving both claims of this step unobservable').toBe(0);

      const anchorPick = await pickIntersectingCategory(page);
      const anchorColumn = anchorPick?.column
        ?? (await page.evaluate((molCol) => {
          const c = [...grok.shell.tv.dataFrame.columns].find((x: any) => x.name !== molCol && x.type === 'string');
          return c ? c.name : null;
        }, CHEM_COL));
      expect(anchorColumn, 'spgi-100 exposes no non-molecular categorical column to seat an anchor card on, so ' +
        '"the new card landed at the TOP" could not be distinguished from "it is the only card"').toBeTruthy();
      if (!(await fp.cardCaptions(page)).includes(anchorColumn!))
        await fp.addCardViaColumnSelector(page, anchorColumn!);

      const before = await panelCardOrder(page);
      expect(before.length, 'the panel must already carry at least one card before the drag — on an empty panel the ' +
        'dropped card is first by default and the "added at the top" claim is satisfied by arithmetic, not by the ' +
        'product placing it there').toBeGreaterThan(0);
      const countBefore = await sustainedCount(page, 'the row set before the header drag');

      await dragColumnHeaderToPanel(page, CHEM_COL);
      await waitForCardCount(page, before.length + 1);

      const after = await panelCardOrder(page);
      expect(after.length, `the drop must add exactly one card (panel held ${before.length})`).toBe(before.length + 1);
      expect(after[0].column,
        `the dropped ${CHEM_COL} filter must be placed at the TOP of the panel — it came back at index ` +
        `${after.findIndex((c) => c.column === CHEM_COL && c.source === SUBSTRUCTURE_SOURCE)} of ` +
        `${after.length} (order: ${JSON.stringify(after.map((c) => c.column))})`).toBe(CHEM_COL);
      expect(after[0].source,
        'the card at the top of the panel must be the Chem substructure filter the molecular column drop creates, ' +
        `not some other filter type on the same column (read: "${after[0].source}")`).toBe(SUBSTRUCTURE_SOURCE);

      const dropped = await substructureStates(page);
      expect(dropped.length,
        `${CHEM_COL} must carry exactly one substructure state after the drop`).toBe(1);
      expect(hasStructure(dropped[0].molBlock),
        'the dropped card must arrive EMPTY — a card that comes up already carrying a structure has filtered the ' +
        'table before the user drew anything, which is the regression this step exists for (state read: ' +
        `${JSON.stringify(dropped[0].molBlock)?.slice(0, 120)})`).toBe(false);
      expect(await sustainedCount(page, 'the row set after the header drag'),
        'an empty dropped card must not narrow the table — a row count below the pre-drag one means the card ' +
        'arrived with a criterion already applied').toBe(countBefore);
      console.log(`Step 11: panel ${before.length} -> ${after.length} cards, top = ${after[0].column}/${after[0].source}`);
    });

  await softStep('Scenario 1, Step 12: the substructure filter and a non-molecular filter intersect, and each toggles independently',
    async () => {
      // Step 11 leaves exactly one substructure card and it is empty, with no sketcher open —
      // it clears every card first, then drops a fresh one. commitSmiles only fills a dialog
      // that is already up, so without this the fill has no target.
      await openSketcherDialog(page);
      await commitSmiles(page, CHEM_PROBE);
      const substructureOnly = await sustainedCount(page, 'the row set with only the substructure filter armed');
      expect(substructureOnly, 'the substructure filter must narrow the table before anything is crossed with it')
        .toBeLessThan(CHEM_FULL);
      expect(substructureOnly, 'the substructure filter must keep at least one row, or the intersection below is ' +
        'empty for a reason that has nothing to do with crossing two filters').toBeGreaterThan(0);

      const pick = await pickIntersectingCategory(page);
      expect(pick, 'no non-molecular category on spgi-100 both splits the substructure hits and reaches rows the ' +
        'substructure filter rejects, so no proper intersection can be built on this dataset and the step would ' +
        'otherwise assert something the data cannot violate').not.toBeNull();
      const {column, category, categoryRows, intersection} = pick!;
      expect(intersection, `derived expectation is degenerate: intersection ${intersection} is not strictly below ` +
        `the substructure-only count ${substructureOnly}`).toBeLessThan(substructureOnly);
      expect(intersection, `derived expectation is degenerate: intersection ${intersection} is not strictly below ` +
        `the category-only count ${categoryRows}`).toBeLessThan(categoryRows);

      if (!(await fp.cardCaptions(page)).includes(column))
        await fp.addCardViaColumnSelector(page, column);
      const geomProbe = await cardOverlayGeometry(page, column);
      const pt = await selectCategoryRow(page, column, category);
      expect(pt, `no rendered row of the ${column} card selects "${category}" when clicked, so the second criterion ` +
        `would never be applied through the panel (card overlay top ${Math.round(geomProbe?.top ?? -1)}, ` +
        `height ${Math.round(geomProbe?.height ?? -1)}, derived rows ${geomProbe?.rows ?? -1}, ` +
        `viewport ${geomProbe?.viewportH ?? -1}, on screen: ${geomProbe?.onScreen ?? 'unknown'})`).not.toBeNull();
      expect(await selectedCategories(page, column),
        `the ${column} card must carry exactly the clicked category as its criterion`).toEqual([category]);

      const both = await sustainedCount(page, `the row set with the substructure filter and ${column}="${category}"`);
      expect(both, `both filters armed must pass exactly the rows satisfying BOTH — the substructure hits that also ` +
        `carry ${column}="${category}" number ${intersection} in spgi-100`).toBe(intersection);
      expect(both, 'the intersection must be strictly smaller than the substructure filter alone, otherwise the ' +
        'second criterion changed nothing and the two are not intersecting').toBeLessThan(substructureOnly);
      expect(both, 'the intersection must be strictly smaller than the category filter alone, otherwise the ' +
        'substructure criterion changed nothing').toBeLessThan(categoryRows);

      expect(await setCardEnabled(page, false),
        'the substructure card must switch off through its own enable checkbox').toBe(true);
      expect(await sustainedCount(page, `the row set with only ${column}="${category}" armed`),
        `switching the substructure card off must leave the ${column} criterion filtering ALONE — the row count ` +
        `must rise to that category's own ${categoryRows} rows, and each filter must therefore act independently`)
        .toBe(categoryRows);

      expect(await setCardEnabled(page, true),
        'the substructure card must switch back on through its own enable checkbox').toBe(true);
      expect(await sustainedCount(page, 'the row set with both filters armed again'),
        're-arming the substructure card must restore exactly the intersection it held before, so the toggle is ' +
        'reversible and neither criterion was destroyed by the other').toBe(intersection);
      console.log(`Step 12: substructure-only ${substructureOnly}, ${column}="${category}" only ${categoryRows}, ` +
        `intersection ${intersection}`);
    });

  await softStep('Scenario 1, Step 13: a DISABLED substructure filter survives collapsing and reopening the panel, still disabled',
    async () => {
      const armed = await sustainedCount(page, 'the row set with the substructure filter armed');
      expect(armed, 'the substructure filter must be filtering before it is switched off, or "the row count went ' +
        'back up" below would be true without the disable doing anything').toBeLessThan(CHEM_FULL);
      const structureBefore = (await substructureStates(page))[0]?.molBlock ?? null;
      expect(hasStructure(structureBefore),
        'the card must carry a drawn structure before the panel is collapsed — the point of the step is that the ' +
        'criterion survives along with the OFF state').toBe(true);

      expect(await setCardEnabled(page, false),
        'the substructure card must switch off through its own enable checkbox').toBe(true);
      const offCount = await sustainedCount(page, 'the row set with the substructure card switched off');

      await closeFilterPanel(page);
      await reopenFilterPanel(page);

      const restored = await readCardEnabled(page);
      expect(restored.found,
        `the ${CHEM_COL} substructure card must come back after the panel is reopened`).toBe(true);
      expect(restored.disabled,
        'the reopened card must still read DISABLED — a card that comes back enabled has silently re-applied a ' +
        'criterion the user switched off, and mere presence of the card does not detect that').toBe(true);
      expect(restored.checked,
        'the reopened card\'s own enable checkbox must still read unchecked, so the OFF state survived on the ' +
        'control the user actually toggles and not only as a CSS class').toBe(false);
      expect(hasStructure((await substructureStates(page))[0]?.molBlock ?? null),
        'the reopened card must still carry the structure that was drawn before the panel was collapsed — a card ' +
        'restored empty has lost the criterion even though it looks right').toBe(true);
      expect(await sustainedCount(page, 'the row set after the panel is reopened'),
        'a still-disabled card must still not be filtering after the reopen — a row count back below the ' +
        'switched-off one means the criterion was re-applied by the reopen').toBe(offCount);

      expect(await setCardEnabled(page, true),
        'the restored card must switch back on through its own enable checkbox').toBe(true);
      expect(await sustainedCount(page, 'the row set with the restored card switched back on'),
        `re-arming the restored card must reproduce exactly the row set it filtered to before the panel was ` +
        `collapsed (${armed}) — anything else means the surviving criterion is not the one that was saved`)
        .toBe(armed);
      console.log(`Step 13: armed ${armed}, disabled ${offCount}, restored-and-rearmed ${armed}`);
    });

  await softStep('Close chem table', async () => {
    await v.cleanupShell(page);
  });

  v.finishSpec('Chem scenario step failures');
});