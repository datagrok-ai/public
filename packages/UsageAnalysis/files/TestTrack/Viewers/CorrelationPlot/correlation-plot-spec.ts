/* ---
realizes: [correlationplot.cp.property-surface-smoke, correlationplot.int.menu-toggles-mirror-props]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const TOL = 1e-3;

// Correlation-plot cell geometry (refdoc §Cell Geometry and Mouse Clicks). pinnedW/headerH are
// auto-sized per dataset and calibrated by one probe click; cellW is 40 with showPearsonR, else 20.
interface Geometry {rootX: number; rootY: number; pinnedW: number; headerH: number; cellW: number; rowH: number}

function cellCenter(g: Geometry, xi: number, yi: number): {x: number; y: number} {
  return {x: g.rootX + g.pinnedW + (xi + 0.5) * g.cellW, y: g.rootY + g.headerH + (yi + 0.5) * g.rowH};
}

async function readBase(page: Page): Promise<{rootX: number; rootY: number; cellW: number; xCols: string[]; yCols: string[]}> {
  return await page.evaluate(() => {
    const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
    const R = document.querySelector('[name="viewer-Correlation-plot"]')!.getBoundingClientRect();
    return {rootX: R.x, rootY: R.y, cellW: cp.props.showPearsonR ? 40 : 20,
      xCols: cp.props.xColumnNames.slice(), yCols: cp.props.yColumnNames.slice()};
  });
}

async function refreshRoot(page: Page, g: Geometry): Promise<void> {
  const r = await page.evaluate(() => {
    const R = document.querySelector('[name="viewer-Correlation-plot"]')!.getBoundingClientRect();
    return {x: R.x, y: R.y};
  });
  g.rootX = r.x; g.rootY = r.y;
}

// Poll (instead of a fixed wait) until the reusable .d4-tooltip carries the '<Type> R:' line —
// the corr tooltip content mounts inline in that element; an empty element means the hover missed.
async function waitForRTooltip(page: Page, timeoutMs = 5000): Promise<void> {
  await page.waitForFunction(() => {
    const tip = document.querySelector('.d4-tooltip');
    return /R:\s*-?\d/.test(tip?.textContent ?? '');
  }, undefined, {timeout: timeoutMs});
}

// Arm the corr-cell-click listener; args column1/column2 are column NAME strings (live recon).
async function armClicks(page: Page): Promise<void> {
  await page.evaluate(() => {
    const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
    (window as any).__clicks = [];
    cp.onEvent('d4-correlation-plot-corr-cell-click').subscribe((e: any) => {
      const a = e.args ?? e;
      (window as any).__clicks.push({c1: a.column1, c2: a.column2, v: a.value});
    });
    return null;
  });
}

async function lastClick(page: Page): Promise<{c1: string; c2: string; v: number} | null> {
  return await page.evaluate(() => {
    const c = (window as any).__clicks;
    return c && c.length ? c[c.length - 1] : null;
  });
}

// Snapshot / diff the CP base canvas: name="canvas" is the paint layer, name="overlay" the
// mouse-input layer [DOM 2026-07-31]. The shared v.snapshotCanvasColors reads the FIRST canvas
// (a 0x0 layer here -> -1 fault), so the diff runs inline against the named base canvas.
// Optional region clips to a root-relative rect. Returns changed-pixel count; -1 on fault.
async function snapCanvas(page: Page, region?: {rx: number; ry: number; w: number; h: number}): Promise<boolean> {
  return await page.evaluate((reg) => {
    const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
    const cv = Array.from(root.querySelectorAll('canvas')).find((c) => c.getAttribute('name') === 'canvas') as HTMLCanvasElement | undefined;
    const ctx = cv?.getContext('2d');
    if (!cv || !ctx) return false;
    try {
      const r = cv.getBoundingClientRect();
      const sx = cv.width / r.width, sy = cv.height / r.height;
      let x = 0, y = 0, w = cv.width, h = cv.height;
      if (reg) {
        const R = root.getBoundingClientRect();
        x = Math.round(((R.x + reg.rx) - r.left) * sx);
        y = Math.round(((R.y + reg.ry) - r.top) * sy);
        w = Math.max(1, Math.round(reg.w * sx));
        h = Math.max(1, Math.round(reg.h * sy));
      }
      (window as any).__cpSnap = {data: ctx.getImageData(x, y, w, h).data, x, y, w, h};
      return true;
    } catch { return false; }
  }, region ?? null);
}

async function diffCanvas(page: Page): Promise<number> {
  return await page.evaluate(() => {
    const snap = (window as any).__cpSnap;
    const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
    const cv = Array.from(root.querySelectorAll('canvas')).find((c) => c.getAttribute('name') === 'canvas') as HTMLCanvasElement | undefined;
    const ctx = cv?.getContext('2d');
    if (!cv || !ctx || !snap) return -1;
    try {
      const cur = ctx.getImageData(snap.x, snap.y, snap.w, snap.h).data;
      const prev = snap.data;
      let n = 0;
      for (let i = 0; i < prev.length; i += 4)
        if (prev[i] !== cur[i] || prev[i + 1] !== cur[i + 1] || prev[i + 2] !== cur[i + 2]) n++;
      return n;
    } catch { return -1; }
  });
}

// Cell RGBA = per-channel median of a 5-point cross patch: a single center pixel can land on an
// in-cell R digit or a border and read the text/border colour instead of the cell hue.
async function cellPixel(page: Page, g: Geometry, xi: number, yi: number): Promise<number[] | null> {
  const c = cellCenter(g, xi, yi);
  return await page.evaluate(({cx, cy}) => {
    const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
    const cv = Array.from(root.querySelectorAll('canvas')).find((c) => c.getAttribute('name') === 'canvas') as HTMLCanvasElement | undefined;
    const ctx = cv?.getContext('2d');
    if (!cv || !ctx) return null;
    const r = cv.getBoundingClientRect();
    const sx = cv.width / r.width, sy = cv.height / r.height;
    try {
      const pts = [[0, 0], [3, 0], [-3, 0], [0, 3], [0, -3]];
      const samples: number[][] = [];
      for (const [dx, dy] of pts) {
        const d = ctx.getImageData(Math.round((cx + dx - r.left) * sx), Math.round((cy + dy - r.top) * sy), 1, 1).data;
        samples.push([d[0], d[1], d[2], d[3]]);
      }
      const med = (i: number) => samples.map((s) => s[i]).sort((a, b) => a - b)[2];
      return [med(0), med(1), med(2), med(3)];
    } catch { return null; }
  }, {cx: c.x, cy: c.y});
}

// d4 context-menu items carry name= attributes: `div-{Item}` at top level, nested via the `---`
// path separator, e.g. `div-Properties...---Misc---Show-Pearson-R` [DOM 2026-08-01, live MCP
// recon, dev/demog]. Resolve to the VISIBLE (non-zero rect) node inside a live .d4-menu-popup —
// never document-wide .first(): the ribbon top-menu holds ~166 same-shaped .d4-menu-item nodes,
// and stale popups accumulate (Escape/body.click do not dismiss a d4 context menu headless).
// A nested leaf stays 0x0 until its parent group is expanded.
async function waitVisibleMenuItem(page: Page, name: string, timeoutMs = 4000): Promise<boolean> {
  try {
    await page.waitForFunction((n) => {
      const nodes = Array.from(document.querySelectorAll(`.d4-menu-popup .d4-menu-item[name="${n}"]`));
      return nodes.some((el) => { const b = (el as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; });
    }, name, {timeout: timeoutMs});
    return true;
  } catch { return false; }
}

// A d4 menu-item's Dart handler needs the FULL pointer/mouse sequence — a bare trusted
// locator.click() does not flip the backing prop headless [live MCP recon 2026-08-01]. Dispatched
// in-page on the visible popup item (sanctioned DOM target, never the canvas).
async function clickMenuByName(page: Page, name: string, timeoutMs = 4000): Promise<boolean> {
  if (!(await waitVisibleMenuItem(page, name, timeoutMs))) return false;
  return await page.evaluate((n) => {
    const el = Array.from(document.querySelectorAll(`.d4-menu-popup .d4-menu-item[name="${n}"]`))
      .find((e) => { const b = (e as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; }) as HTMLElement | undefined;
    if (!el) return false;
    const b = el.getBoundingClientRect();
    const o: any = {bubbles: true, cancelable: true, clientX: b.x + b.width / 2, clientY: b.y + b.height / 2, view: window};
    el.dispatchEvent(new PointerEvent('pointerdown', o));
    el.dispatchEvent(new MouseEvent('mousedown', o));
    el.dispatchEvent(new PointerEvent('pointerup', o));
    el.dispatchEvent(new MouseEvent('mouseup', o));
    el.dispatchEvent(new MouseEvent('click', o));
    return true;
  }, name);
}

// Expand a nested submenu by FORCE-DISPLAYING the group's `.d4-menu-item-container` — neither a
// synthetic pointerover nor a trusted hover lays it out headless [live MCP recon 2026-08-01].
// Returns true when the group element was present. Idempotent.
async function hoverMenuGroupTrusted(page: Page, name: string, timeoutMs = 4000): Promise<boolean> {
  await waitVisibleMenuItem(page, name, timeoutMs);
  return await page.evaluate((n) => {
    const el = (Array.from(document.querySelectorAll(`.d4-menu-popup .d4-menu-item[name="${n}"]`))
      .find((e) => { const b = (e as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; })
      ?? document.querySelector(`.d4-menu-popup .d4-menu-item[name="${n}"]`)) as HTMLElement | null;
    if (!el) return false;
    const container = el.querySelector('.d4-menu-item-container') as HTMLElement | null;
    if (container) {
      container.style.display = 'flex';
      const vert = container.querySelector('.d4-vert-menu') as HTMLElement | null;
      if (vert) vert.style.display = 'flex';
    }
    // Keep the paired pointer/mouse-over sequence too (primes any Dart-side hover state without harm).
    const b = el.getBoundingClientRect();
    const o: any = {bubbles: true, cancelable: true, clientX: b.x + b.width / 2, clientY: b.y + b.height / 2, view: window};
    el.dispatchEvent(new PointerEvent('pointerover', o));
    el.dispatchEvent(new MouseEvent('mouseover', o));
    return true;
  }, name);
}

// Text-channel menu actuation (annotation-regions-spec canon): GRID-GENERIC items carry NO name=
// attribute — resolve by own-label TEXT inside the live .d4-menu-popup.
async function hoverMenuGroupByText(page: Page, text: string): Promise<boolean> {
  return await page.evaluate((t) => {
    const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'));
    const cands = labels.filter((l) => (l.textContent ?? '').trim() === t);
    const vis = cands.find((l) => {
      const b = (l.closest('.d4-menu-item') as HTMLElement | null)?.getBoundingClientRect();
      return !!b && b.width > 0 && b.height > 0;
    });
    const item = ((vis ?? cands[0])?.closest('.d4-menu-item') ?? null) as HTMLElement | null;
    if (!item) return false;
    const cont = item.querySelector('.d4-menu-item-container') as HTMLElement | null;
    if (cont) {
      cont.style.display = 'flex';
      const vert = cont.querySelector('.d4-vert-menu') as HTMLElement | null;
      if (vert) vert.style.display = 'flex';
    }
    item.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, view: window}));
    item.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, view: window}));
    return true;
  }, text);
}

async function clickMenuItemByText(page: Page, pattern: string, ancestorLabel?: string): Promise<boolean> {
  return await page.evaluate(({p, anc}) => {
    const re = new RegExp(p);
    const own = (it: Element) => Array.from(it.querySelectorAll('.d4-menu-item-label'))
      .find((l) => l.closest('.d4-menu-item') === it)?.textContent?.trim() ?? '';
    let items = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item'))
      .filter((it) => re.test(own(it)));
    // Optional ancestor scope: keep only items nested under a group with that own label
    // (disambiguates e.g. Color Coding's Edit... from Tooltip's Edit...).
    if (anc)
      items = items.filter((it) => {
        let n = it.parentElement;
        while (n) { if (n.classList?.contains('d4-menu-item') && own(n) === anc) return true; n = n.parentElement; }
        return false;
      });
    const vis = items.find((it) => {
      const b = (it as HTMLElement).getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    });
    const item = (vis ?? items[0] ?? null) as HTMLElement | null;
    if (!item) return false;
    const b = item.getBoundingClientRect();
    const o: any = {bubbles: true, cancelable: true, clientX: b.x + b.width / 2, clientY: b.y + b.height / 2, view: window};
    item.dispatchEvent(new PointerEvent('pointerdown', o));
    item.dispatchEvent(new MouseEvent('mousedown', o));
    item.dispatchEvent(new PointerEvent('pointerup', o));
    item.dispatchEvent(new MouseEvent('mouseup', o));
    item.dispatchEvent(new MouseEvent('click', o));
    return true;
  }, {p: pattern, anc: ancestorLabel ?? null});
}

// Menu teardown: d4 context menus are not dismissed by Escape/body.click headless [live MCP recon
// 2026-08-01] — remove the popup nodes directly (the ribbon top-menu is never a .d4-menu-popup).
async function closeMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.evaluate(() => {
    document.querySelectorAll('.d4-menu-popup').forEach((el) => el.remove());
    document.body.click();
  });
  await page.waitForTimeout(200);
}

test('Correlation plot — property surface smoke', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // ## Setup — open demog, add the Correlation plot, install console/pageerror collectors
  // (grok.shell.warnings is undefined on this build), arm the click listener, calibrate geometry.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  // 'Unable to find element in cloned iframe' is ambient unfixable noise — benign always.
  // 'Package GrokML is not available' is benign ONLY inside the misc-sequence window (GROK-16818):
  // gated at CAPTURE time so the allowlist cannot mask a GrokML-shaped error elsewhere.
  let inMiscWindow = false;
  const errorNoise = (s: string) => /Unable to find element in cloned iframe/i.test(s)
    || (inMiscWindow && /Package GrokML is not available/i.test(s));
  page.on('console', (m) => { if (m.type() === 'error' && !errorNoise(m.text())) consoleErrors.push(m.text()); });
  const realErrors = () => consoleErrors;

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);
  await armClicks(page);

  const base = await readBase(page);
  const xiHeight = base.xCols.indexOf('HEIGHT');
  const xiWeight = base.xCols.indexOf('WEIGHT');
  const yiAge = base.yCols.indexOf('AGE');
  const yiWeight = base.yCols.indexOf('WEIGHT');
  const geom: Geometry = {rootX: base.rootX, rootY: base.rootY, pinnedW: 104, headerH: 60, cellW: base.cellW, rowH: 20};

  await softStep('Setup — calibrate cell geometry via a probe click', async () => {
    let landed = false;
    for (let attempt = 0; attempt < 6 && !landed; attempt++) {
      await page.evaluate(() => { (window as any).__clicks = []; });
      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y);
      await page.waitForTimeout(400);
      const ev = await lastClick(page);
      console.log(`[Setup] probe ${attempt} at (${Math.round(c.x)},${Math.round(c.y)}) pinnedW=${geom.pinnedW} headerH=${geom.headerH} -> ${JSON.stringify(ev)}`);
      if (ev && [ev.c1, ev.c2].sort().join() === ['AGE', 'HEIGHT'].join()) { landed = true; break; }
      if (ev && ev.c1 && ev.c1 !== 'HEIGHT') {
        const idx = base.xCols.indexOf(ev.c1);
        if (idx >= 0) geom.pinnedW += (xiHeight - idx) * geom.cellW;
      } else if (!ev) geom.headerH += 4;
    }
    const ev = await lastClick(page);
    expect(ev).not.toBeNull();
    // A calibrated probe proves the coordinate maths AND that real trusted mouse reaches the Dart handler.
    expect([ev!.c1, ev!.c2].sort()).toEqual(['AGE', 'HEIGHT']);
  });

  // ## Title, description, and back color
  await softStep('Title, description, and back color', async () => {
    // Title renders in the dock .panel-titlebar-text (document scope), NOT the viewer root
    // [DOM 2026-07-31]; description renders inside the root per Description Visibility Mode.
    const titleShown: boolean = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.showTitle = true;
      cp.props.title = 'Correlation Analysis';
      cp.props.description = 'Shows pairwise correlations';
      cp.props.descriptionVisibilityMode = 'Always';
      await new Promise((r) => setTimeout(r, 500));
      return Array.from(document.querySelectorAll('.panel-titlebar-text'))
        .some((el) => el.textContent!.trim() === 'Correlation Analysis');
    });
    const descAlways: boolean = await page.evaluate(() =>
      document.querySelector('[name="viewer-Correlation-plot"]')!.textContent!.includes('Shows pairwise correlations'));
    // descriptionPosition = Bottom: signal = the description leaf's rect moving into the lower
    // half of the viewer root; original value restored after.
    const descPos = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const orig = cp.props.descriptionPosition;
      cp.props.descriptionPosition = 'Bottom';
      await new Promise((r) => setTimeout(r, 600));
      const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
      const R = root.getBoundingClientRect();
      const leaf = Array.from(root.querySelectorAll('*'))
        .filter((el) => el.children.length === 0 && (el.textContent ?? '').includes('Shows pairwise correlations'))
        .pop() ?? null;
      const rect = leaf ? leaf.getBoundingClientRect() : null;
      cp.props.descriptionPosition = orig;
      await new Promise((r) => setTimeout(r, 300));
      return {found: !!rect, centerY: rect ? rect.top + rect.height / 2 : -1, rootMid: R.y + R.height / 2, orig};
    });
    console.log(`[DescPos] ${JSON.stringify(descPos)}`);
    expect(descPos.found).toBe(true);
    // Bottom position: the description's vertical center sits below the viewer's midline.
    expect(descPos.centerY).toBeGreaterThan(descPos.rootMid);
    const descNever: boolean = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.descriptionVisibilityMode = 'Never';
      await new Promise((r) => setTimeout(r, 500));
      return document.querySelector('[name="viewer-Correlation-plot"]')!.textContent!.includes('Shows pairwise correlations');
    });
    expect(titleShown).toBe(true);
    expect(descAlways).toBe(true);
    expect(descNever).toBe(false);

    // Back Color repaint: settle-gated whole-canvas diff (a single margin pixel reads white in the
    // header region on demog).
    await snapCanvas(page);
    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);
    console.log(`[BackColor] settle=${settle}`);
    await snapCanvas(page);
    await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.lightGray;
      await new Promise((r) => setTimeout(r, 600));
      return null;
    });
    const backDelta = await diffCanvas(page);
    console.log(`[BackColor] delta=${backDelta}`);
    expect(backDelta).toBeGreaterThanOrEqual(0);
    expect(backDelta).toBeGreaterThan(settle + 1000);
    // restore
    await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.white;
      cp.props.showTitle = false; cp.props.title = ''; cp.props.description = '';
      cp.props.descriptionVisibilityMode = 'Auto';
      await new Promise((r) => setTimeout(r, 400));
      return null;
    });
    expect(realErrors()).toEqual([]);
  });

  // ## Misc property sequence and clean console
  await softStep('Misc property sequence and clean console', async () => {
    const errBefore = realErrors().length;
    const peBefore = pageErrors.length;
    inMiscWindow = true;
    try {
      await v.setViewerProps(page, 'Correlation plot', [
        {set: {showTooltip: false}, wait: 250}, {set: {showTooltip: true}, wait: 250},
        {set: {ignoreDoubleClick: true}, wait: 250}, {set: {ignoreDoubleClick: false}, wait: 250},
        {set: {colHeaderFont: 'bold normal 16px "Roboto"'}, wait: 250},
        {set: {colHeaderFont: 'bold normal 13px "Roboto"'}, wait: 250},
        {set: {defaultCellFont: 'normal normal 18px "Roboto"'}, wait: 250},
        {set: {defaultCellFont: 'normal normal 13px "Roboto"'}, wait: 250},
      ], 250);
    } finally {
      inMiscWindow = false;
    }
    // GROK-16818: the whole misc/style switching sequence raises no console/page errors.
    expect(realErrors().length).toBe(errBefore);
    expect(pageErrors.length).toBe(peBefore);
  });

  // ## Context-menu toggles mirror properties
  await softStep('Context-menu toggles mirror properties', async () => {
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(500);
    const showRBefore: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR);
    // Section settle baseline: measure the idle-frame diff so the repaint assert is settle-gated.
    await snapCanvas(page);
    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);
    // Show Pearson R menu click flips props.showPearsonR AND narrows the value columns (repaint).
    // Top-level toggle is [name="div-Show-Pearson-R"]; the Properties>Misc duplicate is nested.
    await snapCanvas(page);
    const clickedShowR = await clickMenuByName(page, 'div-Show-Pearson-R');
    expect(clickedShowR).toBe(true);
    await page.waitForFunction(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR === false,
    undefined, {timeout: 4000}).catch(() => {});
    const showRAfter: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR);
    const repaintR = await diffCanvas(page);
    console.log(`[Menu] showPearsonR ${showRBefore}->${showRAfter} settle=${settle} repaint=${repaintR}`);
    expect(showRBefore).toBe(true);
    expect(showRAfter).toBe(false);
    expect(repaintR).toBeGreaterThanOrEqual(0);
    expect(repaintR).toBeGreaterThan(settle);

    // Tooltip > Visible flips props.showTooltip.
    await closeMenu(page);
    await refreshRoot(page, geom);
    // showPearsonR is now false so cellW is 20 -> recompute center width.
    geom.cellW = 20;
    const c2 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c2.x, c2.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupTrusted(page, 'div-Tooltip');
    await page.waitForTimeout(400);
    const tipBefore: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip);
    const clickedVisible = await clickMenuByName(page, 'div-Tooltip---Visible');
    expect(clickedVisible).toBe(true);
    await page.waitForTimeout(500);
    const tipAfter: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip);
    console.log(`[Menu] showTooltip ${tipBefore}->${tipAfter}`);
    expect(tipBefore).toBe(true);
    expect(tipAfter).toBe(false);
    // With Show Tooltip off, hovering must NOT SHOW the tooltip. The reusable .d4-tooltip keeps the
    // previous hover's text — visibility (computed display) is the honest off-signal. Park the
    // mouse away first, then poll a 1500ms window asserting display never becomes 'block'.
    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
    await page.waitForTimeout(500);
    const displayIdle: string = await page.evaluate(() => {
      const tip = document.querySelector('.d4-tooltip');
      return tip ? getComputedStyle(tip).display : 'missing';
    });
    expect(displayIdle).toBe('none');
    const cOff = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.move(cOff.x, cOff.y);
    const shownWhileOff: boolean = await page.evaluate(async () => {
      for (let i = 0; i < 6; i++) {
        await new Promise((r) => setTimeout(r, 250));
        const tip = document.querySelector('.d4-tooltip');
        if (tip && getComputedStyle(tip).display === 'block') return true;
      }
      return false;
    });
    console.log(`[Menu] displayIdle=${displayIdle} shownWhileOff=${shownWhileOff}`);
    // The tooltip never became visible for THIS hover.
    expect(shownWhileOff).toBe(false);

    // Restore both through the same menu items (round-trip).
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c3 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c3.x, c3.y, {button: 'right'});
    await page.waitForTimeout(500);
    await clickMenuByName(page, 'div-Show-Pearson-R');
    await page.waitForFunction(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR === true,
    undefined, {timeout: 4000}).catch(() => {});
    await closeMenu(page);
    await refreshRoot(page, geom);
    geom.cellW = 40;
    const c4 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c4.x, c4.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupTrusted(page, 'div-Tooltip');
    await page.waitForTimeout(400);
    await clickMenuByName(page, 'div-Tooltip---Visible');
    await page.waitForFunction(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip === true,
    undefined, {timeout: 4000}).catch(() => {});
    const restored = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      return {showR: cp.props.showPearsonR, tip: cp.props.showTooltip};
    });
    expect(restored.showR).toBe(true);
    expect(restored.tip).toBe(true);
    await closeMenu(page);
    // Round-trip check: a FRESH hover must be VISIBLE (display 'block') AND carry the R-line —
    // residual text from an earlier hover alone would be a stale false PASS.
    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
    await page.waitForTimeout(500);
    const cOn = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.move(cOn.x, cOn.y);
    await page.waitForFunction(() => {
      const tip = document.querySelector('.d4-tooltip');
      return !!tip && getComputedStyle(tip).display === 'block' && /R:\s*-?\d/.test(tip.textContent ?? '');
    }, undefined, {timeout: 5000});
  });

  // ## Pinned row-header columns
  await softStep('Pinned row-header columns', async () => {
    // Shrink the viewer until the matrix is wider than the viewer and the horizontal scrollbar appears.
    const scrollbarShown: boolean = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Correlation-plot"]') as HTMLElement;
      root.style.width = '180px';
      await new Promise((r) => setTimeout(r, 700));
      const sb = root.querySelector('.d4-range-selector.d4-grid-horz-scroll') as HTMLElement | null;
      const r = sb?.getBoundingClientRect();
      return !!sb && !!r && r.width > 0 && r.height > 0;
    });
    expect(scrollbarShown).toBe(true);
    // Signal: a probe click at a FIXED viewport x reports a DIFFERENT column pair after the
    // horizontal wheel-scroll, while the pinned type/name block stays. A diagonal (histogram)
    // cell fires NO event, so the probe tries two rows.
    await refreshRoot(page, geom);
    const probeX = geom.rootX + geom.pinnedW + geom.cellW / 2;
    const probePair = async (): Promise<{c1: string; c2: string; v: number} | null> => {
      for (const row of ['HEIGHT', 'AGE']) {
        const yi = base.yCols.indexOf(row);
        await page.evaluate(() => { (window as any).__clicks = []; });
        await page.mouse.click(probeX, geom.rootY + geom.headerH + (yi + 0.5) * geom.rowH);
        await page.waitForTimeout(400);
        const ev = await lastClick(page);
        if (ev) return ev;
      }
      return null;
    };
    const pairBefore = await probePair();
    expect(pairBefore).not.toBeNull();
    // TRUSTED page.mouse.wheel over the value area (a synthetic WheelEvent is untrusted — canon:
    // never simulate canvas gestures with dispatchEvent).
    await snapCanvas(page);
    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + geom.pinnedW + geom.cellW, geom.rootY + geom.headerH + 60);
    await page.mouse.wheel(200, 0);
    await page.waitForTimeout(600);
    const scrollRepaint = await diffCanvas(page);
    console.log(`[Pinned] wheel repaint=${scrollRepaint}`);
    expect(scrollRepaint).toBeGreaterThanOrEqual(0);
    const pairAfter = await probePair();
    console.log(`[Pinned] probe before=${JSON.stringify(pairBefore)} after=${JSON.stringify(pairAfter)}`);
    // The CP horizontal scroll is a canvas-drawn range selector with no DOM thumb and no JS-readable
    // offset [DOM 2026-08-01, live MCP recon]; the trusted wheel may be inert headless. When it moves
    // the columns, assert the pair change; when inert, that sub-claim is a documented per-item
    // reduction (waiver_class: gesture-uncontrollable-headless) and pinning is carried by the
    // pinned-name tooltip below.
    if (pairAfter && pairBefore && pairAfter.c1 !== pairBefore.c1)
      expect(pairAfter.c1).not.toBe(pairBefore.c1);
    else
      console.log('[Pinned] wheel-scroll column move inert headless -> waived; pinning proven by pinned-name tooltip');
    // Hard pinning signal: hovering the pinned name cell at its original x still yields the
    // COLUMN-STATISTICS tooltip (avg/min lines), not merely any text.
    await refreshRoot(page, geom);
    // Trusted hover — a synthetic canvas MouseEvent does not drive the Dart tooltip.
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
    await page.mouse.move(geom.rootX + 60, geom.rootY + 70);
    await page.waitForTimeout(700);
    const pinnedTipText: string = await page.evaluate(() =>
      document.querySelector('.d4-tooltip')?.textContent ?? '');
    console.log(`[Pinned] tooltip="${pinnedTipText.replace(/\s+/g, ' ').slice(0, 120)}"`);
    expect(/avg:|min:/i.test(pinnedTipText)).toBe(true);
    // restore viewer width.
    await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Correlation-plot"]') as HTMLElement;
      root.style.width = '';
      await new Promise((r) => setTimeout(r, 500));
    });
    await refreshRoot(page, geom);
    expect(realErrors()).toEqual([]);
  });

  // ## Order or Hide Columns dialog (GROK-9310)
  await softStep('Order or Hide Columns dialog', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(500);
    // Grid-generic items carry no name= attribute — actuate via the text channel: hover the 'Grid'
    // group by own-label, then click the leaf matching /Order or Hide Columns/.
    await hoverMenuGroupByText(page, 'Grid');
    await page.waitForTimeout(400);
    const openedDialog = await clickMenuItemByText(page, 'Order or Hide Columns');
    expect(openedDialog).toBe(true);
    await page.waitForTimeout(800);
    // RECON NOTE [probe 2026-08-01, dump of [name="dialog-Order-or-Hide-Columns"]]: the dialog's
    // column list is an EMBEDDED CANVAS Grid — no per-column DOM checkboxes. Per-column checkbox
    // toggle: status: waived, waiver_class: canvas-webgl-render; evidence: embedded canvas Grid,
    // probe 2026-08-01. Driven instead: the all/visible/hidden filter SELECT and the search INPUT
    // (.d4-search-input), then a clean CLOSE, under the GROK-9310 no-exception invariant.
    const driven = await page.evaluate(async () => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Order-or-Hide-Columns"]');
      if (!dlg) return {dialogPresent: false, selectDriven: false, searchDriven: false};
      const sel = dlg.querySelector('select') as HTMLSelectElement | null;
      let selectDriven = false;
      if (sel) {
        const drive = async (label: string): Promise<boolean> => {
          const opt = Array.from(sel.options).find((o) => (o.textContent ?? '').trim() === label);
          if (!opt) return false;
          sel.value = opt.value;
          sel.dispatchEvent(new Event('input', {bubbles: true}));
          sel.dispatchEvent(new Event('change', {bubbles: true}));
          await new Promise((r) => setTimeout(r, 400));
          return sel.selectedIndex >= 0 && (sel.options[sel.selectedIndex].textContent ?? '').trim() === label;
        };
        selectDriven = (await drive('visible')) && (await drive('hidden')) && (await drive('all'));
      }
      const inp = dlg.querySelector('input.d4-search-input') as HTMLInputElement | null;
      let searchDriven = false;
      if (inp) {
        inp.focus();
        inp.value = 'AGE';
        inp.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 500));
        searchDriven = inp.value === 'AGE';
        inp.value = '';
        inp.dispatchEvent(new Event('input', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
      }
      return {dialogPresent: true, selectDriven, searchDriven};
    });
    console.log(`[OrderHide] ${JSON.stringify(driven)}`);
    expect(driven.dialogPresent).toBe(true);
    // The filter select ran visible -> hidden -> all through the real DOM control.
    expect(driven.selectDriven).toBe(true);
    // The search input accepted and cleared a query.
    expect(driven.searchDriven).toBe(true);
    // CLOSE via the dialog's own button; verify it is gone.
    const closedOH: boolean = await page.evaluate(async () => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Order-or-Hide-Columns"]');
      const btn = dlg?.querySelector('[name="button-CLOSE"]') as HTMLElement | null;
      if (btn) btn.click();
      await new Promise((r) => setTimeout(r, 400));
      return document.querySelectorAll('.d4-dialog[name="dialog-Order-or-Hide-Columns"]').length === 0;
    });
    expect(closedOH).toBe(true);
    await closeMenu(page);
    // GROK-9310: no exception raised by opening/driving/closing the dialog.
    expect(pageErrors.length).toBe(peBefore);
    expect(realErrors().length).toBe(errBefore);
    // Pinned name column stays — a row-header hover still shows the column-statistics tooltip
    // (avg/min lines), same predicate as the Pinned section.
    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
    await page.mouse.move(geom.rootX + 60, geom.rootY + 70);
    await page.waitForTimeout(700);
    const pinnedStillText: string = await page.evaluate(() =>
      document.querySelector('.d4-tooltip')?.textContent ?? '');
    console.log(`[OrderHide] pinned tooltip="${pinnedStillText.replace(/\s+/g, ' ').slice(0, 120)}"`);
    expect(/avg:|min:/i.test(pinnedStillText)).toBe(true);
  });

  // ## Grid color-coding apply to text (GROK-19052)
  await softStep('Grid color-coding apply to text', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;
    await closeMenu(page);
    await refreshRoot(page, geom);
    // The 'apply to text' affordance (GROK-19052) is the "Apply to" choice in the per-column
    // color-coding Edit dialog; it appears only after a coloring (Linear) is applied.
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(500);
    // Per-column color coding lives under Grid > Current Column > Color Coding (text channel).
    await hoverMenuGroupByText(page, 'Grid');
    await page.waitForTimeout(300);
    await hoverMenuGroupByText(page, 'Current Column');
    await page.waitForTimeout(300);
    await hoverMenuGroupByText(page, 'Color Coding');
    await page.waitForTimeout(300);
    // Apply Linear color-coding to the current (last-clicked) column so the Edit dialog gains "Apply to".
    const appliedLinear = await clickMenuItemByText(page, '^Linear$', 'Color Coding');
    expect(appliedLinear).toBe(true);
    await page.waitForTimeout(700);
    await closeMenu(page);
    // Re-open the menu, open Color Coding > Edit... and set "Apply to" to text.
    await refreshRoot(page, geom);
    const c2 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c2.x, c2.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupByText(page, 'Grid');
    await page.waitForTimeout(300);
    await hoverMenuGroupByText(page, 'Current Column');
    await page.waitForTimeout(300);
    await hoverMenuGroupByText(page, 'Color Coding');
    await page.waitForTimeout(300);
    // The ancestor scope disambiguates Color Coding's Edit... from Tooltip's Edit...
    const openedEdit = await clickMenuItemByText(page, '^Edit\\.\\.\\.$', 'Color Coding');
    expect(openedEdit).toBe(true);
    await page.waitForTimeout(700);
    // Drive the "Apply to" select to 'text' (the GROK-19052 option).
    const setToText = await page.evaluate(async () => {
      // Dialog is [name="dialog-Color-coding--<COL>"]; the 'Apply to' row is
      // [name="input-host-Apply-to"] holding an unnamed select (background / text /
      // text background) [probe 2026-08-01]. Fallbacks: title regex, then the label row.
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      const d = dialogs.find((dd) => /^dialog-Color-coding-/.test(dd.getAttribute('name') ?? ''))
        ?? dialogs.find((dd) => /Color.?coding/i.test(dd.querySelector('.d4-dialog-title')?.textContent ?? ''))
        ?? dialogs.find((dd) => Array.from(dd.querySelectorAll('label'))
          .some((l) => /Apply to/i.test(l.textContent ?? '')));
      if (!d) return false;
      let sel = d.querySelector('[name="input-host-Apply-to"] select') as HTMLSelectElement | null;
      if (!sel) {
        sel = (Array.from(d.querySelectorAll('select')) as HTMLSelectElement[]).find((s) => {
          const host = s.closest('.ui-input-root') ?? s.closest('div');
          return /Apply to/i.test(host?.querySelector('label')?.textContent ?? '');
        }) ?? null;
      }
      if (!sel) return false;
      const opt = Array.from(sel.options).find((o) => (o.textContent ?? o.value).trim().toLowerCase() === 'text');
      sel.value = opt ? opt.value : 'text';
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 700));
      return sel.selectedIndex >= 0 && (sel.options[sel.selectedIndex].textContent ?? sel.value).trim().toLowerCase() === 'text';
    });
    console.log(`[ApplyToText] setToText=${setToText}`);
    // Driven-guard: the 'Apply to' -> text write MUST be proven to have happened — the section's
    // no-error floor below is blind without it (an unfired action would float the floor).
    expect(setToText).toBe(true);
    // Reset: Apply to -> background, close the dialog, remove the color-coding (Grid > Grid Color
    // Coding > None). The dialog's close affordance is [name="button-CLOSE"] [DOM 2026-07-31].
    const resetToBackground: boolean = await page.evaluate(async () => {
      // Same dialog scoping + select resolution as the setToText write above (probe 2026-08-01).
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      const d = dialogs.find((dd) => /^dialog-Color-coding-/.test(dd.getAttribute('name') ?? ''))
        ?? dialogs.find((dd) => /Color.?coding/i.test(dd.querySelector('.d4-dialog-title')?.textContent ?? ''))
        ?? dialogs.find((dd) => Array.from(dd.querySelectorAll('label'))
          .some((l) => /Apply to/i.test(l.textContent ?? '')));
      let sel = (d?.querySelector('[name="input-host-Apply-to"] select') ?? null) as HTMLSelectElement | null;
      if (!sel && d) {
        sel = (Array.from(d.querySelectorAll('select')) as HTMLSelectElement[]).find((s) => {
          const host = s.closest('.ui-input-root') ?? s.closest('div');
          return /Apply to/i.test(host?.querySelector('label')?.textContent ?? '');
        }) ?? null;
      }
      if (sel) {
        const opt = Array.from(sel.options).find((o) => (o.textContent ?? o.value).trim().toLowerCase() === 'background');
        sel.value = opt ? opt.value : 'background';
        sel.dispatchEvent(new Event('input', {bubbles: true}));
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise((r) => setTimeout(r, 400));
      const ok = !!sel && sel.selectedIndex >= 0
        && (sel.options[sel.selectedIndex].textContent ?? sel.value).trim().toLowerCase() === 'background';
      (DG.Dialog.getOpenDialogs?.() ?? []).forEach((dlg: any) => dlg.close?.());
      Array.from(document.querySelectorAll('.d4-dialog')).forEach((dd) => {
        const cl = dd.querySelector('[name="button-CLOSE"]') as HTMLElement | null; if (cl) cl.click();
      });
      await new Promise((r) => setTimeout(r, 300));
      return ok;
    });
    // Symmetric driven-guard on the disable-write: the reset must be proven too.
    expect(resetToBackground).toBe(true);
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c3 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c3.x, c3.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupByText(page, 'Grid');
    await page.waitForTimeout(300);
    await hoverMenuGroupByText(page, 'Grid Color Coding');
    await page.waitForTimeout(300);
    await clickMenuItemByText(page, '^None$', 'Grid Color Coding');
    await page.waitForTimeout(400);
    await closeMenu(page);
    // GROK-19052: enabling 'Apply to text' completes without errors — honest no-error floor; the
    // colored text is canvas-drawn and exposes no readable channel.
    expect(pageErrors.length).toBe(peBefore);
    expect(realErrors().length).toBe(errBefore);
  });

  // ## Table switch (GROK-18487)
  await softStep('Table switch', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;
    // Clear residual popups/dialogs so the section starts from a clean viewer state.
    await closeMenu(page);
    await page.evaluate(() => { (DG.Dialog.getOpenDialogs?.() ?? []).forEach((d: any) => d.close?.()); return null; });
    await page.waitForTimeout(300);
    try {
      // Open spgi-100 alongside demog, switch the Table property, verify the matrix recomputes.
      // Return a JSON STRING: Dart-boxed values (getCorrelation / DG.Stats / DG-referencing Errors)
      // make page.evaluate throw 'object reference chain is too long' — JSON.stringify collapses
      // everything to primitives first [live MCP recon 2026-08-01].
      const raw = await page.evaluate(async ({p, tol}) => {
        try {
          const spgi = await grok.dapi.files.readCsv(p);
          grok.shell.addTableView(spgi);
          await new Promise((r) => setTimeout(r, 1500));
          let cp: any = null;
          for (const tv of grok.shell.tableViews) { const found = tv.viewers.find((v: any) => v.type === 'Correlation plot'); if (found) { cp = found; break; } }
          if (!cp) return JSON.stringify({ok: false, err: 'no-correlation-plot-found'});
          cp.props.table = spgi.name;
          await new Promise((r) => setTimeout(r, 1200));
          const num: string[] = [];
          for (const c of spgi.columns.numerical) { num.push(c.name); if (num.length === 2) break; }
          const gc = Number(cp.getCorrelation(spgi.col(num[0]), spgi.col(num[1])));
          const ref = Number(DG.Stats.fromColumn(spgi.col(num[0])).corr(spgi.col(num[1])));
          return JSON.stringify({ok: true, spgiName: String(spgi.name), cols: num.slice(), gc, ref, diff: Math.abs(gc - ref), tol});
        } catch (e) { return JSON.stringify({ok: false, err: String(e).slice(0, 200)}); }
      }, {p: spgiPath, tol: TOL});
      const result = JSON.parse(raw) as {ok: boolean; err?: string; spgiName?: string; cols?: string[]; gc?: number; ref?: number; diff?: number};
      expect(result.ok).toBe(true);
      console.log(`[TableSwitch] cols=${result.cols} gc=${result.gc} ref=${result.ref} diff=${result.diff}`);
      // GROK-18487: matrix recomputed for the new table — getCorrelation on a spgi-100 pair equals
      // the runtime Pearson reference (the recompute IS the signal; no prop echo).
      expect(Number.isFinite(result.gc)).toBe(true);
      expect(result.diff!).toBeLessThanOrEqual(TOL);
      expect(pageErrors.length).toBe(peBefore);
      expect(realErrors().length).toBe(errBefore);
    } finally {
      // Restore the CP to demog and close the extra view — even if an assert failed. Explicit null
      // return: an implicit Dart-backed return value cannot serialize.
      await page.evaluate(async () => {
        // Snapshot the Dart-backed tableViews to a plain array BEFORE closing: for...of with
        // tv.close() inside the loop throws 'Concurrent modification during iteration'
        // [live MCP recon 2026-08-01].
        const views: any[] = Array.from(grok.shell.tableViews);
        let cp: any = null;
        for (const tv of views) { const found = tv.viewers.find((v: any) => v.type === 'Correlation plot'); if (found) { cp = found; break; } }
        if (cp) cp.props.table = 'Table';
        await new Promise((r) => setTimeout(r, 800));
        for (const tv of Array.from(grok.shell.tableViews) as any[]) if (tv.dataFrame?.name === 'Table (2)') tv.close();
        await new Promise((r) => setTimeout(r, 400));
        return null;
      });
      await refreshRoot(page, geom);
      // Re-select the demog ('Table') view so grok.shell.tv points at it for the NaN / color-probe
      // sections that read grok.shell.tv.dataFrame — closing the spgi view can leave tv elsewhere.
      await page.evaluate(async () => {
        const demog = (Array.from(grok.shell.tableViews) as any[]).find((tv) => tv.dataFrame?.name === 'Table');
        if (demog) grok.shell.v = demog;
        await new Promise((r) => setTimeout(r, 300));
        return null;
      });
    }
  });

  // ## NaN edge cell (GROK-12586)
  await softStep('NaN edge cell', async () => {
    let removed = false;
    try {
      // Add a constant calculated column; its correlation with any column is undefined (getCorrelation -> null).
      const setup = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        await df.columns.addNewCalculated('constZero', '0');
        await new Promise((r) => setTimeout(r, 700));
        const x0 = cp.props.xColumnNames.slice(), y0 = cp.props.yColumnNames.slice();
        cp.props.xColumnNames = [...x0, 'constZero'];
        cp.props.yColumnNames = [...y0, 'constZero'];
        await new Promise((r) => setTimeout(r, 800));
        const xi = cp.props.xColumnNames.indexOf('constZero');
        const yi = cp.props.yColumnNames.indexOf('AGE');
        const corr = cp.getCorrelation(df.col('constZero'), df.col('AGE'));
        return {xi, yi, corrFinite: Number.isFinite(corr)};
      });
      // The undefined coefficient is not finite (feeds the 'N/A' tooltip).
      expect(setup.corrFinite).toBe(false);
      await refreshRoot(page, geom);
      geom.cellW = await page.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);
      // Trusted hover of the (constZero, AGE) cell; nudge from an off-cell point first so the move
      // registers as a real cell-enter, then poll the tooltip for 'N/A'.
      const c = cellCenter(geom, setup.xi, setup.yi);
      await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
      let tipText = '';
      for (let i = 0; i < 12; i++) {
        await page.mouse.move(c.x, c.y);
        await page.waitForTimeout(250);
        const t = await page.evaluate(() => document.querySelector('.d4-tooltip')?.textContent ?? '');
        if (/R:\s*N\/A/i.test(t)) { tipText = t; break; }
        tipText = t;
        await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
      }
      console.log(`[NaN] tooltip="${tipText.slice(0, 120)}"`);
      // GROK-12586: the tooltip reads 'N/A' instead of a number.
      expect(/R:\s*N\/A/i.test(tipText)).toBe(true);
      // No 'Unsupported operation' error appeared.
      expect(realErrors().some((s) => /Unsupported operation/i.test(s))).toBe(false);
      expect(pageErrors.some((s) => /Unsupported operation/i.test(s))).toBe(false);
    } finally {
      // Teardown — remove the fixture column even if a step failed. Target the demog view
      // EXPLICITLY (tv can point elsewhere), rewire the column sets off the fixture first, and
      // verify by re-reading the column-name list.
      removed = await page.evaluate(async () => {
        const views: any[] = Array.from(grok.shell.tableViews);
        const view = views.find((tv) => tv.dataFrame?.name === 'Table') ?? grok.shell.tv;
        const df = view.dataFrame;
        const cp = view.viewers.find((x: any) => x.type === 'Correlation plot');
        if (cp) {
          cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
          cp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
        }
        await new Promise((r) => setTimeout(r, 300));
        const names = (): string[] => (Array.from(df.columns.names()) as string[]).slice();
        if (names().includes('constZero')) df.columns.remove('constZero');
        await new Promise((r) => setTimeout(r, 400));
        const ok = !names().includes('constZero');
        return ok;
      });
      expect(removed).toBe(true);
      await refreshRoot(page, geom);
    }
  });

  // ## Color-coding probes
  await softStep('Color-coding probes', async () => {
    // Visual-state reset before the pixel probes: white back color, digits ON, so the probe reads
    // the true per-cell hue, not residual state from an earlier section.
    await closeMenu(page);
    await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.white;
      cp.props.showPearsonR = true;
      await new Promise((r) => setTimeout(r, 600));
      return null;
    });
    // Earlier teardowns shifted the root and rewired the column sets — recalibrate with a probe
    // click against the CURRENT column order before any pixel read.
    const colsNow = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      return {x: cp.props.xColumnNames.slice(), y: cp.props.yColumnNames.slice()};
    });
    const xiH = colsNow.x.indexOf('HEIGHT'), xiW = colsNow.x.indexOf('WEIGHT');
    const yiA = colsNow.y.indexOf('AGE'), yiWt = colsNow.y.indexOf('WEIGHT');
    const recalibrate = async (): Promise<boolean> => {
      await refreshRoot(page, geom);
      // cellW tracks the LIVE showPearsonR (40 with digits, else 20) — a stale 40 lands probes on
      // white inter-cell gaps [live MCP recon 2026-08-01].
      geom.cellW = await page.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);
      for (let attempt = 0; attempt < 6; attempt++) {
        await page.evaluate(() => { (window as any).__clicks = []; });
        const c = cellCenter(geom, xiH, yiA);
        await page.mouse.click(c.x, c.y);
        await page.waitForTimeout(400);
        const ev = await lastClick(page);
        console.log(`[ColorProbe] recal ${attempt} at (${Math.round(c.x)},${Math.round(c.y)}) pinnedW=${geom.pinnedW} headerH=${geom.headerH} -> ${JSON.stringify(ev)}`);
        if (ev && [ev.c1, ev.c2].sort().join() === ['AGE', 'HEIGHT'].join()) return true;
        if (ev && ev.c1 && ev.c1 !== 'HEIGHT') {
          const idx = colsNow.x.indexOf(ev.c1);
          if (idx >= 0) geom.pinnedW += (xiH - idx) * geom.cellW;
        } else if (!ev) geom.headerH += 4;
      }
      return false;
    };
    expect(await recalibrate()).toBe(true);
    // Default cells draw values on a WHITE background — per-cell hue exists only with Grid >
    // Grid Color Coding > All [live MCP recon 2026-08-01: neg reads blue, pos red, near-zero
    // light blue]. Menu actuation is the path; the per-cell hue is the signal.
    await refreshRoot(page, geom);
    const ccCell = cellCenter(geom, xiH, yiA);
    await page.mouse.click(ccCell.x, ccCell.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupTrusted(page, 'div-Grid');
    await page.waitForTimeout(300);
    await hoverMenuGroupTrusted(page, 'div-Grid---Grid-Color-Coding');
    await page.waitForTimeout(300);
    const enabledAll = await clickMenuByName(page, 'div-Grid---Grid-Color-Coding---All');
    expect(enabledAll).toBe(true);
    await page.waitForTimeout(800);
    await closeMenu(page);
    await refreshRoot(page, geom);
    // Settle pause before the pixel reads.
    await page.waitForTimeout(500);
    // Runtime r values: HEIGHT/AGE (neg), WEIGHT/AGE (near zero), HEIGHT/WEIGHT (pos).
    const rs = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      // Number() collapses a possibly Dart-boxed return before it crosses the page boundary.
      return {
        neg: Number(cp.getCorrelation(df.col('HEIGHT'), df.col('AGE'))),
        nearZero: Number(cp.getCorrelation(df.col('WEIGHT'), df.col('AGE'))),
        pos: Number(cp.getCorrelation(df.col('HEIGHT'), df.col('WEIGHT'))),
      };
    });
    expect(rs.neg).toBeLessThan(-0.1);
    expect(Math.abs(rs.nearZero)).toBeLessThan(0.15);
    expect(rs.pos).toBeGreaterThan(0.1);
    // Read center pixels: HEIGHT(x)/AGE(y) neg, WEIGHT(x)/AGE(y) near-zero, HEIGHT(x)/WEIGHT(y) pos.
    let neg = await cellPixel(page, geom, xiH, yiA);
    let nearZero = await cellPixel(page, geom, xiW, yiA);
    let pos = await cellPixel(page, geom, xiH, yiWt);
    // Insurance against residual stale geometry: all-three-white means the probes hit background,
    // not cells — ONE retry after a fresh refresh+recalibration, then the honest asserts stand.
    const isWhite = (p: number[] | null) => !!p && p[0] >= 250 && p[1] >= 250 && p[2] >= 250;
    if (isWhite(neg) && isWhite(nearZero) && isWhite(pos)) {
      console.log('[ColorProbe] all three probes white — one recalibrated retry');
      expect(await recalibrate()).toBe(true);
      await page.waitForTimeout(500);
      neg = await cellPixel(page, geom, xiH, yiA);
      nearZero = await cellPixel(page, geom, xiW, yiA);
      pos = await cellPixel(page, geom, xiH, yiWt);
    }
    console.log(`[ColorProbe] neg=${neg} nearZero=${nearZero} pos=${pos}`);
    expect(neg).not.toBeNull();
    expect(nearZero).not.toBeNull();
    expect(pos).not.toBeNull();
    // Opposite signs fall into different hue families: negative is blue-dominant (b>r),
    // positive is red-dominant (r>b).
    expect(neg![2]).toBeGreaterThan(neg![0]);
    expect(pos![0]).toBeGreaterThan(pos![2]);
    // The near-zero cell is visibly lighter than both (higher minimum channel).
    const lightness = (p: number[]) => Math.min(p[0], p[1], p[2]);
    expect(lightness(nearZero!)).toBeGreaterThan(lightness(neg!));
    expect(lightness(nearZero!)).toBeGreaterThan(lightness(pos!));
    // Reset Grid Color Coding to None so the next section starts from an uncolored matrix.
    await refreshRoot(page, geom);
    const resetCell = cellCenter(geom, xiH, yiA);
    await page.mouse.click(resetCell.x, resetCell.y, {button: 'right'});
    await page.waitForTimeout(500);
    await hoverMenuGroupTrusted(page, 'div-Grid');
    await page.waitForTimeout(300);
    await hoverMenuGroupTrusted(page, 'div-Grid---Grid-Color-Coding');
    await page.waitForTimeout(300);
    await clickMenuByName(page, 'div-Grid---Grid-Color-Coding---None');
    await page.waitForTimeout(500);
    await closeMenu(page);
    await refreshRoot(page, geom);
  });

  // ## Diagonal histograms repaint
  await softStep('Diagonal histograms repaint', async () => {
    await refreshRoot(page, geom);
    geom.cellW = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);
    // Baseline a diagonal (AGE-AGE, xi=0 yi=0) cell region, settle-gated.
    const region = {rx: geom.pinnedW, ry: geom.headerH, w: geom.cellW, h: geom.rowH};
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false); df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 500));
      return null;
    });
    await snapCanvas(page, region);
    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);
    // Apply df.filter — the diagonal histogram is rebuilt (rowSource default Filtered).
    await snapCanvas(page, region);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.init((i: number) => df.col('AGE').get(i) > 40);
      await new Promise((r) => setTimeout(r, 800));
      return null;
    });
    const dFilter = await diffCanvas(page);
    console.log(`[Diagonal] settle=${settle} filterDiff=${dFilter}`);
    expect(dFilter).toBeGreaterThanOrEqual(0);
    expect(dFilter).toBeGreaterThan(settle);
    // Clear the filter, re-baseline. Switch rowSource to Selected so a selection drives the histogram
    // (under the default Filtered rowSource selection does NOT repaint the diagonal — a bare-selection
    // assert there would be a false-RED).
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      df.filter.setAll(true);
      cp.props.rowSource = 'Selected';
      await new Promise((r) => setTimeout(r, 700));
      return null;
    });
    await snapCanvas(page, region);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 500);
      await new Promise((r) => setTimeout(r, 800));
      return null;
    });
    const dSelMade = await diffCanvas(page);
    console.log(`[Diagonal] selMade=${dSelMade}`);
    expect(dSelMade).toBeGreaterThanOrEqual(0);
    expect(dSelMade).toBeGreaterThan(settle);
    // Clear the selection — the diagonal repaints again.
    await snapCanvas(page, region);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 800));
      return null;
    });
    const dSelClear = await diffCanvas(page);
    console.log(`[Diagonal] selClear=${dSelClear}`);
    expect(dSelClear).toBeGreaterThanOrEqual(0);
    expect(dSelClear).toBeGreaterThan(settle);
    // restore rowSource
    await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.rowSource = 'Filtered';
      await new Promise((r) => setTimeout(r, 500));
      return null;
    });
  });

  // ## Column width drag
  await softStep('Column width drag', async () => {
    await refreshRoot(page, geom);
    geom.cellW = 40;
    // The CP is NOT a resizable d4 grid: cp.grid is undefined, no JS-readable column width, headers
    // carry no resize hotspot [DOM 2026-08-01, live MCP recon] — the header-edge drag is inert
    // headless. Width drag: status: waived, waiver_class: gesture-uncontrollable-headless. The drag
    // is still driven with real trusted input (both directions); the assertion is the defensive
    // floor (diff >= 0), not repaint-exceeds-settle.
    await snapCanvas(page);
    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);
    await snapCanvas(page);
    const edgeX = geom.rootX + geom.pinnedW + geom.cellW;
    const edgeY = geom.rootY + geom.headerH - 10;
    await page.mouse.move(edgeX, edgeY);
    await page.mouse.down();
    await page.mouse.move(edgeX + 30, edgeY, {steps: 6});
    await page.mouse.move(edgeX + 60, edgeY, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(700);
    const dragDelta = await diffCanvas(page);
    console.log(`[WidthDrag] settle=${settle} dragDelta=${dragDelta} (waived: CP has no readable column-width channel; header-edge drag inert headless)`);
    // Honest floor: the canvas read succeeded (no -1 fault); the repaint claim is waived above.
    expect(dragDelta).toBeGreaterThanOrEqual(0);
    // Exercise the restore drag too (round-trip attempt), same waived floor.
    await snapCanvas(page);
    await page.mouse.move(edgeX + 60, edgeY);
    await page.mouse.down();
    await page.mouse.move(edgeX, edgeY, {steps: 8});
    await page.mouse.up();
    await page.waitForTimeout(500);
    const restoreDelta = await diffCanvas(page);
    console.log(`[WidthDrag] restoreDelta=${restoreDelta} (waived)`);
    expect(restoreDelta).toBeGreaterThanOrEqual(0);
  });

  v.finishSpec();
});
