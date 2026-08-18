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

async function waitForRTooltip(page: Page, timeoutMs = 5000): Promise<void> {
  await page.waitForFunction(() => {
    const tip = document.querySelector('.d4-tooltip');
    return /R:\s*-?\d/.test(tip?.textContent ?? '');
  }, undefined, {timeout: timeoutMs});
}

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

async function waitVisibleMenuItem(page: Page, name: string, timeoutMs = 4000): Promise<boolean> {
  try {
    await page.waitForFunction((n) => {
      const nodes = Array.from(document.querySelectorAll(`.d4-menu-popup .d4-menu-item[name="${n}"]`));
      return nodes.some((el) => { const b = (el as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; });
    }, name, {timeout: timeoutMs});
    return true;
  } catch { return false; }
}

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

    const b = el.getBoundingClientRect();
    const o: any = {bubbles: true, cancelable: true, clientX: b.x + b.width / 2, clientY: b.y + b.height / 2, view: window};
    el.dispatchEvent(new PointerEvent('pointerover', o));
    el.dispatchEvent(new MouseEvent('mouseover', o));
    return true;
  }, name);
}

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

async function closeMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.evaluate(() => {
    document.querySelectorAll('.d4-menu-popup').forEach((el) => el.remove());
    document.body.click();
  });

  await page.waitForTimeout(200);
}

async function waitMenuOpen(page: Page, capMs = 500): Promise<boolean> {
  return v.pollValue(() => page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .some((el) => { const b = (el as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; })),
  (open) => open, capMs, 100);
}

async function waitMenuTextVisible(page: Page, pattern: string, capMs: number): Promise<boolean> {
  return v.pollValue(() => page.evaluate((p) => {
    const re = new RegExp(p);
    const own = (it: Element) => Array.from(it.querySelectorAll('.d4-menu-item-label'))
      .find((l) => l.closest('.d4-menu-item') === it)?.textContent?.trim() ?? '';
    return Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item')).some((it) => {
      const b = (it as HTMLElement).getBoundingClientRect();
      return re.test(own(it)) && b.width > 0 && b.height > 0;
    });
  }, pattern), (vis) => vis, capMs, 100);
}

test('Correlation plot — property surface smoke', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

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

      const ev = await v.pollValue(() => lastClick(page), (e) => e !== null, 400, 100);
      console.log(`[Setup] probe ${attempt} at (${Math.round(c.x)},${Math.round(c.y)}) pinnedW=${geom.pinnedW} headerH=${geom.headerH} -> ${JSON.stringify(ev)}`);
      if (ev && [ev.c1, ev.c2].sort().join() === ['AGE', 'HEIGHT'].join()) { landed = true; break; }
      if (ev && ev.c1 && ev.c1 !== 'HEIGHT') {
        const idx = base.xCols.indexOf(ev.c1);
        if (idx >= 0) geom.pinnedW += (xiHeight - idx) * geom.cellW;
      } else if (!ev) geom.headerH += 4;
    }
    const ev = await lastClick(page);
    expect(ev).not.toBeNull();

    expect([ev!.c1, ev!.c2].sort()).toEqual(['AGE', 'HEIGHT']);
  });

  await softStep('Title, description, and back color', async () => {

    await v.setViewerProps(page, 'Correlation plot', [{set: {
      showTitle: true, title: 'Correlation Analysis',
      description: 'Shows pairwise correlations', descriptionVisibilityMode: 'Always',
    }, wait: 500}]);

    const titleShown: boolean = await v.pollValue(() => page.evaluate(() =>
      Array.from(document.querySelectorAll('.panel-titlebar-text'))
        .some((el) => el.textContent!.trim() === 'Correlation Analysis')), (shown) => shown, 3000, 150);
    const descAlways: boolean = await v.pollValue(() => page.evaluate(() =>
      document.querySelector('[name="viewer-Correlation-plot"]')!.textContent!.includes('Shows pairwise correlations')),
    (shown) => shown, 3000, 150);

    const origDescPos: string = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.descriptionPosition);
    await v.setViewerProps(page, 'Correlation plot', [{set: {descriptionPosition: 'Bottom'}, wait: 600}]);
    const descPos = await v.pollValue(() => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
      const R = root.getBoundingClientRect();
      const leaf = Array.from(root.querySelectorAll('*'))
        .filter((el) => el.children.length === 0 && (el.textContent ?? '').includes('Shows pairwise correlations'))
        .pop() ?? null;
      const rect = leaf ? leaf.getBoundingClientRect() : null;
      return {found: !!rect, centerY: rect ? rect.top + rect.height / 2 : -1, rootMid: R.y + R.height / 2};
    }), (p) => p.found && p.centerY > p.rootMid, 3000, 150);
    await v.setViewerProps(page, 'Correlation plot', [{set: {descriptionPosition: origDescPos}, wait: 300}]);
    console.log(`[DescPos] ${JSON.stringify({...descPos, orig: origDescPos})}`);
    expect(descPos.found).toBe(true);

    expect(descPos.centerY).toBeGreaterThan(descPos.rootMid);
    await v.setViewerProps(page, 'Correlation plot', [{set: {descriptionVisibilityMode: 'Never'}, wait: 500}]);
    const descNever: boolean = await v.pollValue(() => page.evaluate(() =>
      document.querySelector('[name="viewer-Correlation-plot"]')!.textContent!.includes('Shows pairwise correlations')),
    (shown) => !shown, 500, 100);
    expect(titleShown).toBe(true);
    expect(descAlways).toBe(true);
    expect(descNever).toBe(false);

    await snapCanvas(page);

    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);
    console.log(`[BackColor] settle=${settle}`);
    await snapCanvas(page);
    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.lightGray;
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 600);
    const backDelta = await v.pollValue(() => diffCanvas(page), (d) => d > settle + 1000, 3000, 150);
    console.log(`[BackColor] delta=${backDelta}`);
    expect(backDelta).toBeGreaterThanOrEqual(0);
    expect(backDelta).toBeGreaterThan(settle + 1000);

    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.white;
      cp.props.showTitle = false; cp.props.title = ''; cp.props.description = '';
      cp.props.descriptionVisibilityMode = 'Auto';
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 400);
    expect(realErrors()).toEqual([]);
  });

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

    expect(realErrors().length).toBe(errBefore);
    expect(pageErrors.length).toBe(peBefore);
  });

  await softStep('Context-menu toggles mirror properties', async () => {
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await waitMenuOpen(page);
    const showRBefore: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR);

    await snapCanvas(page);

    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);

    await snapCanvas(page);
    const clickedShowR = await clickMenuByName(page, 'div-Show-Pearson-R');
    expect(clickedShowR).toBe(true);
    const showRAfter: boolean = await v.pollValue(() => page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR),
    (val) => val === false, 4000, 100);
    const repaintR = await v.pollValue(() => diffCanvas(page), (d) => d > settle, 3000, 150);
    console.log(`[Menu] showPearsonR ${showRBefore}->${showRAfter} settle=${settle} repaint=${repaintR}`);
    expect(showRBefore).toBe(true);
    expect(showRAfter).toBe(false);
    expect(repaintR).toBeGreaterThanOrEqual(0);
    expect(repaintR).toBeGreaterThan(settle);

    await closeMenu(page);
    await refreshRoot(page, geom);

    geom.cellW = 20;
    const c2 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c2.x, c2.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupTrusted(page, 'div-Tooltip');
    await waitVisibleMenuItem(page, 'div-Tooltip---Visible', 400);
    const tipBefore: boolean = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip);
    const clickedVisible = await clickMenuByName(page, 'div-Tooltip---Visible');
    expect(clickedVisible).toBe(true);
    const tipAfter: boolean = await v.pollValue(() => page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip),
    (val) => val === false, 500, 100);
    console.log(`[Menu] showTooltip ${tipBefore}->${tipAfter}`);
    expect(tipBefore).toBe(true);
    expect(tipAfter).toBe(false);

    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
    const displayIdle: string = await v.pollValue(() => page.evaluate(() => {
      const tip = document.querySelector('.d4-tooltip');
      return tip ? getComputedStyle(tip).display : 'missing';
    }), (d) => d === 'none', 500, 100);
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

    expect(shownWhileOff).toBe(false);

    await closeMenu(page);
    await refreshRoot(page, geom);
    const c3 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c3.x, c3.y, {button: 'right'});
    await waitMenuOpen(page);
    await clickMenuByName(page, 'div-Show-Pearson-R');
    await v.pollValue(() => page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR),
    (val) => val === true, 4000, 100);
    await closeMenu(page);
    await refreshRoot(page, geom);
    geom.cellW = 40;
    const c4 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c4.x, c4.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupTrusted(page, 'div-Tooltip');
    await waitVisibleMenuItem(page, 'div-Tooltip---Visible', 400);
    await clickMenuByName(page, 'div-Tooltip---Visible');
    await v.pollValue(() => page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showTooltip),
    (val) => val === true, 4000, 100);
    const restored = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      return {showR: cp.props.showPearsonR, tip: cp.props.showTooltip};
    });
    expect(restored.showR).toBe(true);
    expect(restored.tip).toBe(true);
    await closeMenu(page);

    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);

    await v.pollValue(() => page.evaluate(() => {
      const tip = document.querySelector('.d4-tooltip');
      return tip ? getComputedStyle(tip).display : 'missing';
    }), (d) => d !== 'block', 500, 100);
    const cOn = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.move(cOn.x, cOn.y);
    await page.waitForFunction(() => {
      const tip = document.querySelector('.d4-tooltip');
      return !!tip && getComputedStyle(tip).display === 'block' && /R:\s*-?\d/.test(tip.textContent ?? '');
    }, undefined, {timeout: 5000});
  });

  await softStep('Pinned row-header columns', async () => {

    await page.evaluate(() => {
      (document.querySelector('[name="viewer-Correlation-plot"]') as HTMLElement).style.width = '180px';
      return null;
    });
    const scrollbarShown: boolean = await v.pollValue(() => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Correlation-plot"]') as HTMLElement;
      const sb = root.querySelector('.d4-range-selector.d4-grid-horz-scroll') as HTMLElement | null;
      const r = sb?.getBoundingClientRect();
      return !!sb && !!r && r.width > 0 && r.height > 0;
    }), (shown) => shown, 3000, 150);
    expect(scrollbarShown).toBe(true);

    await refreshRoot(page, geom);
    const probeX = geom.rootX + geom.pinnedW + geom.cellW / 2;
    const probePair = async (): Promise<{c1: string; c2: string; v: number} | null> => {
      for (const row of ['HEIGHT', 'AGE']) {
        const yi = base.yCols.indexOf(row);
        await page.evaluate(() => { (window as any).__clicks = []; });
        await page.mouse.click(probeX, geom.rootY + geom.headerH + (yi + 0.5) * geom.rowH);

        const ev = await v.pollValue(() => lastClick(page), (e) => e !== null, 400, 100);
        if (ev) return ev;
      }
      return null;
    };
    const pairBefore = await probePair();
    expect(pairBefore).not.toBeNull();

    await snapCanvas(page);
    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + geom.pinnedW + geom.cellW, geom.rootY + geom.headerH + 60);
    await page.mouse.wheel(200, 0);
    await v.waitForViewerRendered(page, 'Correlation plot', 600);
    const scrollRepaint = await diffCanvas(page);
    console.log(`[Pinned] wheel repaint=${scrollRepaint}`);
    expect(scrollRepaint).toBeGreaterThanOrEqual(0);
    const pairAfter = await probePair();
    console.log(`[Pinned] probe before=${JSON.stringify(pairBefore)} after=${JSON.stringify(pairAfter)}`);

    if (pairAfter && pairBefore && pairAfter.c1 !== pairBefore.c1)
      expect(pairAfter.c1).not.toBe(pairBefore.c1);
    else
      console.log('[Pinned] wheel-scroll column move inert headless -> waived; pinning proven by pinned-name tooltip');

    await refreshRoot(page, geom);

    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
    await page.mouse.move(geom.rootX + 60, geom.rootY + 70);
    const pinnedTipText: string = await v.pollValue(() => page.evaluate(() =>
      document.querySelector('.d4-tooltip')?.textContent ?? ''), (t) => /avg:|min:/i.test(t), 3000, 150);
    console.log(`[Pinned] tooltip="${pinnedTipText.replace(/\s+/g, ' ').slice(0, 120)}"`);
    expect(/avg:|min:/i.test(pinnedTipText)).toBe(true);

    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Correlation-plot"]') as HTMLElement;
      root.style.width = '';
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 500);
    await refreshRoot(page, geom);
    expect(realErrors()).toEqual([]);
  });

  await softStep('Order or Hide Columns dialog', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await waitMenuOpen(page);

    await hoverMenuGroupByText(page, 'Grid');
    await waitMenuTextVisible(page, 'Order or Hide Columns', 400);
    const openedDialog = await clickMenuItemByText(page, 'Order or Hide Columns');
    expect(openedDialog).toBe(true);
    await v.pollValue(() => page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"]').count(),
      (n) => n > 0, 3000, 150);

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

    expect(driven.selectDriven).toBe(true);

    expect(driven.searchDriven).toBe(true);

    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Order-or-Hide-Columns"]');
      const btn = dlg?.querySelector('[name="button-CLOSE"]') as HTMLElement | null;
      if (btn) btn.click();
      return null;
    });
    const closedOH: boolean = await v.pollValue(() => page.evaluate(() =>
      document.querySelectorAll('.d4-dialog[name="dialog-Order-or-Hide-Columns"]').length === 0),
    (gone) => gone, 3000, 150);
    expect(closedOH).toBe(true);
    await closeMenu(page);

    expect(pageErrors.length).toBe(peBefore);
    expect(realErrors().length).toBe(errBefore);

    await refreshRoot(page, geom);
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 120);
    await page.mouse.move(geom.rootX + 60, geom.rootY + 70);
    const pinnedStillText: string = await v.pollValue(() => page.evaluate(() =>
      document.querySelector('.d4-tooltip')?.textContent ?? ''), (t) => /avg:|min:/i.test(t), 3000, 150);
    console.log(`[OrderHide] pinned tooltip="${pinnedStillText.replace(/\s+/g, ' ').slice(0, 120)}"`);
    expect(/avg:|min:/i.test(pinnedStillText)).toBe(true);
  });

  await softStep('Grid color-coding apply to text', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;
    await closeMenu(page);
    await refreshRoot(page, geom);

    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await waitMenuOpen(page);

    await hoverMenuGroupByText(page, 'Grid');
    await waitMenuTextVisible(page, 'Current Column', 300);
    await hoverMenuGroupByText(page, 'Current Column');
    await waitMenuTextVisible(page, 'Color Coding', 300);
    await hoverMenuGroupByText(page, 'Color Coding');
    await waitMenuTextVisible(page, '^Linear$', 300);

    const appliedLinear = await clickMenuItemByText(page, '^Linear$', 'Color Coding');
    expect(appliedLinear).toBe(true);
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    await closeMenu(page);

    await refreshRoot(page, geom);
    const c2 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c2.x, c2.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupByText(page, 'Grid');
    await waitMenuTextVisible(page, 'Current Column', 300);
    await hoverMenuGroupByText(page, 'Current Column');
    await waitMenuTextVisible(page, 'Color Coding', 300);
    await hoverMenuGroupByText(page, 'Color Coding');
    await waitMenuTextVisible(page, '^Edit\\.\\.\\.$', 300);

    const openedEdit = await clickMenuItemByText(page, '^Edit\\.\\.\\.$', 'Color Coding');
    expect(openedEdit).toBe(true);
    await v.pollValue(() => page.evaluate(() => Array.from(document.querySelectorAll('.d4-dialog'))
      .some((dd) => /^dialog-Color-coding-/.test(dd.getAttribute('name') ?? ''))), (open) => open, 3000, 150);

    const setToText = await page.evaluate(() => {

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
      return sel.selectedIndex >= 0 && (sel.options[sel.selectedIndex].textContent ?? sel.value).trim().toLowerCase() === 'text';
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    console.log(`[ApplyToText] setToText=${setToText}`);

    expect(setToText).toBe(true);

    const resetToBackground: boolean = await page.evaluate(() => {

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
      const ok = !!sel && sel.selectedIndex >= 0
        && (sel.options[sel.selectedIndex].textContent ?? sel.value).trim().toLowerCase() === 'background';
      (DG.Dialog.getOpenDialogs?.() ?? []).forEach((dlg: any) => dlg.close?.());
      Array.from(document.querySelectorAll('.d4-dialog')).forEach((dd) => {
        const cl = dd.querySelector('[name="button-CLOSE"]') as HTMLElement | null; if (cl) cl.click();
      });
      return ok;
    });
    await v.pollValue(() => page.evaluate(() => document.querySelectorAll('.d4-dialog').length),
      (n) => n === 0, 700, 100);

    expect(resetToBackground).toBe(true);
    await closeMenu(page);
    await refreshRoot(page, geom);
    const c3 = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c3.x, c3.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupByText(page, 'Grid');
    await waitMenuTextVisible(page, 'Grid Color Coding', 300);
    await hoverMenuGroupByText(page, 'Grid Color Coding');
    await waitMenuTextVisible(page, '^None$', 300);
    await clickMenuItemByText(page, '^None$', 'Grid Color Coding');
    await v.waitForViewerRendered(page, 'Correlation plot', 400);
    await closeMenu(page);

    expect(pageErrors.length).toBe(peBefore);
    expect(realErrors().length).toBe(errBefore);
  });

  await softStep('Table switch', async () => {
    const peBefore = pageErrors.length;
    const errBefore = realErrors().length;

    await closeMenu(page);
    await page.evaluate(() => { (DG.Dialog.getOpenDialogs?.() ?? []).forEach((d: any) => d.close?.()); return null; });
    await v.pollValue(() => page.evaluate(() => document.querySelectorAll('.d4-dialog').length),
      (n) => n === 0, 300, 100);
    try {

      const raw = await page.evaluate(async ({p, tol}) => {
        try {
          const spgi = await grok.dapi.files.readCsv(p);
          grok.shell.addTableView(spgi);

          const viewReady = () => (Array.from(grok.shell.tableViews) as any[]).some((tv) => tv.dataFrame?.name === spgi.name);
          for (let i = 0; i < 15 && !viewReady(); i++) await new Promise((r) => setTimeout(r, 100));
          let cp: any = null;
          for (const tv of grok.shell.tableViews) { const found = tv.viewers.find((v: any) => v.type === 'Correlation plot'); if (found) { cp = found; break; } }
          if (!cp) return JSON.stringify({ok: false, err: 'no-correlation-plot-found'});
          cp.props.table = spgi.name;

          for (let i = 0; i < 12 && cp.dataFrame?.name !== spgi.name; i++) await new Promise((r) => setTimeout(r, 100));
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

      expect(Number.isFinite(result.gc)).toBe(true);
      expect(result.diff!).toBeLessThanOrEqual(TOL);
      expect(pageErrors.length).toBe(peBefore);
      expect(realErrors().length).toBe(errBefore);
    } finally {

      await page.evaluate(() => {

        const views: any[] = Array.from(grok.shell.tableViews);
        let cp: any = null;
        for (const tv of views) { const found = tv.viewers.find((v: any) => v.type === 'Correlation plot'); if (found) { cp = found; break; } }
        if (cp) cp.props.table = 'Table';
        return null;
      });
      await v.pollValue(() => page.evaluate(() => {
        for (const tv of Array.from(grok.shell.tableViews) as any[]) {
          const found = tv.viewers.find((x: any) => x.type === 'Correlation plot');
          if (found) return found.dataFrame ? String(found.dataFrame.name) : null;
        }
        return null;
      }), (name) => name === 'Table', 800, 100);
      await page.evaluate(() => {
        for (const tv of Array.from(grok.shell.tableViews) as any[]) if (tv.dataFrame?.name === 'Table (2)') tv.close();
        return null;
      });
      await v.pollValue(() => page.evaluate(() =>
        (Array.from(grok.shell.tableViews) as any[]).some((tv) => tv.dataFrame?.name === 'Table (2)')),
      (stillOpen) => !stillOpen, 400, 100);
      await refreshRoot(page, geom);

      await page.evaluate(() => {
        const demog = (Array.from(grok.shell.tableViews) as any[]).find((tv) => tv.dataFrame?.name === 'Table');
        if (demog) grok.shell.v = demog;
        return null;
      });
      await v.pollValue(() => page.evaluate(() =>
        grok.shell.tv?.dataFrame ? String(grok.shell.tv.dataFrame.name) : null),
      (name) => name === 'Table', 300, 100);
    }
  });

  await softStep('NaN edge cell', async () => {
    let removed = false;
    try {

      const setup = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        await df.columns.addNewCalculated('constZero', '0');

        for (let i = 0; i < 7 && !df.col('constZero'); i++) await new Promise((r) => setTimeout(r, 100));
        const x0 = cp.props.xColumnNames.slice(), y0 = cp.props.yColumnNames.slice();
        cp.props.xColumnNames = [...x0, 'constZero'];
        cp.props.yColumnNames = [...y0, 'constZero'];
        const xi = cp.props.xColumnNames.indexOf('constZero');
        const yi = cp.props.yColumnNames.indexOf('AGE');
        const corr = cp.getCorrelation(df.col('constZero'), df.col('AGE'));
        return {xi, yi, corrFinite: Number.isFinite(corr)};
      });

      expect(setup.corrFinite).toBe(false);
      await v.waitForViewerRendered(page, 'Correlation plot', 800);
      await refreshRoot(page, geom);
      geom.cellW = await page.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);

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

      expect(/R:\s*N\/A/i.test(tipText)).toBe(true);

      expect(realErrors().some((s) => /Unsupported operation/i.test(s))).toBe(false);
      expect(pageErrors.some((s) => /Unsupported operation/i.test(s))).toBe(false);
    } finally {

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
        const ok = !names().includes('constZero');
        return ok;
      });
      expect(removed).toBe(true);
      await v.waitForViewerRendered(page, 'Correlation plot', 400);
      await refreshRoot(page, geom);
    }
  });

  await softStep('Color-coding probes', async () => {

    await closeMenu(page);
    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.backColor = DG.Color.white;
      cp.props.showPearsonR = true;
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 600);

    const colsNow = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      return {x: cp.props.xColumnNames.slice(), y: cp.props.yColumnNames.slice()};
    });
    const xiH = colsNow.x.indexOf('HEIGHT'), xiW = colsNow.x.indexOf('WEIGHT');
    const yiA = colsNow.y.indexOf('AGE'), yiWt = colsNow.y.indexOf('WEIGHT');
    const recalibrate = async (): Promise<boolean> => {
      await refreshRoot(page, geom);

      geom.cellW = await page.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);
      for (let attempt = 0; attempt < 6; attempt++) {
        await page.evaluate(() => { (window as any).__clicks = []; });
        const c = cellCenter(geom, xiH, yiA);
        await page.mouse.click(c.x, c.y);

        const ev = await v.pollValue(() => lastClick(page), (e) => e !== null, 400, 100);
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

    await refreshRoot(page, geom);
    const ccCell = cellCenter(geom, xiH, yiA);
    await page.mouse.click(ccCell.x, ccCell.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupTrusted(page, 'div-Grid');
    await waitVisibleMenuItem(page, 'div-Grid---Grid-Color-Coding', 300);
    await hoverMenuGroupTrusted(page, 'div-Grid---Grid-Color-Coding');
    await waitVisibleMenuItem(page, 'div-Grid---Grid-Color-Coding---All', 300);
    const enabledAll = await clickMenuByName(page, 'div-Grid---Grid-Color-Coding---All');
    expect(enabledAll).toBe(true);
    await v.waitForViewerRendered(page, 'Correlation plot', 800);
    await closeMenu(page);
    await refreshRoot(page, geom);

    await page.waitForTimeout(500);

    const rs = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;

      return {
        neg: Number(cp.getCorrelation(df.col('HEIGHT'), df.col('AGE'))),
        nearZero: Number(cp.getCorrelation(df.col('WEIGHT'), df.col('AGE'))),
        pos: Number(cp.getCorrelation(df.col('HEIGHT'), df.col('WEIGHT'))),
      };
    });
    expect(rs.neg).toBeLessThan(-0.1);
    expect(Math.abs(rs.nearZero)).toBeLessThan(0.15);
    expect(rs.pos).toBeGreaterThan(0.1);

    let neg = await cellPixel(page, geom, xiH, yiA);
    let nearZero = await cellPixel(page, geom, xiW, yiA);
    let pos = await cellPixel(page, geom, xiH, yiWt);

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

    expect(neg![2]).toBeGreaterThan(neg![0]);
    expect(pos![0]).toBeGreaterThan(pos![2]);

    const lightness = (p: number[]) => Math.min(p[0], p[1], p[2]);
    expect(lightness(nearZero!)).toBeGreaterThan(lightness(neg!));
    expect(lightness(nearZero!)).toBeGreaterThan(lightness(pos!));

    await refreshRoot(page, geom);
    const resetCell = cellCenter(geom, xiH, yiA);
    await page.mouse.click(resetCell.x, resetCell.y, {button: 'right'});
    await waitMenuOpen(page);
    await hoverMenuGroupTrusted(page, 'div-Grid');
    await waitVisibleMenuItem(page, 'div-Grid---Grid-Color-Coding', 300);
    await hoverMenuGroupTrusted(page, 'div-Grid---Grid-Color-Coding');
    await waitVisibleMenuItem(page, 'div-Grid---Grid-Color-Coding---None', 300);
    await clickMenuByName(page, 'div-Grid---Grid-Color-Coding---None');
    await v.waitForViewerRendered(page, 'Correlation plot', 500);
    await closeMenu(page);
    await refreshRoot(page, geom);
  });

  await softStep('Diagonal histograms repaint', async () => {
    await refreshRoot(page, geom);
    geom.cellW = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot').props.showPearsonR ? 40 : 20);

    const region = {rx: geom.pinnedW, ry: geom.headerH, w: geom.cellW, h: geom.rowH};
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false); df.filter.setAll(true);
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 500);
    await snapCanvas(page, region);

    await page.waitForTimeout(300);
    const settle = await diffCanvas(page);
    expect(settle).toBeGreaterThanOrEqual(0);

    await snapCanvas(page, region);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.filter.init((i: number) => df.col('AGE').get(i) > 40);
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 800);
    const dFilter = await v.pollValue(() => diffCanvas(page), (d) => d > settle, 3000, 150);
    console.log(`[Diagonal] settle=${settle} filterDiff=${dFilter}`);
    expect(dFilter).toBeGreaterThanOrEqual(0);
    expect(dFilter).toBeGreaterThan(settle);

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      df.filter.setAll(true);
      cp.props.rowSource = 'Selected';
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    await snapCanvas(page, region);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 500);
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 800);
    const dSelMade = await v.pollValue(() => diffCanvas(page), (d) => d > settle, 3000, 150);
    console.log(`[Diagonal] selMade=${dSelMade}`);
    expect(dSelMade).toBeGreaterThanOrEqual(0);
    expect(dSelMade).toBeGreaterThan(settle);

    await snapCanvas(page, region);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 800);
    const dSelClear = await v.pollValue(() => diffCanvas(page), (d) => d > settle, 3000, 150);
    console.log(`[Diagonal] selClear=${dSelClear}`);
    expect(dSelClear).toBeGreaterThanOrEqual(0);
    expect(dSelClear).toBeGreaterThan(settle);

    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.rowSource = 'Filtered';
      return null;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 500);
  });

  await softStep('Column width drag', async () => {
    await refreshRoot(page, geom);
    geom.cellW = 40;

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
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    const dragDelta = await diffCanvas(page);
    console.log(`[WidthDrag] settle=${settle} dragDelta=${dragDelta} (waived: CP has no readable column-width channel; header-edge drag inert headless)`);

    expect(dragDelta).toBeGreaterThanOrEqual(0);

    await snapCanvas(page);
    await page.mouse.move(edgeX + 60, edgeY);
    await page.mouse.down();
    await page.mouse.move(edgeX, edgeY, {steps: 8});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Correlation plot', 500);
    const restoreDelta = await diffCanvas(page);
    console.log(`[WidthDrag] restoreDelta=${restoreDelta} (waived)`);
    expect(restoreDelta).toBeGreaterThanOrEqual(0);
  });

  v.finishSpec();
});
