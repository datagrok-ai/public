/* ---
realizes: [grid.cp.cell-appearance, grid.int.color-resolution-order]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, 0).documentBounds;
    const headerY = db.y - grid.colHeaderHeight / 2;
    return {x: db.x + db.width / 2, y: headerY};
  }, col);
}

async function cellCenter(page: Page, col: string, row: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, r}) => {
    const db = grok.shell.tv.grid.cell(c, r).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, r: row});
}

async function laidOutRect(
  page: Page, name: string,
): Promise<{x: number; y: number; w: number; h: number} | null> {
  return page.evaluate((n) => {
    const els = Array.from(document.querySelectorAll('[name="' + n + '"]')) as HTMLElement[];
    for (const el of els) {
      if (el.offsetParent === null) continue;
      const r = el.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) return {x: r.x, y: r.y, w: r.width, h: r.height};
    }
    return null;
  }, name);
}

async function clickMenuLeaf(
  page: Page, at: {x: number; y: number}, groupNames: string[], leafName: string,
): Promise<boolean> {
  for (let open = 0; open < 5; open++) {
    await page.evaluate(({x, y}) => {
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', cm));
      overlay.dispatchEvent(new MouseEvent('mouseup', cm));
      overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
    }, {x: at.x, y: at.y});
    await page.waitForTimeout(550); 

    let ok = true;
    for (const g of groupNames) {
      const gr = await laidOutRect(page, g);
      if (!gr) { ok = false; break; }
      const gcx = gr.x + gr.w / 2; const gcy = gr.y + gr.h / 2;
      await page.mouse.move(gr.x + 4, gcy, {steps: 4});
      await page.mouse.move(gcx, gcy, {steps: 8});
      await page.waitForTimeout(450); 
    }

    const leaf = ok ? await laidOutRect(page, leafName) : null;
    if (!leaf) {
      await page.evaluate(() => (document.body as HTMLElement).click());
      await page.waitForTimeout(250); 
      continue;
    }
    await page.mouse.move(leaf.x + leaf.w / 2, leaf.y + leaf.h / 2, {steps: 4});
    await page.mouse.click(leaf.x + leaf.w / 2, leaf.y + leaf.h / 2);
    await page.waitForTimeout(800); 
    return true;
  }
  return false;
}

async function openGridSettings(page: Page): Promise<boolean> {
  const rows = page.locator('[name="prop-color-coding"]');
  if (await rows.count() > 0) return true;
  const gestures: Array<() => Promise<void>> = [
    async () => {
      const box = await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return null;
        const r = gear.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      });
      if (box === null) return;
      await page.mouse.move(box.x - 40, box.y + 20);
      await page.mouse.move(box.x, box.y, {steps: 6});
      await page.waitForTimeout(250); 
      await page.mouse.click(box.x, box.y);
    },
    async () => {
      await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return;
        const r = gear.getBoundingClientRect();
        const cx = r.x + r.width / 2; const cy = r.y + r.height / 2;
        const host = gear.parentElement ?? gear;
        host.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: cx, clientY: cy}));
        gear.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: cx, clientY: cy}));
        for (const type of ['mousedown', 'mouseup', 'click'])
          gear.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0}));
      });
    },
    async () => { await v.openViewerGear(page, 'Grid'); },
    async () => {
      await page.locator('[name="viewer-Grid"] canvas[name="overlay"]').first().click({position: {x: 60, y: 60}, force: true});
      await page.keyboard.press('F4');
    },
  ];
  for (const gesture of gestures) {
    try {
      await gesture();
    } catch {
      continue;
    }
    try {
      await rows.first().waitFor({state: 'attached', timeout: 6000});
      return true;
    } catch {

    }
  }
  return false;
}

async function snapGridCanvas(page: Page): Promise<boolean> {
  return page.evaluate(() => {
    const cv = document.querySelector('[name="viewer-Grid"] canvas[name="canvas"]') as HTMLCanvasElement | null;
    if (!cv) return false;
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4)
        colors.set((img[i] << 16) | (img[i + 1] << 8) | img[i + 2], 1);
      (window as any).__gridCanvasSnap = img;
      (window as any).__gridCanvasW = cv.width; (window as any).__gridCanvasH = cv.height;
      return true;
    } catch { return false; }
  });
}

async function diffGridCanvas(page: Page): Promise<number> {
  return page.evaluate(() => {
    const w = window as any;
    const cv = document.querySelector('[name="viewer-Grid"] canvas[name="canvas"]') as HTMLCanvasElement | null;
    const prev = w.__gridCanvasSnap as Uint8ClampedArray | undefined;
    if (!cv || !prev) return -1;
    try {
      if (cv.width !== w.__gridCanvasW || cv.height !== w.__gridCanvasH) return -1;
      const cur = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      let delta = 0;
      for (let i = 0; i < cur.length; i += 4) {
        if (cur[i] !== prev[i] || cur[i + 1] !== prev[i + 1] || cur[i + 2] !== prev[i + 2]) delta++;
      }
      return delta;
    } catch { return -1; }
  });
}

test('Grid — Cell Appearance and Color Resolution Order', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  const WHITE = 0xffffffff;
  const consoleErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const baselineErrors = consoleErrors.length;

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const ac = df.col('AGE'); const hc = df.col('HEIGHT');
    let minAgeRow = -1; let maxAgeRow = -1; let minV = Infinity; let maxV = -Infinity; let nullAgeRow = -1;
    for (let i = 0; i < df.rowCount; i++) {
      if (ac.isNone(i)) { if (nullAgeRow < 0) nullAgeRow = i; continue; }
      const val = ac.get(i);
      if (val < minV) { minV = val; minAgeRow = i; }
      if (val > maxV) { maxV = val; maxAgeRow = i; }
    }
    let minHRow = -1; let maxHRow = -1; let hMinV = Infinity; let hMaxV = -Infinity;
    for (let i = 0; i < df.rowCount; i++) {
      if (hc.isNone(i)) continue;
      const val = hc.get(i);
      if (val < hMinV) { hMinV = val; minHRow = i; }
      if (val > hMaxV) { hMaxV = val; maxHRow = i; }
    }

    let targetRow = -1;
    for (let i = 0; i < df.rowCount; i++) if (!ac.isNone(i)) { targetRow = i; break; }
    return {minAgeRow, maxAgeRow, nullAgeRow, minHRow, maxHRow, targetRow};
  });
  expect(setup.minAgeRow).not.toBe(setup.maxAgeRow); 
  expect(setup.nullAgeRow).toBeGreaterThanOrEqual(0); 
  expect(setup.minHRow).not.toBe(setup.maxHRow); 

  await softStep('Step 1-2 — Per-column Linear on AGE via the header menu: min/max cells differ, both differ from background', async () => {
    const c = await headerCenter(page, 'AGE');
    const applied = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Linear');
    expect(applied).toBe(true); 
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
      return {
        ccType: df.col('AGE').getTag('.color-coding-type'),
        minColor: grid.cell('AGE', s.minAgeRow).color >>> 0,
        maxColor: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.ccType).toBe('Linear'); 
    expect(r.minColor).not.toBe(r.maxColor); 
    expect(r.minColor).not.toBe(WHITE); 
    expect(r.maxColor).not.toBe(WHITE);
  });

  await softStep('Step 3-4 — Grid Color Coding None: AGE cells revert to the plain background', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---None');
    expect(applied).toBe(true); 
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        gridCC: grid.props.colorCoding,
        ageMin: grid.cell('AGE', s.minAgeRow).color >>> 0,
        ageMax: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.gridCC).toBe('None'); 

    expect(r.ageMin).toBe(WHITE);
    expect(r.ageMax).toBe(WHITE);
  });

  await softStep('Step 5-6 — Grid Color Coding All: HEIGHT (previously uncolored) now auto-colours; AGE gradient returns', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All');
    expect(applied).toBe(true); 
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        gridCC: grid.props.colorCoding,
        heightMin: grid.cell('HEIGHT', s.minHRow).color >>> 0,
        heightMax: grid.cell('HEIGHT', s.maxHRow).color >>> 0,
        ageMin: grid.cell('AGE', s.minAgeRow).color >>> 0,
        ageMax: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.gridCC).toBe('All'); 

    expect(r.heightMin).not.toBe(r.heightMax);
    expect(r.heightMin).not.toBe(WHITE);
    expect(r.heightMax).not.toBe(WHITE);

    expect(r.ageMin).not.toBe(r.ageMax);
    expect(r.ageMin).not.toBe(WHITE);
  });

  await softStep('Step 7-8 — Grid Color Coding Auto: per-column AGE coding wins; HEIGHT reverts to background', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---Auto');
    expect(applied).toBe(true); 
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        gridCC: grid.props.colorCoding,
        ageMin: grid.cell('AGE', s.minAgeRow).color >>> 0,
        ageMax: grid.cell('AGE', s.maxAgeRow).color >>> 0,
        heightMin: grid.cell('HEIGHT', s.minHRow).color >>> 0,
        heightMax: grid.cell('HEIGHT', s.maxHRow).color >>> 0,
      };
    }, setup);
    expect(r.gridCC).toBe('Auto'); 

    expect(r.ageMin).not.toBe(r.ageMax);
    expect(r.ageMin).not.toBe(WHITE);

    expect(r.heightMin).toBe(WHITE);
    expect(r.heightMax).toBe(WHITE);
  });

  await softStep('Step 9-14 — Coding-application half of the style-vs-coding order (GROK-18638): Linear coding on AGE resolves to the coding colour', async () => {

    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
      return {
        ccType: df.col('AGE').getTag('.color-coding-type'),
        minColor: grid.cell('AGE', s.minAgeRow).color >>> 0,
        maxColor: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.ccType).toBe('Linear'); 
    expect(r.minColor).not.toBe(WHITE); 
    expect(r.minColor).not.toBe(r.maxColor);
  });

  await softStep('Step 15-17 — Narrow the AGE column sharply under Grid Coding All: the resolved cell colour is unchanged', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const toAll = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All');
    expect(toAll).toBe(true);
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const before = await page.evaluate((s) => ({
      color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
      width: grok.shell.tv.grid.columns.byName('AGE').width,
      valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
    }), setup);

    const drag = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const db = grid.cell('AGE', 0).documentBounds;
      const headerY = db.y - grid.colHeaderHeight / 2;
      return {rightBorderX: db.x + db.width, y: headerY, targetX: db.x + 12};
    });
    await page.mouse.move(drag.rightBorderX, drag.y, {steps: 2});
    await page.mouse.down();
    await page.mouse.move(drag.targetX, drag.y, {steps: 12});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const after = await page.evaluate((s) => ({
      color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
      width: grok.shell.tv.grid.columns.byName('AGE').width,
      valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
    }), setup);
    expect(after.width).toBeLessThan(before.width); 
    expect(after.color).toBe(before.color); 

    expect(after.valueString).toBe(before.valueString);
  });

  await softStep('Step 18-20 — Apply a numeric format on AGE via the header menu: the format tag is set and valueString honours it', async () => {

    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').width = 60; grok.shell.tv.grid.invalidate(); });
    await v.waitForViewerRendered(page, 'Grid', 800);
    const c = await headerCenter(page, 'AGE');
    const opened = await page.evaluate(async (coord) => {
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const cm = {bubbles: true, cancelable: true, clientX: coord.x, clientY: coord.y, button: 2, buttons: 2} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', cm));
      overlay.dispatchEvent(new MouseEvent('mouseup', cm));
      overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
      return true;
    }, c);
    expect(opened).toBe(true);
    await page.waitForTimeout(550); 
    const grp = await laidOutRect(page, 'div-Format');
    if (grp) {
      await page.mouse.move(grp.x + 4, grp.y + grp.h / 2, {steps: 4});
      await page.mouse.move(grp.x + grp.w / 2, grp.y + grp.h / 2, {steps: 8});
      await page.waitForTimeout(450); 
    }
    const custom = await laidOutRect(page, 'div-Format---Custom...');
    expect(custom).not.toBeNull(); 
    await page.mouse.move(custom!.x + custom!.w / 2, custom!.y + custom!.h / 2, {steps: 4});
    await page.mouse.click(custom!.x + custom!.w / 2, custom!.y + custom!.h / 2);

    await page.locator('[name="dialog-Format-AGE"] [name="input-Custom"]').waitFor({state: 'attached', timeout: 6000});
    await page.evaluate(() => {
      const inp = document.querySelector('[name="dialog-Format-AGE"] [name="input-Custom"]') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(inp, '0.00');
      inp.dispatchEvent(new Event('input', {bubbles: true}));
      inp.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('[name="dialog-Format-AGE"] [name="button-OK"]').click();

    await page.locator('[name="dialog-Format-AGE"]').waitFor({state: 'detached', timeout: 6000});
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;

      grid.columns.byName('AGE').width = 26; grid.invalidate();
      return {
        formatTag: df.col('AGE').getTag('format'),
        valueString: grid.cell('AGE', s.targetRow).cell.valueString,
      };
    }, setup);
    expect(r.formatTag).toBe('0.00'); 
    expect(r.valueString).toMatch(/\.\d{2}$/); 

    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').width = 60; grok.shell.tv.grid.invalidate(); });
  });

  await softStep('Step 21-24a — Missing Value Color via the gear panel: a null AGE cell resolves the configured colour', async () => {
    const defaultNullColor = await page.evaluate((s) =>
      grok.shell.tv.grid.cell('AGE', s.nullAgeRow).color >>> 0, setup);
    expect(await openGridSettings(page)).toBe(true);
    await page.locator('[name="prop-missing-value-color"]').waitFor({state: 'attached', timeout: 8000});
    await page.evaluate(() => {
      const view = document.querySelector('[name="prop-view-missing-value-color"]') as HTMLElement | null;
      if (!view) return;
      const r = view.getBoundingClientRect();
      const at = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      for (const type of ['mouseover', 'mousedown', 'mouseup', 'click']) view.dispatchEvent(new MouseEvent(type, at));
    });
    await page.locator('.property-grid-item-editor-color-picker-host').waitFor({state: 'attached', timeout: 5000});
    await page.evaluate(() => {
      const host = document.querySelector('.property-grid-item-editor-color-picker-host') as HTMLElement;
      const hex = host.querySelector('input.ui-input-editor[type="text"]') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(hex, '#FFFF00');
      hex.dispatchEvent(new Event('input', {bubbles: true}));
      hex.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        nullColor: grid.cell('AGE', s.nullAgeRow).color >>> 0,
        configured: grid.props.missingValueColor >>> 0,
      };
    }, setup);
    expect(r.nullColor).toBe(r.configured); 
    expect(r.nullColor).toBe(0xffffff00); 
    expect(r.nullColor).not.toBe(defaultNullColor); 
  });

  await softStep('Step 25-27 — Font size via the gear panel produces a settle-gated canvas render delta (GROK-17767)', async () => {
    expect(await openGridSettings(page)).toBe(true);
    await page.locator('[name="prop-default-cell-font"]').waitFor({state: 'attached', timeout: 8000});

    await page.waitForTimeout(900); 
    await snapGridCanvas(page);
    await page.waitForTimeout(700); 
    const idleDelta = await diffGridCanvas(page);
    expect(idleDelta).toBeGreaterThanOrEqual(0); 
    console.log('idle grid canvas delta (px):', idleDelta);
    expect(idleDelta).toBeLessThan(2000); 
    await snapGridCanvas(page);
    const beforeFont = await page.evaluate(() => grok.shell.tv.grid.props.defaultCellFont);
    await page.evaluate(() => {
      const size = document.querySelector('[name="prop-default-cell-font"] .d4-font-size-input') as HTMLInputElement | null;
      if (size) {
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(size, '20');
        size.dispatchEvent(new Event('input', {bubbles: true}));
        size.dispatchEvent(new Event('change', {bubbles: true}));
        size.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter'}));
      }
    });
    await v.waitForViewerRendered(page, 'Grid', 1500);
    await page.waitForTimeout(400); 
    const changeDelta = await diffGridCanvas(page);
    console.log('font-change grid canvas delta (px):', changeDelta);
    expect(changeDelta).toBeGreaterThanOrEqual(0); 
    expect(changeDelta).toBeGreaterThan(3000); 
    const afterFont = await page.evaluate(() => grok.shell.tv.grid.props.defaultCellFont);
    expect(afterFont).not.toBe(beforeFont); 
  });

  const gridErrors = consoleErrors.slice(baselineErrors).filter((e) => /grid|column|index|color/i.test(e));
  expect(gridErrors).toEqual([]); 

  v.finishSpec();
});
