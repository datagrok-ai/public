/* ---
realizes: [grid.cp.cell-appearance, grid.int.color-resolution-order]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;

test.use(specTestOptions);

async function snapGridCanvas(page: Page): Promise<boolean> {
  return page.evaluate(() => {
    const cv = document.querySelector('[name="viewer-Grid"] canvas[name="canvas"]') as HTMLCanvasElement | null;
    if (!cv) return false;
    try {
      const w = window as any;
      w.__gridCanvasSnap = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      w.__gridCanvasW = cv.width; w.__gridCanvasH = cv.height;
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
      for (let i = 0; i < cur.length; i += 4)
        if (cur[i] !== prev[i] || cur[i + 1] !== prev[i + 1] || cur[i + 2] !== prev[i + 2]) delta++;
      return delta;
    } catch { return -1; }
  });
}

interface Setup { minAgeRow: number; maxAgeRow: number; nullAgeRow: number; minHRow: number; maxHRow: number; targetRow: number; }

interface Colors { gridCC: string; ageCC: string; ageMin: number; ageMax: number; heightMin: number; heightMax: number; }

function readColors(page: Page, s: Setup): Promise<Colors> {
  return page.evaluate((s) => {
    const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
    return {
      gridCC: grid.props.colorCoding,
      ageCC: df.col('AGE').getTag('.color-coding-type'),
      ageMin: grid.cell('AGE', s.minAgeRow).color >>> 0,
      ageMax: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      heightMin: grid.cell('HEIGHT', s.minHRow).color >>> 0,
      heightMax: grid.cell('HEIGHT', s.maxHRow).color >>> 0,
    };
  }, s);
}

const WHITE = 0xffffffff;

test('Grid — Cell Appearance and Color Resolution Order', async ({page}) => {
  test.setTimeout(240_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page);
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
    const baselineErrors = errors.count();

    const setup: Setup = await page.evaluate(() => {
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

    const colorsAfter = (ok: (c: Colors) => boolean) => v.pollValue(() => readColors(page, setup), ok, 1500, 50);

    await softStep('Step 1-2 — Per-column Linear on AGE via the header menu: min/max cells differ, both differ from background', async () => {
      const c = await g.headerCenter(page, 'AGE');
      expect(await g.clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Linear')).toBe(true);
      const r = await colorsAfter((x) => x.ageCC === 'Linear' && x.ageMin !== x.ageMax && x.ageMin !== WHITE && x.ageMax !== WHITE);
      expect(r.ageCC).toBe('Linear');
      expect(r.ageMin).not.toBe(r.ageMax);
      expect(r.ageMin).not.toBe(WHITE);
      expect(r.ageMax).not.toBe(WHITE);
    });

    await softStep('Step 3-4 — Grid Color Coding None: AGE cells revert to the plain background', async () => {
      const c = await g.cellCenter(page, 'AGE', 0);
      expect(await g.clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---None')).toBe(true);
      const r = await colorsAfter((x) => x.gridCC === 'None' && x.ageMin === WHITE && x.ageMax === WHITE);
      expect(r.gridCC).toBe('None');
      expect(r.ageMin).toBe(WHITE);
      expect(r.ageMax).toBe(WHITE);
    });

    await softStep('Step 5-6 — Grid Color Coding All: HEIGHT (previously uncolored) now auto-colours; AGE gradient returns', async () => {
      const c = await g.cellCenter(page, 'AGE', 0);
      expect(await g.clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All')).toBe(true);
      const r = await colorsAfter((x) => x.gridCC === 'All' && x.heightMin !== x.heightMax && x.heightMin !== WHITE &&
        x.heightMax !== WHITE && x.ageMin !== x.ageMax && x.ageMin !== WHITE);
      expect(r.gridCC).toBe('All');
      expect(r.heightMin).not.toBe(r.heightMax);
      expect(r.heightMin).not.toBe(WHITE);
      expect(r.heightMax).not.toBe(WHITE);
      expect(r.ageMin).not.toBe(r.ageMax);
      expect(r.ageMin).not.toBe(WHITE);
    });

    await softStep('Step 7-8 — Grid Color Coding Auto: per-column AGE coding wins; HEIGHT reverts to background', async () => {
      const c = await g.cellCenter(page, 'AGE', 0);
      expect(await g.clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---Auto')).toBe(true);
      const r = await colorsAfter((x) => x.gridCC === 'Auto' && x.ageMin !== x.ageMax && x.ageMin !== WHITE &&
        x.heightMin === WHITE && x.heightMax === WHITE);
      expect(r.gridCC).toBe('Auto');
      expect(r.ageMin).not.toBe(r.ageMax);
      expect(r.ageMin).not.toBe(WHITE);
      expect(r.heightMin).toBe(WHITE);
      expect(r.heightMax).toBe(WHITE);
    });

    await softStep('Step 9-14 — Coding-application half of the style-vs-coding order (GROK-18638): Linear coding on AGE resolves to the coding colour', async () => {
      const r = await readColors(page, setup);
      expect(r.ageCC).toBe('Linear');
      expect(r.ageMin).not.toBe(WHITE);
      expect(r.ageMin).not.toBe(r.ageMax);
    });

    await softStep('Step 15-17 — Narrow the AGE column sharply under Grid Coding All: the resolved cell colour is unchanged', async () => {
      const c = await g.cellCenter(page, 'AGE', 0);
      expect(await g.clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All')).toBe(true);
      expect((await colorsAfter((x) => x.gridCC === 'All')).gridCC).toBe('All');
      const readAge = () => page.evaluate((s) => ({
        color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
        width: grok.shell.tv.grid.columns.byName('AGE').width,
        valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
      }), setup);
      const before = await readAge();

      const drag = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const db = grid.cell('AGE', 0).documentBounds;
        return {rightBorderX: db.x + db.width, y: db.y - grid.colHeaderHeight / 2, targetX: db.x + 12};
      });
      await page.mouse.move(drag.rightBorderX, drag.y, {steps: 2});
      await page.mouse.down();
      await page.mouse.move(drag.targetX, drag.y, {steps: 12});
      await page.mouse.up();
      const after = await v.pollValue(readAge, (x) => x.width < before.width, 1500, 50);
      expect(after.width).toBeLessThan(before.width);
      expect(after.color).toBe(before.color);
      expect(after.valueString).toBe(before.valueString);
    });

    await softStep('Step 18-20 — Apply a numeric format on AGE via the header menu: the format tag is set and valueString honours it', async () => {
      await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').width = 60; grok.shell.tv.grid.invalidate(); });
      const c = await g.headerCenter(page, 'AGE');
      expect(await g.clickMenuLeaf(page, c, ['div-Format'], 'div-Format---Custom...')).toBe(true);

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
      expect(await g.openGridSettings(page, 'prop-color-coding')).toBe(true);
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
      const r = await v.pollValue(() => page.evaluate((s) => {
        const grid = grok.shell.tv.grid;
        return {
          nullColor: grid.cell('AGE', s.nullAgeRow).color >>> 0,
          configured: grid.props.missingValueColor >>> 0,
        };
      }, setup), (x) => x.nullColor === 0xffffff00, 1500, 50);
      expect(r.nullColor).toBe(r.configured);
      expect(r.nullColor).toBe(0xffffff00);
      expect(r.nullColor).not.toBe(defaultNullColor);
    });

    await softStep('Step 25-27 — Font size via the gear panel produces a settle-gated canvas render delta (GROK-17767)', async () => {
      expect(await g.openGridSettings(page, 'prop-color-coding')).toBe(true);
      await page.locator('[name="prop-default-cell-font"]').waitFor({state: 'attached', timeout: 8000});

      await v.waitForGridPainted(page, {gapMs: 250, capMs: 900});
      await snapGridCanvas(page);
      // the idle window IS the assertion: nothing may repaint the canvas while no property changes
      await page.waitForTimeout(700);
      const idleDelta = await diffGridCanvas(page);
      expect(idleDelta).toBeGreaterThanOrEqual(0);
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
      const changeDelta = await v.pollValue(() => diffGridCanvas(page), (d) => d > 3000, 1900, 100);
      expect(changeDelta).toBeGreaterThanOrEqual(0);
      expect(changeDelta).toBeGreaterThan(3000);
      const afterFont = await page.evaluate(() => grok.shell.tv.grid.props.defaultCellFont);
      expect(afterFont).not.toBe(beforeFont);
    });

    const gridErrors = errors.list.slice(baselineErrors).filter((e) => /grid|column|index|color/i.test(e));
    expect(gridErrors).toEqual([]);
  } finally {
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});
