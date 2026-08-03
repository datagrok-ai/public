/* ---
realizes: [grid.cp.cell-appearance, grid.int.color-resolution-order]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;


test.use(specTestOptions);

// Page-coordinate center of a column header. X comes from the first data cell's documentBounds,
// which is already in page coordinates; `overlay.rect.x + gridColumn.left` would land on a
// DIFFERENT column, because gridColumn.left is an offset in the grid's virtual coordinate space
// and does not align with the overlay canvas' page-left.
async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, 0).documentBounds;
    const headerY = db.y - grid.colHeaderHeight / 2;
    return {x: db.x + db.width / 2, y: headerY};
  }, col);
}

// Page-coordinate center of a data cell on the overlay canvas.
async function cellCenter(page: Page, col: string, row: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, r}) => {
    const db = grok.shell.tv.grid.cell(c, r).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, r: row});
}

// Rendered rect (page coords) of a menu element whose `name` matches, or null while no such
// copy exists. A d4 nested submenu keeps a detached zero-rect TEMPLATE copy of every leaf
// alongside the laid-out copy in the open popup, and a query by name hits the template first —
// hence the filter on a non-null offsetParent plus a real bounding box.
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

// Open a grid context menu at (clientX, clientY) on the overlay canvas, then walk a chain of
// nested menu groups to one of their leaves and click it.
//
// The nested submenus (Color-Coding, Grid-Color-Coding) use slope-based hover protection: a leaf
// stays a zero-rect detached template until its parent group receives a TRAJECTORY-BEARING,
// TRUSTED pointer movement, which no synthetic MouseEvent chain satisfies. Only the submenu
// expansions need that trusted input — the ROOT menu still opens on synthetic events.
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

// Open the grid's property panel and wait for its rows to exist. The gear's CSS visibility
// flickers with viewer hover/focus, so a plain Playwright click times out on actionability; four
// real gestures are tried in order and the first that reveals the rows wins.
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
      // panel did not appear — fall through to the next gesture
    }
  }
  return false;
}

// Snapshot of the grid's DATA canvas (canvas[name="canvas"], the one that renders cell text).
// The viewers-helper snapshot reads the grid's first, unnamed canvas instead, which carries no
// cell text and so shows no font-size delta.
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

// Summed per-pixel RGB delta of the grid data canvas vs the stored snapshot.
// Returns -1 on fault (no canvas / size mismatch / getImageData failure); the
// caller guards -1 before any ceiling assert.
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

  // Setup: baseline console-error count plus the min/max/null AGE and min/max HEIGHT rows used
  // for the per-cell colour reads. grid.cell(col, tableRowIdx).color resolves the rendered
  // colour without scrolling; the plain background is 0xffffffff.
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
    // A general non-null AGE target row for the format / resize reads.
    let targetRow = -1;
    for (let i = 0; i < df.rowCount; i++) if (!ac.isNone(i)) { targetRow = i; break; }
    return {minAgeRow, maxAgeRow, nullAgeRow, minHRow, maxHRow, targetRow};
  });
  expect(setup.minAgeRow).not.toBe(setup.maxAgeRow); // distinct min/max AGE rows
  expect(setup.nullAgeRow).toBeGreaterThanOrEqual(0); // a null AGE cell exists
  expect(setup.minHRow).not.toBe(setup.maxHRow); // distinct min/max HEIGHT rows

  // --- Scenario 1: Grid-wide colour coding overrides per-column coding ----------
  // The core of grid.color-resolution-order: per-column Linear is applied from the AGE header
  // menu, then None / All / Auto from the cell menu, reading the resolved colours at each step.

  await softStep('Step 1-2 — Per-column Linear on AGE via the header menu: min/max cells differ, both differ from background', async () => {
    const c = await headerCenter(page, 'AGE');
    const applied = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Linear');
    expect(applied).toBe(true); // the Linear leaf was reached and clicked via the header menu
    await page.waitForTimeout(600);
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
      return {
        ccType: df.col('AGE').getTag('.color-coding-type'),
        minColor: grid.cell('AGE', s.minAgeRow).color >>> 0,
        maxColor: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.ccType).toBe('Linear'); // Linear per-column coding is active on AGE
    expect(r.minColor).not.toBe(r.maxColor); // min-AGE and max-AGE cells differ
    expect(r.minColor).not.toBe(WHITE); // and both differ from the plain background
    expect(r.maxColor).not.toBe(WHITE);
  });

  await softStep('Step 3-4 — Grid Color Coding None: AGE cells revert to the plain background', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---None');
    expect(applied).toBe(true); // the None radio was reached and clicked
    await page.waitForTimeout(600);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        gridCC: grid.props.colorCoding,
        ageMin: grid.cell('AGE', s.minAgeRow).color >>> 0,
        ageMax: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.gridCC).toBe('None'); // grid-wide None is set
    // The per-column Linear coding on AGE is suppressed by the grid-wide None.
    expect(r.ageMin).toBe(WHITE);
    expect(r.ageMax).toBe(WHITE);
  });

  await softStep('Step 5-6 — Grid Color Coding All: HEIGHT (previously uncolored) now auto-colours; AGE gradient returns', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All');
    expect(applied).toBe(true); // the All radio was reached and clicked
    await page.waitForTimeout(600);
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
    expect(r.gridCC).toBe('All'); // grid-wide All is set
    // HEIGHT had no per-column coding — grid-wide All auto-colours it.
    expect(r.heightMin).not.toBe(r.heightMax);
    expect(r.heightMin).not.toBe(WHITE);
    expect(r.heightMax).not.toBe(WHITE);
    // AGE gradient is present again too.
    expect(r.ageMin).not.toBe(r.ageMax);
    expect(r.ageMin).not.toBe(WHITE);
  });

  await softStep('Step 7-8 — Grid Color Coding Auto: per-column AGE coding wins; HEIGHT reverts to background', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const applied = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---Auto');
    expect(applied).toBe(true); // the Auto radio was reached and clicked
    await page.waitForTimeout(600);
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
    expect(r.gridCC).toBe('Auto'); // grid-wide Auto is set
    // Explicit per-column Linear coding on AGE wins under Auto.
    expect(r.ageMin).not.toBe(r.ageMax);
    expect(r.ageMin).not.toBe(WHITE);
    // HEIGHT has no per-column coding — Auto does not force it, so it reverts.
    expect(r.heightMin).toBe(WHITE);
    expect(r.heightMax).toBe(WHITE);
  });

  // --- Scenario 2: Column coding overrides explicit style colour (GROK-18638) ---
  // The explicit content-style colour lives in Column Properties > Style > Content, a
  // Dart-internal editor that does not lay out under headless Playwright and surfaces no readable
  // tag. The "coding beats explicit style" half is therefore waived
  // (gesture-uncontrollable-headless); only the coding-application half is driven here.

  await softStep('Step 9-14 — Coding-application half of the style-vs-coding order (GROK-18638): Linear coding on AGE resolves to the coding colour', async () => {
    // AGE already carries Linear from Scenario 1.
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
      return {
        ccType: df.col('AGE').getTag('.color-coding-type'),
        minColor: grid.cell('AGE', s.minAgeRow).color >>> 0,
        maxColor: grid.cell('AGE', s.maxAgeRow).color >>> 0,
      };
    }, setup);
    expect(r.ccType).toBe('Linear'); // per-column coding is enabled on AGE
    expect(r.minColor).not.toBe(WHITE); // and the resolved cell colour is the coding colour
    expect(r.minColor).not.toBe(r.maxColor);
  });

  // --- Scenario 3: Adaptive rendering does not change resolved colour (GROK-19113)

  await softStep('Step 15-17 — Narrow the AGE column sharply under Grid Coding All: the resolved cell colour is unchanged', async () => {
    const c = await cellCenter(page, 'AGE', 0);
    const toAll = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All');
    expect(toAll).toBe(true);
    await page.waitForTimeout(500);
    const before = await page.evaluate((s) => ({
      color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
      width: grok.shell.tv.grid.columns.byName('AGE').width,
      valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
    }), setup);
    // Trusted resize drag from the AGE header right border. The geometry comes from the cell's
    // documentBounds; `overlay.rect.x + gridColumn.left` would grab another column's border.
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
    await page.waitForTimeout(700);
    const after = await page.evaluate((s) => ({
      color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
      width: grok.shell.tv.grid.columns.byName('AGE').width,
      valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
    }), setup);
    expect(after.width).toBeLessThan(before.width); // the column was actually narrowed
    expect(after.color).toBe(before.color); // GROK-19113: the resolved cell colour is unchanged by adaptive shrinking
    // and the full-precision value is still reported despite the narrow display.
    expect(after.valueString).toBe(before.valueString);
  });

  // --- Scenario 4: Column format tag propagates to valueString ------------------

  await softStep('Step 18-20 — Apply a numeric format on AGE via the header menu: the format tag is set and valueString honours it', async () => {
    // AGE is an int column, and its preset list carries no plain two-decimal leaf (the money
    // presets merely embed ".00"), so the format goes through Format > Custom... instead.
    // The width is restored first so the header menu geometry is clean.
    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').width = 60; grok.shell.tv.grid.invalidate(); });
    await page.waitForTimeout(300);
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
    expect(custom).not.toBeNull(); // the Format > Custom... leaf is laid out
    await page.mouse.move(custom!.x + custom!.w / 2, custom!.y + custom!.h / 2, {steps: 4});
    await page.mouse.click(custom!.x + custom!.w / 2, custom!.y + custom!.h / 2);
    await page.waitForTimeout(700);
    // The format is applied only through the real dialog, never a JS-API meta.format write.
    // The dialog root's name follows the pattern dialog-Format-<COLNAME>; see
    // .claude/skills/grok-browser/references/viewers/grid.md §Format > Custom... dialog.
    await page.locator('[name="dialog-Format-AGE"] [name="input-Custom"]').waitFor({state: 'attached', timeout: 6000});
    await page.evaluate(() => {
      const inp = document.querySelector('[name="dialog-Format-AGE"] [name="input-Custom"]') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(inp, '0.00');
      inp.dispatchEvent(new Event('input', {bubbles: true}));
      inp.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('[name="dialog-Format-AGE"] [name="button-OK"]').click();
    await page.waitForTimeout(700);
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
      // Narrow AGE so the displayed text loses characters — valueString must stay full-precision.
      grid.columns.byName('AGE').width = 26; grid.invalidate();
      return {
        formatTag: df.col('AGE').getTag('format'),
        valueString: grid.cell('AGE', s.targetRow).cell.valueString,
      };
    }, setup);
    expect(r.formatTag).toBe('0.00'); // the chosen numeric format is stored on the df column
    expect(r.valueString).toMatch(/\.\d{2}$/); // valueString honours the format (2 decimals)
    // restore width for later geometry
    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').width = 60; grok.shell.tv.grid.invalidate(); });
  });

  // --- Scenario 5: Style settings apply per cell read (GROK-19813) --------------

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
    await page.waitForTimeout(800);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        nullColor: grid.cell('AGE', s.nullAgeRow).color >>> 0,
        configured: grid.props.missingValueColor >>> 0,
      };
    }, setup);
    expect(r.nullColor).toBe(r.configured); // the null AGE cell resolves the configured missing-value colour
    expect(r.nullColor).toBe(0xffffff00); // which is the #FFFF00 that was set
    expect(r.nullColor).not.toBe(defaultNullColor); // and differs from the default background
  });

  await softStep('Step 25-27 — Font size via the gear panel produces a settle-gated canvas render delta (GROK-17767)', async () => {
    expect(await openGridSettings(page)).toBe(true);
    await page.locator('[name="prop-default-cell-font"]').waitFor({state: 'attached', timeout: 8000});
    // Settle-precheck: a small residual idle delta is required so that a stale repaint cannot
    // satisfy the font-change threshold on its own.
    await page.waitForTimeout(900);
    await snapGridCanvas(page);
    await page.waitForTimeout(700);
    const idleDelta = await diffGridCanvas(page);
    expect(idleDelta).toBeGreaterThanOrEqual(0); // fault guard before the ceiling
    console.log('idle grid canvas delta (px):', idleDelta);
    expect(idleDelta).toBeLessThan(2000); // threshold validated on dev
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
    await page.waitForTimeout(1200);
    const changeDelta = await diffGridCanvas(page);
    console.log('font-change grid canvas delta (px):', changeDelta);
    expect(changeDelta).toBeGreaterThanOrEqual(0); // fault guard before the ceiling
    expect(changeDelta).toBeGreaterThan(3000); // GROK-17767: the larger font re-renders the cell text — threshold validated on dev
    const afterFont = await page.evaluate(() => grok.shell.tv.grid.props.defaultCellFont);
    expect(afterFont).not.toBe(beforeFont); // and the default cell font actually changed
  });

  // The Selected Rows Colour case stays in the manual companion grid-ui.md: the selection
  // overlay is a canvas paint that grid.cell.color does not surface. Do not re-add it here.

  const gridErrors = consoleErrors.slice(baselineErrors).filter((e) => /grid|column|index|color/i.test(e));
  expect(gridErrors).toEqual([]); // the appearance flows raised no grid console error

  v.finishSpec();
});
