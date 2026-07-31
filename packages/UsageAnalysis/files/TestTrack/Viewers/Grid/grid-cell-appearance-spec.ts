/* ---
realizes: [grid.cp.cell-appearance, grid.int.color-resolution-order]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;


test.use(specTestOptions);

// Page-coordinate center of a column header from the grid geometry.
// X is derived from the first data cell's documentBounds (true page coords),
// NOT from `overlay.rect.x + gridColumn.left`: gridColumn.left is an offset in the
// grid's virtual coordinate space and does NOT align with the overlay canvas'
// page-left, so `rc.x + gc.left` lands on a DIFFERENT column
//: AGE gc.left=189 → rc.x+gc.left+w/2=523 hitTests RACE, while the AGE
// cell documentBounds.x=348 → center 378 hitTests AGE). documentBounds returns page
// coords directly and mirrors cellCenter, which was always correct.
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

// Rendered rect (page coords) of a menu element whose `name` matches, filtered to
// the actually-laid-out copy: a d4 nested submenu keeps a detached zero-rect
// TEMPLATE copy of every leaf inside the parent group's
// `.d4-menu-item-container-fixed` (display:none until opened) AND, once opened, a
// laid-out copy in a separate `.d4-menu-popup`; querying by name alone hits the
// template first, so filter to the copy with a real bounding box + non-null
// offsetParent. Returns null when no laid-out copy exists yet.
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

// Open a grid context menu at (clientX, clientY) on the overlay canvas, then walk
// a chain of nested menu groups to one of their leaves and click it.
//
// The nested vertical submenus (Color-Coding, Grid-Color-Coding) use slope-based
// hover protection: a leaf stays a zero-rect detached template until its parent
// group receives a TRAJECTORY-BEARING, TRUSTED pointer movement. Synthetic
// MouseEvent dispatch never satisfies the slope tracker, so the
// submenu never renders and the leaf click is silently dropped (recon-verified
// live. Real trusted input via page.mouse.move with steps DOES lay the
// leaves out, and a trusted page.mouse.click on the rendered leaf actuates the
// real menu command — the effect under test is produced ONLY by real user input.
//
// The ROOT context menu is opened synthetically on the overlay canvas (that path
// is trusted-input-independent); only the submenu expansions need trusted input.
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

// Open the grid's property panel and wait for its rows to exist. The gear's CSS
// visibility flickers with viewer hover/focus, so a plain Playwright click times
// out on actionability. Trusted hover-then-click on the gear is tried first, then
// the viewer-scoped settings icon, then F4 on the focused grid. No JS-API
// substitution — every gesture drives the real settings-panel open.
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

// Per-color snapshot + diff of the grid's DATA canvas (canvas[name="canvas"] — the
// one that renders cell text). The viewers-helper snapshot reads root.querySelector
// ('canvas') which is the grid's first (unnamed) canvas, not the data canvas, so a
// grid-specific reader is used for the font-size render delta.
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

  // Setup: baseline console-error count + the min/max/null AGE rows and the min/max
  // HEIGHT rows used throughout for per-cell colour reads. grid.cell(col, tableRowIdx)
  //.color resolves the per-cell rendered colour directly (recon-verified — no scroll
  // needed). White background is 0xffffffff.
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
  // The core of grid.color-resolution-order. Per-column Linear on AGE is applied
  // through the AGE header context menu (trusted submenu driver). Then Grid Color
  // Coding None / All / Auto is switched through the CELL context menu's Grid-Color-
  // Coding submenu, and the resolved per-cell colours are read at each transition.

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
  // The explicit content-style colour is set through the AGE column's Column
  // Properties > Style > Content colour picker — a Dart-internal editor that does
  // not render under headless Playwright (recon-verified: the Style >
  // Content editors do not lay out; no readable tag surfaces the explicit style
  // colour). The "coding beats explicit style" invariant therefore cannot be driven
  // headless in the presence of a set explicit style, so it is waived
  // (gesture-uncontrollable-headless) rather than falsely asserted. What IS driven
  // here is the coding-application half through the real header menu.

  await softStep('Step 9-14 — Coding-application half of the style-vs-coding order (GROK-18638): Linear coding on AGE resolves to the coding colour', async () => {
    // AGE already carries Linear from Scenario 1; confirm the resolved colour is the
    // coding colour (not the plain background), which is the readable half of the
    // resolution-order invariant. The explicit-style precondition is un-drivable
    // headless and is waived below.
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
    // Put AGE under grid-wide All auto-colouring, record the resolved colour, then
    // shrink the column with a REAL trusted resize drag on the header right border,
    // and re-read the colour at the SAME row.
    const c = await cellCenter(page, 'AGE', 0);
    const toAll = await clickMenuLeaf(page, c, ['div-Grid-Color-Coding'], 'div-Grid-Color-Coding---All');
    expect(toAll).toBe(true);
    await page.waitForTimeout(500);
    const before = await page.evaluate((s) => ({
      color: grok.shell.tv.grid.cell('AGE', s.targetRow).color >>> 0,
      width: grok.shell.tv.grid.columns.byName('AGE').width,
      valueString: grok.shell.tv.grid.cell('AGE', s.targetRow).cell.valueString,
    }), setup);
    // Trusted resize drag: press at the AGE header right border, release far left.
    // Geometry is derived from the AGE cell's documentBounds (page coords) — the
    // right border is db.x + db.width; targetX pulls it far left. Using
    // `overlay.rect.x + gridColumn.left` would land on the WRONG column's border
    // (recon-verified: rc.x+gc.left+width=553 hitTests DIS_POP, while the
    // AGE cell right border db.x+db.width=408 is the real AGE resize handle).
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
    // Restore AGE width so the header menu geometry is clean, then apply a 2-decimal
    // format through the header context menu's Format submenu. AGE is an int column,
    // whose preset list carries no plain two-decimal leaf (only int / money / compact /
    // percent / thousand-separator / full-precision / money(<sym>) variants — the money
    // presets embed ".00", which is why a text-match on ".00" wrongly picks "money").
    // The two-decimal format is applied through Format > Custom..., typing 0.00 into
    // the dialog's Custom field — the real menu-plus-dialog actuation of record.
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
    // Trusted glide onto the Format group to lay out its leaves, then click Custom...
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
    // Fill the Format dialog's Custom field with 0.00 and confirm — the format is applied
    // ONLY through the real dialog, never via a JS-API meta.format write.
    // dev.datagrok.ai server 1.28.0):
    //   [name="dialog-Format-AGE"]  — dialog root; name pattern is dialog-Format-<COLNAME> [DOM]
    //   [name="input-Custom"]       — free-text format-pattern <INPUT type="text"> [DOM]
    //   [name="button-OK"]          — applies format, closes dialog; writes df col 'format' tag [DOM]
    // Reference:.claude/skills/grok-browser/references/viewers/grid.md §Format > Custom... dialog
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
      // Narrow AGE so the displayed text loses characters, then confirm valueString
      // is still the full-precision formatted value.
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
    // Settle-precheck: snapshot, let late renders drain, and assert the residual
    // idle delta is small so a stale repaint cannot satisfy the font-change ceiling.
    await page.waitForTimeout(900);
    await snapGridCanvas(page);
    await page.waitForTimeout(700);
    const idleDelta = await diffGridCanvas(page);
    expect(idleDelta).toBeGreaterThanOrEqual(0); // fault guard before the ceiling
    console.log('idle grid canvas delta (px):', idleDelta);
    expect(idleDelta).toBeLessThan(2000); // threshold validated on dev
    // Baseline for the real change, then increase the default cell font size.
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

  // The Selected Rows Colour case is owned by grid-ui.md (a manual [TESTED] section):
  // the only mouse gesture that selects rows is a strip drag that is refuted for
  // automation, and the selection overlay is a canvas paint that grid.cell.color does
  // not surface. Asserting a row-selection driven by that gesture here would be a
  // false-RED (assert-then-can't-drive), so it stays in the companion. Do not re-add it.

  // Whole-scenario console-error floor: no console error was raised across the run.
  const gridErrors = consoleErrors.slice(baselineErrors).filter((e) => /grid|column|index|color/i.test(e));
  expect(gridErrors).toEqual([]); // the appearance flows raised no grid console error

  v.finishSpec();
});
