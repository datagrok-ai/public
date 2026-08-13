/* ---
realizes: [trellisplot.cp.scroll-categories, trellisplot.int.pack-categories-vs-viewport-scroll]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Ambient noise the packed dev build emits regardless of what is actuated. Nothing
// product-specific is filtered, so a NullError or a Dart stack trace still counts.
const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

// The X '+'/'-' icons carry no aria/disabled attribute — disabled is a CSS class setting
// `pointer-events: none` (trellis_plot.css:125), so enablement is read off computed style.
// Opacity carries no signal: the icon container sits at opacity 0 until CSS :hover.
async function xAxisState(page: Page): Promise<{
  cells: number; plusEnabled: boolean; minusEnabled: boolean;
}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const icons = root.querySelector('[name="x-axis-icons"]');
    const plus = icons?.querySelector('[name="icon-plus"]') as HTMLElement | null;
    const minus = icons?.querySelector('[name="icon-minus"]') as HTMLElement | null;
    const enabled = (el: HTMLElement | null) =>
      !!el && getComputedStyle(el).pointerEvents !== 'none';
    return {
      cells: root.querySelectorAll('.d4-trellis-plot-cell').length,
      plusEnabled: enabled(plus),
      minusEnabled: enabled(minus),
    };
  });
}

// Paging workhorse. A DOM click() fires the handler whatever the computed pointer-events are,
// acceptable HERE only because every call site drives an icon already asserted enabled. It
// cannot show that a DISABLED icon does nothing — that is realClickXIcon's job.
async function clickXIcon(page: Page, which: 'plus' | 'minus'): Promise<void> {
  await page.evaluate((w) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    (root.querySelector(`[name="x-axis-icons"] [name="icon-${w}"]`) as HTMLElement)?.click();
  }, which);
  await page.waitForTimeout(800);
}

// Real user-channel click: Playwright hit-tests the element first, so an icon whose disabled
// class sets `pointer-events: none` never receives it and the call times out. Returns whether
// the click landed, so the extremes are graded on BEHAVIOURAL inertness, not on a CSS read.
async function realClickXIcon(page: Page, which: 'plus' | 'minus', timeout: number): Promise<boolean> {
  try {
    await page.locator(`[name="viewer-Trellis-plot"] [name="x-axis-icons"] [name="icon-${which}"]`)
      .click({timeout});
  }
  catch {
    return false;
  }
  await page.waitForTimeout(800);
  return true;
}

// Menu items carry no name=, so "Reset X columns" is matched by label under two rules: scoped to
// `.d4-menu-popup` (unscoped also returns the MAIN menu), and gated on a non-zero box (a collapsed
// child is clickable through display:none). Refdoc trellis_plot.md pitfall 22 [DOM 2026-08-11, 2026-08-12].
async function resetXColumnsViaMenu(page: Page): Promise<boolean> {
  const clicked = await page.evaluate(async () => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const combos = Array.from(root.querySelectorAll('[name="div-column-combobox-"]'));
    // X selectors are horizontal (width > height); Y selectors are vertical.
    const xCombo = combos.find((c) => {
      const r = c.getBoundingClientRect();
      return r.width > r.height;
    }) as HTMLElement | undefined;
    if (!xCombo) return false;
    const r = xCombo.getBoundingClientRect();
    xCombo.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 2,
    }));
    await new Promise((res) => setTimeout(res, 600));
    const label = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .find((el) => {
        if ((el.textContent || '').trim() !== 'Reset X columns') return false;
        const box = el.getBoundingClientRect();
        return box.width > 0 && box.height > 0;
      });
    if (!label) { document.body.click(); return false; }
    (label.closest('.d4-menu-item') as HTMLElement).click();
    return true;
  });
  await page.waitForTimeout(2000);
  return clicked;
}

test('Trellis plot: category viewport paging (+/- icons and pack-categories coupling)', async ({page}) => {
  test.setTimeout(600_000);
  page.setDefaultTimeout(120_000);

  // Attached before any action, so the floors assert against an already-listening collector.
  // grok.shell.warnings is undefined on this build [DOM 2026-08-06], so pageerror/console is the
  // substitute. Each floor is a DELTA: a global empty-collector check would be false-red.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);
  await page.waitForTimeout(5000);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // DIS_POP = 6, RACE = 4 on the live DemoFiles demog [DOM 2026-08-06].
  const cardinalities = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {disPop: df.col('DIS_POP').categories.length, race: df.col('RACE').categories.length};
  });
  expect(cardinalities).toEqual({disPop: 6, race: 4});

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

  await page.evaluate(async () => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['DIS_POP'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.packCategories = true;
    await new Promise((r) => setTimeout(r, 2000));
  });

  // Entry state, against the scenario's "default shows a subset" assumption: all 6
  // DIS_POP categories already fit, so '+' starts DISABLED and '-' enabled
  // [DOM 2026-08-06]. yCategoriesCount = 4 is the unit every paging step moves by.
  const yCatCount = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.yCategoriesCount as number;
  });
  expect(yCatCount).toBe(4);

  // #### Scenario 1 Step 3 (entry state): read BEFORE the first click, so the preamble's
  // "'+' starts disabled, '-' enabled" is a graded fact and not prose the code assumes.
  const entry = await xAxisState(page);
  expect(entry.plusEnabled).toBe(false);
  expect(entry.minusEnabled).toBe(true);

  // Opening edge of the window the Step 7 no-error floor covers: the paging in Steps
  // 4-6 only, nothing login or dataset load emitted before it.
  const pagingWindow = {pageErrors: pageErrors.length, consoleErrors: consoleErrors.length};

  // #### Scenario 1 Step 4: '+' pages one X-category row of cells into the viewport
  await softStep('Scenario 1 Step 4', async () => {
    // Shrink by one first, so '+' has room to page a category back in.
    await clickXIcon(page, 'minus');
    const paged = await xAxisState(page);
    expect(paged.plusEnabled).toBe(true);
    const cellsBeforePlus = paged.cells;

    await clickXIcon(page, 'plus');
    const afterPlus = await xAxisState(page);
    expect(afterPlus.cells).toBe(cellsBeforePlus + yCatCount);
  });

  // #### Scenario 1 Step 5: '-' pages one X-category row back out (DOM round-trip)
  await softStep('Scenario 1 Step 5', async () => {
    // Round-trip probe: shrink first so '+' is enabled, then assert '-' is symmetric
    // with '+' rather than merely landing somewhere smaller.
    await clickXIcon(page, 'minus');
    const preExpansion = (await xAxisState(page)).cells;
    await clickXIcon(page, 'plus');
    const expanded = (await xAxisState(page)).cells;
    expect(expanded).toBe(preExpansion + yCatCount);
    await clickXIcon(page, 'minus');
    const afterMinus = (await xAxisState(page)).cells;
    expect(afterMinus).toBe(preExpansion);
  });

  // #### Scenario 1 Step 6: icon disabled styling at both viewport extremes
  await softStep('Scenario 1 Step 6', async () => {
    // 8 iterations comfortably exceed the 6 DIS_POP categories, so the loop always
    // reaches the minimum viewport.
    for (let i = 0; i < 8; i++) {
      const s = await xAxisState(page);
      if (!s.minusEnabled) break;
      await clickXIcon(page, 'minus');
    }
    const atMin = await xAxisState(page);
    expect(atMin.cells).toBe(yCatCount); // 1 X category x 4 RACE
    expect(atMin.minusEnabled).toBe(false);
    expect(atMin.plusEnabled).toBe(true);

    // Behavioural inertness at the floor, as a DIFFERENTIAL over one channel: the real
    // click must FAIL on the disabled '-' and LAND on the enabled '+' beside it. Without
    // the second half, a universally intercepted click would "prove" the first for free.
    expect(await realClickXIcon(page, 'minus', 5000)).toBe(false);
    expect((await xAxisState(page)).cells).toBe(atMin.cells);
    expect(await realClickXIcon(page, 'plus', 20000)).toBe(true);
    expect((await xAxisState(page)).cells).toBe(atMin.cells + yCatCount);

    // Grow back to the ceiling, same bound.
    for (let i = 0; i < 8; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const atMax = await xAxisState(page);
    expect(atMax.cells).toBe(cardinalities.disPop * yCatCount); // 6 x 4 = 24
    expect(atMax.plusEnabled).toBe(false);
    expect(atMax.minusEnabled).toBe(true);

    // Mirror differential at the ceiling: '+' is the inert one here, '-' the live one.
    expect(await realClickXIcon(page, 'plus', 5000)).toBe(false);
    expect((await xAxisState(page)).cells).toBe(atMax.cells);
    expect(await realClickXIcon(page, 'minus', 20000)).toBe(true);
    expect((await xAxisState(page)).cells).toBe(atMax.cells - yCatCount);
  });

  // #### Scenario 1 Step 7: no uncaught console error during the paging in steps 4-6
  await softStep('Scenario 1 Step 7', async () => {
    // Slices, not whole collectors: this floor owns the Steps 4-6 window only.
    expect(pageErrors.slice(pagingWindow.pageErrors)).toEqual([]);
    expect(consoleErrors.slice(pagingWindow.consoleErrors)).toEqual([]);
  });

  // #### Scenario 1 Step 8: resetting the X column selection via its context menu while paged
  await softStep('Scenario 1 Step 8', async () => {
    // 12 SEX x DIS_POP combinations exceed the viewport clamp, so the reset is exercised
    // from a genuinely paged state — and driven through the real context menu, not a
    // prop write, since the menu path is what the step guards. [DOM 2026-08-06]
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX', 'DIS_POP'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 2000));
    });
    const pagedCells = (await xAxisState(page)).cells;
    const pagedXCat = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return tp.xCategoriesCount as number;
    });
    expect(pagedCells).toBeLessThan(pagedXCat * yCatCount);

    const errorsBefore = pageErrors.length + consoleErrors.length;
    const reset = await resetXColumnsViaMenu(page);
    expect(reset).toBe(true);

    const after = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      return {
        xNames: tp.props.xColumnNames as string[],
        xCat: tp.xCategoriesCount as number,
        cells: root.querySelectorAll('.d4-trellis-plot-cell').length,
      };
    });
    // Reset clears the X column assignment -> one category row of yCategories cells.
    expect(after.xNames).toEqual([]);
    expect(after.cells).toBe(yCatCount);
    const errorsAfter = pageErrors.length + consoleErrors.length;
    // Resetting while paged is the crash-prone path this floor watches.
    expect(errorsAfter).toBe(errorsBefore);

    // Restore X = DIS_POP for downstream steps' independence.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['DIS_POP'];
      await new Promise((r) => setTimeout(r, 1200));
    });
  });

  // Carried from Scenario 2 Step 4 into Step 6: the packed and unpacked cell counts of the
  // SAME fully opened viewport are the two ends of the delta the scenario promises.
  let packedMaxCells = -1;
  let unpackedMaxCells = -1;

  // #### Scenario 2 Step 4: Pack Categories ON collapses the structurally empty X columns
  await softStep('Scenario 2 Step 4', async () => {
    // The X pair must carry structurally empty combinations, or packing has nothing to collapse
    // and every assertion below holds even when packCategories does nothing. The emptiness is
    // derived from the live frame, so a demog that ever fills every combination fails loudly.
    const truth = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['RACE', 'SEVERITY'];
      tp.props.yColumnNames = ['SEX'];
      tp.props.packCategories = true;
      await new Promise((r) => setTimeout(r, 2500));
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE'), severity = df.col('SEVERITY');
      const seen = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) seen.add(`${race.get(i)}|${severity.get(i)}`);
      return {
        full: race.categories.length * severity.categories.length,
        populated: seen.size,
        sexCats: df.col('SEX').categories.length as number,
        xCat: tp.xCategoriesCount as number,
        yCat: tp.yCategoriesCount as number,
      };
    });
    // xCategoriesCount is the raw category product on both sides of the toggle
    // (grok_api.dart:539, which packing never touches), so it is blind to packing — the
    // rendered cell count is the only channel that can witness it.
    expect(truth.xCat).toBe(truth.full);
    expect(truth.yCat).toBe(truth.sexCats);
    expect(truth.populated).toBeLessThan(truth.full);

    const paged = await xAxisState(page);
    expect(paged.plusEnabled).toBe(true);
    expect(paged.cells % truth.yCat).toBe(0); // whole rows of yCategories
    await clickXIcon(page, 'plus');
    const afterPlus = await xAxisState(page);
    // Paging still advances exactly one X-combination column while packing is on.
    expect(afterPlus.cells).toBe(paged.cells + truth.yCat);

    // Only at the fully opened viewport does the rendered count equal
    // (visible X categories) x (Y categories) — the state where packing is measurable.
    for (let i = 0; i < 30; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const atMax = await xAxisState(page);
    expect(atMax.plusEnabled).toBe(false);
    expect(atMax.cells % truth.yCat).toBe(0);
    expect(atMax.cells).toBeGreaterThan(0);
    // The populated-pair count bounds the packed grid from ABOVE — an upper bound, not an
    // equality, because packing additionally drops combinations whose inner-viewer value
    // column is empty on every one of their rows (trellis_plot_core.dart:185-190).
    expect(atMax.cells).toBeLessThanOrEqual(truth.populated * truth.yCat);
    packedMaxCells = atMax.cells;
    unpackedMaxCells = truth.full * truth.yCat;
  });

  // #### Scenario 2 Step 6: Pack Categories OFF restores the empty-combination cells
  await softStep('Scenario 2 Step 6', async () => {
    // Guards the carry-over: without Step 4's measurement the comparison below would be
    // against a sentinel and would pass on nothing.
    expect(packedMaxCells).toBeGreaterThan(0);
    const errorsBefore = pageErrors.length + consoleErrors.length;

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.packCategories = false;
      await new Promise((r) => setTimeout(r, 2500));
    });
    // The toggle resets the category scroll window, so the viewport is opened again
    // before the two counts are compared at the same viewport state.
    for (let i = 0; i < 30; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const offMax = await xAxisState(page);
    expect(offMax.plusEnabled).toBe(false);
    // Nothing is collapsed now: the whole category product occupies the grid.
    expect(offMax.cells).toBe(unpackedMaxCells);
    // The difference against Step 4 IS the pack-categories effect, measured on the same
    // fully opened viewport.
    expect(offMax.cells).toBeGreaterThan(packedMaxCells);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.packCategories = true;
      await new Promise((r) => setTimeout(r, 1500));
    });
    const errorsAfter = pageErrors.length + consoleErrors.length;
    expect(errorsAfter).toBe(errorsBefore);
  });

  v.finishSpec();
});
