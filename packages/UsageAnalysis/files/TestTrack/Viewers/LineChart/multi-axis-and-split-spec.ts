/* ---
realizes: [linechart.cp.multi-axis-and-split]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import {countCanvasPixels} from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await page.evaluate((p) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    for (const [k, val] of Object.entries(p)) (lc.props as any)[k] = val;
  }, props);
  await page.waitForTimeout(500);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function chartCanvasNonEmpty(page: Page): Promise<boolean> {
  // Content-level check (github-2904): a blank chart still has a non-zero
  // bounding rect, so measure drawn pixels instead of geometry. The threshold
  // sits ABOVE the axes/labels-only floor so an axes-only "blank" state fails it.
  // -1 (no canvas / getImageData fault) fails the threshold too.
  return (await countCanvasPixels(page, 'Line chart')).total > 38000;
}

test('Line Chart — Multi-Axis and Split', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X', 'Chemical Space Y']});

  // The 38000-px blank threshold in chartCanvasNonEmpty assumes the canvas size the
  // fixed viewport produces — fail loudly if a layout or viewport change resizes
  // it, instead of silently devaluing the threshold.
  const cvDims = await page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
    const cv = lc?.root?.querySelector('canvas') as HTMLCanvasElement | null;
    return cv ? {w: cv.width, h: cv.height} : {w: -1, h: -1};
  });
  expect(cvDims.w, 'canvas width changed — recalibrate the 38000-px blank threshold in chartCanvasNonEmpty').toBeGreaterThan(450);
  expect(cvDims.h, 'canvas height changed — recalibrate the 38000-px blank threshold in chartCanvasNonEmpty').toBeGreaterThan(200);

  await softStep('S1: enable Multi Axis with 2 Y columns', async () => {
    const before = realErrors().length;
    await setProps(page, {multiAxis: true});
    expect((await getProps(page, 'multiAxis')).multiAxis).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  // github-2904: a split under multi-axis must keep the chart rendering.
  await softStep('S1: first split column, chart stays non-empty', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toEqual(['Stereo Category']);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  // github-2904: the second split column is where the chart went blank.
  await softStep('S1: second split column, chart does not go blank', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category', 'Series']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(2);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    const followup = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 100);
    expect(followup).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: 3 Y columns survive an edit, not reset to 1', async () => {
    await setProps(page, {splitColumnNames: [], yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(3);
    // GROK-18484: EDITING the Y list through the property-panel multi-select
    // dialog must not silently reset it. Drive the real UI: the Y row of the
    // context panel's property grid, its [...] editor, the Select-columns
    // dialog's Search input, Escape. The dialog's column LIST is canvas-rendered
    // (no DOM rows/checkboxes), so a checkbox-toggle+OK mutation is not reachable
    // in automation — open/search/close is the deepest scriptable editor
    // interaction; the wrong-dialog and reset failure modes are held by the
    // '3 checked' identity assert and the exact-set assert below.
    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
    // The context panel usually auto-opens for a selected viewer; click the
    // gear only when the Y property row is not present (a blind gear click
    // toggles an open panel CLOSED).
    const clickYDots = () => page.evaluate(() => {
      // Viewer props render as table.property-grid rows: TR.property-grid-item
      // with the prop name in TD.property-grid-item-name. The Y section header
      // also carries the text 'Y' but has .property-grid-category.
      const rows = Array.from(document.querySelectorAll('table.property-grid tr.property-grid-item')) as HTMLElement[];
      const yRow = rows.find((tr) => !tr.classList.contains('property-grid-category') &&
        (tr.querySelector('td.property-grid-item-name')?.textContent ?? '').trim() === 'Y');
      if (!yRow) return false;
      const dots = (Array.from(yRow.querySelectorAll('*')) as HTMLElement[])
        .find((el) => el.childElementCount === 0 && /^(\.\.\.|…)$/.test((el.textContent ?? '').trim()));
      if (!dots) return false;
      dots.click();
      return true;
    });
    let dotsClicked = await clickYDots();
    if (!dotsClicked) {
      await page.evaluate(() => {
        const el = document.querySelector('[name="viewer-Line-chart"]') as HTMLElement;
        const panelBase = el.closest('.panel-base') as HTMLElement;
        (panelBase?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement)?.click();
      });
      await page.waitForTimeout(900);
      dotsClicked = await clickYDots();
    }
    expect(dotsClicked).toBe(true);
    const search = page.locator('.d4-dialog .d4-column-grid input[placeholder="Search"]');
    await search.waitFor({timeout: 5000});
    // Identity: this dialog belongs to the Y row — exactly the 3 Y columns are
    // pre-checked (Split's dialog would read '0 checked' here).
    const checkedLabel = await page.evaluate(() =>
      (Array.from(document.querySelectorAll('.d4-dialog *')) as HTMLElement[])
        .filter((el) => el.childElementCount === 0)
        .map((el) => (el.textContent ?? '').trim())
        .find((t) => /^\d+ checked$/.test(t)) ?? null);
    expect(checkedLabel).toBe('3 checked');
    await search.fill('Chemical');
    await page.waitForTimeout(600);
    // Close deterministically via the CANCEL button (an Escape can be consumed
    // by the focused search input), then prove THIS dialog is gone. The dialog
    // root is position:fixed (offsetParent is null even when shown) — check
    // closure by rect. The exact-set assert below makes a wrong-row click or a
    // silent reset fail loudly, not pass on a length check.
    await page.locator('.d4-dialog button', {hasText: 'CANCEL'}).first().click();
    await page.waitForTimeout(800);
    // Same marker as the open-detect above (.d4-column-grid), not the
    // title-derived dlg-select-columns class — a dialog-title change would
    // silently turn a slug-based check vacuous.
    const openDialogs = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-dialog'))
        .filter((d) => d.getBoundingClientRect().width > 0 &&
          d.querySelector('.d4-column-grid')).length);
    expect(openDialogs).toBe(0);
    const y = (await getProps(page, 'yColumnNames')).yColumnNames;
    expect(y).toEqual(['Chemical Space X', 'Chemical Space Y', 'TPSA']);
  });

  // GROK-20033: the column-selector search input must stay inside the selector.
  await softStep('S2: search input is inside the selector bounds', async () => {
    const before = realErrors().length;
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      const combos = Array.from(lc.root.querySelectorAll('[name^="div-column-combobox"]')) as HTMLElement[];
      const yCombo = combos.find((c) => c.getBoundingClientRect().x > 800);
      (yCombo ?? combos[0])?.click();
    });
    await page.waitForTimeout(800);
    const result = await page.evaluate(() => {
      const input = Array.from(document.querySelectorAll('input.ui-input-editor'))
        .find((i) => (i as HTMLElement).offsetParent !== null &&
          (i as HTMLInputElement).placeholder?.includes('Search')) as HTMLInputElement | undefined;
      if (!input) return {present: false, contained: false};
      const container = input.closest('.ui-input-root');
      if (!container) return {present: true, contained: false};
      const ir = input.getBoundingClientRect();
      const cr = container.getBoundingClientRect();
      const contained = ir.top >= cr.top - 1 && ir.bottom <= cr.bottom + 1;
      return {present: true, contained};
    });
    if (result.present)
      expect(result.contained).toBe(true); // GROK-20033 invariant
    else
      expect(realErrors().length).toBe(before); // no-error floor fallback (search input absent)
  });

  await softStep('S3: disable Multi Axis', async () => {
    const before = realErrors().length;
    await setProps(page, {multiAxis: false});
    expect((await getProps(page, 'multiAxis')).multiAxis).toBe(false);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S3: clear all split columns', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: []});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(0);
    expect(realErrors().length).toBe(before);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart multi-axis-and-split failures:\n${stepErrors.join('\n')}`);
});
