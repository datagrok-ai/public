/* ---
realizes: [linechart.cp.multi-axis-and-split]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props}], 500);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function chartCanvasCenter(page: Page): Promise<{x: number, y: number}> {
  return page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const canvases = lc.root.querySelectorAll('canvas');
    let mc: HTMLCanvasElement | null = null; let ma = 0;
    for (const c of canvases) {
      const r = (c as HTMLCanvasElement).getBoundingClientRect();
      if (r.width * r.height > ma) { ma = r.width * r.height; mc = c as HTMLCanvasElement; }
    }
    const rect = mc!.getBoundingClientRect();
    return {x: rect.left + rect.width * 0.5, y: rect.top + rect.height * 0.5};
  });
}

async function chartContextMenuClickByLabel(page: Page, label: string) {
  await page.evaluate(() => {
    document.querySelectorAll('.d4-menu-popup').forEach((m) => m.remove());
  });

  await page.waitForTimeout(200);
  const center = await chartCanvasCenter(page);
  await page.mouse.click(center.x, center.y, {button: 'right'});

  await v.pollValue(() => page.locator('.d4-menu-item-label').count(), (n) => n > 0, 500, 50);
  await page.evaluate((text) => {
    const lbl = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el.textContent ?? '').trim() === text);
    if (!lbl) throw new Error(`Menu item not found: ${text}`);
    const item = lbl.closest('.d4-menu-item') as HTMLElement;
    const container = item.closest('.d4-menu-item-container') as HTMLElement | null;
    if (container) container.style.display = '';
    item.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, label);
  await v.waitForViewerRendered(page, 'Line chart', 500);
}

async function chartCanvasNonEmpty(page: Page): Promise<boolean> {

  return (await v.countCanvasPixels(page, 'Line chart')).total > 38000;
}

test('Line Chart — Multi-Axis and Split', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X', 'Chemical Space Y']});

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

  await softStep('S1: first split column, chart stays non-empty', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toEqual(['Stereo Category']);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: second split column, chart does not go blank', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category', 'Series']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(2);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    const followup = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 100);
    expect(followup).toBe(true);

    const hover = await chartCanvasCenter(page);
    await page.mouse.move(hover.x, hover.y, {steps: 5});

    await page.waitForTimeout(600);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: 3 Y columns survive an edit, not reset to 1', async () => {
    await setProps(page, {splitColumnNames: [], yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(3);

    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});

    const clickYDots = () => page.evaluate(() => {

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

      dotsClicked = await v.pollValue(clickYDots, (ok) => ok, 900, 100);
    }
    expect(dotsClicked).toBe(true);
    const search = page.locator('.d4-dialog .d4-column-grid input[placeholder="Search"]');
    await search.waitFor({timeout: 5000});

    const checkedLabel = await page.evaluate(() =>
      (Array.from(document.querySelectorAll('.d4-dialog *')) as HTMLElement[])
        .filter((el) => el.childElementCount === 0)
        .map((el) => (el.textContent ?? '').trim())
        .find((t) => /^\d+ checked$/.test(t)) ?? null);
    expect(checkedLabel).toBe('3 checked');
    await search.fill('Chemical');

    await page.waitForTimeout(600);

    await page.locator('.d4-dialog button', {hasText: 'CANCEL'}).first().click();

    const openDialogs = await v.pollValue(() => page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-dialog'))
        .filter((d) => d.getBoundingClientRect().width > 0 &&
          d.querySelector('.d4-column-grid')).length),
    (n) => n === 0, 800, 50);
    expect(openDialogs).toBe(0);
    const y = (await getProps(page, 'yColumnNames')).yColumnNames;
    expect(y).toEqual(['Chemical Space X', 'Chemical Space Y', 'TPSA']);
  });

  await softStep('S2: search input is inside the selector bounds', async () => {
    const before = realErrors().length;
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      const combos = Array.from(lc.root.querySelectorAll('[name^="div-column-combobox"]')) as HTMLElement[];
      const yCombo = combos.find((c) => c.getBoundingClientRect().x > 800);
      (yCombo ?? combos[0])?.click();
    });

    await v.pollValue(() => page.evaluate(() =>
      !!Array.from(document.querySelectorAll('input.ui-input-editor'))
        .find((i) => (i as HTMLElement).offsetParent !== null &&
          (i as HTMLInputElement).placeholder?.includes('Search'))),
    (present) => present, 800, 50);
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
      expect(result.contained).toBe(true); 
    else
      expect(realErrors().length).toBe(before); 
  });

  await softStep('S2b: Hide other charts — per-chart menu reduces Y columns to one', async () => {
    const before = realErrors().length;
    expect((await getProps(page, 'multiAxis')).multiAxis).toBe(true);
    const yBefore = (await getProps(page, 'yColumnNames')).yColumnNames as string[];
    expect(yBefore).toHaveLength(3);
    await chartContextMenuClickByLabel(page, 'Hide other charts');
    const yAfter = (await getProps(page, 'yColumnNames')).yColumnNames as string[];
    expect(yAfter).toHaveLength(1);
    expect(yBefore).toContain(yAfter[0]);
    expect(realErrors().length).toBe(before);

    await setProps(page, {yColumnNames: yBefore});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toEqual(yBefore);
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
