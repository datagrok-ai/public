/* ---
realizes: [matrixplot.cp.axes-layout-and-style, matrixplot-auto-layout-overrides-axis-toggles, viewers.matrix-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Stack trace [A-Za-z]+/.test(text);

/** Shared axis tick-strip canvas width (0 when the strip is collapsed). */
const stripWidth = (page: Page, side: 'top' | 'left') => page.evaluate((s: string) => {
  const c = document.querySelector(`[name="viewer-Matrix-plot"] .d4-layout-${s} canvas`) as HTMLCanvasElement | null;
  return c ? c.width : -1;
}, side);

/** Computed font of a matrix column-label div (re-queried — the divs are
 * re-created on a font change, so a captured reference goes stale). */
const labelFont = (page: Page) => page.evaluate(() => {
  const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
  const d = [...root.querySelectorAll('div')]
    .find((e) => e.children.length === 0 && e.textContent!.trim() === 'AGE');
  return d ? getComputedStyle(d).font : null;
});

const readProp = (page: Page, prop: string) => page.evaluate((p: string) => {
  const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
  return mp?.props[p];
}, prop);

/** Toggle a property-grid checkbox via a synthetic `.click()` — the input is
 * custom-rendered (visually hidden) but the click routes through its onChange
 * handler, keeping the drive on the checkbox path rather than the JS-API prop. */
async function clickCheckbox(page: Page, name: string) {
  await page.evaluate((n: string) => {
    const cb = document.querySelector(`[name="${n}"]`) as HTMLElement | null;
    if (!cb) throw new Error(`checkbox ${n} not found`);
    cb.scrollIntoView({block: 'center'});
    cb.click();
  }, name);
}

/** Drive the settings Font row's size textbox and commit (GROK-18736 path). */
async function setFontSize(page: Page, size: number) {
  await page.evaluate((n: number) => {
    const inp = document.querySelector('[name="prop-view-font"] input.d4-font-size-input') as HTMLInputElement;
    inp.focus();
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(inp, String(n));
    inp.dispatchEvent(new Event('input', {bubbles: true}));
    inp.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter', keyCode: 13} as any));
    inp.dispatchEvent(new Event('change', {bubbles: true}));
  }, size);
  await page.waitForTimeout(700);
}

test('Matrix Plot — Axes Visibility, Auto Layout, Font', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');
  await v.openViewerGear(page, 'Matrix plot');

  // Expand the Axes category if collapsed so the real Show X/Y Axes
  // checkboxes are interactable.
  await page.evaluate(() => {
    const cat = document.querySelector('[name="prop-category-axes"]') as HTMLElement | null;
    if (cat && cat.querySelector('.property-grid-icon-plus')) cat.click();
  });
  await page.waitForTimeout(400);

  await softStep('Scenario 1 — Show X Axes checkbox with Auto Layout on (GROK-19106)', async () => {
    // GROK-19106 reproduces only under the default Auto Layout ON, driven
    // through the real checkbox (the JS-API prop path never reproduced it).
    expect(await readProp(page, 'autoLayout')).toBe(true);
    const wShown = await stripWidth(page, 'top');
    expect(wShown).toBeGreaterThan(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'top')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'top')).toBeGreaterThan(0);
  });

  await softStep('Scenario 2 — explicit axes flags with Auto Layout off', async () => {
    await clickCheckbox(page, 'prop-view-auto-layout');
    await page.waitForTimeout(400);
    expect(await readProp(page, 'autoLayout')).toBe(false);

    // X strip
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'top')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'top')).toBeGreaterThan(0);

    // Y strip
    await clickCheckbox(page, 'prop-view-show-y-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'left')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-y-axes');
    await page.waitForTimeout(700);
    expect(await stripWidth(page, 'left')).toBeGreaterThan(0);

    // Round-trip Auto Layout back to the default.
    await clickCheckbox(page, 'prop-view-auto-layout');
    await page.waitForTimeout(400);
    expect(await readProp(page, 'autoLayout')).toBe(true);
  });

  await softStep('Scenario 3 — Font change reaches the labels (GROK-18736)', async () => {
    expect(await labelFont(page)).toMatch(/(^|\s)10px/);
    await setFontSize(page, 24);
    expect(await labelFont(page)).toMatch(/(^|\s)24px/);
    await setFontSize(page, 10);
    expect(await labelFont(page)).toMatch(/(^|\s)10px/);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
