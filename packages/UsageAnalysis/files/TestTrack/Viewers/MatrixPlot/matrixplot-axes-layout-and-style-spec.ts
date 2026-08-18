/* ---
realizes: [matrixplot.cp.axes-layout-and-style, matrixplot.int.auto-layout-overrides-axis-toggles, viewers.matrix-plot]
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

const stripWidth = (page: Page, side: 'top' | 'left') => page.evaluate((s: string) => {
  const c = document.querySelector(`[name="viewer-Matrix-plot"] .d4-layout-${s} canvas`) as HTMLCanvasElement | null;
  return c ? c.width : -1;
}, side);

const waitForStrip = (page: Page, side: 'top' | 'left', ok: (w: number) => boolean, capMs: number) =>
  v.pollValue(() => stripWidth(page, side), ok, capMs, 50);

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

async function clickCheckbox(page: Page, name: string) {
  await page.evaluate((n: string) => {
    const cb = document.querySelector(`[name="${n}"]`) as HTMLElement | null;
    if (!cb) throw new Error(`checkbox ${n} not found`);
    cb.scrollIntoView({block: 'center'});
    cb.click();
  }, name);
}

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

  await v.pollValue(() => labelFont(page),
    (f) => new RegExp(`(^|\\s)${size}px`).test(f ?? ''), 700, 50);
}

test('Matrix Plot — Axes Visibility, Auto Layout, Font', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');
  await v.openViewerGear(page, 'Matrix plot');

  await page.evaluate(() => {
    const cat = document.querySelector('[name="prop-category-axes"]') as HTMLElement | null;
    if (cat && cat.querySelector('.property-grid-icon-plus')) cat.click();
  });

  await v.pollValue(() => page.evaluate(() => {
    const cat = document.querySelector('[name="prop-category-axes"]');
    return !!cat && !cat.querySelector('.property-grid-icon-plus');
  }), (expanded) => expanded, 400, 25);

  await softStep('Scenario 1 — Show X Axes checkbox with Auto Layout on (GROK-19106)', async () => {

    expect(await readProp(page, 'autoLayout')).toBe(true);
    const wShown = await stripWidth(page, 'top');
    expect(wShown).toBeGreaterThan(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await waitForStrip(page, 'top', (w) => w === 0, 700);
    expect(await stripWidth(page, 'top')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await waitForStrip(page, 'top', (w) => w > 0, 700);
    expect(await stripWidth(page, 'top')).toBeGreaterThan(0);
  });

  await softStep('Scenario 2 — explicit axes flags with Auto Layout off', async () => {
    await clickCheckbox(page, 'prop-view-auto-layout');
    await v.pollValue(() => readProp(page, 'autoLayout'), (on) => on === false, 400, 25);
    expect(await readProp(page, 'autoLayout')).toBe(false);

    await clickCheckbox(page, 'prop-view-show-x-axes');
    await waitForStrip(page, 'top', (w) => w === 0, 700);
    expect(await stripWidth(page, 'top')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-x-axes');
    await waitForStrip(page, 'top', (w) => w > 0, 700);
    expect(await stripWidth(page, 'top')).toBeGreaterThan(0);

    await clickCheckbox(page, 'prop-view-show-y-axes');
    await waitForStrip(page, 'left', (w) => w === 0, 700);
    expect(await stripWidth(page, 'left')).toBe(0);
    await clickCheckbox(page, 'prop-view-show-y-axes');
    await waitForStrip(page, 'left', (w) => w > 0, 700);
    expect(await stripWidth(page, 'left')).toBeGreaterThan(0);

    await clickCheckbox(page, 'prop-view-auto-layout');
    await v.pollValue(() => readProp(page, 'autoLayout'), (on) => on === true, 400, 25);
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
