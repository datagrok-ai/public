/* ---
realizes: [scatterplot.cp.axes-and-encode, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import * as sp from './scatterplot-shared';

declare const grok: any;

// The server lane of the axes-and-encode scenario: Scenario 2 (its formula line comes from the
// PowerPack dialog) and Scenario 5, the layout and project round-trips. Scenarios 1, 3 and 4
// are scatterplot-axes-and-encode-spec.ts on the local lane; the peak state is set up here
// directly through the API.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const PEAK = {x: 'AGE', y: 'HEIGHT', color: 'RACE', size: 'WEIGHT', markers: 'SEX'};

let inProjectSaveWindow = false;
const isBenignError = (text: string) => sp.isAmbientError(text) || (inProjectSaveWindow && (
  /Unable to find element in cloned iframe/.test(text) || /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text)));

const readConfig = (page: Page) => page.evaluate(() => {
  let s: any = null;
  for (const view of grok.shell.tableViews)
    for (const vw of view.viewers)
      if (vw.type === 'Scatter plot') s = vw;
  return {
    found: !!s,
    x: s?.props.xColumnName, y: s?.props.yColumnName, color: s?.props.colorColumnName,
    size: s?.props.sizeColumnName, markers: s?.props.markersColumnName,
  };
});

const selectorLabels = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  const label = (r: string) => (s.root.querySelector(
    `[name="div-column-combobox-${r}"] .d4-column-selector-column`)?.textContent ?? '').trim();
  return {x: label('x'), y: label('y')};
});

async function renameColumnViaGrid(page: Page, current: string, next: string): Promise<void> {
  const pt = await page.evaluate((name: string) => {
    const grid = grok.shell.tv.grid;
    const r = grid.root.getBoundingClientRect();
    const gc = grid.columns.byName(name);
    return {x: r.x + gc.left + gc.width / 2, y: r.y + (grid.colHeaderH || 20) / 2};
  }, current);
  await page.mouse.click(pt.x, pt.y, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.locator('.d4-menu-popup [name="div-Column-Properties..."]').last().click();
  const nameInput = page.locator('.d4-dialog input[name="input-New-name--"]');
  await nameInput.waitFor({timeout: 8000});
  await nameInput.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(next);
  await v.pollValue(() => nameInput.inputValue(), (val) => val === next, 1500, 50);
  await page.locator('.d4-dialog [name="button-OK"]').last().click();
  await v.pollValue(() => page.evaluate((n: string) =>
    document.querySelectorAll('.d4-dialog').length === 0 &&
    grok.shell.tv.dataFrame.columns.names().includes(n), next), (ok) => ok, 5000, 50);
}

const FORMULA_DIALOG = '[name="dialog-Formula-Lines"]';

async function openFormulaLines(page: Page): Promise<void> {
  await sp.openPlotContextMenu(page);
  await sp.navigatePopup(page, ['div-Tools', 'div-Tools---Formula-Lines...']);
  await page.locator(FORMULA_DIALOG).waitFor({timeout: 10000});
  await page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`).waitFor({state: 'visible', timeout: 5000});
}

const formulaEditorValues = (page: Page) => page.evaluate((dlg: string) =>
  ([...document.querySelectorAll(`${dlg} textarea`)] as HTMLTextAreaElement[])
    .map((a) => a.value).join('|'), FORMULA_DIALOG);

const formulaLineCount = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  return JSON.parse(s.props.formulaLines || '[]').length;
});

test('Scatter Plot — Axes and Encodings Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(300_000);

  const errors = sp.trackErrors(page, isBenignError);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Setup — a scatter plot at the peak configuration', async () => {
    await page.evaluate((peak) => grok.shell.tv.addViewer('Scatter plot', {
      xColumnName: peak.x, yColumnName: peak.y, colorColumnName: peak.color,
      sizeColumnName: peak.size, markersColumnName: peak.markers,
    }), PEAK);
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 15_000});
    await v.installEventWaits(page);
    expect(await readConfig(page)).toEqual({found: true, ...PEAK});
  });

  await softStep('Scenario 2 — One column on both axes, then renamed (GROK-19334)', async () => {
    const errBefore = errCount();
    const probeName = 'AGE_PROBE';

    await sp.pickPanelColumn(page, 'prop-y', 'div-column-combobox-y', 'y', 'AGE');
    const both = await readConfig(page);
    expect(both.x).toBe('AGE');
    expect(both.y).toBe('AGE');

    await openFormulaLines(page);
    const noLines = await formulaEditorValues(page);
    await page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`).click();
    await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
    await sp.navigatePopup(page, ['div-Line']);
    await v.pollValue(() => formulaEditorValues(page), (cur) => cur !== noLines, 4000, 50);
    await page.locator(`${FORMULA_DIALOG} [name="button-OK"]`).click();
    expect(await v.pollValue(() => formulaLineCount(page), (n) => n === 1, 4000, 50)).toBe(1);

    try {
      await renameColumnViaGrid(page, 'AGE', probeName);
      const renamed = await v.pollValue(async () => ({...await readConfig(page), labels: await selectorLabels(page)}),
        (c) => c.x === probeName && c.labels.x === probeName && c.labels.y === probeName, 3000, 50);
      expect(renamed.x).toBe(probeName);
      expect(renamed.y).toBe(probeName);
      expect(renamed.labels.x).toBe(probeName);
      expect(renamed.labels.y).toBe(probeName);
      expect(errors.all().slice(errBefore)).toEqual([]);
    } finally {
      const names = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
      if (names.includes(probeName)) await renameColumnViaGrid(page, probeName, 'AGE');
      await openFormulaLines(page);
      const withLine = await formulaEditorValues(page);
      await page.locator(`${FORMULA_DIALOG} [name="button-Delete"]`).click();
      await v.pollValue(() => formulaEditorValues(page), (cur) => cur !== withLine, 2500, 50);
      await page.locator(`${FORMULA_DIALOG} [name="button-OK"]`).click();
      await v.pollValue(() => formulaLineCount(page), (n) => n === 0, 4000, 50);
      await sp.pickPanelColumn(page, 'prop-y', 'div-column-combobox-y', 'y', 'HEIGHT');
    }
    const reverted = await readConfig(page);
    expect(reverted.x).toBe('AGE');
    expect(reverted.y).toBe('HEIGHT');
    expect(await formulaLineCount(page)).toBe(0);
  });

  await softStep('Scenario 5a — Layout and project persistence at peak configuration — layout round-trip (GROK-18945)', async () => {
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-scatterplot-axes-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await page.evaluate(() => grok.shell.tv.addViewer('Histogram'));
      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 10_000}).toBe(1);

      const [clearedColor] = await v.setViewerProps(page, sp.SP_TYPE,
        [{set: {colorColumnName: ''}, wait: 800, read: 'colorColumnName'}]);
      expect(clearedColor).toBe('');
      await page.evaluate(async (id) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);
      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const restored = tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
        return {
          hasScatter: !!restored,
          hasHistogram: tv.viewers.some((vw: any) => vw.type === 'Histogram'),
          x: restored?.props.xColumnName, y: restored?.props.yColumnName,
          color: restored?.props.colorColumnName, size: restored?.props.sizeColumnName,
          markers: restored?.props.markersColumnName,
        };
      }), (r) => r.hasScatter && !r.hasHistogram && r.color === PEAK.color, 6000, 100);
      expect(result).toEqual({hasScatter: true, hasHistogram: false, ...PEAK});
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Scenario 5b — Layout and project persistence at peak configuration — project save / Close All / reopen', async () => {
    let projectId: string | null = null;
    inProjectSaveWindow = true;
    try {
      projectId = (await saveProjectViaApi(page, 'zz-scatterplot-axes-encode-probe-' + Date.now())).projectId;
      expect(projectId).toBeTruthy();

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => (await grok.dapi.projects.find(id)).open(), projectId);
      await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 30_000});

      const result = await v.pollValue(() => readConfig(page), (r) => r.found && r.color === PEAK.color, 10_000, 200);
      expect(result).toEqual({found: true, ...PEAK});
    } finally {
      inProjectSaveWindow = false;
      if (projectId) await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
