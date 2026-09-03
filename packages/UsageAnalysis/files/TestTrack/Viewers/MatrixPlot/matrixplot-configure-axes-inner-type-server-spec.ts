/* ---
realizes: [matrixplot.cp.configure-axes-inner-type, viewers.matrix-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

// The server lane of the configure scenario: Scenario 5, whose subject is the configured
// Matrix plot surviving a layout and a project round-trip. The peak configuration is set up
// through the API here; the dialog and property-panel paths that reach it are proven in
// matrixplot-configure-axes-inner-type-spec.ts on the local lane.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const PEAK = {xColumnNames: ['AGE', 'HEIGHT'], yColumnNames: ['AGE', 'HEIGHT', 'WEIGHT'], cellPlotType: 'Scatter plot'};

const cellCount = (page: Page) => page.evaluate(() =>
  document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer').length);

const viewerTypes = (page: Page) => page.evaluate(() => {
  const types: string[] = [];
  for (const view of grok.shell.tableViews)
    for (const vw of view.viewers) types.push(vw.type);
  return types;
});

const readSets = (page: Page) => page.evaluate(() => {
  let mp: any = null;
  for (const view of grok.shell.tableViews)
    for (const vw of view.viewers)
      if (vw.type === 'Matrix plot') mp = vw;
  return mp ? {x: mp.props.xColumnNames, y: mp.props.yColumnNames, cellPlotType: mp.props.cellPlotType} : null;
});

test('Matrix Plot — layout and project persistence', async ({page}: {page: Page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Setup — Matrix plot at the peak configuration: X AGE, HEIGHT; Y AGE, HEIGHT, WEIGHT; Scatter plot cells', async () => {
    await page.evaluate((peak) => { grok.shell.tv.addViewer('Matrix plot', peak); }, PEAK);
    await page.locator('[name="viewer-Matrix-plot"]').waitFor({timeout: 15_000});
    await expect.poll(() => cellCount(page), {timeout: 10_000}).toBe(6);
    expect(await readSets(page)).toEqual({x: PEAK.xColumnNames, y: PEAK.yColumnNames, cellPlotType: PEAK.cellPlotType});
  });

  await softStep('Scenario 5a — layout round-trip restores the saved viewer set and config', async () => {
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-matrixplot-configure-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await page.evaluate(() => { grok.shell.tv.addViewer('Scatter plot'); });
      await expect.poll(() => viewerTypes(page), {timeout: 5_000}).toContain('Scatter plot');

      await page.evaluate(async (id) => { grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id)); }, layoutId);
      await expect.poll(async () => {
        const types = await viewerTypes(page);
        return types.includes('Matrix plot') && !types.includes('Scatter plot');
      }, {timeout: 10_000}).toBe(true);
      await expect.poll(() => readSets(page), {timeout: 5_000})
        .toEqual({x: PEAK.xColumnNames, y: PEAK.yColumnNames, cellPlotType: PEAK.cellPlotType});
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Scenario 5b — project save / Close All / reopen restores the Matrix plot (GROK-10925)', async () => {
    let projectId: string | null = null;
    try {
      projectId = (await saveProjectViaApi(page, 'zz-matrixplot-persistence-probe-' + Date.now())).projectId;
      expect(projectId).toBeTruthy();

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => { await (await grok.dapi.projects.find(id)).open(); }, projectId);
      await expect.poll(() => viewerTypes(page), {timeout: 20_000}).toContain('Matrix plot');
      await expect.poll(() => readSets(page), {timeout: 5_000})
        .toEqual({x: PEAK.xColumnNames, y: PEAK.yColumnNames, cellPlotType: PEAK.cellPlotType});
    } finally {
      if (projectId) await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
