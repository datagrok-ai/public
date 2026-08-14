import {test, expect} from '@playwright/test';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

// The legend model is two sections (main = colour categories, extra = markers), and the view
// reconciles by key rather than rebuilding. Changing one section must leave the other's DOM
// elements untouched — that is what makes a rebuild cheap and non-jumpy.

test('Legend sections — colour and markers are independent sections', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv'});
  await v.addLegendViewers(page, {column: 'RACE', viewers: ['Scatter plot']});
  await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: 'SEX'}, wait: 1500}]);
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('Two sections, each with its own items', async () => {
    const state = await v.getLegendState(page, 'Scatter plot');
    const ids = state!.sections.map((s) => s.id).sort();
    expect(ids).toEqual(['extra', 'main']);
    expect(state!.sections.find((s) => s.id === 'main')!.items).toBeGreaterThan(0);
    expect(state!.sections.find((s) => s.id === 'extra')!.items).toBeGreaterThan(0);
  });

  await softStep('Changing markers leaves the main section keys untouched', async () => {
    const mainBefore = await v.getLegendItemKeys(page, 'Scatter plot', 'main');
    const extraBefore = await v.getLegendItemKeys(page, 'Scatter plot', 'extra');

    await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: 'DIS_POP'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');

    const mainAfter = await v.getLegendItemKeys(page, 'Scatter plot', 'main');
    const extraAfter = await v.getLegendItemKeys(page, 'Scatter plot', 'extra');

    expect(mainAfter).toEqual(mainBefore);
    expect(extraAfter).not.toEqual(extraBefore);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });
  v.finishSpec();
});

test('Legend sections — a large category count never enters corner mode', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv'});
  // USUBJID is a per-row id, so the category count is far above the 100-category corner gate.
  await v.addLegendViewers(page, {column: 'USUBJID', viewers: ['Scatter plot']});
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('Mode is docked, not corner', async () => {
    const state = await v.getLegendState(page, 'Scatter plot');
    expect(state!.mode).not.toBe('corner');
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });
  v.finishSpec();
});
