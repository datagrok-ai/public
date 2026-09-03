/* ---
realizes: [boxplot.cp.value-category-axes-persist]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {deleteProjectWithCleanup} from '../../helpers/projects';
import {BOX, bpProp, viewportRect, readLadder, dragTopHandle} from './boxplot-helpers';

declare const grok: any;
declare const DG: any;

// saveProjectViaApi with the view state taken from saveLayout({saveWithData: true}): a plain
// getInfo() strips the viewport (viewer_base.dart removeDefaultValues) unless the view was flagged
// in-project, which only the ribbon's Save dialog does.
async function saveProjectWithViewport(page: Page, name: string): Promise<{projectId: string}> {
  return page.evaluate(async (n) => {
    const tv = grok.shell.tv;
    const project = DG.Project.create();
    project.name = n;
    const tableInfo = tv.dataFrame.getTableInfo();
    const viewInfo = tv.getInfo();
    viewInfo.viewState = tv.saveLayout({saveWithData: true}).viewState;
    project.addChild(tableInfo);
    project.addChild(viewInfo);
    await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.views.save(viewInfo);
    await grok.dapi.projects.save(project);
    return {projectId: String(project.id)};
  }, name);
}

// Scenario 5 of the settings ladder: the configured box plot saved as a project, closed and
// reopened. The ladder that proves each setting takes is boxplot-settings-ladder-spec.ts on the
// local lane; the state is set up here directly through the properties.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Box Plot settings ladder — project round-trip', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const bp = grok.shell.tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
  });
  await page.locator(BOX).waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);

  let projectIds: {projectId: string} | null = null;
  try {
    await softStep('Scenario 5 setup: the full ladder with group comparison on, then a range-slider zoom', async () => {
      // The ladder's own order. A category change re-derives the marker colour on a later tick, and
      // an explicit colour set during that tick does not count as explicit (box_plot_core.dart
      // _colorAutoResetting), so each category step waits for its derived colour to land first.
      const markerColor = () => bpProp(page, 'markerColorColumnName');
      await v.setViewerProps(page, 'Box plot', [{set: {category1ColumnName: 'SEX'}, wait: 1500}]);
      expect(await v.pollValue(markerColor, (c) => c === 'SEX', 2000, 50)).toBe('SEX');
      await v.setViewerProps(page, 'Box plot', [
        {set: {markerColorColumnName: 'HEIGHT', invertColorScheme: true, colorMin: 20, colorMax: 80}, wait: 800},
        {set: {category2ColumnName: 'RACE', showMinorCategories: true, showAllCategories: true}, wait: 900},
      ]);
      expect(await v.pollValue(markerColor, (c) => c !== 'HEIGHT', 1000, 50)).toBe('HEIGHT');
      await v.setViewerProps(page, 'Box plot', [
        {set: {valueColumnName: 'WEIGHT'}, wait: 1200},
        {set: {axisType: 'logarithmic', invertYAxis: true, plotStyle: 'violin'}, wait: 1400},
        {set: {markerColorColumnName: 'SEX'}, wait: 1200},
        {set: {showGroupComparison: true, controlComparisons: true, controlGroup: 'F'}, wait: 1200},
        {set: {covariateColumnName: 'HEIGHT'}, wait: 1500},
      ]);
      expect(await readLadder(page)).toMatchObject({
        valueColumnName: 'WEIGHT', category1ColumnName: 'SEX', category2ColumnName: 'RACE',
        markerColorColumnName: 'SEX', invertColorScheme: true, colorMin: 20, colorMax: 80,
        axisType: 'logarithmic', invertYAxis: true, plotStyle: 'violin',
      });
      expect(await bpProp(page, 'showGroupComparison')).toBe(true);
      expect(await bpProp(page, 'covariateColumnName')).toBe('HEIGHT');
      await v.waitForCanvasQuiet(page, 'Box plot');

      const before = await viewportRect(page);
      await dragTopHandle(page, 0.30);
      const zoomed = await v.pollValue(() => viewportRect(page),
        (vp) => vp.height < before.height * 0.9, 1000, 150);
      console.log('Scenario 5 viewport height before/after zoom:', before.height, zoomed.height);

      expect(zoomed.height).toBeLessThan(before.height * 0.9);
      await page.evaluate((vp) => (window as any).__preSaveViewport = vp, zoomed);
      const res = await saveProjectWithViewport(page, 'zz-boxplot-ladder-' + Date.now());
      projectIds = {projectId: res.projectId};
      expect(res.projectId.length).toBeGreaterThan(0);
    });

    await softStep('[anchor: PROJECT-ROUND-TRIP] Scenario 5 Step 4-5: reopen the project restores the ladder, zoom, and group-comparison props', async () => {
      const preSave = await page.evaluate(() => (window as any).__preSaveViewport);
      await v.closeAllAndWait(page);
      await page.evaluate(async (pid) => {
        const p = await grok.dapi.projects.find(pid);
        await p.open();
      }, projectIds!.projectId);
      await page.locator(BOX).waitFor({timeout: 30000});

      const ladder = await v.pollValue(() => readLadder(page),
        (l) => l.valueColumnName === 'WEIGHT' && l.axisType === 'logarithmic' && l.plotStyle === 'violin',
        3000, 150);
      console.log('Scenario 5 restored ladder:', JSON.stringify(ladder));
      expect(ladder.valueColumnName).toBe('WEIGHT');
      expect(ladder.category1ColumnName).toBe('SEX');
      expect(ladder.category2ColumnName).toBe('RACE');
      expect(ladder.markerColorColumnName).toBe('SEX');
      expect(ladder.invertColorScheme).toBe(true);
      expect(ladder.colorMin).toBe(20);
      expect(ladder.colorMax).toBe(80);
      expect(ladder.axisType).toBe('logarithmic');
      expect(ladder.invertYAxis).toBe(true);
      expect(ladder.plotStyle).toBe('violin');
      const restoredVp = await v.pollValue(() => viewportRect(page),
        (vp) => Math.abs(vp.height - preSave.height) < 0.5 && Math.abs(vp.top - preSave.top) < 0.5,
        3000, 150);
      console.log('Scenario 5 viewport preSave/restored:', JSON.stringify(preSave), JSON.stringify(restoredVp));
      expect(restoredVp.height).toBeCloseTo(preSave.height, 0);
      expect(restoredVp.top).toBeCloseTo(preSave.top, 0);
      expect(await bpProp(page, 'showGroupComparison')).toBe(true);
      expect(await bpProp(page, 'controlGroup')).toBe('F');
      expect(await bpProp(page, 'covariateColumnName')).toBe('HEIGHT');
    });
  } finally {
    if (projectIds) await deleteProjectWithCleanup(page, projectIds);
  }

  await v.closeAllAndWait(page);
  v.finishSpec();
});
