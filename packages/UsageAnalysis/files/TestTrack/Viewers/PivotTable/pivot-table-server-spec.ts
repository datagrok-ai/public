/* ---
realizes: [pivottable.cp.configure-crosstab-values, pivottable.cp.inner-grid-look-viewers]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as proj from '../../helpers/projects';
import {PIVOT, INNER_CANVAS, openPivot, pivotProps, gridLookColumn, pivotTitleDom, applyLinearColorCoding} from './pivot-helpers';

declare const grok: any;

// The server lane of the pivot section: the steps whose subject is a layout or a project surviving
// a round-trip through the server. The configuration under test is set up directly through the
// viewer API; the ladders that prove the configuration UI itself are the local-lane specs.
test.use(specTestOptions);

const MED_CROSSTAB = {
  groupByColumnNames: ['DIS_POP'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['med'], pivotColumnNames: ['SEVERITY'],
};

const pivotPresent = () => Array.from(grok.shell.tv?.viewers ?? []).some((x: any) => x.type === 'Pivot table');

test('Pivot Table — layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await openPivot(page);

  await softStep('Configure crosstab Scenario 5 Step 9: a non-default aggregation survives a layout saved to the gallery and re-applied via the layout API (github-2535)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: MED_CROSSTAB}], 500);
    expect((await pivotProps(page)).aggTypes).toEqual(['med']);

    // awaiting layouts.save IS the completion signal, the same call the Save to Gallery menu makes
    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-pivot-crosstab-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await v.setViewerProps(page, 'Pivot table',
        [{set: {aggregateAggTypes: ['avg'], groupByColumnNames: ['SEX']}}], 400);
      expect((await pivotProps(page)).aggTypes).toEqual(['avg']);

      await page.evaluate(async (id) => grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id)), layoutId);
      const restored = await v.pollValue(() => page.evaluate(() => {
        const pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
        return {
          member: !!pv2,
          groupBy: pv2?.props.groupByColumnNames, pivot: pv2?.props.pivotColumnNames,
          agg: pv2?.props.aggregateColumnNames, aggTypes: pv2?.props.aggregateAggTypes,
        };
      }), (r) => !!r.groupBy?.includes('DIS_POP'), 5000, 100);
      expect(restored.member).toBe(true);
      expect(restored.groupBy).toContain('DIS_POP');
      expect(restored.pivot).toContain('SEVERITY');
      expect(restored.agg).toContain('AGE');
      expect(restored.aggTypes).toContain('med');
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Configure crosstab Scenario 5 Step 13: the configuration survives a project save / close / reopen round-trip (github-2535)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: MED_CROSSTAB}], 500);
    expect((await pivotProps(page)).aggTypes).toEqual(['med']);

    let projectId: string | null = null;
    try {
      projectId = (await proj.saveProjectViaApi(page, `PivotCrosstabProj${Date.now()}`)).projectId;
      expect(projectId).toBeTruthy();

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => (await grok.dapi.projects.find(id)).open(), projectId);
      await page.waitForFunction(pivotPresent, null, {timeout: 30_000});

      const result = await page.evaluate(() => {
        const pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
        return {
          groupBy: pv2.props.groupByColumnNames, pivot: pv2.props.pivotColumnNames,
          agg: pv2.props.aggregateColumnNames, aggTypes: pv2.props.aggregateAggTypes,
        };
      });
      expect(result.groupBy).toContain('DIS_POP');
      expect(result.pivot).toContain('SEVERITY');
      expect(result.aggTypes).toContain('med');
    } finally {
      if (projectId) await proj.deleteProjectWithCleanup(page, {projectId});
    }
  });

  const SPGI_VALUE_X = 584;
  const SPGI_VALUE_COL = 'avg(CAST Idea ID)';

  await softStep('Inner grid look Scenario 6 Step 9: after re-applying the saved layout, the pivot is present, titled "Pivot Overview", and the coloured column keeps its colour', async () => {
    await openPivot(page, 'System:AppData/Chem/tests/spgi-100.csv');
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i));
      const cats = cols.filter((c: any) => c.type === 'string');
      const nums = cols.filter((c: any) => c.type === 'double' || c.type === 'int');
      pv.props.groupByColumnNames = [cats[0].name, cats[1].name];
      pv.props.pivotColumnNames = [];
      pv.props.aggregateColumnNames = [nums[0].name];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.title = 'Pivot Overview';
    });
    await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});

    await applyLinearColorCoding(page, SPGI_VALUE_X);
    const spgiCoded = await v.pollValue(() => gridLookColumn(page, SPGI_VALUE_COL),
      (c) => c?.colorCodingType === 'Linear', 3000, 100);
    expect(spgiCoded?.colorCodingType).toBe('Linear');

    const layoutId: string = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      layout.name = 'zz-pivot-look-' + Date.now();
      return String((await grok.dapi.layouts.save(layout)).id);
    });
    try {
      await v.setViewerProps(page, 'Pivot table', [{set: {title: 'Perturbed'}}], 500);
      await page.evaluate(async (id) => grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id)), layoutId);
      await page.waitForFunction(pivotPresent, null, {timeout: 15_000});
      await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});

      const restoredTitleDom = await v.pollValue(() => pivotTitleDom(page),
        (t) => t.includes('Pivot Overview'), 3000, 100);
      expect(restoredTitleDom).toContain('Pivot Overview');

      const restoredCol = await v.pollValue(() => gridLookColumn(page, SPGI_VALUE_COL),
        (c) => c?.colorCodingType === 'Linear', 3000, 100);
      expect(restoredCol?.colorCodingType).toBe('Linear');
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Inner grid look Scenario 6 Step 10: driving the ribbon Save opens the Save-project dialog with no error (project-persistence entry point)', async () => {
    const consoleErrors: string[] = [];
    const onConsole = (m: any) => { if (m.type() === 'error' && !/cloned iframe/i.test(m.text())) consoleErrors.push(m.text()); };
    page.on('console', onConsole);
    try {
      await page.locator('[name="button-Save"]').first().click();
      const saveDialog = page.locator('.d4-dialog').filter({hasText: 'Save project'});
      await saveDialog.first().waitFor({timeout: 8000});
      expect(await saveDialog.count()).toBeGreaterThan(0);
      expect(consoleErrors).toEqual([]);

      await page.locator('[name="button-CANCEL"]').first().click().catch(() => page.keyboard.press('Escape'));
      await saveDialog.first().waitFor({state: 'detached', timeout: 5000}).catch(() => {});
    } finally {
      page.off('console', onConsole);
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
