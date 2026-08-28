/* ---
realizes: [densityplot.cp.setup-binning-color, viewers.density-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// The server lane of the mixed setup-binning-color scenario: Scenarios 6a and 6b persist the
// configured viewer through dapi.layouts and a saved project, which local mode does not serve.
// Scenarios 1-5 are client behaviour and stay in densityplot-setup-binning-color-spec.ts.

test('Density Plot — layout and project persistence', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'density-plot', 'Density-plot');

  await softStep('Scenario 6a — layout round-trip restores the saved viewer set and config', async () => {

    await v.setViewerProps(page, 'Density plot', [{
      set: {
        xColumnName: 'AGE', yColumnName: 'HEIGHT', bins: 200,
        binShape: 'rectangle', invertColorScheme: true,
      },
      wait: 800,
    }]);
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      await page.waitForTimeout(1000); 
      await page.evaluate(() => { grok.shell.tv.addViewer('Scatter plot'); });
      await v.pollValue(
        () => page.evaluate(() => grok.shell.tv.viewers.some((vw: any) => vw.type === 'Scatter plot')),
        (present) => present, 600, 100);
      await page.evaluate(async (id) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);
      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const dp = tv.viewers.find((vw: any) => vw.type === 'Density plot');
        return {
          hasScatter: tv.viewers.some((vw: any) => vw.type === 'Scatter plot'),
          hasDensity: tv.viewers.some((vw: any) => vw.type === 'Density plot'),
          x: dp?.props.xColumnName, y: dp?.props.yColumnName,
          bins: dp?.props.bins, binShape: dp?.props.binShape,
          invert: dp?.props.invertColorScheme,
        };
      }), (r) => r.hasDensity && !r.hasScatter, 3000, 150);

      expect(result.hasDensity).toBe(true);
      expect(result.hasScatter).toBe(false);

      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.bins).toBe(200);
      expect(result.binShape).toBe('rectangle');
      expect(result.invert).toBe(true);
    } finally {
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Scenario 6b — project save / Close All / reopen restores the Density Plot', async () => {
    const projName = 'zz-densityplot-persistence-probe-' + Date.now();
    let projectId: string | null = null;
    try {
      const saved = await saveProjectViaApi(page, projName);
      projectId = saved.projectId;
      expect(projectId).toBeTruthy();

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => {
        const full = await grok.dapi.projects.find(id);
        await full.open();
      }, projectId);

      const result = await v.pollValue(() => page.evaluate(() => {
        const types: string[] = [];
        let dp: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers) {
            types.push(vw.type);
            if (vw.type === 'Density plot') dp = vw;
          }
        return {
          types,
          x: dp?.props.xColumnName, y: dp?.props.yColumnName,
          bins: dp?.props.bins, binShape: dp?.props.binShape,
          invert: dp?.props.invertColorScheme,
        };
      }), (r) => r.types.includes('Density plot'), 20000, 1000);

      expect(result.types).toContain('Density plot');
      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.bins).toBe(200);
      expect(result.binShape).toBe('rectangle');
      expect(result.invert).toBe(true);
    } finally {
      if (projectId)
        await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  v.finishSpec();
});
