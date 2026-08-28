/* ---
realizes: [pcplot.cp.layout-project-persistence]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// The server lane of the mixed setup-color-filter scenario: the layout round-trip and the
// project save / Close All / reopen are the only steps whose subject IS server state.
// Everything else about the PC plot is client behaviour and stays in pcplot-setup-color-filter-spec.ts.

test('PC Plot — layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pc-plot', 'PC-Plot', 15000);

  await v.installEventWaits(page);

  await softStep('Layout round-trip — saved layout restores the configured viewer set and props', async () => {

    await v.setViewerProps(page, 'PC Plot', [{
      set: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT'], colorColumnName: 'RACE', title: 'PC Persistence Probe'},
      wait: 800,
    }]);
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      await page.evaluate(() => { grok.shell.tv.addViewer('Scatter plot'); });
      await v.pollValue(
        () => page.evaluate(() => grok.shell.tv.viewers.some((vw: any) => vw.type === 'Scatter plot')),
        (present) => present, 500, 100);
      await page.evaluate(async (id) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);
      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot');
        return {
          hasScatter: tv.viewers.some((vw: any) => vw.type === 'Scatter plot'),
          hasPc: tv.viewers.some((vw: any) => vw.type === 'PC Plot'),
          cols: pc?.props.columnNames?.slice(),
          color: pc?.props.colorColumnName,
          title: pc?.props.title,
        };
      }), (r) => r.hasPc && !r.hasScatter, 3000, 150);

      expect(result.hasScatter).toBe(false);
      expect(result.hasPc).toBe(true);

      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {

      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved)
            await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Project save / Close All / reopen — project restores the configured viewer', async () => {
    const projName = 'zz-pcplot-persistence-probe-' + Date.now();
    let projectId: string | undefined;
    try {
      // every assertion below is ordinary viewer @Prop state, so the API save path applies
      // (see helpers-registry.yaml for the saveProjectViaApi / saveProjectViaUI boundary)
      const saved = await saveProjectViaApi(page, projName);
      projectId = saved.projectId;
      expect(projectId).toBeTruthy();

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => {
        const full = await grok.dapi.projects.find(id);
        await full.open();
      }, projectId);

      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const pc = tv ? Array.from(tv.viewers).find((x: any) => x.type === 'PC Plot') as any : null;
        return {
          pcRestored: (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'PC Plot'),
          cols: pc?.props?.columnNames?.slice(),
          color: pc?.props?.colorColumnName,
          title: pc?.props?.title,
        };
      }), (r) => r.pcRestored, 4500, 150);

      expect(result.pcRestored).toBe(true);

      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {
      await deleteProjectWithCleanup(page, {projectId});
    }
  });

  v.finishSpec();
});
