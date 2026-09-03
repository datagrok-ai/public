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

// The server lane of the two mixed PC plot scenarios, in one test so the section pays one
// server-lane table open: To Script renders its output through a platform function that local
// mode does not serve, and the layout / project round-trips are server state. Everything else
// about the PC plot is client behaviour and stays in pc-plot-spec.ts and
// pcplot-setup-color-filter-spec.ts.

test('PC Plot — To Script, layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const pcPlotPresent = () => page.evaluate(() => !!grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot'));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pc-plot', 'PC-Plot', 15000);

  await v.installEventWaits(page);

  await softStep('Menu Ribbon and To Script', async () => {
    expect(await pcPlotPresent()).toBe(true);

    await page.evaluate(async () => {
      const w = window as any;
      const canvas = document.querySelector('[name="viewer-PC-Plot"] canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2}));

      await w.__poll(() => document.querySelectorAll('.d4-menu-item-label').length,
        (n: number) => n > 0, 500);
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await w.__poll(() => document.querySelectorAll('.d4-menu-popup').length,
        (n: number) => n === 0, 800);
    });

    const balloon = await page.evaluate(async () => {
      const w = window as any;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();

      const reset = async () => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        await w.__poll(() => document.querySelectorAll('.d4-menu-popup').length,
          (n: number) => n === 0, 500);
      };
      const attempt = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        const opened = await w.__menuLeaf('To Script', 'To JavaScript').catch(() => false);
        if (!opened) { await reset(); return ''; }
        const balloonText = await w.__poll(
          () => (document.querySelector('.d4-balloon') as HTMLElement | null)?.innerText ?? '',
          (t: string) => t.length > 0, 6000, 250);
        if (balloonText.length > 0) return balloonText;
        await reset();
        return '';
      };

      let text = '';
      for (let a = 0; a < 5 && !text; a++) text = await attempt();
      return {present: text.length > 0, text};
    });

    expect(balloon.present).toBe(true);
    expect(balloon.text).toContain('addViewer');

    await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot');
      if (pc) pc.close();
    });
    await v.pollValue(() => page.evaluate(() => !grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')),
      (gone) => gone, 500, 100);
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-pc-plot"]');
      if (icon) (icon as HTMLElement).click();
    });
    const reopened = await v.pollValue(pcPlotPresent, (present) => present, 1000, 100);
    expect(reopened).toBe(true);
  });

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

  expect(pageErrors).toEqual([]);

  v.finishSpec();
});
