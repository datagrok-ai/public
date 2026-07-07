/* ---
sub_features_covered: [legend.column, legend.refresh.on-data-change]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

test('Line chart legend', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Setup + Sc1 steps 1-2: add Line chart, Split=Series → categorical legend', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await new Promise((r) => setTimeout(r, 1000));
      const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.props.splitColumnName = 'Series';
      try { lc.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(items).toBeGreaterThan(0);
  });

  await softStep('Sc1 verification: legend hover surface available (real DOM)', async () => {
    const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]').first();
    if (await legend.count() === 0) {
      const alt = page.locator('[name="viewer-Line chart"] [name="legend"]').first();
      await alt.waitFor({timeout: 10000});
      await alt.hover();
    } else {
      await legend.waitFor({timeout: 10000});
      await legend.hover();
    }
  });

  await softStep('Sc2 steps 1-2: enable Multi Axis (multiAxis=true)', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      try { lc.props.multiAxis = true; } catch (e) { return {multiAxis: false, err: String(e)}; }
      await new Promise((r) => setTimeout(r, 1500));
      return {multiAxis: lc.props.multiAxis};
    });
    expect(res.multiAxis).toBe(true);
  });

  let layoutId1: string | null = null;
  await softStep('Sc3 steps 1-3: save+reapply layout (multiAxis+split persist)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LineChart_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      const lc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      return {layoutId: saved.id, multiAxis: lc?.props.multiAxis, split: lc?.props.splitColumnName};
    });
    layoutId1 = res.layoutId;
    expect(res.multiAxis).toBe(true);
    expect(res.split).toBe('Series');
  });

  await softStep('Sc4 steps 1-2: yColumnNames = [Average Mass, TPSA]', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.props.yColumnNames = ['Average Mass', 'TPSA'];
      await new Promise((r) => setTimeout(r, 2000));
      const items = lc.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {yCols: lc.props.yColumnNames, totalItems: items.length};
    });
    expect(res.yCols).toEqual(['Average Mass', 'TPSA']);
    expect(res.totalItems).toBeGreaterThan(0);
  });

  await softStep('Sc4 step 3: replace Y column → NIBR logP', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.props.yColumnNames = ['Average Mass', 'NIBR logP'];
      await new Promise((r) => setTimeout(r, 1800));
      return {yCols: lc.props.yColumnNames, items: lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length};
    });
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  let layoutId2: string | null = null;
  await softStep('Sc4 steps 4-5: save+reapply layout (new Y persists)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LineChart2_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      const lc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      return {layoutId: saved.id, yCols: lc?.props.yColumnNames};
    });
    layoutId2 = res.layoutId;
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  let projectId: string | null = null;
  await softStep('Sc3 steps 4-5 / Sc4 steps 6-7: project save+close+reopen (FK graceful-degrade)', async () => {
    const res = await page.evaluate(async () => {
      let pid: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'LineChartProj_' + Date.now();
        const __df = (window as any).grok.shell.tv.dataFrame;
        const __ti = __df.getTableInfo();
        proj.addChild(__ti);
        proj.addChild((window as any).grok.shell.tv.saveLayout());
        await (window as any).grok.dapi.tables.uploadDataFrame(__df);
        await (window as any).grok.dapi.tables.save(__ti);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        pid = saved.id;
      } catch (e: any) {
        return {phase: 'save', ok: false, error: String(e).slice(0, 200)};
      }
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1200));
      try {
        const reopened = await (window as any).grok.dapi.projects.find(pid);
        await reopened.open();
      } catch (e: any) {
        return {phase: 'reopen', ok: false, error: String(e).slice(0, 200), projectId: pid};
      }
      await new Promise((r) => setTimeout(r, 3500));
      const tv = (window as any).grok.shell.tv;
      if (!tv) return {phase: 'reopen', ok: false, error: 'no tv after reopen', projectId: pid};
      const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
      return {
        phase: 'verified', ok: true, projectId: pid,
        multiAxis: lc?.props.multiAxis, split: lc?.props.splitColumnName,
        yCols: lc?.props.yColumnNames,
      };
    });
    expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
    projectId = res.projectId ?? null;
    expect(res.multiAxis).toBe(true);
    expect(res.split).toBe('Series');
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async ([lid1, lid2, pid]: [string | null, string | null, string | null]) => {
      for (const id of [lid1, lid2]) {
        if (id) try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch (_) {}
      }
      if (pid) try { await (window as any).grok.dapi.projects.delete(await (window as any).grok.dapi.projects.find(pid)); } catch (_) {}
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    }, [layoutId1, layoutId2, projectId]);
  });

  v.finishSpec();
});
