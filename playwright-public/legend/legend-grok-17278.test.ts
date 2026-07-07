/* ---
sub_features_covered: [legend.column, legend.item.color-picker]
--- */
// GROK-17278: legend color customizations serialize into both layout and project state.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('GROK-17278: line chart legend color persists across layout + project round-trips', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Steps 2-3: add Line chart, Split = Stereo Category', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await new Promise((r) => setTimeout(r, 1000));
      const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.props.splitColumnName = 'Stereo Category';
      try { lc.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(items).toBeGreaterThan(0);
  });

  await softStep('Step 4: change R_ONE colour via legend picker (Line chart)', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Line chart',
      category: 'R_ONE',
      rgb: [31, 119, 180],
      hex: '#1f77b4',
      column: 'Stereo Category',
    });
  });

  let layoutId: string | null = null;
  await softStep('Step 5 + Step 9 invariant: save + re-apply layout, R_ONE remains blue', async () => {
    const res = await page.evaluate(async () => {
      const withTimeout = <T>(p: Promise<T>, ms: number, label: string): Promise<T> => {
        let t: any;
        return Promise.race([
          p,
          new Promise<T>((_, rej) => { t = setTimeout(() => rej(new Error(`Timeout ${ms}ms: ${label}`)), ms); }),
        ]).finally(() => clearTimeout(t));
      };
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'GROK17278_' + Date.now();
      try {
        const saved = await withTimeout((window as any).grok.dapi.layouts.save(layout), 30000, 'layouts.save');
        await new Promise((r) => setTimeout(r, 1000));
        const found = await withTimeout((window as any).grok.dapi.layouts.find(saved.id), 15000, 'layouts.find');
        tv.loadLayout(found);
        await new Promise((r) => setTimeout(r, 3500));
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
        return {layoutId: saved.id, ok: true, rOneAfterReload: String(t['R_ONE'] ?? '').toLowerCase()};
      } catch (e: any) {
        return {layoutId: null, ok: false, error: String(e?.message ?? e).slice(0, 200)};
      }
    });
    expect(res.ok, res.ok ? '' : `layout save+reapply failed: ${res.error}`).toBe(true);
    layoutId = res.layoutId;
    expect(res.rOneAfterReload, 'R_ONE retains blue after layout round-trip (GROK-17278)').toBe('#1f77b4');
  });

  let projectId: string | null = null;
  await softStep('Steps 6-8 + Step 8 invariant: project save+closeAll+reopen, R_ONE remains blue', async () => {
    const res = await page.evaluate(async () => {
      let pid: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'GROK17278Proj_' + Date.now();
        const __df = (window as any).grok.shell.tv.dataFrame;
        const __ti = __df.getTableInfo();
        proj.addChild(__ti);
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
      const col = tv.dataFrame.col('Stereo Category');
      const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return {
        phase: 'verified', ok: true, projectId: pid,
        rOneAfter: String(t['R_ONE'] ?? '').toLowerCase(),
      };
    });
    expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
    projectId = res.projectId ?? null;
    expect(res.rOneAfter, 'R_ONE retains blue after project save+close+reopen (GROK-17278)').toBe('#1f77b4');
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async ([lid, pid]: [string | null, string | null]) => {
      if (lid) try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(lid)); } catch (_) {}
      if (pid) try { await (window as any).grok.dapi.projects.delete(await (window as any).grok.dapi.projects.find(pid)); } catch (_) {}
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    }, [layoutId, projectId]);
  });

  v.finishSpec();
});
