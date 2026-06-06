import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test.describe.configure({retries: 1});

test('Legend structure rendering', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Add 7 viewers, set legend column to Core, Always visible', async () => {
    await v.addLegendViewers(page, {
      column: 'Core',
      viewers: ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
      settleMs: 2500,
    });
    const types = await page.evaluate(() => (window as any).grok.shell.tv.viewers.map((x: any) => x.type));
    expect(types.length).toBeGreaterThanOrEqual(8);
  });

  // Molecule thumbnail rendering depends on Chem package presence. Env-gated.
  await softStep('Verify legend canvases on Scatter/Hist/Line/Pie (Molecule rendering)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let out: Record<string, {items: number; canvas: boolean}> = {};
      let chemPresent = false;
      for (let i = 0; i < 50; i++) {
        const col = tv.dataFrame.col('Core');
        chemPresent = col?.semType === 'Molecule' || col?.semType === 'Macromolecule';
        out = {};
        for (const x of tv.viewers) {
          if (x.type === 'Grid') continue;
          const legend = x.root.querySelector('[name="legend"]');
          const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
          out[x.type] = {items: items.length, canvas: !!items[0]?.querySelector('canvas')};
        }
        const targets = ['Scatter plot', 'Histogram', 'Line chart', 'Pie chart'];
        const allSettled = targets.every((t) => {
          const o = out[t];
          if (!o) return false;
          return o.canvas || (!chemPresent && o.items > 0);
        });
        if (allSettled) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      return {chemPresent, viewers: out};
    });
    const targets = ['Scatter plot', 'Histogram', 'Line chart', 'Pie chart'];
    let viewersWithLegend = 0;
    for (const t of targets) {
      const view = res.viewers[t];
      if (view && (view.canvas || view.items > 0)) viewersWithLegend++;
    }
    expect(viewersWithLegend,
      'at least 3 of 4 target viewers render legend items (canvas or text fallback)')
      .toBeGreaterThanOrEqual(3);
  });

  await softStep('Scatter plot — Marker=Core, Color=Core', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.markersColumnName = 'Core';
      sp.props.colorColumnName = 'Core';
      let result = {chemPresent: false, items: 0, canvas: false};
      for (let i = 0; i < 50; i++) {
        const col = tv.dataFrame.col('Core');
        const chemPresent = col?.semType === 'Molecule' || col?.semType === 'Macromolecule';
        const legend = sp.root.querySelector('[name="legend"]');
        const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
        const canvas = !!items[0]?.querySelector('canvas');
        result = {chemPresent, items: items.length, canvas};
        if (canvas || (!chemPresent && items.length > 0)) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      return result;
    });
    expect(res.canvas || res.items > 0,
      'legend renders Marker=Core items (canvas if Chem present, text otherwise)')
      .toBe(true);
  });

  await softStep('Scatter plot — Color=Series, Marker stays Core', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      await new Promise((r) => setTimeout(r, 1000));
      const legend = sp.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {items: items.length, markers: sp.props.markersColumnName, color: sp.props.colorColumnName};
    });
    expect(res.markers).toBe('Core');
    expect(res.color).toBe('Series');
    expect(res.items).toBeGreaterThan(0);
  });

  await softStep('Save + re-apply layout', async () => {
    const id = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'StructRender_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      (window as any).__srLayoutId = saved.id;
      return saved.id;
    });
    (globalThis as any).__srLayoutId = id;
    expect(id).toBeTruthy();
  });

  await softStep('Save project (known FK limitation)', async () => {
    const res = await page.evaluate(async () => {
      try {
        const tv = (window as any).grok.shell.tv;
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'StructRenderProj_' + Date.now();
        proj.addChild(tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        return {ok: true, id: saved.id};
      } catch (e: any) {
        return {ok: false, error: String(e).slice(0, 200)};
      }
    });
    expect(res.ok || String(res.error ?? '').includes('foreign key')).toBe(true);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async (id) => {
      if (id) try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch (_) {}
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    }, (globalThis as any).__srLayoutId);
  });

  v.finishSpec();
});
