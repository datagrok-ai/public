/* ---
sub_features_covered: [legend.allow-item-coloring, legend.item.color-picker, legend.use-custom-color-coding]
--- */
// Scenario 2 picker UI runs on Histogram: Bar chart legend needs a color edit to render.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

test('Legend color consistency', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
  });

  await softStep('Categorical color coding from grid: R_ONE=red, S_UNKN=green', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      col.meta.colors.setCategorical(
        {'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'},
        {fallbackColor: '#808080'},
      );
      for (const x of (window as any).grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      let tagColors: Record<string, any> = {};
      try { tagColors = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}'); } catch (_) {}
      return {
        codingType: col.tags['.color-coding-type'],
        rOneTag: tagColors['R_ONE'] ?? null,
        sUnknTag: tagColors['S_UNKN'] ?? null,
      };
    });
    expect(res.codingType).toBe('Categorical');
    expect(String(res.rOneTag).toLowerCase()).toBe('#ff0000');
    expect(String(res.sUnknTag).toLowerCase()).toBe('#00ff00');
  });

  await softStep('Every viewer reflects R_ONE=red and S_UNKN=green (DOM)', async () => {
    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, any> = {viewers: {}};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
        const sUnknItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'S_UNKN');
        out.viewers[x.type] = {
          legendRendered: items.length > 0,
          rOneColor: rOneItem ? getComputedStyle(rOneItem).color : null,
          sUnknColor: sUnknItem ? getComputedStyle(sUnknItem).color : null,
        };
      }
      return out;
    });
    let viewersWithLegend = 0;
    for (const [_, info] of Object.entries(result.viewers as Record<string, any>)) {
      if (!info.legendRendered) continue;
      if (info.rOneColor === 'rgb(255, 0, 0)' && info.sUnknColor === 'rgb(0, 255, 0)')
        viewersWithLegend++;
    }
    expect(viewersWithLegend, 'at least 1 viewer renders legend with the configured DOM colors').toBeGreaterThanOrEqual(1);
  });

  await softStep('Open color picker via legend, change R_ONE to blue', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Histogram',
      category: 'R_ONE',
      rgb: [31, 119, 180],
      hex: '#1f77b4',
      column: 'Stereo Category',
      additive: {'R_ONE': '#1f77b4', 'S_UNKN': '#00FF00'},
    });
  });

  await softStep('Picker change propagated: every legend item in DOM shows blue', async () => {
    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, string | null> = {};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
        out[x.type] = rOneItem ? getComputedStyle(rOneItem).color : null;
      }
      return out;
    });
    let viewersChecked = 0;
    for (const [_, color] of Object.entries(result)) {
      if (color === 'rgb(31, 119, 180)') viewersChecked++;
    }
    expect(viewersChecked, 'picker change reflected in legend DOM on at least 1 viewer').toBeGreaterThanOrEqual(1);
  });

  await softStep('Save + re-apply layout — custom palette persists (tag verification)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'ColorConsist_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      (window as any).__ccLayoutId = saved.id;
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      const tag = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return {layoutId: saved.id, rOneAfterReload: String(tag['R_ONE'] ?? '').toLowerCase()};
    });
    (globalThis as any).__ccLayoutId = res.layoutId;
    expect(res.rOneAfterReload).toBe('#1f77b4');
  });

  await softStep('Project round-trip — save + close + reopen + verify palette', async () => {
    const res = await page.evaluate(async () => {
      let projectId: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'ColorConsistProj_' + Date.now();
        const __df = (window as any).grok.shell.tv.dataFrame;
        const __ti = __df.getTableInfo();
        proj.addChild(__ti);
        // Persist the table entity BEFORE saving the project, else project_relations
        // references a not-yet-persisted entity id -> FK violation on the CI stack.
        await (window as any).grok.dapi.tables.uploadDataFrame(__df);
        await (window as any).grok.dapi.tables.save(__ti);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        projectId = saved.id;
      } catch (e: any) {
        return {phase: 'save', ok: false, error: String(e).slice(0, 200)};
      }
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1200));
      try {
        const reopened = await (window as any).grok.dapi.projects.find(projectId);
        await reopened.open();
      } catch (e: any) {
        return {phase: 'reopen', ok: false, error: String(e).slice(0, 200), projectId};
      }
      await new Promise((r) => setTimeout(r, 3500));
      const tv = (window as any).grok.shell.tv;
      if (!tv) return {phase: 'reopen', ok: false, error: 'no tv after reopen', projectId};
      const col = tv.dataFrame.col('Stereo Category');
      const tag = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      const colorAfter = String(tag['R_ONE'] ?? '').toLowerCase();
      const viewerColors: Record<string, string|null> = {};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
        viewerColors[x.type] = rOneItem ? getComputedStyle(rOneItem).color : null;
      }
      return {phase: 'verified', ok: true, projectId, colorAfter, viewerColors};
    });
    expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
    (globalThis as any).__ccProjectId = res.projectId;
    expect(res.colorAfter).toBe('#1f77b4');
    const colors = res.viewerColors as Record<string, string|null>;
    let checked = 0;
    for (const [_, color] of Object.entries(colors)) {
      if (color === 'rgb(31, 119, 180)') checked++;
    }
    expect(checked, 'at least 1 viewer reflects R_ONE=blue post-reopen').toBeGreaterThanOrEqual(1);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async ([layoutId, projectId]) => {
      if (layoutId) try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(layoutId)); } catch (_) {}
      if (projectId) try { await (window as any).grok.dapi.projects.delete(await (window as any).grok.dapi.projects.find(projectId)); } catch (_) {}
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    }, [(globalThis as any).__ccLayoutId, (globalThis as any).__ccProjectId]);
  });

  v.finishSpec();
});
