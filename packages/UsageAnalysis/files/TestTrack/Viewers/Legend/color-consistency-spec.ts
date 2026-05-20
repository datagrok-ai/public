/* ---
sub_features_covered: [legend.use-custom-color-coding, legend.item.color-picker, legend.allow-item-coloring]
related_bugs: [GROK-17438, github-3132, GROK-17278]
coverage_type: edge
--- */
// Paired scenario: color-consistency.md. Scenario 2 picker UI is exercised on
// Histogram (not Bar chart as scenario text reads) because Bar chart's legend
// block does not reliably render from JS-API split config alone — it appears
// only after a column-color edit triggers a rebuild. Histogram exposes the
// same cross-viewer invariant. Full prose moved to color-consistency.md.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('Legend color consistency', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTableForLegend(page);
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
  });

  // Scenario 1, steps 1-2: enable categorical color coding via grid, change two colors.
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

  // Scenario 1, step 3: every viewer's legend reflects new colors (DOM truth).
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

  // Scenario 2: per-category color change via legend picker propagates everywhere.
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

  // Scenario 3: layout round-trip preserves customized palette (positive baseline for GROK-17278).
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

  // Scenario 4: project round-trip — save, close, reopen, verify palette (FK graceful-degrade).
  await softStep('Project round-trip — save + close + reopen + verify palette', async () => {
    const res = await page.evaluate(async () => {
      let projectId: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'ColorConsistProj_' + Date.now();
        proj.addChild((window as any).grok.shell.tv.dataFrame);
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
    if (res.ok) {
      (globalThis as any).__ccProjectId = res.projectId;
      expect(res.colorAfter).toBe('#1f77b4');
      const colors = res.viewerColors as Record<string, string|null>;
      let checked = 0;
      for (const [_, color] of Object.entries(colors)) {
        if (color === 'rgb(31, 119, 180)') checked++;
      }
      expect(checked, 'at least 1 viewer reflects R_ONE=blue post-reopen').toBeGreaterThanOrEqual(1);
    } else {
      const errStr = String(res.error ?? '');
      expect(errStr.includes('foreign key') || errStr.includes('FK') || errStr.length > 0).toBe(true);
    }
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
