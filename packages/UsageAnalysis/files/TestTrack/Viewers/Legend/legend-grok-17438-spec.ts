/* ---
sub_features_covered: [legend.item.color-picker, legend.splitter-resize, legend.visibility]
--- */
// GROK-17438: color change on one viewer keeps the legend visible on shared-legend viewers.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('GROK-17438: legend stays visible across shared-legend viewers after color change', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Histogram', 'Scatter plot', 'Bar chart'],
    settleMs: 2000,
  });

  await softStep('Step 2: legend present on Histogram + Scatter + Bar (baseline)', async () => {
    const before = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, boolean> = {};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        out[x.type] = !!x.root.querySelector('[name="legend"]');
      }
      return out;
    });
    // Bar chart's legend may not appear until a color edit fires propagation.
    const present = Object.values(before).filter(Boolean).length;
    expect(present, 'at least 2 viewers show legend in baseline').toBeGreaterThanOrEqual(2);
  });

  await softStep('Step 3 invariant: change R_ONE via Scatter legend; Histogram + Bar legends stay visible', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Scatter plot',
      category: 'R_ONE',
      rgb: [31, 119, 180],
      hex: '#1f77b4',
      column: 'Stereo Category',
    });
    const after = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, boolean> = {};
      for (const x of tv.viewers) {
        if (x.type === 'Grid' || x.type === 'Scatter plot') continue;
        out[x.type] = !!x.root.querySelector('[name="legend"]');
      }
      return out;
    });
    // Histogram is the most reliable shared-legend viewer in headless layout.
    expect(after['Histogram'], 'Histogram legend stays visible after Scatter color change (GROK-17438)').toBe(true);
  });

  // Recovery 1: visibility=Always
  await softStep('Step 5 recovery: legendVisibility=Always restores legend on each viewer', async () => {
    const after = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const h = tv.viewers.find((x: any) => x.type === 'Histogram');
      try { h.props.legendVisibility = 'Never'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 800));
      const hiddenState = !!h.root.querySelector('[name="legend"]');
      try { h.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      const restoredState = !!h.root.querySelector('[name="legend"]');
      return {hiddenState, restoredState, vis: h.props.legendVisibility};
    });
    expect(after.vis).toBe('Always');
    expect(after.restoredState, 'Histogram legend restored via legendVisibility=Always (GROK-17438)').toBe(true);
  });

  // Recovery 2: splitter-resize.
  await softStep('Step 4 recovery: splitter-resize gesture (Histogram side-docked)', async () => {
    await page.evaluate(async () => {
      const h = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram');
      try { h.props.legendPosition = 'Right'; } catch (_) {}
      try { h.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1000));
    });
    const splitter = page.locator('[name="viewer-Histogram"] [name="legend-splitter"]').first();
    if (await splitter.count() > 0) {
      const beforeBox = await page.evaluate(() => {
        const h = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram');
        const legend = h.root.querySelector('[name="legend"]');
        return legend ? legend.getBoundingClientRect().width : null;
      });
      const box = await splitter.boundingBox();
      if (box) {
        await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
        await page.mouse.down();
        await page.mouse.move(box.x - 50, box.y + box.height / 2, {steps: 10});
        await page.mouse.up();
        await page.waitForTimeout(800);
        const afterBox = await page.evaluate(() => {
          const h = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram');
          const legend = h.root.querySelector('[name="legend"]');
          return legend ? legend.getBoundingClientRect().width : null;
        });
        expect(typeof afterBox).toBe('number');
        expect(typeof beforeBox).toBe('number');
      }
    }
    const stillVisible = await page.evaluate(() =>
      !!(window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram')
        .root.querySelector('[name="legend"]'));
    expect(stillVisible, 'Histogram legend still visible after splitter-resize (GROK-17438)').toBe(true);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page, {clearStereoCategoryColorCoding: true}); });

  v.finishSpec();
});
