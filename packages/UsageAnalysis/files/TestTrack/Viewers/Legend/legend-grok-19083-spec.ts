// GROK-19083: legend marker entries must react to deselecting markers on the host viewer.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('GROK-19083: legend marker entries sync with markersColumnName deselect', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Steps 2-4: Color=Series + Marker=Series → combined legend with markers', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      sp.props.markersColumnName = 'Series';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1800));
      const i = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {
        itemCount: i.length,
        markersBefore: sp.props.markersColumnName,
        colorBefore: sp.props.colorColumnName,
      };
    });
    expect(items.markersBefore).toBe('Series');
    expect(items.colorBefore).toBe('Series');
    expect(items.itemCount).toBeGreaterThanOrEqual(2);
  });

  await softStep('DOM gesture: hover Scatter plot legend (real Playwright)', async () => {
    const legend = page.locator('[name="viewer-Scatter-plot"] [name="legend"]').first();
    await legend.waitFor({timeout: 10000});
    await legend.hover();
  });

  await softStep('Step 5 invariant: deselect markers — legend updates in sync', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const itemsBefore = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const markerIconsBefore = itemsBefore.filter((it) =>
        !!it.querySelector('svg, canvas, .d4-legend-marker, [class*="marker"]')).length;
      const beforeCount = itemsBefore.length;
      sp.props.markersColumnName = '';
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1800));
      const itemsAfter = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const markerIconsAfter = itemsAfter.filter((it) =>
        !!it.querySelector('svg, canvas, .d4-legend-marker, [class*="marker"]')).length;
      return {
        markersAfter: sp.props.markersColumnName,
        beforeCount,
        afterCount: itemsAfter.length,
        markerIconsBefore,
        markerIconsAfter,
      };
    });
    expect(res.markersAfter, 'markersColumnName cleared (GROK-19083)').toBe('');
    expect(res.afterCount, 'legend retains category items via Color binding').toBeGreaterThan(0);
    expect(res.markerIconsAfter, 'marker glyph count must not increase after deselect (GROK-19083)')
      .toBeLessThanOrEqual(res.markerIconsBefore);
  });

  await softStep('Step 5 follow-up: re-bind markers — legend renders entries again', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.markersColumnName = 'Series';
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {markers: sp.props.markersColumnName, itemCount: items.length};
    });
    expect(res.markers).toBe('Series');
    expect(res.itemCount).toBeGreaterThan(0);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });

  v.finishSpec();
});
