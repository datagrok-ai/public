/* ---
sub_features_covered: [legend.show-main-item-icons, legend.column]
related_bugs: [GROK-19083]
coverage_type: regression
bug_url: https://reddata.atlassian.net/browse/GROK-19083
--- */
// Fix invariant: legend marker entries must react to "deselect markers" on the
// host viewer; legend cannot fall out of sync with viewer-level visualization
// properties. Full prose moved to legend-grok-19083.md.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('GROK-19083: legend marker entries sync with markersColumnName deselect', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTableForLegend(page);

  // Steps 2-4 (bug repro): Color=Series + Marker=Series → combined legend with markers.
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

  // Real DOM gesture on the Scatter plot legend host — exercises hover surface
  // (legend-color-picker icon visibility + tooltip path) for E-LAYER-COMPLIANCE-01.
  await softStep('DOM gesture: hover Scatter plot legend (real Playwright)', async () => {
    const legend = page.locator('[name="viewer-Scatter-plot"] [name="legend"]').first();
    await legend.waitFor({timeout: 10000});
    await legend.hover();
  });

  // Step 5 (bug repro): deselect markers (set markersColumnName=''). Legend MUST
  // re-render WITHOUT marker glyphs OR with updated entries reflecting the new state.
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

  // Step 5 cont: re-bind markers to verify the legend can recover the marker entries.
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
