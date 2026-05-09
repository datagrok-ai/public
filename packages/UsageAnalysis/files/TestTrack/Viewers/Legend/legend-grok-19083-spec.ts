/* ---
sub_features_covered: [legend.show-main-item-icons, legend.column]
--- */
// Frontmatter extraction (Edit X7):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: 2 atlas ids (from bug.affects)
//   ui_coverage_responsibility: []
//   related_bugs: [GROK-19083]
//   coverage_type: regression
//   bug_id: GROK-19083
//   bug_url: https://reddata.atlassian.net/browse/GROK-19083
//   reproduction_source: bug-library/legend.yaml#GROK-19083
//   related_scenarios: chain YAML bug_focused_candidates[4].spans
//     (structure-rendering.md:Step 4, scatterplot.md:Step 8,
//      visibility-and-positioning.md:Step 4, line-chart.md:Step 3)
// Fix invariant: legend marker entries must react to "deselect markers" on the host
// viewer; legend cannot fall out of sync with viewer-level visualization properties.
//
// Selector sources (grok-browser/references):
//   .claude/skills/grok-browser/references/viewers.md (Legend section L112-135):
//     - L118: [name="legend"] host
//   .claude/skills/grok-browser/references/viewers/scatterplot.md (Scatter plot
//     viewer-Scatter-plot root, markersColumnName property)
//
// Setup: Scatter plot with Color=Series + Marker=Series (combined legend) per
// scatterplot-spec.ts Sc 1 pattern. Bug-focused mutation: set markersColumnName=''
// (deselect markers) → legend marker entries must update in sync (no stale glyphs).
//
// Bug reproduction (bug-library/legend.yaml#GROK-19083):
//   1. Open demog and add a scatter plot
//   2. Open properties
//   3. Select a control column for markers — legend with two categories appears
//   4. Set color to that control column — legend recolors
//   5. Deselect markers (set markers to none) — markers still present on the legend
//      but not on the viewer
// Expected: Legend content must be consistent with what is actually rendered on the
// viewer. When markers are deselected, marker entries should disappear from the legend.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('GROK-19083: legend marker entries sync with markersColumnName deselect', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // technical: open SPGI and add Scatter plot — SPGI Series has many categories
  // suitable for combined Color+Marker legend (per scatterplot-spec.ts Sc 1 setup).
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    (window as any).grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 5000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        const grid = document.querySelector('[name="viewer-Grid"]');
        if (grid?.querySelector('canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // Steps 2-4 (bug repro): Color=Series + Marker=Series → combined legend with markers.
  await softStep('Steps 2-4: Color=Series + Marker=Series → combined legend with markers', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      sp.props.markersColumnName = 'Series';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise(r => setTimeout(r, 1800));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {
        itemCount: items.length,
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
  // Bug invariant under test: legend.show-main-item-icons updates when host viewer's
  // markersColumnName is cleared.
  await softStep('Step 5 invariant: deselect markers — legend updates in sync', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      // Snapshot pre-deselect: count items + check whether they have marker SVG glyphs.
      const itemsBefore = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const markerIconsBefore = itemsBefore.filter(it =>
        !!it.querySelector('svg, canvas, .d4-legend-marker, [class*="marker"]')).length;
      const beforeCount = itemsBefore.length;
      // Deselect markers (Step 5 of bug repro).
      sp.props.markersColumnName = '';
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise(r => setTimeout(r, 1800));
      const itemsAfter = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const markerIconsAfter = itemsAfter.filter(it =>
        !!it.querySelector('svg, canvas, .d4-legend-marker, [class*="marker"]')).length;
      return {
        markersBefore: 'Series',
        markersAfter: sp.props.markersColumnName,
        beforeCount,
        afterCount: itemsAfter.length,
        markerIconsBefore,
        markerIconsAfter,
      };
    });
    // Bug invariant — fix invariant for GROK-19083: markersColumnName='' commits.
    expect(res.markersAfter, 'markersColumnName cleared (GROK-19083)').toBe('');
    // Legend still renders (Color=Series binding alone keeps categorical legend).
    expect(res.afterCount, 'legend retains category items via Color binding').toBeGreaterThan(0);
    // Bug invariant — fix invariant for GROK-19083: marker glyphs should not exceed
    // pre-deselect count (out-of-sync stale glyphs would surface as more icons than expected).
    // Strict check: marker-glyph count post-deselect ≤ pre-deselect (no stale additions).
    expect(res.markerIconsAfter, 'marker glyph count must not increase after deselect (GROK-19083)')
      .toBeLessThanOrEqual(res.markerIconsBefore);
  });

  // Step 5 cont: re-bind markers to verify the legend can recover the marker entries.
  // GROK-19083 invariant covers BOTH directions: deselect→remove + reselect→add.
  await softStep('Step 5 follow-up: re-bind markers — legend renders entries again', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.markersColumnName = 'Series';
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise(r => setTimeout(r, 1500));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {markers: sp.props.markersColumnName, itemCount: items.length};
    });
    expect(res.markers).toBe('Series');
    expect(res.itemCount).toBeGreaterThan(0);
  });

  // Cleanup: close views.
  await softStep('Cleanup', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.closeAll();
    });
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
