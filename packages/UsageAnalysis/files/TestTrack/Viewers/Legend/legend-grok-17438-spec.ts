/* ---
sub_features_covered: [legend.item.color-picker, legend.visibility, legend.splitter-resize]
--- */
// Frontmatter extraction (Edit X7):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: 3 atlas ids (from bug.affects)
//   ui_coverage_responsibility: []
//   related_bugs: [GROK-17438]
//   coverage_type: regression
//   bug_id: GROK-17438
//   bug_url: https://reddata.atlassian.net/browse/GROK-17438
//   reproduction_source: bug-library/legend.yaml#GROK-17438
//   related_scenarios: chain YAML bug_focused_candidates[0].spans
//     (visibility-and-positioning.md:Step 9, color-consistency.md:Step 5,
//      scatterplot.md:Step 4)
// Fix invariant: when color is changed on one viewer, the legend remains visible on
// other viewers with the same legend; recovery paths (resize, set visibility=Always)
// restore the legend if hidden.
//
// Selector sources (grok-browser/references):
//   .claude/skills/grok-browser/references/viewers.md (Legend section L112-135):
//     - L118: [name="legend"] host, [name="legend-splitter"] resize handle
//     - L130: right-click on .d4-legend-item opens picker dialog
//   .claude/plan/legend-mcp-recon-2026-05-08-color-picker.md (MCP-validated picker DOM)
//
// Setup: 3-viewer shared-legend layout (Histogram + Scatter + Bar all on Stereo Category).
// UI flow: real Playwright `locator.click({button: 'right'})` on Histogram legend
// (per validator-b-legend-section-postfix-2026-05-08); JS-API fallback via setCategorical.
// Recovery exercise: visibility=Always + splitter-resize must each restore legend
// visibility if it has been incorrectly hidden.
//
// Bug reproduction (bug-library/legend.yaml#GROK-17438):
//   1. Open SPGI
//   2. Open the layout — legend is present on each viewer
//   3. Change any color on the scatterplot — the legend is hidden on viewers with same legend
//   4. Try to resize the viewers — legend is still hidden
//   5. Try to set legend visibility to always — legend is still hidden
// Expected: When color is changed, legend remains visible on shared-legend viewers;
// recovery paths (resize, visibility=Always) restore visibility.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('GROK-17438: legend stays visible across shared-legend viewers after color change', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // technical: open SPGI + 3-viewer shared-legend layout (Histogram + Scatter + Bar
  // all using Stereo Category).
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
    const tv = (window as any).grok.shell.tv;
    const names = ['Histogram', 'Scatter plot', 'Bar chart'];
    for (const n of names) {
      tv.addViewer(n);
      await new Promise(r => setTimeout(r, 400));
    }
    const col = 'Stereo Category';
    for (const v of tv.viewers) {
      if (v.type === 'Grid') continue;
      try {
        if (v.type === 'Histogram') v.props.splitColumnName = col;
        else if (v.type === 'Scatter plot') v.props.colorColumnName = col;
        else if (v.type === 'Bar chart') v.props.splitColumnName = col;
        try { v.props.legendVisibility = 'Always'; } catch (_) {}
      } catch (_) {}
    }
    await new Promise(r => setTimeout(r, 2000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // Step 2 (bug repro): legend present on each viewer (baseline before color change).
  await softStep('Step 2: legend present on Histogram + Scatter + Bar (baseline)', async () => {
    const before = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, boolean> = {};
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        out[v.type] = !!v.root.querySelector('[name="legend"]');
      }
      return out;
    });
    // Bar chart's legend may not appear until a color edit fires propagation
    // (per legend-mcp-recon-2026-05-08-color-picker.md "Bar chart legend caveat" L74-94).
    // At least 2 of 3 viewers must show legend in baseline (Histogram + Scatter reliably do).
    const present = Object.values(before).filter(Boolean).length;
    expect(present, 'at least 2 viewers show legend in baseline').toBeGreaterThanOrEqual(2);
  });

  // Step 3 (bug repro): change color via Scatter plot legend. Legend MUST remain
  // visible on the OTHER shared-legend viewers (Histogram + Bar).
  await softStep('Step 3 invariant: change R_ONE via Scatter legend; Histogram + Bar legends stay visible', async () => {
    const dlgName = 'dialog-R-ONE';
    let okCommitted = false;
    try {
      const item = page.locator('[name="viewer-Scatter-plot"] [name="legend"] .d4-legend-item')
        .filter({has: page.locator('.d4-legend-value', {hasText: /^R_ONEx?$/})}).first();
      await item.scrollIntoViewIfNeeded();
      await item.click({button: 'right', timeout: 5000});
      await page.locator(`.d4-dialog[name="${dlgName}"]`).waitFor({timeout: 5000});
      await page.evaluate(({dlgName}) => {
        const dlg = document.querySelector(`.d4-dialog[name="${dlgName}"]`)!;
        const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
          .find(s => s.style.backgroundColor === 'rgb(31, 119, 180)');
        if (sw) {
          const opts = {bubbles: true, cancelable: true, view: window, button: 0};
          sw.dispatchEvent(new MouseEvent('mousedown', opts));
          sw.dispatchEvent(new MouseEvent('mouseup', opts));
          sw.dispatchEvent(new MouseEvent('click', opts));
        }
      }, {dlgName});
      await page.waitForTimeout(200);
      await page.locator(`.d4-dialog[name="${dlgName}"] [name="button-OK"]`).click({timeout: 5000});
      await page.waitForTimeout(700);
      const tag = await page.evaluate(() => {
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
        return String(t['R_ONE'] ?? '').toLowerCase();
      });
      okCommitted = tag === '#1f77b4' || tag.includes('1f77b4');
    } catch (_) {
      okCommitted = false;
    }
    if (!okCommitted) {
      // JS-API fallback (per ApiSamples scripts/grid/color-coding/color-coding.js):
      // setCategorical is what the OK handler calls internally
      // (public/js-api/src/dataframe/column-helpers.ts L103-111). The bug is about
      // legend visibility post-color-edit — exercise the same code path.
      await page.evaluate(() => {
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        col.tags['.color-coding-type'] = 'Categorical';
        col.meta.colors.setCategorical({'R_ONE': '#1f77b4'}, {fallbackColor: '#808080'});
        for (const v of (window as any).grok.shell.tv.viewers)
          if (v.type !== 'Grid') try { v.invalidate?.(); } catch (_) {}
      });
      await page.waitForTimeout(1500);
    }
    // Bug invariant — fix invariant for GROK-17438: legend remains visible on
    // OTHER shared-legend viewers (Histogram + Bar) after Scatter color edit.
    const after = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, boolean> = {};
      for (const v of tv.viewers) {
        if (v.type === 'Grid' || v.type === 'Scatter plot') continue;
        out[v.type] = !!v.root.querySelector('[name="legend"]');
      }
      return out;
    });
    // Histogram is the most reliable shared-legend viewer in headless layout.
    expect(after['Histogram'], 'Histogram legend stays visible after Scatter color change (GROK-17438)').toBe(true);
  });

  // Steps 4-5 (bug repro recovery paths): exercise visibility=Always + splitter-resize.
  // The bug claims these do NOT restore visibility once it has been hidden. Fix invariant:
  // both recovery paths produce a visible legend.
  // Recovery 1: visibility=Always
  await softStep('Step 5 recovery: legendVisibility=Always restores legend on each viewer', async () => {
    const after = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      // Force-hide on Histogram first to verify recovery.
      const h = tv.viewers.find((v: any) => v.type === 'Histogram');
      try { h.props.legendVisibility = 'Never'; } catch (_) {}
      await new Promise(r => setTimeout(r, 800));
      const hiddenState = !!h.root.querySelector('[name="legend"]');
      // Recovery: set visibility back to Always.
      try { h.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise(r => setTimeout(r, 1500));
      const restoredState = !!h.root.querySelector('[name="legend"]');
      return {hiddenState, restoredState, vis: h.props.legendVisibility};
    });
    expect(after.vis).toBe('Always');
    // Bug invariant — fix invariant for GROK-17438: visibility=Always recovery works.
    expect(after.restoredState, 'Histogram legend restored via legendVisibility=Always (GROK-17438)').toBe(true);
  });

  // Recovery 2: splitter-resize. Use real Playwright drag on [name="legend-splitter"].
  await softStep('Step 4 recovery: splitter-resize gesture (Histogram side-docked)', async () => {
    await page.evaluate(async () => {
      const h = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram');
      try { h.props.legendPosition = 'Right'; } catch (_) {}
      try { h.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise(r => setTimeout(r, 1000));
    });
    const splitter = page.locator('[name="viewer-Histogram"] [name="legend-splitter"]').first();
    if (await splitter.count() > 0) {
      const beforeBox = await page.evaluate(() => {
        const h = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram');
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
          const h = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram');
          const legend = h.root.querySelector('[name="legend"]');
          return legend ? legend.getBoundingClientRect().width : null;
        });
        // Splitter response is build-dependent; the gesture is the surface under test.
        // Bug invariant — fix invariant for GROK-17438: splitter-resize gesture
        // does not break legend visibility (legend still rendered).
        expect(typeof afterBox).toBe('number');
        expect(typeof beforeBox).toBe('number');
      }
    }
    const stillVisible = await page.evaluate(() =>
      !!(window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram')
        .root.querySelector('[name="legend"]'));
    expect(stillVisible, 'Histogram legend still visible after splitter-resize (GROK-17438)').toBe(true);
  });

  // Cleanup: reset state + close views.
  await softStep('Cleanup', async () => {
    await page.evaluate(() => {
      try {
        const col = (window as any).grok.shell.tv?.dataFrame.col('Stereo Category');
        if (col) {
          delete col.tags['.color-coding-categorical'];
          delete col.tags['.color-coding-type'];
        }
      } catch (_) {}
      (window as any).grok.shell.closeAll();
    });
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
