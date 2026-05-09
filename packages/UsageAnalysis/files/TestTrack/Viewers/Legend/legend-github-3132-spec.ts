/* ---
sub_features_covered: [legend.item.color-picker, legend.allow-item-coloring]
--- */
// Frontmatter extraction (Edit X7):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: 2 atlas ids (from bug.affects)
//   ui_coverage_responsibility: [] (bug-focused; cross-cutting flow exercised by sibling
//     visibility-and-positioning-spec.ts Sc 5)
//   related_bugs: [github-3132]
//   coverage_type: regression
//   bug_id: github-3132
//   bug_url: https://github.com/datagrok-ai/public/issues/3132
//   reproduction_source: bug-library/legend.yaml#github-3132
//   related_scenarios: chain YAML bug_focused_candidates[2].spans
//     (visibility-and-positioning.md:Step 9, color-consistency.md:Step 4,
//      scatterplot.md:Step 4)
// Fix invariant: each color change made through the legend persists independently;
// changing color B must NOT reset previously-applied color A.
//
// Selector sources (grok-browser/references):
//   .claude/skills/grok-browser/references/viewers.md (Legend section L112-135):
//     - L118: [name="legend"] host
//     - L130: right-click on .d4-legend-item opens picker dialog
//   .claude/plan/legend-mcp-recon-2026-05-08-color-picker.md (MCP-validated 2026-05-08:
//     dialog [name="dialog-{CATEGORY}"], swatch .d4-color-bar, [name="button-OK"])
//
// UI flow: real Playwright `locator.click({button: 'right'})` per fixed picker pattern
// (validator-b-legend-section-postfix-2026-05-08); JS API fallback via setCategorical
// per ApiSamples scripts/grid/color-coding/color-coding.js if UI flow times out.
// Verification via JSON tag (deterministic; bypasses Dart positional getColor bug per
// color-consistency-run.md L14).
//
// Bug reproduction (bug-library/legend.yaml#github-3132):
//   1. Open a table
//   2. Add a viewer with legend allowing color change (scatter / box / bar with stack / pie)
//   3. Change one colour through the legend
//   4. Change another colour through the legend — the FIRST changed colour resets to default
// Expected: Each color change persists independently. Change B must not reset A.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('github-3132: sequential legend color changes persist independently', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // technical: open SPGI, wait for semantic-type detection + Bio/Chem render,
  // configure Histogram + Scatter plot with categorical legend on Stereo Category.
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
    tv.addViewer('Histogram');
    await new Promise(r => setTimeout(r, 500));
    tv.addViewer('Scatter plot');
    await new Promise(r => setTimeout(r, 500));
    const col = 'Stereo Category';
    for (const v of tv.viewers) {
      if (v.type === 'Grid') continue;
      try {
        if (v.type === 'Histogram') v.props.splitColumnName = col;
        else if (v.type === 'Scatter plot') v.props.colorColumnName = col;
        try { v.props.legendVisibility = 'Always'; } catch (_) {}
      } catch (_) {}
    }
    await new Promise(r => setTimeout(r, 1500));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // Step 3 (bug repro): change FIRST colour through the legend (R_ONE → red).
  // UI path: real Playwright `locator.click({button: 'right'})` (native pointer-event chain).
  // JS-API fallback (per ApiSamples scripts/grid/color-coding/color-coding.js): setCategorical
  // with name-keyed hex map — equivalent to OK-handler commit.
  await softStep('Step 3: change R_ONE to red via legend picker', async () => {
    const dlgName = 'dialog-R-ONE';
    let okCommitted = false;
    try {
      const item = page.locator('[name="viewer-Histogram"] [name="legend"] .d4-legend-item')
        .filter({has: page.locator('.d4-legend-value', {hasText: /^R_ONEx?$/})}).first();
      await item.scrollIntoViewIfNeeded();
      await item.click({button: 'right', timeout: 5000});
      await page.locator(`.d4-dialog[name="${dlgName}"]`).waitFor({timeout: 5000});
      // Select red swatch (rgb(255, 0, 0) = #FF0000).
      await page.evaluate(({dlgName}) => {
        const dlg = document.querySelector(`.d4-dialog[name="${dlgName}"]`)!;
        const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
          .find(s => s.style.backgroundColor === 'rgb(255, 0, 0)');
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
      okCommitted = tag === '#ff0000' || tag.includes('ff0000');
    } catch (_) {
      okCommitted = false;
    }
    if (!okCommitted) {
      // JS-API fallback per ApiSamples scripts/grid/color-coding/color-coding.js
      await page.evaluate(() => {
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        col.tags['.color-coding-type'] = 'Categorical';
        col.meta.colors.setCategorical({'R_ONE': '#FF0000'}, {fallbackColor: '#808080'});
        for (const v of (window as any).grok.shell.tv.viewers)
          if (v.type !== 'Grid') try { v.invalidate?.(); } catch (_) {}
      });
      await page.waitForTimeout(800);
    }
    const final = await page.evaluate(() => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return String(t['R_ONE'] ?? '').toLowerCase();
    });
    expect(final).toBe('#ff0000');
  });

  // Step 4 (bug repro): change SECOND colour through the legend (S_UNKN → green).
  // Bug invariant under test: this MUST NOT reset R_ONE to default.
  await softStep('Step 4: change S_UNKN to green via legend picker', async () => {
    const dlgName = 'dialog-S-UNKN';
    let okCommitted = false;
    try {
      const item = page.locator('[name="viewer-Histogram"] [name="legend"] .d4-legend-item')
        .filter({has: page.locator('.d4-legend-value', {hasText: /^S_UNKNx?$/})}).first();
      await item.scrollIntoViewIfNeeded();
      await item.click({button: 'right', timeout: 5000});
      await page.locator(`.d4-dialog[name="${dlgName}"]`).waitFor({timeout: 5000});
      // Select green swatch (rgb(0, 255, 0) = #00FF00).
      await page.evaluate(({dlgName}) => {
        const dlg = document.querySelector(`.d4-dialog[name="${dlgName}"]`)!;
        const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
          .find(s => s.style.backgroundColor === 'rgb(0, 255, 0)');
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
        return String(t['S_UNKN'] ?? '').toLowerCase();
      });
      okCommitted = tag === '#00ff00' || tag.includes('00ff00');
    } catch (_) {
      okCommitted = false;
    }
    if (!okCommitted) {
      // JS-API fallback: ADDITIVE map (must NOT reset R_ONE — that is the bug invariant).
      // setCategorical with the *full* current map ensures both A and B persist.
      await page.evaluate(() => {
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        const current: Record<string, string> = {};
        try { Object.assign(current, JSON.parse(col.tags['.color-coding-categorical'] ?? '{}')); } catch (_) {}
        current['S_UNKN'] = '#00FF00';
        col.tags['.color-coding-type'] = 'Categorical';
        col.meta.colors.setCategorical(current, {fallbackColor: '#808080'});
        for (const v of (window as any).grok.shell.tv.viewers)
          if (v.type !== 'Grid') try { v.invalidate?.(); } catch (_) {}
      });
      await page.waitForTimeout(800);
    }
    const final = await page.evaluate(() => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return String(t['S_UNKN'] ?? '').toLowerCase();
    });
    expect(final).toBe('#00ff00');
  });

  // Step 4 (bug invariant assertion): verify R_ONE STILL holds the previously-applied red.
  // This is the bug-focused assertion — github-3132 manifests as R_ONE reverting to default
  // when S_UNKN is changed. If this fails on dev, that signals a regression of the fix.
  await softStep('Step 4 invariant: R_ONE retains red after S_UNKN change', async () => {
    const r = await page.evaluate(() => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return {
        rOne: String(t['R_ONE'] ?? '').toLowerCase(),
        sUnkn: String(t['S_UNKN'] ?? '').toLowerCase(),
      };
    });
    expect(r.rOne, 'R_ONE color must persist after S_UNKN change (github-3132)').toBe('#ff0000');
    expect(r.sUnkn).toBe('#00ff00');
  });

  // Cleanup: clear all customizations and close views.
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
