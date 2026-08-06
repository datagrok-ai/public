/* ---
realizes: [boxplot.cp.value-category-axes-persist, boxplot.int.category-sets-marker-color]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function setBpProp(page: Page, prop: string, value: any, settleMs = 900): Promise<void> {
  await page.evaluate(({p, val}) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    bp.props[p] = val;
  }, {p: prop, val: value});
  await page.waitForTimeout(settleMs);
}

// bp.viewport is a Rect getter: top/bottom are the visible value bounds, height the span.
async function viewportRect(page: Page): Promise<{top: number; bottom: number; height: number}> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const vp = bp.viewport;
    return {top: vp.top, bottom: vp.bottom, height: vp.height};
  });
}

// Handle centers of the vertical value-axis range slider (svg[type="range-slider"]
// with height > width). A freshly (re)created box plot — e.g. after loadLayout —
// keeps STALE handle geometry (circles collapsed near the slider origin) until
// pointer traffic over the slider strip forces a re-layout; wiggle the trusted
// pointer along the strip between polls until the handles separate.
async function verticalSliderHandles(
  page: Page,
): Promise<{top: {x: number; y: number}; bottom: {x: number; y: number}}> {
  const readSlider = () => page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const slider = Array.from(bp.root.querySelectorAll('svg[type="range-slider"]'))
      .find((s: any) => {
        const r = s.getBoundingClientRect();
        return r.height > r.width;
      }) as SVGElement | undefined;
    if (!slider) return null;
    const sr = slider.getBoundingClientRect();
    const circles = Array.from(slider.querySelectorAll('circle'));
    if (circles.length < 2) return null;
    const centers = circles.slice(0, 2)
      .map((c) => {
        const r = c.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      })
      .sort((a, b) => a.y - b.y);
    return {rect: {x: sr.x, y: sr.y, w: sr.width, h: sr.height},
      top: centers[0], bottom: centers[1]};
  });
  for (let i = 0; i < 20; i++) {
    const s = await readSlider();
    if (s && s.bottom.y - s.top.y > 50) return {top: s.top, bottom: s.bottom};
    if (s) {
      const fy = 0.15 + (i % 5) * 0.15;
      await page.mouse.move(s.rect.x + s.rect.w / 2, s.rect.y + s.rect.h * fy, {steps: 4});
    }
    await page.waitForTimeout(300);
  }
  throw new Error('vertical range slider handles did not lay out to a usable span');
}

// The group-comparison icon is hover-gated: it appears only after trusted
// pointer traffic over the canvas top-left reveal zone.
async function revealIcon(page: Page, iconName: string): Promise<{x: number; y: number}> {
  const origin = await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y};
  });
  for (const [dx, dy] of [[35, 15], [40, 17], [30, 14], [45, 16]]) {
    await page.mouse.move(origin.x + dx, origin.y + dy);
    await page.waitForTimeout(150);
  }
  return page.evaluate((name) => {
    const el = document.querySelector(`[name="${name}"]`) as HTMLElement;
    const r = el.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, iconName);
}

// Excludes viewport (asserted separately, on the project round-trip only).
async function readLadder(page: Page): Promise<Record<string, any>> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const p = bp.props;
    return {
      valueColumnName: p.valueColumnName,
      category1ColumnName: p.category1ColumnName,
      category2ColumnName: p.category2ColumnName,
      showMinorCategories: p.showMinorCategories,
      showAllCategories: p.showAllCategories,
      markerColorColumnName: p.markerColorColumnName,
      invertColorScheme: p.invertColorScheme,
      colorMin: p.colorMin,
      colorMax: p.colorMax,
      axisType: p.axisType,
      invertYAxis: p.invertYAxis,
      plotStyle: p.plotStyle,
    };
  });
}

test('Box Plot settings ladder and persistence round-trip', async ({page}) => {
  test.setTimeout(600_000);

  // Collected run-wide; the log-axis step reads the delta across its own boundary.
  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await page.waitForTimeout(1500);

  // ==================================================================
  // Scenario 1: Pre-ladder datetime axis-type gate probe (GROK-20395)
  // ==================================================================
  await softStep('[anchor: PRE-LADDER] Scenario 1: datetime Value disables the property-panel Axis Type (GROK-20395)', async () => {
    await setBpProp(page, 'valueColumnName', 'STARTED', 1200);
    // Route the box plot into the Context Panel so its property grid renders.
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      grok.shell.o = bp;
    });
    // The prop-axis-type row sits in a collapsed (display:none) property category,
    // so a visibility wait times out; poll the row's presence and read its
    // opacity via getComputedStyle instead.
    const axisTypeOpacity = await page.evaluate(async () => {
      for (let i = 0; i < 40; i++) {
        const row = document.querySelector('.property-grid tr[name="prop-axis-type"]') as HTMLElement | null;
        if (row) return getComputedStyle(row).opacity;
        await new Promise((r) => setTimeout(r, 250));
      }
      return null;
    });
    console.log('Scenario 1 Axis Type row opacity (datetime value):', axisTypeOpacity);
    expect(axisTypeOpacity).not.toBeNull();
    expect(parseFloat(axisTypeOpacity as string)).toBeLessThan(1);
    await setBpProp(page, 'valueColumnName', 'AGE', 1000);
    expect(await bpProp(page, 'valueColumnName')).toBe('AGE');
  });

  // ==================================================================
  // Scenario 2: Settings ladder — accumulate, do not revert
  // ==================================================================

  await softStep('[anchor: Step 4] Scenario 2 Step 3: Category 1 = SEX auto-sets Marker Color to SEX (category-sets-marker-color)', async () => {
    // Marker Color starts TRACKING the auto-selected Category 1 — do NOT clear it
    // first: an explicit clear to '' de-links the tracking and suppresses the
    // auto-follow side effect.
    expect(await bpProp(page, 'markerColorColumnName')).toBe(await bpProp(page, 'category1ColumnName'));
    await setBpProp(page, 'category1ColumnName', 'SEX', 1500);
    expect(await bpProp(page, 'category1ColumnName')).toBe('SEX');
    expect(await bpProp(page, 'markerColorColumnName')).toBe('SEX');
  });

  await softStep('Scenario 2 Step 4-5: explicit Marker Color HEIGHT, invert scheme, color min/max', async () => {
    await setBpProp(page, 'markerColorColumnName', 'HEIGHT', 800);
    await setBpProp(page, 'invertColorScheme', true, 500);
    await setBpProp(page, 'colorMin', 20, 500);
    await setBpProp(page, 'colorMax', 80, 700);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('HEIGHT');
    expect(await bpProp(page, 'invertColorScheme')).toBe(true);
    expect(await bpProp(page, 'colorMin')).toBe(20);
    expect(await bpProp(page, 'colorMax')).toBe(80);
  });

  await softStep('Scenario 2 Step 6-8: Category 2 = RACE, Show Minor + Show All Categories', async () => {
    await setBpProp(page, 'category2ColumnName', 'RACE', 900);
    await setBpProp(page, 'showMinorCategories', true, 500);
    await setBpProp(page, 'showAllCategories', true, 700);
    expect(await bpProp(page, 'category2ColumnName')).toBe('RACE');
    expect(await bpProp(page, 'showMinorCategories')).toBe(true);
    expect(await bpProp(page, 'showAllCategories')).toBe(true);
  });

  await softStep('[anchor: Step 5] Scenario 2 Step 10: Value=WEIGHT keeps Marker Color HEIGHT + color scheme (GROK-18876)', async () => {
    await setBpProp(page, 'valueColumnName', 'WEIGHT', 1200);
    expect(await bpProp(page, 'valueColumnName')).toBe('WEIGHT');
    expect(await bpProp(page, 'markerColorColumnName')).toBe('HEIGHT');
    expect(await bpProp(page, 'invertColorScheme')).toBe(true);
    expect(await bpProp(page, 'colorMin')).toBe(20);
    expect(await bpProp(page, 'colorMax')).toBe(80);
  });

  await softStep('Scenario 2 Step 11: Value Min 20, Value Max 60', async () => {
    await setBpProp(page, 'valueMin', 20, 500);
    await setBpProp(page, 'valueMax', 60, 700);
    expect(await bpProp(page, 'valueMin')).toBe(20);
    expect(await bpProp(page, 'valueMax')).toBe(60);
  });

  await softStep('[anchor: Step 8] Scenario 2 Step 13: Axis Type Log throws no console error, range within data bounds (GROK-18515, GROK-20397)', async () => {
    // Clear the manual bounds so the log axis governs the viewport directly.
    await setBpProp(page, 'valueMin', null, 400);
    await setBpProp(page, 'valueMax', null, 700);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'axisType', 'logarithmic', 1400);
    expect(await bpProp(page, 'axisType')).toBe('logarithmic');
    // grok.shell.warnings is not exposed to JS — the console/pageerror channels
    // are the no-throw signal.
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Step 13 console error delta:', JSON.stringify(errDelta), 'pageerror delta:', JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    const {top, bottom} = await viewportRect(page);
    const weightMax = await page.evaluate(() => {
      const col = grok.shell.t.col('WEIGHT');
      return col.stats.max;
    });
    console.log('Step 13 viewport top/bottom:', top, bottom, 'WEIGHT max:', weightMax);
    // Positive lower bound (log axis) and no stretch far beyond the data max.
    expect(top).toBeGreaterThan(0);
    expect(bottom).toBeLessThanOrEqual(weightMax * 1.1);
  });

  await softStep('Scenario 2 Step 14-15: Invert Y Axis on, Plot Style violin', async () => {
    await setBpProp(page, 'invertYAxis', true, 600);
    await setBpProp(page, 'plotStyle', 'violin', 800);
    expect(await bpProp(page, 'invertYAxis')).toBe(true);
    expect(await bpProp(page, 'plotStyle')).toBe('violin');
  });

  await softStep('Scenario 2 Step 17: range-slider drag narrows the visible value range', async () => {
    await v.waitForCanvasQuiet(page, 'Box plot');
    const before = await viewportRect(page);
    const handles = await verticalSliderHandles(page);
    await page.mouse.move(handles.top.x, handles.top.y);
    await page.mouse.down();
    const dropY = handles.top.y + (handles.bottom.y - handles.top.y) * 0.35;
    await page.mouse.move(handles.top.x, dropY, {steps: 16});
    await page.mouse.up();
    await page.waitForTimeout(1000);
    const after = await viewportRect(page);
    console.log('Step 17 viewport height before/after zoom:', before.height, after.height);
    // Real-margin bar: at least 10% narrower, not a pixel-level wobble.
    expect(after.height).toBeLessThan(before.height * 0.9);
  });

  await softStep('[anchor: Step 10] Scenario 2 Step 19: changing Marker Color to SEX leaves the zoomed viewport unchanged (GROK-20469)', async () => {
    const before = await viewportRect(page);
    await setBpProp(page, 'markerColorColumnName', 'SEX', 1200);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('SEX');
    const after = await viewportRect(page);
    console.log('Step 19 viewport before/after marker-color change:', JSON.stringify(before), JSON.stringify(after));
    expect(after.top).toBeCloseTo(before.top, 1);
    expect(after.bottom).toBeCloseTo(before.bottom, 1);
    expect(after.height).toBeCloseTo(before.height, 1);
  });

  // ==================================================================
  // Scenario 3: Group-comparison persistence tail
  // ==================================================================
  await softStep('[anchor: Step 11] Scenario 3: enable Group Comparison, pick control group, set Adjust By covariate', async () => {
    if (await bpProp(page, 'showGroupComparison') !== true) {
      const pt = await revealIcon(page, 'show-group-stats');
      await page.mouse.click(pt.x, pt.y);
      await page.waitForTimeout(1500);
    }
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);
    // The on-chart control-group <select> renders with an empty option list and a
    // display:none host (options are not populated headless) — the props path
    // (controlComparisons + controlGroup) is the sanctioned setter.
    await setBpProp(page, 'controlComparisons', true, 900);
    await setBpProp(page, 'controlGroup', 'F', 1200);
    expect(await bpProp(page, 'controlGroup')).toBe('F');
    // Adjust By has no simple on-chart set-selector — the props path is the
    // setter, the DOM caption/selector the observable.
    await setBpProp(page, 'covariateColumnName', 'HEIGHT', 1500);
    expect(await bpProp(page, 'covariateColumnName')).toBe('HEIGHT');
    const dom = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const caption = Array.from(root.querySelectorAll('.d4-column-selector-caption'))
        .some((c) => /Adjust by:/i.test(c.textContent ?? ''));
      const adjustSelector = Array.from(root.querySelectorAll('.d4-column-selector'))
        .map((s) => (s.textContent ?? '').replace(/\s+/g, ' ').trim());
      return {caption, hasHeightSelector: adjustSelector.some((t) => /Adjust by:\s*HEIGHT/i.test(t))};
    });
    console.log('Scenario 3 Adjust-by caption/selector:', JSON.stringify(dom));
    expect(dom.caption).toBe(true);
    expect(dom.hasHeightSelector).toBe(true);
  });

  // ==================================================================
  // Scenario 4: Layout round-trip (without viewport)
  // ==================================================================
  let savedLayout: any = null;
  let layoutSaved = false;
  await softStep('[anchor: LAYOUT-ROUND-TRIP] Scenario 4: layout round-trip restores the Box Plot ladder, Scatter re-arm absent', async () => {
    // Turn off the group-comparison extras before the snapshot so the saved
    // layout carries the plain accumulated ladder.
    await setBpProp(page, 'covariateColumnName', '', 700);
    await setBpProp(page, 'controlComparisons', false, 500);
    await setBpProp(page, 'showGroupComparison', false, 800);
    const ladderBefore = await readLadder(page);
    try {
      const saved = await page.evaluate(() => {
        const layout = grok.shell.tv.saveLayout();
        (window as any).__probeLayout = layout;
        return {name: layout.name ?? null, id: layout.id ?? null};
      });
      savedLayout = saved;
      layoutSaved = true;
    } catch (e) {
      throw new Error(`saveLayout failed: ${e}`);
    }
    await page.waitForTimeout(400);
    // Re-arm: prove the restore changes the viewer set, not a no-op.
    await page.evaluate(() => {
      const tv = grok.shell.tv;
      const bp = tv.viewers.find((x: any) => x.type === 'Box plot');
      if (bp) bp.close();
      tv.addViewer('Scatter plot');
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});
    await page.waitForTimeout(800);
    const armed = await page.evaluate(() => ({
      hasBox: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot'),
      hasScatter: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot'),
    }));
    expect(armed.hasBox).toBe(false);
    expect(armed.hasScatter).toBe(true);
    await page.evaluate(() => grok.shell.tv.loadLayout((window as any).__probeLayout));
    await page.waitForTimeout(3000);
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
    const restored = await page.evaluate(() => ({
      hasBox: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot'),
      hasScatter: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot'),
    }));
    expect(restored.hasBox).toBe(true);
    expect(restored.hasScatter).toBe(false);
    const ladderAfter = await readLadder(page);
    console.log('Scenario 4 ladder before/after layout round-trip:',
      JSON.stringify(ladderBefore), JSON.stringify(ladderAfter));
    expect(ladderAfter).toEqual(ladderBefore);
    expect(ladderAfter.valueColumnName).toBe('WEIGHT');
    expect(ladderAfter.category1ColumnName).toBe('SEX');
    expect(ladderAfter.category2ColumnName).toBe('RACE');
    expect(ladderAfter.axisType).toBe('logarithmic');
    expect(ladderAfter.plotStyle).toBe('violin');
    expect(ladderAfter.invertYAxis).toBe(true);
  });

  // ==================================================================
  // Scenario 5: Project round-trip (with viewport and group-comparison props)
  // ==================================================================
  let projectIds: {projectId: string} | null = null;
  const projectName = 'zz-boxplot-ladder-' + Date.now();
  try {
    await softStep('Scenario 5 setup: re-arm group comparison then zoom, save project via the ribbon', async () => {
      await v.waitForCanvasQuiet(page, 'Box plot');
      // Group comparison FIRST: turning it on (and its recompute) resets the
      // value-axis viewport to full — a zoom applied before it would be
      // silently discarded.
      if (await bpProp(page, 'showGroupComparison') !== true) {
        const pt = await revealIcon(page, 'show-group-stats');
        await page.mouse.click(pt.x, pt.y);
        await page.waitForTimeout(1500);
      }
      await setBpProp(page, 'controlComparisons', true, 700);
      await setBpProp(page, 'controlGroup', 'F', 1000);
      await setBpProp(page, 'covariateColumnName', 'HEIGHT', 1200);
      await v.waitForCanvasQuiet(page, 'Box plot');
      // Zoom so a non-full viewport is saved — LAST state-changing action
      // before the save.
      const before = await viewportRect(page);
      const handles = await verticalSliderHandles(page);
      await page.mouse.move(handles.top.x, handles.top.y);
      await page.mouse.down();
      await page.mouse.move(handles.top.x, handles.top.y + (handles.bottom.y - handles.top.y) * 0.30, {steps: 16});
      await page.mouse.up();
      await page.waitForTimeout(1000);
      const zoomed = await viewportRect(page);
      console.log('Scenario 5 viewport height before/after zoom re-arm:', before.height, zoomed.height);
      // Real-margin bar: at least 10% narrower.
      expect(zoomed.height).toBeLessThan(before.height * 0.9);
      await page.evaluate((vp) => (window as any).__preSaveViewport = vp, zoomed);
      const res = await saveProjectViaUI(page, projectName);
      projectIds = {projectId: res.projectId};
      expect(res.projectId.length).toBeGreaterThan(0);
    });

    await softStep('[anchor: PROJECT-ROUND-TRIP] Scenario 5 Step 4-5: reopen the project restores the ladder, zoom, and group-comparison props', async () => {
      const preSave = await page.evaluate(() => (window as any).__preSaveViewport);
      await page.evaluate(() => grok.shell.closeAll());
      await page.waitForTimeout(1500);
      await page.evaluate(async (pid) => {
        const p = await grok.dapi.projects.find(pid);
        await p.open();
      }, projectIds!.projectId);
      await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 30000});
      await page.waitForTimeout(3000);
      const ladder = await readLadder(page);
      console.log('Scenario 5 restored ladder:', JSON.stringify(ladder));
      expect(ladder.valueColumnName).toBe('WEIGHT');
      expect(ladder.category1ColumnName).toBe('SEX');
      expect(ladder.category2ColumnName).toBe('RACE');
      expect(ladder.markerColorColumnName).toBe('SEX');
      expect(ladder.invertColorScheme).toBe(true);
      expect(ladder.colorMin).toBe(20);
      expect(ladder.colorMax).toBe(80);
      expect(ladder.axisType).toBe('logarithmic');
      expect(ladder.invertYAxis).toBe(true);
      expect(ladder.plotStyle).toBe('violin');
      const restoredVp = await viewportRect(page);
      console.log('Scenario 5 viewport preSave/restored:', JSON.stringify(preSave), JSON.stringify(restoredVp));
      expect(restoredVp.height).toBeCloseTo(preSave.height, 0);
      expect(restoredVp.top).toBeCloseTo(preSave.top, 0);
      expect(await bpProp(page, 'showGroupComparison')).toBe(true);
      expect(await bpProp(page, 'controlGroup')).toBe('F');
      expect(await bpProp(page, 'covariateColumnName')).toBe('HEIGHT');
    });
  } finally {
    if (projectIds) await deleteProjectWithCleanup(page, projectIds);
    if (layoutSaved) {
      await page.evaluate(async () => {
        try {
          const layout = (window as any).__probeLayout;
          if (layout?.id) {
            const saved = await grok.dapi.layouts.find(layout.id).catch(() => null);
            if (saved) await grok.dapi.layouts.delete(saved);
          }
        } catch (_) { /* best effort */ }
      }).catch(() => {});
    }
  }

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
