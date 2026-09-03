/* ---
realizes: [boxplot.cp.value-category-axes-persist, boxplot.int.category-sets-marker-color]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {BOX, bpProp, setBpProp, viewportRect, readLadder, clickToggleIcon, dragTopHandle} from './boxplot-helpers';

declare const grok: any;

// Scenarios 1-4 of the settings ladder; the layout round-trip is in memory (tv.saveLayout /
// tv.loadLayout). Scenario 5, whose subject is a project surviving the server, is
// boxplot-settings-ladder-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Box Plot settings ladder and layout round-trip', async ({page}) => {
  test.setTimeout(300_000);

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { if (!isLocalBootNoise(String(e))) pageErrors.push(String(e)); });

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const bp = grok.shell.tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
  });
  await page.locator(BOX).waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);

  await softStep('[anchor: PRE-LADDER] Scenario 1: datetime Value disables the property-panel Axis Type (GROK-20395)', async () => {
    await setBpProp(page, 'valueColumnName', 'STARTED', 1200);

    await page.evaluate(() => {
      grok.shell.o = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    });

    const axisTypeOpacity = await v.pollValue(() => page.evaluate(() => {
      const row = document.querySelector('.property-grid tr[name="prop-axis-type"]') as HTMLElement | null;
      return row ? getComputedStyle(row).opacity : null;
    }), (o) => o !== null, 10000, 250);
    console.log('Scenario 1 Axis Type row opacity (datetime value):', axisTypeOpacity);
    expect(axisTypeOpacity).not.toBeNull();
    expect(parseFloat(axisTypeOpacity as string)).toBeLessThan(1);
    await setBpProp(page, 'valueColumnName', 'AGE', 1000);
    expect(await bpProp(page, 'valueColumnName')).toBe('AGE');
  });

  await softStep('[anchor: Step 4] Scenario 2 Step 3: Category 1 = SEX auto-sets Marker Color to SEX (category-sets-marker-color)', async () => {
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
    await setBpProp(page, 'valueMin', null, 400);
    await setBpProp(page, 'valueMax', null, 700);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'axisType', 'logarithmic', 1400);
    expect(await bpProp(page, 'axisType')).toBe('logarithmic');

    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Step 13 console error delta:', JSON.stringify(errDelta), 'pageerror delta:', JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    const {top, bottom} = await viewportRect(page);
    const weightMax = await page.evaluate(() => grok.shell.t.col('WEIGHT').stats.max);
    console.log('Step 13 viewport top/bottom:', top, bottom, 'WEIGHT max:', weightMax);

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
    await dragTopHandle(page, 0.35);
    const after = await v.pollValue(() => viewportRect(page),
      (vp) => vp.height < before.height * 0.9, 1000, 150);
    console.log('Step 17 viewport height before/after zoom:', before.height, after.height);

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

  await softStep('[anchor: Step 11] Scenario 3: enable Group Comparison, pick control group, set Adjust By covariate', async () => {
    if (await bpProp(page, 'showGroupComparison') !== true)
      await clickToggleIcon(page, 'show-group-stats');
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);

    await setBpProp(page, 'controlComparisons', true, 900);
    await setBpProp(page, 'controlGroup', 'F', 1200);
    expect(await bpProp(page, 'controlGroup')).toBe('F');

    await setBpProp(page, 'covariateColumnName', 'HEIGHT', 1500);
    expect(await bpProp(page, 'covariateColumnName')).toBe('HEIGHT');
    const dom = await page.evaluate((sel) => {
      const root = document.querySelector(sel)!;
      const caption = Array.from(root.querySelectorAll('.d4-column-selector-caption'))
        .some((c) => /Adjust by:/i.test(c.textContent ?? ''));
      const adjustSelector = Array.from(root.querySelectorAll('.d4-column-selector'))
        .map((s) => (s.textContent ?? '').replace(/\s+/g, ' ').trim());
      return {caption, hasHeightSelector: adjustSelector.some((t) => /Adjust by:\s*HEIGHT/i.test(t))};
    }, BOX);
    console.log('Scenario 3 Adjust-by caption/selector:', JSON.stringify(dom));
    expect(dom.caption).toBe(true);
    expect(dom.hasHeightSelector).toBe(true);
  });

  await softStep('[anchor: LAYOUT-ROUND-TRIP] Scenario 4: layout round-trip restores the Box Plot ladder, Scatter re-arm absent', async () => {
    await setBpProp(page, 'covariateColumnName', '', 700);
    await setBpProp(page, 'controlComparisons', false, 500);
    await setBpProp(page, 'showGroupComparison', false, 800);
    const ladderBefore = await readLadder(page);
    await page.evaluate(() => { (window as any).__probeLayout = grok.shell.tv.saveLayout(); });

    await v.waitForViewerRendered(page, 'Box plot', 400);

    await page.evaluate(() => {
      const tv = grok.shell.tv;
      const bp = tv.viewers.find((x: any) => x.type === 'Box plot');
      if (bp) bp.close();
      tv.addViewer('Scatter plot');
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});
    const viewerSet = () => page.evaluate(() => ({
      hasBox: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot'),
      hasScatter: !!grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot'),
    }));
    const armed = await v.pollValue(viewerSet, (s) => !s.hasBox && s.hasScatter, 800, 100);
    expect(armed.hasBox).toBe(false);
    expect(armed.hasScatter).toBe(true);
    await page.evaluate(() => grok.shell.tv.loadLayout((window as any).__probeLayout));
    await page.locator(BOX).waitFor({timeout: 10000});
    const restored = await v.pollValue(viewerSet, (s) => s.hasBox && !s.hasScatter, 3000, 150);
    expect(restored.hasBox).toBe(true);
    expect(restored.hasScatter).toBe(false);

    const ladderAfter = await v.pollValue(() => readLadder(page),
      (l) => l.valueColumnName === 'WEIGHT' && l.axisType === 'logarithmic' && l.plotStyle === 'violin',
      3000, 150);
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

  await v.closeAllAndWait(page);
  v.finishSpec();
});
