/* ---
realizes: [trellisplot.cp.global-scale-inner-axes, trellisplot.int.global-scale-range-slider-sync]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Ambient noise the packed dev build emits regardless of what is actuated. Nothing
// product-specific is filtered, so a NullError or a Dart stack trace still counts.
const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

// Per-cell canvas hash. Returns null for a missing or unreadable canvas, so EVERY
// caller must null-guard: a vanished canvas hashes null on both sides of a diff and
// would make the assertion pass vacuously.
async function cellHashes(page: Page, idxs: number[]): Promise<(number | null)[]> {
  return page.evaluate((idxs) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    function hash(cellIdx: number): number | null {
      const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
      const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
      if (!cv) return null;
      try {
        const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        let h = 0;
        for (let i = 0; i < img.length; i += 4)
          h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
        return h;
      } catch { return null; }
    }
    return idxs.map(hash);
  }, idxs);
}

// Context menu on the charts-grid background, NOT inside a cell. Menu-open works
// synthetically — unlike the Scenario 2 slider drag, which no-ops without trusted input.
async function trellisMenuLabels(page: Page): Promise<string[]> {
  const labels = await page.evaluate(async () => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const grid = (root.querySelector('.d4-trellis-plot-charts-grid') as HTMLElement) ?? root;
    const gr = grid.getBoundingClientRect();
    grid.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true,
      clientX: gr.left + 4, clientY: gr.top + 4}));
    await new Promise((r) => setTimeout(r, 800));
    // Scoped to `.d4-menu-popup` ON PURPOSE: an unscoped `.d4-menu-item-label` sweep also
    // returns the APPLICATION MAIN menu, so the label is graded against whichever surface
    // happens to own the string (trellis_plot.md pitfall 22 [DOM 2026-08-11]).
    const out = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .map((e) => (e as HTMLElement).innerText.trim());
    // close the menu (outside mousedown + Escape) so it does not swallow later gestures.
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
    return out;
  });
  await page.waitForTimeout(400);
  return labels;
}

// [type="range-slider"] matches BOTH the inner-axis sliders and the bare category
// scroll sliders, so an axis toggle can only be asserted as a count DELTA.
async function rangeSliderCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    return root.querySelectorAll('[type="range-slider"]').length;
  });
}

// The wheel-target point: DPR is 1 on dev, and the cell div and its canvas share a box.
async function cellCenter(page: Page, idx: number): Promise<{x: number; y: number} | null> {
  return page.evaluate((idx) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const cell = root.querySelectorAll('.d4-trellis-plot-cell')[idx] as HTMLElement | undefined;
    if (!cell) return null;
    const r = cell.getBoundingClientRect();
    return {x: r.left + r.width / 2, y: r.top + r.height / 2};
  }, idx);
}

// The type switch is a SETTING, not the tested action, so driving it through the
// prop is legitimate; the wheel below is the tested action and must stay trusted.
async function setInnerType(page: Page, viewerType: string): Promise<void> {
  await page.evaluate(async (vt) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.viewerType = vt;
    await new Promise((r) => setTimeout(r, 2500));
  }, viewerType);
}

// Prove an inner-type switch LANDED (beforeSwitch = cell-0 hash from just before it). The cell holds a BARE
// CANVAS, not a named viewer root (viewer_base.dart:158; bar_chart_meta.dart:26, box_plot_meta.dart:18,
// scatterplot_meta.dart:28), so identity is a weak viewerType read-back plus the repaint delta. [DOM 2026-08-11]
async function assertInnerTypeSwitched(page: Page, viewerType: string,
  beforeSwitch: number | null): Promise<void> {
  const applied = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.props.viewerType as string;
  });
  const [afterSwitch] = await cellHashes(page, [0]);
  console.log(`[inner type -> ${viewerType}] readBack=${applied} hash ${beforeSwitch} -> ${afterSwitch}`);
  expect(applied, `inner viewer type did not switch to ${viewerType}`).toBe(viewerType);
  expect(beforeSwitch).not.toBeNull();
  expect(afterSwitch).not.toBeNull();
  expect(afterSwitch,
    `cells were not repainted by the switch to ${viewerType} — the arm below would probe the previous inner type`)
    .not.toBe(beforeSwitch);
}

const innerTypeTabs = ['Scatter plot', 'Bar chart', 'Box plot', 'Histogram', 'Line chart', 'Pie chart'];

// PURE GLUE: opens the panel, selects the tab, expands the collapsed "Misc" category.
// It NEVER writes a property value — that is what lets Step 2 assert the untouched
// default.
async function openInnerViewerTab(page: Page, tabName?: string): Promise<void> {
  // The gear lives on the dock `.panel-titlebar`, a SIBLING of the viewer root.
  await v.openViewerGear(page, 'Trellis plot');
  // Fall back to any viewer-type-labelled tab, so the cleanup path still lands on the
  // inner-viewer tab.
  await page.evaluate(({name, names}) => {
    const headers = Array.from(document.querySelectorAll('.d4-tab-header')) as HTMLElement[];
    const tab = (name ? headers.find((h) => h.innerText.trim() === name) : undefined) ??
      headers.find((h) => names.includes(h.innerText.trim()));
    tab?.click();
  }, {name: tabName ?? '', names: innerTypeTabs});
  await page.waitForTimeout(800);
  // The Allow Zoom row sits under the COLLAPSED "Misc" category — in the DOM at height
  // 0, so a visibility-gated click times out. The expand is guarded: clicking an already
  // expanded category would collapse it again. Expanding is layout, not a value write.
  const row = page.locator('.property-grid tr[name="prop-allow-zoom"]');
  if (await row.count() > 0 && !(await row.isVisible())) {
    await page.locator('.property-grid tr[name="prop-category-misc"]').first().click();
    await page.waitForTimeout(500);
  }
}

// READ ONLY — no click, no write. Returns null when the row is absent: `allowZoom` is
// declared on exactly two looks (scatterplot_look.dart:294, density_plot_look.dart), so a
// Bar chart or Box plot inner viewer has no such row. Call openInnerViewerTab first.
async function allowZoomState(page: Page): Promise<boolean | null> {
  return page.evaluate(() => {
    const row = document.querySelector('.property-grid tr[name="prop-allow-zoom"]');
    const cb = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
    return cb ? cb.checked : null;
  });
}

// Drives the checkbox through the real property-grid editor. Only for the managed Step 5/6
// transitions. Never call it on a path that then asserts the product default: the write
// would repair the very regression the step exists to expose.
async function setAllowZoom(page: Page, desired: boolean, tabName = 'Scatter plot'): Promise<boolean | null> {
  await openInnerViewerTab(page, tabName);
  const row = page.locator('.property-grid tr[name="prop-allow-zoom"]');
  if (await row.count() === 0) return null;
  await row.waitFor({state: 'visible', timeout: 10000});
  const box = row.locator('input[type="checkbox"]').first();
  if (await box.isChecked() !== desired) {
    await box.click();
    await page.waitForTimeout(1500);
  }
  return box.isChecked();
}

// Trusted wheel over a page point: real CDP input the inner-viewer canvas
// tracks (synthetic wheel events do not drive the d4 mouse-wheel handler).
async function wheelOver(page: Page, pt: {x: number; y: number}, steps = 5): Promise<void> {
  await page.mouse.move(pt.x, pt.y);
  await page.waitForTimeout(200);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, 120);
    await page.waitForTimeout(60);
  }
  await page.waitForTimeout(600);
}

test('Trellis plot: global scale, inner axes, range slider reset', async ({page}) => {
  test.setTimeout(900_000);
  page.setDefaultTimeout(120_000);

  // Attached before login, so it covers every actuation below. The client exposes no
  // `grok.shell.warnings` symbol, so a floor built on it could never fail. Every step asserts a
  // DELTA over its own actuation — a global zero floor would be false-red on ambient noise.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);
  await page.waitForTimeout(5000);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {rowCount: df.rowCount, sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
  });
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await page.evaluate(async () => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
    tp.props.showRangeSliders = true;
    tp.props.showXAxes = 'Always';
    tp.props.showYAxes = 'Always';
    await new Promise((r) => setTimeout(r, 2000));
  });
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
  await expect(cellLocator).toHaveCount(8);

  // ================= Scenario 1: Global Scale toggle and axis-driven slider presence =================

  // #### Scenario 1 Step 5: enabling Global Scale re-renders the probed cells
  await softStep('Scenario 1 Step 5', async () => {
    // Probe cells row0col0 and row1col0 — the flat index is yi*xCat+xi, hence 0 and 2.
    const probes = [0, 2];
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    const before = await cellHashes(page, probes);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.globalScale = true;
      await new Promise((r) => setTimeout(r, 2500));
    });
    const after = await cellHashes(page, probes);
    // Two cells still rendering different pixels is what proves the per-cell binding
    // survived the switch to a shared scale.
    expect(before.every((h) => h !== null)).toBe(true);
    expect(after.every((h) => h !== null)).toBe(true);
    console.log(`[Scenario 1 Step 5] before=${JSON.stringify(before)} after=${JSON.stringify(after)}`);
    expect(after[0]).not.toBe(before[0]);
    expect(after[1]).not.toBe(before[1]);
    expect(after[0]).not.toBe(after[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
  });

  // #### Scenario 1 Step 7: the context menu offers 'Reset Inner Range Sliders'
  await softStep('Scenario 1 Step 7', async () => {
    const labels = await trellisMenuLabels(page);
    expect(labels).toContain('Reset Inner Range Sliders');
  });

  // #### Scenario 1 Step 9: hiding BOTH inner axes removes the reset item
  await softStep('Scenario 1 Step 9', async () => {
    // [DOM 2026-08-06]: hiding only the X axis does NOT remove the item — it exists
    // while ANY inner-axis slider does (inner_viewer_axes.dart:243-245), so BOTH axes
    // must be hidden.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Never';
      tp.props.showYAxes = 'Never';
      await new Promise((r) => setTimeout(r, 2000));
    });
    const labels = await trellisMenuLabels(page);
    // Gate the negative before making it: a menu that never opened yields [], and []
    // satisfies not.toContain. 'Properties...' sits on the trellis menu regardless of
    // axis state (trellis_plot.md, Context Menu > Top-level items), so it witnesses a real open menu.
    expect(labels).toContain('Properties...');
    expect(labels).not.toContain('Reset Inner Range Sliders');
  });

  // #### Scenario 1 Step 10: re-enabling both axes brings the reset item back
  await softStep('Scenario 1 Step 10', async () => {
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      await new Promise((r) => setTimeout(r, 2000));
    });
    const labels = await trellisMenuLabels(page);
    expect(labels).toContain('Reset Inner Range Sliders');
  });

  // #### Scenario 1 Step 11: the [type="range-slider"] count drops when an axis is hidden
  await softStep('Scenario 1 Step 11', async () => {
    const countAxesShown = await rangeSliderCount(page);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Never';
      await new Promise((r) => setTimeout(r, 2000));
    });
    const countXHidden = await rangeSliderCount(page);
    console.log(`[Scenario 1 Step 11] count_axes_shown=${countAxesShown} count_x_hidden=${countXHidden}`);
    expect(countAxesShown).toBeGreaterThan(countXHidden);
    // Restore both axes shown for Scenario 2.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      await new Promise((r) => setTimeout(r, 2000));
    });
    await expect(cellLocator).toHaveCount(8);
  });

  // ================= Scenario 2: Shared range slider re-bounds all cells; Reset restores =================

  // #### Scenario 2 Step 4: a trusted drag on the shared X-axis slider re-bounds every cell
  await softStep('Scenario 2 Step 4', async () => {
    // With Global Scale on, the single shared X-axis slider re-bounds ALL cells at once
    // (features/inner_viewer_axes.dart:202-215). The drag needs TRUSTED input —
    // synthetic Mouse/Pointer events silently no-op [DOM 2026-08-06].
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.globalScale = true;
      await new Promise((r) => setTimeout(r, 1500));
    });
    const probes = [0, 2];
    // Window opens AFTER the carried-over globalScale write, so the delta grades the
    // hover-reveal, the drag and its settle — the actuation this step owns.
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    const baseline = await cellHashes(page, probes);
    expect(baseline.every((h) => h !== null)).toBe(true);
    expect(baseline[0]).not.toBe(baseline[1]);

    // The inner sliders are visibility:hidden until the axis area is hovered, hence the
    // mousemove before the track geometry is read.
    const track = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const sliders = Array.from(root.querySelectorAll(
        '.d4-range-selector > svg[type="range-slider"][name="x-slider"]')) as SVGElement[];
      if (sliders.length === 0) return null;
      const s = sliders[0];
      const r = s.getBoundingClientRect();
      s.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: r.left + r.width / 2, clientY: r.top + 5}));
      return {left: r.left, top: r.top, width: r.width, height: r.height};
    });
    expect(track).not.toBeNull();
    const t = track!;
    const y = t.top + t.height / 2;
    await page.mouse.move(t.left + t.width - 6, y);
    await page.mouse.down();
    await page.mouse.move(t.left + t.width * 0.75, y, {steps: 6});
    await page.mouse.move(t.left + t.width * 0.5, y, {steps: 8});
    await page.mouse.up();
    await page.waitForTimeout(2000);

    const after = await cellHashes(page, probes);
    console.log(`[Scenario 2 Step 4] baseline=${JSON.stringify(baseline)} after=${JSON.stringify(after)}`);
    expect(after.every((h) => h !== null)).toBe(true);
    // Both probed cells moved off baseline: the re-bound reached every cell, not only
    // the one under the slider.
    expect(after[0]).not.toBe(baseline[0]);
    expect(after[1]).not.toBe(baseline[1]);
    expect(after[0]).not.toBe(after[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    // Stash BOTH states for Step 5: the narrowed one it must move off, and the baseline
    // it must land back on exactly.
    await page.evaluate(({narrowed, base}) => {
      (window as any).__narrowed = narrowed;
      (window as any).__baseline = base;
    }, {narrowed: after, base: baseline});
  });

  // #### Scenario 2 Step 5: 'Reset Inner Range Sliders' restores every cell
  await softStep('Scenario 2 Step 5', async () => {
    // Reset is graded as an EXACT return to the full-range baseline, the live-verified behaviour
    // (trellis_plot.md [DOM 2026-08-05]). "Moved off the narrowed state" alone would be satisfied
    // by any repaint, so it is kept only as the companion anti-vacuity check.
    const probes = [0, 2];
    const narrowed = await page.evaluate(() => (window as any).__narrowed as (number | null)[]);
    const baseline = await page.evaluate(() => (window as any).__baseline as (number | null)[]);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const grid = (root.querySelector('.d4-trellis-plot-charts-grid') as HTMLElement) ?? root;
      const gr = grid.getBoundingClientRect();
      grid.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true,
        clientX: gr.left + 4, clientY: gr.top + 4}));
      await new Promise((r) => setTimeout(r, 800));
      // Same `.d4-menu-popup` scope as trellisMenuLabels, for the same reason.
      const target = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .find((e) => (e as HTMLElement).innerText.trim() === 'Reset Inner Range Sliders');
      (target?.closest('.d4-menu-item') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 2000));
    });
    const after = await cellHashes(page, probes);
    console.log(`[Scenario 2 Step 5] baseline=${JSON.stringify(baseline)} narrowed=${JSON.stringify(narrowed)} ` +
      `after=${JSON.stringify(after)}`);
    expect(after.every((h) => h !== null)).toBe(true);
    expect(after[0]).not.toBe(narrowed[0]);
    expect(after[1]).not.toBe(narrowed[1]);
    expect(after[0], 'cell 0 did not return exactly to its pre-drag baseline').toBe(baseline[0]);
    expect(after[1], 'cell 2 did not return exactly to its pre-drag baseline').toBe(baseline[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
  });

  // ================= Scenario 3: inner viewers do not zoom on wheel (GROK-14587) =================

  // Only Scatter has a gate: `allowZoom` is declared on the scatter look
  // (scatterplot_look.dart:294, default true) and turned off by the trellis() preset
  // (:498,:514). Bar chart and Box plot have no such property and do not zoom at all.

  // Scenario 3 owns its start state: the Scenario 1-2 ladder leaves Show Range Sliders
  // on, and those sliders are hover-revealed — every wheel gesture below hovers the cell
  // first, so the reveal would move the very pixels the per-cell diff is watching.
  await page.evaluate(async () => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.globalScale = false;
    tp.props.showRangeSliders = false;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
    await new Promise((r) => setTimeout(r, 2500));
  });
  await expect(cellLocator).toHaveCount(8);

  try {
    // #### Scenario 3 Step 2: scatter plot cell does not zoom on wheel (Allow Zoom off by default)
    await softStep('Scenario 3 Step 2', async () => {
      await setInnerType(page, 'Scatter plot');
      // GROK-14587 invariant, read UNTOUCHED: the trellis() preset that turns allowZoom off IS
      // the regression surface, so it must be observed before anything writes it — a setter
      // would silently repair a regressed preset instead of failing here.
      await openInnerViewerTab(page, 'Scatter plot');
      const defaultAllowZoom = await allowZoomState(page);
      console.log(`[Scenario 3 Step 2] untouched default allowZoom=${defaultAllowZoom}`);
      expect(defaultAllowZoom).toBe(false);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 2] scatter before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    // #### Scenario 3 Step 3: bar chart cell does not zoom on wheel
    await softStep('Scenario 3 Step 3', async () => {
      // No Allow Zoom read here by design: a Bar chart inner viewer has no such property,
      // so this arm asserts the absence of zoom itself, not the state of a gate.
      const [beforeSwitch] = await cellHashes(page, [0]);
      await setInnerType(page, 'Bar chart');
      // Structural settle only: 8 cells is the category product for EVERY inner type, so
      // it says nothing about which type the cells hold — that is the witness below.
      await expect(cellLocator).toHaveCount(8);
      await assertInnerTypeSwitched(page, 'Bar chart', beforeSwitch);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 3] bar before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    // #### Scenario 3 Step 4: box plot cell does not zoom on wheel
    await softStep('Scenario 3 Step 4', async () => {
      // Same as Step 3: Box plot has no `allowZoom` property either.
      const [beforeSwitch] = await cellHashes(page, [0]);
      await setInnerType(page, 'Box plot');
      await expect(cellLocator).toHaveCount(8);
      await assertInnerTypeSwitched(page, 'Box plot', beforeSwitch);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 4] box before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    // #### Scenario 3 Step 5: with Allow Zoom = true, wheel zooms the scatter cell (positive control)
    await softStep('Scenario 3 Step 5', async () => {
      // Scatter is the one inner type with the gate, so this managed false -> true transition is
      // what setAllowZoom is for. No separate switch witness is needed: `prop-allow-zoom` exists
      // for Scatter ONLY, so setAllowZoom returns null if the cells are still on Box plot.
      await setInnerType(page, 'Scatter plot');
      const on = await setAllowZoom(page, true);
      expect(on).toBe(true);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 5] allowZoom=true before=${before} after=${after}`);
      expect(after).not.toBeNull();
      // Positive control: the same gesture that left the cell unchanged in Steps 2-4 now
      // moves it, so those "unchanged" results cannot be a wheel that never arrived.
      expect(after).not.toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    // #### Scenario 3 Step 6: resetting Allow Zoom to false restores the no-zoom default
    await softStep('Scenario 3 Step 6', async () => {
      // The other half of the managed round-trip (true -> false): it proves that the
      // GATE, not the wheel plumbing, is what changed between Step 5 and here.
      const off = await setAllowZoom(page, false);
      expect(off).toBe(false);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 6] allowZoom=false before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });
  } finally {
    // Allow Zoom is a persisted inner-viewer property: never leave the viewer opted in,
    // even when a step above threw before the round-trip reached Step 6.
    await setAllowZoom(page, false).catch(() => {});
  }

  v.finishSpec();
});
