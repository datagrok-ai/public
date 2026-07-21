import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

test('PC Plot tests', async ({page}) => {
  test.setTimeout(600_000);

  // The canvas-only steps below change nothing but how the plot is painted, so the
  // check is that driving them raises nothing and leaves the viewer alive.
  // grok.shell.warnings is undefined on this build, hence the page/console baseline.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const viewerAlive = () => page.evaluate(() =>
    !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot')
    && !!document.querySelector('[name="viewer-PC-Plot"]'));

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 10000});

  // Both adds go through the Toolbox icon ([name="icon-pc-plot"]); the Menu
  // Ribbon "Add viewer" gallery is a canvas-rendered dialog with no headless
  // handles (see pc-plot.md Actuation note). To Script emits a code balloon.
  await softStep('Menu Ribbon and To Script', async () => {
    const viewerExists = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(viewerExists).toBe(true);

    // Warm up the context menu — the first menu on a freshly attached viewer
    // builds cold and the To Script action does not fire on it.
    await page.evaluate(async () => {
      const canvas = document.querySelector('[name="viewer-PC-Plot"] canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2}));
      await new Promise(r => setTimeout(r, 500));
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise(r => setTimeout(r, 800));
    });

    const balloon = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const reset = async () => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        await new Promise(r => setTimeout(r, 500));
      };
      const attempt = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 700));
        const toScript = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(el => el.textContent!.trim() === 'To Script');
        if (!toScript) { await reset(); return ''; }
        const parent = toScript.closest('.d4-menu-item')!;
        parent.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        parent.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 500));
        const js = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(el => el.textContent!.trim() === 'To JavaScript');
        if (!js) { await reset(); return ''; }
        js.closest('.d4-menu-item')!.click();
        // The generated script is delivered in a `.d4-balloon`; poll for it.
        for (let i = 0; i < 24; i++) {
          const b = document.querySelector('.d4-balloon');
          if (b && ((b as HTMLElement).innerText || '').length > 0) return (b as HTMLElement).innerText;
          await new Promise(r => setTimeout(r, 250));
        }
        await reset();
        return '';
      };
      // The first context menu on a fresh viewer can build cold — retry, always
      // resetting the menu state between tries so a stuck-open menu can't jam
      // the next attempt.
      let text = '';
      for (let a = 0; a < 5 && !text; a++) text = await attempt();
      return {present: text.length > 0, text};
    });
    // A non-empty balloon carrying the generated viewer-creation script.
    expect(balloon.present).toBe(true);
    expect(balloon.text).toContain('addViewer');

    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot');
      if (pc) pc.close();
      await new Promise(r => setTimeout(r, 500));
      const icon = document.querySelector('[name="icon-pc-plot"]');
      if (icon) (icon as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
    });
    const reopened = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(reopened).toBe(true);
  });

  await softStep('Axis scale via the context menu', async () => {
    const menuResult = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
      };
      const clickSub = async (parent: string, child: string) => {
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return false;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 300));
        return !!c;
      };
      await openMenu();
      await clickSub('Y Axis', 'Global');
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const afterGlobal = pc.props.normalizeEachColumn;
      await openMenu();
      await clickSub('Y Axis', 'Normalized');
      const afterNorm = pc.props.normalizeEachColumn;
      return { afterGlobal, afterNorm };
    });
    expect(menuResult.afterGlobal).toBe(false);
    expect(menuResult.afterNorm).toBe(true);
  });

  // The current/mouse-over/all painted lines are a canvas outcome with no DOM
  // counterpart, so the prop-drive block is a no-error floor. The Selection
  // context-menu items DO flip readable props, so the menu toggles for Show
  // Current Line / Show All Lines are asserted as a menu -> state round-trip.
  await softStep('Selection & line display', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.showCurrentLine = false; await wait();
      pc.props.showCurrentLine = true; await wait();
      pc.props.showMouseOverLine = false; await wait();
      pc.props.showMouseOverLine = true; await wait();
      pc.props.showMouseOverRowGroup = true; await wait();
      pc.props.showAllLines = false; await wait();
      pc.props.showAllLines = true; await wait();
      pc.props.showMouseOverRowGroup = false;
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);

    // Context-menu Selection toggles flip the same props the Context Panel does.
    const menu = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const clickSub = async (parent: string, child: string) => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
        const p = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === parent);
        if (!p) return false;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 350));
        const c = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === child);
        if (c) { c.closest('.d4-menu-item')!.click(); await new Promise(r => setTimeout(r, 350)); }
        return !!c;
      };
      const curBefore = pc.props.showCurrentLine;
      await clickSub('Selection', 'Show Current Line');
      const curToggled = pc.props.showCurrentLine;
      await clickSub('Selection', 'Show Current Line');
      const curRestored = pc.props.showCurrentLine;
      const allBefore = pc.props.showAllLines;
      await clickSub('Selection', 'Show All Lines');
      const allToggled = pc.props.showAllLines;
      await clickSub('Selection', 'Show All Lines');
      const allRestored = pc.props.showAllLines;
      return {curBefore, curToggled, curRestored, allBefore, allToggled, allRestored};
    });
    expect(menu.curToggled).toBe(!menu.curBefore);
    expect(menu.curRestored).toBe(menu.curBefore);
    expect(menu.allToggled).toBe(!menu.allBefore);
    expect(menu.allRestored).toBe(menu.allBefore);
    expect(errorCount()).toBe(errBefore);
  });

  // Line widths, label orientation and margins are pure painting (the Lines
  // context-menu Line Width slider writes the same `lineWidth` prop driven
  // here — see pc-plot.md). The axis sliders are read afterwards to confirm the
  // layout pass rebuilt the plot.
  await softStep('Style & layout', async () => {
    const errBefore = errorCount();
    const sliders = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.lineWidth = 3; await wait();
      pc.props.currentLineWidth = 5; await wait();
      pc.props.mouseOverLineWidth = 5; await wait();
      pc.props.labelsOrientation = 'Vert'; await wait();
      pc.props.minMaxOrientation = 'Vert'; await wait();
      pc.props.horzMargin = 60; await wait();
      pc.props.autoLayout = false; await wait();
      pc.props.lineWidth = 0.5; pc.props.currentLineWidth = 2; pc.props.mouseOverLineWidth = 2;
      pc.props.labelsOrientation = 'Auto'; pc.props.minMaxOrientation = 'Auto';
      pc.props.horzMargin = 40; pc.props.autoLayout = true;
      await new Promise(r => setTimeout(r, 400));
      return document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]').length;
    });
    expect(sliders).toBeGreaterThan(0);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Narrow two per-axis range sliders (real DOM handles), watch the shared
  // df.filter drop, restore it with Reset View, then round-trip the Show
  // Filters state through the context menu (the slider DOM persists; the
  // toggled state lives on the `showFilters` prop).
  await softStep('Reset and filter visibility from the context menu', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const df = grok.shell.tv.dataFrame;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const vr = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
      await new Promise(r => setTimeout(r, 400));
      const fullCount = df.filter.trueCount;

      const dragMax = async (axis: string) => {
        const svg = document.querySelector(`[name="axis-slider-${axis}"]`);
        if (!svg) return false;
        const maxHandle = svg.querySelector('[name="max-handle"]')!;
        const hr = maxHandle.getBoundingClientRect();
        const cx = hr.x + hr.width / 2, cy = hr.y + hr.height / 2;
        const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
        maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await new Promise(r => setTimeout(r, 50));
        for (let dy = 20; dy <= 200; dy += 30) {
          document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          await new Promise(r => setTimeout(r, 20));
        }
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 200)));
        await new Promise(r => setTimeout(r, 500));
        return true;
      };
      await dragMax('AGE');
      const afterAge = df.filter.trueCount;
      await dragMax('HEIGHT');
      const afterBoth = df.filter.trueCount;

      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();
      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
        await new Promise(r => setTimeout(r, 500));
      };
      await openMenu();
      const rv = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === 'Reset View');
      if (rv) rv.closest('.d4-menu-item')!.click();
      await new Promise(r => setTimeout(r, 700));
      const afterReset = df.filter.trueCount;

      const clickFilterSub = async (child: string) => {
        await openMenu();
        const p = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === 'Filter');
        if (!p) return false;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 350));
        const c = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === child);
        if (c) { c.closest('.d4-menu-item')!.click(); await new Promise(r => setTimeout(r, 400)); }
        return !!c;
      };
      const showBefore = pc.props.showFilters;
      await clickFilterSub('Show Filters');
      const showToggled = pc.props.showFilters;
      await clickFilterSub('Show Filters');
      const showRestored = pc.props.showFilters;
      return {fullCount, afterAge, afterBoth, afterReset, showBefore, showToggled, showRestored};
    });
    // Two-axis narrowing drops the filter progressively, Reset View fully restores it.
    expect(result.afterAge).toBeLessThan(result.fullCount);
    expect(result.afterBoth).toBeLessThan(result.afterAge);
    expect(result.afterReset).toBe(result.fullCount);
    // Show Filters menu item round-trips the showFilters state.
    expect(result.showToggled).toBe(!result.showBefore);
    expect(result.showRestored).toBe(result.showBefore);
  });

  // Both the Filter Panel and the in-chart range sliders write the shared
  // df.filter (AND-combined). Reset View clears ONLY the in-chart part; the
  // Filter Panel "Reset filters" button clears everything back to full.
  await softStep('Filter panel interaction', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const df = grok.shell.tv.dataFrame;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const vr = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
      await new Promise(r => setTimeout(r, 400));
      const fullCount = df.filter.trueCount;

      // Filter Panel histogram narrows AGE.
      grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise(r => setTimeout(r, 700));
      const afterPanel = df.filter.trueCount;

      // In-chart range slider narrows HEIGHT on top of the Filter Panel filter.
      const dragHeight = async () => {
        const svg = document.querySelector('[name="axis-slider-HEIGHT"]');
        if (!svg) return false;
        const maxHandle = svg.querySelector('[name="max-handle"]')!;
        const hr = maxHandle.getBoundingClientRect();
        const cx = hr.x + hr.width / 2, cy = hr.y + hr.height / 2;
        const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
        maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await new Promise(r => setTimeout(r, 50));
        for (let dy = 20; dy <= 200; dy += 30) {
          document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          await new Promise(r => setTimeout(r, 20));
        }
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 200)));
        await new Promise(r => setTimeout(r, 500));
        return true;
      };
      await dragHeight();
      const afterBoth = df.filter.trueCount;

      // Reset View clears only the in-chart slider; the Filter Panel filter survives.
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
      await new Promise(r => setTimeout(r, 500));
      const rv = Array.from(document.querySelectorAll('.d4-menu-item-label')).find(el => el.textContent!.trim() === 'Reset View');
      if (rv) rv.closest('.d4-menu-item')!.click();
      await new Promise(r => setTimeout(r, 700));
      const afterResetView = df.filter.trueCount;

      // Re-narrow, then the Filter Panel "Reset filters" button clears everything.
      await dragHeight();
      const afterReDrag = df.filter.trueCount;
      const btn = document.querySelector('.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement | null;
      if (btn) btn.click();
      await new Promise(r => setTimeout(r, 800));
      const afterPanelReset = df.filter.trueCount;
      return {fullCount, afterPanel, afterBoth, afterResetView, afterReDrag, afterPanelReset};
    });
    // Filter Panel filter takes effect, and the in-chart slider narrows further.
    expect(result.afterPanel).toBeLessThan(result.fullCount);
    expect(result.afterBoth).toBeLessThan(result.afterPanel);
    // Reset View resets ONLY the in-chart part: the Filter Panel filter remains.
    expect(result.afterResetView).toBe(result.afterPanel);
    // Re-narrowing works again, and the Filter Panel reset restores the full count.
    expect(result.afterReDrag).toBeLessThan(result.afterPanel);
    expect(result.afterPanelReset).toBe(result.fullCount);
  });

  // Each axis carries its own range slider named after its column, and the sliders
  // sit in painted order, so their sequence is the rendered axis order.
  await softStep('Column reordering from the Context Panel list', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const sliderOrder = () =>
        Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
          .map((e) => e.getAttribute('name')!.replace('axis-slider-', ''));
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const before = sliderOrder();
      pc.props.columnNames = ['WEIGHT', 'AGE', 'HEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const after = sliderOrder();
      return { before, after };
    });
    expect(result.before).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(result.after).toEqual(['WEIGHT', 'AGE', 'HEIGHT']);
  });

  // densityStyle defaults to 'circles' (and showDensity defaults on). Every
  // box-plot component toggle draws to canvas, so the rest is the no-error floor.
  await softStep('Density component toggles', async () => {
    const errBefore = errorCount();
    const defaults = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      return {densityStyle: pc.props.densityStyle, showDensity: pc.props.showDensity};
    });
    expect(defaults.densityStyle).toBe('circles');
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.densityStyle = 'box plot'; await wait();
      pc.props.showInterquartileRange = false; await wait();
      pc.props.showInterquartileRange = true;
      pc.props.showUpperDash = false; pc.props.showUpperDash = true; await wait();
      pc.props.showLowerDash = false; pc.props.showLowerDash = true; await wait();
      pc.props.showMeanCross = false; pc.props.showMeanCross = true; await wait();
      pc.props.showMedian = false; pc.props.showMedian = true; await wait();
      pc.props.showCircles = true; await wait();
      pc.props.densityStyle = 'violin plot'; await wait();
      pc.props.bins = 200; await wait();
      pc.props.whiskerLineWidth = 5; await wait();
      pc.props.densityStyle = 'circles';
      pc.props.bins = 100; pc.props.whiskerLineWidth = 2;
      pc.props.showDensity = false;
      await new Promise(r => setTimeout(r, 300));
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Legend visibility is a real DOM signal; the gradient options (log axis,
  // inversion, min/max clamps) are canvas-only and are driven under the floor.
  // Conditional grid color coding renders a DOM legend listing its bins, while
  // linear/numeric coloring has no DOM legend — that contrast is the readable
  // signal that the plot picked up the grid column's color-coding change.
  await softStep('Color coding, legend & grid coloring', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = (ms = 600) => new Promise(r => setTimeout(r, ms));
      const legend = () => {
        const el = document.querySelector('[name="viewer-PC-Plot"] .d4-legend') as HTMLElement | null;
        return {present: !!el, labels: el ? (el.innerText || '').split('\n').map(s => s.trim()).filter(Boolean) : [],
          text: el ? (el.innerText || '').replace(/\s+/g, ' ').trim() : ''};
      };

      // Numeric colouring: gradient options are canvas-only, just drive them.
      pc.props.colorColumnName = 'AGE'; await wait();
      pc.props.colorAxisType = 'logarithmic'; await wait(300);
      pc.props.invertColorScheme = true; await wait(300);
      pc.props.invertColorScheme = false;
      pc.props.colorMin = 30; pc.props.colorMax = 60; await wait(300);
      pc.props.colorMin = null; pc.props.colorMax = null; pc.props.colorAxisType = 'linear';

      // Categorical colouring: the legend lists the column's categories.
      pc.props.colorColumnName = 'RACE'; await wait();
      const categorical = legend();
      pc.props.legendPosition = 'Left'; await wait(300);
      pc.props.legendPosition = 'Right'; pc.props.legendPosition = 'Top';
      pc.props.legendPosition = 'Bottom'; await wait(300);
      pc.props.legendVisibility = 'Never'; await wait();
      const hidden = legend();
      pc.props.legendVisibility = 'Auto'; await wait();
      const restored = legend();

      // Colour coding set on the grid column, read back through the plot legend.
      pc.props.colorColumnName = 'HEIGHT';
      pc.props.legendPosition = 'Auto'; pc.props.legendVisibility = 'Auto';
      const df = grok.shell.tv.dataFrame;
      df.col('HEIGHT').meta.colors.setConditional({'20-150': DG.Color.green, '150-250': DG.Color.orange});
      await wait(800);
      const conditional = legend();
      df.col('HEIGHT').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
      await wait(800);
      const linear = legend();
      df.col('HEIGHT').meta.colors.setLinear();
      pc.props.colorColumnName = '';
      pc.props.legendPosition = 'Auto';
      await wait(300);
      return {categorical, hidden, restored, conditional, linear};
    });
    expect(result.categorical.present).toBe(true);
    expect(result.hidden.present).toBe(false);
    expect(result.restored.labels).toEqual(result.categorical.labels);
    // Conditional coding on the grid column surfaces its bins in the plot legend.
    expect(result.conditional.present).toBe(true);
    expect(result.conditional.text).toContain('20-150');
    expect(result.conditional.text).toContain('150-250');
    // Switching to a linear/numeric scheme drops the DOM legend (gradient is canvas).
    expect(result.linear.present).toBe(false);
    expect(errorCount()).toBe(errBefore);
  });

  // The description is rendered inside the viewer element and can be read back;
  // the title renders in the panel titlebar (.panel-titlebar-text).
  await softStep('Title and description', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 600));
      const root = pc.root as HTMLElement;
      const panel = (root.closest('.panel-base') as HTMLElement) ?? root;
      const titlebarText = () => (panel.querySelector('.panel-titlebar-text')?.textContent ?? '').trim();
      const shownText = () => (root.innerText || '').replace(/\s+/g, ' ').trim();
      pc.props.title = 'My PC Plot'; await wait();
      const titleShown = titlebarText();
      pc.props.description = 'Test description'; await wait();
      const withDescription = shownText();
      pc.props.descriptionPosition = 'Bottom'; await wait();
      const moved = shownText();
      pc.props.title = ''; pc.props.description = ''; await wait();
      const titleCleared = titlebarText();
      const cleared = shownText();
      return {titleShown, withDescription, moved, titleCleared, cleared};
    });
    expect(result.titleShown).toContain('My PC Plot');
    expect(result.withDescription).toContain('Test description');
    expect(result.moved).toContain('Test description');
    expect(result.titleCleared).not.toContain('My PC Plot');
    expect(result.cleared).not.toContain('Test description');
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Pick Up / Apply', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pc1 = tv.viewers.find(v => v.type === 'PC Plot')!;
      pc1.props.columnNames = ['AGE', 'WEIGHT', 'STARTED'];
      pc1.props.logColumnsColumnNames = ['AGE'];
      pc1.props.colorColumnName = 'RACE';
      pc1.props.legendPosition = 'Left';
      pc1.props.title = 'Source Plot';
      tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      const pcs = () => tv.viewers.filter(v => v.type === 'PC Plot');
      // Address the target plot through tv.viewers (same ordering the reads
      // below use) so a DOM-vs-viewers order mismatch can't pass vacuously.
      const clickSub = async (idx: number, parent: string, child: string) => {
        const canvas = pcs()[idx].root.querySelector('canvas[name="canvas"]')!;
        const rect = canvas.getBoundingClientRect();
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 500));
      };
      await clickSub(0, 'Pick Up / Apply', 'Pick Up');
      await clickSub(1, 'Pick Up / Apply', 'Apply');
      await new Promise(r => setTimeout(r, 500));
      const applied = {
        cols: pcs()[1]?.props.columnNames?.slice(),
        color: pcs()[1]?.props.colorColumnName,
        log: pcs()[1]?.props.logColumnsColumnNames?.slice(),
        legend: pcs()[1]?.props.legendPosition,
        title: pcs()[1]?.props.title,
      };

      // Step 8: changing the first plot's axes must NOT touch the second plot.
      pcs()[0].props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      await new Promise(r => setTimeout(r, 600));
      const pc2ColsAfterPc1Change = pcs()[1]?.props.columnNames?.slice();
      const pc1ColsAfterChange = pcs()[0]?.props.columnNames?.slice();

      // Step 9: a range slider on the second plot filters the shared DataFrame,
      // so the first plot updates too (both observe the same df.filter).
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      await new Promise(r => setTimeout(r, 400));
      const fullBefore = df.filter.trueCount;
      const viewers = document.querySelectorAll('[name="viewer-PC-Plot"]');
      const pc2El = viewers[1] as HTMLElement;
      const vr = pc2El.getBoundingClientRect();
      pc2El.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
      await new Promise(r => setTimeout(r, 400));
      const svg = pc2El.querySelector('[name="axis-slider-AGE"]');
      let draggedPc2 = false;
      if (svg) {
        const maxHandle = svg.querySelector('[name="max-handle"]')!;
        const hr = maxHandle.getBoundingClientRect();
        const cx = hr.x + hr.width / 2, cy = hr.y + hr.height / 2;
        const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
        maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await new Promise(r => setTimeout(r, 50));
        for (let dy = 20; dy <= 200; dy += 30) {
          document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          await new Promise(r => setTimeout(r, 20));
        }
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 200)));
        await new Promise(r => setTimeout(r, 500));
        draggedPc2 = true;
      }
      const filteredByPc2 = df.filter.trueCount;

      pcs()[1]?.close();
      const pcF = tv.viewers.find(v => v.type === 'PC Plot')!;
      df.filter.setAll(true);
      pcF.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      pcF.props.logColumnsColumnNames = []; pcF.props.colorColumnName = '';
      pcF.props.legendPosition = 'Auto'; pcF.props.title = '';
      await new Promise(r => setTimeout(r, 300));
      return {applied, pc2ColsAfterPc1Change, pc1ColsAfterChange, fullBefore, filteredByPc2, draggedPc2};
    });
    // Step 7: the second plot matches the first's picked-up state.
    expect(result.applied.cols).toEqual(['AGE', 'WEIGHT', 'STARTED']);
    expect(result.applied.color).toBe('RACE');
    expect(result.applied.log).toEqual(['AGE']);
    expect(result.applied.legend).toBe('Left');
    expect(result.applied.title).toBe('Source Plot');
    // Step 8: the second plot is independent of the first's later axis change.
    expect(result.pc1ColsAfterChange).toEqual(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);
    expect(result.pc2ColsAfterPc1Change).toEqual(['AGE', 'WEIGHT', 'STARTED']);
    // Step 9: a slider on the second plot filters the shared DataFrame.
    expect(result.draggedPc2).toBe(true);
    expect(result.filteredByPc2).toBeLessThan(result.fullBefore);
  });

  await softStep('Layout round-trip', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const hasScatter = tv.viewers.some(v => v.type === 'Scatter plot');
      const hasPc = tv.viewers.some(v => v.type === 'PC Plot');
      await grok.dapi.layouts.delete(saved);
      return { hasScatter, hasPc };
    });
    expect(result.hasScatter).toBe(false);
    expect(result.hasPc).toBe(true);
  });

  // Only the UI Save button captures the VIEW LAYOUT into a project (a JS-API
  // Project.create().addChild(saveLayout()) throws "Unable to add entity"), so
  // the project is saved through the real ribbon Save button, then closeAll +
  // reopen restores the PC plot. Pattern proven in
  // LineChart/legend-color-and-persistence-spec.ts.
  await softStep('Project save / Close All / reopen', async () => {
    const projName = 'zz-pcplot-save-restore-' + Date.now();
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.title = 'Persisted PC Plot';
      await new Promise(r => setTimeout(r, 400));
    });

    await page.locator('[name="button-Save"]').first().click();
    await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
    await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
    await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});
    await page.waitForTimeout(3000);
    // A "Share <project>" dialog pops up after a successful save — dismiss it.
    const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
    if (await cancel.count() > 0) await cancel.click({force: true});
    await page.waitForTimeout(800);

    const result = await page.evaluate(async (name) => {
      let proj = null;
      for (let a = 0; a < 6 && !proj; a++) {
        try { proj = await grok.dapi.projects.filter('name = "' + name + '"').first(); } catch (e) {}
        if (!proj) await new Promise(r => setTimeout(r, 1200));
      }
      if (!proj) return {found: false};
      try {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(proj.id);
        await full.open();
        await new Promise(r => setTimeout(r, 4500));
        const tv = grok.shell.tv;
        const pcRestored = (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'PC Plot');
        const pc = tv ? Array.from(tv.viewers).find((x: any) => x.type === 'PC Plot') as any : null;
        const titleRestored = pc?.props?.title;
        return {found: true, pcRestored, titleRestored};
      } finally {
        await grok.dapi.projects.delete(proj); // never leak the probe project
      }
    }, projName);

    expect(result.found).toBe(true);
    expect(result.pcRestored).toBe(true);
    expect(result.titleRestored).toBe('Persisted PC Plot');
  });

  await softStep('Table switching and transformation', async () => {
    const result = await page.evaluate(async (path) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      const df2 = await grok.dapi.files.readCsv(path);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 3000);
      });
      const hasBioChem = Array.from({length: df2.columns.length}, (_, i) => df2.columns.byIndex(i))
        .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      // The PC plot is added on the spgi view; setting the Table property binds
      // it to spgi (the axes below are spgi columns, confirmed by the pivot).
      const pc = grok.shell.tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      pc.props.table = df2.name;
      await new Promise(r => setTimeout(r, 500));
      const tableSet = pc.dataFrame?.name;
      const sliderAxes = () =>
        Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
          .map((e) => e.getAttribute('name')!.replace('axis-slider-', ''));
      await new Promise(r => setTimeout(r, 1500));
      const axesBefore = sliderAxes();
      pc.props.transformation = '[{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]';
      await new Promise(r => setTimeout(r, 3000));
      // The pivot replaces the axes with one generated column per Series value, so
      // the slider names show whether the aggregation was applied.
      const axesAfter = sliderAxes();
      pc.props.transformation = '';
      await new Promise(r => setTimeout(r, 1500));
      const axesReverted = sliderAxes();
      pc.close();
      return { spgiRows: df2.rowCount, tableSet, axesBefore, axesAfter, axesReverted };
    }, spgiPath);
    expect(result.spgiRows).toBe(100);
    expect(result.tableSet).toBeTruthy();
    // Pivoted axes are the Series categories, not the raw numeric columns.
    expect(result.axesAfter).not.toEqual(result.axesBefore);
    expect(result.axesAfter).toContain('Triazoles');
    expect(result.axesReverted).toEqual(result.axesBefore);
  });

  v.finishSpec();
});
