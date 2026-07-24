/* ---
realizes: [pcplot.cp.setup-columns-color-filter, pcplot.cp.layout-project-persistence, pcplot-color-column-legend-coding, pcplot-range-filter-cross-viewer]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


declare const grok: any;
declare const DG: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';


test('PC Plot — Setup, Column Selection, Color, In-Chart Range Filter, Log Scale', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error')
      consoleErrors.push(m.text());
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  await softStep('Column setup — select AGE, HEIGHT, WEIGHT (axis count = 3)', async () => {
    // Cross-channel signal: columnNames is WRITTEN via props but read back from
    // the RENDERED DOM (one axis-slider element per column), so a broken
    // re-render fails instead of echoing the prop value. Setting the columns
    // through the Context Panel dialog is not scriptable headless (the
    // Select-columns list is canvas-rendered).
    const axes = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise((r) => setTimeout(r, 900));
      return Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
        .map((e) => e.getAttribute('name')!.replace('axis-slider-', ''));
    });
    expect(axes).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
  });

  await softStep('In-chart range-filter drop + Reset View restore (PRIMARY SIGNAL)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();
      // reveal the per-axis range sliders (hidden until hover)
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise((r) => setTimeout(r, 400));
      const fullCount = df.filter.trueCount;
      // The AGE axis slider is real DOM — <svg name="axis-slider-AGE"> with
      // min-handle / max-handle circles that track standard mouse events. Dragging
      // the max handle downward narrows the AGE window.
      const svg = document.querySelector('[name="axis-slider-AGE"]')!;
      const maxHandle = svg.querySelector('[name="max-handle"]')!;
      const hr = maxHandle.getBoundingClientRect();
      const cx = hr.x + hr.width / 2;
      const cy = hr.y + hr.height / 2;
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
      await new Promise((r) => setTimeout(r, 50));
      for (let dy = 20; dy <= 300; dy += 30) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
        await new Promise((r) => setTimeout(r, 20));
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 300)));
      await new Promise((r) => setTimeout(r, 600));
      const filteredCount = df.filter.trueCount;
      // Reset View via the context menu (fully restores the filter).
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const crect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: crect.left + crect.width / 2, clientY: crect.top + crect.height / 2,
      }));
      await new Promise((r) => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const rv = items.find((el) => el.textContent!.trim() === 'Reset View');
      if (rv)
        (rv.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 600));
      const restoredCount = df.filter.trueCount;
      return {fullCount, filteredCount, restoredCount};
    });
    expect(result.filteredCount).toBeLessThan(result.fullCount);
    expect(result.restoredCount).toBe(result.fullCount);
  });

  await softStep('GROK-18000 — add then remove a column, axes update immediately (DOM axis-slider count 3 → 4 → 3)', async () => {
    // GROK-18000: the axes must update immediately with no manual refresh. The
    // read-back is the RENDERED axis-slider set, not the prop echo: adding STARTED
    // (a valid DateTime axis — it renders axis-slider-STARTED)
    // grows the slider set to four, removing it returns to three.
    const errBefore = pageErrors.length + consoleErrors.length;
    const names = (): Promise<string[]> => page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
        .map((e) => e.getAttribute('name')!.replace('axis-slider-', '')));
    const base = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      return pc.props.columnNames.slice();
    });
    await page.evaluate(async (b) => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.columnNames = [...b, 'STARTED'];
      await new Promise((r) => setTimeout(r, 500));
    }, base);
    const afterAdd = await names();
    await page.evaluate(async (b) => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.columnNames = b;
      await new Promise((r) => setTimeout(r, 500));
    }, base);
    const afterRemove = await names();
    expect(afterAdd.length).toBe(base.length + 1);
    expect(afterAdd).toContain('STARTED');
    expect(afterRemove.length).toBe(base.length);
    expect(afterRemove).not.toContain('STARTED');
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('GROK-17754 — color by HEIGHT, switch coloring type, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const df = grok.shell.tv.dataFrame;
      pc.props.colorColumnName = 'HEIGHT';
      await new Promise((r) => setTimeout(r, 400));
      df.col('HEIGHT').meta.colors.setCategorical();
      await new Promise((r) => setTimeout(r, 400));
      df.col('HEIGHT').meta.colors.setLinear();
      await new Promise((r) => setTimeout(r, 400));
      pc.props.colorColumnName = '';
      await new Promise((r) => setTimeout(r, 400));
    });
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Show Filtered Out Lines toggle, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.showFilteredOutLines = true;
      await new Promise((r) => setTimeout(r, 400));
      pc.props.showFilteredOutLines = false;
      await new Promise((r) => setTimeout(r, 300));
    });
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Per-column logarithmic scale for AGE, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.logColumnsColumnNames = ['AGE'];
      await new Promise((r) => setTimeout(r, 400));
      pc.props.logColumnsColumnNames = [];
      await new Promise((r) => setTimeout(r, 300));
    });
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Categorical coloring renders a legend listing RACE categories', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const root = document.querySelector('[name="viewer-PC-Plot"]')!;
      const df = grok.shell.tv.dataFrame;
      const raceCats = df.col('RACE').categories.slice();
      pc.props.colorColumnName = '';
      await new Promise((r) => setTimeout(r, 400));
      const legendBefore = root.querySelectorAll('.d4-legend').length;
      pc.props.colorColumnName = 'RACE';
      await new Promise((r) => setTimeout(r, 800));
      const legends = root.querySelectorAll('.d4-legend');
      const legendText = legends.length ? legends[0].textContent : '';
      // Per-item labels from `.d4-legend-value` — the clean category strings, so a
      // set comparison catches EXTRA entries the `toContain` loop would miss.
      const legendValues = legends.length
        ? Array.from(legends[0].querySelectorAll('.d4-legend-value')).map((e) => (e.textContent ?? '').trim())
        : [];
      pc.props.colorColumnName = '';
      await new Promise((r) => setTimeout(r, 800));
      const legendAfterClear = root.querySelectorAll('.d4-legend').length;
      return {raceCats, legendBefore, legendAfterCount: legends.length, legendText, legendValues, legendAfterClear};
    });
    expect(result.legendBefore).toBe(0);
    expect(result.legendAfterCount).toBeGreaterThan(0);
    for (const cat of result.raceCats)
      expect(result.legendText).toContain(cat);
    // Exactly the RACE categories, no extras — set equality both ways.
    expect([...result.legendValues].sort()).toEqual([...result.raceCats].sort());
    expect(result.legendAfterClear).toBe(0);
  });

  await softStep('Numeric coloring gradient drive — invert via the context menu, log axis, min/max clamp', async () => {
    // With a numeric color column the right-click menu carries a Color Scheme
    // group — the menu-reachable part of the gradient surface. Its Invert
    // Color Scheme item flips the invertColorScheme state (asserted as
    // menu -> state) and Edit... opens the color-coding dialog. The label
    // "Color Scheme" also names the scheme picker inside the group, so the
    // group is resolved by its d4-menu-group class and the child clicks are
    // scoped to it. The remaining gradient options (log axis, min/max clamps)
    // repaint the canvas only and stay a no-error floor over the prop drive.
    const errBefore = pageErrors.length + consoleErrors.length;
    const menu = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const wait = (ms = 300) => new Promise((r) => setTimeout(r, ms));
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2}));
        await wait(500);
      };
      const closeMenu = async () => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        await wait(300);
      };
      const findColorGroup = () => {
        for (const lbl of Array.from(document.querySelectorAll('.d4-menu-item-label'))) {
          if (lbl.textContent!.trim() !== 'Color Scheme') continue;
          const item = lbl.closest('.d4-menu-item')!;
          if (item.classList.contains('d4-menu-group')) return item;
        }
        return null;
      };
      const clickColorSub = async (child: string) => {
        await openMenu();
        const group = findColorGroup();
        if (!group) { await closeMenu(); return false; }
        group.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
        group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await wait(350);
        const c = Array.from(group.querySelectorAll('.d4-menu-item-label'))
          .find((el) => el.textContent!.trim() === child);
        if (c) { (c.closest('.d4-menu-item') as HTMLElement).click(); await wait(400); }
        else await closeMenu();
        return !!c;
      };
      pc.props.colorColumnName = 'AGE';
      await wait(600);
      const invertBefore = pc.props.invertColorScheme;
      const invClicked1 = await clickColorSub('Invert Color Scheme');
      const invertToggled = pc.props.invertColorScheme;
      const invClicked2 = await clickColorSub('Invert Color Scheme');
      const invertRestored = pc.props.invertColorScheme;
      const editClicked = await clickColorSub('Edit...');
      await wait(500);
      const dlg = document.querySelector('.d4-dialog');
      const dialogHeader = dlg ? ((dlg.querySelector('.d4-dialog-header') as HTMLElement)?.innerText ?? '').trim() : '';
      const closeBtn = dlg?.querySelector('button[name="button-CLOSE"]') as HTMLElement | null;
      if (closeBtn)
        closeBtn.click();
      await wait(500);
      const dialogClosed = !document.querySelector('.d4-dialog');
      // Remaining gradient options: canvas repaint only — no-error floor.
      pc.props.colorAxisType = 'logarithmic'; await wait();
      pc.props.colorMin = 30; pc.props.colorMax = 60; await wait();
      pc.props.colorMin = null; pc.props.colorMax = null;
      pc.props.colorAxisType = 'linear';
      pc.props.colorColumnName = '';
      await wait();
      return {invertBefore, invClicked1, invertToggled, invClicked2, invertRestored,
        editClicked, dialogHeader, dialogClosed};
    });
    // Menu -> state: the read prop is the one the MENU item flipped, with a
    // second click restoring it (round-trip).
    expect(menu.invClicked1).toBe(true);
    expect(menu.invClicked2).toBe(true);
    expect(menu.invertToggled).toBe(!menu.invertBefore);
    expect(menu.invertRestored).toBe(menu.invertBefore);
    // Edit... opens the color-coding dialog for the color column.
    expect(menu.editClicked).toBe(true);
    expect(menu.dialogHeader).toContain('Color-coding: AGE');
    expect(menu.dialogClosed).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Legend position cycle and visibility round-trip — Never removes the legend, Auto restores the same labels', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const root = document.querySelector('[name="viewer-PC-Plot"]')!;
      const wait = (ms = 300) => new Promise((r) => setTimeout(r, ms));
      const legendLabels = (): string[] | null => {
        const el = root.querySelector('.d4-legend');
        return el
          ? Array.from(el.querySelectorAll('.d4-legend-value')).map((e) => (e.textContent ?? '').trim())
          : null;
      };
      pc.props.colorColumnName = 'RACE';
      await wait(800);
      const initial = legendLabels();
      pc.props.legendPosition = 'Left'; await wait();
      pc.props.legendPosition = 'Right'; pc.props.legendPosition = 'Top';
      pc.props.legendPosition = 'Bottom'; await wait();
      const afterCycle = legendLabels();
      pc.props.legendVisibility = 'Never';
      await wait(600);
      const hiddenCount = root.querySelectorAll('.d4-legend').length;
      pc.props.legendVisibility = 'Auto';
      await wait(600);
      const restored = legendLabels();
      pc.props.legendPosition = 'Auto';
      pc.props.colorColumnName = '';
      await wait(600);
      const clearedCount = root.querySelectorAll('.d4-legend').length;
      return {initial, afterCycle, hiddenCount, restored, clearedCount};
    });
    expect(result.initial).not.toBeNull();
    expect(result.initial!.length).toBeGreaterThan(0);
    // The legend survives the position cycle with the same labels.
    expect(result.afterCycle).toEqual(result.initial);
    // Never removes the legend element from the DOM; Auto brings it back
    // with the SAME labels (round-trip).
    expect(result.hiddenCount).toBe(0);
    expect(result.restored).toEqual(result.initial);
    expect(result.clearedCount).toBe(0);
  });

  await softStep('Grid conditional color coding surfaces its bins in the plot legend; linear drops it', async () => {
    // Conditional color coding set on the grid column renders a DOM legend
    // listing its bins, while a linear/numeric gradient has no DOM legend —
    // that contrast is the readable signal that the plot picked up the grid
    // column's color-coding change.
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const df = grok.shell.tv.dataFrame;
      const root = document.querySelector('[name="viewer-PC-Plot"]')!;
      const wait = (ms = 800) => new Promise((r) => setTimeout(r, ms));
      const legend = () => {
        const el = root.querySelector('.d4-legend') as HTMLElement | null;
        return {present: !!el, text: el ? (el.innerText || '').replace(/\s+/g, ' ').trim() : ''};
      };
      pc.props.colorColumnName = 'HEIGHT';
      df.col('HEIGHT').meta.colors.setConditional({'20-150': DG.Color.green, '150-250': DG.Color.orange});
      await wait();
      const conditional = legend();
      df.col('HEIGHT').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
      await wait();
      const linear = legend();
      // Clean up: reset the grid column coloring and clear the color column.
      df.col('HEIGHT').meta.colors.setLinear();
      pc.props.colorColumnName = '';
      await wait(400);
      const cleared = legend();
      return {conditional, linear, cleared};
    });
    expect(result.conditional.present).toBe(true);
    expect(result.conditional.text).toContain('20-150');
    expect(result.conditional.text).toContain('150-250');
    expect(result.linear.present).toBe(false);
    expect(result.cleared.present).toBe(false);
  });

  // Trailing persistence check: the two round-trips below persist the state
  // this spec just configured (axes AGE/HEIGHT/WEIGHT, color RACE, probe
  // title) — restored-set membership plus restored props, not a saved-layout
  // prop echo.
  await softStep('Layout round-trip — saved layout restores the configured viewer set and props', async () => {
    // Save the layout first and capture its id on the node side, so the
    // finally teardown can delete the probe layout even when a later
    // evaluate or assertion fails.
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      pc.props.colorColumnName = 'RACE';
      pc.props.title = 'PC Persistence Probe';
      await new Promise((r) => setTimeout(r, 800));
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        await new Promise((r) => setTimeout(r, 1000));
        tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 500));
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 3000));
        const hasScatter = tv.viewers.some((vw: any) => vw.type === 'Scatter plot');
        const hasPc = tv.viewers.some((vw: any) => vw.type === 'PC Plot');
        const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot');
        return {
          hasScatter, hasPc,
          cols: pc?.props.columnNames?.slice(),
          color: pc?.props.colorColumnName,
          title: pc?.props.title,
        };
      }, layoutId);
      // The restored viewer set equals the SAVED set…
      expect(result.hasScatter).toBe(false);
      expect(result.hasPc).toBe(true);
      // …and the restored PC Plot carries the configured props.
      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {
      // Never leak the probe layout, even when an assertion above failed.
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved)
            await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  // The project is saved through the real ribbon Save button (only the UI Save
  // captures the view layout), then closeAll + reopen restores the PC plot.
  await softStep('Project save / Close All / reopen — project restores the configured viewer', async () => {
    const projName = 'zz-pcplot-persistence-probe-' + Date.now();
    try {
      await page.locator('[name="button-Save"]').first().click();
      await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
      await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
      await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});
      await page.waitForTimeout(3000);
      // A "Share <project>" dialog pops up after a successful save — dismiss it.
      const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
      if (await cancel.count() > 0)
        await cancel.click({force: true});
      await page.waitForTimeout(800);

      const result = await page.evaluate(async (name) => {
        let proj = null;
        for (let a = 0; a < 6 && !proj; a++) {
          try {
            proj = await grok.dapi.projects.filter('name = "' + name + '"').first();
          } catch (e) {}
          if (!proj)
            await new Promise((r) => setTimeout(r, 1200));
        }
        if (!proj)
          return {found: false};
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(proj.id);
        await full.open();
        await new Promise((r) => setTimeout(r, 4500));
        const tv = grok.shell.tv;
        const pcRestored = (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'PC Plot');
        const pc = tv ? Array.from(tv.viewers).find((x: any) => x.type === 'PC Plot') as any : null;
        return {
          found: true, pcRestored,
          cols: pc?.props?.columnNames?.slice(),
          color: pc?.props?.colorColumnName,
          title: pc?.props?.title,
        };
      }, projName);

      expect(result.found).toBe(true);
      expect(result.pcRestored).toBe(true);
      // The reopened project restores the configured props — a cross-session
      // round-trip, distinct from setting the props.
      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {
      // Never leak the probe project, even when the save/reopen flow failed.
      await page.evaluate(async (name) => {
        try {
          const leftover = await grok.dapi.projects.filter('name = "' + name + '"').first();
          if (leftover)
            await grok.dapi.projects.delete(leftover);
        } catch (_) {}
      }, projName);
    }
  });

  // Cleanup: strip the probe styling so the spec ends on the baseline state.
  await page.evaluate(async () => {
    const pc = grok.shell.tv?.viewers?.find((vw: any) => vw.type === 'PC Plot');
    if (pc) {
      pc.props.colorColumnName = '';
      pc.props.title = '';
    }
    await new Promise((r) => setTimeout(r, 300));
  });

  v.finishSpec();
});
