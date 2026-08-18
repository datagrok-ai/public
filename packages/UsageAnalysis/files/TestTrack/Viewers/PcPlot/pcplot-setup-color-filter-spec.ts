/* ---
realizes: [pcplot.cp.setup-columns-color-filter, pcplot.cp.layout-project-persistence, pcplot.int.color-column-legend-coding, pcplot.int.range-filter-cross-viewer]
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

  await v.addViewerByIcon(page, 'pc-plot', 'PC-Plot', 15000);

  const axisNames = (): Promise<string[]> => page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
      .map((e) => e.getAttribute('name')!.replace('axis-slider-', '')));
  const legendCount = (): Promise<number> => page.evaluate(() =>
    document.querySelectorAll('[name="viewer-PC-Plot"] .d4-legend').length);
  const legendLabels = (): Promise<string[] | null> => page.evaluate(() => {
    const el = document.querySelector('[name="viewer-PC-Plot"] .d4-legend');
    return el
      ? Array.from(el.querySelectorAll('.d4-legend-value')).map((e) => (e.textContent ?? '').trim())
      : null;
  });

  await softStep('Column setup — select AGE, HEIGHT, WEIGHT (axis count = 3)', async () => {

    await v.setViewerProps(page, 'PC Plot', [{set: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}, wait: 900}]);

    const axes = await v.pollValue(axisNames, (a) => a.length === 3, 3000, 150);
    expect(axes).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
  });

  await softStep('In-chart range-filter drop + Reset View restore (PRIMARY SIGNAL)', async () => {
    const filterCount = (): Promise<number> =>
      page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();

      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
    });

    await v.pollValue(() => page.evaluate(() => {
      const h = document.querySelector('[name="axis-slider-AGE"] [name="max-handle"]');
      if (!h) return 0;
      const r = h.getBoundingClientRect();
      return Math.min(r.width, r.height);
    }), (size) => size > 0, 400, 100);
    const fullCount = await filterCount();

    await page.evaluate(async () => {

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
    });
    const filteredCount = await v.pollValue(filterCount, (n) => n < fullCount, 600, 100);

    await page.evaluate(() => {
      const canvas = document.querySelector('[name="viewer-PC-Plot"] canvas[name="canvas"]')!;
      const crect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: crect.left + crect.width / 2, clientY: crect.top + crect.height / 2,
      }));
    });
    await v.pollValue(
      () => page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .some((el) => el.textContent!.trim() === 'Reset View')),
      (present) => present, 500, 100);
    await page.evaluate(() => {
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const rv = items.find((el) => el.textContent!.trim() === 'Reset View');
      if (rv)
        (rv.closest('.d4-menu-item') as HTMLElement).click();
    });
    const restoredCount = await v.pollValue(filterCount, (n) => n === fullCount, 600, 100);

    expect(filteredCount).toBeLessThan(fullCount);
    expect(restoredCount).toBe(fullCount);
  });

  await softStep('GROK-18000 — add then remove a column, axes update immediately (DOM axis-slider count 3 → 4 → 3)', async () => {

    const errBefore = pageErrors.length + consoleErrors.length;
    const base: string[] = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      return pc.props.columnNames.slice();
    });
    await v.setViewerProps(page, 'PC Plot', [{set: {columnNames: [...base, 'STARTED']}, wait: 500}]);
    const afterAdd = await v.pollValue(axisNames, (a) => a.length === base.length + 1, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {columnNames: base}, wait: 500}]);
    const afterRemove = await v.pollValue(axisNames, (a) => a.length === base.length, 3000, 150);
    expect(afterAdd.length).toBe(base.length + 1);
    expect(afterAdd).toContain('STARTED');
    expect(afterRemove.length).toBe(base.length);
    expect(afterRemove).not.toContain('STARTED');
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('GROK-17754 — color by HEIGHT, switch coloring type, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: 'HEIGHT'}, wait: 400}]);
    await page.evaluate(() => { grok.shell.tv.dataFrame.col('HEIGHT').meta.colors.setCategorical(); });
    await v.waitForViewerRendered(page, 'PC Plot', 400);
    await page.evaluate(() => { grok.shell.tv.dataFrame.col('HEIGHT').meta.colors.setLinear(); });
    await v.waitForViewerRendered(page, 'PC Plot', 400);
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: ''}, wait: 400}]);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Show Filtered Out Lines toggle, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'PC Plot', [
      {set: {showFilteredOutLines: true}, wait: 400},
      {set: {showFilteredOutLines: false}, wait: 300},
    ]);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Per-column logarithmic scale for AGE, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'PC Plot', [
      {set: {logColumnsColumnNames: ['AGE']}, wait: 400},
      {set: {logColumnsColumnNames: []}, wait: 300},
    ]);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Categorical coloring renders a legend listing RACE categories', async () => {
    const raceCats: string[] = await page.evaluate(() => grok.shell.tv.dataFrame.col('RACE').categories.slice());
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: ''}, wait: 400}]);
    const legendBefore = await v.pollValue(legendCount, (n) => n === 0, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: 'RACE'}, wait: 800}]);

    const shown = await v.pollValue(() => page.evaluate(() => {
      const legends = document.querySelectorAll('[name="viewer-PC-Plot"] .d4-legend');
      return {
        legendAfterCount: legends.length,
        legendText: legends.length ? legends[0].textContent : '',
        legendValues: legends.length
          ? Array.from(legends[0].querySelectorAll('.d4-legend-value')).map((e) => (e.textContent ?? '').trim())
          : [],
      };
    }), (s) => s.legendValues.length === raceCats.length, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: ''}, wait: 800}]);
    const legendAfterClear = await v.pollValue(legendCount, (n) => n === 0, 3000, 150);

    expect(legendBefore).toBe(0);
    expect(shown.legendAfterCount).toBeGreaterThan(0);
    for (const cat of raceCats)
      expect(shown.legendText).toContain(cat);

    expect([...shown.legendValues].sort()).toEqual([...raceCats].sort());
    expect(legendAfterClear).toBe(0);
  });

  await softStep('Numeric coloring gradient drive — invert via the context menu, log axis, min/max clamp', async () => {

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

      pc.props.colorAxisType = 'logarithmic'; await wait();
      pc.props.colorMin = 30; pc.props.colorMax = 60; await wait();
      pc.props.colorMin = null; pc.props.colorMax = null;
      pc.props.colorAxisType = 'linear';
      pc.props.colorColumnName = '';
      await wait();
      return {invertBefore, invClicked1, invertToggled, invClicked2, invertRestored,
        editClicked, dialogHeader, dialogClosed};
    });

    expect(menu.invClicked1).toBe(true);
    expect(menu.invClicked2).toBe(true);
    expect(menu.invertToggled).toBe(!menu.invertBefore);
    expect(menu.invertRestored).toBe(menu.invertBefore);

    expect(menu.editClicked).toBe(true);
    expect(menu.dialogHeader).toContain('Color-coding: AGE');
    expect(menu.dialogClosed).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Legend position cycle and visibility round-trip — Never removes the legend, Auto restores the same labels', async () => {
    const labelled = (l: string[] | null) => l !== null && l.length > 0;
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: 'RACE'}, wait: 800}]);
    const initial = await v.pollValue(legendLabels, labelled, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [
      {set: {legendPosition: 'Left'}, wait: 300},
      {set: {legendPosition: 'Right'}, wait: 0},
      {set: {legendPosition: 'Top'}, wait: 0},
      {set: {legendPosition: 'Bottom'}, wait: 300},
    ]);
    const afterCycle = await v.pollValue(legendLabels, labelled, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {legendVisibility: 'Never'}, wait: 600}]);
    const hiddenCount = await v.pollValue(legendCount, (n) => n === 0, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {legendVisibility: 'Auto'}, wait: 600}]);
    const restored = await v.pollValue(legendLabels, labelled, 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {legendPosition: 'Auto', colorColumnName: ''}, wait: 600}]);
    const clearedCount = await v.pollValue(legendCount, (n) => n === 0, 3000, 150);

    expect(initial).not.toBeNull();
    expect(initial!.length).toBeGreaterThan(0);

    expect(afterCycle).toEqual(initial);

    expect(hiddenCount).toBe(0);
    expect(restored).toEqual(initial);
    expect(clearedCount).toBe(0);
  });

  await softStep('Grid conditional color coding surfaces its bins in the plot legend; linear drops it', async () => {

    const legend = () => page.evaluate(() => {
      const el = document.querySelector('[name="viewer-PC-Plot"] .d4-legend') as HTMLElement | null;
      return {present: !!el, text: el ? (el.innerText || '').replace(/\s+/g, ' ').trim() : ''};
    });
    await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.colorColumnName = 'HEIGHT';
      grok.shell.tv.dataFrame.col('HEIGHT').meta.colors
        .setConditional({'20-150': DG.Color.green, '150-250': DG.Color.orange});
    });
    await v.waitForViewerRendered(page, 'PC Plot', 800);
    const conditional = await v.pollValue(legend, (l) => l.present && l.text.includes('150-250'), 3000, 150);
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.col('HEIGHT').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
    });
    await v.waitForViewerRendered(page, 'PC Plot', 800);
    const linear = await v.pollValue(legend, (l) => !l.present, 3000, 150);

    await page.evaluate(() => { grok.shell.tv.dataFrame.col('HEIGHT').meta.colors.setLinear(); });
    await v.setViewerProps(page, 'PC Plot', [{set: {colorColumnName: ''}, wait: 400}]);
    const cleared = await v.pollValue(legend, (l) => !l.present, 3000, 150);

    expect(conditional.present).toBe(true);
    expect(conditional.text).toContain('20-150');
    expect(conditional.text).toContain('150-250');
    expect(linear.present).toBe(false);
    expect(cleared.present).toBe(false);
  });

  await softStep('Layout round-trip — saved layout restores the configured viewer set and props', async () => {

    await v.setViewerProps(page, 'PC Plot', [{
      set: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT'], colorColumnName: 'RACE', title: 'PC Persistence Probe'},
      wait: 800,
    }]);
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      await page.waitForTimeout(1000); 
      await page.evaluate(() => { grok.shell.tv.addViewer('Scatter plot'); });
      await v.pollValue(
        () => page.evaluate(() => grok.shell.tv.viewers.some((vw: any) => vw.type === 'Scatter plot')),
        (present) => present, 500, 100);
      await page.evaluate(async (id) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);
      const result = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot');
        return {
          hasScatter: tv.viewers.some((vw: any) => vw.type === 'Scatter plot'),
          hasPc: tv.viewers.some((vw: any) => vw.type === 'PC Plot'),
          cols: pc?.props.columnNames?.slice(),
          color: pc?.props.colorColumnName,
          title: pc?.props.title,
        };
      }), (r) => r.hasPc && !r.hasScatter, 3000, 150);

      expect(result.hasScatter).toBe(false);
      expect(result.hasPc).toBe(true);

      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {

      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved)
            await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Project save / Close All / reopen — project restores the configured viewer', async () => {
    const projName = 'zz-pcplot-persistence-probe-' + Date.now();
    try {
      await page.locator('[name="button-Save"]').first().click();
      await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
      await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
      await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});

      const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
      await v.pollValue(() => cancel.count(), (n) => n > 0, 3000, 150);
      if (await cancel.count() > 0)
        await cancel.click({force: true});
      await v.pollValue(() => page.locator('.d4-dialog').count(), (n) => n === 0, 800, 100);

      const projId = await v.pollValue(() => page.evaluate(async (name) => {
        try {
          return (await grok.dapi.projects.filter('name = "' + name + '"').first())?.id ?? null;
        } catch (e) {
          return null;
        }
      }, projName), (id) => id !== null, 7200, 1200);

      const result: any = {found: projId !== null};
      if (projId !== null) {
        await v.closeAllAndWait(page);
        await page.evaluate(async (id) => {
          const full = await grok.dapi.projects.find(id);
          await full.open();
        }, projId);
        Object.assign(result, await v.pollValue(() => page.evaluate(() => {
          const tv = grok.shell.tv;
          const pc = tv ? Array.from(tv.viewers).find((x: any) => x.type === 'PC Plot') as any : null;
          return {
            pcRestored: (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'PC Plot'),
            cols: pc?.props?.columnNames?.slice(),
            color: pc?.props?.colorColumnName,
            title: pc?.props?.title,
          };
        }), (r) => r.pcRestored, 4500, 150));
      }

      expect(result.found).toBe(true);
      expect(result.pcRestored).toBe(true);

      expect(result.cols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.color).toBe('RACE');
      expect(result.title).toBe('PC Persistence Probe');
    } finally {

      await page.evaluate(async (name) => {
        try {
          const leftover = await grok.dapi.projects.filter('name = "' + name + '"').first();
          if (leftover)
            await grok.dapi.projects.delete(leftover);
        } catch (_) {}
      }, projName);
    }
  });

  await page.evaluate(() => {
    const pc = grok.shell.tv?.viewers?.find((vw: any) => vw.type === 'PC Plot');
    if (pc) {
      pc.props.colorColumnName = '';
      pc.props.title = '';
    }
  });
  await v.waitForViewerRendered(page, 'PC Plot', 300);

  v.finishSpec();
});
