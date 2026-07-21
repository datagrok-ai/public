/* ---
realizes: [pcplot.cp.setup-columns-color-filter]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


declare const grok: any;


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
    // Select-columns list is canvas-rendered — 2026-07-21 recon).
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

  await softStep('GROK-18000 — add then remove a column, no error (no-error floor)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const base = pc.props.columnNames.slice();
      pc.props.columnNames = [...base, 'STARTED'];
      await new Promise((r) => setTimeout(r, 400));
      pc.props.columnNames = base;
      await new Promise((r) => setTimeout(r, 400));
    });
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
      pc.props.colorColumnName = '';
      await new Promise((r) => setTimeout(r, 800));
      const legendAfterClear = root.querySelectorAll('.d4-legend').length;
      return {raceCats, legendBefore, legendAfterCount: legends.length, legendText, legendAfterClear};
    });
    expect(result.legendBefore).toBe(0);
    expect(result.legendAfterCount).toBeGreaterThan(0);
    for (const cat of result.raceCats)
      expect(result.legendText).toContain(cat);
    expect(result.legendAfterClear).toBe(0);
  });

  v.finishSpec();
});
