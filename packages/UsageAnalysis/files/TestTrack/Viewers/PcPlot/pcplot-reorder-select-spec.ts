/* ---
realizes: [pcplot.cp.reorder-and-select]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';


test('PC Plot — Axis Reorder, Polyline Selection, and Current-Row Sync', async ({page}) => {
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

  await page.evaluate(async () => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
    await new Promise((r) => setTimeout(r, 800));
  });

  await softStep('Scenario 1 Step 4 — shift+drag rectangle selects polylines (selection rises above zero)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 200));
      const before = df.selection.trueCount;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number, extra: any) => Object.assign(
        {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0}, extra || {});
      // shift+drag a rectangle in the band between the first two axes (AGE/HEIGHT)
      const x1 = r0.x + r0.width * 0.32, y1 = r0.y + r0.height * 0.40;
      const x2 = r0.x + r0.width * 0.45, y2 = r0.y + r0.height * 0.55;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(x1, y1, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 30));
      for (let t = 0; t <= 1.0001; t += 0.25) {
        const x = x1 + (x2 - x1) * t, y = y1 + (y2 - y1) * t;
        overlay.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        await new Promise((r) => setTimeout(r, 25));
      }
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      document.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 400));
      return {before, after: df.selection.trueCount};
    });
    expect(result.before).toBe(0);
    expect(result.after).toBeGreaterThan(0);
  });

  await softStep('Scenario 1 Step 5 — additive second shift+drag band (selection rises again, not replaced)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const beforeSecond = df.selection.trueCount;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number, extra: any) => Object.assign(
        {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0}, extra || {});
      // A second band, this one in the HEIGHT/WEIGHT region.
      const x1 = r0.x + r0.width * 0.58, y1 = r0.y + r0.height * 0.42;
      const x2 = r0.x + r0.width * 0.72, y2 = r0.y + r0.height * 0.58;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(x1, y1, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 30));
      for (let t = 0; t <= 1.0001; t += 0.25) {
        const x = x1 + (x2 - x1) * t, y = y1 + (y2 - y1) * t;
        overlay.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        await new Promise((r) => setTimeout(r, 25));
      }
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      document.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 400));
      return {beforeSecond, afterSecond: df.selection.trueCount};
    });
    expect(result.afterSecond).toBeGreaterThan(result.beforeSecond);
  });

  await softStep('Scenario 1 Step 6 — click empty space clears the selection (round-trip to zero)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const beforeClear = df.selection.trueCount;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      // Click the empty top margin, above the polylines, to clear the selection.
      const ex = r0.x + r0.width * 0.5, ey = r0.y + r0.height * 0.02;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('click', mk(ex, ey)));
      await new Promise((r) => setTimeout(r, 400));
      return {beforeClear, afterClear: df.selection.trueCount};
    });
    expect(result.beforeClear).toBeGreaterThan(0);
    expect(result.afterClear).toBe(0);
  });

  await softStep('Scenario 1 Step 7 — click a polyline sets current row off -1', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = -1;
      await new Promise((r) => setTimeout(r, 200));
      const before = df.currentRowIdx;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      // Click a point in the mid band where polylines cross between the first two axes.
      const cx = r0.x + r0.width * 0.42, cy = r0.y + r0.height * 0.45;
      overlay.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy)));
      await new Promise((r) => setTimeout(r, 60));
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy)));
      overlay.dispatchEvent(new MouseEvent('click', mk(cx, cy)));
      await new Promise((r) => setTimeout(r, 300));
      return {before, after: df.currentRowIdx};
    });
    expect(result.before).toBe(-1);
    expect(result.after).toBeGreaterThanOrEqual(0);
  });

  await softStep('Scenario 2 Step 4 — drag a column label reorders the axes (columnNames order changes, still 3)', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const before = pc.getOptions().look.columnNames.slice();
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      // Column-name labels are canvas-rendered near the top; dragging the rightmost
      // label to the leftmost position reorders the axes.
      const labelY = r0.y + 8;
      const leftX = r0.x + r0.width * 0.08;
      const rightX = r0.x + r0.width * 0.88;
      overlay.dispatchEvent(new MouseEvent('mousemove', mk(rightX, labelY)));
      await new Promise((r) => setTimeout(r, 60));
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(rightX, labelY)));
      await new Promise((r) => setTimeout(r, 40));
      for (let t = 0; t <= 1.0001; t += 0.1) {
        const x = rightX + (leftX - rightX) * t;
        overlay.dispatchEvent(new MouseEvent('mousemove', mk(x, labelY)));
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, labelY)));
        await new Promise((r) => setTimeout(r, 30));
      }
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(leftX, labelY)));
      document.dispatchEvent(new MouseEvent('mouseup', mk(leftX, labelY)));
      await new Promise((r) => setTimeout(r, 500));
      const after = pc.getOptions().look.columnNames.slice();
      return {before, after};
    });
    expect(result.after.length).toBe(3);
    expect(result.after).not.toEqual(result.before);
    expect([...result.after].sort()).toEqual([...result.before].sort());
  });

  await softStep('Scenario 2 Step 7 — shift+drag on the reordered chart selects (selection rises above zero)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 200));
      const before = df.selection.trueCount;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number, extra: any) => Object.assign(
        {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0}, extra || {});
      const x1 = r0.x + r0.width * 0.32, y1 = r0.y + r0.height * 0.40;
      const x2 = r0.x + r0.width * 0.45, y2 = r0.y + r0.height * 0.55;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(x1, y1, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 30));
      for (let t = 0; t <= 1.0001; t += 0.25) {
        const x = x1 + (x2 - x1) * t, y = y1 + (y2 - y1) * t;
        overlay.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, y, {shiftKey: true})));
        await new Promise((r) => setTimeout(r, 25));
      }
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      document.dispatchEvent(new MouseEvent('mouseup', mk(x2, y2, {shiftKey: true})));
      await new Promise((r) => setTimeout(r, 400));
      return {before, after: df.selection.trueCount};
    });
    expect(result.before).toBe(0);
    expect(result.after).toBeGreaterThan(0);
  });

  await softStep('Scenario 2 Step 8 — click empty space clears on the reordered chart (round-trip to zero)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const beforeClear = df.selection.trueCount;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      const ex = r0.x + r0.width * 0.5, ey = r0.y + r0.height * 0.02;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('click', mk(ex, ey)));
      await new Promise((r) => setTimeout(r, 400));
      return {beforeClear, afterClear: df.selection.trueCount};
    });
    expect(result.beforeClear).toBeGreaterThan(0);
    expect(result.afterClear).toBe(0);
  });

  v.finishSpec();
});
