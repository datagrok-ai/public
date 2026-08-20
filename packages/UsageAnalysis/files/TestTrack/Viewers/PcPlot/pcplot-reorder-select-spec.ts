/* ---
realizes: [pcplot.cp.reorder-and-select, pcplot.int.area-select-cross-viewer, pcplot.int.current-row-sync]
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

  await v.installEventWaits(page);

  await page.evaluate(() => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
  });

  const readAxisNames = () => page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
      .map((e) => e.getAttribute('name')!.replace('axis-slider-', '')));

  const selectionCount = () => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
  const columnNames = () => page.evaluate(() =>
    grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!.props.columnNames.slice());

  await softStep('Setup — confirm three axes AGE, HEIGHT, WEIGHT (DOM axis-slider names)', async () => {

    const names = await v.pollValue(readAxisNames, (n) => n.length === 3, 800, 100);
    expect(names).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
  });

  await softStep('Scenario 1 Step 5 — shift+drag rectangle selects polylines (selection rises above zero)', async () => {
    const before = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cleared = new Promise((res) => {
        const sub = df.onSelectionChanged.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 200);
      });
      df.selection.setAll(false);
      await cleared;
      return df.selection.trueCount;
    });
    await page.evaluate(async () => {
      const w = window as any;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]') as HTMLElement;
      const r0 = overlay.getBoundingClientRect();
      await w.__drag(overlay,
        {x: r0.x + r0.width * 0.32, y: r0.y + r0.height * 0.40},
        {x: r0.x + r0.width * 0.45, y: r0.y + r0.height * 0.55},
        {modifiers: {shiftKey: true}});
    });
    const after = await v.pollValue(selectionCount, (n) => n > 0, 400, 50);
    expect(before).toBe(0);
    expect(after).toBeGreaterThan(0);
  });

  await softStep('Scenario 1 Step 7 — additive second shift+drag band (selection rises again, not replaced)', async () => {
    const beforeSecond = await selectionCount();
    await page.evaluate(async () => {
      const w = window as any;
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]') as HTMLElement;
      const r0 = overlay.getBoundingClientRect();
      await w.__drag(overlay,
        {x: r0.x + r0.width * 0.58, y: r0.y + r0.height * 0.42},
        {x: r0.x + r0.width * 0.72, y: r0.y + r0.height * 0.58},
        {modifiers: {shiftKey: true}});
    });
    const afterSecond = await v.pollValue(selectionCount, (n) => n > beforeSecond, 400, 50);
    expect(afterSecond).toBeGreaterThan(beforeSecond);
  });

  await softStep('Scenario 1 Step 9 — click empty space clears the selection (round-trip to zero)', async () => {
    const beforeClear = await selectionCount();
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});

      const ex = r0.x + r0.width * 0.5, ey = r0.y + r0.height * 0.02;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('click', mk(ex, ey)));
    });
    const afterClear = await v.pollValue(selectionCount, (n) => n === 0, 400, 50);
    expect(beforeClear).toBeGreaterThan(0);
    expect(afterClear).toBe(0);
  });

  await softStep('Scenario 1 Step 11 — click a polyline sets current row off -1', async () => {
    const before = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const settled = new Promise((res) => {
        const sub = df.onCurrentRowChanged.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 200);
      });
      df.currentRowIdx = -1;
      await settled;
      return df.currentRowIdx;
    });
    await page.evaluate(async () => {
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const w = window as any;
      const at = {x: r0.x + r0.width * 0.42, y: r0.y + r0.height * 0.45};
      await w.__drag(overlay, at, at, {hoverFirst: true, holdMs: 60, steps: 0, stepMs: 0});
      overlay.dispatchEvent(new MouseEvent('click',
        {bubbles: true, cancelable: true, clientX: at.x, clientY: at.y, button: 0}));
    });
    const after = await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx), (i) => i >= 0, 300, 50);
    expect(before).toBe(-1);
    expect(after).toBeGreaterThanOrEqual(0);
  });

  await softStep('Scenario 2 Step 5 — drag a column label reorders the axes (columnNames order changes, still 3)', async () => {
    const before = await columnNames();
    await page.evaluate(async () => {
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const w = window as any;
      const labelY = r0.y + 8;
      await w.__drag(overlay,
        {x: r0.x + r0.width * 0.88, y: labelY},
        {x: r0.x + r0.width * 0.08, y: labelY},
        {hoverFirst: true, holdMs: 50, steps: 10, stepMs: 30});
    });
    const after = await v.pollValue(
      columnNames, (a: string[]) => a.join('|') !== before.join('|'), 500, 50);
    expect(before).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(after.length).toBe(3);
    expect(after).not.toEqual(before);
    expect([...after].sort()).toEqual([...before].sort());
  });

  await softStep('Scenario 2 Step 7 — shift+drag on the reordered chart selects (selection rises above zero)', async () => {
    const before = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cleared = new Promise((res) => {
        const sub = df.onSelectionChanged.subscribe(() => { sub.unsubscribe(); res(undefined); });
        setTimeout(res, 200);
      });
      df.selection.setAll(false);
      await cleared;
      return df.selection.trueCount;
    });
    await page.evaluate(async () => {
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const w = window as any;
      await w.__drag(overlay,
        {x: r0.x + r0.width * 0.32, y: r0.y + r0.height * 0.40},
        {x: r0.x + r0.width * 0.45, y: r0.y + r0.height * 0.55},
        {modifiers: {shiftKey: true}});
    });
    const after = await v.pollValue(selectionCount, (n) => n > 0, 400, 50);
    expect(before).toBe(0);
    expect(after).toBeGreaterThan(0);
  });

  await softStep('Scenario 2 Step 9 — click empty space clears on the reordered chart (round-trip to zero)', async () => {
    const beforeClear = await selectionCount();
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="viewer-PC-Plot"] canvas[name="overlay"]')!;
      const r0 = overlay.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      const ex = r0.x + r0.width * 0.5, ey = r0.y + r0.height * 0.02;
      overlay.dispatchEvent(new MouseEvent('mousedown', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('mouseup', mk(ex, ey)));
      overlay.dispatchEvent(new MouseEvent('click', mk(ex, ey)));
    });
    const afterClear = await v.pollValue(selectionCount, (n) => n === 0, 400, 50);
    expect(beforeClear).toBeGreaterThan(0);
    expect(afterClear).toBe(0);
  });

  v.finishSpec();
});
