/* ---
realizes: []
--- */

import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

test('PC Plot tests', async ({page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => {
    if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text());
  });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const viewerAlive = () => page.evaluate(() =>
    !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot')
    && !!document.querySelector('[name="viewer-PC-Plot"]'));
  const pcPlotPresent = () => page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 10000});

  await v.installEventWaits(page);

  await softStep('Axis scale via the context menu', async () => {
    const menuResult = await page.evaluate(async () => {
      const w = window as any;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();

      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await w.__poll(() => document.querySelectorAll('.d4-menu-item-label').length,
          (n: number) => n > 0, 500);
      };
      const clickSub = (parent: string, child: string) =>
        w.__menuLeaf(parent, child).then(() => true, () => false);
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

  await softStep('Selection & line display', async () => {
    const errBefore = errorCount();
    await v.setViewerProps(page, 'PC Plot', [
      {set: {showCurrentLine: false}},
      {set: {showCurrentLine: true}},
      {set: {showMouseOverLine: false}},
      {set: {showMouseOverLine: true}},
      {set: {showMouseOverRowGroup: true}},
      {set: {showAllLines: false}},
      {set: {showAllLines: true}},
      {set: {showMouseOverRowGroup: false}},
    ], 150);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);

    const menu = await page.evaluate(async () => {
      const w = window as any;
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();

      const clickSub = async (parent: string, child: string) => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        return w.__menuLeaf(parent, child).then(() => true, () => false);
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

    const settledPx = async () => {
      await v.waitForCanvasQuiet(page, 'PC Plot', {timeoutMs: 1500, optional: true});
      let prev = (await v.countCanvasPixels(page, 'PC Plot')).total;
      let cur = prev;
      for (let i = 0; i < 1; i++) {
        cur = (await v.countCanvasPixels(page, 'PC Plot')).total;
        if (Math.abs(cur - prev) < 200) break;
        prev = cur;
      }
      return cur;
    };
    const setState = async (allLines: boolean, selectFirst: number) => {
      await page.evaluate((s) => {
        const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
        const df = grok.shell.tv.dataFrame;
        if (s.selectFirst === 0)
          df.selection.setAll(false);
        else
          for (let i = 0; i < s.selectFirst; i++) df.selection.set(i, true);
        pc.props.showAllLines = s.allLines;
      }, {allLines, selectFirst});
      await v.waitForViewerRendered(page, 'PC Plot', 400);
    };

    await setState(true, 0);
    const allPx = await settledPx();
    await setState(false, 0);
    const hiddenPx = await settledPx();
    await setState(false, 40);
    const selectedPx = await settledPx();
    await setState(true, 40);
    const restoredPx = await settledPx();

    await setState(true, 0);

    console.log(`Selection & line display px: allPx=${allPx} hiddenPx=${hiddenPx} selectedPx=${selectedPx} restoredPx=${restoredPx}`);

    expect(hiddenPx).toBeGreaterThanOrEqual(0);
    expect(allPx - hiddenPx).toBeGreaterThan(2000);

    expect(selectedPx - hiddenPx).toBeGreaterThan(500);
    expect(allPx - selectedPx).toBeGreaterThan(1000);

    expect(restoredPx - selectedPx).toBeGreaterThan(1000);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Style & layout', async () => {
    const errBefore = errorCount();
    await v.setViewerProps(page, 'PC Plot', [
      {set: {lineWidth: 3}},
      {set: {currentLineWidth: 5}},
      {set: {mouseOverLineWidth: 5}},
      {set: {labelsOrientation: 'Vert'}},
      {set: {minMaxOrientation: 'Vert'}},
      {set: {horzMargin: 60}},
      {set: {autoLayout: false}},
      {set: {
        lineWidth: 0.5, currentLineWidth: 2, mouseOverLineWidth: 2,
        labelsOrientation: 'Auto', minMaxOrientation: 'Auto',
        horzMargin: 40, autoLayout: true,
      }, wait: 400},
    ], 150);
    const sliders = await v.pollValue(
      () => page.locator('[name="viewer-PC-Plot"] [name^="axis-slider-"]').count(),
      (n) => n > 0, 3000, 150);
    expect(sliders).toBeGreaterThan(0);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Show Filters from the context menu', async () => {
    await v.setViewerProps(page, 'PC Plot', [{set: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}}], 800);
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();

      const clickFilterSub = async (child: string) => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
        return w.__menuLeaf('Filter', child).then(() => true, () => false);
      };
      const showBefore = pc.props.showFilters;
      await clickFilterSub('Show Filters');
      const showToggled = pc.props.showFilters;
      await clickFilterSub('Show Filters');
      const showRestored = pc.props.showFilters;
      return {showBefore, showToggled, showRestored};
    });

    expect(result.showToggled).toBe(!result.showBefore);
    expect(result.showRestored).toBe(result.showBefore);
  });

  await softStep('Title and description', async () => {
    const errBefore = errorCount();
    const readTexts = () => page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const root = pc.root as HTMLElement;
      const panel = (root.closest('.panel-base') as HTMLElement) ?? root;
      return {
        titlebar: (panel.querySelector('.panel-titlebar-text')?.textContent ?? '').trim(),
        shown: (root.innerText || '').replace(/\s+/g, ' ').trim(),
      };
    });
    await v.setViewerProps(page, 'PC Plot', [{set: {title: 'My PC Plot'}}], 600);
    const titleShown = (await v.pollValue(readTexts,
      (t) => t.titlebar.includes('My PC Plot'), 3000, 150)).titlebar;
    await v.setViewerProps(page, 'PC Plot', [{set: {description: 'Test description'}}], 600);
    const withDescription = (await v.pollValue(readTexts,
      (t) => t.shown.includes('Test description'), 3000, 150)).shown;
    await v.setViewerProps(page, 'PC Plot', [{set: {descriptionPosition: 'Bottom'}}], 600);
    const moved = (await v.pollValue(readTexts,
      (t) => t.shown.includes('Test description'), 3000, 150)).shown;
    await v.setViewerProps(page, 'PC Plot', [{set: {title: '', description: ''}}], 600);
    const after = await v.pollValue(readTexts,
      (t) => !t.titlebar.includes('My PC Plot') && !t.shown.includes('Test description'), 3000, 150);
    expect(titleShown).toContain('My PC Plot');
    expect(withDescription).toContain('Test description');
    expect(moved).toContain('Test description');
    expect(after.titlebar).not.toContain('My PC Plot');
    expect(after.shown).not.toContain('Test description');
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Pick Up / Apply', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.tv;
      const pc1 = tv.viewers.find(v => v.type === 'PC Plot')!;
      pc1.props.columnNames = ['AGE', 'WEIGHT', 'STARTED'];
      pc1.props.logColumnsColumnNames = ['AGE'];
      pc1.props.colorColumnName = 'RACE';
      pc1.props.legendPosition = 'Left';
      pc1.props.title = 'Source Plot';
      tv.addViewer('PC Plot');
    });
    await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.viewers.filter(v => v.type === 'PC Plot').length),
      (n) => n >= 2, 500, 100);

    await page.evaluate(async () => {
      const w = window as any;
      const pcs = () => grok.shell.tv.viewers.filter(v => v.type === 'PC Plot');

      const clickSub = async (idx: number, parent: string, child: string) => {
        const canvas = pcs()[idx].root.querySelector('canvas[name="canvas"]')!;
        const rect = canvas.getBoundingClientRect();
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await w.__menuLeaf(parent, child).catch(() => {});
      };
      await clickSub(0, 'Pick Up / Apply', 'Pick Up');
      await clickSub(1, 'Pick Up / Apply', 'Apply');
    });

    const applied = await v.pollValue(() => page.evaluate(() => {
      const pcs = grok.shell.tv.viewers.filter(v => v.type === 'PC Plot');
      return {
        cols: pcs[1]?.props.columnNames?.slice(),
        color: pcs[1]?.props.colorColumnName,
        log: pcs[1]?.props.logColumnsColumnNames?.slice(),
        legend: pcs[1]?.props.legendPosition,
        title: pcs[1]?.props.title,
      };
    }), (a) => a.title === 'Source Plot', 500, 100);

    await page.evaluate(() => {
      grok.shell.tv.viewers.filter(v => v.type === 'PC Plot')[0]
        .props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
    });
    await v.waitForViewerRendered(page, 'PC Plot', 600);
    const cols = await page.evaluate(() => {
      const pcs = grok.shell.tv.viewers.filter(v => v.type === 'PC Plot');
      return {pc1: pcs[0]?.props.columnNames?.slice(), pc2: pcs[1]?.props.columnNames?.slice()};
    });

    await page.evaluate(() => grok.shell.tv.dataFrame.filter.setAll(true));
    await v.waitForViewerRendered(page, 'PC Plot', 400);
    const fullBefore = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.evaluate(() => {
      const pc2El = document.querySelectorAll('[name="viewer-PC-Plot"]')[1] as HTMLElement;
      const vr = pc2El.getBoundingClientRect();
      pc2El.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
    });
    await v.pollValue(() => page.evaluate(() =>
      !!document.querySelectorAll('[name="viewer-PC-Plot"]')[1]?.querySelector('[name="axis-slider-AGE"]')),
    (present) => present, 400, 100);
    const draggedPc2 = await page.evaluate(async () => {
      const w = window as any;
      const pc2El = document.querySelectorAll('[name="viewer-PC-Plot"]')[1] as HTMLElement;
      const svg = pc2El.querySelector('[name="axis-slider-AGE"]');
      if (!svg) return false;
      const maxHandle = svg.querySelector('[name="max-handle"]')!;
      const hr = maxHandle.getBoundingClientRect();
      const cx = hr.x + hr.width / 2, cy = hr.y + hr.height / 2;
      const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
      await w.__drag(svg as HTMLElement, {x: cx, y: cy + 20}, {x: cx, y: cy + 200},
        {steps: 6, stepMs: 20, holdMs: 50});
      return true;
    });
    const filteredByPc2 = await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
      (n) => n < fullBefore, 500, 100);

    await page.evaluate(() => {
      const tv = grok.shell.tv;
      tv.viewers.filter(v => v.type === 'PC Plot')[1]?.close();
      const pcF = tv.viewers.find(v => v.type === 'PC Plot')!;
      tv.dataFrame.filter.setAll(true);
      pcF.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      pcF.props.logColumnsColumnNames = []; pcF.props.colorColumnName = '';
      pcF.props.legendPosition = 'Auto'; pcF.props.title = '';
    });
    await v.waitForViewerRendered(page, 'PC Plot', 300);

    expect(applied.cols).toEqual(['AGE', 'WEIGHT', 'STARTED']);
    expect(applied.color).toBe('RACE');
    expect(applied.log).toEqual(['AGE']);
    expect(applied.legend).toBe('Left');
    expect(applied.title).toBe('Source Plot');

    expect(cols.pc1).toEqual(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);
    expect(cols.pc2).toEqual(['AGE', 'WEIGHT', 'STARTED']);

    expect(draggedPc2).toBe(true);
    expect(filteredByPc2).toBeLessThan(fullBefore);
  });

  await softStep('Table switching and transformation', async () => {
    await v.closeAllAndWait(page);
    await v.openTable(page, {path: spgiPath, semTypeTimeoutMs: 3000});
    const spgi = await page.evaluate(() => ({
      rows: grok.shell.tv.dataFrame.rowCount, name: grok.shell.tv.dataFrame.name}));

    await page.evaluate(() => grok.shell.tv.addViewer('PC Plot'));
    await v.pollValue(pcPlotPresent, (present) => present, 500, 100);
    await v.setViewerProps(page, 'PC Plot', [{set: {table: spgi.name}}], 500);
    const tableSet = await page.evaluate(() =>
      grok.shell.tv.viewers.find(v => v.type === 'PC Plot')?.dataFrame?.name);

    const sliderAxes = () => page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
        .map((e) => e.getAttribute('name')!.replace('axis-slider-', '')));
    const axesBefore = await v.pollValue(sliderAxes, (a) => a.length > 0, 1500, 150);

    const pivot = '[{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]';
    await v.setViewerProps(page, 'PC Plot', [{set: {transformation: pivot}}], 3000);

    const axesAfter = await v.pollValue(sliderAxes, (a) => a.includes('Triazoles'), 3000, 150);
    await v.setViewerProps(page, 'PC Plot', [{set: {transformation: ''}}], 1500);
    const axesReverted = await v.pollValue(sliderAxes,
      (a) => JSON.stringify(a) === JSON.stringify(axesBefore), 1500, 150);
    await page.evaluate(() => grok.shell.tv.viewers.find(v => v.type === 'PC Plot')?.close());

    expect(spgi.rows).toBe(100);
    expect(tableSet).toBeTruthy();

    expect(axesAfter).not.toEqual(axesBefore);
    expect(axesAfter).toContain('Triazoles');
    expect(axesReverted).toEqual(axesBefore);
  });

  v.finishSpec();
});
