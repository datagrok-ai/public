/* ---
realizes: []
--- */

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
  // grok.shell.warnings is not exposed to JS here, hence the page/console baseline.
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
  // On top of that, the Show All Lines OFF + grid selection path carries a real
  // canvas signal: fewer painted lines = less ink, so the block at the end
  // measures settle-gated pixel counts across the toggle in both directions.
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

    // Canvas-delta: with Show All Lines off, only the selected lines are
    // painted. Each measurement is settle-gated (two consecutive counts must
    // agree) so a delta is the toggle's effect, not a render tail.
    const settledPx = async () => {
      let prev = (await v.countCanvasPixels(page, 'PC Plot')).total;
      let cur = prev;
      for (let i = 0; i < 5; i++) {
        await page.waitForTimeout(300);
        cur = (await v.countCanvasPixels(page, 'PC Plot')).total;
        if (Math.abs(cur - prev) < 200) break;
        prev = cur;
      }
      return cur;
    };
    const setState = (allLines: boolean, selectFirst: number) => page.evaluate(async (s) => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const df = grok.shell.tv.dataFrame;
      if (s.selectFirst === 0)
        df.selection.setAll(false);
      else
        for (let i = 0; i < s.selectFirst; i++) df.selection.set(i, true);
      pc.props.showAllLines = s.allLines;
      await new Promise(r => setTimeout(r, 400));
    }, {allLines, selectFirst});

    await setState(true, 0);
    const allPx = await settledPx();
    await setState(false, 0);
    const hiddenPx = await settledPx();
    await setState(false, 40);
    const selectedPx = await settledPx();
    await setState(true, 40);
    const restoredPx = await settledPx();
    // Round-trip: clear the selection, showAllLines is back at its default.
    await setState(true, 0);

    // Keep the measured ink values visible on green runs so the fixed
    // thresholds below can be audited against live numbers.
    console.log(`Selection & line display px: allPx=${allPx} hiddenPx=${hiddenPx} selectedPx=${selectedPx} restoredPx=${restoredPx}`);

    // Precheck on the empty-selection baseline: a valid measurement (>= 0, no
    // canvas fault) that sits far below the all-lines ink — hiding actually
    // removed the polylines, so the deltas below measure the toggle itself.
    expect(hiddenPx).toBeGreaterThanOrEqual(0);
    expect(allPx - hiddenPx).toBeGreaterThan(2000);
    // Selecting rows paints ONLY those lines: ink rises off the hidden floor
    // yet stays well below the all-lines total.
    expect(selectedPx - hiddenPx).toBeGreaterThan(500);
    expect(allPx - selectedPx).toBeGreaterThan(1000);
    // Re-enabling Show All Lines paints all lines again.
    expect(restoredPx - selectedPx).toBeGreaterThan(1000);
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

  // The axis-slider DOM elements persist regardless of the Show Filters state
  // (the range-handle visuals are canvas-drawn), so the assertable signal is the
  // `showFilters` prop the context-menu Filter > Show Filters item flips.
  await softStep('Show Filters from the context menu', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();
      const clickFilterSub = async (child: string) => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
        await new Promise(r => setTimeout(r, 500));
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
      return {showBefore, showToggled, showRestored};
    });
    // Show Filters menu item round-trips the showFilters state.
    expect(result.showToggled).toBe(!result.showBefore);
    expect(result.showRestored).toBe(result.showBefore);
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
