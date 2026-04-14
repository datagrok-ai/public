import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Tree map viewer tests', async () => {
  test.setTimeout(600_000);

  // Reuse existing logged-in Chrome session
  const browser = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
  }
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined'
        && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi
        && grok.dapi.files;
    } catch { return false; }
  }, {timeout: 60000});
  await page.waitForFunction(() => {
    try { grok.shell.closeAll(); return true; }
    catch { return false; }
  }, {timeout: 60000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // #### Add Tree map viewer via toolbox icon
  await softStep('Add Tree map viewer', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-tree-map"]') as HTMLElement;
      icon.click();
    });
    await page.locator('[name="viewer-Tree-map"]').waitFor({timeout: 5000});
    const state = await page.evaluate(() => {
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      return {found: !!tm, splits: tm?.props.splitByColumnNames};
    });
    expect(state.found).toBe(true);
    expect(state.splits.length).toBe(1);
  });

  // #### Set first split column to RACE
  await softStep('Set first split column to RACE', async () => {
    const result = await page.evaluate(async () => {
      const selectors = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      selectors[0].value = 'RACE';
      selectors[0].dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      const sel2 = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      return {splits: tm.props.splitByColumnNames, values: Array.from(sel2).map(s => s.value)};
    });
    expect(result.splits).toEqual(['RACE']);
    expect(result.values).toEqual(['RACE', '']);
  });

  // #### Add second split level (SEX) via trailing empty selector
  await softStep('Add second split level (SEX)', async () => {
    const result = await page.evaluate(async () => {
      const selectors = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      const last = selectors[selectors.length - 1];
      last.value = 'SEX';
      last.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      const sel2 = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      return {splits: tm.props.splitByColumnNames, count: sel2.length, last: sel2[sel2.length - 1].value};
    });
    expect(result.splits).toEqual(['RACE', 'SEX']);
    expect(result.count).toBe(3);
    expect(result.last).toBe('');
  });

  // #### Truncate split chain by clearing 2nd selector
  await softStep('Truncate split chain', async () => {
    const result = await page.evaluate(async () => {
      const selectors = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      selectors[1].value = '';
      selectors[1].dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      const sel2 = document.querySelectorAll(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as NodeListOf<HTMLSelectElement>;
      return {splits: tm.props.splitByColumnNames, count: sel2.length};
    });
    expect(result.splits).toEqual(['RACE']);
    expect(result.count).toBe(2);
  });

  // #### Color column via UI popup, then aggregation via <select>
  // Actual DOM name is lowercase `div-column-combobox-color` (no trailing dash).
  // The aggregation <select> needs a two-step value change — a single `change` event is
  // swallowed when the DOM value equals the current logical value.
  await softStep('Color column and aggregation', async () => {
    // JS API fallback: in headless CDP, the column-combobox popup does not open from
    // synthesized mousedown/click events (it works in interactive sessions).
    const result = await page.evaluate(async () => {
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      tm.props.colorColumnName = 'AGE';
      await new Promise(r => setTimeout(r, 400));
      const r1 = {col: tm.props.colorColumnName, aggr1: tm.props.colorAggrType};
      // Two-step aggregation change to force propagation
      const sel = document.querySelector(
        '[name="viewer-Tree-map"] [name="div-column-combobox-color"] select.d4-column-selector-aggregation') as HTMLSelectElement;
      sel.value = 'avg';
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      sel.value = 'max';
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      return {...r1, aggr2: tm.props.colorAggrType};
    });
    expect(result.col).toBe('AGE');
    expect(result.aggr1).toBe('avg');
    expect(result.aggr2).toBe('max');
  });

  // #### Click, Shift+click, Ctrl+click selection on canvas
  await softStep('Selection with modifiers', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.dataFrame.selection.setAll(false);
      const canvas = document.querySelector(
        '[name="viewer-Tree-map"] canvas[name="canvas"]') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const fire = (type: string, x: number, y: number, mods: any = {}) =>
        canvas.dispatchEvent(new MouseEvent(type, {
          bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0,
          shiftKey: !!mods.shift, ctrlKey: !!mods.ctrl,
        }));
      const p1 = {x: r.left + 100, y: r.top + r.height / 2};
      const p3 = {x: r.left + r.width - 100, y: r.top + r.height / 2};
      fire('mousedown', p1.x, p1.y); fire('mouseup', p1.x, p1.y); fire('click', p1.x, p1.y);
      await new Promise(res => setTimeout(res, 300));
      const sel1 = grok.shell.tv.dataFrame.selection.trueCount;
      fire('click', p3.x, p3.y, {ctrl: true});
      await new Promise(res => setTimeout(res, 300));
      const sel2 = grok.shell.tv.dataFrame.selection.trueCount;
      return {sel1, sel2};
    });
    expect(result.sel1).toBeGreaterThan(0);
    // After ctrl-click the third rectangle, selection should have changed
    expect(result.sel2).not.toBe(result.sel1);
  });

  // #### Hover tooltip
  await softStep('Hover tooltip shows leaf info', async () => {
    const result = await page.evaluate(async () => {
      const canvas = document.querySelector(
        '[name="viewer-Tree-map"] canvas[name="canvas"]') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, cancelable: true,
        clientX: r.left + r.width / 2, clientY: r.top + r.height / 2,
      }));
      await new Promise(res => setTimeout(res, 500));
      const tt = document.querySelector('.d4-tooltip, [class*="tooltip"]') as HTMLElement;
      return {found: !!tt, text: tt?.textContent ?? ''};
    });
    expect(result.found).toBe(true);
    expect(result.text.length).toBeGreaterThan(0);
  });

  // #### Property pane modifications — JS API fallback
  // Title-bar gear icon is absent on Tree map even with body.selenium; use viewer.props.*.
  await softStep('Property pane modifications', async () => {
    const result = await page.evaluate(async () => {
      const tm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tree map') as any;
      const r: any = {};

      tm.props.sizeColumnName = 'AGE';
      tm.props.sizeAggrType = 'avg';
      await new Promise(res => setTimeout(res, 300));
      r.size = {col: tm.props.sizeColumnName, aggr: tm.props.sizeAggrType};

      tm.props.outerMarginLeft = 10;
      tm.props.outerMarginRight = 10;
      tm.props.outerMarginTop = 10;
      tm.props.outerMarginBottom = 10;
      await new Promise(res => setTimeout(res, 300));
      r.margins = {l: tm.props.outerMarginLeft, r: tm.props.outerMarginRight,
                   t: tm.props.outerMarginTop, b: tm.props.outerMarginBottom};

      tm.props.showColumnSelectionPanel = false;
      await new Promise(res => setTimeout(res, 300));
      const splitSel = document.querySelector(
        '[name="viewer-Tree-map"] select.d4-column-selector-tree-map') as HTMLSelectElement;
      r.panelHidden = !splitSel || splitSel.offsetParent === null;

      tm.props.showColumnSelectionPanel = true;
      await new Promise(res => setTimeout(res, 300));

      tm.props.rowSource = 'Filtered';
      r.rowSource = tm.props.rowSource;

      tm.props.filter = '${AGE} > 40';
      await new Promise(res => setTimeout(res, 500));
      r.filter = tm.props.filter;

      // Reset
      tm.props.filter = '';
      tm.props.outerMarginLeft = 0; tm.props.outerMarginRight = 0;
      tm.props.outerMarginTop = 0; tm.props.outerMarginBottom = 0;
      return r;
    });
    expect(result.size).toEqual({col: 'AGE', aggr: 'avg'});
    expect(result.margins).toEqual({l: 10, r: 10, t: 10, b: 10});
    expect(result.panelHidden).toBe(true);
    expect(result.rowSource).toBe('Filtered');
    expect(result.filter).toBe('${AGE} > 40');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
