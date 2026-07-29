/* ---
realizes: [scatterplot.cp.axes-and-encode, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// NullError / "Stack trace" / cloned-iframe messages are benign ONLY during the
// ribbon project save; elsewhere that same Dart class is the regression signal.
let inProjectSaveWindow = false;
const isBenignError = (text: string) => {
  if (isAmbientError(text)) return true;
  if (inProjectSaveWindow)
    return /Unable to find element in cloned iframe/.test(text) ||
      /Stack trace [A-Za-z]+/.test(text) ||
      /NullError: method not found: '\w+' on null/.test(text);
  return false;
};

/** Viewer root is disambiguated: the Formula Lines dialog embeds its own preview
 * plot under the same name. */
const canvasRect = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** Poll `cond` until it holds. The caller keeps its own assert, so a poll that
 * runs out of time leaves the observed state to be graded as it is. */
async function waitUntil(cond: () => Promise<boolean>, timeoutMs: number, stepMs = 100): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  for (;;) {
    if (await cond()) return true;
    if (Date.now() >= deadline) return false;
    await new Promise((r) => setTimeout(r, stepMs));
  }
}

const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  return sp?.props[n] ?? null;
}, name);

const propIs = (page: Page, propName: string, value: unknown, timeoutMs: number) =>
  waitUntil(async () => await readProp(page, propName) === value, timeoutMs);

/** Assign a column through the on-viewer selector. The shared helper verifies
 * the property itself, so the short commit settle is retried through a poll on
 * that property before the pick is repeated at the helper's own pace. */
async function pickOnViewer(page: Page, role: string, column: string): Promise<void> {
  try {
    await v.pickColumnViaSelectorTrusted(page, {role, columnName: column, commitSettleMs: 250});
    return;
  } catch (_) { /* fall through to the poll, then to a full retry */ }
  if (await propIs(page, `${role}ColumnName`, column, 2500)) return;
  await v.pickColumnViaSelectorTrusted(page, {role, columnName: column});
}

async function waitBackdrop(page: Page, timeout = 5000): Promise<boolean> {
  return await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout}).then(() => true).catch(() => false);
}

/** The popup grid is canvas-rendered, so the name goes in with real keys; the
 * first key is separated because the popup focuses its search input a tick
 * later. */
async function commitColumn(page: Page, column: string): Promise<void> {
  const text = column.toLowerCase();
  await page.keyboard.press(text[0]);
  await page.waitForTimeout(150);
  if (text.length > 1) await page.keyboard.type(text.slice(1));
  await page.waitForTimeout(200);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
}

/** A row inside a collapsed category keeps its DOM node with an empty box, so
 * readiness is the editor's own rectangle, not the row's presence. The category
 * header is a toggle, hence the retry. */
async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<void> {
  const ready = () => page.evaluate((sel: string) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el || !el.offsetParent) return false;
    const b = el.getBoundingClientRect();
    return b.width > 0 && b.height > 0;
  }, editorSelector);
  for (let i = 0; i < 8; i++) {
    if (await ready()) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    if (await waitUntil(ready, 1000)) return;
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

/** The value cell turns into a select only once clicked. */
async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string, propName: string,
): Promise<void> {
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const select = page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`);
  await waitUntil(async () => await select.count() > 0 && await select.isVisible(), 2000);
  await select.selectOption(value);
  await propIs(page, propName, value, 2500);
}

async function pickInPanel(
  page: Page, rowName: string, comboName: string, column: string, propName: string,
): Promise<void> {
  const sel = page.locator(`[name="${rowName}"] [name="${comboName}"]`);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  const opened = await waitBackdrop(page);
  if (!opened) throw new Error(`${comboName} popup did not open`);
  await commitColumn(page, column);
  await propIs(page, propName, column, 2500);
}

const readConfig = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const label = (r: string) => (root?.querySelector(
    `[name="div-column-combobox-${r}"] .d4-column-selector-column`)?.textContent ?? '').trim();
  return {
    x: sp.props.xColumnName, y: sp.props.yColumnName, color: sp.props.colorColumnName,
    size: sp.props.sizeColumnName, markers: sp.props.markersColumnName,
    labels: {x: label('x'), y: label('y'), color: label('color'), size: label('size')},
  };
});

/** The grid's header context menu has no Rename entry; the New name field lives
 * in the dialog behind Column Properties. The header cell is addressed through
 * the grid's own column geometry. */
async function renameColumnViaGrid(page: Page, current: string, next: string): Promise<void> {
  const pt = await page.evaluate((name: string) => {
    const grid = grok.shell.tv.grid;
    const r = grid.root.getBoundingClientRect();
    const gc = grid.columns.byName(name);
    const headerH = grid.colHeaderH || 20;
    return {x: r.x + gc.left + gc.width / 2, y: r.y + headerH / 2};
  }, current);
  await page.mouse.click(pt.x, pt.y, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.locator('.d4-menu-popup [name="div-Column-Properties..."]').last().click();
  const nameInput = page.locator('.d4-dialog input[name="input-New-name--"]');
  await nameInput.waitFor({timeout: 8000});
  await nameInput.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(next);
  await waitUntil(async () => await nameInput.inputValue() === next, 1500);
  await page.locator('.d4-dialog [name="button-OK"]').last().click();
  await waitUntil(() => page.evaluate((n: string) =>
    document.querySelectorAll('.d4-dialog').length === 0 &&
    grok.shell.tv.dataFrame.columns.names().includes(n), next), 5000);
}

/** Submenu leaves stay hidden until their parent group is hovered. */
async function openFormulaLines(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.locator('.d4-menu-popup [name="div-Tools"]').last().hover();
  const leaf = page.locator('.d4-menu-popup [name="div-Tools---Formula-Lines..."]');
  await waitUntil(async () => await leaf.count() > 0, 3000);
  await leaf.last().click();
  await page.locator('[name="dialog-Formula-Lines"]').waitFor({timeout: 10000});
  const addNew = page.locator('[name="dialog-Formula-Lines"] [name="button-Add-new"]');
  await waitUntil(async () => await addNew.count() > 0 && await addNew.isVisible(), 5000);
}

/** Values of the dialog's formula editors — the DOM signal that its list of
 * annotations gained or lost an entry. */
const formulaEditorValues = (page: Page) => page.evaluate(() =>
  ([...document.querySelectorAll('[name="dialog-Formula-Lines"] textarea')] as HTMLTextAreaElement[])
    .map((a) => a.value).join('|'));

const formulaLineCount = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
  return JSON.parse(sp.props.formulaLines || '[]').length;
});

/** The X and Y rows share `prop-min` / `prop-max`, so an axis bound is addressed
 * through its unique view cell. */
async function setAxisBound(
  page: Page, cellName: string, value: string | null, propName: string,
): Promise<void> {
  const cell = page.locator(`[name="${cellName}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  // The cell turns into an editor that takes the focus one tick after the click.
  await waitUntil(() => page.evaluate(() =>
    document.activeElement instanceof HTMLInputElement), 2000);
  await page.keyboard.press('Control+A');
  if (value === null) await page.keyboard.press('Delete');
  else await page.keyboard.type(value);
  await page.keyboard.press('Enter');
  await propIs(page, propName, value === null ? null : Number(value), 2500);
}

const rowOpacity = (page: Page, rowName: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`) as HTMLElement | null)?.style.opacity ?? null, rowName);

test('Scatter Plot — Axes, Encodings, Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(900_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;
  const allErrors = () => [...pageErrors, ...consoleErrors];

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await waitUntil(() => page.evaluate(() => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'));
    const c = root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    return !!c && c.getBoundingClientRect().width > 0;
  }), 15000);

  await softStep('Scenario 1 — Set the axes and the encodings through the on-viewer selectors (GROK-18411)', async () => {
    const errBefore = errCount();

    // The viewer's auto-pick on attach is not contractual, so every column used
    // later is set explicitly.
    await pickOnViewer(page, 'x', 'AGE');
    await pickOnViewer(page, 'y', 'HEIGHT');

    // GROK-18411: the Color selector's label text is itself a hit target, not
    // only the dropdown triangle. It carries no text until a column is assigned,
    // so the first pick goes through the caption and the label hit target is
    // exercised right after, on the populated label.
    await pickOnViewer(page, 'color', 'RACE');
    const labelHit = await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'RACE', target: 'column', requirePopup: false,
    });
    expect(labelHit.popupOpened).toBe(true);

    await pickOnViewer(page, 'size', 'WEIGHT');

    await v.openViewerGear(page, 'Scatter plot');
    await revealPropEditor(page, '[name="prop-markers"] [name="div-column-combobox-markers"]', 'marker');
    await pickInPanel(page, 'prop-markers', 'div-column-combobox-markers', 'SEX', 'markersColumnName');

    // Switch X away and back, so the read-back is not the immediate echo of the
    // last pick.
    await pickOnViewer(page, 'x', 'WEIGHT');
    await pickOnViewer(page, 'x', 'AGE');

    const cfg = await readConfig(page);
    expect(cfg.x).toBe('AGE');
    expect(cfg.y).toBe('HEIGHT');
    expect(cfg.color).toBe('RACE');
    expect(cfg.size).toBe('WEIGHT');
    expect(cfg.markers).toBe('SEX');
    expect(cfg.labels.x).toBe('AGE');
    expect(cfg.labels.y).toBe('HEIGHT');
    expect(cfg.labels.color).toBe('RACE');
    expect(cfg.labels.size).toBe('WEIGHT');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 — One column on both axes, then renamed (GROK-19334)', async () => {
    const errBefore = errCount();
    const probeName = 'AGE_PROBE';

    await pickOnViewer(page, 'y', 'AGE');
    const both = await readConfig(page);
    expect(both.x).toBe('AGE');
    expect(both.y).toBe('AGE');

    await openFormulaLines(page);
    const noLines = await formulaEditorValues(page);
    await page.locator('[name="dialog-Formula-Lines"] [name="button-Add-new"]').click();
    await page.locator('.d4-menu-popup [name="div-Line"]').last().click();
    await waitUntil(async () => await formulaEditorValues(page) !== noLines, 4000);
    await page.locator('[name="dialog-Formula-Lines"] [name="button-OK"]').click();
    await waitUntil(async () => await formulaLineCount(page) === 1, 4000);
    expect(await formulaLineCount(page)).toBe(1);

    try {
      await renameColumnViaGrid(page, 'AGE', probeName);
      const renamed = await readConfig(page);
      expect(renamed.x).toBe(probeName);
      expect(renamed.y).toBe(probeName);
      expect(renamed.labels.x).toBe(probeName);
      expect(renamed.labels.y).toBe(probeName);
      expect(errCount()).toBe(errBefore);
    } finally {
      const names = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
      if (names.includes(probeName)) await renameColumnViaGrid(page, probeName, 'AGE');
      await openFormulaLines(page);
      const withLine = await formulaEditorValues(page);
      await page.locator('[name="dialog-Formula-Lines"] [name="button-Delete"]').click();
      await waitUntil(async () => await formulaEditorValues(page) !== withLine, 2500);
      await page.locator('[name="dialog-Formula-Lines"] [name="button-OK"]').click();
      await waitUntil(async () => await formulaLineCount(page) === 0, 4000);
      await pickOnViewer(page, 'y', 'HEIGHT');
    }
    const reverted = await readConfig(page);
    expect(reverted.x).toBe('AGE');
    expect(reverted.y).toBe('HEIGHT');
    expect(await formulaLineCount(page)).toBe(0);
  });

  await softStep('Scenario 3 — Logarithmic and inverted axis with a reversed range window (GROK-13110)', async () => {
    const errBefore = errCount();
    const errIndexBefore = allErrors().length;
    await v.openViewerGear(page, 'Scatter plot');
    await setChoiceProp(page, 'prop-x-axis-type', 'prop-view-x-axis-type', 'x', 'logarithmic', 'xAxisType');
    await page.locator('[name="prop-invert-x-axis"] input[type="checkbox"]').click();
    await propIs(page, 'invertXAxis', true, 2500);

    // GROK-13110: minimum above maximum, on an inverted logarithmic axis.
    await setAxisBound(page, 'prop-view-x-min', '60', 'xMin');
    await setAxisBound(page, 'prop-view-x-max', '20', 'xMax');
    await waitUntil(() => page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      const vp = sp?.viewport;
      return !!vp && Number.isFinite(vp.width) && vp.width > 0 &&
        Number.isFinite(vp.height) && vp.height > 0;
    }), 3000);

    const state = await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      const vp = sp.viewport;
      return {
        axisType: sp.props.xAxisType, invert: sp.props.invertXAxis,
        xMin: sp.props.xMin, xMax: sp.props.xMax,
        attached: document.body.contains(sp.root),
        vp: vp && {x: vp.x, y: vp.y, width: vp.width, height: vp.height},
      };
    });
    expect(state.axisType).toBe('logarithmic');
    expect(state.invert).toBe(true);
    expect(state.xMin).toBe(60);
    expect(state.xMax).toBe(20);
    expect(state.attached).toBe(true);
    expect(state.vp).not.toBeNull();
    for (const k of ['x', 'y', 'width', 'height'] as const)
      expect(Number.isFinite((state.vp as any)[k])).toBe(true);
    expect(state.vp!.width).toBeGreaterThan(0);
    expect(state.vp!.height).toBeGreaterThan(0);
    const raised = allErrors().slice(errIndexBefore);
    expect(raised.filter((e) => /Wrong range/i.test(e))).toEqual([]);
    expect(errCount()).toBe(errBefore);

    await setAxisBound(page, 'prop-view-x-min', null, 'xMin');
    await setAxisBound(page, 'prop-view-x-max', null, 'xMax');
    await page.locator('[name="prop-invert-x-axis"] input[type="checkbox"]').click();
    await propIs(page, 'invertXAxis', false, 2500);
    await setChoiceProp(page, 'prop-x-axis-type', 'prop-view-x-axis-type', 'x', 'linear', 'xAxisType');
    const back = await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
      return {axisType: sp.props.xAxisType, invert: sp.props.invertXAxis,
        xMin: sp.props.xMin ?? null, xMax: sp.props.xMax ?? null};
    });
    expect(back.axisType).toBe('linear');
    expect(back.invert).toBe(false);
    expect(back.xMin).toBeNull();
    expect(back.xMax).toBeNull();
  });

  await softStep('Scenario 4 — Axis type control disabled for a datetime axis (GROK-20395)', async () => {
    const errBefore = errCount();
    await revealPropEditor(page, '[name="prop-view-x-axis-type"]', 'x');

    await pickOnViewer(page, 'x', 'STARTED');
    await waitUntil(async () => await rowOpacity(page, 'prop-x-axis-type') === '0.5', 3000);
    // An inline opacity of 0.5 on the row is the only disabled signal the
    // property grid exposes; the context menu does not mirror it.
    expect(await rowOpacity(page, 'prop-x-axis-type')).toBe('0.5');
    expect(await rowOpacity(page, 'prop-x-map')).toBe('1');

    await pickOnViewer(page, 'x', 'AGE');
    await waitUntil(async () => await rowOpacity(page, 'prop-x-axis-type') === '1', 3000);
    expect(await rowOpacity(page, 'prop-x-axis-type')).toBe('1');
    expect(await rowOpacity(page, 'prop-x-map')).toBe('0.5');

    const cfg = await readConfig(page);
    expect(cfg.x).toBe('AGE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5a — Layout and project persistence at peak configuration — layout round-trip (GROK-18945)', async () => {
    const before = await readConfig(page);
    expect(before).toMatchObject({x: 'AGE', y: 'HEIGHT', color: 'RACE', size: 'WEIGHT', markers: 'SEX'});

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1500));
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        tv.addViewer('Histogram');
        await new Promise((r) => setTimeout(r, 1000));
        const sp = tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
        // Clearing a column sets the prop to the empty string, not null.
        sp.props.colorColumnName = '';
        await new Promise((r) => setTimeout(r, 800));
        const clearedColor = sp.props.colorColumnName;
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 4000));
        const restored = tv.viewers.find((vw: any) => vw.type === 'Scatter plot') as any;
        return {
          clearedColor,
          hasScatter: tv.viewers.some((vw: any) => vw.type === 'Scatter plot'),
          hasHistogram: tv.viewers.some((vw: any) => vw.type === 'Histogram'),
          x: restored?.props.xColumnName, y: restored?.props.yColumnName,
          color: restored?.props.colorColumnName, size: restored?.props.sizeColumnName,
          markers: restored?.props.markersColumnName,
        };
      }, layoutId);
      expect(result.clearedColor).toBe('');
      expect(result.hasScatter).toBe(true);
      expect(result.hasHistogram).toBe(false);
      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.color).toBe('RACE');
      expect(result.size).toBe('WEIGHT');
      expect(result.markers).toBe('SEX');
    } finally {
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Scenario 5b — Layout and project persistence at peak configuration — project save / Close All / reopen', async () => {
    const projName = 'zz-scatterplot-axes-encode-probe-' + Date.now();
    let projectId: string | null = null;
    inProjectSaveWindow = true;
    try {
      const saved = await saveProjectViaUI(page, projName);
      projectId = saved.projectId;
      expect(projectId).toBeTruthy();

      const result = await page.evaluate(async (id) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(id);
        await full.open();
        // Viewers re-materialize asynchronously — poll for the restored plot.
        let types: string[] = [];
        for (let t = 0; t < 20; t++) {
          await new Promise((r) => setTimeout(r, 1000));
          types = [];
          for (const view of grok.shell.tableViews)
            for (const vw of view.viewers) types.push(vw.type);
          if (types.includes('Scatter plot')) break;
        }
        let sp: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers)
            if (vw.type === 'Scatter plot') sp = vw;
        return {
          types,
          x: sp?.props.xColumnName, y: sp?.props.yColumnName, color: sp?.props.colorColumnName,
          size: sp?.props.sizeColumnName, markers: sp?.props.markersColumnName,
        };
      }, projectId);

      expect(result.types).toContain('Scatter plot');
      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.color).toBe('RACE');
      expect(result.size).toBe('WEIGHT');
      expect(result.markers).toBe('SEX');
    } finally {
      inProjectSaveWindow = false;
      if (projectId) await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
