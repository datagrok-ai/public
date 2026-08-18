/* ---
realizes: [pivottable.cp.inner-grid-look-viewers, pivottable.int.viewer-columns-per-pivot-category]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const PIVOT = '[name="viewer-Pivot-table"]';

const INNER_CANVAS = `${PIVOT} .grok-pivot-grid [name="viewer-Grid"] canvas[name="canvas"]`;

async function pivotProps(page: Page) {
  return page.evaluate(() => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    return {
      groupBy: pv.props.groupByColumnNames,
      pivot: pv.props.pivotColumnNames,
      agg: pv.props.aggregateColumnNames,
      aggTypes: pv.props.aggregateAggTypes,
    };
  });
}

async function gridLookColumns(page: Page): Promise<{
  count: number;
  cols: {columnName: string; visible: boolean; colorCodingType: string; isColorCoded: boolean; width: number}[];
}> {
  return page.evaluate(() => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    const gl = pv.getOptions(true).look.gridLook;
    const cols = (gl.columns || []).map((c: any) => ({
      columnName: c.columnName, visible: c.visible,
      colorCodingType: c.colorCodingType, isColorCoded: c.isColorCoded, width: c.width,
    }));
    return {count: cols.length, cols};
  });
}

async function gridLookColumn(page: Page, columnName: string) {
  const {cols} = await gridLookColumns(page);
  return cols.find((c) => c.columnName === columnName);
}

async function aggregateChipCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const scope = document.querySelector('[name="viewer-Pivot-table"]');
    const row = scope && [...scope.querySelectorAll('.grok-pivot-column-panel')]
      .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Aggregate"]'));
    return row ? row.querySelectorAll('.d4-tag').length : 0;
  });
}

async function pivotTitleDom(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.panel-titlebar-tabhost .panel-titlebar-text'))
      .map((n) => (n.textContent || '').trim())
      .filter(Boolean));
}

async function innerRect(page: Page) {
  await page.locator(INNER_CANVAS).first().scrollIntoViewIfNeeded().catch(() => {});
  const deadline = Date.now() + 10000;
  let box = await page.locator(INNER_CANVAS).first().boundingBox();
  while (Date.now() < deadline) {
    const vp = page.viewportSize() ?? {width: 1920, height: 1080};
    if (box && box.width > 20 && box.height > 20 &&
        box.x >= 0 && box.y >= 0 && box.x + 40 <= vp.width && box.y + 20 <= vp.height)
      return box;
    await page.waitForTimeout(300);   
    box = await page.locator(INNER_CANVAS).first().boundingBox();
  }
  if (!box) throw new Error('inner pivot grid canvas not visible');
  return box;
}

const KEY_X = 62;                            
const VALUE_X = 166;                         
const VALUE_COL = 'Critical avg(AGE)';       
const HEADER_Y = 12;                         
const rowY = (r: number) => 24 + r * 24 + 12;

function headerPoint(box: {x: number; y: number; width: number; height: number},
    vp: {width: number; height: number}, xLocal: number): {px: number; py: number} {
  const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));
  const px = clamp(clamp(box.x + xLocal, box.x + 2, box.x + box.width - 2), 1, vp.width - 1);
  const py = clamp(clamp(box.y + HEADER_Y, box.y + 2, box.y + box.height - 2), 1, vp.height - 1);
  return {px, py};
}

async function rightClickHeader(page: Page, xLocal: number) {
  const vp = page.viewportSize() ?? {width: 1920, height: 1080};
  let lastErr: unknown = null;
  for (let attempt = 0; attempt < 3; attempt++) {
    try {
      const box = await innerRect(page);   
      const {px, py} = headerPoint(box, vp, xLocal);
      await page.evaluate(({sel, cx, cy}) => {
        const cv = (document.querySelector(sel.replace('canvas[name="canvas"]', 'canvas[name="overlay"]'))
          ?? document.querySelector(sel)) as HTMLCanvasElement | null;
        if (!cv) throw new Error('inner-grid overlay canvas not found for contextmenu dispatch');
        cv.dispatchEvent(new MouseEvent('contextmenu',
          {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 2, buttons: 2}));
      }, {sel: INNER_CANVAS, cx: px, cy: py});
      await page.locator('.d4-menu-popup').first().waitFor({timeout: 4000});
      return;
    } catch (e) {
      lastErr = e;
      await page.keyboard.press('Escape').catch(() => {});
      await page.waitForTimeout(400);      
    }
  }
  throw new Error(`rightClickHeader: could not open the inner-grid header menu at xLocal=${xLocal} — ${String(lastErr)}`);
}

async function expandGroup(page: Page, groupName: string, childName: string) {
  const box = await page.locator(`.d4-menu-popup [name="${groupName}"]`).first().boundingBox();
  if (!box) throw new Error(`grid-menu group ${groupName} not found`);
  const px = box.x + box.width / 2, py = box.y + box.height / 2;
  const childReady = () => page.evaluate((name) => {
    const el = document.querySelector(`.d4-menu-popup [name="${name}"]`) as HTMLElement | null;
    if (!el) return false;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  }, childName);
  for (let i = 0; i < 25; i++) {
    await page.mouse.move(px + (i % 2), py);   
    await page.waitForTimeout(150);            
    if (await childReady()) return;
  }
  throw new Error(`grid-menu ${groupName} flyout did not reveal ${childName}`);
}

async function clickLeaf(page: Page, leafName: string) {
  const box = await page.locator(`.d4-menu-popup [name="${leafName}"]`).first().boundingBox();
  if (!box) throw new Error(`grid-menu leaf ${leafName} not laid out`);
  const lx = box.x + box.width / 2, ly = box.y + box.height / 2;
  await page.mouse.move(lx, ly);
  await page.waitForTimeout(120);   
  await page.mouse.down();
  await page.mouse.up();
  await v.waitForViewerRendered(page, 'Pivot table', 400);
}

async function clickGridLeaf(page: Page, leafName: string) {
  await expandGroup(page, 'div-Grid', leafName);
  await clickLeaf(page, leafName);
}

async function applyLinearColorCoding(page: Page, xLocal: number = VALUE_X) {
  await rightClickHeader(page, xLocal);            
  await expandGroup(page, 'div-Grid', 'div-Grid---Color-Coding');
  await expandGroup(page, 'div-Grid---Color-Coding', 'div-Grid---Color-Coding---Linear');
  await clickLeaf(page, 'div-Grid---Color-Coding---Linear');
  await page.keyboard.press('Escape');
  await v.waitForViewerRendered(page, 'Pivot table', 600);
}

const VIEWER_PICKER = `${PIVOT} .grok-pivot-column-tags-title[d4-name="Aggregate"] .d4-combo-popup`;
async function pickViewerColumn(page: Page) {
  const box = await page.locator(VIEWER_PICKER).first().boundingBox();
  if (!box) throw new Error('viewer-picker combo-popup not visible');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);   
  await page.locator(`${VIEWER_PICKER}.d4-combo-popup-expanded .d4-list-item`).first().waitFor({timeout: 5000});
  await page.locator(`${VIEWER_PICKER}.d4-combo-popup-expanded .d4-list-item [name="icon-scatter-plot"]`).first().click();
  await v.waitForViewerRendered(page, 'Pivot table', 1200);
}

test('Pivot Table — Inner grid look and viewer columns', async ({page}) => {
  test.setTimeout(360_000);

  const consoleErrors: string[] = [];
  const IGNORE = /Unable to find element in cloned iframe/i;
  page.on('console', (m) => { if (m.type() === 'error' && !IGNORE.test(m.text())) consoleErrors.push(m.text()); });
  const errCount = () => consoleErrors.length;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
  await v.waitForViewerRendered(page, 'Pivot table', 1000);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows with the auto cross-tab', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));
    const props = await pivotProps(page);
    expect(props.groupBy).toEqual(['DIS_POP']);
    expect(props.agg).toEqual(['AGE']);
    expect(props.aggTypes).toEqual(['avg']);
    expect(props.pivot).toEqual(['SEVERITY']);
    await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});
  });

  await softStep('Scenario 1 Step 4: the value column is colour-coded Linear after the header-menu action, no console error', async () => {
    const errBefore = errCount();

    const baseline = await gridLookColumn(page, VALUE_COL);
    expect(baseline?.colorCodingType).toBe('Off');
    await applyLinearColorCoding(page);

    const coded = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.colorCodingType === 'Linear', 3000, 150);
    expect(coded?.colorCodingType).toBe('Linear');
    expect(coded?.isColorCoded).toBe(true);
    expect(errCount()).toBe(errBefore);   
  });

  await softStep('Scenario 3 Step 4: Grid > Hide hides the value column in the inner grid (gridLook visible flips false)', async () => {
    const errBefore = errCount();

    const beforeCol = await gridLookColumn(page, VALUE_COL);
    expect(beforeCol?.visible).toBe(true);

    const rc = await innerRect(page);
    const vpc = page.viewportSize() ?? {width: 1920, height: 1080};
    const cell = headerPoint(rc, vpc, VALUE_X);
    await page.mouse.click(cell.px, rc.y + rowY(0));   
    await v.waitForViewerRendered(page, 'Pivot table', 300);
    await rightClickHeader(page, VALUE_X);
    await clickGridLeaf(page, 'div-Grid---Hide');
    await page.keyboard.press('Escape');
    await v.waitForViewerRendered(page, 'Pivot table', 700);

    const afterCol = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.visible === false, 3000, 150);
    expect(afterCol?.visible).toBe(false);             
  });

  await softStep('Scenario 3 Step 6: the select-all checkbox in Order or Hide Columns restores the hidden column', async () => {

    const errBefore = errCount();
    await rightClickHeader(page, KEY_X);
    await clickGridLeaf(page, 'div-Grid---Order-or-Hide-Columns...');
    await page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 8000});
    await page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"] input[type="checkbox"]').first().click();
    await v.waitForViewerRendered(page, 'Pivot table', 500);
    const close = page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"] [name="button-CLOSE"]');
    if (await close.count() > 0) await close.first().click();
    else await page.keyboard.press('Escape');
    await v.pollValue(() => page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"]').count(),
      (n) => n === 0, 600, 100);
    const restored = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.visible === true, 3000, 150);
    expect(restored?.visible).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 Step 3: with one pivot column, the viewer picker adds a viewer column, no console error', async () => {
    const errBefore = errCount();

    await v.setViewerProps(page, 'Pivot table', [{set: {pivotColumnNames: ['SEVERITY']}}], 700);
    expect((await pivotProps(page)).pivot).toEqual(['SEVERITY']);

    const chipsBefore = await aggregateChipCount(page);

    await pickViewerColumn(page);

    const chipsAfter = await v.pollValue(() => aggregateChipCount(page),
      (n) => n === chipsBefore + 1, 3000, 150);
    expect(chipsAfter).toBe(chipsBefore + 1);   
    expect(errCount()).toBe(errBefore);         
  });

  await softStep('Scenario 4 Step 6: with two pivot columns configured, the viewer picker runs error-free', async () => {

    const errBefore = errCount();
    await v.setViewerProps(page, 'Pivot table', [{set: {pivotColumnNames: ['SEVERITY', 'SEX']}}], 900);
    expect((await pivotProps(page)).pivot).toEqual(['SEVERITY', 'SEX']);
    await pickViewerColumn(page);
    expect(errCount()).toBe(errBefore);
  });

  const LAYOUT_NAME = `pivot-look-test-${Date.now()}`;
  let layoutId: string | null = null;

  try {

    const SPGI_VALUE_X = 584;
    const SPGI_VALUE_COL = 'avg(CAST Idea ID)';
    await softStep('Scenario 6 Step 9: after re-applying the saved layout, the pivot is present, titled "Pivot Overview", and the coloured column keeps its colour', async () => {

      await page.evaluate(async () => {
        grok.shell.closeAll();
        const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
        grok.shell.addTableView(df);
        await new Promise((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
          setTimeout(resolve, 3000);
        });
      });
      await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
      await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
      await v.waitForViewerRendered(page, 'Pivot table', 1000);

      await page.evaluate(() => {
        const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
        const df = grok.shell.tv.dataFrame;
        const cats = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
          .filter((c: any) => c.type === 'string');
        const nums = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
          .filter((c: any) => c.type === 'double' || c.type === 'int');
        if (cats.length >= 2) pv.props.groupByColumnNames = [cats[0].name, cats[1].name];
        pv.props.pivotColumnNames = [];
        if (nums.length >= 1) { pv.props.aggregateColumnNames = [nums[0].name]; pv.props.aggregateAggTypes = ['avg']; }
        pv.props.title = 'Pivot Overview';
      });
      await v.waitForViewerRendered(page, 'Pivot table', 900);
      await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});

      await applyLinearColorCoding(page, SPGI_VALUE_X);
      const spgiCoded = await v.pollValue(() => gridLookColumn(page, SPGI_VALUE_COL),
        (c) => c?.colorCodingType === 'Linear', 3000, 150);
      expect(spgiCoded?.colorCodingType).toBe('Linear');

      layoutId = await page.evaluate(async (name: string) => {
        const tv = grok.shell.tv;
        const layout = tv.saveLayout();
        layout.name = name;
        await grok.dapi.layouts.save(layout);

        const t0 = Date.now();
        while (Date.now() - t0 < 1000) {
          if (await grok.dapi.layouts.find(layout.id)) break;
          await new Promise((r) => setTimeout(r, 100));
        }
        return layout.id;
      }, LAYOUT_NAME);

      await v.setViewerProps(page, 'Pivot table', [{set: {title: 'Perturbed'}}], 500);
      await page.evaluate(async (id: string) => {
        const saved = await grok.dapi.layouts.find(id);

        try { grok.shell.tv.loadLayout(saved); } catch (_) {}

        const t0 = Date.now();
        while (Date.now() - t0 < 3000) {
          const pv = Array.from(grok.shell.tv?.viewers ?? []).find((x: any) => x.type === 'Pivot table') as any;
          if (pv?.props.title === 'Pivot Overview') break;
          await new Promise((r) => setTimeout(r, 100));
        }
      }, layoutId!);
      await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});
      const restored = await page.evaluate(() => {
        const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
        return {present: !!pv};
      });
      expect(restored.present).toBe(true);

      const restoredTitleDom = await v.pollValue(() => pivotTitleDom(page),
        (t) => t.includes('Pivot Overview'), 3000, 150);
      expect(restoredTitleDom).toContain('Pivot Overview');

      const restoredCol = await v.pollValue(() => gridLookColumn(page, SPGI_VALUE_COL),
        (c) => c?.colorCodingType === 'Linear', 3000, 150);
      expect(restoredCol?.colorCodingType).toBe('Linear');
    });

    await softStep('Scenario 6 Step 10: driving the ribbon Save opens the Save-project dialog with no error (project-persistence entry point)', async () => {

      const errBefore = errCount();
      await page.locator('[name="button-Save"]').first().click();   
      const saveDialog = page.locator('.d4-dialog').filter({hasText: 'Save project'});
      await saveDialog.first().waitFor({timeout: 8000});
      expect(await saveDialog.count()).toBeGreaterThan(0);   
      expect(errCount()).toBe(errBefore);                    

      await page.locator('[name="button-CANCEL"]').first().click().catch(() => page.keyboard.press('Escape'));
      await v.pollValue(() => saveDialog.count(), (n) => n === 0, 400, 100);
    });
  } finally {
    await page.evaluate(async (lid: string | null) => {
      try { if (lid) { const l = await grok.dapi.layouts.find(lid); if (l) await grok.dapi.layouts.delete(l); } } catch (_) {}
    }, layoutId);
  }

  v.finishSpec();
});
