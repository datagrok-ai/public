/* ---
realizes: [grid.cp.appearance-summary-persist]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

const isBenignSaveWindowError = (text: string): boolean =>
  /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);

test.use(specTestOptions);

async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(c);
    const dataTop = grid.cell(c, 0).documentBounds.y;
    const headerY = dataTop - grid.colHeaderHeight / 2;
    return {x: rc.x + gc.left + gc.width / 2, y: headerY};
  }, col);
}

async function laidOutRect(
  page: Page, name: string,
): Promise<{x: number; y: number; w: number; h: number} | null> {
  return page.evaluate((n) => {
    const els = Array.from(document.querySelectorAll('[name="' + n + '"]')) as HTMLElement[];
    for (const el of els) {
      if (el.offsetParent === null) continue;
      const r = el.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) return {x: r.x, y: r.y, w: r.width, h: r.height};
    }
    return null;
  }, name);
}

async function clickMenuLeaf(
  page: Page, at: {x: number; y: number}, groupNames: string[], leafName: string,
): Promise<boolean> {
  for (let open = 0; open < 5; open++) {
    await page.evaluate(({x, y}) => {
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', cm));
      overlay.dispatchEvent(new MouseEvent('mouseup', cm));
      overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
    }, {x: at.x, y: at.y});
    await page.waitForTimeout(550);

    let ok = true;
    for (const g of groupNames) {
      const gr = await laidOutRect(page, g);
      if (!gr) { ok = false; break; }
      const gcx = gr.x + gr.w / 2; const gcy = gr.y + gr.h / 2;
      await page.mouse.move(gr.x + 4, gcy, {steps: 4});
      await page.mouse.move(gcx, gcy, {steps: 8});
      await page.waitForTimeout(450);
    }

    const leaf = ok ? await laidOutRect(page, leafName) : null;
    if (!leaf) {

      await page.evaluate(() => (document.body as HTMLElement).click());
      await page.waitForTimeout(250);
      continue;
    }
    await page.mouse.move(leaf.x + leaf.w / 2, leaf.y + leaf.h / 2, {steps: 4});
    await page.mouse.click(leaf.x + leaf.w / 2, leaf.y + leaf.h / 2);
    await page.waitForTimeout(800);
    return true;
  }
  return false;
}

async function openGridSettings(page: Page): Promise<boolean> {
  const rows = page.locator('[name="prop-row-height"]');
  if (await rows.count() > 0) return true;
  const gestures: Array<() => Promise<void>> = [

    async () => {
      const box = await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return null;
        const r = gear.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      });
      if (box === null) return;
      await page.mouse.move(box.x - 40, box.y + 20);
      await page.mouse.move(box.x, box.y, {steps: 6});
      await page.waitForTimeout(250);
      await page.mouse.click(box.x, box.y);
    },
    async () => {
      await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return;
        const r = gear.getBoundingClientRect();
        const cx = r.x + r.width / 2; const cy = r.y + r.height / 2;
        const host = gear.parentElement ?? gear;
        host.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: cx, clientY: cy}));
        gear.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: cx, clientY: cy}));
        for (const type of ['mousedown', 'mouseup', 'click'])
          gear.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0}));
      });
    },

    async () => { await v.openViewerGear(page, 'Grid'); },
    async () => {
      await page.locator('[name="viewer-Grid"] canvas[name="overlay"]').first().click({position: {x: 60, y: 60}, force: true});
      await page.keyboard.press('F4');
    },
  ];
  for (const gesture of gestures) {
    try {
      await gesture();
    } catch {
      continue;
    }
    try {

      await rows.first().waitFor({state: 'attached', timeout: 6000});
      return true;
    } catch {

    }
  }
  return false;
}

test('Grid — Appearance, Summary Columns, and Persistence', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const grid = grok.shell.tv.grid;
    const hc = df.col('HEIGHT');
    let nullHeightRow = -1;
    for (let i = 0; i < df.rowCount; i++) if (hc.isNone(i)) { nullHeightRow = i; break; }
    const ac = df.col('AGE');
    let minAgeRow = -1; let maxAgeRow = -1; let minV = Infinity; let maxV = -Infinity;
    for (let i = 0; i < df.rowCount; i++) {
      if (ac.isNone(i)) continue;
      const val = ac.get(i);
      if (val < minV) { minV = val; minAgeRow = i; }
      if (val > maxV) { maxV = val; maxAgeRow = i; }
    }
    return {
      nullHeightRow,
      minAgeRow, maxAgeRow, minV, maxV,
      defaultRowHeight: grid.props.rowHeight,
      defaultBg: grid.cell('AGE', 0).color, 
    };
  });
  expect(setup.nullHeightRow).toBeGreaterThanOrEqual(0); 
  expect(setup.minAgeRow).not.toBe(setup.maxAgeRow); 

  await softStep('Step 4 — Linear colour coding on AGE: min/max cells differ, both differ from background', async () => {
    const c = await headerCenter(page, 'AGE');
    const applied = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Linear');
    expect(applied).toBe(true); 
    await page.waitForTimeout(600); 
    const r = await page.evaluate((s) => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {
        ccType: df.col('AGE').getTag('.color-coding-type'),
        minColor: grid.cell('AGE', s.minAgeRow).color,
        maxColor: grid.cell('AGE', s.maxAgeRow).color,

        bgMinRow: grid.cell('DEMOG', s.minAgeRow).color,
        bgMaxRow: grid.cell('DEMOG', s.maxAgeRow).color,
      };
    }, setup);
    expect(r.ccType).toBe('Linear'); 

    expect(r.minColor).not.toBe(r.maxColor);
    expect(r.minColor).not.toBe(r.bgMinRow);
    expect(r.maxColor).not.toBe(r.bgMaxRow);
  });

  const condRanges = {'<160': '#0000FF', '>180': '#FF0000'};

  await softStep('Step 5 — Conditional colour coding on HEIGHT: each in-range cell resolves its configured colour', async () => {
    const c = await headerCenter(page, 'HEIGHT');
    const applied = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Conditional');
    expect(applied).toBe(true); 
    await page.waitForTimeout(600); 
    const r = await page.evaluate((ranges) => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const hc = df.col('HEIGHT');
      hc.meta.colors.setConditional(ranges);
      grid.invalidate();

      let lowRow = -1; let highRow = -1; let midRow = -1;
      for (let i = 0; i < df.rowCount; i++) {
        if (hc.isNone(i)) continue;
        const val = hc.get(i);
        if (lowRow < 0 && val < 160) lowRow = i;
        if (highRow < 0 && val > 180) highRow = i;
        if (midRow < 0 && val >= 160 && val <= 180) midRow = i;
        if (lowRow >= 0 && highRow >= 0 && midRow >= 0) break;
      }
      return {
        ccType: hc.getTag('.color-coding-type'),
        cond: hc.getTag('.color-coding-conditional'),
        lowColor: grid.cell('HEIGHT', lowRow).color >>> 0,
        highColor: grid.cell('HEIGHT', highRow).color >>> 0,
        midRow,
        midColor: midRow >= 0 ? grid.cell('HEIGHT', midRow).color >>> 0 : -1,
        bg: grid.cell('DEMOG', midRow >= 0 ? midRow : lowRow).color >>> 0,
      };
    }, condRanges);
    expect(r.ccType).toBe('Conditional'); 
    expect(r.cond).toBe(JSON.stringify(condRanges)); 
    expect(r.lowColor).toBe(0xff0000ff); 
    expect(r.highColor).toBe(0xffff0000); 
    expect(r.lowColor).not.toBe(r.highColor); 
    if (r.midRow >= 0) {
      expect(r.midColor).not.toBe(r.lowColor); 
      expect(r.midColor).not.toBe(r.highColor);
    }
  });

  await softStep('Step 6 — Categorical colour coding on SEX: distinct SEX values get distinct colours', async () => {
    const c = await headerCenter(page, 'SEX');
    const applied = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Categorical');
    expect(applied).toBe(true); 
    await page.waitForTimeout(600); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const sc = df.col('SEX');
      let mRow = -1; let fRow = -1;
      for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
        const val = sc.get(i);
        if (val === 'M' && mRow < 0) mRow = i;
        if (val === 'F' && fRow < 0) fRow = i;
      }
      return {
        ccType: sc.getTag('.color-coding-type'),
        mColor: grid.cell('SEX', mRow).color,
        fColor: grid.cell('SEX', fRow).color,
      };
    });
    expect(r.ccType).toBe('Categorical'); 
    expect(r.mColor).not.toBe(r.fColor); 
  });

  await softStep('Step 7 — Linked colour coding on WEIGHT (source SEX): WEIGHT cell colour equals SEX cell colour', async () => {
    const c = await headerCenter(page, 'WEIGHT');

    const opened = await clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Edit...');
    expect(opened).toBe(true); 
    await page.locator('.d4-dialog[name="dialog-Color-coding--WEIGHT"]').waitFor({timeout: 8000});

    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Color-coding--WEIGHT"]') as HTMLElement;
      const typeSel = dlg.querySelector('[name="input-Type"]') as HTMLSelectElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(typeSel, 'Linked');
      typeSel.dispatchEvent(new Event('change', {bubbles: true}));
      typeSel.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.locator('.d4-dialog[name="dialog-Color-coding--WEIGHT"] [name="input-Source-column"]')
      .waitFor({timeout: 5000});
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Color-coding--WEIGHT"]') as HTMLElement;
      const src = dlg.querySelector('[name="input-Source-column"]') as HTMLElement;
      const label = (src.querySelector('.d4-column-selector-column') as HTMLElement) ?? src;
      const r = label.getBoundingClientRect();
      label.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2}));
    });
    await page.waitForTimeout(600);
    await page.keyboard.press('s');
    await page.waitForTimeout(120);
    await page.keyboard.type('ex');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(500);

    await page.locator('.d4-dialog[name="dialog-Color-coding--WEIGHT"] [name="button-CLOSE"]')
      .first().click({timeout: 5000}).catch(() => {});
    await page.waitForTimeout(700);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const sc = df.col('SEX');

      let mRow = -1; let fRow = -1;
      for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
        const val = sc.get(i);
        if (val === 'M' && mRow < 0) mRow = i;
        if (val === 'F' && fRow < 0) fRow = i;
      }
      return {
        ccType: df.col('WEIGHT').getTag('.color-coding-type'),
        src: df.col('WEIGHT').getTag('.%color-coding-linked-column-name'),
        mWeight: grid.cell('WEIGHT', mRow).color, mSex: grid.cell('SEX', mRow).color,
        fWeight: grid.cell('WEIGHT', fRow).color, fSex: grid.cell('SEX', fRow).color,
      };
    });
    expect(r.ccType).toBe('Linked'); 
    expect(r.src).toBe('SEX'); 
    expect(r.mWeight).toBe(r.mSex); 
    expect(r.fWeight).toBe(r.fSex); 
    expect(r.mWeight).not.toBe(r.fWeight); 
  });

  await softStep('Step 9 — Row height via the gear panel: cell bounds height reflects the new value', async () => {
    const before = await page.evaluate(() => grok.shell.tv.grid.cell('AGE', 0).bounds.height);
    expect(await openGridSettings(page)).toBe(true);
    await page.evaluate(() => {
      const rh = document.querySelector('[name="prop-row-height"] input.property-grid-slider-textbox') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(rh, '48');
      rh.dispatchEvent(new Event('input', {bubbles: true}));
      rh.dispatchEvent(new Event('change', {bubbles: true}));
      rh.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter'}));
    });
    await v.waitForViewerRendered(page, 'Grid', 800);
    const after = await page.evaluate(() => grok.shell.tv.grid.cell('AGE', 0).bounds.height);
    expect(after).not.toBe(before); 
    expect(after).toBe(48); 
  });

  await softStep('Step 10 — Missing value colour via the gear panel: a null cell resolves the configured colour', async () => {
    const defaultNullColor = await page.evaluate((s) =>
      grok.shell.tv.grid.cell('HEIGHT', s.nullHeightRow).color, setup);
    expect(await openGridSettings(page)).toBe(true);
    await page.locator('[name="prop-missing-value-color"]').waitFor({state: 'attached', timeout: 8000});

    await page.evaluate(() => {
      const view = document.querySelector('[name="prop-view-missing-value-color"]') as HTMLElement | null;
      if (!view) return;
      const r = view.getBoundingClientRect();
      const at = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      for (const type of ['mouseover', 'mousedown', 'mouseup', 'click'])
        view.dispatchEvent(new MouseEvent(type, at));
    });
    await page.locator('.property-grid-item-editor-color-picker-host').waitFor({state: 'attached', timeout: 5000});
    await page.evaluate(() => {
      const host = document.querySelector('.property-grid-item-editor-color-picker-host') as HTMLElement;
      const hex = host.querySelector('input.ui-input-editor[type="text"]') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(hex, '#FFAAAA');
      hex.dispatchEvent(new Event('input', {bubbles: true}));
      hex.dispatchEvent(new Event('change', {bubbles: true}));
    });

    await page.waitForTimeout(800);
    const r = await page.evaluate((s) => {
      const grid = grok.shell.tv.grid;
      return {
        nullColor: grid.cell('HEIGHT', s.nullHeightRow).color,
        configured: grid.props.missingValueColor,
      };
    }, setup);
    expect(r.nullColor).toBe(r.configured); 
    expect(r.nullColor >>> 0).toBe(0xffffaaaa); 
    expect(r.nullColor).not.toBe(defaultNullColor); 
  });

  const summaryTypes: {leaf: string; cellType: string}[] = [
    {leaf: 'Sparklines', cellType: 'sparkline'},
    {leaf: 'Bar-Chart', cellType: 'barchart'},
    {leaf: 'Pie-Chart', cellType: 'piechart'},
    {leaf: 'Radar', cellType: 'radar'},
    {leaf: 'Smart-Form', cellType: 'smartform'},
    {leaf: 'Tags', cellType: 'tags'},
    {leaf: 'Confidence-Interval', cellType: 'confidenceinterval'},
  ];
  let colsAfterSummary = 0;

  await softStep('Step 12 — Summary columns: add all seven one-click types; column count and cellType track each add', async () => {
    const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    let count = before;
    const results: {leaf: string; added: boolean; cellType: string}[] = [];
    for (const t of summaryTypes) {
      const cellCoord = await page.evaluate(() => {
        const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
        return {x: db.x + db.width / 2, y: db.y + db.height / 2};
      });
      const added = await clickMenuLeaf(page, cellCoord, ['div-Add', 'div-Add---Summary-Columns'],
        `div-Add---Summary-Columns---${t.leaf}`);
      await page.waitForTimeout(500);
      const state = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        return {len: grid.columns.length, lastType: grid.columns.byIndex(grid.columns.length - 1).cellType};
      });
      results.push({leaf: t.leaf, added, cellType: state.lastType});
      if (added && state.len === count + 1) count = state.len;
    }
    colsAfterSummary = count;

    for (let i = 0; i < summaryTypes.length; i++) {
      expect(results[i].added).toBe(true); 
      expect(results[i].cellType).toBe(summaryTypes[i].cellType); 
    }
    expect(colsAfterSummary).toBe(before + summaryTypes.length); 
  });

  await softStep('Step 13 — Stats rows: add min and max; summary columns survive and no console error (GROK-19809)', async () => {

    const consoleErrors: string[] = [];
    const onErr = (msg: any) => {
      if (msg.type() === 'error' && !isBenignSaveWindowError(msg.text())) consoleErrors.push(msg.text());
    };
    page.on('console', onErr);
    for (const stat of ['min', 'max']) {

      const cellCoord = await page.evaluate(() => {
        const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
        return {x: db.x + db.width / 2, y: db.y + db.height / 2};
      });
      const added = await clickMenuLeaf(page, cellCoord, ['div-Add', 'div-Add---Column-Stats'],
        `div-Add---Column-Stats---${stat}`);
      expect(added).toBe(true); 
      await page.waitForTimeout(500);
    }
    page.off('console', onErr);
    const colsAfter = await page.evaluate(() => grok.shell.tv.grid.columns.length);

    expect(colsAfter).toBe(colsAfterSummary); 
    expect(consoleErrors).toEqual([]); 
  });

  await softStep('Step 15 — Persistence: save layout, add a foreign viewer, re-apply the layout (GROK-19769)', async () => {
    const r = await page.evaluate(async (s) => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const grid = tv.grid;
      const rows = (() => {
        const hc = df.col('HEIGHT'); const ac = df.col('AGE'); const scol = df.col('SEX');
        let lowRow = -1; let highRow = -1; let midRow = -1;
        for (let i = 0; i < df.rowCount; i++) {
          if (hc.isNone(i)) continue;
          const val = hc.get(i);
          if (lowRow < 0 && val < 160) lowRow = i;
          if (highRow < 0 && val > 180) highRow = i;
          if (midRow < 0 && val >= 160 && val <= 180) midRow = i;
          if (lowRow >= 0 && highRow >= 0 && midRow >= 0) break;
        }
        let minAgeRow = -1; let maxAgeRow = -1; let minV = Infinity; let maxV = -Infinity;
        for (let i = 0; i < df.rowCount; i++) {
          if (ac.isNone(i)) continue;
          const val = ac.get(i);
          if (val < minV) { minV = val; minAgeRow = i; }
          if (val > maxV) { maxV = val; maxAgeRow = i; }
        }
        let mRow = -1; let fRow = -1;
        for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
          const val = scol.get(i);
          if (val === 'M' && mRow < 0) mRow = i;
          if (val === 'F' && fRow < 0) fRow = i;
        }
        return {lowRow, highRow, midRow, minAgeRow, maxAgeRow, mRow, fRow};
      })();
      const battery = (g: any) => ({

        ageMinColor: g.cell('AGE', rows.minAgeRow).color >>> 0,
        ageMaxColor: g.cell('AGE', rows.maxAgeRow).color >>> 0,
        agePlainColor: g.cell('DEMOG', rows.minAgeRow).color >>> 0,

        heightLowColor: g.cell('HEIGHT', rows.lowRow).color >>> 0,
        heightHighColor: g.cell('HEIGHT', rows.highRow).color >>> 0,
        heightMidColor: rows.midRow >= 0 ? g.cell('HEIGHT', rows.midRow).color >>> 0 : -1,

        sexMColor: g.cell('SEX', rows.mRow).color >>> 0,
        sexFColor: g.cell('SEX', rows.fRow).color >>> 0,

        weightMEqualsSexM: g.cell('WEIGHT', rows.mRow).color === g.cell('SEX', rows.mRow).color,
        weightFEqualsSexF: g.cell('WEIGHT', rows.fRow).color === g.cell('SEX', rows.fRow).color,
      });
      const before = {
        cols: grid.columns.length,
        ageCC: df.col('AGE').getTag('.color-coding-type'),
        heightCC: df.col('HEIGHT').getTag('.color-coding-type'),
        heightCond: df.col('HEIGHT').getTag('.color-coding-conditional'),
        sexCC: df.col('SEX').getTag('.color-coding-type'),
        weightCC: df.col('WEIGHT').getTag('.color-coding-type'),
        weightSrc: df.col('WEIGHT').getTag('.%color-coding-linked-column-name'),
        ageBounds: grid.cell('AGE', 0).bounds.height,
        nullColor: grid.cell('HEIGHT', s.nullHeightRow).color,
        cellTypes: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).cellType),
        battery: battery(grid),
      };
      const layout = await grok.dapi.layouts.save(tv.saveLayout());
      await new Promise((res) => setTimeout(res, 800));
      tv.addViewer('Scatter plot');
      await new Promise((res) => setTimeout(res, 900));
      const hadScatter = tv.viewers.some((x: any) => x.type === 'Scatter plot');
      tv.loadLayout(layout);
      await new Promise((res) => setTimeout(res, 2800));
      const g2 = grok.shell.tv.grid;
      const df2 = grok.shell.tv.dataFrame;
      const after = {
        cols: g2.columns.length,
        ageCC: df2.col('AGE').getTag('.color-coding-type'),
        heightCC: df2.col('HEIGHT').getTag('.color-coding-type'),
        heightCond: df2.col('HEIGHT').getTag('.color-coding-conditional'),
        sexCC: df2.col('SEX').getTag('.color-coding-type'),
        weightCC: df2.col('WEIGHT').getTag('.color-coding-type'),
        weightSrc: df2.col('WEIGHT').getTag('.%color-coding-linked-column-name'),
        ageBounds: g2.cell('AGE', 0).bounds.height,
        nullColor: g2.cell('HEIGHT', s.nullHeightRow).color,
        cellTypes: Array.from({length: g2.columns.length}, (_: any, i: number) => g2.columns.byIndex(i).cellType),
        battery: battery(g2),
        scatterGone: !grok.shell.tv.viewers.some((x: any) => x.type === 'Scatter plot'),
      };
      await grok.dapi.layouts.delete(layout);
      return {before, after, hadScatter};
    }, setup);
    expect(r.hadScatter).toBe(true); 
    expect(r.after.scatterGone).toBe(true); 

    expect(r.after.ageCC).toBe('Linear');
    expect(r.after.battery.ageMinColor).not.toBe(r.after.battery.ageMaxColor);
    expect(r.after.battery.ageMinColor).not.toBe(r.after.battery.agePlainColor);
    expect(r.after.battery.ageMaxColor).not.toBe(r.after.battery.agePlainColor);
    expect(r.after.heightCC).toBe('Conditional');
    expect(r.after.heightCond).toBe(r.before.heightCond); 
    expect(r.after.battery.heightLowColor).toBe(0xff0000ff); 
    expect(r.after.battery.heightHighColor).toBe(0xffff0000); 

    expect(r.after.battery.heightMidColor).not.toBe(0xff0000ff);
    expect(r.after.battery.heightMidColor).not.toBe(0xffff0000);
    expect(r.after.sexCC).toBe('Categorical');

    expect(r.after.battery.sexMColor).not.toBe(r.after.battery.sexFColor);
    expect(r.after.weightCC).toBe('Linked');
    expect(r.after.weightSrc).toBe('SEX');

    expect(r.after.battery.weightMEqualsSexM).toBe(true);
    expect(r.after.battery.weightFEqualsSexF).toBe(true);
    expect(r.after.ageBounds).toBe(r.before.ageBounds); 
    expect(r.after.nullColor).toBe(r.before.nullColor); 

    expect(r.after.cols).toBe(r.before.cols); 
    expect(r.after.cellTypes).toEqual(r.before.cellTypes); 
  });

  const projectName = 'grid-appearance-summary-persist-' + Date.now();
  let savedProjectId: string | null = null;

  let peakCellHeight: number | null = null;

  await softStep('Step 16 — Persistence: save the view as a project via the ribbon Save button', async () => {
    peakCellHeight = await page.evaluate(() => grok.shell.tv.grid.cell('AGE', 0).bounds.height);

    const saved = await saveProjectViaUI(page, projectName);
    savedProjectId = saved.projectId;
    expect(savedProjectId).not.toBeNull();

    const saveBalloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim()));
    const realSaveBalloons = saveBalloons.filter((t) => !isBenignSaveWindowError(t));
    expect(realSaveBalloons).toEqual([]); 
  });

  await softStep('Step 17 — Persistence: Close All and reopen the project; the full battery holds, error delta 0', async () => {

    const consoleErrors: string[] = [];
    const onErr = (msg: any) => {
      if (msg.type() === 'error' && !isBenignSaveWindowError(msg.text())) consoleErrors.push(msg.text());
    };
    page.on('console', onErr);
    const r = await page.evaluate(async (args) => {
      grok.shell.closeAll();
      await new Promise((res) => setTimeout(res, 1500));
      const proj = await grok.dapi.projects.find(args.pid);
      await proj.open();
      await new Promise((res) => setTimeout(res, 4500));
      const grid = grok.shell.tv?.grid;
      const df = grok.shell.tv?.dataFrame;
      const loadFailureBalloons = Array.from(
        document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim())
        .filter((t) => /error loading/i.test(t));

      const battery = (() => {
        if (!grid) return null;
        const hc = df.col('HEIGHT'); const ac = df.col('AGE'); const scol = df.col('SEX');
        let lowRow = -1; let highRow = -1; let midRow = -1;
        for (let i = 0; i < df.rowCount; i++) {
          if (hc.isNone(i)) continue;
          const val = hc.get(i);
          if (lowRow < 0 && val < 160) lowRow = i;
          if (highRow < 0 && val > 180) highRow = i;
          if (midRow < 0 && val >= 160 && val <= 180) midRow = i;
          if (lowRow >= 0 && highRow >= 0 && midRow >= 0) break;
        }
        let minAgeRow = -1; let maxAgeRow = -1; let minV = Infinity; let maxV = -Infinity;
        for (let i = 0; i < df.rowCount; i++) {
          if (ac.isNone(i)) continue;
          const val = ac.get(i);
          if (val < minV) { minV = val; minAgeRow = i; }
          if (val > maxV) { maxV = val; maxAgeRow = i; }
        }
        let mRow = -1; let fRow = -1;
        for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
          const val = scol.get(i);
          if (val === 'M' && mRow < 0) mRow = i;
          if (val === 'F' && fRow < 0) fRow = i;
        }
        return {
          ageMinColor: grid.cell('AGE', minAgeRow).color >>> 0,
          ageMaxColor: grid.cell('AGE', maxAgeRow).color >>> 0,
          agePlainColor: grid.cell('DEMOG', minAgeRow).color >>> 0,
          heightLowColor: grid.cell('HEIGHT', lowRow).color >>> 0,
          heightHighColor: grid.cell('HEIGHT', highRow).color >>> 0,
          heightMidColor: midRow >= 0 ? grid.cell('HEIGHT', midRow).color >>> 0 : -1,
          sexMColor: grid.cell('SEX', mRow).color >>> 0,
          sexFColor: grid.cell('SEX', fRow).color >>> 0,
          weightMEqualsSexM: grid.cell('WEIGHT', mRow).color === grid.cell('SEX', mRow).color,
          weightFEqualsSexF: grid.cell('WEIGHT', fRow).color === grid.cell('SEX', fRow).color,
        };
      })();
      return {
        reopened: !!grid,
        loadFailureBalloons,
        ageCC: grid ? df.col('AGE').getTag('.color-coding-type') : null,
        heightCC: grid ? df.col('HEIGHT').getTag('.color-coding-type') : null,
        sexCC: grid ? df.col('SEX').getTag('.color-coding-type') : null,
        weightCC: grid ? df.col('WEIGHT').getTag('.color-coding-type') : null,
        weightSrc: grid ? df.col('WEIGHT').getTag('.%color-coding-linked-column-name') : null,
        battery,
        ageBounds: grid ? grid.cell('AGE', 0).bounds.height : -1,
        nullColor: grid ? grid.cell('HEIGHT', args.nullHeightRow).color : -1,
        cols: grid ? grid.columns.length : -1,
        cellTypes: grid ? Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).cellType) : [],
      };
    }, {pid: savedProjectId, nullHeightRow: setup.nullHeightRow});
    page.off('console', onErr);
    expect(r.reopened).toBe(true); 
    expect(r.loadFailureBalloons).toEqual([]); 
    expect(consoleErrors).toEqual([]); 
    expect(r.ageCC).toBe('Linear');
    expect(r.battery!.ageMinColor).not.toBe(r.battery!.ageMaxColor); 
    expect(r.battery!.ageMinColor).not.toBe(r.battery!.agePlainColor); 
    expect(r.battery!.ageMaxColor).not.toBe(r.battery!.agePlainColor);
    expect(r.heightCC).toBe('Conditional');
    expect(r.battery!.heightLowColor).toBe(0xff0000ff); 
    expect(r.battery!.heightHighColor).toBe(0xffff0000); 
    expect(r.battery!.heightMidColor).not.toBe(0xff0000ff); 
    expect(r.battery!.heightMidColor).not.toBe(0xffff0000);
    expect(r.sexCC).toBe('Categorical');
    expect(r.battery!.sexMColor).not.toBe(r.battery!.sexFColor); 
    expect(r.weightCC).toBe('Linked');
    expect(r.weightSrc).toBe('SEX');
    expect(r.battery!.weightMEqualsSexM).toBe(true); 
    expect(r.battery!.weightFEqualsSexF).toBe(true);
    expect(r.ageBounds).toBe(peakCellHeight); 
    expect(r.nullColor >>> 0).toBe(0xffffaaaa); 
    expect(r.cols).toBe(colsAfterSummary); 

    expect(r.cellTypes).toEqual(expect.arrayContaining(summaryTypes.map((t) => t.cellType)));
  });

  await softStep('Step 18 — Teardown: delete the probe project', async () => {
    if (savedProjectId)
      await deleteProjectWithCleanup(page, {projectId: savedProjectId});
  });

  v.finishSpec();
});
