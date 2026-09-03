/* ---
realizes: [grid.cp.appearance-summary-persist]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;

// The colour-coding and gear-panel steps of the appearance scenario. The summary columns are
// package cell renderers (PowerGrid) and the layout/project round-trips are server state, so
// Steps 12-18 run on the server lane in grid-server-spec.ts.
test.use(specTestOptions);

test('Grid — Appearance, Summary Columns, and Persistence', async ({page}) => {
  test.setTimeout(240_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

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
      return {nullHeightRow, minAgeRow, maxAgeRow, minV, maxV, defaultRowHeight: grid.props.rowHeight, defaultBg: grid.cell('AGE', 0).color};
    });
    expect(setup.nullHeightRow).toBeGreaterThanOrEqual(0);
    expect(setup.minAgeRow).not.toBe(setup.maxAgeRow);

    await softStep('Step 4 — Linear colour coding on AGE: min/max cells differ, both differ from background', async () => {
      const c = await g.headerCenter(page, 'AGE');
      expect(await g.clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Linear')).toBe(true);
      const r = await v.pollValue(() => page.evaluate((s) => {
        const df = grok.shell.tv.dataFrame;
        const grid = grok.shell.tv.grid;
        return {
          ccType: df.col('AGE').getTag('.color-coding-type'),
          minColor: grid.cell('AGE', s.minAgeRow).color,
          maxColor: grid.cell('AGE', s.maxAgeRow).color,
          bgMinRow: grid.cell('DEMOG', s.minAgeRow).color,
          bgMaxRow: grid.cell('DEMOG', s.maxAgeRow).color,
        };
      }, setup), (x) => x.ccType === 'Linear' && x.minColor !== x.maxColor && x.minColor !== x.bgMinRow && x.maxColor !== x.bgMaxRow, 1500, 50);
      expect(r.ccType).toBe('Linear');
      expect(r.minColor).not.toBe(r.maxColor);
      expect(r.minColor).not.toBe(r.bgMinRow);
      expect(r.maxColor).not.toBe(r.bgMaxRow);
    });

    const condRanges = {'<160': '#0000FF', '>180': '#FF0000'};

    await softStep('Step 5 — Conditional colour coding on HEIGHT: each in-range cell resolves its configured colour', async () => {
      const c = await g.headerCenter(page, 'HEIGHT');
      expect(await g.clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Conditional')).toBe(true);
      await v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.col('HEIGHT').getTag('.color-coding-type')),
        (t) => t === 'Conditional', 1500, 50);
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
      const c = await g.headerCenter(page, 'SEX');
      expect(await g.clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Categorical')).toBe(true);
      const r = await v.pollValue(() => page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const grid = grok.shell.tv.grid;
        const sc = df.col('SEX');
        let mRow = -1; let fRow = -1;
        for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
          const val = sc.get(i);
          if (val === 'M' && mRow < 0) mRow = i;
          if (val === 'F' && fRow < 0) fRow = i;
        }
        return {ccType: sc.getTag('.color-coding-type'), mColor: grid.cell('SEX', mRow).color, fColor: grid.cell('SEX', fRow).color};
      }), (x) => x.ccType === 'Categorical' && x.mColor !== x.fColor, 1500, 50);
      expect(r.ccType).toBe('Categorical');
      expect(r.mColor).not.toBe(r.fColor);
    });

    await softStep('Step 7 — Linked colour coding on WEIGHT (source SEX): WEIGHT cell colour equals SEX cell colour', async () => {
      const c = await g.headerCenter(page, 'WEIGHT');
      expect(await g.clickMenuLeaf(page, c, ['div-Color-Coding'], 'div-Color-Coding---Edit...')).toBe(true);
      const dlg = '.d4-dialog[name="dialog-Color-coding--WEIGHT"]';
      await page.locator(dlg).waitFor({timeout: 8000});

      await page.evaluate((sel) => {
        const typeSel = document.querySelector(`${sel} [name="input-Type"]`) as HTMLSelectElement;
        const setter = Object.getOwnPropertyDescriptor(window.HTMLSelectElement.prototype, 'value')!.set!;
        setter.call(typeSel, 'Linked');
        typeSel.dispatchEvent(new Event('change', {bubbles: true}));
        typeSel.dispatchEvent(new Event('input', {bubbles: true}));
      }, dlg);
      await page.locator(`${dlg} [name="input-Source-column"]`).waitFor({timeout: 5000});
      expect(await g.pickColumnInCombo(page, `${dlg} [name="input-Source-column"]`, 'SEX')).toBe(true);
      await v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.col('WEIGHT').getTag('.%color-coding-linked-column-name')),
        (t) => t === 'SEX', 1500, 50);

      await page.locator(`${dlg} [name="button-CLOSE"]`).first().click({timeout: 5000}).catch(() => {});
      await page.locator(dlg).waitFor({state: 'detached', timeout: 5000}).catch(() => {});
      const r = await v.pollValue(() => page.evaluate(() => {
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
      }), (x) => x.ccType === 'Linked' && x.mWeight === x.mSex && x.fWeight === x.fSex && x.mWeight !== x.fWeight, 1500, 50);
      expect(r.ccType).toBe('Linked');
      expect(r.src).toBe('SEX');
      expect(r.mWeight).toBe(r.mSex);
      expect(r.fWeight).toBe(r.fSex);
      expect(r.mWeight).not.toBe(r.fWeight);
    });

    await softStep('Step 9 — Row height via the gear panel: cell bounds height reflects the new value', async () => {
      const readHeight = () => page.evaluate(() => grok.shell.tv.grid.cell('AGE', 0).bounds.height);
      const before = await readHeight();
      expect(await g.openGridSettings(page)).toBe(true);
      await page.evaluate(() => {
        const rh = document.querySelector('[name="prop-row-height"] input.property-grid-slider-textbox') as HTMLInputElement;
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(rh, '48');
        rh.dispatchEvent(new Event('input', {bubbles: true}));
        rh.dispatchEvent(new Event('change', {bubbles: true}));
        rh.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter'}));
      });
      const after = await v.pollValue(readHeight, (h) => h !== before, 1000, 50);
      expect(after).not.toBe(before);
      expect(after).toBe(48);
    });

    await softStep('Step 10 — Missing value colour via the gear panel: a null cell resolves the configured colour', async () => {
      const defaultNullColor = await page.evaluate((s) => grok.shell.tv.grid.cell('HEIGHT', s.nullHeightRow).color, setup);
      expect(await g.openGridSettings(page)).toBe(true);
      await page.locator('[name="prop-missing-value-color"]').waitFor({state: 'attached', timeout: 8000});

      await page.evaluate(() => {
        const view = document.querySelector('[name="prop-view-missing-value-color"]') as HTMLElement | null;
        if (!view) return;
        const r = view.getBoundingClientRect();
        const at = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
        for (const type of ['mouseover', 'mousedown', 'mouseup', 'click']) view.dispatchEvent(new MouseEvent(type, at));
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
      const r = await v.pollValue(() => page.evaluate((s) => {
        const grid = grok.shell.tv.grid;
        return {nullColor: grid.cell('HEIGHT', s.nullHeightRow).color, configured: grid.props.missingValueColor};
      }, setup), (x) => (x.nullColor >>> 0) === 0xffffaaaa, 1000, 50);
      expect(r.nullColor).toBe(r.configured);
      expect(r.nullColor >>> 0).toBe(0xffffaaaa);
      expect(r.nullColor).not.toBe(defaultNullColor);
    });
  } finally {
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});
