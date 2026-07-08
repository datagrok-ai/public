import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Word Cloud tests', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Phase 2: Open dataset
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // Phase 3: Add Word Cloud viewer by clicking toolbox icon
  await v.addViewerByIcon(page, 'Word-cloud', 'Word-cloud', 10000);
  await page.waitForTimeout(1000);

  // Open viewer settings (gear icon on panel-base title bar)
  await page.evaluate(() => {
    const wcEl = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
    const panelBase = wcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  // #### Column assignment
  await softStep('Column assignment', async () => {
    const result = await page.evaluate(async () => {
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      const df = grok.shell.tv.dataFrame;
      const r: any[] = [];

      wc.props.columnColumnName = 'RACE';
      await new Promise(res => setTimeout(res, 600));
      r.push({col: wc.props.columnColumnName, error: !!document.querySelector('[name="viewer-Word-cloud"] .d4-viewer-error')});

      wc.props.columnColumnName = 'DIS_POP';
      await new Promise(res => setTimeout(res, 600));
      r.push({col: wc.props.columnColumnName, error: !!document.querySelector('[name="viewer-Word-cloud"] .d4-viewer-error')});

      // SITE is not in demog.csv — skip but record available columns
      const hasSite = df.columns.names().includes('SITE');
      r.push({siteAvailable: hasSite});

      wc.props.columnColumnName = 'SEX';
      await new Promise(res => setTimeout(res, 600));
      const sexCats = df.getCol('SEX').categories.length;
      r.push({col: wc.props.columnColumnName, sexCats});

      // Restore to RACE
      wc.props.columnColumnName = 'RACE';
      await new Promise(res => setTimeout(res, 400));

      return r;
    });
    expect(result[0].col).toBe('RACE');
    expect(result[0].error).toBe(false);
    expect(result[1].col).toBe('DIS_POP');
    expect(result[1].error).toBe(false);
    // SITE column is not present in demog.csv; skip its check
    expect(result[3].col).toBe('SEX');
    expect(result[3].sexCats).toBe(2);
  });

  // #### Shape
  await softStep('Shape', async () => {
    const result = await page.evaluate(async () => {
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      const shapes = ['circle', 'diamond', 'triangle-forward', 'triangle', 'pentagon', 'star'];
      const r: string[] = [];
      for (const s of shapes) {
        wc.props.shape = s;
        await new Promise(res => setTimeout(res, 400));
        r.push(wc.props.shape);
      }
      return r;
    });
    expect(result).toEqual(['circle', 'diamond', 'triangle-forward', 'triangle', 'pentagon', 'star']);
  });

  // #### Font size range
  await softStep('Font size range', async () => {
    const result = await v.setViewerProps(page, 'Word cloud', [
      {set: {minTextSize: 10, maxTextSize: 60}, wait: 400, read: ['minTextSize', 'maxTextSize']},
      {set: {minTextSize: 20, maxTextSize: 20}, wait: 400, read: ['minTextSize', 'maxTextSize']},
      {set: {minTextSize: 12, maxTextSize: 48}, wait: 400, read: ['minTextSize', 'maxTextSize']},
    ]);
    expect(result[0]).toEqual({minTextSize: 10, maxTextSize: 60});
    expect(result[1]).toEqual({minTextSize: 20, maxTextSize: 20});
    expect(result[2]).toEqual({minTextSize: 12, maxTextSize: 48});
  });

  // #### Text rotation
  await softStep('Text rotation', async () => {
    const result = await v.setViewerProps(page, 'Word cloud', [
      {set: {minRotationDegree: 0, maxRotationDegree: 0}, wait: 400, read: ['minRotationDegree', 'maxRotationDegree']},
      {set: {minRotationDegree: -90, maxRotationDegree: 90}, wait: 400, read: ['minRotationDegree', 'maxRotationDegree']},
      {set: {rotationStep: 45}, wait: 400, read: 'rotationStep'},
      {
        set: {minRotationDegree: -30, maxRotationDegree: 30, rotationStep: 5},
        wait: 400,
        read: ['minRotationDegree', 'maxRotationDegree', 'rotationStep'],
      },
    ]);
    expect(result[0]).toEqual({minRotationDegree: 0, maxRotationDegree: 0});
    expect(result[1]).toEqual({minRotationDegree: -90, maxRotationDegree: 90});
    expect(result[2]).toBe(45);
    expect(result[3]).toEqual({minRotationDegree: -30, maxRotationDegree: 30, rotationStep: 5});
  });

  // #### Grid size
  await softStep('Grid size', async () => {
    const result = await page.evaluate(async () => {
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      const r: number[] = [];
      for (const g of [2, 20, 8]) {
        wc.props.gridSize = g;
        await new Promise(res => setTimeout(res, 400));
        r.push(wc.props.gridSize);
      }
      return r;
    });
    expect(result).toEqual([2, 20, 8]);
  });

  // #### Draw out of bound
  await softStep('Draw out of bound', async () => {
    const result = await v.setViewerProps(page, 'Word cloud', [
      {set: {drawOutOfBound: true}, wait: 400, read: 'drawOutOfBound'},
      {set: {drawOutOfBound: false}, wait: 400, read: 'drawOutOfBound'},
      {set: {drawOutOfBound: true}, wait: 200},
    ]);
    expect(result).toEqual([true, false]);
  });

  // #### Font family
  await softStep('Font family', async () => {
    const result = await page.evaluate(async () => {
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      const r: string[] = [];
      for (const f of ['serif', 'monospace', 'sans-serif']) {
        wc.props.fontFamily = f;
        await new Promise(res => setTimeout(res, 400));
        r.push(wc.props.fontFamily);
      }
      return r;
    });
    expect(result).toEqual(['serif', 'monospace', 'sans-serif']);
  });

  // #### Bold
  await softStep('Bold', async () => {
    const result = await v.setViewerProps(page, 'Word cloud', [
      {set: {bold: true}, read: 'bold'},
      {set: {bold: false}, read: 'bold'},
      {set: {bold: true}, wait: 200},
    ]);
    expect(result).toEqual([true, false]);
  });

  // #### Tooltip on hover (canvas-based: assert a tooltip element appears on mousemove)
  await softStep('Tooltip on hover', async () => {
    const result = await page.evaluate(async () => {
      const viewerEl = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const rect = viewerEl.getBoundingClientRect();
      viewerEl.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(res => setTimeout(res, 600));
      const tip1 = !!document.querySelector('.d4-tooltip');

      viewerEl.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width * 0.3, clientY: rect.top + rect.height * 0.3,
      }));
      await new Promise(res => setTimeout(res, 600));
      const tip2 = !!document.querySelector('.d4-tooltip');

      viewerEl.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
      return {tip1, tip2};
    });
    expect(result.tip1 || result.tip2).toBe(true);
  });

  // #### Word click and selection — omitted (echarts-wordcloud canvas, no DOM-dispatchable events)

  // #### Filter interaction
  await softStep('Filter interaction', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      wc.props.columnColumnName = 'RACE';
      await new Promise(res => setTimeout(res, 300));

      grok.shell.tv.getFiltersGroup({createDefaultFilters: true});
      await new Promise(res => setTimeout(res, 1000));

      const total = df.rowCount;
      const sexCol = df.getCol('SEX');
      const maleBits = df.filter.clone();
      maleBits.init((i: number) => sexCol.get(i) === 'M');
      df.filter.copyFrom(maleBits);
      await new Promise(res => setTimeout(res, 1200));
      const filtered = df.filter.trueCount;
      const hasCanvas = !!document.querySelector('[name="viewer-Word-cloud"] canvas');

      df.filter.setAll(true);
      await new Promise(res => setTimeout(res, 800));
      const restored = df.filter.trueCount;

      return {total, filtered, restored, hasCanvas};
    });
    expect(result.filtered).toBeLessThan(result.total);
    expect(result.filtered).toBeGreaterThan(0);
    expect(result.hasCanvas).toBe(true);
    expect(result.restored).toBe(result.total);
  });

  // #### Viewer resize
  await softStep('Viewer resize', async () => {
    const result = await page.evaluate(async () => {
      const viewerEl = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const before = viewerEl.getBoundingClientRect();

      const panelBase = viewerEl.closest('.panel-base') as HTMLElement;
      const maxIcon = panelBase.querySelector('[name="icon-expand-arrows"]') as HTMLElement;
      maxIcon.click();
      await new Promise(res => setTimeout(res, 800));
      const maxed = document.querySelector('[name="viewer-Word-cloud"]')!.getBoundingClientRect();

      const panelBase2 = document.querySelector('[name="viewer-Word-cloud"]')!.closest('.panel-base') as HTMLElement;
      const maxIcon2 = panelBase2.querySelector('[name="icon-expand-arrows"]') as HTMLElement;
      maxIcon2.click();
      await new Promise(res => setTimeout(res, 800));
      const restored = document.querySelector('[name="viewer-Word-cloud"]')!.getBoundingClientRect();

      return {before: {w: before.width, h: before.height}, maxed: {w: maxed.width, h: maxed.height}, restored: {w: restored.width, h: restored.height}};
    });
    expect(result.maxed.w).toBeGreaterThan(result.before.w);
    expect(Math.abs(result.restored.w - result.before.w)).toBeLessThan(100);
  });

  // #### Context menu
  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const viewerEl = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const rect = viewerEl.getBoundingClientRect();
      viewerEl.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(res => setTimeout(res, 600));
      const ctxMenu = document.querySelector('.d4-menu-popup');
      const items = ctxMenu ? Array.from(ctxMenu.querySelectorAll('.d4-menu-item-label')).map(el => el.textContent!.trim()) : [];
      // close the menu
      document.body.click();
      await new Promise(res => setTimeout(res, 200));
      return {
        menuShown: !!ctxMenu,
        itemCount: items.length,
        hasProperties: items.some(i => /Properties/i.test(i)),
        hasSaveOrExport: items.some(i => /Save|Export/i.test(i)),
      };
    });
    expect(result.menuShown).toBe(true);
    expect(result.hasProperties).toBe(true);
    expect(result.hasSaveOrExport).toBe(true);
  });

  // #### Error state (too many categories)
  await softStep('Error state (too many categories)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const wc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Word cloud') as any;
      const r: any[] = [];

      // USUBJID has 5850 unique values
      const hiCardCol = df.columns.names().find((n: string) =>
        df.getCol(n).type === DG.TYPE.STRING && df.getCol(n).categories.length > 500);
      wc.props.columnColumnName = hiCardCol!;
      await new Promise(res => setTimeout(res, 1000));
      const err = document.querySelector('[name="viewer-Word-cloud"] .d4-viewer-error');
      r.push({hasError: !!err, errorText: err?.textContent ?? null});

      wc.props.columnColumnName = 'RACE';
      await new Promise(res => setTimeout(res, 1000));
      const err2 = document.querySelector('[name="viewer-Word-cloud"] .d4-viewer-error');
      const canvas = document.querySelector('[name="viewer-Word-cloud"] canvas');
      r.push({errorCleared: !err2, hasCanvas: !!canvas});

      return r;
    });
    expect(result[0].hasError).toBe(true);
    expect(result[0].errorText).toContain('500 or fewer unique categories');
    expect(result[1].errorCleared).toBe(true);
    expect(result[1].hasCanvas).toBe(true);
  });

  v.finishSpec();
});
