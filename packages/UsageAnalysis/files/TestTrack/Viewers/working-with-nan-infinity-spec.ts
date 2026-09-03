/* ---
realizes: [viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const spgiInfinityPath = 'System:AppData/ApiTests/datasets/SPGI_v2_infinity.csv';

async function saveCurrentLayout(page: Page): Promise<string> {
  return page.evaluate(async () => {
    const layout = grok.shell.tv.saveLayout();
    await grok.dapi.layouts.save(layout);
    return layout.id;
  });
}

test('Working with NaN and Infinity values in viewers', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: spgiInfinityPath, semTypeTimeoutMs: 3000});

  const layoutIds: string[] = [];
  try {
    await softStep('1. Open SPGI_v2_infinity.csv — dataset loads with infinity values', async () => {
      const result = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const infCols: string[] = [];
        for (let i = 0; i < df.columns.length; i++) {
          const col = df.columns.byIndex(i);
          if (col.type === 'double' || col.type === 'int') {
            for (let r = 0; r < Math.min(col.length, 100); r++) {
              const val = col.get(r);
              if (val !== null && !isNaN(val) && !isFinite(val)) { infCols.push(col.name); break; }
            }
          }
        }
        return {rows: df.rowCount, cols: df.columns.length, infCols};
      });
      expect(result.rows).toBe(3624);
      expect(result.cols).toBe(88);
      expect(result.infCols.length).toBeGreaterThan(0);
    });

    await softStep('2. Add Scatter plot with Chemical Space X/Y (infinity columns)', async () => {
      await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot', 8000, 'Scatter plot');
      const [props] = await v.setViewerProps(page, 'Scatter plot', [
        {set: {xColumnName: 'Chemical Space X', yColumnName: 'Chemical Space Y'}, read: ['xColumnName', 'yColumnName']},
      ]);
      expect(props.xColumnName).toBe('Chemical Space X');
      expect(props.yColumnName).toBe('Chemical Space Y');
    });

    await softStep('3. Add Histogram with Chemical Space X (infinity column)', async () => {
      await v.addViewerByIcon(page, 'histogram', 'Histogram', 8000, 'Histogram');
      const [col] = await v.setViewerProps(page, 'Histogram', [
        {set: {valueColumnName: 'Chemical Space X'}, read: 'valueColumnName'},
      ]);
      expect(col).toBe('Chemical Space X');
    });

    await softStep('4. Save layout for SPGI dataset', async () => {
      const id = await saveCurrentLayout(page);
      if (id) layoutIds.push(id);
      expect(id).toBeTruthy();
    });

    await softStep('5. Run console script — demog with NaN height[0] and Infinity weight[0]', async () => {
      const rows = await page.evaluate(() => {
        const t = grok.data.demo.demog();
        const view = grok.shell.addTableView(t);
        view.scatterPlot({
          x: 'height', y: 'weight',
          size: 'age',
          color: 'race',
        });
        t.col('height').set(0, NaN);
        t.col('weight').set(0, Infinity);
        return t.rowCount;
      });
      await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 15000});
      const [sp] = await v.setViewerProps(page, 'Scatter plot', [{
        set: {showRegressionLine: true, markerType: 'square'}, wait: 1500,
        read: ['xColumnName', 'yColumnName', 'sizeColumnName', 'colorColumnName', 'showRegressionLine', 'markerType'],
      }]);
      expect(rows).toBe(10000);
      expect(sp.xColumnName).toBe('height');
      expect(sp.yColumnName).toBe('weight');
      expect(sp.sizeColumnName).toBe('age');
      expect(sp.colorColumnName).toBe('race');
      expect(sp.showRegressionLine).toBe(true);
      expect(sp.markerType).toBe('square');
    });

    await softStep('6. Scatter plot renders correctly with NaN/Infinity rows excluded', async () => {
      const sp = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const vw = tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return vw ? {x: vw.props.xColumnName, y: vw.props.yColumnName, regression: vw.props.showRegressionLine} : null;
      });
      expect(sp).not.toBeNull();
      expect(sp!.x).toBe('height');
      expect(sp!.y).toBe('weight');
      expect(sp!.regression).toBe(true);

      const canvas = page.locator('[name="viewer-Scatter-plot"] canvas').first();
      await expect(canvas).toBeVisible();
    });

    await softStep('7. Add Histogram on height column (NaN column)', async () => {
      const histCountBefore = await page.evaluate(() =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram').length,
      );
      await page.evaluate(() => {
        const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
        icon?.click();
      });

      await page.waitForFunction((before: number) =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram').length > before,
      histCountBefore, {timeout: 8000},
      );
      const col = await page.evaluate(() => {
        const hists = grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram') as any[];
        const hist = hists[hists.length - 1];
        hist.props.valueColumnName = 'height';
        return hist.props.valueColumnName;
      });
      expect(col).toBe('height');
    });

    await softStep('8. Add Histogram on weight column (Infinity column)', async () => {
      const histCountBefore = await page.evaluate(() =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram').length,
      );
      await page.evaluate(() => {
        const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
        icon?.click();
      });
      await page.waitForFunction((before: number) =>
        grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram').length > before,
      histCountBefore, {timeout: 8000},
      );
      const col = await page.evaluate(() => {
        const hists = grok.shell.tv.viewers.filter((x: any) => x.type === 'Histogram') as any[];
        const hist = hists[hists.length - 1];
        hist.props.valueColumnName = 'weight';
        return hist.props.valueColumnName;
      });
      expect(col).toBe('weight');
    });

    await softStep('9. Save layout for demog dataset', async () => {
      const id = await saveCurrentLayout(page);
      if (id) layoutIds.push(id);
      expect(id).toBeTruthy();
    });
  } finally {
    // dapi.save resolves before the entity always reads back: __findSaved retries the find
    await page.evaluate(async (ids: string[]) => {
      const w = window as any;
      for (const id of ids) {
        const layout = await w.__findSaved(() => grok.dapi.layouts.find(id), 2000);
        if (layout) await grok.dapi.layouts.delete(layout);
      }
    }, layoutIds);
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
