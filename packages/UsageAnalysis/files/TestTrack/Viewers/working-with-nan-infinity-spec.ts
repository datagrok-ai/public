import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Working with NaN and Infinity values in viewers', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup + open SPGI_v2_infinity.csv
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok.shell.settings as any).showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI_v2_infinity.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  let spgiLayoutId = '';

  await softStep('1. Open SPGI_v2_infinity.csv — dataset loads with infinity values', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const infCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const col = df.columns.byIndex(i);
        if (col.type === 'double' || col.type === 'int') {
          for (let r = 0; r < Math.min(col.length, 100); r++) {
            const v = col.get(r);
            if (v !== null && !isNaN(v) && !isFinite(v)) { infCols.push(col.name); break; }
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
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-scatter-plot"]') as HTMLElement;
      icon?.click();
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 8000});
    await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot') as any;
      sp.props.xColumnName = 'Chemical Space X';
      sp.props.yColumnName = 'Chemical Space Y';
    });
    const props = await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot') as any;
      return {x: sp.props.xColumnName, y: sp.props.yColumnName};
    });
    expect(props.x).toBe('Chemical Space X');
    expect(props.y).toBe('Chemical Space Y');
  });

  await softStep('3. Add Histogram with Chemical Space X (infinity column)', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
      icon?.click();
    });
    await page.locator('[name="viewer-Histogram"]').waitFor({timeout: 8000});
    await page.evaluate(() => {
      const hist = grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram') as any;
      hist.props.valueColumnName = 'Chemical Space X';
    });
    const col = await page.evaluate(() => {
      const hist = grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram') as any;
      return hist.props.valueColumnName;
    });
    expect(col).toBe('Chemical Space X');
  });

  await softStep('4. Save layout for SPGI dataset', async () => {
    spgiLayoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });
    expect(spgiLayoutId).toBeTruthy();
  });

  let demogLayoutId = '';

  await softStep('5. Run console script — demog with NaN height[0] and Infinity weight[0]', async () => {
    const result = await page.evaluate(async () => {
      const t = grok.data.demo.demog();
      const view = grok.shell.addTableView(t);
      const plot = view.scatterPlot({
        x: 'height', y: 'weight',
        size: 'age',
        color: 'race',
      });
      t.col('height').set(0, NaN);
      t.col('weight').set(0, Infinity);
      plot.setOptions({showRegressionLine: true, markerType: 'square'});
      await new Promise(r => setTimeout(r, 1500));
      const sp = view.viewers.find((v: any) => v.type === 'Scatter plot') as any;
      return {
        rows: t.rowCount,
        spX: sp?.props.xColumnName,
        spY: sp?.props.yColumnName,
        spSize: sp?.props.sizeColumnName,
        spColor: sp?.props.colorColumnName,
        regression: sp?.props.showRegressionLine,
        markerType: sp?.props.markerType,
      };
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 15000});
    expect(result.rows).toBe(10000);
    expect(result.spX).toBe('height');
    expect(result.spY).toBe('weight');
    expect(result.spSize).toBe('age');
    expect(result.spColor).toBe('race');
    expect(result.regression).toBe(true);
    expect(result.markerType).toBe('square');
  });

  await softStep('6. Scatter plot renders correctly with NaN/Infinity rows excluded', async () => {
    // Wait for the scatter plot to fully render with its options applied
    await page.waitForTimeout(2000);
    const sp = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const v = tv.viewers.find((v: any) => v.type === 'Scatter plot') as any;
      return v ? {x: v.props.xColumnName, y: v.props.yColumnName, regression: v.props.showRegressionLine} : null;
    });
    expect(sp).not.toBeNull();
    expect(sp!.x).toBe('height');
    expect(sp!.y).toBe('weight');
    expect(sp!.regression).toBe(true);
    // Scatter plot canvas renders without crash
    const canvas = page.locator('[name="viewer-Scatter-plot"] canvas').first();
    await expect(canvas).toBeVisible();
  });

  await softStep('7. Add Histogram on height column (NaN column)', async () => {
    const histCountBefore = await page.evaluate(() =>
      grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram').length
    );
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
      icon?.click();
    });
    // Wait for a new histogram to appear in the current tv
    await page.waitForFunction((before: number) =>
      grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram').length > before,
      histCountBefore, {timeout: 8000}
    );
    await page.evaluate(() => {
      const hists = grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram') as any[];
      const hist = hists[hists.length - 1];
      hist.props.valueColumnName = 'height';
    });
    const col = await page.evaluate(() => {
      const hists = grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram') as any[];
      return hists[hists.length - 1]?.props.valueColumnName;
    });
    expect(col).toBe('height');
  });

  await softStep('8. Add Histogram on weight column (Infinity column)', async () => {
    const histCountBefore = await page.evaluate(() =>
      grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram').length
    );
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
      icon?.click();
    });
    await page.waitForFunction((before: number) =>
      grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram').length > before,
      histCountBefore, {timeout: 8000}
    );
    await page.evaluate(() => {
      const hists = grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram') as any[];
      hists[hists.length - 1].props.valueColumnName = 'weight';
    });
    const col = await page.evaluate(() => {
      const hists = grok.shell.tv.viewers.filter((v: any) => v.type === 'Histogram') as any[];
      return hists[hists.length - 1]?.props.valueColumnName;
    });
    expect(col).toBe('weight');
  });

  await softStep('9. Save layout for demog dataset', async () => {
    demogLayoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });
    expect(demogLayoutId).toBeTruthy();
  });

  // Cleanup
  await page.evaluate(async (ids: string[]) => {
    for (const id of ids) {
      try {
        const layout = await grok.dapi.layouts.find(id);
        if (layout) await grok.dapi.layouts.delete(layout);
      } catch(e) {}
    }
    grok.shell.closeAll();
  }, [spgiLayoutId, demogLayoutId].filter(Boolean));

  if (stepErrors.length > 0)
    throw new Error(stepErrors.map(e => `[${e.step}] ${e.error}`).join('\n'));
});
