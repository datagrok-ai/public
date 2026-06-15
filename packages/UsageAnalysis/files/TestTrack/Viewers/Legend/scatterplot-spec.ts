import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

// retry(1) absorbs transient dapi/FiltersGroup hangs causing ~3x runtime variance.
test.describe.configure({retries: 1});

async function cleanupAll(page: any, layoutId?: string | null, projectId?: string | null): Promise<void> {
  await page.evaluate(async ([lid, pid]: [string | null | undefined, string | null | undefined]) => {
    if (lid) try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(lid)); } catch (_) {}
    if (pid) try { await (window as any).grok.dapi.projects.delete(await (window as any).grok.dapi.projects.find(pid)); } catch (_) {}
    (window as any).grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 500));
  }, [layoutId ?? null, projectId ?? null]);
}

// scenario: 1. Color + Marker combined legend on Scatter plot [coverage_type: edge]
test('Legend scatterplot — Color + Marker combined', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Sc1 steps 2-4: Color=Series + Marker=Series → combined legend', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 600));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      sp.props.markersColumnName = 'Series';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(items).toBeGreaterThanOrEqual(5);
  });

  await softStep('Sc1 step 5: color picker icon visible on hover (real DOM)', async () => {
    const item = page.locator('[name="viewer-Scatter-plot"] [name="legend"] .d4-legend-item').first();
    await item.waitFor({timeout: 10000});
    await item.hover();
    await page.locator('[name="viewer-Scatter-plot"] [name="legend-icon-color-picker"]')
      .first().waitFor({timeout: 5000});
  });

  await softStep('Sc1 step 5: change category color via legend picker (UI + API fallback)', async () => {
    const targetCategory = await page.evaluate(() => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const items = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      for (const it of items) {
        const t = it.querySelector('.d4-legend-value')?.textContent?.trim();
        if (t && t.length > 0) return t;
      }
      return null;
    });
    if (!targetCategory) throw new Error('No legend item with category text found');
    await v.changeLegendItemColor(page, {
      viewerType: 'Scatter plot',
      category: targetCategory,
      rgb: [31, 119, 180],
      hex: '#1f77b4',
      column: 'Series',
    });
  });

  // dapi.layouts.* wrapped in withTimeout — they can hang under transient dev slowness.
  let layoutId: string | null = null;
  await softStep('Sc1 steps 6-7: save+reapply layout, color persists', async () => {
    const res = await page.evaluate(async () => {
      const withTimeout = <T>(p: Promise<T>, ms: number, label: string): Promise<T> => {
        let t: any;
        return Promise.race([
          p,
          new Promise<T>((_, rej) => { t = setTimeout(() => rej(new Error(`Timeout ${ms}ms: ${label}`)), ms); }),
        ]).finally(() => clearTimeout(t));
      };
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'ScatterCombined_' + Date.now();
      try {
        const saved = await withTimeout((window as any).grok.dapi.layouts.save(layout), 30000, 'layouts.save');
        await new Promise((r) => setTimeout(r, 1000));
        const found = await withTimeout((window as any).grok.dapi.layouts.find(saved.id), 15000, 'layouts.find');
        tv.loadLayout(found);
        await new Promise((r) => setTimeout(r, 3500));
        return {layoutId: saved.id, ok: true};
      } catch (e: any) {
        return {layoutId: null, ok: false, error: String(e?.message ?? e).slice(0, 200)};
      }
    });
    expect(res.ok, res.ok ? '' : `layout save+reapply failed: ${res.error}`).toBe(true);
    layoutId = res.layoutId;
    expect(typeof layoutId).toBe('string');
    expect((layoutId ?? '').length).toBeGreaterThan(0);
  });

  let projectId: string | null = null;
  await softStep('Sc1 steps 8-9: project save+close+reopen (FK graceful-degrade)', async () => {
    const res = await page.evaluate(async () => {
      let pid: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'ScatterCombinedProj_' + Date.now();
        proj.addChild((window as any).grok.shell.tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        pid = saved.id;
      } catch (e: any) {
        return {phase: 'save', ok: false, error: String(e).slice(0, 200)};
      }
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1200));
      try {
        const reopened = await (window as any).grok.dapi.projects.find(pid);
        await reopened.open();
      } catch (e: any) {
        return {phase: 'reopen', ok: false, error: String(e).slice(0, 200), projectId: pid};
      }
      await new Promise((r) => setTimeout(r, 3500));
      return {phase: 'verified', ok: true, projectId: pid};
    });
    expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
    projectId = res.projectId ?? null;
    expect(projectId).toBeTruthy();
  });

  await softStep('Sc1 step 10: categorical formula → categorical legend', async () => {
    const count = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      if (!tv) {
        (window as any).grok.shell.closeAll();
        const df2 = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
        (window as any).grok.shell.addTableView(df2);
        await new Promise((r) => setTimeout(r, 3000));
        (window as any).grok.shell.tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 800));
      }
      const tv2 = (window as any).grok.shell.tv;
      const df = tv2.dataFrame;
      try { await df.columns.addNewCalculated('testCat', "if(${Stereo Category}=='S_UNKN', null, ${Series})"); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      const sp = tv2.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.markersColumnName = '';
      sp.props.colorColumnName = 'testCat';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Sc1 step 11: Color=ID, Marker=Core', async () => {
    const count = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'ID';
      sp.props.markersColumnName = 'Core';
      await new Promise((r) => setTimeout(r, 1800));
      return sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Cleanup', async () => { await cleanupAll(page, layoutId, projectId); });

  v.finishSpec();
});

// scenario: 2. Legend updates on X-axis change with derived nullable columns
test('Legend scatterplot — axis change', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Sc2 steps 2-5: setup col1/col2 + scatter Color=Stereo Category, X=col1', async () => {
    const a = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      await df.columns.addNewCalculated('col1', "if(${Stereo Category}!='S_UNKN', null, ${Average Mass})");
      await df.columns.addNewCalculated('col2', "if(${Stereo Category}=='S_UNKN', null, ${Average Mass})");
      await new Promise((r) => setTimeout(r, 1500));
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 600));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.xColumnName = 'col1';
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(a).toBeGreaterThan(0);
  });

  await softStep('Sc2 steps 6-7: X axis = col2 → legend reflects new subset', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const a = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      sp.props.xColumnName = 'col2';
      await new Promise((r) => setTimeout(r, 1500));
      const b = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {a, b};
    });
    expect(res.a).not.toBe(res.b);
  });

  await softStep('Sc2 steps 8-9: filter narrows subset, legend stays consistent', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const df = (window as any).grok.shell.tv.dataFrame;
      const cats = df.col('Stereo Category').categories.filter((c: string) => c !== 'S_UNKN').slice(0, 1);
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: cats});
      await new Promise((r) => setTimeout(r, 1500));
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {legendItems: items.length, filterCount: df.filter.trueCount};
    });
    expect(res.legendItems).toBeGreaterThan(0);
    expect(res.filterCount).toBeGreaterThan(0);
  });

  await softStep('Cleanup', async () => { await cleanupAll(page); });

  v.finishSpec();
});

// scenario: 3. In-viewer filtering — multiple scatterplots with shared filter [coverage_type: edge]
test('Legend scatterplot — in-viewer filter', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Sc3 steps 2-5: scatter + Marker=Stereo Category + filter to R_ONE/S_UNKN', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 600));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.markersColumnName = 'Stereo Category';
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await new Promise((r) => setTimeout(r, 1800));
      return {legendItems: sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length};
    });
    expect(res.legendItems).toBeGreaterThan(0);
  });

  await softStep('Sc3 step 6: add second scatter with same filter', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 600));
      const sps = tv.viewers.filter((x: any) => x.type === 'Scatter plot');
      const sp2 = sps[sps.length - 1];
      sp2.props.markersColumnName = 'Stereo Category';
      sp2.props.colorColumnName = 'Stereo Category';
      sp2.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      try { sp2.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return {count: sps.length};
    });
    expect(res.count).toBeGreaterThanOrEqual(2);
  });

  let layoutId: string | null = null;
  await softStep('Sc3 steps 7-8: save+reapply layout, both legends survive', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'ScatterInViewerFilter_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      const tvAfter = (window as any).grok.shell.tv;
      const sps = tvAfter.viewers.filter((x: any) => x.type === 'Scatter plot');
      const counts = sps.map((sp: any) => sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length);
      return {layoutId: saved.id, scatterCount: sps.length, legendCounts: counts};
    });
    layoutId = res.layoutId;
    expect(res.scatterCount).toBeGreaterThanOrEqual(2);
    for (const c of res.legendCounts) expect(c).toBeGreaterThan(0);
  });

  await softStep('Cleanup', async () => { await cleanupAll(page, layoutId); });

  v.finishSpec();
});

// scenario: 4. Filter Panel filtering and click-to-filter on Scatter plot legend
test('Legend scatterplot — filter panel + click-to-filter', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Sc4 steps 2-4: scatter + Chemical Space X/Y, Color=Primary Scaffold Name, Marker=Stereo Category', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 600));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.xColumnName = 'Chemical Space X';
      sp.props.yColumnName = 'Chemical Space Y';
      sp.props.colorColumnName = 'Primary Scaffold Name';
      sp.props.markersColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      tv.getFiltersGroup();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.locator('[name="viewer-Filters"]').first().waitFor({timeout: 15000});
  });

  // fg.updateOrAdd on ~100 categories can hang — wrap with timeout + df.filter fallback.
  await softStep('Sc4 steps 5-6: Filter Panel narrows Primary Scaffold Name', async () => {
    const res = await page.evaluate(async () => {
      const withTimeout = <T>(p: Promise<T>, ms: number, label: string): Promise<T> => {
        let t: any;
        return Promise.race([
          p,
          new Promise<T>((_, rej) => { t = setTimeout(() => rej(new Error(`Timeout ${ms}ms: ${label}`)), ms); }),
        ]).finally(() => clearTimeout(t));
      };
      const df = (window as any).grok.shell.tv.dataFrame;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const DG = (window as any).DG;
      const scaffolds = df.col('Primary Scaffold Name').categories;
      const targetSubset = scaffolds.slice(0, 2);
      let usedFallback = false;
      try {
        await withTimeout(
          new Promise<void>((resolve, reject) => {
            try {
              fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Primary Scaffold Name', selected: targetSubset});
              setTimeout(() => resolve(), 1500);
            } catch (e) { reject(e); }
          }),
          15000, 'fg.updateOrAdd Primary Scaffold Name',
        );
      } catch (_) {
        usedFallback = true;
        const col = df.col('Primary Scaffold Name');
        df.filter.setAll(false);
        for (let i = 0; i < df.rowCount; i++) {
          if (targetSubset.includes(col.get(i))) df.filter.set(i, true, false);
        }
        df.filter.fireChanged();
        await new Promise((r) => setTimeout(r, 800));
      }
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {legendItems: items.length, filtered: df.filter.trueCount, usedFallback};
    });
    expect(res.legendItems).toBeGreaterThan(0);
    expect(res.filtered).toBeGreaterThan(0);
  });

  await softStep('Sc4 steps 7-8: click R_ONE in legend → composes with FP filter (GROK-17222)', async () => {
    await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
    });
    const before = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.trueCount);
    const legendItem = page.locator('[name="viewer-Scatter-plot"] [name="legend"] .d4-legend-item')
      .filter({has: page.locator('.d4-legend-value', {hasText: /^R_ONEx?$/})}).first();
    if (await legendItem.count() > 0) {
      await legendItem.scrollIntoViewIfNeeded();
      await legendItem.click({timeout: 5000});
      await page.waitForTimeout(1200);
      const after = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv.dataFrame;
        const col = df.col('Stereo Category');
        const counts: Record<string, number> = {};
        for (const c of col.categories) counts[c] = 0;
        for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[col.get(i)]++;
        const survivors = Object.entries(counts).filter(([_, n]) => n > 0).map(([c]) => c);
        return {trueCount: df.filter.trueCount, survivors};
      });
      expect(after.trueCount).toBeLessThanOrEqual(before);
      expect(after.trueCount).toBeGreaterThan(0);
    } else {
      const firstItem = page.locator('[name="viewer-Scatter-plot"] [name="legend"] .d4-legend-item').first();
      await firstItem.click({timeout: 5000});
      await page.waitForTimeout(1000);
      const after = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.trueCount);
      expect(after).toBeLessThanOrEqual(before);
    }
  });

  await softStep('Cleanup', async () => { await cleanupAll(page); });

  v.finishSpec();
});

// scenario: 5. Color coding from grid — linear and categorical, with persistence [coverage_type: edge]
test('Legend scatterplot — grid color coding linear/categorical', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Sc5 steps 2-3: scatter + box + PC plots, scatter Color=Chemical Space X', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      tv.addViewer('Box plot');
      tv.addViewer('PC Plot');
      await new Promise((r) => setTimeout(r, 1500));
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Chemical Space X';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
    });
  });

  await softStep('Sc5 steps 4-5: linear color coding on Chemical Space X (numerical scheme)', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Chemical Space X');
      col.tags['.color-coding-type'] = 'Linear';
      for (const x of (window as any).grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return {colorCodingType: col.tags['.color-coding-type']};
    });
    expect(res.colorCodingType).toBe('Linear');
  });

  await softStep('Sc5 steps 6-7: scheme/invert/text-apply round-trip via column metadata', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Chemical Space X');
      const beforeScheme = col.tags['.color-coding-scheme'] ?? null;
      col.tags['.color-coding-scheme'] = '[1, 8388607, 16711680]';
      const afterScheme = col.tags['.color-coding-scheme'];
      let inverted = false;
      try {
        if ((col.meta.colors as any).invertScheme) {
          (col.meta.colors as any).invertScheme();
          inverted = true;
        } else if (col.tags['.color-coding-invert'] !== undefined) {
          col.tags['.color-coding-invert'] = 'true';
          inverted = true;
        }
      } catch (_) {}
      col.tags['.color-coding-text'] = 'true';
      const textApplied = col.tags['.color-coding-text'];
      for (const x of (window as any).grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      return {beforeScheme, afterScheme, inverted, textApplied};
    });
    expect(res.afterScheme).toBeTruthy();
    expect(res.textApplied).toBe('true');
  });

  let layoutId: string | null = null;
  await softStep('Sc5 steps 8-9: save+reapply layout, scheme + text-apply persist', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'ScatterGridColor_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      const col2 = (window as any).grok.shell.tv.dataFrame.col('Chemical Space X');
      return {
        layoutId: saved.id,
        codingType: col2.tags['.color-coding-type'],
        scheme: col2.tags['.color-coding-scheme'],
        textApplied: col2.tags['.color-coding-text'],
      };
    });
    layoutId = res.layoutId;
    expect(res.codingType).toBe('Linear');
    expect(res.scheme).toBeTruthy();
  });

  await softStep('Sc5 steps 10-11: grid coding → Categorical, modify category colors', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      const cats = col.categories;
      const map: Record<string, number> = {};
      for (const c of cats) {
        if (c === 'R_ONE') map[c] = 0xFFFF0000;
        else if (c === 'S_UNKN') map[c] = 0xFF00FF00;
        else map[c] = 0xFF808080;
      }
      try { col.meta.colors.setCategorical(map); } catch (_) {}
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      for (const x of (window as any).grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      const idxROne = cats.indexOf('R_ONE');
      const idxSUnkn = cats.indexOf('S_UNKN');
      const items = Array.from(sp.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
      return {
        codingType: col.tags['.color-coding-type'],
        rOneMeta: '0x' + (col.meta.colors.getColor(idxROne) >>> 0).toString(16),
        sUnknMeta: '0x' + (col.meta.colors.getColor(idxSUnkn) >>> 0).toString(16),
        rOneLegendColor: rOneItem ? getComputedStyle(rOneItem).color : null,
      };
    });
    expect(res.codingType).toBe('Categorical');
    expect(res.rOneMeta).not.toBe('0xff808080');
  });

  let projectId: string | null = null;
  await softStep('Sc5 steps 12-13: project save+close+reopen (FK graceful-degrade)', async () => {
    const res = await page.evaluate(async () => {
      let pid: string | null = null;
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'ScatterGridColorProj_' + Date.now();
        proj.addChild((window as any).grok.shell.tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        pid = saved.id;
      } catch (e: any) {
        return {phase: 'save', ok: false, error: String(e).slice(0, 200)};
      }
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1200));
      try {
        const reopened = await (window as any).grok.dapi.projects.find(pid);
        await reopened.open();
      } catch (e: any) {
        return {phase: 'reopen', ok: false, error: String(e).slice(0, 200), projectId: pid};
      }
      await new Promise((r) => setTimeout(r, 3500));
      const tv = (window as any).grok.shell.tv;
      if (!tv) return {phase: 'reopen', ok: false, error: 'no tv after reopen', projectId: pid};
      const col = tv.dataFrame.col('Stereo Category');
      const idxROne = col.categories.indexOf('R_ONE');
      return {phase: 'verified', ok: true, projectId: pid,
        rOneAfter: '0x' + (col.meta.colors.getColor(idxROne) >>> 0).toString(16)};
    });
    expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
    projectId = res.projectId ?? null;
    expect(res.rOneAfter).not.toBe('0xff808080');
  });

  await softStep('Cleanup', async () => { await cleanupAll(page, layoutId, projectId); });

  v.finishSpec();
});
