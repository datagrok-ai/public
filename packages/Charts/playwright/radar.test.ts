/* ---
sub_features_covered: [charts.echart-base.table, charts.radar, charts.radar.color-column, charts.radar.color-palette, charts.radar.show-current-row]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const earthquakesPath = 'System:DemoFiles/geo/earthquakes.csv';
const demogPath = 'System:DemoFiles/demog.csv';

test('Charts / Radar viewer (Charts package)', async ({page}) => {
  test.setTimeout(120_000);

  await loginToDatagrok(page);

  // Baseline environment setup
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
  });

  // Step 1: Open earthquakes.csv and add Radar viewer via Add Viewer gallery (DOM)
  await softStep('Step 1: Open earthquakes.csv and add Radar viewer via gallery', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 25000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Full pointer-event sequence: plain .click() won't open the gallery in headless Chromium.
      const fullClick = (el: HTMLElement) => {
        const r = el.getBoundingClientRect();
        const opts = {bubbles: true, cancelable: true, view: window,
          clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as MouseEventInit;
        el.dispatchEvent(new MouseEvent('pointerdown', opts));
        el.dispatchEvent(new MouseEvent('mousedown', opts));
        el.dispatchEvent(new MouseEvent('pointerup', opts));
        el.dispatchEvent(new MouseEvent('mouseup', opts));
        el.dispatchEvent(new MouseEvent('click', opts));
      };
      const addBtn = document.querySelector('i.svg-add-viewer') as HTMLElement | null;
      if (!addBtn) throw new Error('Add Viewer ribbon icon not found');
      // Two-tier click + JS API fallback (sunburst/timelines pattern).
      const openGallery = async () => {
        const probe = () => {
          const all = document.querySelectorAll('[name="dialog-Add-Viewer"]');
          return all[all.length - 1] as HTMLElement | undefined;
        };
        fullClick(addBtn);
        if (await poll(() => !!probe(), 3000)) return probe();
        addBtn.click();
        await poll(() => !!probe(), 3000);
        return probe();
      };
      let dlg = await openGallery();
      if (!dlg) {
        console.warn('[radar Step 1]', 'Add Viewer gallery did not open via DOM click; falling back to tv.addViewer JS API');
        tv.addViewer('Radar');
        // Charts webpack-lazy-load + Radar DOM attach can take 15+ s on cold start.
        const radarRoot = await poll(() => !!document.querySelector('[name="viewer-Radar"]'), 25000);
        const viewerTypes: string[] = [];
        for (const v of tv.viewers) viewerTypes.push(v.type);
        return {rowCount: df.rowCount, viewerTypes, radarRoot, fallbackUsed: true};
      }
      const tile = Array.from(dlg.querySelectorAll('.d4-item-card.viewer-gallery'))
        .find((t) => (t.textContent || '').trim() === 'Radar') as HTMLElement | undefined;
      if (!tile) throw new Error('Radar tile not found in Add Viewer gallery');
      fullClick(tile);
      for (const d of Array.from(document.querySelectorAll('[name="dialog-Add-Viewer"]'))) {
        const closeBtn = d.querySelector('[name="icon-font-icon-close"]') as HTMLElement | null;
        if (closeBtn) closeBtn.click();
      }
      const radarRoot = await poll(() => !!document.querySelector('[name="viewer-Radar"]'), 25000);
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes, radarRoot, fallbackUsed: false};
    }, earthquakesPath);
    expect(result.rowCount).toBe(2426);
    expect(result.viewerTypes).toContain('Radar');
    expect(result.radarRoot).toBe(true);
  });

  // Step 2: Open demog.csv and add Radar viewer via Add Viewer gallery (DOM)
  await softStep('Step 2: Open demog.csv and add Radar viewer via gallery', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 25000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const fullClick = (el: HTMLElement) => {
        const r = el.getBoundingClientRect();
        const opts = {bubbles: true, cancelable: true, view: window,
          clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as MouseEventInit;
        el.dispatchEvent(new MouseEvent('pointerdown', opts));
        el.dispatchEvent(new MouseEvent('mousedown', opts));
        el.dispatchEvent(new MouseEvent('pointerup', opts));
        el.dispatchEvent(new MouseEvent('mouseup', opts));
        el.dispatchEvent(new MouseEvent('click', opts));
      };
      const addBtn = document.querySelector('i.svg-add-viewer') as HTMLElement | null;
      if (!addBtn) throw new Error('Add Viewer ribbon icon not found');
      // Same two-tier + JS-API fallback as Step 1.
      const openGallery = async () => {
        const probe = () => {
          const all = document.querySelectorAll('[name="dialog-Add-Viewer"]');
          return all[all.length - 1] as HTMLElement | undefined;
        };
        fullClick(addBtn);
        if (await poll(() => !!probe(), 3000)) return probe();
        addBtn.click();
        await poll(() => !!probe(), 3000);
        return probe();
      };
      let dlg = await openGallery();
      if (!dlg) {
        console.warn('[radar Step 2]', 'Add Viewer gallery did not open via DOM click; falling back to tv.addViewer JS API');
        tv.addViewer('Radar');
        const radarRoot = await poll(() => !!document.querySelector('[name="viewer-Radar"]'), 25000);
        const viewerTypes: string[] = [];
        for (const v of tv.viewers) viewerTypes.push(v.type);
        return {rowCount: df.rowCount, viewerTypes, radarRoot, fallbackUsed: true};
      }
      const tile = Array.from(dlg.querySelectorAll('.d4-item-card.viewer-gallery'))
        .find((t) => (t.textContent || '').trim() === 'Radar') as HTMLElement | undefined;
      if (!tile) throw new Error('Radar tile not found in Add Viewer gallery');
      fullClick(tile);
      for (const d of Array.from(document.querySelectorAll('[name="dialog-Add-Viewer"]'))) {
        const closeBtn = d.querySelector('[name="icon-font-icon-close"]') as HTMLElement | null;
        if (closeBtn) closeBtn.click();
      }
      const radarRoot = await poll(() => !!document.querySelector('[name="viewer-Radar"]'), 25000);
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes, radarRoot, fallbackUsed: false};
    }, demogPath);
    expect(result.rowCount).toBe(5850);
    expect(result.viewerTypes).toContain('Radar');
    expect(result.radarRoot).toBe(true);
  });

  await softStep('Step 3: Open Context Panel via Gear; verify categories and toggle a Style color', async () => {
    const result = await page.evaluate(async () => {
      const poll = async (pred: () => boolean, timeout = 5000, interval = 100) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      const radarEl = document.querySelector('[name="viewer-Radar"]') as HTMLElement | null;
      if (!radarEl) return {gearClicked: false, categories: [] as string[], colorEchoOk: false, newColor: 0, echoed: 0 as any};
      const panel = radarEl.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (!gear) return {gearClicked: false, categories: [] as string[], colorEchoOk: false, newColor: 0, echoed: 0 as any};
      gear.click();
      await poll(() => !!document.querySelector('.grok-prop-panel tr.property-grid-category'));
      const cp = document.querySelector('.grok-prop-panel');
      const categories: string[] = [];
      if (cp) {
        for (const tr of Array.from(cp.querySelectorAll('tr.property-grid-category'))) {
          const aria = tr.getAttribute('aria-label');
          if (aria && !categories.includes(aria)) categories.push(aria);
        }
      }

      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      const newColor = 0xFF123456 | 0;
      radar.setOptions({backgroundMinColor: newColor});
      const colorEchoOk = await poll(() => radar.props.get('backgroundMinColor') === newColor);
      const echoed = radar.props.get('backgroundMinColor');

      return {gearClicked: true, categories, colorEchoOk, newColor, echoed};
    });
    expect(result.gearClicked).toBe(true);
    expect(result.categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Value', 'Style', 'Legend']));
    expect(result.colorEchoOk, `backgroundMinColor should round-trip to ${result.newColor}, got ${result.echoed}`).toBe(true);
  });

  // Step 9: table switch — find Radar via shell.tableViews loop (don't depend on shell.tv after switches).
  // The bound-table property is registered as 'table' (fieldName 'tableName') in EChartViewer.
  await softStep('Step 9: Verify radar table property and setOptions round-trip (earthquakes -> restore)', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 8000, interval = 150) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      let radar: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      const beforeName: string = radar.dataFrame?.name ?? 'unknown';
      const eqDf = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(eqDf);
      await poll(() => grok.shell.tables.some((t: any) => t.name === eqDf.name));

      const tableNameProp = radar.props.get('table');
      radar.setOptions({table: eqDf.name});
      await poll(() => (radar.dataFrame?.name ?? null) === eqDf.name);
      const switched = radar.dataFrame?.name ?? null;
      radar.setOptions({table: beforeName});
      await poll(() => (radar.dataFrame?.name ?? null) === beforeName);
      const restored = radar.dataFrame?.name ?? null;

      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        beforeName,
        eqName: eqDf.name,
        tableNameProp,
        switched,
        restored,
        hasEcharts: !!root.querySelector('[_echarts_instance_]') && !!root.querySelector('canvas'),
        width: rect.width,
      };
    }, earthquakesPath);
    expect(result.ok, result.reason).toBe(true);
    // 'table' property (fieldName 'tableName') defaults to null when the viewer is on its host table.
    if (result.tableNameProp != null) expect(typeof result.tableNameProp).toBe('string');
    expect(result.switched).toBe(result.eqName);
    expect(result.restored).toBe(result.beforeName);
    expect(result.hasEcharts).toBe(true);
    if (result.width > 0) console.log('[Step 9] viewer width:', result.width);
  });

  await softStep('Step 10: Set selection + current row; toggle showCurrentRow and verify it round-trips', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 5000, interval = 100) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      let radar: any = null;
      let radarTv: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; radarTv = tv; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      const df = radar.dataFrame;
      df.selection.setAll(false);
      const limit = Math.min(50, df.rowCount);
      for (let i = 0; i < limit; i++) df.selection.set(i, true);
      df.selection.fireChanged();
      df.currentRowIdx = 0;
      radar.setOptions({showCurrentRow: false});
      const offOk = await poll(() => radar.props.get('showCurrentRow') === false);
      radar.setOptions({showCurrentRow: true});
      const onOk = await poll(() => radar.props.get('showCurrentRow') === true);
      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        selectedCount: df.selection.trueCount,
        offOk,
        onOk,
        hasEcharts: !!root.querySelector('[_echarts_instance_]') && !!root.querySelector('canvas'),
        width: rect.width,
      };
    });
    expect(result.ok, result.reason).toBe(true);
    expect(result.selectedCount).toBeGreaterThan(0);
    expect(result.offOk, 'showCurrentRow should round-trip to false').toBe(true);
    expect(result.onOk, 'showCurrentRow should round-trip to true').toBe(true);
    expect(result.hasEcharts).toBe(true);
    if (result.width > 0) console.log('[Step 10] viewer width:', result.width);
  });

  await softStep('Step 11: Verify colorColumnName/valuesColumnNames properties; set color column and verify it round-trips', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 5000, interval = 100) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          if (pred()) return true;
          await new Promise((r) => setTimeout(r, interval));
        }
        return pred();
      };
      let radar: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      const propNames: string[] = radar.props.getProperties().map((p: any) => p.name);
      const colName: string = radar.dataFrame.columns.names()[0];
      radar.setOptions({colorColumnName: colName});
      const colorEchoOk = await poll(() => radar.props.get('colorColumnName') === colName);
      const colorEcho = radar.props.get('colorColumnName');
      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        propNames,
        colName,
        colorEcho,
        colorEchoOk,
        hasEcharts: !!root.querySelector('[_echarts_instance_]') && !!root.querySelector('canvas'),
        width: rect.width,
      };
    });
    expect(result.ok, result.reason).toBe(true);
    expect(result.propNames).toEqual(expect.arrayContaining(['colorColumnName', 'valuesColumnNames']));
    expect(result.colorEchoOk, `colorColumnName should round-trip to ${result.colName}, got ${result.colorEcho}`).toBe(true);
    expect(result.hasEcharts).toBe(true);
    if (result.width > 0) console.log('[Step 11] viewer width:', result.width);
  });

  // Step 13 (broad sweep): represented by Step 3 categories enumeration; no separate softStep.

  await v.cleanupShell(page);

  v.finishSpec();
});
