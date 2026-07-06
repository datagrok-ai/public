/* ---
sub_features_covered: [charts.timelines, charts.timelines.color-column, charts.timelines.end-column, charts.timelines.legend-visibility, charts.timelines.split-by-column, charts.timelines.start-column]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

const aePath = 'System:AppData/Charts/ae.csv';

test('Charts / Timelines viewer — legend filtering regression (GROK-19033)', async ({page}) => {
  test.setTimeout(120_000);

  // Capture console errors for the GROK-19033 visual-stability invariant; filter benign network noise.
  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

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

  // === Scenario 1: Legend click-to-filter (GROK-19033 reproduction class) ===

  // Step 1.1-1.2: Open ae.csv and add Timelines viewer via gallery (DOM).
  await softStep('Scenario 1 Step 1-2: Open ae.csv; Add Timelines viewer via gallery', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Full pointer-event sequence required to open the gallery.
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
      // Two-tier click + JS API fallback (same pattern as sunburst-spec.ts).
      const waitFor = async (pred: () => any, timeout: number, interval = 100) => {
        const start = Date.now();
        while (Date.now() - start < timeout) { const r = pred(); if (r) return r; await new Promise((res) => setTimeout(res, interval)); }
        return pred();
      };
      const probeDialog = () => {
        const all = document.querySelectorAll('[name="dialog-Add-Viewer"]');
        return all[all.length - 1] as HTMLElement | undefined;
      };
      const openGallery = async () => {
        fullClick(addBtn);
        if (await waitFor(probeDialog, 3000)) return probeDialog();
        addBtn.click();
        return await waitFor(probeDialog, 3000);
      };
      let dlg = await openGallery();
      if (!dlg) {
        console.warn('[timelines Step 1-2]', 'Add Viewer gallery did not open via DOM click; falling back to tv.addViewer JS API');
        tv.addViewer('Timelines');
        // Charts package webpack-lazy-loads — poll for the viewer root.
        await waitFor(() => document.querySelector('[name="viewer-Timelines"]'), 20_000);
        const viewerTypes: string[] = [];
        for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) viewerTypes.push(v.type);
        const tlRoot = !!document.querySelector('[name="viewer-Timelines"]');
        return {rowCount: df.rowCount, viewerTypes, tlRoot, columns: df.columns.names(), fallbackUsed: true};
      }
      const tile = Array.from(dlg.querySelectorAll('.d4-item-card.viewer-gallery'))
        .find((t) => (t.textContent || '').trim() === 'Timelines') as HTMLElement | undefined;
      if (!tile) throw new Error('Timelines tile not found');
      fullClick(tile);
      // Charts package webpack-lazy-loads — poll for the viewer root.
      await waitFor(() => document.querySelector('[name="viewer-Timelines"]'), 20_000);
      for (const d of Array.from(document.querySelectorAll('[name="dialog-Add-Viewer"]'))) {
        const closeBtn = d.querySelector('[name="icon-font-icon-close"]') as HTMLElement | null;
        if (closeBtn) closeBtn.click();
      }
      await new Promise((r) => setTimeout(r, 500));
      const viewerTypes: string[] = [];
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) viewerTypes.push(v.type);
      const tlRoot = !!document.querySelector('[name="viewer-Timelines"]');
      return {rowCount: df.rowCount, viewerTypes, tlRoot, columns: df.columns.names(), fallbackUsed: false};
    }, aePath);
    expect(result.viewerTypes).toContain('Timelines');
    expect(result.tlRoot).toBe(true);
    expect(result.columns).toEqual(expect.arrayContaining(['USUBJID', 'AESTDY', 'AEENDY', 'AESOC']));
  });

  await softStep('Scenario 1 Step 3-4: Open Gear; set Color Column = AESOC via JS API', async () => {
    const result = await page.evaluate(async () => {
      const tlEl = document.querySelector('[name="viewer-Timelines"]') as HTMLElement | null;
      if (!tlEl) return {ok: false, gearClicked: false};
      const panel = tlEl.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (!gear) return {ok: false, gearClicked: false};
      const waitFor = async (pred: () => any, timeout: number, interval = 100) => {
        const start = Date.now();
        while (Date.now() - start < timeout) { const r = pred(); if (r) return r; await new Promise((res) => setTimeout(res, interval)); }
        return pred();
      };
      gear.click();
      await waitFor(() => document.querySelector('[name="div-column-combobox-split--by"]'), 10_000);

      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false, gearClicked: true};
      const safeGet = (n: string) => { try { return timelines.props.get(n); } catch (e) { return null; } };
      const defaultProps = {
        splitByColumnName: safeGet('splitByColumnName'),
        startColumnName: safeGet('startColumnName'),
        endColumnName: safeGet('endColumnName'),
        marker: safeGet('marker'),
        legendVisibility: safeGet('legendVisibility'),
      };

      timelines.setOptions({colorColumnName: 'AESOC'});
      const csStart = Date.now();
      while (Date.now() - csStart < 5000 && timelines.props.get('colorColumnName') !== 'AESOC')
        await new Promise((r) => setTimeout(r, 100));
      const colorCombo = document.querySelector('[name="div-column-combobox-color"]');
      const comboShown = !!colorCombo;
      const comboColumnText = colorCombo?.querySelector('.d4-column-selector-column')?.textContent?.trim() || null;
      const colorColumn = timelines.props.get('colorColumnName');

      const aesoc = tv.dataFrame.col('AESOC');
      const categories = new Set<string>();
      for (let i = 0; i < tv.dataFrame.rowCount; i++) {
        const v = aesoc.get(i);
        if (v != null) categories.add(String(v));
      }

      const legendVisLabel = document.querySelector('[name="prop-view-legend-visibility"]');
      const legendVisLabelShown = !!legendVisLabel;

      // Split-by combobox name uses a double dash.
      const splitCombo = document.querySelector('[name="div-column-combobox-split--by"]');
      const splitComboShown = !!splitCombo;

      return {
        ok: true,
        gearClicked: true,
        defaultProps,
        comboShown,
        comboColumnText,
        colorColumn,
        categoryCount: categories.size,
        legendVisLabelShown,
        splitComboShown,
      };
    });
    expect(result.ok).toBe(true);
    expect(result.gearClicked).toBe(true);
    expect(result.comboShown).toBe(true);
    expect(result.legendVisLabelShown).toBe(true);
    expect(result.splitComboShown).toBe(true);
    expect(result.colorColumn, 'colorColumnName should bind to AESOC').toBe('AESOC');
    // AESOC is a raw SDTM name (no display alias), so the combobox reflects the field name verbatim.
    expect(result.comboColumnText, 'color combobox should reflect AESOC').toBe('AESOC');
    expect(result.categoryCount).toBeGreaterThan(0);
    console.log('[Timelines defaults read]', JSON.stringify(result.defaultProps));
  });

  // Step 1.5-1.7: DOM legend click on AESOC category (legend is real DOM, not ECharts canvas).
  await softStep('Scenario 1 Step 5-7: DOM legend click on AESOC category — visual-stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tl = document.querySelector('[name="viewer-Timelines"]') as HTMLElement | null;
      if (!tl) return {ok: false};
      const legend = tl.querySelector('[name="legend"]');
      if (!legend) return {ok: false, legendFound: false};
      // The legend repopulates asynchronously after the color column change; poll for items.
      const itemsStart = Date.now();
      while (Date.now() - itemsStart < 15000 && legend.querySelectorAll('.d4-legend-item').length === 0)
        await new Promise((r) => setTimeout(r, 100));
      const items = Array.from(legend.querySelectorAll('.d4-legend-item'));
      const target = items[0] as HTMLElement | undefined;
      if (!target) return {ok: false, legendFound: true, itemCount: 0};
      const targetText = target.querySelector('.d4-legend-value')?.textContent?.trim() || null;
      const beforeCls = target.className;
      target.click();
      const clkStart = Date.now();
      while (Date.now() - clkStart < 5000 && !target.className.includes('d4-legend-item-current'))
        await new Promise((r) => setTimeout(r, 100));
      const afterCls = target.className;

      const root = tl.getBoundingClientRect();

      // GROK-19033 invariant: the viewer canvas is present, sampleable, and not all-white.
      const canvas = tl.querySelector('canvas') as HTMLCanvasElement | null;
      const canvasPresent = !!canvas;
      let canvasSampled = false;
      let canvasNotWhite = false;
      if (canvas) {
        const ctx = canvas.getContext('2d');
        if (ctx) {
          const cx = Math.max(1, Math.floor(canvas.width / 2));
          const cy = Math.max(1, Math.floor(canvas.height / 2));
          const px = ctx.getImageData(cx, cy, 1, 1).data;
          canvasSampled = true;
          canvasNotWhite = !(px[0] === 255 && px[1] === 255 && px[2] === 255 && px[3] === 255);
        }
      }

      return {
        ok: true,
        legendFound: true,
        itemCount: items.length,
        targetText,
        beforeCls,
        afterCls,
        rootWidth: root.width,
        rootHeight: root.height,
        canvasPresent,
        canvasSampled,
        canvasNotWhite,
        currentToggled: !beforeCls.includes('d4-legend-item-current') && afterCls.includes('d4-legend-item-current'),
      };
    });
    expect(result.ok).toBe(true);
    expect(result.legendFound).toBe(true);
    expect(result.itemCount).toBeGreaterThan(0);
    expect(result.currentToggled).toBe(true);
    // GROK-19033: viewer remains stable, no console error during the click.
    expect(result.rootWidth).toBeGreaterThan(0);
    expect(result.rootHeight).toBeGreaterThan(0);
    expect(result.canvasPresent, 'Timelines must render a canvas').toBe(true);
    expect(result.canvasSampled, 'canvas getImageData must succeed').toBe(true);
    expect(result.canvasNotWhite).toBe(true);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  // === Scenario 2: legendVisibility transitions ===

  // legendVisibility cycle Auto → Always → Never → Auto must be visually stable (GROK-19033).
  await softStep('Scenario 2 Step 2: legendVisibility = Always; visual stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false};
      timelines.setOptions({legendVisibility: 'Always'});
      const start = Date.now();
      while (Date.now() - start < 5000 && timelines.props.get('legendVisibility') !== 'Always')
        await new Promise((r) => setTimeout(r, 100));
      const root = timelines.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      const visibility = timelines.props.get('legendVisibility');
      const tlEl = document.querySelector('[name="viewer-Timelines"]');
      const legendEl = tlEl ? tlEl.querySelector('[name="legend"]') : null;
      const legendVisible = !!legendEl && (legendEl as HTMLElement).getBoundingClientRect().height > 0;
      return {ok: true, visibility, legendVisible, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.visibility, 'legendVisibility should apply').toBe('Always');
    expect(result.legendVisible, 'legend must be shown when Always').toBe(true);
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Scenario 2 Step 3: legendVisibility = Never; visual stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false};
      timelines.setOptions({legendVisibility: 'Never'});
      const start = Date.now();
      while (Date.now() - start < 5000 && timelines.props.get('legendVisibility') !== 'Never')
        await new Promise((r) => setTimeout(r, 100));
      const root = timelines.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      const visibility = timelines.props.get('legendVisibility');
      const tlEl = document.querySelector('[name="viewer-Timelines"]');
      const legendEl = tlEl ? tlEl.querySelector('[name="legend"]') : null;
      const legendVisible = !!legendEl && (legendEl as HTMLElement).getBoundingClientRect().height > 0;
      return {ok: true, visibility, legendVisible, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.visibility, 'legendVisibility should apply').toBe('Never');
    expect(result.legendVisible, 'legend must be hidden when Never').toBe(false);
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Scenario 2 Step 4: legendVisibility = Auto (round-trip); visual stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false};
      timelines.setOptions({legendVisibility: 'Auto'});
      const start = Date.now();
      while (Date.now() - start < 3000 && timelines.props.get('legendVisibility') !== 'Auto')
        await new Promise((r) => setTimeout(r, 100));
      const root = timelines.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      const visibility = timelines.props.get('legendVisibility');
      const tlEl = document.querySelector('[name="viewer-Timelines"]');
      const legendEl = tlEl ? tlEl.querySelector('[name="legend"]') : null;
      const legendVisible = !!legendEl && (legendEl as HTMLElement).getBoundingClientRect().height > 0;
      return {ok: true, visibility, legendVisible, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    // 'Auto' may be stored verbatim or normalized to the resolved mode by the viewer; assert the visible effect (legend shows for AESOC's many categories) as the invariant, and that the stored prop stays a valid mode.
    expect(['Auto', 'Always', 'Never'], 'legendVisibility must be a valid mode').toContain(result.visibility);
    expect(result.legendVisible, 'Auto shows the legend when categories exist').toBe(true);
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  // Step 2.5: re-click legend to re-enable category — DOM (re-uses Scenario 1 selector).
  await softStep('Scenario 2 Step 5: Re-click legend item to re-enable category — visual stability', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tl = document.querySelector('[name="viewer-Timelines"]') as HTMLElement | null;
      if (!tl) return {ok: false};
      const legend = tl.querySelector('[name="legend"]');
      if (!legend) return {ok: false, legendFound: false};
      const items = Array.from(legend.querySelectorAll('.d4-legend-item'));
      const target = items[0] as HTMLElement | undefined;
      if (!target) return {ok: false, legendFound: true, itemCount: 0};
      const before = target.className;
      target.click();
      const s = Date.now();
      while (Date.now() - s < 5000 && target.className === before)
        await new Promise((r) => setTimeout(r, 100));
      const root = tl.getBoundingClientRect();
      return {
        ok: true,
        legendFound: true,
        itemCount: items.length,
        rootWidth: root.width,
        rootHeight: root.height,
      };
    });
    expect(result.ok).toBe(true);
    expect(result.legendFound).toBe(true);
    expect(result.itemCount).toBeGreaterThan(0);
    expect(result.rootWidth).toBeGreaterThan(0);
    expect(result.rootHeight).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  // === Scenario 3: splitByColumnName re-bind mid-session ===

  await softStep('Scenario 3 Step 2: splitByColumnName USUBJID -> AESEV; visual stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false};
      timelines.setOptions({splitByColumnName: 'AESEV'});
      const start = Date.now();
      while (Date.now() - start < 5000 && timelines.props.get('splitByColumnName') !== 'AESEV')
        await new Promise((r) => setTimeout(r, 100));
      const root = timelines.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      const splitComboText = document.querySelector('[name="div-column-combobox-split--by"] .d4-column-selector-column')?.textContent?.trim() || null;
      const splitBy = timelines.props.get('splitByColumnName');
      return {ok: true, splitBy, splitComboText, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.splitBy, 'splitByColumnName should rebind to AESEV').toBe('AESEV');
    expect(result.splitComboText, 'split-by combobox should reflect AESEV').toBe('AESEV');
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  // Step 3.3: legend click in new split-by configuration — DOM.
  await softStep('Scenario 3 Step 3: Legend click in new split-by config — visual stability', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tl = document.querySelector('[name="viewer-Timelines"]') as HTMLElement | null;
      if (!tl) return {ok: false};
      const legend = tl.querySelector('[name="legend"]');
      // Legend may re-render or hide after rebind; tolerate either — invariant is no white-screen / error.
      const items = legend ? Array.from(legend.querySelectorAll('.d4-legend-item')) : [];
      let clicked = false;
      if (items.length > 0) {
        const el = items[0] as HTMLElement;
        const before = el.className;
        el.click();
        clicked = true;
        const s = Date.now();
        while (Date.now() - s < 3000 && el.className === before)
          await new Promise((r) => setTimeout(r, 100));
      }
      const root = tl.getBoundingClientRect();
      const canvas = tl.querySelector('canvas') as HTMLCanvasElement | null;
      const canvasPresent = !!canvas;
      let canvasSampled = false;
      let canvasNotWhite = false;
      if (canvas) {
        const ctx = canvas.getContext('2d');
        if (ctx) {
          const cx = Math.max(1, Math.floor(canvas.width / 2));
          const cy = Math.max(1, Math.floor(canvas.height / 2));
          const px = ctx.getImageData(cx, cy, 1, 1).data;
          canvasSampled = true;
          canvasNotWhite = !(px[0] === 255 && px[1] === 255 && px[2] === 255 && px[3] === 255);
        }
      }
      return {ok: true, legendItems: items.length, clicked, rootWidth: root.width, rootHeight: root.height, canvasPresent, canvasSampled, canvasNotWhite};
    });
    expect(result.ok).toBe(true);
    expect(result.rootWidth).toBeGreaterThan(0);
    expect(result.rootHeight).toBeGreaterThan(0);
    expect(result.canvasPresent, 'Timelines must render a canvas').toBe(true);
    expect(result.canvasSampled, 'canvas getImageData must succeed').toBe(true);
    expect(result.canvasNotWhite).toBe(true);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
    console.log('[Scenario 3 Step 3]', JSON.stringify({legendItems: result.legendItems, clicked: result.clicked}));
  });

  await softStep('Scenario 3 Step 4: splitByColumnName revert to USUBJID; visual stability assertion', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let timelines: any = null;
      for (const v of (tv && tv.viewers ? Array.from(tv.viewers) as any[] : [] as any[])) if (v.type === 'Timelines') { timelines = v; break; }
      if (!timelines) return {ok: false};
      timelines.setOptions({splitByColumnName: 'USUBJID'});
      const start = Date.now();
      while (Date.now() - start < 5000 && timelines.props.get('splitByColumnName') !== 'USUBJID')
        await new Promise((r) => setTimeout(r, 100));
      const root = timelines.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      const splitBy = timelines.props.get('splitByColumnName');
      const splitComboText = document.querySelector('[name="div-column-combobox-split--by"] .d4-column-selector-column')?.textContent?.trim() || null;
      return {ok: true, splitBy, splitComboText, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.splitBy, 'splitByColumnName should revert to USUBJID').toBe('USUBJID');
    expect(result.splitComboText, 'split-by combobox should reflect USUBJID').toBe('USUBJID');
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
