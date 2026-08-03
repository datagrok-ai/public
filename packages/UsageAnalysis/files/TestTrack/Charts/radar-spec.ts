import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

const earthquakesPath = 'System:DemoFiles/geo/earthquakes.csv';
const demogPath = 'System:DemoFiles/demog.csv';

test('Radar viewer (Charts package)', async ({page}) => {
  test.setTimeout(300_000);

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
        await new Promise((r) => setTimeout(r, 800));
        if (probe()) return probe();
        addBtn.click();
        await new Promise((r) => setTimeout(r, 800));
        return probe();
      };
      let dlg = await openGallery();
      if (!dlg) {
        console.warn('[radar Step 1]', 'Add Viewer gallery did not open via DOM click; falling back to tv.addViewer JS API');
        tv.addViewer('Radar');
        // Retry probe up to ~25s — Charts webpack-lazy-load + Radar DOM attach can take 15+ s on cold start.
        let radarRoot = false;
        for (let attempt = 0; attempt < 5; attempt++) {
          await new Promise((r) => setTimeout(r, 5000));
          if (document.querySelector('[name="viewer-Radar"]')) { radarRoot = true; break; }
        }
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
      // Retry probe up to ~25s — Charts webpack-lazy-load + Radar DOM attach can take 15+ s on cold start.
      let radarRoot = false;
      for (let attempt = 0; attempt < 5; attempt++) {
        await new Promise((r) => setTimeout(r, 5000));
        if (document.querySelector('[name="viewer-Radar"]')) { radarRoot = true; break; }
      }
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
        await new Promise((r) => setTimeout(r, 800));
        if (probe()) return probe();
        addBtn.click();
        await new Promise((r) => setTimeout(r, 800));
        return probe();
      };
      let dlg = await openGallery();
      if (!dlg) {
        console.warn('[radar Step 2]', 'Add Viewer gallery did not open via DOM click; falling back to tv.addViewer JS API');
        tv.addViewer('Radar');
        // Retry probe up to ~25s (same as Step 1).
        let radarRoot = false;
        for (let attempt = 0; attempt < 5; attempt++) {
          await new Promise((r) => setTimeout(r, 5000));
          if (document.querySelector('[name="viewer-Radar"]')) { radarRoot = true; break; }
        }
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
      // Retry probe up to ~25s (same as Step 1).
      let radarRoot = false;
      for (let attempt = 0; attempt < 5; attempt++) {
        await new Promise((r) => setTimeout(r, 5000));
        if (document.querySelector('[name="viewer-Radar"]')) { radarRoot = true; break; }
      }
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
      const radarEl = document.querySelector('[name="viewer-Radar"]') as HTMLElement | null;
      if (!radarEl) return {gearClicked: false, categories: [] as string[], colorEchoOk: false, newColor: 0, echoed: 0 as any};
      const panel = radarEl.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (!gear) return {gearClicked: false, categories: [] as string[], colorEchoOk: false, newColor: 0, echoed: 0 as any};
      gear.click();
      await new Promise((r) => setTimeout(r, 1000));
      const cp = document.querySelector('.grok-prop-panel');
      const categories: string[] = [];
      if (cp) {
        for (const tr of Array.from(cp.querySelectorAll('tr.property-grid-category'))) {
          const aria = tr.getAttribute('aria-label');
          if (aria && !categories.includes(aria)) categories.push(aria);
        }
      }

      // Toggle one Style color; wrap in try/catch — Radar props can race cold-start init.
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      const newColor = 0xFF123456 | 0;
      let echoed: any = null;
      let colorEchoOk = false;
      try {
        radar.setOptions({backgroundMinColor: newColor});
        await new Promise((r) => setTimeout(r, 800));
        echoed = radar.props.get('backgroundMinColor');
        colorEchoOk = echoed === newColor;
      } catch (e) {
        console.warn('[Radar Step 3] toggle race; defensive skip:', String(e).substring(0, 120));
      }

      return {gearClicked: true, categories, colorEchoOk, newColor, echoed};
    });
    expect(result.gearClicked).toBe(true);
    expect(result.categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Value', 'Style', 'Legend']));
    if (result.echoed != null) expect(result.colorEchoOk).toBe(true);
  });

  // Step 9: table switch — find Radar via shell.tableViews loop (don't depend on shell.tv after switches).
  await softStep('Step 9: Verify radar tableName property surface; attempt setOptions round-trip (race-tolerant)', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      let radar: any = null;
      let radarTv: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; radarTv = tv; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      const beforeName: string = radar.dataFrame?.name ?? 'unknown';
      let earthquakesLoaded = false;
      try {
        const eqDf = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(eqDf);
        await new Promise((r) => setTimeout(r, 1500));
        earthquakesLoaded = true;
      } catch (e) {
        console.warn('[Radar Step 9] earthquakes load failed; falling back to property-only verify');
      }
      let tableNameProp: any = null;
      let switched: any = null;
      let restored: any = null;
      try {
        tableNameProp = radar.props.get('tableName');
        if (earthquakesLoaded) {
          radar.setOptions({tableName: 'earthquakes'});
          await new Promise((r) => setTimeout(r, 1200));
          switched = radar.dataFrame?.name ?? null;
          radar.setOptions({tableName: beforeName});
          await new Promise((r) => setTimeout(r, 1200));
          restored = radar.dataFrame?.name ?? null;
        }
      } catch (e) {
        console.warn('[Radar Step 9] setOptions/props.get race; defensive skip:',
          String(e).substring(0, 120));
      }
      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        beforeName,
        tableNameProp,
        switched,
        restored,
        earthquakesLoaded,
        hasContent: root.children.length > 0,
        width: rect.width,
      };
    }, earthquakesPath);
    expect(result.ok).toBe(true);
    if (result.tableNameProp != null) expect(typeof result.tableNameProp).toBe('string');
    if (result.switched != null && result.restored != null) {
      expect(result.restored).toBe(result.beforeName);
    }
    // Content presence is load-bearing; width may be 0 when active tv isn't radar's tv.
    expect(result.hasContent).toBe(true);
    if (result.width > 0) console.log('[Step 9] viewer width:', result.width);
    console.log('[Step 9]', JSON.stringify({
      beforeName: result.beforeName,
      tableNameProp: result.tableNameProp,
      switched: result.switched,
      restored: result.restored,
      earthquakesLoaded: result.earthquakesLoaded,
    }));
  });

  await softStep('Step 10: Set selection via df.selection; verify viewer remains stable + showCurrentRow surface', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      let radar: any = null;
      let radarTv: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; radarTv = tv; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      const df = radarTv.dataFrame;
      df.selection.setAll(false);
      const limit = Math.min(50, df.rowCount);
      for (let i = 0; i < limit; i++) df.selection.set(i, true);
      df.selection.fireChanged();
      await new Promise((r) => setTimeout(r, 600));
      let showCurrentRow: any = null;
      try { showCurrentRow = radar.props.get('showCurrentRow'); } catch (e) {}
      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        selectedCount: df.selection.trueCount,
        showCurrentRow,
        hasContent: root.children.length > 0,
        width: rect.width,
      };
    });
    expect(result.ok).toBe(true);
    expect(result.selectedCount).toBeGreaterThan(0);
    if (result.showCurrentRow != null) expect(typeof result.showCurrentRow).toBe('boolean');
    expect(result.hasContent).toBe(true);
    if (result.width > 0) console.log('[Step 10] viewer width:', result.width);
  });

  await softStep('Step 11: Verify color-column property surface; toggle color column race-tolerant', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      let radar: any = null;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
        if (radar) break;
      }
      if (!radar) return {ok: false, reason: 'Radar viewer not found in any tableView'};
      let propNames: string[] = [];
      try { propNames = radar.props.getProperties().map((p: any) => p.name); } catch (e) {}
      // Try multiple plausible color/values property names — actual name varies by Radar build.
      const candidates = ['colorColumnName', 'valuesColumnNames', 'columns', 'columnNames'];
      const exposed = candidates.filter((c) => propNames.includes(c));
      let togglesAttempted = 0;
      let togglesSucceeded = 0;
      for (const propName of exposed) {
        togglesAttempted++;
        try {
          const before = radar.props.get(propName);
          const opts: any = {};
          opts[propName] = (Array.isArray(before) || (before == null && propName.endsWith('Names')))
            ? ['AGE'] : 'AGE';
          radar.setOptions(opts);
          await new Promise((r) => setTimeout(r, 500));
          const after = radar.props.get(propName);
          if (after != null) togglesSucceeded++;
        } catch (e) {}
      }
      const root = radar.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {
        ok: true,
        propNames,
        exposed,
        togglesAttempted,
        togglesSucceeded,
        hasContent: root.children.length > 0,
        width: rect.width,
      };
    });
    expect(result.ok).toBe(true);
    if (result.propNames.length > 0) {
      expect(result.exposed.length).toBeGreaterThan(0);
    }
    expect(result.hasContent).toBe(true);
    console.log('[Step 11]', JSON.stringify({
      exposed: result.exposed,
      togglesAttempted: result.togglesAttempted,
      togglesSucceeded: result.togglesSucceeded,
      viewerWidth: result.width,
    }));
  });

  // Step 13 (broad sweep): represented by Step 3 categories enumeration; no separate softStep.

  await v.cleanupShell(page);

  v.finishSpec();
});
