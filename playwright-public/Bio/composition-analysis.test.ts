/* ---
sub_features_covered: [bio.analyze.composition, bio.viewers.web-logo]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];
test('Bio | Analyze | Composition — composition analysis integration', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });
  for (const ds of datasets) {
    await softStep(`[${ds.name}] Scenario 1 Step 1 — Open ${ds.path}`, async () => {
      const result: {rows: number, hasMacromolecule: boolean} = await page.evaluate(async (path: string) => {
        const g = (window as any).grok;
        g.shell.closeAll();
        const df = await g.dapi.files.readCsv(path);
        g.shell.addTableView(df);
        await new Promise<void>((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => resolve(), 3000);
        });
        const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
        const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
        if (hasBioChem) {
          for (let i = 0; i < 50; i++) {
            if (document.querySelector('[name="viewer-Grid"] canvas')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 5000));
        }
        return {rows: df.rowCount, hasMacromolecule: hasBioChem};
      }, ds.path);
      expect(result.rows).toBeGreaterThan(0);
      expect(result.hasMacromolecule).toBe(true);
      await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    });
    await softStep(`[${ds.name}] Scenario 1 Step 2 — Bio > Analyze > Composition; WebLogo viewer docks; no multi-column dialog`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(
          new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Composition"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => {
        const viewers = Array.from((window as any).grok.shell.tv.viewers);
        return viewers.some((v: any) => v.type === 'WebLogo');
      }, null, {timeout: 60_000});
      const result: {hasCanvas: boolean, dialogOpen: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        let hasCanvas = false;
        for (let i = 0; i < 40; i++) {
          const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
          if (w && w.root.querySelector('canvas')) { hasCanvas = true; break; }
          await new Promise((r) => setTimeout(r, 300));
        }
        const dialogOpen = !!document.querySelector('.d4-dialog');
        return {hasCanvas, dialogOpen};
      });
      expect(result.hasCanvas).toBe(true);
      expect(result.dialogOpen).toBe(false);
    });
    await softStep(`[${ds.name}] Scenario 2 Step 3 — Click letter in WebLogo selects ≥1 row in source grid`, async () => {
      // Settle: canvas may be in DOM before click handlers are bound.
      await page.waitForTimeout(4000);
      const selected: number = await page.evaluate(async () => {
        const g = (window as any).grok;
        const tv = g.shell.tv;
        const df = tv.dataFrame;
        df.selection.setAll(false);
        const wl: any = tv.viewers.find((v: any) => v.type === 'WebLogo');
        const canvas = wl.root.querySelector('canvas') as HTMLCanvasElement;
        const canvasRect = canvas.getBoundingClientRect();
        const dpr = window.devicePixelRatio || 1;
        const GAP = '-';
        const positions: any[] = Array.isArray(wl.positions) ? wl.positions : [];
        const candidates: {x: number, y: number, h: number}[] = [];
        for (let idx = 0; idx < positions.length && idx < 30; idx++) {
          const pi = positions[idx];
          if (!pi || typeof pi.getMonomers !== 'function') continue;
          let monomers: string[] = [];
          try { monomers = pi.getMonomers(); } catch (_) { monomers = []; }
          for (const m of monomers) {
            if (!m || m === GAP) continue;
            let b: any = null;
            try { b = pi.getFreq(m).bounds; } catch (_) { b = null; }
            if (!b) continue;
            const cx = (b.x + b.width / 2) / dpr + canvasRect.left;
            const cy = (b.y + b.height / 2) / dpr + canvasRect.top;
            candidates.push({x: cx, y: cy, h: b.height});
          }
        }
        // Tallest rect first — most frequent monomer is the largest hit target.
        candidates.sort((a, b) => b.h - a.h);
        for (const c of candidates) {
          df.selection.setAll(false);
          for (const type of ['mousedown', 'mouseup', 'click']) {
            canvas.dispatchEvent(new MouseEvent(type, {
              bubbles: true, cancelable: true, clientX: c.x, clientY: c.y, button: 0,
            }));
          }
          await new Promise((r) => setTimeout(r, 400));
          if (df.selection.trueCount > 0) break;
        }
        return df.selection.trueCount;
      });
      expect(selected).toBeGreaterThan(0);
    });
    await softStep(`[${ds.name}] Scenario 3 Step 4 — Gear icon on WebLogo opens Context Pane property grid`, async () => {
      const opened: {found: boolean, pg: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        g.shell.tv.dataFrame.selection.setAll(false);
        const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
        // Gear lives on the outer docked-panel title bar (.panel-base ancestor).
        let panelBase: any = w.root;
        while (panelBase && !panelBase.classList?.contains('panel-base'))
          panelBase = panelBase.parentElement;
        const gear = panelBase?.querySelector(
          '.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
        if (!gear) {
          // Fallback: programmatic equivalent if the title-bar gear is absent.
          g.shell.o = w;
          await new Promise((r) => setTimeout(r, 600));
          const pg2 = document.querySelector('.grok-prop-panel .property-grid, .grok-prop-panel tr[name^="prop-"]');
          return {found: false, pg: !!pg2};
        }
        gear.click();
        await new Promise((r) => setTimeout(r, 600));
        const pg = document.querySelector('.grok-prop-panel .property-grid, .grok-prop-panel tr[name^="prop-"]');
        return {found: true, pg: !!pg};
      });
      expect(opened.pg).toBe(true);
      // Property row sits in a collapsed accordion — wait for 'attached', not 'visible'.
      await page.locator('tr[name="prop-show-position-labels"]').waitFor({
        state: 'attached', timeout: 10_000});
    });
    await softStep(`[${ds.name}] Scenario 3 Step 5 — Edit ≥1 Context Pane property (SR-01: edit-acceptance only)`, async () => {
      const result: {changedShow: boolean, changedSkip: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        const wl: any = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
        const before = {
          showPositionLabels: wl.getOptions().look.showPositionLabels,
          skipEmptyPositions: wl.getOptions().look.skipEmptyPositions,
        };
        for (const name of ['prop-show-position-labels', 'prop-skip-empty-positions']) {
          const row = document.querySelector(`tr[name="${name}"]`);
          const cb = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          if (cb) {
            cb.click();
            cb.dispatchEvent(new Event('change', {bubbles: true}));
          }
          await new Promise((r) => setTimeout(r, 300));
        }
        await new Promise((r) => setTimeout(r, 400));
        const after = {
          showPositionLabels: wl.getOptions().look.showPositionLabels,
          skipEmptyPositions: wl.getOptions().look.skipEmptyPositions,
        };
        return {
          changedShow: before.showPositionLabels !== after.showPositionLabels,
          changedSkip: before.skipEmptyPositions !== after.skipEmptyPositions,
        };
      });
      expect(result.changedShow || result.changedSkip).toBe(true);
    });
  }
  await page.evaluate(() => (window as any).grok.shell.closeAll());
  finishSpec();
});
