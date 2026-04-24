import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];

test('Composition Analysis manual test', async ({page}) => {
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
    await softStep(`[${ds.name}] Open ${ds.path}`, async () => {
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

    await softStep(`[${ds.name}] Open Bio > Analyze > Composition — WebLogo viewer opens`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Composition"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => {
        const viewers = Array.from((window as any).grok.shell.tv.viewers);
        return viewers.some((v: any) => v.type === 'WebLogo');
      }, null, {timeout: 60_000});
      const hasCanvas: boolean = await page.evaluate(async () => {
        const g = (window as any).grok;
        for (let i = 0; i < 40; i++) {
          const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
          if (w && w.root.querySelector('canvas')) return true;
          await new Promise((r) => setTimeout(r, 300));
        }
        return false;
      });
      expect(hasCanvas).toBe(true);
    });

    await softStep(`[${ds.name}] Click a letter in the WebLogo — rows are selected`, async () => {
      // Give WebLogo extra time to finish rendering and wire event handlers
      // (subsequent datasets fail without this settle — first render works fast, later ones need ~3s).
      await page.waitForTimeout(3000);
      const selected: number = await page.evaluate(async () => {
        const g = (window as any).grok;
        const tv = g.shell.tv;
        const df = tv.dataFrame;
        df.selection.setAll(false);
        const w = tv.viewers.find((v: any) => v.type === 'WebLogo');
        const canvas = w.root.querySelector('canvas') as HTMLCanvasElement;
        const rect = canvas.getBoundingClientRect();
        // Probe a few x-offsets so narrower letter columns (HELM/MSA) still get a hit.
        const offsets = [18, 30, 50, 80, 120, 160];
        for (const dx of offsets) {
          const x = rect.x + dx;
          const y = rect.y + rect.height / 2;
          for (const type of ['mousedown', 'mouseup', 'click']) {
            canvas.dispatchEvent(new MouseEvent(type, {
              bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0,
            }));
          }
          await new Promise((r) => setTimeout(r, 400));
          if (df.selection.trueCount > 0) break;
        }
        return df.selection.trueCount;
      });
      expect(selected).toBeGreaterThan(0);
    });

    await softStep(`[${ds.name}] Click Gear on WebLogo title bar — property grid appears`, async () => {
      const opened: {found: boolean, pg: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        g.shell.tv.dataFrame.selection.setAll(false);
        const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
        let panelBase: any = w.root;
        while (panelBase && !panelBase.classList?.contains('panel-base')) panelBase = panelBase.parentElement;
        const gear = panelBase?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
        if (!gear) return {found: false, pg: false};
        gear.click();
        await new Promise((r) => setTimeout(r, 600));
        const pg = document.querySelector('.grok-prop-panel .property-grid');
        return {found: true, pg: !!pg};
      });
      expect(opened.found).toBe(true);
      expect(opened.pg).toBe(true);
      // Property row may be in a collapsed accordion section — only check it is attached, not visible.
      await page.locator('tr[name="prop-show-position-labels"]').waitFor({state: 'attached', timeout: 10_000});
    });

    await softStep(`[${ds.name}] Change arbitrary properties in the Context Pane`, async () => {
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

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
