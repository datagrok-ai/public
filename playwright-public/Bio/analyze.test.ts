/* ---
sub_features_covered: [bio.analyze.activity-cliffs, bio.analyze.activity-cliffs.editor, bio.analyze.activity-cliffs.top-menu, bio.analyze.composition, bio.analyze.sequence-space, bio.analyze.sequence-space.editor, bio.analyze.sequence-space.top-menu]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as bio from '../helpers/bio';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];
for (const ds of datasets) {
  test(`Bio Analyze umbrella on ${ds.name}`, async ({page}) => {
    test.setTimeout(180_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      const detected = () => Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some((c: any) => c.semType === 'Macromolecule');
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { if (detected()) { sub.unsubscribe(); resolve(); } });
        const deadline = Date.now() + 30_000;
        const poll = () => { if (detected() || Date.now() > deadline) { sub.unsubscribe(); resolve(); } else setTimeout(poll, 200); };
        poll();
      });
      if (detected()) {
        for (let i = 0; i < 60; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
      }
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      throw new Error('Bio package did not warm up: no probe function callable');
    });
    await page.evaluate(async () => {
      const candidates = [
        ['Bio:sequenceSpaceTopMenu', 'Bio:sequenceSpace'],
        ['Bio:activityCliffsTopMenu', 'Bio:activityCliffs', 'Bio:macromoleculeActivityCliffs'],
        ['Bio:compositionAnalysisWidget', 'Bio:composition', 'Bio:compositionAnalysis'],
      ];
      const findAny = (names: string[]): boolean => {
        for (const n of names) {
          try {
            if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
          } catch { /* try next */ }
        }
        return false;
      };
      const deadline = Date.now() + 15_000;
      while (Date.now() < deadline) {
        if (candidates.every(findAny)) return;
        await new Promise((r) => setTimeout(r, 300));
      }
      throw new Error('Bio analyze functions did not register within 15s');
    });
    await page.waitForFunction(() => !!(grok as any).functions.find?.('Bio:sequenceSpace'),
      null, {timeout: 15_000});
    await softStep(`${ds.name}: Bio > Analyze > Sequence Space — run with defaults`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Sequence-Space...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Sequence Space');
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      const baseScatter: number = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        ({base, baseSc}) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length > baseSc,
        {base: baseCols, baseSc: baseScatter}, {timeout: 60_000});
      const hasEmbedCols = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names()
        .map((s: string) => s.toLowerCase()).some((n: string) => n.includes('embed') || n.startsWith('sequence space')));
      expect(hasEmbedCols).toBe(true);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    await softStep(`${ds.name}: Bio > Analyze > Activity Cliffs — run with defaults`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Activity-Cliffs...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Activity Cliffs');
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      const baseScatter: number = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        ({base, baseSc}) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length > baseSc,
        {base: baseCols, baseSc: baseScatter}, {timeout: 60_000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    await softStep(`${ds.name}: Bio > Analyze > Composition — WebLogo docks (no dialog)`, async () => {
      // Composition has no "..." suffix — opens directly with no dialog; docks a WebLogo viewer.
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Composition');
      await page.waitForFunction(
        () => Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'WebLogo'),
        null, {timeout: 60_000});
      const hasWebLogo = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'WebLogo'));
      expect(hasWebLogo).toBe(true);
    });
    // Scenario 3 — Composition Gear → Context Panel wiring (FASTA only).
    if (ds.name === 'FASTA') {
      await softStep(`${ds.name}: Composition Gear → Context Panel binds WebLogo property surface`, async () => {
        // WebLogo title bar surfaces no Gear under body.selenium; grok.shell.o opens the same Context Panel.
        // The property grid mounts asynchronously after the object is selected, so poll for it instead
        // of a single fixed-delay read (an 800ms snapshot races the panel build under load).
        const result: {hasPropPanel: boolean, hasPropertyGrid: boolean, selectedIsWebLogo: boolean, panelText: string} =
          await page.evaluate(async () => {
          // Force the Context Panel open — under simpleMode it can stay collapsed,
          // so selecting the object alone would not mount the property grid.
          try { (grok.shell as any).windows.showContextPanel = true; } catch (e) {}
          let propPanel: Element | null = null;
          let propertyGrid: Element | null = null;
          let panelText = '';
          for (let i = 0; i < 40; i++) {
            // Re-select each round: the WebLogo viewer ref and the panel build can
            // race, and re-assigning grok.shell.o re-triggers the property surface.
            const wl = (grok.shell.tv as any).viewers.find((v: any) => v.type === 'WebLogo');
            if (wl) (grok as any).shell.o = wl;
            await new Promise((r) => setTimeout(r, 500));
            propPanel = document.querySelector('.grok-prop-panel');
            propertyGrid = document.querySelector('.grok-prop-panel .property-grid');
            panelText = propPanel?.textContent ?? '';
            if (propPanel && propertyGrid && /Data|Layout|Behavior|Style/.test(panelText)) break;
          }
          return {hasPropPanel: !!propPanel, hasPropertyGrid: !!propertyGrid,
            selectedIsWebLogo: (grok as any).shell.o?.type === 'WebLogo', panelText};
        });
        expect(result.hasPropPanel).toBe(true);
        expect(result.hasPropertyGrid).toBe(true);
        // Bind-check: the mounted surface belongs to the WebLogo viewer and shows its property sections.
        expect(result.selectedIsWebLogo).toBe(true);
        expect(result.panelText).toMatch(/Data|Layout|Behavior|Style/);
      });
    }
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    finishSpec();
  });
}
