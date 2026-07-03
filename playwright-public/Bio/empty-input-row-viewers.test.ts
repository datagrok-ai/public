/* ---
sub_features_covered: [bio.analyze.activity-cliffs, bio.analyze.activity-cliffs.top-menu, bio.search.diversity, bio.search.diversity.top-menu, bio.search.similarity, bio.search.similarity.top-menu, bio.viewers.diversity-search, bio.viewers.similarity-search]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
type ViewerCase = {
  label: string;
  groupSelector: string;
  leafSelector: string;
  hasEditorDialog: boolean;
  viewerSelector: string | null;
  subFeatures: string[];
};
const viewerCases: ViewerCase[] = [
  {
    label: 'Sequence Similarity Search',
    groupSelector: '[name="div-Bio---Search"]',
    leafSelector: '[name="div-Bio---Search---Similarity-Search"]',
    hasEditorDialog: false,
    viewerSelector: '[name="viewer-Sequence-Similarity-Search"]',
    subFeatures: [
      'bio.viewers.similarity-search',
      'bio.search.similarity',
      'bio.search.similarity.top-menu',
    ],
  },
  {
    label: 'Sequence Diversity Search',
    groupSelector: '[name="div-Bio---Search"]',
    leafSelector: '[name="div-Bio---Search---Diversity-Search"]',
    hasEditorDialog: false,
    viewerSelector: '[name="viewer-Sequence-Diversity-Search"]',
    subFeatures: [
      'bio.viewers.diversity-search',
      'bio.search.diversity',
      'bio.search.diversity.top-menu',
    ],
  },
  {
    label: 'Activity Cliffs',
    groupSelector: '[name="div-Bio---Analyze"]',
    leafSelector: '[name="div-Bio---Analyze---Activity-Cliffs..."]',
    hasEditorDialog: true,
    viewerSelector: null,
    subFeatures: [
      'bio.analyze.activity-cliffs',
      'bio.analyze.activity-cliffs.top-menu',
    ],
  },
];
for (const vc of viewerCases) {
  test(`Bio ${vc.label} rejects empty current-row input with balloon`, async ({page}) => {
    test.setTimeout(180_000);
    test.fail(true, 'GROK-16111: Bio current-row viewers silently KNN on empty input, no rejection balloon — remove test.fail when fixed');
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_FASTA.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => resolve(), 4000);
      });
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
      if (hasMacromolecule) {
        for (let i = 0; i < 75; i++) {
          const canvas = document.querySelector('[name="viewer-Grid"] canvas') as HTMLCanvasElement;
          if (canvas && canvas.width > 0) break;
          await new Promise((r) => setTimeout(r, 200));
        }
      }
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (let i = 0; i < 30; i++) {
        for (const fn of probes) {
          try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
        }
        await new Promise((r) => setTimeout(r, 500));
      }
    });
    await softStep('Setup: open filter_FASTA, empty row 0, hook balloon, set current row', async () => {
      const setup = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
        if (!macro) throw new Error('No Macromolecule column detected on filter_FASTA.csv');
        if (df.rowCount < 2) throw new Error(`filter_FASTA.csv must have >=2 rows; got ${df.rowCount}`);
        macro.set(0, '');
        df.currentRowIdx = 0;
        (window as any).__balloonCalls = [];
        const origWarn = grok.shell.warning.bind(grok.shell);
        const origErr = grok.shell.error.bind(grok.shell);
        grok.shell.warning = ((msg: any, opts?: any) => {
          try { (window as any).__balloonCalls.push({type: 'warning', msg: String(msg).slice(0, 300)}); } catch { /* ignore */ }
          return origWarn(msg, opts);
        }) as any;
        grok.shell.error = ((msg: any, opts?: any) => {
          try { (window as any).__balloonCalls.push({type: 'error', msg: String(msg).slice(0, 300)}); } catch { /* ignore */ }
          return origErr(msg, opts);
        }) as any;
        return {
          rowCount: df.rowCount,
          macroName: macro.name,
          row0value: macro.get(0),
          row0empty: macro.get(0) === '' || macro.get(0) == null,
          currentRowIdx: df.currentRowIdx,
        };
      });
      expect(setup.rowCount).toBeGreaterThanOrEqual(2);
      expect(setup.macroName).not.toBeNull();
      expect(setup.row0empty).toBe(true);
      expect(setup.currentRowIdx).toBe(0);
    });
    const baseRowCount: number = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    await softStep(`${vc.label}: open Bio > ${vc.groupSelector.includes('Search') ? 'Search' : 'Analyze'} > ${vc.label}${vc.hasEditorDialog ? ' (defaults OK)' : ''}`, async () => {
      await page.evaluate(async ({groupSel, leafSel}) => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        let group: Element | null = null;
        for (let i = 0; i < 40; i++) {
          group = document.querySelector(groupSel);
          if (group) break;
          await new Promise((r) => setTimeout(r, 100));
        }
        if (!group) throw new Error(`Bio group ${groupSel} not found in top menu`);
        group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        let leaf: Element | null = null;
        for (let i = 0; i < 40; i++) {
          leaf = document.querySelector(leafSel);
          if (leaf) break;
          await new Promise((r) => setTimeout(r, 100));
        }
        if (!leaf) throw new Error(`Bio leaf ${leafSel} not found under ${groupSel}`);
        (leaf as HTMLElement).click();
      }, {groupSel: vc.groupSelector, leafSel: vc.leafSelector});
      if (vc.hasEditorDialog) {
        await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
        await page.locator('.d4-dialog [name="button-OK"]').click();
      }
    });
    await softStep(`${vc.label}: balloon surfaces, no silent zero-row result`, async () => {
      await page.waitForFunction(() => {
        return Array.isArray((window as any).__balloonCalls)
          && (window as any).__balloonCalls.length > 0;
      }, null, {timeout: 30_000}).catch(() => { /* bounded readiness gate; the balloon assertion below is the real check */ });
      const probe = await page.evaluate((viewerSel) => {
        const df = grok.shell.tv.dataFrame;
        const calls = ((window as any).__balloonCalls || []) as Array<{type: string; msg: string}>;
        const viewerTypes = Array.from((grok.shell.tv as any).viewers).map((v: any) => v.type);
        const docked = viewerSel
          ? !!document.querySelector(viewerSel)
          : viewerTypes.includes('Scatter plot');
        return {
          balloonCount: calls.length,
          balloonTypes: calls.map((c) => c.type),
          balloonMsgs: calls.map((c) => c.msg),
          rowCount: df.rowCount,
          viewerTypes,
          docked,
        };
      }, vc.viewerSelector);
      // Invariant 2 (.md): no silent zero-row result — the source table is not rewritten and the viewer reacts (docks or rejects), never a silent no-op.
      expect(probe.rowCount, `${vc.label}: source table must not be silently rewritten on empty input`).toBe(baseRowCount);
      expect(probe.docked || probe.balloonCount > 0, `${vc.label}: viewer must react on empty input (dock or reject), not silently no-op`).toBe(true);
      // Invariant 1 (.md): empty current-row input must surface a rejection balloon. GROK-16111 keeps this failing today; test.fail() tracks the xfail and auto-flips to red when fixed.
      expect(probe.balloonCount, 'GROK-16111: empty current-row input must surface a rejection balloon').toBeGreaterThan(0);
    });
    finishSpec();
  });
}
