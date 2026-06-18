import {test, expect} from '@playwright/test';
import {loginAndOpenFile, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

// Paired scenario: helm-bio-menu-integration.md (cross-feature smoke).
// Selector sources (grok-browser/references):
//   - .claude/skills/grok-browser/references/bio.md   (§ "Top-menu Bio entries" — div-Bio--- vectors)
//   - .claude/skills/grok-browser/references/helm.md   (HELM column tags / renderer)

test.use(specTestOptions);
test.use({timeout: 600_000});

// Dot-namespaced for the direct file-browse URL (see loginAndOpenFile).
const DATASET_PATH = 'System.AppData/Helm/samples/helm-showcase.csv';

// Bio errors that indicate a Helm breakage (the only failures this cross-feature
// smoke owns). Generic Bio compute artifacts (e.g. UMAP "Dimensionality
// reduction failed" on the tiny showcase) are tolerated and logged.
const HELM_ERR = /helm|jsdraw|pseudo-?molfile|getmolfiles|cell\s*render|monomer.{0,24}(render|null|fail)/i;

// HELM-applicable Bio leaves, each with the concrete OUTCOME to wait for
// (verified live on the showcase). `expect`:
//   'viewer' — a viewer of `type` MUST dock
//   'cols'   — colCount MUST grow by >= `min`
//   'table'  — a new result table MUST open (MSA)
//   'run'    — runs cleanly but produces no deterministic artifact with defaults
//              on this dataset (still must raise no Helm error)
// Skipped leaves are documented in the paired .md.
type BioFn = {leaf: string, label: string} &
  ({expect: 'viewer', type: string} | {expect: 'cols', min: number} | {expect: 'table'} | {expect: 'run'});

const BIO_FUNCS: BioFn[] = [
  {leaf: 'div-Bio---Analyze---Composition', label: 'Composition', expect: 'viewer', type: 'WebLogo'},
  {leaf: 'div-Bio---Analyze---Sequence-Space...', label: 'Sequence Space', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Analyze---Activity-Cliffs...', label: 'Activity Cliffs', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Analyze---MSA...', label: 'MSA', expect: 'table'},
  {leaf: 'div-Bio---Analyze---Hierarchical-Clustering...', label: 'Hierarchical Clustering', expect: 'run'},
  {leaf: 'div-Bio---Transform---Convert-Sequence-Notation...', label: 'Convert Sequence Notation', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Transform---Split-to-Monomers...', label: 'Split to Monomers', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Calculate---Extract-Region...', label: 'Extract Region', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Annotate---Apply-Numbering-Scheme...', label: 'Apply Numbering Scheme', expect: 'run'},
  {leaf: 'div-Bio---Annotate---Scan-Liabilities...', label: 'Scan Liabilities', expect: 'cols', min: 1},
  {leaf: 'div-Bio---Annotate---Manage-Annotations...', label: 'Manage Annotations', expect: 'run'},
  {leaf: 'div-Bio---Search---Similarity-Search', label: 'Similarity Search', expect: 'viewer', type: 'Sequence Similarity Search'},
  {leaf: 'div-Bio---Search---Diversity-Search', label: 'Diversity Search', expect: 'viewer', type: 'Sequence Diversity Search'},
  {leaf: 'div-Bio---Search---Subsequence-Search-...', label: 'Subsequence Search', expect: 'viewer', type: 'Filters'},
  {leaf: 'div-Bio---PolyTool---Convert...', label: 'PolyTool Convert', expect: 'run'},
  {leaf: 'div-Bio---PolyTool---Enumerate-HELM...', label: 'PolyTool Enumerate HELM', expect: 'run'},
];

test('Helm — Bio menu cross-feature integration on a HELM column', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Open the showcase DIRECTLY via its instance-derived file URL.
  await loginAndOpenFile(page, DATASET_PATH);
  await page.evaluate(() => {
    const g = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.windows.simpleMode = true;
    g.shell.windows.showContextPanel = true;
  });
  await page.waitForFunction(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    return df?.col('HELM')?.semType === 'Macromolecule';
  }, null, {timeout: 45_000});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForTimeout(1500);

  // Baseline + capture the showcase view so we can return to it after a leaf
  // (MSA) opens its own result table. Make HELM the current column so the Bio
  // top menu attaches.
  const base = await page.evaluate(() => {
    const g = (window as any).grok;
    (window as any).__sc = g.shell.tv;             // showcase TableView handle
    const df = g.shell.tv.dataFrame;
    df.currentCol = df.col('HELM');
    return {
      cols: df.columns.length,
      viewers: Array.from(g.shell.tv.viewers).map((v: any) => v.type),
      tables: Array.from(g.shell.tables).length,
    };
  });

  const bioPresent = await page.evaluate(() => !!document.querySelector('[name="div-Bio"]'));
  expect(bioPresent, 'Bio top menu MUST attach for a Macromolecule HELM column').toBe(true);

  // Open a Bio leaf via the dispatchEvent vector (per bio.md; works in overflow).
  const openBioLeaf = async (leaf: string): Promise<boolean> => {
    await page.evaluate(() => document.body.dispatchEvent(new MouseEvent('click', {bubbles: true})));
    await page.waitForTimeout(150);
    await page.evaluate(() => document.querySelector('[name="div-Bio"]')
      ?.dispatchEvent(new MouseEvent('click', {bubbles: true})));
    await page.waitForTimeout(400);
    await page.evaluate((l) => document.querySelector(`[name="${l.split('---').slice(0, 2).join('---')}"]`)
      ?.dispatchEvent(new MouseEvent('mouseover', {bubbles: true})), leaf);
    await page.waitForTimeout(400);
    return page.evaluate((l) => {
      const el = document.querySelector(`[name="${l}"]`);
      if (!el) return false;
      el.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      return true;
    }, leaf);
  };

  // Read the live state used for outcome polling + Helm-error detection.
  const readState = () => page.evaluate(() => {
    const g = (window as any).grok;
    const tv = g.shell.tv;
    return {
      cols: tv?.dataFrame?.columns?.length ?? 0,
      viewers: tv ? Array.from(tv.viewers).map((v: any) => v.type) : [],
      tables: Array.from(g.shell.tables).length,
      balloons: Array.from(document.querySelectorAll('.d4-balloon'))
        .filter((b: any) => /error/i.test(b.className)).map((b: any) => (b.textContent || '').trim()),
    };
  });

  // Return to the showcase view, cancel dialogs, close added viewers, drop added
  // columns, re-select HELM.
  const cleanup = async () => {
    await page.evaluate((baseViewers: string[]) => {
      const g = (window as any).grok;
      document.querySelectorAll('.d4-dialog [name="button-CANCEL"]').forEach((b: any) => { try { b.click(); } catch (_) {/* */} });
      // MSA (and any leaf opening its own table) switches the active view — go back.
      if ((window as any).__sc && g.shell.v !== (window as any).__sc) {
        try { g.shell.v = (window as any).__sc; } catch (_) {/* */}
      }
      const tv = g.shell.tv;
      if (tv) {
        for (const v of Array.from(tv.viewers) as any[])
          if (v.type !== 'Grid' && !baseViewers.includes(v.type)) { try { v.close(); } catch (_) {/* */} }
      }
    }, base.viewers);
    await page.waitForTimeout(300);
    await page.evaluate((baseCols: number) => {
      const g = (window as any).grok;
      const df = g.shell.tv?.dataFrame;
      if (df) { while (df.columns.length > baseCols) df.columns.remove(df.columns.byIndex(df.columns.length - 1).name); }
      if (df?.col('HELM')) df.currentCol = df.col('HELM');
    }, base.cols);
    await page.waitForTimeout(200);
  };

  for (const fn of BIO_FUNCS) {
    const expectDesc = fn.expect === 'viewer' ? `docks the "${fn.type}" viewer`
      : fn.expect === 'cols' ? `adds >= ${fn.min} column(s)`
        : fn.expect === 'table' ? 'opens a result table' : 'runs cleanly';
    await softStep(`Bio | ${fn.label}: ${expectDesc}; no Helm-related error`, async () => {
      const opened = await openBioLeaf(fn.leaf);
      expect(opened, `Bio leaf "${fn.leaf}" MUST be locatable in the menu`).toBe(true);
      await page.waitForTimeout(2000);

      // Accept defaults if a dialog with an enabled OK appeared.
      if (await page.locator('.d4-dialog').count() > 0) {
        await page.evaluate(() => {
          const ok = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement | null;
          if (ok && !ok.classList.contains('disabled')) ok.click();
        });
      }

      // Poll for the concrete outcome (and capture errors). Viewers dock fast;
      // column/table production (converting/splitting 55 long sequences) is slow
      // and async, so give it a generous window.
      const maxPolls = fn.expect === 'viewer' ? 24 : 60; // 12s vs 30s
      let met = fn.expect === 'run';
      let last = await readState();
      for (let i = 0; i < maxPolls; i++) {
        last = await readState();
        if (fn.expect === 'viewer' && last.viewers.includes(fn.type)) met = true;
        else if (fn.expect === 'cols' && last.cols >= base.cols + fn.min) met = true;
        else if (fn.expect === 'table' && last.tables > base.tables) met = true;
        if (met) break;
        await page.waitForTimeout(500);
      }
      if (fn.expect === 'cols')
        console.warn(`Bio | ${fn.label}: column delta = ${last.cols - base.cols} (base=${base.cols}, observed=${last.cols}).`);

      // Positive outcome assertion (where deterministic).
      if (fn.expect === 'viewer')
        expect(met, `Bio | ${fn.label} MUST dock the "${fn.type}" viewer (observed: ${JSON.stringify(last.viewers)})`).toBe(true);
      else if (fn.expect === 'cols')
        expect(met, `Bio | ${fn.label} MUST add >= ${fn.min} column(s) (base=${base.cols}, observed=${last.cols})`).toBe(true);
      else if (fn.expect === 'table')
        expect(met, `Bio | ${fn.label} MUST open a result table (base=${base.tables}, observed=${last.tables})`).toBe(true);

      // Helm-integration assertion: Bio consuming the HELM column must not break Helm.
      const helmErrs = last.balloons.filter((t) => HELM_ERR.test(t));
      const otherErrs = last.balloons.filter((t) => !HELM_ERR.test(t));
      if (otherErrs.length > 0)
        console.warn(`Bio | ${fn.label}: tolerated non-Helm balloon(s): ${JSON.stringify(otherErrs).slice(0, 200)}`);
      expect(helmErrs.length,
        `Bio | ${fn.label} MUST NOT raise a Helm-related error balloon: ${JSON.stringify(helmErrs)}`).toBe(0);

      await cleanup();
    });
  }

  await softStep('Post-sweep: the HELM column + renderer survived all Bio operations', async () => {
    const integ = await page.evaluate(() => {
      const g = (window as any).grok;
      if ((window as any).__sc && g.shell.v !== (window as any).__sc) { try { g.shell.v = (window as any).__sc; } catch (_) {/* */} }
      const df = g.shell.tv?.dataFrame;
      const c = df?.col('HELM');
      return {
        semType: c?.semType ?? null,
        renderer: c?.getTag ? c.getTag('cell.renderer') : null,
        rowCount: df?.rowCount ?? 0,
        canvas: !!document.querySelector('[name="viewer-Grid"] canvas'),
      };
    });
    expect(integ.semType, 'post-sweep: HELM column MUST remain semType=Macromolecule').toBe('Macromolecule');
    expect(integ.renderer, 'post-sweep: HELM column MUST remain cell.renderer=helm').toBe('helm');
    expect(integ.rowCount, 'post-sweep: showcase row count MUST be intact (55)').toBe(55);
    expect(integ.canvas, 'post-sweep: the grid MUST still render').toBe(true);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());
  finishSpec();
});
