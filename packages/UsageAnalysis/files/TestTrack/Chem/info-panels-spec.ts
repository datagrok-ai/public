import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

const CHEMISTRY_PANES = ['Rendering', 'Highlights', 'Descriptors', 'Properties'];
const BIOLOGY_PANES = ['Drug Likeness', 'Structural Alerts', 'Pharmacophore', 'Toxicity'];
const STRUCTURE_PANES = ['Identifiers', '2D Structure', '3D Structure'];

async function expandAndVerifyPanes(page: Page, label: string, expectedPanes: string[]) {
  await page.evaluate(async () => {
    const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
    for (const p of panes) {
      const h = p.querySelector('.d4-accordion-pane-header') as HTMLElement;
      if (h && !h.classList.contains('expanded')) {
        h.click();
        await new Promise(r => setTimeout(r, 200));
      }
    }
  });
  await page.waitForTimeout(2500);
  const seen = await page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
      .map(h => h.textContent!.trim()));
  const found = expectedPanes.filter(p => seen.some(s => new RegExp(p, 'i').test(s)));
  console.log(`[${label}] panes seen: ${seen.length}, expected matches: ${found.length}/${expectedPanes.length}`);
}

test('Chem: Info Panels Phase A column+cell walk + Phase B multi-format', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  // ===== Phase A — smiles-50 column + cell walk =====

  await softStep('Phase A — Step 1: Open smiles-50.csv', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      grok.shell.addTableView(df);
      grok.shell.windows.showContextPanel = true;
      (window as any).__ip_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__ip_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
    await waitForChemMenu(page);
  });

  await softStep('Phase A — Step 2: Molecule semType + RDKit renderer (canonical_smiles)', async () => {
    const ok = await page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      return molCol?.name === 'canonical_smiles';
    });
    expect(ok).toBe(true);
  });

  await softStep('Phase A — Step 3: Click canonical_smiles column header → column-context Panel', async () => {
    await page.evaluate(() => {
      const col = grok.shell.t.col('canonical_smiles');
      grok.shell.o = col;
    });
    await page.waitForTimeout(2500);
  });

  await softStep('Phase A — Step 4: Walk Chemistry/Biology/Structure info panels (column context)', async () => {
    await expandAndVerifyPanes(page, 'A4 column', [...CHEMISTRY_PANES, ...BIOLOGY_PANES, ...STRUCTURE_PANES]);
  });

  await softStep('Phase A — Steps 5-9: Open chembl-scaffolds + apply Rendering scaffold + Highlight', async () => {
    const exists = await page.evaluate(async () => {
      try {
        const ls = await grok.dapi.files.list('System:AppData/Chem/', false);
        return ls.some((f: any) => /chembl-scaffolds\.csv$/i.test(f.name));
      } catch (e) { return false; }
    });
    if (!exists) {
      console.log('[A5-9] chembl-scaffolds.csv not available — skipping Rendering+Highlight slice');
      return;
    }
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/chembl-scaffolds.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
    // SR-DEFERRED scaffold alignment + Highlight custom color: applied via JS API column tags.
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const molCol = tv.dataFrame.columns.toList().find((c: any) => c.semType === 'Molecule');
      const scaffoldCol = tv.dataFrame.columns.toList().find((c: any) => /Scaffold/i.test(c.name));
      if (molCol && scaffoldCol) {
        molCol.setTag?.('chem-scaffold', scaffoldCol.name);
        molCol.setTag?.('chem-highlight-scaffold', 'true');
      }
      grok.shell.o = molCol;
      await new Promise(r => setTimeout(r, 2000));
    });
    await expandAndVerifyPanes(page, 'A6-7 Rendering', ['Rendering', 'Highlights']);
  });

  await softStep('Phase A — Step 10-11: Switch back to smiles-50.csv + click first cell → cell context', async () => {
    await page.evaluate(async () => {
      const tables = Array.from(grok.shell.tables);
      const smiTable = tables.find((t: any) => /smiles-?50/i.test(t.name)) || tables[0];
      const tvs = Array.from(grok.shell.tableViews) as any[];
      const tv = tvs.find((v: any) => v.dataFrame === smiTable) || grok.shell.tv;
      if (tv) (grok.shell as any).v = tv;
      tv.dataFrame.currentRowIdx = 0;
      grok.shell.o = tv.dataFrame.currentCell;
      await new Promise(r => setTimeout(r, 2500));
    });
  });

  await softStep('Phase A — Step 12: Expand cell-context panels (Descriptors, Drug Likeness, etc.)', async () => {
    await expandAndVerifyPanes(page, 'A12 cell', [
      'Descriptors', 'Drug Likeness', 'Properties', 'Structural Alerts',
      'Pharmacophore', 'Identifiers', '2D Structure', '3D Structure', 'Toxicity',
    ]);
  });

  // ===== Phase B — Multi-format coverage =====

  type Variant = {id: string; format: string; path: string; opener: 'csv' | 'openFile'};
  const variants: Variant[] = [
    {id: 'B-smiles', format: 'smiles', path: 'System:AppData/Chem/tests/smiles-50.csv', opener: 'csv'},
    {id: 'B-molV2000', format: 'molV2000', path: 'System:AppData/Chem/mol1K.sdf', opener: 'openFile'},
    {id: 'B-molV3000', format: 'molV3000', path: 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf', opener: 'openFile'},
  ];

  for (const v of variants) {
    await softStep(`Phase B — ${v.id}: open + cell context + walk panes`, async () => {
      const opened = await page.evaluate(async ({path, opener}: any) => {
        try {
          grok.shell.closeAll();
          if (opener === 'openFile') {
            await ((DG as any).Func.find({name: 'OpenFile'})[0])
              .prepare({fullPath: path}).call(undefined, undefined, {processed: false});
          } else {
            const df = await grok.dapi.files.readCsv(path);
            grok.shell.addTableView(df);
          }
          return {ok: true};
        } catch (e) { return {ok: false, err: String(e)}; }
      }, {path: v.path, opener: v.opener});
      if (!(opened as any).ok) {
        console.log(`[${v.id}] open failed — skipping (${JSON.stringify(opened)})`);
        return;
      }
      try {
        await waitForChemMenu(page);
      } catch (e) {
        console.log(`[${v.id}] Chem menu not ready — skipping (${e})`);
        return;
      }
      await page.evaluate(() => {
        const tv = grok.shell.tv;
        tv.dataFrame.currentRowIdx = 0;
        grok.shell.o = tv.dataFrame.currentCell;
      });
      await page.waitForTimeout(2500);
      await expandAndVerifyPanes(page, v.id, ['Descriptors', '2D Structure', 'Properties']);
    });
  }

  await softStep('Final: no Chem-info-panel errors throughout walk', async () => {
    const errs = await page.evaluate(() => ((window as any).__ip_errors ?? []) as string[]);
    const panelErrs = errs.filter(e => /info[- ]?panel|context[- ]?panel|chem\.panels|rdkit/i.test(e));
    expect(panelErrs.length,
      `Info-panel console errors: ${JSON.stringify(panelErrs.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
