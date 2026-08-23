/* ---
realizes: [chem.calculate.descriptors]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

const SMILES_CSV = 'System:DemoFiles/chem/smiles.csv';

async function openSmilesAndWaitForChem(page: Page) {
  await softStep('Open smiles.csv + wait for Molecule semType and Chem menu', async () => {
    const info = await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1000)); 
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 5000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200)); 
      }
      await new Promise((r) => setTimeout(r, 4000)); 
      return {rowCount: df.rowCount};
    }, SMILES_CSV);
    await waitForChemMenu(page);
    await waitForMolecule(page);
    // Read the semType back only after waitForMolecule returns: the open evaluate above races
    // detection against its own 5 s timeout, so a snapshot taken there can predate the barrier.
    const molCols = await page.evaluate(() => grok.shell.t.columns.toList()
      .filter((c: any) => c.semType === 'Molecule').map((c: any) => c.name));
    expect(info.rowCount, 'smiles.csv should open with 1000 rows').toBe(1000);
    expect(molCols.length, 'smiles.csv should carry a Molecule column').toBeGreaterThan(0);
  });
}

async function openDescriptorsDialog(page: Page) {
  await softStep('Chem > Calculate > Descriptors... opens the Chemical Descriptors dialog', async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .some((m) => m.textContent!.trim() === 'Calculate'), null, {timeout: 15000});
    await page.evaluate(() => {
      const calc = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((m) => m.textContent!.trim() === 'Calculate') as HTMLElement | undefined;
      if (!calc) throw new Error('"Calculate" menu group not found');
      const calcItem = calc.closest('.d4-menu-item') as HTMLElement;
      calcItem.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      calcItem.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .some((m) => m.textContent!.trim() === 'Descriptors...'), null, {timeout: 15000});
    await page.evaluate(() => {
      const descr = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((m) => m.textContent!.trim() === 'Descriptors...') as HTMLElement | undefined;
      if (!descr) throw new Error('"Descriptors..." leaf not found under Calculate');
      (descr.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="dialog-Chemical-Descriptors"]').waitFor({timeout: 15000});
    // The waitFor above resolves when the dialog root attaches; the Dart model fills the
    // Molecules selector after that, so the auto-selection is polled rather than snapshotted.
    await expect.poll(() => page.evaluate(() =>
      document.querySelector('[name="dialog-Chemical-Descriptors"] [name="input-host-Molecules"] .d4-column-selector-column')?.textContent?.trim() ?? ''),
    {message: 'Descriptors dialog should auto-select the SMILES molecule column (last observed selector text is shown as Received)', timeout: 10000})
      .toBe('canonical_smiles');
  });
}

const DLG = '[name="dialog-Chemical-Descriptors"]';

async function checkedSummaryCount(page: Page): Promise<number> {
  return page.evaluate((dlgSel) => {
    const dlg = document.querySelector(dlgSel)!;
    const txt = Array.from(dlg.querySelectorAll('*'))
      .map((e) => (e.textContent || '').trim())
      .find((t) => /^\d+ checked$/.test(t));
    return txt ? parseInt(txt, 10) : -1;
  }, DLG);
}

async function selectMolWtAndLogP(page: Page) {
  // Tree expanders and the None link take real clicks; the descriptor-node checkboxes do not —
  // a Playwright click leaves the dialog's own "N checked" summary unmoved (recon 2026-08-19,
  // 3/3), so they are toggled by setting .checked and dispatching input/change. The summary
  // poll below is the check that the dialog model actually registered the selection.
  await softStep('Select MolWt and MolLogP descriptors in the tree (checkbox state + dispatched events; tree clicks do not reach the dialog model)', async () => {
    const expanders = [
      `${DLG} [name="tree-expander-Descriptors"]`,
      `${DLG} [name="tree-expander-Crippen"]`,
    ];
    for (const sel of expanders) {
      const exp = page.locator(sel);
      const expanded = await exp.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded'));
      if (!expanded) await exp.click();
    }
    await page.locator(`${DLG} [name="tree-Descriptors---MolWt"] input[type="checkbox"]`).waitFor({state: 'visible', timeout: 10000});
    await page.locator(`${DLG} [name="tree-Crippen---MolLogP"] input[type="checkbox"]`).waitFor({state: 'visible', timeout: 10000});

    await page.locator(DLG).getByText('None', {exact: true}).click();
    await expect.poll(() => checkedSummaryCount(page),
      {message: 'None link should clear all persisted descriptor defaults', timeout: 5000}).toBe(0);

    await page.evaluate((dlgSel) => {
      const dlg = document.querySelector(dlgSel)!;
      for (const sel of ['[name="tree-Descriptors---MolWt"]', '[name="tree-Crippen---MolLogP"]']) {
        const cb = dlg.querySelector(`${sel} input[type="checkbox"]`) as HTMLInputElement;
        cb.checked = true;
        cb.dispatchEvent(new Event('input', {bubbles: true}));
        cb.dispatchEvent(new Event('change', {bubbles: true}));
      }
    }, DLG);

    await expect.poll(() => checkedSummaryCount(page),
      {message: 'exactly MolWt + MolLogP must be selected (Dart model registers the checkbox events)', timeout: 5000}).toBe(2);
  });
}

async function okAndVerifyColumnsAppended(page: Page, consoleErrors: string[]) {

  await softStep('Scenario 1 Step 5: OK submits to chem-chem container; MolWt + MolLogP columns append, row count unchanged', async () => {
    const before = await page.evaluate(() => {
      const df = grok.shell.t;
      (window as any).__descrBeforeCols = df.columns.names();
      return {beforeCols: df.columns.names(), beforeRows: df.rowCount};
    });
    await page.evaluate(() => {
      const df = grok.shell.t;
      const beforeCols: string[] = (window as any).__descrBeforeCols;
      (window as any).__descrAppended = new Promise<void>((resolve) => {
        const sub = df.onColumnsAdded.subscribe(() => {
          const newOnes = df.columns.names().filter((n: string) => !beforeCols.includes(n));
          if (newOnes.some((n: string) => /molwt/i.test(n)) && newOnes.some((n: string) => /logp/i.test(n))) {
            sub.unsubscribe();
            resolve();
          }
        });
      });
    });
    await page.locator('[name="dialog-Chemical-Descriptors"] [name="button-OK"]').click();
    const result = await page.evaluate(async (beforeRows) => {
      const df = grok.shell.t;
      const beforeCols: string[] = (window as any).__descrBeforeCols;
      await Promise.race([
        (window as any).__descrAppended,
        new Promise((r) => setTimeout(r, 120000)),
      ]);
      const newCols = df.columns.names().filter((n: string) => !beforeCols.includes(n));
      const ok = newCols.some((n: string) => /molwt/i.test(n)) && newCols.some((n: string) => /logp/i.test(n));
      return {ok, beforeRows, afterRows: df.rowCount, beforeColCount: beforeCols.length, newCols};
    }, before.beforeRows);
    console.log(`S1.5 appended ${result.newCols.length} columns ${JSON.stringify(result.newCols)} to a table of ` +
      `${result.beforeColCount} columns; rows ${result.beforeRows} -> ${result.afterRows}`);
    expect(result.ok, `MolWt + MolLogP columns did not append within 120s. newCols=${JSON.stringify(result.newCols)}`).toBe(true);
    // The dialog's selection is sticky and opens with five descriptors already checked, so a run that
    // failed to clear it appends extras: the appended set — not the presence of two names — is the claim.
    expect([...result.newCols].sort(),
      `exactly MolWt + MolLogP must be appended; extra columns mean the sticky selection was not cleared. newCols=${JSON.stringify(result.newCols)}`)
      .toEqual(['MolLogP', 'MolWt']);
    expect(result.afterRows, 'row count must be unchanged by descriptor calculation').toBe(result.beforeRows);
  });

  await softStep('Scenario 1 Step 4: the Descriptors dialog closes once the calculation is submitted', async () => {
    // Recon 2026-08-20 (chem.md) saw an OK click that issued no request at all; the dialog staying
    // open was what discriminated that state from a working submit.
    await expect(page.locator(DLG),
      'the Chemical Descriptors dialog must close after OK; a dialog still open means OK did not submit')
      .toHaveCount(0, {timeout: 30_000});
  });

  await softStep('Scenario 1 Step 6: MolWt and MolLogP carry non-empty numeric values (not all-null / all-zero)', async () => {
    const vals = await page.evaluate(() => {
      const df = grok.shell.t;
      const out: Record<string, any> = {};
      for (const cn of ['MolWt', 'MolLogP']) {
        const col = df.col(cn);
        if (!col) { out[cn] = null; continue; }
        let nonNull = 0, nonZero = 0;
        const n = Math.min(50, df.rowCount);
        for (let r = 0; r < n; r++) {
          const v = col.get(r);
          if (v !== null && v !== undefined && !Number.isNaN(v)) { nonNull++; if (v !== 0) nonZero++; }
        }
        out[cn] = {type: col.type, sampled: n, nonNull, nonZero, first: col.get(0), last: col.get(df.rowCount - 1)};
      }
      return out;
    });
    console.log(`S1.6 ${JSON.stringify(vals)}`);
    for (const cn of ['MolWt', 'MolLogP']) {
      const v = (vals as any)[cn];
      expect(v, `${cn} column missing after calculation`).not.toBeNull();
      expect(['double', 'float', 'int'], `${cn} must be a numeric column, got type "${v.type}"`).toContain(v.type);
      expect(v.nonNull, `${cn} is all-null over ${v.sampled} sampled rows — container returned no results`).toBe(v.sampled);
      expect(v.nonZero, `${cn} is all-zero over ${v.sampled} sampled rows — container returned no results`).toBeGreaterThan(0);
      expect(Number.isFinite(v.first), `${cn} first row is not a finite number (got ${v.first})`).toBe(true);
      expect(Number.isFinite(v.last), `${cn} last row is not a finite number (got ${v.last})`).toBe(true);
    }
  });

  await softStep('No console errors during dialog open, OK submit, and column-append', async () => {
    expect(consoleErrors, `console errors during Descriptors flow: ${JSON.stringify(consoleErrors)}`).toEqual([]);
  });
}

test('Chem: Descriptors via Docker — columns appended (happy path)', async ({page}) => {
  test.setTimeout(600_000);

  // Attached before the first dialog opens: the step below claims "no console errors during
  // dialog open", and a listener registered after the open cannot see them.
  const consoleErrors: string[] = [];
  const AMBIENT = [/favicon/i, /ResizeObserver loop/i, /Permissions policy violation/i,
    /Unable to find element in cloned iframe/i];
  const isAmbient = (text: string) => AMBIENT.some((re) => re.test(text));
  page.on('pageerror', (e) => { if (!isAmbient(String(e))) consoleErrors.push(String(e)); });
  page.on('console', (m) => { if (m.type() === 'error' && !isAmbient(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await openSmilesAndWaitForChem(page);
  await openDescriptorsDialog(page);
  await selectMolWtAndLogP(page);
  await okAndVerifyColumnsAppended(page, consoleErrors);

  finishSpec();
});

// Scenario 2 (container-unavailable bounded error, GROK-17621) is not automated: the chem-chem
// Docker container is provisioned and reachable on the dev validation stand, and the Playwright
// harness cannot stop a platform-managed container to establish the "unavailable" precondition.
// Recon 2026-08-19: the descriptor call succeeded in ~4.5 s with no error balloon, so the
// error-notification and unchanged-column-count observables are not producible here. The manual
// coverage lives in chem-descriptors-docker-ui.md, which realizes chem.int.descriptors-docker-bounded
// as a manual-only scenario.
