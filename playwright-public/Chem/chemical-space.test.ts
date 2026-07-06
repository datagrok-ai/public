/* ---
sub_features_covered: [chem.analyze.chemical-space, chem.analyze.chemical-space.editor, chem.analyze.chemical-space.embeddings, chem.analyze.chemical-space.top-menu, chem.analyze.chemical-space.transform]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

async function openDatasetAndWaitForMolecule(page: Page, label: string, datasetPath: string) {
  await softStep(`[${label}] Open ${datasetPath} + wait for Chem menu (Molecule semType)`, async () => {
    const isSdf = datasetPath.toLowerCase().endsWith('.sdf');
    await page.evaluate(async ({path, isSdf}) => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      for (let i = 0; i < 50 && grok.shell.tv != null; i++) await new Promise(r => setTimeout(r, 100));
      if (isSdf) {
        await ((DG as any).Func.find({name: 'OpenFile'})[0])
          .prepare({fullPath: path}).call(undefined, undefined, {processed: false});
      } else {
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
      }
    }, {path: datasetPath, isSdf});
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await waitForChemMenu(page);
    await waitForMolecule(page);
  });
}

async function openChemicalSpaceDialog(page: Page, label: string) {
  await softStep(`[${label}] Open Chem → Analyze → Chemical Space dialog`, async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .some(m => m.textContent!.trim() === 'Chemical Space...'), null, {timeout: 5000});
    await page.evaluate(() => {
      const cs = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Chemical Space...') as HTMLElement | undefined;
      if (!cs) throw new Error('"Chemical Space..." sub-menu item not found');
      (cs.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title, `Dialog title expected /Chem(ical)? Space/, got "${title}"`).toMatch(/Chem(ical)?\s*Space/i);
    // .md step 2: editor must surface the method (UMAP/t-SNE) selector + Cluster MCS toggle.
    await expect(page.locator('.d4-dialog [name="input-Method"]'), `[${label}] method selector missing`).toHaveCount(1);
    await expect(page.locator('.d4-dialog [name="input-Cluster-MCS"]'), `[${label}] Cluster MCS toggle missing`).toHaveCount(1);
  });
}

async function clickOkAndWaitForEmbedding(page: Page, label: string, minSuffix: number,
  expectClusterCol = false) {
  let foundSuffix = -1;
  await softStep(`[${label}] OK → embedding scatter plot + Embed_X_N/Embed_Y_N (N > ${minSuffix})`, async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async (lastSuffix: number) => {
      const hasTag = (col: any, re: RegExp) => col != null &&
        Array.from(col.tags.keys()).some((k: any) => re.test(k));
      for (let i = 0; i < 45; i++) {
        const df = grok.shell.tv?.dataFrame;
        if (!df) { await new Promise(r => setTimeout(r, 2000)); continue; }
        const embedCols = df.columns.toList().filter((c: any) => /^Embed_X_(\d+)$/.test(c.name));
        const hasScatter = Array.from(grok.shell.tv.viewers).some((v: any) => v.type === 'Scatter plot');
        if (embedCols.length > 0 && hasScatter) {
          const maxSuffix = Math.max(...embedCols.map((c: any) => parseInt(c.name.match(/^Embed_X_(\d+)$/)![1])));
          if (maxSuffix > lastSuffix) {
            const xName = `Embed_X_${maxSuffix}`, yName = `Embed_Y_${maxSuffix}`;
            // Tags are set asynchronously ~5-8s after column creation; poll up to 20s for both axes.
            let xTag = false, yTag = false;
            for (let j = 0; j < 10; j++) {
              const cols = grok.shell.tv?.dataFrame;
              xTag = hasTag(cols?.col(xName), /chem-space-embedding-col/);
              yTag = hasTag(cols?.col(yName), /chem-space-embedding-col/);
              if (xTag && yTag) break;
              await new Promise(r => setTimeout(r, 2000));
            }
            const df2 = grok.shell.tv!.dataFrame;
            const xCol = df2.col(xName), yCol = df2.col(yName);
            const scatterBound = Array.from(grok.shell.tv!.viewers).some((v: any) =>
              v.type === 'Scatter plot' && v.props.xColumnName === xName && v.props.yColumnName === yName);
            const clusterCol = df2.columns.toList().some((c: any) => hasTag(c, /chem-space-cluster-col/));
            return {ok: true, maxSuffix, hasYCol: yCol != null, xTag, yTag, scatterBound, clusterCol,
              xMissing: xCol?.stats?.missingValueCount, yMissing: yCol?.stats?.missingValueCount};
          }
        }
        await new Promise(r => setTimeout(r, 2000));
      }
      const df = grok.shell.tv?.dataFrame;
      return {ok: false, dfCols: df?.columns?.toList()?.map((c: any) => c.name)?.slice(0, 30),
        viewerTypes: Array.from(grok.shell.tv?.viewers ?? []).map((v: any) => v.type)};
    }, minSuffix);
    const r = result as any;
    expect(r.ok,
      `[${label}] Chemical Space did not produce a fresh Embed_X_N (N > ${minSuffix}) within 90s. result=${JSON.stringify(result)}`,
    ).toBe(true);
    expect(r.hasYCol, `[${label}] paired Embed_Y_${r.maxSuffix} column not created`).toBe(true);
    expect(r.xTag && r.yTag,
      `[${label}] Embed_X/Embed_Y missing '.%chem-space-embedding-col' tag (x=${r.xTag}, y=${r.yTag})`).toBe(true);
    expect(r.xMissing, `[${label}] Embed_X_${r.maxSuffix} has ${r.xMissing} missing values (embedding not computed)`).toBe(0);
    expect(r.yMissing, `[${label}] Embed_Y_${r.maxSuffix} has ${r.yMissing} missing values (embedding not computed)`).toBe(0);
    expect(r.scatterBound,
      `[${label}] no Scatter plot bound to Embed_X_${r.maxSuffix}/Embed_Y_${r.maxSuffix}`).toBe(true);
    if (expectClusterCol)
      expect(r.clusterCol,
        `[${label}] Cluster MCS run did not add a '.%chem-space-cluster-col'-tagged column`).toBe(true);
    foundSuffix = r.maxSuffix;
  });
  return foundSuffix;
}

async function runChemicalSpaceWalk(page: Page, label: string, datasetPath: string,
  customParam: 'method' | 'cluster-mcs') {

  await openDatasetAndWaitForMolecule(page, label, datasetPath);
  await openChemicalSpaceDialog(page, label);
  const defaultSuffix = await clickOkAndWaitForEmbedding(page, `${label}/defaults`, 0);
  await openChemicalSpaceDialog(page, `${label}/reopen`);

  await softStep(`[${label}] Change param: ${customParam}`, async () => {
    if (customParam === 'method') {
      await page.locator('.d4-dialog [name="input-Method"]').selectOption('t-SNE');
      expect(await page.locator('.d4-dialog [name="input-Method"]').inputValue()).toBe('t-SNE');
    } else {
      const toggle = page.locator('.d4-dialog [name="input-Cluster-MCS"]');
      await toggle.check();
      await expect(toggle).toBeChecked();
    }
  });

  // .md step 6: cluster-mcs run must materialise a chem-space-cluster-col; t-SNE only re-runs the embedding
  // (deterministic method-vs-method difference is stochastic — not asserted, see deferred).
  await clickOkAndWaitForEmbedding(page, `${label}/custom`, defaultSuffix, customParam === 'cluster-mcs');

  await softStep(`[${label}] Close active view`, async () => {
    await page.evaluate(() => {
      // Dismiss any lingering editor dialog first — a still-open modal keeps the table view
      // alive, so closeAll() alone won't null grok.shell.tv (the cause of the 5s timeout here).
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        cancel?.click();
      });
      grok.shell.closeAll();
    });
    await page.waitForFunction(() => grok.shell.tv == null, null, {timeout: 15_000});
  });
}

test('Chem: Chemical Space multi-format walk (smiles-50 / molV2000 / molV3000)', async ({page}) => {
  test.setTimeout(600_000); // 6 dim-reduction runs (2 × 3 datasets) @ ~45-90s each + margin

  const consoleErrors: string[] = [];
  const isBenignError = (t: string) => /Failed to load resource/.test(t) || /404 \(\)/.test(t) || /favicon/.test(t);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);

  await runChemicalSpaceWalk(page, 'D1 smiles-50', 'System:AppData/Chem/tests/smiles-50.csv', 'method');
  await runChemicalSpaceWalk(page, 'D2 molV2000', 'System:AppData/Chem/mol1K.sdf', 'method');
  await runChemicalSpaceWalk(page, 'D3 molV3000', 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf', 'cluster-mcs');

  await softStep('No console errors across the multi-format walk', async () =>
    expect(consoleErrors, `console errors: ${consoleErrors.join(' | ')}`).toEqual([]));

  finishSpec();
});
