/* ---
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.add, chem.analyze.scaffold-tree.viewer]
--- */
// github-3004: Scaffold Tree from active TableView must bind to the active TableView, not first-opened (fixed 1.21.0).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Chem: github-3004 Scaffold Tree from active TableView binds to active TableView (multi-table)', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
    });
    await page.waitForFunction(() => Array.from((window as any).grok.shell.tableViews).length === 0, null, {timeout: 30000});
  });

  await softStep('Open TableView A (tableA) + TableView B (tableB)', async () => {
    await page.evaluate(async () => {
      const dfA = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      dfA.name = 'tableA';
      grok.shell.addTableView(dfA);
      const dfB = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      dfB.name = 'tableB';
      grok.shell.addTableView(dfB);
    });
  });

  await softStep('Wait for Chem menu registration (semType + cellRenderer + scaffold-tree.add)', async () => {
    await waitForChemMenu(page);
  });

  await softStep('Verify both tables have Molecule semType (poll up to 60s for both)', async () => {
    await page.waitForFunction(() => {
      const tables = Array.from((window as any).grok.shell.tables) as any[];
      return tables.length >= 2 && tables.every((t: any) => t.columns.toList().some((c: any) => c.semType === 'Molecule'));
    }, null, {timeout: 60000});
    const states = await page.evaluate(() => (Array.from(grok.shell.tables) as any[]).map((t: any) => ({
      name: t.name,
      hasMol: t.columns.toList().some((c: any) => c.semType === 'Molecule'),
    })));
    expect(states.length >= 2 && states.every((s: any) => s.hasMol),
      `Tables missing Molecule semType after 60s poll: ${JSON.stringify(states)}`).toBe(true);
  });

  await softStep('Force tableB to be active', async () => {
    await page.evaluate(async () => {
      const tvs = Array.from(grok.shell.tableViews) as any[];
      const tvB = tvs.find((tv: any) => tv.dataFrame?.name === 'tableB');
      if (tvB) (grok.shell as any).v = tvB;
    });
    await page.waitForFunction(() => (window as any).grok.shell.tv?.dataFrame?.name === 'tableB', null, {timeout: 10000});
  });

  await softStep('Pre-condition: active table is tableB after step 2', async () => {
    const state = await page.evaluate(() => ({
      shell_t_name: grok.shell.t?.name,
      shell_tv_df_name: grok.shell.tv?.dataFrame?.name,
      tables_open: grok.shell.tables.map((t: any) => t.name),
    }));
    expect(
      state.shell_t_name,
      `Pre-condition failed: active table expected 'tableB', got '${state.shell_t_name}'. Setup invariant broken — Datagrok default newly-added-view-becomes-active changed? tables_open=${JSON.stringify(state.tables_open)}`,
    ).toBe('tableB');
    expect(state.shell_tv_df_name).toBe('tableB');
  });

  await softStep('Invoke top-menu Chem | Analyze | Scaffold Tree', async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-menu-item-label', {hasText: 'Scaffold Tree'}).first().waitFor({state: 'attached', timeout: 30000});
    await page.evaluate(() => {
      const stItem = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(el => el.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      if (!stItem) throw new Error('Top-menu "Scaffold Tree" sub-menu item not found');
      (stItem.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 10000});
  });

  await softStep('Assert regression-lock: Scaffold Tree viewer binds to tableB (active), NOT tableA', async () => {
    const result = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const stViewer = Array.from(tv.viewers).find((v: any) => v.type === 'Scaffold Tree');
      return {
        viewer_present: !!stViewer,
        viewer_df_name: stViewer ? (stViewer as any).dataFrame?.name : null,
        active_tv_df_name: tv?.dataFrame?.name,
      };
    });
    expect(result.viewer_present, 'Scaffold Tree viewer not attached to active TableView').toBe(true);
    expect(
      result.viewer_df_name,
      `github-3004 regression: Scaffold Tree viewer bound to '${result.viewer_df_name}' instead of 'tableB' (active). Multi-table state-tracking bug landed again. active_tv_df=${result.active_tv_df_name}`,
    ).toBe('tableB');
  });

  await softStep('Verify TableView A has no Scaffold Tree viewer (no spurious binding)', async () => {
    const result = await page.evaluate(() => {
      const tvs = Array.from(grok.shell.tableViews);
      const tvA = tvs.find((tv: any) => tv.dataFrame?.name === 'tableA');
      if (!tvA) return {tvA_found: false};
      const stOnA = Array.from(tvA.viewers).filter((v: any) => v.type === 'Scaffold Tree').length;
      return {tvA_found: true, scaffold_tree_count_on_A: stOnA};
    });
    expect(result.tvA_found, 'TableView A not found in shell.tableViews — multi-table state setup broken').toBe(true);
    expect(
      result.scaffold_tree_count_on_A,
      `github-3004 regression: Scaffold Tree viewer spuriously attached to TableView A — multi-table state-tracking error. count=${result.scaffold_tree_count_on_A}`,
    ).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
