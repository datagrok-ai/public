/* ---
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.filter, chem.analyze.scaffold-tree.viewer, chem.search.substructure, chem.search.substructure.api, chem.sketcher.cell-editor]
--- */
// GROK-12758: Scaffold Tree node Edit→CANCEL then per-node filter checkbox must not corrupt substructure-search state.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: GROK-12758 Scaffold Tree node Edit-then-Filter does not corrupt substructure-search state', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await softStep('Open spgi-100.csv + wait for Chem menu', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
    await waitForMolecule(page);
    await page.evaluate(() => {
      (window as any).__df = grok.shell.t;
      (window as any).__grok12758_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__grok12758_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
  });

  await softStep('Verify Molecule semType + Grid renders', async () => {
    const result = await page.evaluate(() => {
      const df = (window as any).__df;
      const cols = df.columns.toList().map((c: any) => ({name: c.name, semType: c.semType}));
      const hasMol = cols.some((c: any) => c.semType === 'Molecule');
      return {hasMol, rowCount: df.rowCount, cols: cols.slice(0, 5)};
    });
    if (!result.hasMol)
      throw new Error(`Setup failed: no Molecule column on spgi-100.csv after 30s settle. cols=${JSON.stringify(result.cols)}`);
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 10000});
  });

  await softStep('Open Scaffold Tree via top-menu Chem | Analyze | Scaffold Tree', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      if (!chemMenu) throw new Error('Top-menu Chem entry not found');
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 700));
      const stItem = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(el => el.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      if (!stItem) throw new Error('Top-menu "Scaffold Tree" sub-menu item not found');
      (stItem.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 10000});
  });

  await softStep('Generate scaffold tree via viewer.generateTree() JS API', async () => {
    const ok = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const viewer: any = Array.from(tv.viewers).find((v: any) => /Scaffold Tree/i.test(v.type || ''));
      if (!viewer) return {ok: false, reason: 'Scaffold Tree viewer not in tv.viewers'};
      try {
        if (typeof viewer.generateTree === 'function') {
          await viewer.generateTree();
          return {ok: true, method: 'viewer.generateTree()'};
        }
        await (window as any).grok.functions.call('Chem:GenerateScaffoldTree', {table: tv.dataFrame});
        return {ok: true, method: 'Chem:GenerateScaffoldTree'};
      } catch (e) {
        return {ok: false, reason: String(e)};
      }
    });
    expect((ok as any).ok, `Tree generation failed: ${JSON.stringify(ok)}`).toBe(true);
    await page.waitForTimeout(3000);
  });

  await softStep('Wait for scaffold tree to populate with ≥4 visible nodes (poll up to 90s)', async () => {
    const finalState = await page.evaluate(async () => {
      for (let i = 0; i < 18; i++) {
        const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
        if (viewer) {
          const allNodes = Array.from(viewer.querySelectorAll('.d4-tree-view-node'));
          const visibleNodes = allNodes.filter(n => n.querySelector('canvas'));
          const toolbar = viewer.querySelector('.chem-scaffold-tree-toolbar');
          const emptyState = (toolbar?.className ?? '').includes('empty-tree');
          if (!emptyState && visibleNodes.length >= 4) {
            return {ok: true, visibleNodes: visibleNodes.length, totalNodes: allNodes.length, toolbarCls: toolbar?.className, atMs: (i+1) * 5000};
          }
          if (i === 17) {
            return {ok: false, visibleNodes: visibleNodes.length, totalNodes: allNodes.length, toolbarCls: toolbar?.className, emptyState};
          }
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      return {ok: false, reason: 'viewer disappeared'};
    });
    if (!(finalState as any).ok)
      throw new Error(`Scaffold tree did not populate ≥4 visible nodes within 90s poll. state=${JSON.stringify(finalState)}`);
  });

  await softStep('Open 4th scaffold node in Sketcher (Edit pencil icon)', async () => {
    const opened = await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
      if (!viewer) return {ok: false, reason: 'no viewer'};
      const visibleNodes = Array.from(viewer.querySelectorAll('.d4-tree-view-node'))
        .filter(n => n.querySelector('canvas'));
      const fourth = visibleNodes[3];
      if (!fourth) return {ok: false, reason: `only ${visibleNodes.length} visible nodes`};
      const editBtn = fourth.querySelector('[aria-label="Edit scaffold"]') as HTMLElement | null;
      if (!editBtn) return {ok: false, reason: 'no Edit pencil icon on 4th node'};
      editBtn.click();
      return {ok: true};
    });
    expect((opened as any).ok, `Edit-scaffold click failed: ${JSON.stringify(opened)}`).toBe(true);
    await page.locator('.d4-dialog').waitFor({timeout: 8000});
    const dialogTitle = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      const title = dlg?.querySelector('.d4-dialog-header, .d4-dialog-title');
      return title?.textContent?.trim() ?? '';
    });
    expect(dialogTitle, `Expected "Edit Scaffold..." dialog; got "${dialogTitle}"`).toMatch(/Edit Scaffold/i);
  });

  await softStep('CANCEL the Edit Scaffold dialog (no edits applied)', async () => {
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
    // Wait for the dialog to actually detach rather than sleeping a fixed interval before the count.
    await page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 8000});
    const dialogCount = await page.evaluate(() => document.querySelectorAll('.d4-dialog').length);
    expect(dialogCount, 'Dialog did not close after CANCEL').toBe(0);
  });

  await softStep('Click filter checkbox on the SAME 4th scaffold node (bug-trigger step)', async () => {
    const clicked = await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
      if (!viewer) return {ok: false};
      const visibleNodes = Array.from(viewer.querySelectorAll('.d4-tree-view-node'))
        .filter(n => n.querySelector('canvas'));
      const fourth = visibleNodes[3];
      const checkbox = fourth?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
      if (!checkbox) return {ok: false, reason: 'no checkbox on 4th node'};
      const beforeChecked = checkbox.checked;
      checkbox.click();
      return {ok: true, beforeChecked, afterChecked: checkbox.checked};
    });
    expect((clicked as any).ok, `Checkbox click failed: ${JSON.stringify(clicked)}`).toBe(true);
    // Substructure filter applies asynchronously — poll for it to settle (some rows masked, not all)
    // instead of a blind sleep before the next step reads df.filter.
    await expect.poll(async () => page.evaluate(() => {
      const df = (window as any).__df;
      return df.filter.trueCount > 0 && df.filter.trueCount < df.rowCount;
    }), {timeout: 12_000, intervals: [250, 500, 1000]}).toBe(true);
  });

  await softStep('Assert clean filter state — no searchSubstructure errors, BitSet honored, selection untouched', async () => {
    const state = await page.evaluate(() => {
      const df = (window as any).__df;
      const errs = ((window as any).__grok12758_errors as string[] | undefined) ?? [];
      const searchErrors = errs.filter(e => /Search pattern cannot be set|searchSubstructure/i.test(e));
      return {
        rowCount: df.rowCount,
        filterTrue: df.filter.trueCount,
        selectionTrue: df.selection.trueCount,
        searchErrors: searchErrors.slice(0, 5),
        totalErrors: errs.length,
      };
    });
    // (A1) No searchSubstructure-class errors
    expect(
      state.searchErrors.length,
      `GROK-12758 regression: searchSubstructure-class errors fired after Edit-Scaffold + checkbox sequence. errors=${JSON.stringify(state.searchErrors)}`,
    ).toBe(0);
    // (A2) Filter applied cleanly: some rows masked, but not all-zero
    expect(
      state.filterTrue,
      `GROK-12758 regression: filter.trueCount=${state.filterTrue} of ${state.rowCount} — expected strictly between 0 and rowCount.`,
    ).toBeGreaterThan(0);
    expect(state.filterTrue).toBeLessThan(state.rowCount);
    // (A3) Selection untouched (the bug surface includes "crossed-out" molecule rendering — selection BitSet pollution)
    expect(
      state.selectionTrue,
      `GROK-12758 regression: selection.trueCount=${state.selectionTrue} after checkbox click — should be 0 (selection BitSet must NOT be modified by filter action).`,
    ).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
