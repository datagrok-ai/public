// Paired scenario: Advanced/scaffold-tree-functions.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../../spec-login';
import {finishSpec} from '../../helpers/viewers';

test.use(specTestOptions);

test('Chem: Scaffold Tree add + generate + node-click filter + toolbox + property-panel', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await softStep('Step 1: Open smiles-50.csv', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      grok.shell.addTableView(df);
      (window as any).__stf_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__stf_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
    await waitForChemMenu(page);
  });

  await softStep('Step 2-3: Chem → Analyze → Scaffold Tree → viewer mounted (empty state)', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 600));
      const item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      (item.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 5000));
    });
    await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 15000});
    const emptyState = await page.evaluate(() => {
      const toolbar = document.querySelector('.chem-scaffold-tree-toolbar');
      return toolbar?.className?.includes('empty-tree') ?? false;
    });
    expect(emptyState).toBe(true);
  });

  await softStep('Step 4-5: Click magic-wand → scaffold tree generates → nodes appear', async () => {
    await page.evaluate(async () => {
      const wand = document.querySelector('[name="viewer-Scaffold-Tree"] [aria-label="Generate"]') as HTMLElement;
      wand?.click();
    });
    // Wait up to 60s for nodes
    const ready = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
        const nodes = viewer?.querySelectorAll('.d4-tree-view-node') ?? [];
        const visible = Array.from(nodes).filter(n => n.querySelector('canvas.chem-canvas'));
        if (visible.length >= 1) return {ok: true, count: visible.length};
        await new Promise(r => setTimeout(r, 2000));
      }
      return {ok: false};
    });
    expect((ready as any).ok, `Scaffold tree did not populate within 60s`).toBe(true);
  });

  await softStep('Step 6-7: Click first scaffold node → table filters', async () => {
    const before = await page.evaluate(() => grok.shell.t.rowCount);
    const click = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
      const nodes = Array.from(viewer!.querySelectorAll('.d4-tree-view-node'))
        .filter(n => n.querySelector('canvas.chem-canvas'));
      if (nodes.length === 0) return {ok: false, reason: 'no visible nodes'};
      const checkbox = nodes[0].querySelector('input[type="checkbox"]') as HTMLInputElement;
      if (!checkbox) return {ok: false, reason: 'no checkbox on first node'};
      checkbox.click();
      await new Promise(r => setTimeout(r, 3000));
      return {ok: true};
    });
    expect((click as any).ok).toBe(true);
    const after = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(after).toBeLessThan(before);
  });

  await softStep('Step 8: Viewer toolbox renders without errors', async () => {
    const ok = await page.evaluate(() => {
      const toolbar = document.querySelector('[name="viewer-Scaffold-Tree"] .chem-scaffold-tree-toolbar');
      return !!toolbar;
    });
    expect(ok).toBe(true);
  });

  await softStep('Step 9: Open property panel via viewer click + Properties...', async () => {
    await page.evaluate(async () => {
      const viewer = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Scaffold Tree');
      grok.shell.o = viewer;
      await new Promise(r => setTimeout(r, 2000));
    });
    // Verify Context Panel has expandable property panes for the viewer
    const hasPanes = await page.evaluate(() =>
      document.querySelectorAll('.d4-accordion-pane-header').length > 0);
    expect(hasPanes).toBe(true);
  });

  await softStep('Final: no console errors during Scaffold Tree walk', async () => {
    const errs = await page.evaluate(() => ((window as any).__stf_errors ?? []) as string[]);
    const stfErrs = errs.filter(e => /scaffold|searchSubstructure/i.test(e));
    expect(stfErrs.length, `Console errors: ${JSON.stringify(stfErrs.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
