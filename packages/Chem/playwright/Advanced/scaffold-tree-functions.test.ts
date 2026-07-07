/* ---
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.add, chem.analyze.scaffold-tree.filter, chem.analyze.scaffold-tree.generate, chem.analyze.scaffold-tree.viewer]
--- */
// Paired scenario: Advanced/scaffold-tree-functions.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Chem: Scaffold Tree add + generate + node-click filter + toolbox + property-panel', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

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
    // Wait for the Molecule semType to actually land before mounting the viewer —
    // a Scaffold Tree mounted before detection has no molecule column and Generate
    // produces nothing (the menu attaches earlier than the async detector).
    await waitForMolecule(page);
  });

  await softStep('Step 2-3: Chem → Analyze → Scaffold Tree → viewer mounted (empty state)', async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .some(m => m.textContent!.trim() === 'Scaffold Tree'), null, {timeout: 15000});
    await page.evaluate(() => {
      const item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Scaffold Tree') as HTMLElement;
      (item.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 15000});
    const emptyState = await page.evaluate(() => {
      const toolbar = document.querySelector('.chem-scaffold-tree-toolbar');
      return toolbar?.className?.includes('empty-tree') ?? false;
    });
    expect(emptyState).toBe(true);
  });

  await softStep('Step 4-5: Click magic-wand → scaffold tree generates → nodes appear', async () => {
    // Fire the magic-wand "Generate" once via the JS-API viewer root (the icon
    // there is interactive immediately; the named-DOM node can lag). Do NOT
    // re-click — a second click interrupts an in-flight generation and it never
    // settles within the poll window.
    // Deliver a REAL mouse click via a Playwright locator (the magic-wand d4
    // handler does not respond to a synthetic Element.click() under the test
    // runner). Close any lingering menu overlay first.
    await page.keyboard.press('Escape').catch(() => {});
    await page.locator('[name="viewer-Scaffold-Tree"] [aria-label="Generate"]').first().click({timeout: 15000});
    // Generation can take 30-90s (cold RDKit / server load); poll from the
    // Playwright side (short evaluates) rather than one long-running in-page loop —
    // a single 120s page.evaluate blocking the main thread starves the async
    // generation work and the tree never fills. Count molecule nodes (each a
    // <canvas>, excluding the hidden root).
    const visibleNodeCount = async () => page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
      const nodes = viewer?.querySelectorAll('.d4-tree-view-node') ?? [];
      return Array.from(nodes).filter(n => n.querySelector('canvas')).length;
    });
    await expect.poll(visibleNodeCount,
      {timeout: 120_000, intervals: [2000], message: 'Scaffold tree did not populate within 120s'})
      .toBeGreaterThanOrEqual(1);
  });

  await softStep('Step 6-7: Click first scaffold node → table filters', async () => {
    const beforeFiltered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    const click = await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Scaffold-Tree"]');
      const nodes = Array.from(viewer!.querySelectorAll('.d4-tree-view-node'))
        .filter(n => n.querySelector('canvas'));
      if (nodes.length === 0) return {ok: false, reason: 'no visible nodes'};
      const checkbox = nodes[0].querySelector('input[type="checkbox"]') as HTMLInputElement;
      if (!checkbox) return {ok: false, reason: 'no checkbox on first node'};
      checkbox.click();
      return {ok: true};
    });
    expect((click as any).ok, (click as any).reason).toBe(true);
    await expect.poll(() => page.evaluate(() => grok.shell.t.filter.trueCount),
      {timeout: 20000, intervals: [500], message: 'Scaffold filter did not apply'})
      .toBeLessThan(beforeFiltered);
    const after = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(after, 'Filter hid every row — scaffold match is broken').toBeGreaterThan(0);
  });

  await softStep('Step 8: Viewer toolbox renders without errors', async () => {
    const res = await page.evaluate(() => {
      const toolbar = document.querySelector('[name="viewer-Scaffold-Tree"] .chem-scaffold-tree-toolbar');
      if (!toolbar) return {ok: false, icons: 0};
      const icons = toolbar.querySelectorAll('i, button, [class*="fa-"], [class*="grok-icon"], .d4-ribbon-item').length;
      return {ok: true, icons};
    });
    expect(res.ok, 'Scaffold Tree toolbar not rendered').toBe(true);
    expect(res.icons, 'Scaffold Tree toolbar rendered no action icons').toBeGreaterThan(0);
  });

  await softStep('Step 9: Open property panel via viewer click + Properties...', async () => {
    const focused = await page.evaluate(() => {
      const viewer = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Scaffold Tree');
      grok.shell.o = viewer;
      return grok.shell.o != null && (grok.shell.o as any).type === 'Scaffold Tree';
    });
    expect(focused, 'Context Panel not focused on the Scaffold Tree viewer').toBe(true);
    // Context Panel switched to the viewer's property set — accordion panes render
    await expect.poll(() => page.evaluate(() =>
      document.querySelectorAll('.d4-accordion-pane-header').length),
      {timeout: 15000, intervals: [500], message: 'Viewer property panel did not render'})
      .toBeGreaterThan(0);
    // Expand any collapsed property group and confirm interactable controls render
    await page.evaluate(() => {
      for (const pane of Array.from(document.querySelectorAll('.d4-accordion-pane'))) {
        if (!pane.querySelector('input, select, textarea'))
          (pane.querySelector('.d4-accordion-pane-header') as HTMLElement)?.click();
      }
    });
    await expect.poll(() => page.evaluate(() =>
      document.querySelectorAll('.d4-accordion-pane input, .d4-accordion-pane select, .d4-accordion-pane textarea').length),
      {timeout: 10000, intervals: [500], message: 'Property panel exposed no interactable controls'})
      .toBeGreaterThan(0);
  });

  await softStep('Final: no console errors during Scaffold Tree walk', async () => {
    const errs = await page.evaluate(() => ((window as any).__stf_errors ?? []) as string[]);
    const stfErrs = errs.filter(e => /scaffold|searchSubstructure/i.test(e));
    expect(stfErrs.length, `Console errors: ${JSON.stringify(stfErrs.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
