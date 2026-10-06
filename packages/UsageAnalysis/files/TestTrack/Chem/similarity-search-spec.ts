import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: Similarity Search', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 100));
    }
    // The flat settle stood in for the Molecule semType landing on the opened frame.
    const semDeadline = Date.now() + 5000;
    while (Date.now() < semDeadline && !df.columns.toList().some((c: any) => c.semType === 'Molecule'))
      await new Promise(r => setTimeout(r, 100));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Chem → Search → Similarity Search → viewer appears', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      // Labels from a previously opened menu stay in the document, so only a node that was not
      // already there is this menu's leaf; clicking a stale one actuates nothing.
      const stale = new Set(Array.from(document.querySelectorAll('.d4-menu-item-label')));
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      const find = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => !stale.has(m) && m.textContent!.trim() === 'Similarity Search...') as HTMLElement | undefined;
      const menuDeadline = Date.now() + 500;
      let sim = find();
      while (!sim && Date.now() < menuDeadline) { await new Promise(r => setTimeout(r, 25)); sim = find(); }
      if (!sim)
        sim = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => m.textContent!.trim() === 'Similarity Search...') as HTMLElement | undefined;
      (sim!.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      // Wait for the viewer the assertion below reads, not for a flat interval.
      const attachDeadline = Date.now() + 6000;
      while (Date.now() < attachDeadline &&
             !Array.from(grok.shell.tv.viewers).some((v: any) => /Similarity/i.test(v.type || '')))
        await new Promise(r => setTimeout(r, 100));
    });
    const hasSim = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => /Similarity/i.test(v.type || '')));
    expect(hasSim).toBe(true);
  });

  await softStep('Modify viewer options (fingerprint/limit/metric/cutoff) without error', async () => {
    const results = await page.evaluate(async () => {
      const simViewer: any = Array.from(grok.shell.tv.viewers).find((v: any) => /Similarity/i.test(v.type || ''));
      if (!simViewer) return {error: 'viewer not found'};
      const res: Record<string, boolean> = {};
      try { simViewer.setOptions({fingerprint: 'Pattern'}); await new Promise(r => setTimeout(r, 1500)); res.fingerprint = true; } catch (e) { res.fingerprint = false; }
      try { simViewer.setOptions({limit: 5}); await new Promise(r => setTimeout(r, 1500)); res.limit = true; } catch (e) { res.limit = false; }
      try { simViewer.setOptions({distanceMetric: 'Dice'}); await new Promise(r => setTimeout(r, 1500)); res.metric = true; } catch (e) { res.metric = false; }
      try { simViewer.setOptions({cutoff: 1.0}); await new Promise(r => setTimeout(r, 1500)); res.cutoff = true; } catch (e) { res.cutoff = false; }
      simViewer.setOptions({fingerprint: 'Morgan', limit: 12, distanceMetric: 'Tanimoto', cutoff: 0.01});
      return res;
    });
    expect(results.fingerprint).toBe(true);
    expect(results.limit).toBe(true);
    expect(results.metric).toBe(true);
    expect(results.cutoff).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
