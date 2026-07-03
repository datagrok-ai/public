import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {finishSpec} from '../../helpers/viewers';

test.use(specTestOptions);

const tiles = '.chem-viewer-grid > div';

test('Chem: Similarity Search', async ({page}) => {
  test.setTimeout(120_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
    grok.shell.addTableView(df);
    for (let i = 0; i < 100; i++) {
      if (df.col('canonical_smiles')?.semType === 'Molecule') break;
      await new Promise(r => setTimeout(r, 100));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Chem → Search → Similarity Search → viewer appears with hits', async () => {
    await page.evaluate(() => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-menu-item-label').filter({hasText: 'Similarity Search...'}).first().waitFor({state: 'attached', timeout: 30000});
    await page.evaluate(() => {
      const sim = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Similarity Search...') as HTMLElement;
      (sim.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.waitForFunction(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => /Similarity/i.test(v.type || '')), null, {timeout: 30000});
    const hasSim = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => /Similarity/i.test(v.type || '')));
    expect(hasSim, 'Chem Similarity Search viewer added').toBe(true);
    await expect.poll(() => page.locator(tiles).count(),
      {timeout: 30000, message: 'top-N similar-molecule tiles rendered'}).toBeGreaterThan(0);
  });

  await softStep('Property knobs (fingerprint/limit/metric/cutoff) re-query without errors', async () => {
    const warnBefore = await page.evaluate(() => (grok.shell.warnings || []).length);
    const errBefore = pageErrors.length;

    const enums = await page.evaluate(() => {
      const v: any = Array.from(grok.shell.tv.viewers).find((x: any) => /Similarity/i.test(x.type || ''));
      const fps: string[] = (v.getProperty('fingerprint')?.choices ?? []).slice();
      const metrics: string[] = (v.getProperty('distanceMetric')?.choices ?? []).slice();
      const look = v.getOptions().look;
      return {
        cur: {fp: look.fingerprint, metric: look.distanceMetric},
        fps: fps.length ? fps : ['Morgan', 'RDKit', 'MACCS', 'AtomPair', 'TopologicalTorsion', 'Pattern'],
        metrics: metrics.length ? metrics : ['Tanimoto', 'Asymmetric', 'Cosine', 'Sokal'],
      };
    });

    const setOpt = async (opt: Record<string, any>) => {
      await page.evaluate((o) => {
        const v: any = Array.from(grok.shell.tv.viewers).find((x: any) => /Similarity/i.test(x.type || ''));
        v.setOptions(o);
      }, opt);
    };

    // Fingerprint knob — exercise every available type; visit the default last so each step differs.
    const fpWalk = [...enums.fps.filter((f) => f !== enums.cur.fp), enums.cur.fp];
    for (const fp of fpWalk) {
      await setOpt({fingerprint: fp});
      await expect(page.locator('.chem-similarity-header'), `fingerprint ${fp} applied`)
        .toContainText(fp, {timeout: 30000});
      expect(await page.locator(tiles).count(), `hits render for fingerprint ${fp}`).toBeGreaterThan(0);
    }

    // Distance metric knob — exercise every available metric; default visited last.
    const metricWalk = [...enums.metrics.filter((m) => m !== enums.cur.metric), enums.cur.metric];
    for (const m of metricWalk) {
      await setOpt({distanceMetric: m});
      await expect(page.locator('.chem-similarity-header'), `metric ${m} applied`)
        .toContainText(m, {timeout: 30000});
      expect(await page.locator(tiles).count(), `hits render for metric ${m}`).toBeGreaterThan(0);
    }

    // Limit knob — decrease then increase; the rendered tile count must track the limit.
    await setOpt({limit: 5});
    await expect.poll(() => page.locator(tiles).count(),
      {timeout: 30000, message: 'tile count reflects limit=5'}).toBe(5);
    await setOpt({limit: 12});
    await expect.poll(() => page.locator(tiles).count(),
      {timeout: 30000, message: 'tile count reflects limit=12'}).toBe(12);

    // Cutoff knob — cutoff=1 keeps only the exact-match self-hit (smiles.csv is deduplicated).
    await setOpt({cutoff: 1});
    await expect.poll(() => page.locator(tiles).count(),
      {timeout: 30000, message: 'cutoff=1 leaves only the exact-match self-hit'}).toBe(1);
    await setOpt({cutoff: 0.01});

    const warnDelta = await page.evaluate((n) => (grok.shell.warnings || []).length - n, warnBefore);
    expect(warnDelta, 'no new warnings during property-knob walk').toBe(0);
    expect(pageErrors.slice(errBefore), 'no uncaught errors during property-knob walk').toEqual([]);
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
