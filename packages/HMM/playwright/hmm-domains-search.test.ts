/* ---
sub_features_covered: [hmm.domains.top-menu, hmm.search.top-menu, hmm.panels.anarci, hmm.panels.pfam]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const ROWS = 40;

async function openBioMenu(page: Page, group: string, leaf: string): Promise<void> {
  const groupLoc = page.locator(`[name="div-Bio---${group}"]`);
  if (!await groupLoc.isVisible().catch(() => false))
    await page.locator('[name="div-Bio"]').click();
  await groupLoc.waitFor({state: 'visible', timeout: 15_000});
  const box = (await groupLoc.boundingBox())!;
  await page.mouse.move(box.x + 5, box.y + box.height / 2, {steps: 12});
  await page.mouse.move(box.x + box.width - 15, box.y + box.height / 2, {steps: 12});
  const leafLoc = page.locator(`[name="div-Bio---${group}---${leaf}"]`);
  await leafLoc.waitFor({state: 'visible', timeout: 15_000});
  await leafLoc.click();
}

test('HMM: Find Domains, Profile HMM Search and the sequence panels', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async ({path, rows}) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const src = await grok.dapi.files.readCsv(path);
    const df = src.clone(DG.BitSet.create(src.rowCount, (i) => i < rows));
    df.name = 'antibodies';
    grok.shell.addTableView(df);
    for (let i = 0; i < 200; i++) {
      if (df.col('AntibodyHC')?.semType === 'Macromolecule' && document.querySelector('[name="viewer-Grid"] canvas'))
        break;
      await new Promise((r) => setTimeout(r, 200));
    }
  }, {path: ANTIBODIES, rows: ROWS});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 60_000});

  await softStep('S1: Find Domains (built-in Pfam library, gathering cutoffs)', async () => {
    await openBioMenu(page, 'Annotate', 'Find-Domains-(HMMER)...');
    const dlg = page.locator('[name="dialog-Find-Domains-(HMMER)"]');
    await dlg.waitFor({timeout: 60_000});
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC domains'), null,
      {timeout: 120_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const summary = Array.from({length: df.rowCount}, (_, i) => df.get('AntibodyHC domains', i) as string);
      const definitions = JSON.parse(df.col('AntibodyHC')!.getTag('.annotations') || '[]').map((a: any) => a.id);
      return {rows: df.rowCount, vset: summary.filter((s) => s.includes('V-set')).length,
        c1: summary.filter((s) => s.includes('C1-set')).length, definitions};
    });
    expect(info.vset).toBeGreaterThan(info.rows * 0.8);
    expect(info.c1).toBeGreaterThan(info.rows * 0.5);
    expect(info.definitions.some((id: string) => id.startsWith('hmm-PF07686'))).toBe(true);
  });
  await page.evaluate(() => {
    const grid = grok.shell.tv.grid;
    grid.scrollToCell('AntibodyHC', 0);
    grid.col('AntibodyHC')!.width = 900;
  });
  await page.waitForTimeout(1000);
  await page.screenshot({path: 'playwright/test-output/hmm-domains.png'});

  await softStep('S2: Profile HMM Search with the built-in V-set family', async () => {
    await openBioMenu(page, 'Search', 'Profile-HMM-Search...');
    const dlg = page.locator('[name="dialog-Profile-HMM-Search"]');
    await dlg.waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const select = document.querySelector('[name="dialog-Profile-HMM-Search"] [name="input-host-Family"] select') as
        HTMLSelectElement | null;
      return !!select && select.value.startsWith('V-set');
    }, null, {timeout: 60_000});
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('V-set significant'), null,
      {timeout: 120_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const significant = Array.from({length: df.rowCount}, (_, i) => df.get('V-set significant', i)).filter((v) => v);
      return {rows: df.rowCount, significant: significant.length, format: df.col('V-set E-value')!.meta.format};
    });
    expect(info.significant).toBeGreaterThan(info.rows * 0.8);
    expect(info.format).toBe('scientific');
  });

  await softStep('S3: ANARCI and Pfam Domains panels for a heavy chain', async () => {
    const text = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const value = DG.SemanticValue.fromTableCell(df.cell(1, 'AntibodyHC'));
      const render = async (name: string) => {
        const widget: DG.Widget = await grok.functions.call(`HMM:${name}`, {sequence: value});
        for (let i = 0; i < 100 && widget.root.querySelector('.grok-loader'); i++)
          await new Promise((r) => setTimeout(r, 100));
        return widget.root.textContent ?? '';
      };
      return {anarci: await render('anarciPanel'), pfam: await render('pfamDomainsPanel')};
    });
    expect(text.anarci).toContain('Heavy');
    expect(text.anarci).toContain('IGHV');
    expect(text.anarci).toContain('CDR3 (IMGT)');
    expect(text.pfam).toContain('V-set');
  });
  finishSpec(stepErrors);
});
