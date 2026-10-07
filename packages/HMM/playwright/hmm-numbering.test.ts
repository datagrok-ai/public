/* ---
sub_features_covered: [hmm.numbering.engine, hmm.numbering.bio-dialog, hmm.germlines.top-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const ROWS = 60;

/** Walks the pointer into a Bio top-menu group (d4 expands submenus off the mouse path). */
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

async function setChoice(page: Page, dialog: string, input: string, value: string): Promise<void> {
  const select = page.locator(`[name="dialog-${dialog}"] [name="input-host-${input}"] select`);
  await select.selectOption({label: value});
}

test('HMM: ANARCI numbering in Bio\'s numbering dialog, germlines from the Bio menu', async ({page}) => {
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

  await softStep('S1: the HMM engine is offered with all six ANARCI schemes', async () => {
    const info = await page.evaluate(() => {
      const engines = DG.Func.find({meta: {role: 'antibodyNumbering'}});
      const hmm = engines.find((f) => f.friendlyName === 'ANARCI (HMMER)');
      const scheme = hmm?.inputs.find((p) => p.name === 'scheme');
      return {labels: engines.map((f) => f.friendlyName), choices: scheme?.choices ?? [], packageName: hmm?.package?.name};
    });
    expect(info.labels).toContain('ANARCI (HMMER)');
    expect(info.choices).toEqual(['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy']);
  });

  await softStep('S2: Apply Numbering Scheme with ANARCI (HMMER), Chothia', async () => {
    await openBioMenu(page, 'Annotate', 'Apply-Numbering-Scheme...');
    const dlg = page.locator('[name="dialog-Apply-Antibody-Numbering"]');
    await dlg.waitFor({timeout: 60_000});
    // Sequence defaults to the first Macromolecule column, AntibodyHC.
    await setChoice(page, 'Apply-Antibody-Numbering', 'Engine', 'ANARCI (HMMER)');
    await setChoice(page, 'Apply-Antibody-Numbering', 'Scheme', 'chothia');
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC (aligned)'), null,
      {timeout: 120_000});
  });

  await softStep('S3: aligned column, Chothia insertions and FR/CDR annotations', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const aligned = df.col('AntibodyHC (aligned)')!;
      const source = df.col('AntibodyHC')!;
      const positions = (aligned.getTag('.positionNames') ?? '').split(', ');
      const regions = JSON.parse(aligned.getTag('.annotations') || '[]').map((a: any) => a.name);
      const annotationCol = df.col(source.getTag('.annotationColumnName') ?? '');
      let rowsWithRegions = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const hits = JSON.parse(annotationCol?.get(i) || '[]');
        if (hits.some((h: any) => h.endPositionIndex != null)) rowsWithRegions++;
      }
      const lengths = new Set(Array.from({length: df.rowCount}, (_, i) => (aligned.get(i) ?? '').length));
      return {scheme: aligned.getTag('.numberingScheme'), positions, regions, rowsWithRegions,
        rows: df.rowCount, lengths: [...lengths]};
    });
    expect(info.scheme).toBe('chothia');
    // Chothia places heavy-chain insertions at 31A, 52A and 82A-C.
    for (const p of ['31A', '52A', '82A', '82B', '82C']) expect(info.positions).toContain(p);
    expect(info.regions).toEqual(['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']);
    expect(info.rowsWithRegions).toBe(info.rows);
    expect(info.lengths.length).toBe(1);
  });
  await page.screenshot({path: 'playwright/test-output/hmm-numbering-chothia.png'});

  await softStep('S3b: IMGT aligned rows keep residue order (CDR3 insertions 111A.. then ..112A, 112)', async () => {
    await openBioMenu(page, 'Annotate', 'Apply-Numbering-Scheme...');
    const dlg = page.locator('[name="dialog-Apply-Antibody-Numbering"]');
    await dlg.waitFor({timeout: 60_000});
    await setChoice(page, 'Apply-Antibody-Numbering', 'Engine', 'ANARCI (HMMER)');
    await setChoice(page, 'Apply-Antibody-Numbering', 'Scheme', 'imgt');
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC (aligned) (2)'), null,
      {timeout: 120_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const aligned = df.col('AntibodyHC (aligned) (2)')!;
      const source = df.col('AntibodyHC')!;
      const positions = (aligned.getTag('.positionNames') ?? '').split(', ');
      const scrambled: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if ((aligned.get(i) ?? '').replace(/-/g, '') !== (source.get(i) ?? '').replace(/-/g, '')) scrambled.push(i);
      }
      const at = (p: string) => positions.indexOf(p);
      return {scrambled, has112A: at('112A') >= 0, order112: at('112A') < at('112') && at('111A') < at('112A')};
    });
    expect(info.has112A).toBe(true);
    expect(info.order112).toBe(true);
    expect(info.scrambled).toEqual([]);
  });

  await softStep('S4: Germlines and Species (ANARCI) adds annotated columns', async () => {
    await openBioMenu(page, 'Annotate', 'Germlines-and-Species-(ANARCI)...');
    const dlg = page.locator('[name^="dialog-Antibody-Germlines"]');
    await dlg.waitFor({timeout: 60_000});
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC v gene'), null,
      {timeout: 120_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const values = (name: string) => Array.from({length: df.rowCount}, (_, i) => df.get(name, i));
      return {chains: [...new Set(values('AntibodyHC chain'))], species: [...new Set(values('AntibodyHC species'))],
        genes: values('AntibodyHC v gene').filter((g) => g).length, rows: df.rowCount,
        evalues: values('AntibodyHC evalue').filter((e) => e > 0 && e < 1e-20).length};
    });
    // Some AntibodyHC entries are scFv (VL-linker-VH): ANARCI reports the N-terminal domain first.
    expect(info.chains).toContain('Heavy');
    expect(info.chains.every((c) => ['Heavy', 'Kappa', 'Lambda'].includes(c))).toBe(true);
    expect(info.species.every((s) => s === 'human' || s === 'mouse')).toBe(true);
    expect(info.genes).toBe(info.rows);
    expect(info.evalues).toBe(info.rows);
  });
  await page.screenshot({path: 'playwright/test-output/hmm-germlines.png'});
  finishSpec(stepErrors);
});
