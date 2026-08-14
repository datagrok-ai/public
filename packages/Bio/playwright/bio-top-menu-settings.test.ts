/* ---
sub_features_covered: [bio.analyze.composition.column-choice, bio.annotate.numbering-scheme.kabat, bio.calculate.get-region.custom-range, bio.calculate.identity.reference, bio.calculate.similarity.reference, bio.transform.convert-notation.helm, bio.transform.to-atomic-level.options]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {openBioMenu, dialog, setChoice, setText, setBool, loadBioTable, closeExtraViewers} from './helpers';

test.use(specTestOptions);

// Every menu entry below is already exercised with its defaults elsewhere in this suite; here each
// one runs with a changed setting, so a setting that is read but never applied shows up.
const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const PEPTIDES = 'System:AppData/Bio/samples/FASTA_PT_activity.csv';

const colNames = (page: Page): Promise<string[]> =>
  page.evaluate(() => grok.shell.tv.dataFrame.columns.names());

/** Runs `body` and returns the columns it appended to the table. */
async function addedColumns(page: Page, body: () => Promise<void>): Promise<string[]> {
  const before = await colNames(page);
  await body();
  await page.waitForFunction((n) => grok.shell.tv.dataFrame.columns.length > n, before.length,
    {timeout: 180_000});
  return (await colNames(page)).filter((n) => !before.includes(n));
}

test('Bio top menu on antibodies — Composition, numbering scheme and region run with changed settings',
  async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await loadBioTable(page, ANTIBODIES, 40, 'antibodies');

    await softStep('Composition on the non-default column docks a WebLogo bound to that column', async () => {
      await openBioMenu(page, 'Analyze', 'Composition');
      const dlg = dialog(page, 'Composition-Analysis');
      await dlg.waitFor({timeout: 60_000});
      await setChoice(page, 'Composition-Analysis', 'Column', 'AntibodyLC');
      await dlg.locator('[name="button-OK"]').click();
      await dlg.waitFor({state: 'detached', timeout: 30_000});

      await page.waitForFunction(() => Array.from(((grok.shell.tv as any)?.viewers ?? []) as any[])
        .some((v: any) => v.type === 'WebLogo'), null, {timeout: 120_000});
      const boundTo = await page.evaluate(() => {
        const wl = Array.from(((grok.shell.tv as any)?.viewers ?? []) as any[])
          .find((v: any) => v.type === 'WebLogo');
        return wl?.props?.sequenceColumnName ?? wl?.sequenceColumnName ?? null;
      });
      expect(boundTo).toBe('AntibodyLC');
      await closeExtraViewers(page);
    });

    await softStep('Apply Numbering Scheme with Kabat tags the result as Kabat, not the default IMGT',
      async () => {
        const added = await addedColumns(page, async () => {
          await openBioMenu(page, 'Annotate', 'Apply-Numbering-Scheme...');
          const dlg = dialog(page, 'Apply-Antibody-Numbering');
          await dlg.waitFor({timeout: 60_000});
          expect(await dlg.locator('[name="input-host-Scheme"] select').inputValue()).toBe('imgt');
          await setChoice(page, 'Apply-Antibody-Numbering', 'Scheme', 'kabat');
          await dlg.locator('[name="button-OK"]').click();
          await dlg.waitFor({state: 'detached', timeout: 60_000});
        });
        const aligned = added.find((n) => n.includes('aligned'));
        expect(aligned).toBeTruthy();
        const scheme = await page.evaluate((n) =>
          grok.shell.tv.dataFrame.col(n)!.getTag('.numberingScheme'), aligned!);
        expect(scheme).toBe('kabat');
      });

    await softStep('Extract Region honours a custom range and column name', async () => {
      const added = await addedColumns(page, async () => {
        await openBioMenu(page, 'Calculate', 'Extract-Region...');
        const dlg = dialog(page, 'Get-Sequence-Region');
        await dlg.waitFor({timeout: 60_000});
        await setChoice(page, 'Get-Sequence-Region', 'Start', '20');
        await setChoice(page, 'Get-Sequence-Region', 'End', '60');
        await setText(page, 'Get-Sequence-Region', 'Column-name', 'HC 20-60');
        await dlg.locator('[name="button-OK"]').click();
        await dlg.waitFor({state: 'detached', timeout: 60_000});
      });
      expect(added).toContain('HC 20-60');
      const lengths = await page.evaluate(() => {
        const col = grok.shell.tv.dataFrame.col('HC 20-60')!;
        const out: number[] = [];
        for (let i = 0; i < col.length; i++) out.push((col.get(i) ?? '').length);
        return out;
      });
      // Positions 20..60 inclusive — sequences shorter than 60 residues yield a shorter slice.
      expect(Math.max(...lengths)).toBe(41);
      expect(Math.min(...lengths)).toBeGreaterThan(0);
    });

    finishSpec();
  });

test('Bio top menu on peptides — scoring references, HELM conversion and atomic-level options',
  async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await loadBioTable(page, PEPTIDES, 20, 'peptides');

    const reference = await page.evaluate(() => grok.shell.tv.dataFrame.col('sequence')!.get(0)!);

    await softStep('Identity against an explicit reference scores that very row as a perfect match',
      async () => {
        const added = await addedColumns(page, async () => {
          await openBioMenu(page, 'Calculate', 'Identity...');
          const dlg = dialog(page, 'Identity');
          await dlg.waitFor({timeout: 60_000});
          await setText(page, 'Identity', 'Reference', reference);
          await dlg.locator('[name="button-OK"]').click();
          await dlg.waitFor({state: 'detached', timeout: 60_000});
        });
        expect(added).toHaveLength(1);
        expect(await page.evaluate((n) => grok.shell.tv.dataFrame.col(n)!.get(0), added[0])).toBe(1);
      });

    await softStep('Similarity against the same reference peaks on the reference row', async () => {
      const added = await addedColumns(page, async () => {
        await openBioMenu(page, 'Calculate', 'Similarity...');
        const dlg = dialog(page, 'Similarity');
        await dlg.waitFor({timeout: 60_000});
        await setText(page, 'Similarity', 'Reference', reference);
        await dlg.locator('[name="button-OK"]').click();
        await dlg.waitFor({state: 'detached', timeout: 60_000});
      });
      expect(added).toHaveLength(1);
      const stats = await page.evaluate((n) => {
        const col = grok.shell.tv.dataFrame.col(n)!;
        const vals: number[] = [];
        for (let i = 0; i < col.length; i++) vals.push(col.get(i));
        return {first: vals[0], max: Math.max(...vals)};
      }, added[0]);
      expect(stats.first).toBe(stats.max);
    });

    await softStep('Convert Sequence Notation to HELM instead of the default separator', async () => {
      const added = await addedColumns(page, async () => {
        await openBioMenu(page, 'Transform', 'Convert-Sequence-Notation...');
        const dlg = dialog(page, 'Convert-Sequence-Notation');
        await dlg.waitFor({timeout: 60_000});
        expect(await dlg.locator('[name="input-host-Convert-to"] select').inputValue()).toBe('separator');
        await setChoice(page, 'Convert-Sequence-Notation', 'Convert-to', 'helm');
        await dlg.locator('[name="button-OK"]').click();
        await dlg.waitFor({state: 'detached', timeout: 60_000});
      });
      expect(added).toHaveLength(1);
      const meta = await page.evaluate((n) => {
        const col = grok.shell.tv.dataFrame.col(n)!;
        return {units: col.meta.units, semType: col.semType, sample: col.get(0)};
      }, added[0]);
      expect(meta.units).toBe('helm');
      expect(meta.semType).toBe('Macromolecule');
      expect(meta.sample).toContain('PEPTIDE1{');
    });

    await softStep('To Atomic Level with monomer highlighting on and non-linear off', async () => {
      const added = await addedColumns(page, async () => {
        await openBioMenu(page, 'Transform', 'To-Atomic-Level...');
        const dlg = dialog(page, 'To-Atomic-Level');
        await dlg.waitFor({timeout: 60_000});
        await setBool(page, 'To-Atomic-Level', 'Non-linear', false);
        await setBool(page, 'To-Atomic-Level', 'Highlight-monomers', true);
        await dlg.locator('[name="button-OK"]').click();
        await dlg.waitFor({state: 'detached', timeout: 60_000});
      });
      const semTypes = await page.evaluate((names) => names
        .map((n: string) => grok.shell.tv.dataFrame.col(n)!.semType), added);
      expect(semTypes).toContain('Molecule');
    });

    finishSpec();
  });
