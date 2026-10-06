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
      // A docked-but-blank WebLogo would still satisfy the binding check, so prove it painted.
      const painted = await page.evaluate(async () => {
        for (let i = 0; i < 60; i++) {
          const cnv = document.querySelector('[name="viewer-WebLogo"] canvas') as HTMLCanvasElement | null;
          if (cnv && cnv.width > 0) {
            const px = cnv.getContext('2d')!.getImageData(0, 0, cnv.width, cnv.height).data;
            for (let j = 3; j < px.length; j += 4) if (px[j] !== 0) return true;
          }
          await new Promise((r) => setTimeout(r, 500));
        }
        return false;
      });
      expect(painted).toBe(true);
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
        const alignedCol = await page.evaluate((n) => {
          const col = grok.shell.tv.dataFrame.col(n)!;
          const vals: string[] = [];
          for (let i = 0; i < col.length; i++) vals.push(col.get(i) ?? '');
          const filled = vals.filter((v) => v.replace(/-/g, '').length > 0);
          return {scheme: col.getTag('.numberingScheme'), positionNames: col.getTag('.positionNames') ?? '',
            lengths: Array.from(new Set(vals.map((v) => v.length))), filled: filled.length, rows: vals.length};
        }, aligned!);
        expect(alignedCol.scheme).toBe('kabat');
        // The point of an aligned column: every row padded to one width, and actually numbered.
        expect(alignedCol.lengths).toHaveLength(1);
        expect(alignedCol.filled).toBe(alignedCol.rows);
        expect(alignedCol.positionNames.split(',').length).toBeGreaterThan(10);
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
      const slicesMatchSource = await page.evaluate(() => {
        const src = grok.shell.tv.dataFrame.col('AntibodyHC')!;
        const reg = grok.shell.tv.dataFrame.col('HC 20-60')!;
        for (let i = 0; i < reg.length; i++)
          if ((reg.get(i) ?? '') !== (src.get(i) ?? '').substring(19, 60)) return i;
        return -1;
      });
      expect(slicesMatchSource).toBe(-1);
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
        const scores = await page.evaluate((n) => {
          const col = grok.shell.tv.dataFrame.col(n)!;
          const vals: number[] = [];
          for (let i = 0; i < col.length; i++) vals.push(col.get(i));
          return {type: col.type, first: vals[0], min: Math.min(...vals), max: Math.max(...vals),
            nulls: vals.filter((v) => v == null || Number.isNaN(v)).length};
        }, added[0]);
        expect(scores.first).toBe(1);
        expect(scores.nulls).toBe(0);
        expect(scores.min).toBeGreaterThanOrEqual(0);
        expect(scores.max).toBeLessThanOrEqual(1);
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
        return {first: vals[0], min: Math.min(...vals), max: Math.max(...vals),
          nulls: vals.filter((v) => v == null || Number.isNaN(v)).length};
      }, added[0]);
      expect(stats.first).toBe(stats.max);
      expect(stats.nulls).toBe(0);
      expect(stats.min).toBeGreaterThanOrEqual(0);
      expect(stats.max).toBeLessThanOrEqual(1);
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
        const src = grok.shell.tv.dataFrame.col('sequence')!;
        const vals: string[] = [];
        for (let i = 0; i < col.length; i++) vals.push(col.get(i) ?? '');
        return {units: col.meta.units, semType: col.semType, sample: vals[0],
          allHelm: vals.every((v) => v.startsWith('PEPTIDE1{') && v.includes('}$')),
          monomers: vals[0].substring(vals[0].indexOf('{') + 1, vals[0].indexOf('}')).split('.').length,
          srcLength: (src.get(0) ?? '').length};
      }, added[0]);
      expect(meta.units).toBe('helm');
      expect(meta.semType).toBe('Macromolecule');
      expect(meta.allHelm).toBe(true);
      // One HELM monomer per residue of the source peptide.
      expect(meta.monomers).toBe(meta.srcLength);
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
      const mol = await page.evaluate((names) => {
        const cols = names.map((n: string) => grok.shell.tv.dataFrame.col(n)!);
        const molCol = cols.find((c: any) => c.semType === 'Molecule');
        return molCol ? {
          // Set from the Highlight monomers checkbox — proves the changed setting reached the run.
          highlight: molCol.getTag('.sequence-src-highlight-monomers'),
          molblock: (molCol.get(0) ?? '').slice(0, 200),
          empty: (() => { let n = 0; for (let i = 0; i < molCol.length; i++) if (!molCol.get(i)) n++; return n; })(),
        } : null;
      }, added);
      expect(mol).not.toBeNull();
      expect(mol!.highlight).toBe('true');
      expect(mol!.molblock).toContain('V3000');
      expect(mol!.empty).toBe(0);
    });

    finishSpec();
  });
