/* ---
sub_features_covered: [bio.analyze.compare-sequences, bio.analyze.compare-sequences.top-menu, bio.annotate.manage-annotations, bio.annotate.manage-annotations.top-menu, bio.annotate.scan-liabilities, bio.annotate.scan-liabilities.rules, bio.annotate.scan-liabilities.top-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {openBioMenu, dialog, cancelDialog, setChoice, setText, setBool, loadBioTable} from './helpers';

test.use(specTestOptions);

// Antibody variable regions: two Macromolecule columns in one notation and alphabet, so the pair
// drives Compare Sequences, and real liability motifs (NG/NS/M/NxS) make the scanner produce hits.
const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const ROWS = 40;

const colNames = (page: Page): Promise<string[]> =>
  page.evaluate(() => grok.shell.tv.dataFrame.columns.names());

// Clearing the last annotation leaves the tag as an empty string, not '[]'.
const storedAnnotations = (page: Page): Promise<number> => page.evaluate(() => {
  const tag = grok.shell.tv.dataFrame.col('AntibodyHC')!.getTag('.annotations');
  return tag ? JSON.parse(tag).length : 0;
});

test('Bio Analyze | Compare sequences — defaults, custom result name, identical-column guard', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await loadBioTable(page, ANTIBODIES, ROWS, 'antibodies');

  await softStep('Compare sequences with defaults appends an auto-named difference column', async () => {
    const before = await colNames(page);
    await openBioMenu(page, 'Analyze', 'Compare-sequences...');
    const dlg = dialog(page, 'Compare-Sequences');
    await dlg.waitFor({timeout: 60_000});
    expect((await dlg.locator('.d4-dialog-title').textContent())?.trim()).toContain('Compare Sequences');
    // Both pickers defaulting to one column would make OK fail by construction.
    const c1 = await dlg.locator('[name="input-host-Sequence-column-1"] select').inputValue();
    const c2 = await dlg.locator('[name="input-host-Sequence-column-2"] select').inputValue();
    expect(c1).not.toBe(c2);
    await dlg.locator('[name="button-OK"]').click();
    await dlg.waitFor({state: 'detached', timeout: 30_000});

    await page.waitForFunction((n) => grok.shell.tv.dataFrame.columns.length > n, before.length,
      {timeout: 60_000});
    const added = (await colNames(page)).filter((n) => !before.includes(n));
    expect(added).toEqual([c1 + ' vs ' + c2]);
    const diff = await page.evaluate((n) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(n)!;
      const hc = df.col('AntibodyHC')!;
      const lc = df.col('AntibodyLC')!;
      let mismatched = 0;
      for (let i = 0; i < col.length; i++)
        if ((col.get(i) ?? '') !== `${hc.get(i)}#${lc.get(i)}`) mismatched++;
      return {renderer: col.getTag('cell.renderer'), mismatched};
    }, added[0]);
    expect(diff.renderer).toBe('MacromoleculeDifference');
    // The column is the two sources paired per row — that pairing is the whole feature.
    expect(diff.mismatched).toBe(0);
  });

  await softStep('A custom result name and swapped columns are both honoured', async () => {
    const before = await colNames(page);
    await openBioMenu(page, 'Analyze', 'Compare-sequences...');
    const dlg = dialog(page, 'Compare-Sequences');
    await dlg.waitFor({timeout: 60_000});
    await setChoice(page, 'Compare-Sequences', 'Sequence-column-1', 'AntibodyLC');
    await setChoice(page, 'Compare-Sequences', 'Sequence-column-2', 'AntibodyHC');
    await setText(page, 'Compare-Sequences', 'Result-column-name', 'LC minus HC');
    await dlg.locator('[name="button-OK"]').click();
    await dlg.waitFor({state: 'detached', timeout: 30_000});

    await page.waitForFunction((n) => grok.shell.tv.dataFrame.columns.length > n, before.length,
      {timeout: 60_000});
    expect((await colNames(page)).filter((n) => !before.includes(n))).toEqual(['LC minus HC']);
  });

  await softStep('Choosing the same column twice is rejected and adds nothing', async () => {
    const before = await colNames(page);
    await page.evaluate(() => {
      for (const b of Array.from(document.querySelectorAll('.d4-balloon, .grok-balloon-error'))) b.remove();
    });
    await openBioMenu(page, 'Analyze', 'Compare-sequences...');
    const dlg = dialog(page, 'Compare-Sequences');
    await dlg.waitFor({timeout: 60_000});
    await setChoice(page, 'Compare-Sequences', 'Sequence-column-2', 'AntibodyHC');
    await setChoice(page, 'Compare-Sequences', 'Sequence-column-1', 'AntibodyHC');
    await dlg.locator('[name="button-OK"]').click();
    await dlg.waitFor({state: 'detached', timeout: 30_000});

    await expect(page.locator('.d4-balloon.error, .grok-balloon-error').first())
      .toContainText(/distinct columns/i, {timeout: 30_000});
    expect(await colNames(page)).toEqual(before);
  });

  finishSpec();
});

test('Bio Annotate | Scan Liabilities — default rules, then a narrowed rule set with a summary column',
  async ({page}) => {
    test.setTimeout(300_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await loadBioTable(page, ANTIBODIES, ROWS, 'antibodies');

    let defaultHits = 0;
    await softStep('Default rule set writes the companion annotation column', async () => {
      await openBioMenu(page, 'Annotate', 'Scan-Liabilities...');
      const dlg = dialog(page, 'Scan-Sequence-Liabilities');
      await dlg.waitFor({timeout: 60_000});
      // Free Cysteine ships off, everything else on — pin the shipped defaults before changing them.
      await expect(dlg.locator('[name="input-host-Deamidation-(NG)"] input')).toBeChecked();
      await expect(dlg.locator('[name="input-host-Free-Cysteine"] input')).not.toBeChecked();
      await expect(dlg.locator('[name="input-host-Create-annotation-column"] input')).toBeChecked();
      await expect(dlg.locator('[name="input-host-Create-summary-count-column"] input')).not.toBeChecked();
      await dlg.locator('[name="button-OK"]').click();
      await dlg.waitFor({state: 'detached', timeout: 60_000});

      await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.names()
        .includes('~AntibodyHC_annotations'), null, {timeout: 60_000});
      defaultHits = await page.evaluate(() => {
        const col = grok.shell.tv.dataFrame.col('~AntibodyHC_annotations')!;
        let n = 0;
        for (let i = 0; i < col.length; i++) {
          const raw = col.get(i);
          if (raw) n += JSON.parse(raw).length;
        }
        return n;
      });
      expect(defaultHits).toBeGreaterThan(0);
      expect(await colNames(page)).not.toContain('AntibodyHC_liability_count');
      // Each hit must point at the motif it claims to have matched, at the position it reports.
      const misplaced = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const seq = df.col('AntibodyHC')!;
        const ann = df.col('~AntibodyHC_annotations')!;
        let bad = 0;
        let checked = 0;
        for (let i = 0; i < ann.length; i++) {
          const raw = ann.get(i);
          if (!raw) continue;
          for (const hit of JSON.parse(raw)) {
            if (!hit.matchedMonomers) continue;
            checked++;
            if ((seq.get(i) ?? '').substr(hit.positionIndex, hit.matchedMonomers.length) !== hit.matchedMonomers)
              bad++;
          }
        }
        return {bad, checked};
      });
      expect(misplaced.checked).toBeGreaterThan(0);
      expect(misplaced.bad).toBe(0);
    });

    await softStep('Only the oxidation rules, output switched to the summary count column', async () => {
      await openBioMenu(page, 'Annotate', 'Scan-Liabilities...');
      const dlg = dialog(page, 'Scan-Sequence-Liabilities');
      await dlg.waitFor({timeout: 60_000});
      for (const rule of ['Deamidation-(NG)', 'Deamidation-(NS)', 'Deamidation-(NA)', 'Deamidation-(ND)',
        'Deamidation-(NT)', 'Isomerization-(DG)', 'Isomerization-(DS)', 'N-glycosylation'])
        await setBool(page, 'Scan-Sequence-Liabilities', rule, false);
      await setBool(page, 'Scan-Sequence-Liabilities', 'Create-annotation-column', false);
      await setBool(page, 'Scan-Sequence-Liabilities', 'Create-summary-count-column', true);
      await dlg.locator('[name="button-OK"]').click();
      await dlg.waitFor({state: 'detached', timeout: 60_000});

      await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.names()
        .includes('AntibodyHC_liability_count'), null, {timeout: 60_000});
      const summary = await page.evaluate(() => {
        const col = grok.shell.tv.dataFrame.col('AntibodyHC_liability_count')!;
        let total = 0;
        for (let i = 0; i < col.length; i++) total += col.get(i) ?? 0;
        return {type: col.type, total};
      });
      expect(summary.type).toBe('int');
      // Two oxidation rules out of ten can only match a subset of what the full set matched.
      expect(summary.total).toBeGreaterThan(0);
      expect(summary.total).toBeLessThan(defaultHits);
    });

    finishSpec();
  });

test('Bio Annotate | Manage Annotations — lists, deletes and clears the scanned annotations',
  async ({page}) => {
    test.setTimeout(300_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await loadBioTable(page, ANTIBODIES, ROWS, 'antibodies');

    await softStep('Seed annotations through Scan Liabilities', async () => {
      await openBioMenu(page, 'Annotate', 'Scan-Liabilities...');
      const dlg = dialog(page, 'Scan-Sequence-Liabilities');
      await dlg.waitFor({timeout: 60_000});
      await dlg.locator('[name="button-OK"]').click();
      await dlg.waitFor({state: 'detached', timeout: 60_000});
      await page.waitForFunction(() => {
        const tag = grok.shell.tv.dataFrame.col('AntibodyHC')!.getTag('.annotations');
        return !!tag && JSON.parse(tag).length > 0;
      }, null, {timeout: 60_000});
    });

    await softStep('The dialog lists one row per annotation and offers Clear All', async () => {
      await openBioMenu(page, 'Annotate', 'Manage-Annotations...');
      const dlg = dialog(page, 'Manage-Annotations');
      await dlg.waitFor({timeout: 60_000});
      const stored = await storedAnnotations(page);
      expect(stored).toBeGreaterThan(0);
      expect(await dlg.locator('i[class*="fa-trash"]').count()).toBe(stored);
      expect(await dlg.locator('[name="button-Clear-All"]').count()).toBe(1);
    });

    await softStep('The trash icon drops exactly one annotation', async () => {
      const before = await storedAnnotations(page);
      await dialog(page, 'Manage-Annotations').locator('i[class*="fa-trash"]').first().click();
      await expect(dialog(page, 'Manage-Annotations').locator('i[class*="fa-trash"]'))
        .toHaveCount(before - 1, {timeout: 30_000});
      expect(await storedAnnotations(page)).toBe(before - 1);
    });

    await softStep('Clear All empties both the list and the column tag', async () => {
      await dialog(page, 'Manage-Annotations').locator('[name="button-Clear-All"]').click();
      await expect(dialog(page, 'Manage-Annotations'))
        .toContainText('No annotations on this column.', {timeout: 30_000});
      expect(await storedAnnotations(page)).toBe(0);
    });

    await cancelDialog(page, 'Manage-Annotations');
    finishSpec();
  });
