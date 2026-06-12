import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const CYCLIZED_PATH = 'System:AppData/SequenceTranslator/samples/cyclized.csv';
const SIRNA_PATH = 'System:AppData/SequenceTranslator/samples/sirna-demo.csv';

async function validateViaHelper(page: Page, sequence: string): Promise<{error: string | null; returnValue: boolean | null}> {
  return page.evaluate(async (seq) => {
    const grok = (window as any).grok;
    let error: string | null = null;
    let returnValue: boolean | null = null;
    try {
      const helper = await grok.functions.call('SequenceTranslator:getTranslationHelper', {});
      const detector = helper.createFormatDetector(seq);
      const format = detector.getFormat();
      const validResult = format === null ? false : helper.createSequenceValidator(seq).isValidSequence(format);
      returnValue = Boolean(validResult);
    } catch (e: any) {
      error = String(e);
    }
    return {error, returnValue};
  }, sequence);
}

async function clickBioPolyToolItem(page: Page, itemName: string): Promise<void> {

  await page.locator('[name="div-Bio"]').click();
  await page.waitForTimeout(400);

  await page.locator('[name="div-Bio---PolyTool"]').hover();
  await page.waitForTimeout(600);

  await page.locator(`[name="${itemName}"]`).click();
}

test('SequenceTranslator — Macromolecule FASTA/Separator column lifecycle: notation-provider, PolyTool, API', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv(path);
    (window as any).grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 8000);
    });
    const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, CYCLIZED_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await page.waitForFunction(() => {
    return (window as any).DG.Func.find({
      name: 'applyNotationProviderForCyclized',
      package: 'SequenceTranslator',
    }).length > 0;
  }, null, {timeout: 60_000});

  await page.evaluate(async () => {
    await (window as any).grok.functions.call('SequenceTranslator:init', {});
    await new Promise((r) => setTimeout(r, 1000));
  });

  await softStep('Scenario 2 step 1: cyclized.csv seqs column detected as Macromolecule with custom units + separator tag', async () => {
    const probe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const seqsCol = df.col('seqs');
      if (!seqsCol) return {found: false};
      const tempKeys: string[] = [];
      if (seqsCol.temp) {
        for (const k in seqsCol.temp) tempKeys.push(k);
      }
      return {
        found: true,
        semType: seqsCol.semType,
        units: seqsCol.meta ? seqsCol.meta.units : null,
        aligned: seqsCol.getTag('aligned'),
        alphabet: seqsCol.getTag('alphabet'),
        alphabetIsMultichar: seqsCol.getTag('.alphabetIsMultichar'),
        separator: seqsCol.getTag('separator'),
        polyToolDataRole: seqsCol.getTag('polytool-data-role'),
        tempKeys: tempKeys,
        rowCount: df.rowCount,
      };
    });
    expect(probe.found, 'seqs column must exist in cyclized.csv').toBe(true);
    expect(probe.semType, 'seqs column semType must be Macromolecule after notation refiner fires').toBe('Macromolecule');
    expect(probe.units, 'seqs column units must be custom (set by applyNotationProviderForCyclized)').toBe('custom');
    expect(probe.aligned, 'seqs column aligned tag must be SEQ').toBe('SEQ');
    expect(probe.alphabet, 'seqs column alphabet tag must be UN').toBe('UN');
    expect(probe.alphabetIsMultichar, 'seqs column alphabetIsMultichar must be true (multi-char monomers in cyclized notation)').toBe('true');
    expect(probe.separator, 'seqs column separator tag must be present (cyclized notation uses separator)').toBeTruthy();
    expect(probe.polyToolDataRole, 'seqs column polytool-data-role tag must be template (set by apply step)').toBe('template');

    expect(probe.tempKeys.includes('seq-handler.notation-provider'),
      'col.temp must include seq-handler.notation-provider key (CyclizedNotationProvider attached)').toBe(true);
    expect(probe.rowCount, 'cyclized.csv must have 14 rows').toBe(14);
  });

  await softStep('Scenario 2 step 2: applyNotationProviderForCyclized JS API call sets expected column tags', async () => {

    const probe = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const seqsCol = df.col('seqs');

      await (window as any).grok.functions.call('SequenceTranslator:applyNotationProviderForCyclized', {
        col: seqsCol,
        separator: '-',
      });
      await new Promise((r) => setTimeout(r, 500));
      return {
        units: seqsCol.meta ? seqsCol.meta.units : null,
        aligned: seqsCol.getTag('aligned'),
        alphabet: seqsCol.getTag('alphabet'),
        alphabetIsMultichar: seqsCol.getTag('.alphabetIsMultichar'),
        polyToolDataRole: seqsCol.getTag('polytool-data-role'),
        separator: seqsCol.getTag('separator'),
      };
    });
    expect(probe.units, 'units must be custom after applyNotationProviderForCyclized').toBe('custom');
    expect(probe.aligned, 'aligned must be SEQ after applyNotationProviderForCyclized').toBe('SEQ');
    expect(probe.alphabet, 'alphabet must be UN after applyNotationProviderForCyclized').toBe('UN');
    expect(probe.alphabetIsMultichar, 'alphabetIsMultichar must be true after applyNotationProviderForCyclized').toBe('true');
    expect(probe.polyToolDataRole, 'polytool-data-role must be template after applyNotationProviderForCyclized').toBe('template');
    expect(probe.separator, 'separator tag must equal the provided separator after call').toBe('-');
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after applyNotationProviderForCyclized API call').toBe(0);
  });

  await softStep('Scenario 1 step 1: Bio | PolyTool | Convert... opens PolyTool Conversion dialog on cyclized column; CANCEL closes', async () => {
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});

    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Convert...');

    await page.locator('[name="dialog-PolyTool-Conversion"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-PolyTool-Conversion"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        columnInputPresent: !!(dlg?.querySelector('[name="input-host-Column"]')),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
        okPresent: !!(dlg?.querySelector('[name="button-OK"]')),
      };
    });
    expect(dlgProbe.title, 'PolyTool Conversion dialog title must be correct').toBe('PolyTool Conversion');
    expect(dlgProbe.columnInputPresent, 'Column input must be in PolyTool Conversion dialog').toBe(true);
    expect(dlgProbe.cancelPresent, 'CANCEL button must be in PolyTool Conversion dialog').toBe(true);
    await page.locator('[name="dialog-PolyTool-Conversion"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-PolyTool-Conversion"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after PolyTool Conversion CANCEL').toBe(0);
  });

  await softStep('Scenario 3 step 1: Bio | PolyTool | Combine Sequences... opens Combine Sequences dialog; CANCEL closes', async () => {

    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Combine-Sequences...');

    await page.locator('[name="dialog-Combine-Sequences"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Combine-Sequences"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        tableInputPresent: !!(dlg?.querySelector('[name="input-host-Table"]')),
        columnInputPresent: !!(dlg?.querySelector('[name="input-host-Column"]')),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
      };
    });
    expect(dlgProbe.title, 'Combine Sequences dialog title must be correct').toBe('Combine Sequences');
    expect(dlgProbe.tableInputPresent, 'Table input must be in Combine Sequences dialog').toBe(true);
    expect(dlgProbe.columnInputPresent, 'Column input must be in Combine Sequences dialog').toBe(true);
    expect(dlgProbe.cancelPresent, 'CANCEL button must be in Combine Sequences dialog').toBe(true);
    await page.locator('[name="dialog-Combine-Sequences"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Combine-Sequences"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Combine Sequences CANCEL').toBe(0);
  });

  await softStep('Scenario 4 step 1: validateSequence returns true for valid lowercase nucleotide sequence (FASTA format)', async () => {

    const validResult = await validateViaHelper(page, 'acgu');
    const alsoValidResult = await validateViaHelper(page, 'AGCT');
    const invalidResult = await validateViaHelper(page, 'notasequence12345');
    expect(validResult.error, 'validateSequence must not throw for valid input').toBeNull();
    expect(invalidResult.error, 'validateSequence must not throw for garbage input').toBeNull();
    expect(validResult.returnValue, "validateSequence('acgu') must return true for valid nucleotide sequence").toBe(true);
    expect(alsoValidResult.returnValue, "validateSequence('AGCT') must return true for valid nucleotide sequence").toBe(true);
    expect(invalidResult.returnValue, 'validateSequence with garbage string must return false').toBe(false);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon from validateSequence calls').toBe(0);
  });

  await softStep('Scenario 4 step 2: translateOligonucleotideSequence converts HELM to Axolabs (HELM-pivot)', async () => {
    const result = await page.evaluate(async (sirnaPath) => {

      const sirnadf = await (window as any).grok.dapi.files.readCsv(sirnaPath);
      const senseHelmVal = sirnadf.col('sense_helm') ? String(sirnadf.col('sense_helm').get(0)) : null;
      if (!senseHelmVal) return {error: 'sense_helm column not found in sirna-demo.csv', result: null};

      let translated = null;
      let error = null;
      try {
        translated = await (window as any).grok.functions.call(
          'SequenceTranslator:translateOligonucleotideSequence',
          {sequence: senseHelmVal, sourceFormat: 'HELM', targetFormat: 'Axolabs'}
        );
      } catch (e: any) {
        error = String(e);
      }
      return {
        result: translated ? String(translated).slice(0, 80) : null,
        error,
      };
    }, SIRNA_PATH);
    expect(result.error, 'translateOligonucleotideSequence HELM->Axolabs must not throw').toBeNull();
    expect(result.result, 'translateOligonucleotideSequence must return a non-empty Axolabs string').toBeTruthy();
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon from translateOligonucleotideSequence').toBe(0);
  });

  await softStep('Cleanup: close all views', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n');
    throw new Error(`Step failures:\n${summary}`);
  }
});
