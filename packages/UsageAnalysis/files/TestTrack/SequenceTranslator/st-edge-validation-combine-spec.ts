/* ---
sub_features_covered: [sequencetranslator.api.validate-sequence, sequencetranslator.detectors.context-menu-macromolecule, sequencetranslator.oligo-renderer.combine-sense-antisense]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

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

test('SequenceTranslator — Edge: validateSequence false-branch + combineSenseAntisense mismatched units', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.evaluate(async (path: string) => {
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
    const hasMacromolecule = (cols as any[]).some((c) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, SIRNA_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await page.waitForFunction(() => {
    return (window as any).DG.Func.find({
      name: 'convertHelmToOligoNucleotide',
      package: 'SequenceTranslator',
    }).length > 0;
  }, null, {timeout: 60_000});

  await page.evaluate(async () => {
    await (window as any).grok.functions.call('SequenceTranslator:init', {});

    await (window as any).grok.functions.call('SequenceTranslator:getTranslationHelper', {});
    await new Promise((r) => setTimeout(r, 1000));
  });

  await page.locator('[name="Browse"]').waitFor({state: 'attached', timeout: 30_000});

  await softStep('Scenario 1 step 1: validateSequence returns false for garbage input', async () => {
    const result = await validateViaHelper(page, 'NOTAVALIDSEQUENCE!@#$');
    expect(result.error,
      'validateSequence logic must not throw for garbage input — should return false, not throw').toBeNull();
    expect(result.returnValue,
      'validateSequence must return false for unrecognized garbage sequence (FormatDetector returns null)').toBe(false);
  });

  await softStep('Scenario 1 step 2: validateSequence returns false for truncated HELM string', async () => {
    const result = await validateViaHelper(page, 'RNA1{r(A)p.r(C)p');
    expect(result.error,
      'validateSequence logic must not throw for truncated HELM — should complete normally').toBeNull();
    expect(result.returnValue,
      'validateSequence must return false for truncated/incomplete HELM string').toBe(false);
  });

  await softStep('Scenario 2 step 1: validateSequence logic confirms garbage input is invalid (translator edge case)', async () => {

    const result = await validateViaHelper(page, 'NOTASEQUENCE!@#$');
    expect(result.error, 'validateSequence logic must not throw for the translator garbage-input sequence').toBeNull();
    expect(result.returnValue,
      'validateSequence returns false for the same garbage input used in Scenario 2 translator test').toBe(false);
  });

  await softStep('Scenario 3 step 1: validateSequence returns false for mixed-format (valid prefix + garbage suffix)', async () => {
    const result = await validateViaHelper(page, 'RNA1{r(A)p.r(C)p}$$$$ extra_garbage');
    expect(result.error,
      'validateSequence logic must not throw for mixed-format input (valid HELM prefix + garbage suffix)').toBeNull();
    expect(result.returnValue,
      'validateSequence must return false for mixed-format string — FormatDetector / SequenceValidator rejects it').toBe(false);
  });

  await softStep('Scenario 4 step 1: create synthetic table with mismatched units (sense=helm, antisense=fasta)', async () => {
    const result = await page.evaluate(async () => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;

      const senseHelmVal = 'RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$';
      const antisenseFastaVal = 'UGCA';

      const senseCol = DG.Column.fromList('string', 'sense_helm', [senseHelmVal, senseHelmVal]);
      senseCol.semType = 'Macromolecule';
      senseCol.setTag('quality', 'Macromolecule');
      senseCol.setTag('units', 'helm');
      senseCol.setTag('alphabet', 'RNA');
      senseCol.setTag('aligned', 'SEQ');

      const antiCol = DG.Column.fromList('string', 'antisense_fasta', [antisenseFastaVal, antisenseFastaVal]);
      antiCol.semType = 'Macromolecule';
      antiCol.setTag('quality', 'Macromolecule');
      antiCol.setTag('units', 'fasta');
      antiCol.setTag('alphabet', 'RNA');
      antiCol.setTag('aligned', 'SEQ');

      const df = DG.DataFrame.fromColumns([senseCol, antiCol]);
      grok.shell.addTableView(df);
      (window as any).__testDf = df;
      await new Promise((r) => setTimeout(r, 1000));

      return {
        rowCount: df.rowCount,
        senseUnits: senseCol.getTag('units'),
        antiUnits: antiCol.getTag('units'),
        senseValue: senseCol.get(0),
        antiValue: antiCol.get(0),
      };
    });
    expect(result.rowCount, 'synthetic DataFrame must have 2 rows').toBe(2);
    expect(result.senseUnits, 'sense column units must be helm').toBe('helm');
    expect(result.antiUnits, 'antisense column units must be fasta (mismatched)').toBe('fasta');
    expect(result.senseValue, 'sense column must contain HELM value').toBeTruthy();
    expect(result.antiValue, 'antisense column must contain FASTA value').toBeTruthy();
  });

  await softStep('Scenario 4 step 2: combineSenseAntisenseToOligoNucleotide with mismatched units — graceful sense-only output', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = (window as any).__testDf;
      if (!df) return {
        error: 'synthetic DataFrame not found in __testDf',
        threwError: false,
        oligoColProduced: false,
        oligoColHasAntisense: false,
        oligoColValue0: null as string | null,
        errorBalloonCount: 0,
      };

      const senseCol = df.col('sense_helm');
      const antiCol = df.col('antisense_fasta');
      if (!senseCol || !antiCol) return {
        error: 'columns not found',
        threwError: false,
        oligoColProduced: false,
        oligoColHasAntisense: false,
        oligoColValue0: null as string | null,
        errorBalloonCount: 0,
      };

      let threwError = false;
      let errorMsg: string | null = null;
      let oligoCol: any = null;
      try {

        oligoCol = await grok.functions.call(
          'SequenceTranslator:combineSenseAntisenseToOligoNucleotide',
          {table: df, senseCol, antiCol},
        );
      } catch (e: any) {
        threwError = true;
        errorMsg = String(e);
      }

      await new Promise((r) => setTimeout(r, 2000));

      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const foundOligoCol = (cols as any[]).find((c) => {
        return c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide';
      });
      const oligoColProduced = !!foundOligoCol;

      let oligoColValue0: string | null = null;
      let oligoColHasAntisense = false;
      if (foundOligoCol) {
        oligoColValue0 = foundOligoCol.get(0) ?? null;
        oligoColHasAntisense = oligoColValue0 !== null && oligoColValue0.includes('|RNA2');
      }

      const errorBalloonCount = document.querySelectorAll('.d4-balloon.error').length;

      return {
        error: errorMsg,
        threwError,
        oligoColProduced,
        oligoColHasAntisense,
        oligoColValue0,
        errorBalloonCount,
      };
    });

    expect(result.threwError,
      `combineSenseAntisenseToOligoNucleotide must NOT throw with mismatched units (helm+fasta) — ` +
      `function gracefully handles non-HELM antisense by skipping it. Error: ${result.error ?? 'none'}`
    ).toBe(false);

    expect(result.oligoColProduced,
      'combineSenseAntisenseToOligoNucleotide must produce an OligoNucleotide column even ' +
      'when antisense is non-HELM (FASTA) — sense-only output is the documented graceful behavior'
    ).toBe(true);

    if (result.oligoColProduced) {
      expect(result.oligoColHasAntisense,
        `OligoNucleotide column produced from mismatched units (helm+fasta) must NOT contain ` +
        `an RNA2 antisense chain — FASTA antisense is silently skipped per converters.ts#L57. ` +
        `Cell value[0]: "${result.oligoColValue0 ?? 'null'}"`
      ).toBe(false);
    }
  });

  await page.evaluate(async () => {
    delete (window as any).__testDf;
    (window as any).grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    throw new Error(
      `${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  [${e.step}] ${e.error}`).join('\n'),
    );
  }
});
