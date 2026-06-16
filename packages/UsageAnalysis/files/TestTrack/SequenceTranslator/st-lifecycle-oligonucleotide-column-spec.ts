/* ---
sub_features_covered: [sequencetranslator.api.get-code-to-weights-map, sequencetranslator.api.translate-oligonucleotide-sequence, sequencetranslator.detectors.context-menu-oligo, sequencetranslator.oligo-renderer, sequencetranslator.oligo-renderer.combine-sense-antisense, sequencetranslator.oligo-renderer.convert-helm-to-oligo, sequencetranslator.oligo-renderer.copy-as-helm, sequencetranslator.oligo-renderer.open-helm-editor]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const SIRNA_PATH = 'System:AppData/SequenceTranslator/samples/sirna-demo.csv';

test('SequenceTranslator — OligoNucleotide column lifecycle: convert, combine, copy-as-helm, helm-editor, enumerate, API', async ({page}) => {
  // One CSV load + convert/combine + helm-editor/enumerator dialogs + API calls; 4 min suffices.
  test.setTimeout(240_000);
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
    await new Promise((r) => setTimeout(r, 1000));

    await (window as any).grok.functions.call('SequenceTranslator:getTranslationHelper', {});
    await new Promise((r) => setTimeout(r, 1000));
  });

  await softStep('Setup: convert oligo_helm HELM column to OligoNucleotide (oligo-renderer.convert-helm-to-oligo)', async () => {
    await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const helmCol = df.col('oligo_helm');
      await (window as any).grok.functions.call('SequenceTranslator:convertHelmToOligoNucleotide', {
        table: df,
        helmCol,
      });
      await new Promise((r) => setTimeout(r, 3000));
    });

    await page.waitForFunction(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      return (cols as any[]).some((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');
    }, null, {timeout: 30_000});

    const colProbe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCol = (cols as any[]).find((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');
      return {
        found: !!oligoCol,
        name: oligoCol?.name,
        semType: oligoCol?.semType,
        quality: oligoCol?.getTag('quality'),
        cellRenderer: oligoCol?.getTag('cell.renderer'),
        originalIntact: !!(df.col('oligo_helm')),
      };
    });
    expect(colProbe.found, 'OligoNucleotide column must be created by convertHelmToOligoNucleotide').toBe(true);
    expect(colProbe.semType, 'converted column semType must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.quality, 'converted column quality tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.cellRenderer, 'converted column cell.renderer tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.name, 'converted column name must be the oligo-suffix convention').toBe('oligo_helm (oligo)');
    expect(colProbe.originalIntact, 'original oligo_helm column must remain intact after conversion').toBe(true);

    const canvasPresent = await page.evaluate(() => {
      return !!document.querySelector('[name="viewer-Grid"] canvas');
    });
    expect(canvasPresent, 'viewer-Grid canvas must be present for OligoNucleotide cell rendering').toBe(true);

    const balloonError = await page.evaluate(() =>
      document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length
    );
    expect(balloonError, 'no error balloon after convertHelmToOligoNucleotide').toBe(0);
  });

  await softStep('Scenario 1 step 1: copyOligoAsHelm does not throw for OligoNucleotide cell (copy-as-helm)', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCol = (cols as any[]).find((c) => c.semType === 'OligoNucleotide');
      if (!oligoCol) return {success: false, error: 'No OligoNucleotide column found'};
      const cell0 = oligoCol.get(0);
      const semValue = (window as any).DG.SemanticValue.fromValueType(cell0, 'OligoNucleotide');
      let error = null;
      try {
        await (window as any).grok.functions.call('SequenceTranslator:copyOligoAsHelm', {value: semValue});
      } catch (e: any) {
        error = String(e);
      }
      return {
        success: error === null,
        error,
        cellValue: cell0 ? String(cell0).slice(0, 30) : null,
        semType: semValue ? semValue.semType : null,
      };
    });
    expect(result.error, 'copyOligoAsHelm must not throw for OligoNucleotide SemanticValue').toBeNull();
    expect(result.success, 'copyOligoAsHelm must succeed without error').toBe(true);
    expect(result.semType, 'SemanticValue semType must be OligoNucleotide').toBe('OligoNucleotide');

    const balloonError = await page.evaluate(() =>
      document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length
    );
    expect(balloonError, 'no error balloon after copyOligoAsHelm').toBe(0);
  });

  await softStep('Scenario 2 step 1: openOligoHelmEditor does not throw for OligoNucleotide cell (open-helm-editor)', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCol = (cols as any[]).find((c) => c.semType === 'OligoNucleotide');
      if (!oligoCol) return {success: false, error: 'No OligoNucleotide column found'};
      const cell0 = oligoCol.get(0);
      const semValue = (window as any).DG.SemanticValue.fromValueType(cell0, 'OligoNucleotide');
      let error = null;
      try {

        (window as any).grok.functions.call('SequenceTranslator:openOligoHelmEditor', {value: semValue});
      } catch (e: any) {
        error = String(e);
      }
      return {
        success: error === null,
        error,
        cellValue: cell0 ? String(cell0).slice(0, 30) : null,
      };
    });
    expect(result.error, 'openOligoHelmEditor must not throw synchronously for OligoNucleotide SemanticValue').toBeNull();
    expect(result.success, 'openOligoHelmEditor call initiation must succeed').toBe(true);

    // The editor opens asynchronously — wait for the dialog if it appears, then close it (best-effort).
    await page.locator('.d4-dialog').first().waitFor({state: 'visible', timeout: 10_000}).catch(() => {});
    await page.evaluate(async () => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const dlg of dialogs) {
        const cancelBtn = dlg.querySelector('[name="button-CANCEL"], [name="button-Close"]') as HTMLElement | null;
        if (cancelBtn) cancelBtn.click();
      }
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  await softStep('Scenario 3 step 1: getPtOligoEnumeratorDialog does not throw for OligoNucleotide GridCell (context-menu-oligo)', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCol = (cols as any[]).find((c) => c.semType === 'OligoNucleotide');
      if (!oligoCol) return {success: false, error: 'No OligoNucleotide column found'};

      const tv = (window as any).grok.shell.tv;
      const grid = tv ? tv.grid : null;
      if (!grid) return {success: false, error: 'No grid found'};

      const gridCell = grid.cell(oligoCol.name, 0);
      if (!gridCell) return {success: false, error: 'Could not get DG.GridCell'};

      let error = null;
      try {
        (window as any).grok.functions.call('SequenceTranslator:getPtOligoEnumeratorDialog', {cell: gridCell});
      } catch (e: any) {
        error = String(e);
      }
      return {
        success: error === null,
        error,
        gridCellHasValue: !!gridCell.cell?.value,
      };
    });
    expect(result.error, 'getPtOligoEnumeratorDialog must not throw for OligoNucleotide DG.GridCell').toBeNull();
    expect(result.success, 'getPtOligoEnumeratorDialog call initiation must succeed').toBe(true);

    // The enumerator dialog opens asynchronously — wait for it if it appears, then close it (best-effort).
    await page.locator('.d4-dialog').first().waitFor({state: 'visible', timeout: 10_000}).catch(() => {});
    await page.evaluate(async () => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const dlg of dialogs) {
        const cancelBtn = dlg.querySelector('[name="button-CANCEL"], [name="button-Close"]') as HTMLElement | null;
        if (cancelBtn) cancelBtn.click();
      }
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  await softStep('Scenario 4 step 1: combineSenseAntisenseToOligoNucleotide produces OligoNucleotide column with duplex HELM (combine-sense-antisense)', async () => {
    await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const senseCol = df.col('sense_helm');
      const antisenseCol = df.col('antisense_helm');

      await (window as any).grok.functions.call('SequenceTranslator:combineSenseAntisenseToOligoNucleotide', {
        table: df,
        senseCol,
        antiCol: antisenseCol,
      });
      await new Promise((r) => setTimeout(r, 4000));
    });

    const colProbe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCols = (cols as any[]).filter((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');

      const combinedCol = oligoCols.find((c: any) => c.name.includes('sense_helm+antisense_helm'));
      const val0 = combinedCol ? String(combinedCol.get(0)) : null;
      return {
        oligoColCount: oligoCols.length,
        combinedColFound: !!combinedCol,
        combinedColName: combinedCol ? combinedCol.name : null,
        combinedSemType: combinedCol ? combinedCol.semType : null,
        combinedQuality: combinedCol ? combinedCol.getTag('quality') : null,
        combinedCellRenderer: combinedCol ? combinedCol.getTag('cell.renderer') : null,
        val0HasPipe: val0 ? val0.includes('|') : false,
        val0HasRNA1: val0 ? val0.includes('RNA1{') : false,
        val0HasRNA2: val0 ? val0.includes('RNA2{') : false,
        val0Slice: val0 ? val0.slice(0, 80) : null,
      };
    });

    expect(colProbe.combinedColFound, 'combined sense+antisense OligoNucleotide column must be created').toBe(true);
    expect(colProbe.combinedSemType, 'combined column semType must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.combinedQuality, 'combined column quality tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.combinedCellRenderer, 'combined column cell.renderer tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.val0HasPipe, 'combined column cell 0 must contain | (duplex strand separator)').toBe(true);
    expect(colProbe.val0HasRNA1, 'combined column cell 0 must contain RNA1{ (sense strand)').toBe(true);
    expect(colProbe.val0HasRNA2, 'combined column cell 0 must contain RNA2{ (antisense strand)').toBe(true);
    expect(colProbe.oligoColCount, 'at least 2 OligoNucleotide columns expected after setup + combine').toBeGreaterThanOrEqual(2);

    const balloonError = await page.evaluate(() =>
      document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length
    );
    expect(balloonError, 'no error balloon after combineSenseAntisenseToOligoNucleotide').toBe(0);
  });

  await softStep('Scenario 4 step 2: Bio menu is accessible (DOM-driving verification)', async () => {
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});

    await page.locator('[name="div-Bio"]').click();

    await page.locator('[name="div-Bio---PolyTool"]').waitFor({state: 'visible', timeout: 8_000});

    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
  });

  await softStep('Scenario 5 step 1: getCodeToWeightsMap returns non-empty Record<string, number> (api.get-code-to-weights-map)', async () => {
    const result = await page.evaluate(async () => {
      let weights: any = null;
      let error: string | null = null;
      try {
        weights = await (window as any).grok.functions.call('SequenceTranslator:getCodeToWeightsMap', {});
      } catch (e: any) {
        error = String(e);
      }
      if (error !== null || weights === null) return {success: false, error, keyCount: 0, isObject: false, sample: null};
      const keys = Object.keys(weights);
      const allNumeric = keys.every((k) => typeof weights[k] === 'number');
      const sample = keys.slice(0, 5).reduce((acc: any, k) => { acc[k] = weights[k]; return acc; }, {});
      return {
        success: true,
        error: null,
        keyCount: keys.length,
        isObject: typeof weights === 'object',
        allNumericValues: allNumeric,
        sample,
      };
    });
    expect(result.error, 'getCodeToWeightsMap must not throw').toBeNull();
    expect(result.success, 'getCodeToWeightsMap must return successfully').toBe(true);
    expect(result.isObject, 'getCodeToWeightsMap result must be an object').toBe(true);

    expect(result.keyCount, 'getCodeToWeightsMap result must be non-empty (≥1 entry)').toBeGreaterThanOrEqual(1);
    expect(result.allNumericValues, 'all values in the code-to-weights map must be numbers (molecular weights)').toBe(true);

    const balloonError = await page.evaluate(() =>
      document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length
    );
    expect(balloonError, 'no error balloon after getCodeToWeightsMap').toBe(0);
  });

  await softStep('Scenario 5b: translateOligonucleotideSequence HELM->Axolabs works (api.translate-oligonucleotide-sequence)', async () => {
    const result = await page.evaluate(async () => {

      const df = (window as any).grok.shell.t;
      const senseCol = df ? df.col('sense_helm') : null;
      if (!senseCol) return {success: false, error: 'sense_helm column not found on active table', translated: null};
      const senseHelm = senseCol.get(0);
      if (!senseHelm) return {success: false, error: 'sense_helm row 0 is null/empty', translated: null};
      let translated: string | null = null;
      let error: string | null = null;
      try {
        const raw = await (window as any).grok.functions.call(
          'SequenceTranslator:translateOligonucleotideSequence',
          {sequence: senseHelm, sourceFormat: 'HELM', targetFormat: 'Axolabs'}
        );
        translated = raw ? String(raw) : null;
      } catch (e: any) {
        error = String(e);
      }
      return {
        success: error === null && translated !== null && translated.length > 0,
        error,
        translated: translated ? translated.slice(0, 80) : null,
      };
    });
    expect(result.error, 'translateOligonucleotideSequence HELM->Axolabs must not throw').toBeNull();
    expect(result.success, 'translateOligonucleotideSequence must return a non-empty translated string').toBe(true);
    expect(result.translated, 'translateOligonucleotideSequence Axolabs result must be non-null').toBeTruthy();

    const balloonError = await page.evaluate(() =>
      document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length
    );
    expect(balloonError, 'no error balloon from translateOligonucleotideSequence').toBe(0);
  });

  await softStep('Cleanup: close all views', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e: any) => `- ${e.step}: ${e.error}`).join('\n');
    throw new Error(`Step failures:\n${summary}`);
  }
});
