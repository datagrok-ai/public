import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const SIRNA_PATH = 'System:AppData/SequenceTranslator/samples/sirna-demo.csv';
const CYCLIZED_PATH = 'System:AppData/SequenceTranslator/samples/cyclized.csv';

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

async function expandDgSubmenu(page: Page, itemSelector: string): Promise<void> {
  await page.evaluate(async (sel) => {
    const item = document.querySelector(sel) as HTMLElement | null;
    if (!item) throw new Error(`expandDgSubmenu: item not found: ${sel}`);
    const rect = item.getBoundingClientRect();
    const cx = rect.left + 10;
    const cy = rect.top + 5;
    const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
    for (const evType of eventSeq) {
      const isPointer = evType.startsWith('pointer');
      const evt = isPointer
        ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
        : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
      item.dispatchEvent(evt);
      await new Promise((r) => setTimeout(r, 100));
    }
    await new Promise((r) => setTimeout(r, 500));
  }, itemSelector);
}

async function clickBioPolyToolItem(page: Page, itemName: string): Promise<void> {

  await page.locator('[name="div-Bio"]').click();
  await page.waitForTimeout(400);

  await page.evaluate(async () => {
    const pt = document.querySelector('[name="div-Bio---PolyTool"]') as HTMLElement | null;
    if (!pt) throw new Error('[name="div-Bio---PolyTool"] not found after Bio click');
    const ptRect = pt.getBoundingClientRect();
    if (ptRect.width === 0)
      throw new Error('[name="div-Bio---PolyTool"] has width=0 after Bio click (menu not open)');
    const cx = ptRect.left + 10;
    const cy = ptRect.top + 5;
    const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
    for (const evType of eventSeq) {
      const isPointer = evType.startsWith('pointer');
      const evt = isPointer
        ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
        : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
      pt.dispatchEvent(evt);
      await new Promise((r) => setTimeout(r, 100));
    }
    await new Promise((r) => setTimeout(r, 500));
  });

  await page.evaluate(async (name) => {
    const item = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
    if (!item) throw new Error(`[name="${name}"] not found in PolyTool submenu`);
    item.click();
  }, itemName);
}

test('SequenceTranslator — Macromolecule HELM column lifecycle: init, converters, PolyTool, API', async ({page}) => {
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
  });

  await softStep('Scenario 1 step 1: sirna-demo.csv opens with Macromolecule HELM columns (lifecycle.init evidence)', async () => {
    const probe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const helmCols = (cols as any[]).filter((c) => c.semType === 'Macromolecule' && c.meta?.units === 'helm');
      return {
        rowCount: df.rowCount,
        helmColCount: helmCols.length,
        helmColNames: helmCols.map((c: any) => c.name),
      };
    });
    expect(probe.rowCount, 'sirna-demo.csv must open with rows').toBeGreaterThanOrEqual(10);
    expect(probe.helmColCount, 'at least 3 Macromolecule/helm columns expected (sense_helm, antisense_helm, oligo_helm)').toBeGreaterThanOrEqual(3);
  });

  await softStep('Scenario 1 step 2: convert oligo_helm HELM column to OligoNucleotide via JS API', async () => {
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
        semType: oligoCol?.semType,
        quality: oligoCol?.getTag('quality'),
        cellRenderer: oligoCol?.getTag('cell.renderer'),
        originalColIntact: !!(df.col('oligo_helm')),
      };
    });
    expect(colProbe.found, 'OligoNucleotide column must be appended after conversion').toBe(true);
    expect(colProbe.semType, 'column semType must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.quality, 'column quality tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.cellRenderer, 'column cell.renderer tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.originalColIntact, 'original oligo_helm column must not be modified').toBe(true);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after HELM to OligoNucleotide conversion').toBe(0);
  });

  await softStep('Scenario 2 step 1: combine sense_helm + antisense_helm to OligoNucleotide via JS API', async () => {
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
      const combinedCols = (cols as any[]).filter((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');
      return {
        oligoColCount: combinedCols.length,
        firstValue: combinedCols.length > 0 ? String(combinedCols[0].get(0)).slice(0, 30) : null,
      };
    });

    expect(colProbe.oligoColCount, 'at least one OligoNucleotide column from combine operation').toBeGreaterThanOrEqual(1);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after combine sense+antisense').toBe(0);
  });

  await softStep('Scenario 3 setup: open cyclized.csv (custom-notation Macromolecule column for PolyTool)', async () => {
    await page.evaluate(async (path) => {
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
    const colProbe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const seqsCol = df.col('seqs');
      return {
        rowCount: df.rowCount,
        seqsSemType: seqsCol?.semType,
        seqsUnits: seqsCol?.meta?.units,
        seqsPolyTool: seqsCol?.getTag('polytool-data-role'),
      };
    });
    expect(colProbe.rowCount, 'cyclized.csv must have rows').toBeGreaterThanOrEqual(1);
    expect(colProbe.seqsSemType, 'seqs column must detect as Macromolecule').toBe('Macromolecule');
    expect(colProbe.seqsUnits, 'seqs column must have units=custom').toBe('custom');
  });

  await softStep('Scenario 3 step 1: Bio | PolyTool | Convert... opens PolyTool Conversion dialog; CANCEL closes', async () => {

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
    expect(dlgProbe.title, 'PolyTool Conversion dialog title').toBe('PolyTool Conversion');
    expect(dlgProbe.columnInputPresent, 'Column input must be in PolyTool Conversion dialog').toBe(true);
    expect(dlgProbe.cancelPresent, 'CANCEL button must be in PolyTool Conversion dialog').toBe(true);
    await page.locator('[name="dialog-PolyTool-Conversion"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-PolyTool-Conversion"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after PolyTool Conversion CANCEL').toBe(0);
  });

  await softStep('Scenario 4 step 1: Bio | PolyTool | Enumerate HELM... opens PolyTool Helm Enumeration dialog; CANCEL closes', async () => {
    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Enumerate-HELM...');
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-PolyTool-Helm-Enumeration"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
      };
    });
    expect(dlgProbe.title, 'PolyTool Helm Enumeration dialog title').toBe('PolyTool Helm Enumeration');
    expect(dlgProbe.cancelPresent, 'CANCEL button must be in PolyTool Helm Enumeration dialog').toBe(true);
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after PolyTool Helm Enumeration CANCEL').toBe(0);
  });

  await softStep('Scenario 5 step 1: Bio | PolyTool | Combine Sequences... opens Combine Sequences dialog; CANCEL closes', async () => {
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
    expect(dlgProbe.title, 'Combine Sequences dialog title').toBe('Combine Sequences');
    expect(dlgProbe.tableInputPresent, 'Table input must be in Combine Sequences dialog').toBe(true);
    expect(dlgProbe.columnInputPresent, 'Column input must be in Combine Sequences dialog').toBe(true);
    await page.locator('[name="dialog-Combine-Sequences"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Combine-Sequences"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Combine Sequences CANCEL').toBe(0);
  });

  await softStep('Scenario 5b: enumerateSingleHelmSequence produces expected row count', async () => {
    const result = await page.evaluate(async () => {
      const helmSequence = 'PEPTIDE1{A.C.G.T}$$$$';
      const positions = [1];
      const monomerLists = [['A', 'G', 'T']];
      const df = await (window as any).grok.functions.call('SequenceTranslator:enumerateSingleHelmSequence', {
        helmSequence,
        positions,
        monomerLists,
        toAtomicLevel: false,
      });
      if (!df) return {success: false, rowCount: 0, colCount: 0};
      return {success: true, rowCount: df.rowCount, colCount: df.columns.length};
    });
    expect(result.success, 'enumerateSingleHelmSequence must return a DataFrame').toBe(true);

    expect(result.rowCount, 'enumerateSingleHelmSequence with 3 monomers must produce 3 rows').toBe(3);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon from enumerateSingleHelmSequence').toBe(0);
  });

  await softStep('Scenario 6 step 1: validateSequence returns true for valid nucleotide sequence', async () => {

    const v1 = await validateViaHelper(page, 'acgu');
    const v2 = await validateViaHelper(page, 'AGCT');
    const v3 = await validateViaHelper(page, 'notasequence12345');
    expect(v1.error, 'validateSequence must not throw for valid input').toBeNull();
    expect(v3.error, 'validateSequence must not throw for garbage input').toBeNull();
    expect(v1.returnValue, "validateSequence('acgu') must return true").toBe(true);
    expect(v2.returnValue, "validateSequence('AGCT') must return true").toBe(true);
    expect(v3.returnValue, 'validateSequence with garbage must return false').toBe(false);
  });

  await softStep('Scenario 6 step 2: translateOligonucleotideSequence converts HELM to Axolabs', async () => {
    const result = await page.evaluate(async () => {

      const senseHelm = await (window as any).grok.dapi.files.readCsv(
        'System:AppData/SequenceTranslator/samples/sirna-demo.csv'
      ).then((df: any) => df.col('sense_helm').get(0));

      let translated = null;
      let error = null;
      try {
        translated = await (window as any).grok.functions.call(
          'SequenceTranslator:translateOligonucleotideSequence',
          {sequence: senseHelm, sourceFormat: 'HELM', targetFormat: 'Axolabs'}
        );
      } catch (e: any) {
        error = String(e);
      }
      return {translated: translated ? String(translated).slice(0, 80) : null, error};
    });
    expect(result.error, 'translateOligonucleotideSequence HELM->Axolabs must not throw').toBeNull();
    expect(result.translated, 'translateOligonucleotideSequence must return a non-empty string').toBeTruthy();
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
