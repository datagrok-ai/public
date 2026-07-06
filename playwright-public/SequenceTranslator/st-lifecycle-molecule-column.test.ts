/* ---
sub_features_covered: [sequencetranslator.detectors.context-menu-molecule, sequencetranslator.polytool.enumerate-markush-top-menu, sequencetranslator.polytool.get-pt-chem-enumerator-dialog]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Datagrok hides (display:none) rather than detaches a dialog on close, so the
// repeated Markush opens across scenarios pile up stale [name="dialog-Markush-Enumerator"]
// nodes. That trips strict-mode locators and makes first-match querySelector probes read a
// stale dialog. Keep only the newest (last in DOM = the just-opened one) and drop the rest.
async function dropStaleMarkush(page: Page): Promise<void> {
  await page.evaluate(() => {
    const dlgs = Array.from(document.querySelectorAll('[name="dialog-Markush-Enumerator"]'));
    for (let i = 0; i < dlgs.length - 1; i++) dlgs[i].remove();
  });
}

async function clickChemTransformItem(page: Page, itemName: string): Promise<void> {
  await page.locator('[name="div-Chem"]').waitFor({state: 'visible', timeout: 15_000});
  await page.locator('[name="div-Chem"]').click();
  await page.waitForTimeout(400);

  await page.evaluate(async () => {
    const transform = document.querySelector('[name="div-Chem---Transform"]');
    if (!transform) throw new Error('[name="div-Chem---Transform"] not found after Chem click');
    const rect = transform.getBoundingClientRect();
    if (rect.width === 0)
      throw new Error('[name="div-Chem---Transform"] has width=0 after Chem click (menu not open)');
    const cx = rect.left + 10;
    const cy = rect.top + 5;
    const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
    for (const evType of eventSeq) {
      const isPointer = evType.startsWith('pointer');
      const evt = isPointer
        ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
        : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
      transform.dispatchEvent(evt);
      await new Promise((r) => setTimeout(r, 100));
    }
    await new Promise((r) => setTimeout(r, 500));
  });

  await page.evaluate(async (name) => {
    const item = document.querySelector(`[name="${name}"]`);
    if (!item) throw new Error(`[name="${name}"] not found in Transform submenu`);
    (item as HTMLElement).click();
  }, itemName);
}

async function importCoresFromTable(page: Page, tableName: string, colName: string): Promise<void> {

  await page.evaluate(() => {
    const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
    if (!dlg) throw new Error('dialog-Markush-Enumerator not found for core import');
    const importBtns = Array.from(dlg.querySelectorAll('[name="button-↓-Import…"]'));
    if (importBtns.length === 0) throw new Error('No import buttons found in Markush Enumerator dialog');
    (importBtns[0] as HTMLElement).click();
  });
  await page.locator('[name="dialog-Import-Cores"]').waitFor({state: 'visible', timeout: 10_000});

  await page.evaluate(async (tblName) => {
    const importDlg = document.querySelector('[name="dialog-Import-Cores"]');
    if (!importDlg) throw new Error('dialog-Import-Cores not found');
    const tableSel = importDlg.querySelector('[name="input-Table"]') as HTMLSelectElement;
    if (!tableSel) throw new Error('input-Table not found in dialog-Import-Cores');
    tableSel.value = tblName;
    tableSel.dispatchEvent(new Event('change', {bubbles: true}));
    tableSel.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 800));
  }, tableName);

  await page.evaluate(async (col) => {
    const importDlg = document.querySelector('[name="dialog-Import-Cores"]');
    if (!importDlg) throw new Error('dialog-Import-Cores not found');
    const colSel = importDlg.querySelector('[name="input-Column"]') as HTMLSelectElement;
    if (!colSel) throw new Error('input-Column not found in dialog-Import-Cores');
    colSel.value = col;
    colSel.dispatchEvent(new Event('change', {bubbles: true}));
    colSel.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 400));
  }, colName);

  await page.locator('[name="dialog-Import-Cores"] [name="button-OK"]').click();
  await page.locator('[name="dialog-Import-Cores"]').waitFor({state: 'hidden', timeout: 10_000});
  await page.waitForTimeout(500);
}

async function importRGroupsFromTable(page: Page, tableName: string, colName: string, targetRNum: number): Promise<void> {

  await page.evaluate(() => {
    const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
    if (!dlg) throw new Error('dialog-Markush-Enumerator not found for R-group import');
    const importBtns = Array.from(dlg.querySelectorAll('[name="button-↓-Import…"]'));
    if (importBtns.length < 2) throw new Error('Second import button not found in Markush Enumerator dialog');
    (importBtns[1] as HTMLElement).click();
  });
  await page.locator('[name="dialog-Import-R-Groups"]').waitFor({state: 'visible', timeout: 10_000});

  await page.evaluate(async (tblName) => {
    const importDlg = document.querySelector('[name="dialog-Import-R-Groups"]');
    if (!importDlg) throw new Error('dialog-Import-R-Groups not found');
    const tableSel = importDlg.querySelector('[name="input-Table"]') as HTMLSelectElement;
    if (!tableSel) throw new Error('input-Table not found in dialog-Import-R-Groups');
    tableSel.value = tblName;
    tableSel.dispatchEvent(new Event('change', {bubbles: true}));
    tableSel.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 800));
  }, tableName);

  await page.evaluate(async (col) => {
    const importDlg = document.querySelector('[name="dialog-Import-R-Groups"]');
    if (!importDlg) throw new Error('dialog-Import-R-Groups not found');
    const colSel = importDlg.querySelector('[name="input-Column"]') as HTMLSelectElement;
    if (!colSel) throw new Error('input-Column not found in dialog-Import-R-Groups');
    colSel.value = col;
    colSel.dispatchEvent(new Event('change', {bubbles: true}));
    colSel.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 400));
  }, colName);

  await page.evaluate(async (rNum) => {
    const importDlg = document.querySelector('[name="dialog-Import-R-Groups"]');
    if (!importDlg) throw new Error('dialog-Import-R-Groups not found');
    const targetInput = importDlg.querySelector('[name="input-Target-R#"]') as HTMLInputElement;
    if (!targetInput) throw new Error('input-Target-R# not found in dialog-Import-R-Groups');
    targetInput.value = String(rNum);
    targetInput.dispatchEvent(new Event('change', {bubbles: true}));
    targetInput.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 300));
  }, targetRNum);

  await page.locator('[name="dialog-Import-R-Groups"] [name="button-OK"]').click();
  await page.locator('[name="dialog-Import-R-Groups"]').waitFor({state: 'hidden', timeout: 10_000});
  await page.waitForTimeout(500);
}

test('SequenceTranslator — Molecule column lifecycle: Markush Enumeration top-menu and context-menu dialog', async ({page}) => {
  // RDKit module load + several Cartesian/Zip Markush enumerations (chem-heavy); 5 min cap.
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();

    const df = (window as any).DG.DataFrame.fromCsv('smiles\nCC([*:1])C([*:2])=O\nFC([*:1])C([*:2])=O\nClC([*:1])C([*:2])=O');
    (window as any).grok.shell.addTableView(df);

    df.col('smiles').semType = 'Molecule';

    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 8000);
    });

    const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
    const hasMolecule = (cols as any[]).some((c) => c.semType === 'Molecule');
    if (hasMolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 3000));
    }

    const coresDf = await (window as any).grok.dapi.files.readCsv('System:AppData/SequenceTranslator/tests/chem_enum_cores.csv');
    (window as any).grok.shell.addTableView(coresDf);
    await new Promise((r) => setTimeout(r, 500));

    const rgroupsDf = await (window as any).grok.dapi.files.readCsv('System:AppData/SequenceTranslator/tests/chem_enum_rgroups.csv');
    (window as any).grok.shell.addTableView(rgroupsDf);
    await new Promise((r) => setTimeout(r, 500));

    await (window as any).grok.data.detectSemanticTypes(coresDf);
    await (window as any).grok.data.detectSemanticTypes(rgroupsDf);
    await new Promise((r) => setTimeout(r, 3000));

    const mainTv = Array.from((window as any).grok.shell.tableViews)
      .find((tv: any) => tv.dataFrame.name === 'Table');
    if (mainTv) (window as any).grok.shell.tv = mainTv;
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await page.waitForFunction(() => {
    return (window as any).DG.Func.find({
      name: 'getPtChemEnumeratorDialog',
      package: 'SequenceTranslator',
    }).length > 0;
  }, null, {timeout: 60_000});

  await page.evaluate(async () => {
    await (window as any).grok.functions.call('SequenceTranslator:init', {});
    await new Promise((r) => setTimeout(r, 1000));
  });

  await page.waitForFunction(() => {
    return (window as any).DG.Func.find({
      package: 'Chem',
      name: 'getRdKitModule',
    }).length > 0;
  }, null, {timeout: 60_000});

  await page.waitForTimeout(2000);

  await softStep('Scenario 1 step 1: synthetic Markush SMILES DataFrame opens with Molecule semType column', async () => {
    const probe = await page.evaluate(() => {
      const tables = (window as any).grok.shell.tableViews;
      const mainTv = Array.from(tables).find((tv: any) => tv.dataFrame.name === 'Table');
      if (!mainTv) return {found: false, rowCount: -1, semType: null, firstValue: null};
      const df = (mainTv as any).dataFrame;
      const smilesCol = df.col('smiles');
      return {
        found: true,
        rowCount: df.rowCount,
        semType: smilesCol ? smilesCol.semType : null,
        firstValue: smilesCol ? String(smilesCol.get(0)).slice(0, 40) : null,
      };
    });
    expect(probe.found, 'main smiles table must be present').toBe(true);
    expect(probe.rowCount, 'Markush SMILES DataFrame must have 3 rows').toBe(3);
    expect(probe.semType, 'smiles column must have Molecule semType').toBe('Molecule');
    expect(probe.firstValue, 'first row must contain [*:1] R-group label').toContain('[*:1]');
  });

  await softStep('Scenario 1 step 2: Chem | Transform | Markush Enumeration... opens Markush Enumerator dialog', async () => {
    await clickChemTransformItem(page, 'div-Chem---Transform---Markush-Enumeration...');
    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 20_000});

    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
      if (!dlg) return {found: false, title: null, enumTypePresent: false, cancelPresent: false, okPresent: false, enumOptions: []};
      const title = dlg.querySelector('.d4-dialog-title')?.textContent?.trim() ?? null;
      const enumTypeSelect = dlg.querySelector('[name="input-Enumerator-type"]');
      const options = enumTypeSelect ? Array.from((enumTypeSelect as HTMLSelectElement).options).map((o) => o.value) : [];
      return {
        found: true,
        title,
        enumTypePresent: !!enumTypeSelect,
        enumOptions: options,
        cancelPresent: !!(dlg.querySelector('[name="button-CANCEL"]')),
        okPresent: !!(dlg.querySelector('[name="button-OK"]')),
      };
    });

    expect(dlgProbe.found, 'Markush Enumerator dialog must open').toBe(true);
    expect(dlgProbe.title, 'dialog title must be "Markush Enumerator"').toBe('Markush Enumerator');
    expect(dlgProbe.enumTypePresent, 'Enumerator type selector must be present').toBe(true);
    expect(dlgProbe.enumOptions, 'Enumerator type options must include Zip and Cartesian').toEqual(expect.arrayContaining(['Zip', 'Cartesian']));
    expect(dlgProbe.cancelPresent, 'CANCEL button must be present').toBe(true);
    expect(dlgProbe.okPresent, 'OK button must be present').toBe(true);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 10_000});

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Markush Enumerator CANCEL').toBe(0);
  });

  await softStep('Scenario 1 steps 3-5: enumerate Markush structures (Cartesian) and verify result table', async () => {

    await page.evaluate(async () => {
      await (window as any).grok.functions.call('SequenceTranslator:chemEnumerateMarkushTopMenu', {});
      await new Promise((r) => setTimeout(r, 1500));
    });
    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    await importCoresFromTable(page, 'Table (2)', 'Core');

    await importRGroupsFromTable(page, 'Table (3)', 'R1', 1);

    await importRGroupsFromTable(page, 'Table (3)', 'R2', 2);

    await importRGroupsFromTable(page, 'Table (3)', 'R3', 3);

    const okEnabled = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
      if (!dlg) return false;
      const okBtn = dlg.querySelector('[name="button-OK"]') as HTMLButtonElement;
      return okBtn ? !okBtn.disabled : false;
    });
    expect(okEnabled, 'OK button must be enabled after importing cores and R-groups').toBe(true);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]').selectOption('Cartesian');

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 30_000});

    await page.waitForTimeout(3000);

    const resultProbe = await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const enumTv = views.find((tv: any) => tv.dataFrame.name === 'Chem Enumeration');
      if (!enumTv) return {found: false, rowCount: -1, enumColPresent: false, nullCount: -1, residualRLabels: -1, semType: null};
      const df = (enumTv as any).dataFrame;
      const enumCol = df.col('Enumerated');
      if (!enumCol) return {found: true, rowCount: df.rowCount, enumColPresent: false, nullCount: -1, residualRLabels: -1, semType: null};
      let nullCount = 0;
      let residualRLabels = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const v = enumCol.get(i);
        if (v == null || String(v).trim() === '') nullCount++;
        else if (String(v).includes('[*:')) residualRLabels++;
      }
      return {
        found: true,
        rowCount: df.rowCount,
        enumColPresent: true,
        nullCount,
        residualRLabels,
        semType: enumCol.semType,
      };
    });

    expect(resultProbe.found, 'Chem Enumeration result table must be present after enumeration').toBe(true);
    expect(resultProbe.enumColPresent, '"Enumerated" column must exist in result table').toBe(true);
    expect(resultProbe.rowCount, 'result table must contain at least 1 enumerated molecule').toBeGreaterThan(0);
    expect(resultProbe.nullCount, 'no null values in Enumerated column').toBe(0);
    expect(resultProbe.residualRLabels, 'no residual R-labels in enumerated SMILES (all valid molecules)').toBe(0);
    expect(resultProbe.semType, 'Enumerated column must have Molecule semType').toBe('Molecule');

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Cartesian enumeration').toBe(0);
  });

  await softStep('Scenario 1 step seeding check: chemEnumerateMarkushTopMenu opens dialog (seeds from current cell)', async () => {

    await page.evaluate(async () => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const mainTv = views.find((tv: any) => tv.dataFrame.name === 'Table');
      if (mainTv) (window as any).grok.shell.tv = mainTv;
    });
    await page.waitForTimeout(500);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.dataFrame.currentRow = 0;
      await (window as any).grok.functions.call('SequenceTranslator:chemEnumerateMarkushTopMenu', {});
      await new Promise((r) => setTimeout(r, 1500));
    });

    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    const dlgTitle = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
      return dlg?.querySelector('.d4-dialog-title')?.textContent?.trim() ?? null;
    });
    expect(dlgTitle, 'Markush Enumerator dialog must open from top-menu function call').toBe('Markush Enumerator');

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 10_000});
  });

  await softStep('Scenario 2 step 1: getPtChemEnumeratorDialog opens Markush Enumerator seeded from Molecule cell', async () => {

    await page.evaluate(async () => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const mainTv = views.find((tv: any) => tv.dataFrame.name === 'Table');
      if (mainTv) (window as any).grok.shell.tv = mainTv;
    });
    await page.waitForTimeout(300);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;

      const cell = tv.dataFrame.cell(0, 'smiles');
      await (window as any).grok.functions.call('SequenceTranslator:getPtChemEnumeratorDialog', {cell});
      await new Promise((r) => setTimeout(r, 1500));
    });

    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
      if (!dlg) return {found: false, title: null};
      return {
        found: true,
        title: dlg.querySelector('.d4-dialog-title')?.textContent?.trim() ?? null,
      };
    });

    expect(dlgProbe.found, 'Markush Enumerator dialog must open from getPtChemEnumeratorDialog (context-menu path)').toBe(true);
    expect(dlgProbe.title, 'context-menu dialog title must be "Markush Enumerator"').toBe('Markush Enumerator');

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 10_000});

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after context-menu dialog CANCEL').toBe(0);
  });

  await softStep('Scenario 2 steps 2-4: getPtChemEnumeratorDialog + enumerate produces valid result column', async () => {

    await page.evaluate(async () => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const mainTv = views.find((tv: any) => tv.dataFrame.name === 'Table');
      if (mainTv) (window as any).grok.shell.tv = mainTv;
      const tv = (window as any).grok.shell.tv;

      const cell = tv.dataFrame.cell(0, 'smiles');
      await (window as any).grok.functions.call('SequenceTranslator:getPtChemEnumeratorDialog', {cell});
      await new Promise((r) => setTimeout(r, 1500));
    });
    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    await importCoresFromTable(page, 'Table (2)', 'Core');
    await importRGroupsFromTable(page, 'Table (3)', 'R1', 1);
    await importRGroupsFromTable(page, 'Table (3)', 'R2', 2);
    await importRGroupsFromTable(page, 'Table (3)', 'R3', 3);

    const okEnabled = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Markush-Enumerator"]');
      if (!dlg) return false;
      const btn = dlg.querySelector('[name="button-OK"]') as HTMLButtonElement;
      return btn ? !btn.disabled : false;
    });
    expect(okEnabled, 'OK must be enabled after core+R-group import (context-menu path)').toBe(true);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]').selectOption('Cartesian');
    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 30_000});
    await page.waitForTimeout(3000);

    const resultProbe = await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const enumViews = views.filter((tv: any) => tv.dataFrame.name === 'Chem Enumeration');

      if (enumViews.length === 0) return {found: false, rowCount: -1};
      const df = (enumViews[enumViews.length - 1] as any).dataFrame;
      const enumCol = df.col('Enumerated');
      let residualR = 0;
      if (enumCol) {
        for (let i = 0; i < df.rowCount; i++) {
          const v = enumCol.get(i);
          if (v && String(v).includes('[*:')) residualR++;
        }
      }
      return {
        found: true,
        rowCount: df.rowCount,
        enumColPresent: !!enumCol,
        residualRLabels: residualR,
        semType: enumCol ? enumCol.semType : null,
      };
    });

    expect(resultProbe.found, 'Chem Enumeration result must be produced from context-menu path').toBe(true);
    expect(resultProbe.rowCount, 'context-menu enumeration must produce at least 1 result row').toBeGreaterThan(0);
    expect(resultProbe.enumColPresent, 'Enumerated column must be present in context-menu result').toBe(true);
    expect(resultProbe.residualRLabels, 'no residual R-labels in context-menu enumeration results').toBe(0);

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after context-menu enumeration').toBe(0);
  });

  await softStep('Scenario 2: detectors context-menu dispatch works for Molecule cells at rows 1 and 2', async () => {
    for (const rowIdx of [1, 2]) {
      await page.evaluate(async (idx) => {
        const views = Array.from((window as any).grok.shell.tableViews);
        const mainTv = views.find((tv: any) => tv.dataFrame.name === 'Table');
        if (mainTv) (window as any).grok.shell.tv = mainTv;
        const tv = (window as any).grok.shell.tv;

        const cell = tv.dataFrame.cell(idx, 'smiles');
        await (window as any).grok.functions.call('SequenceTranslator:getPtChemEnumeratorDialog', {cell});
        await new Promise((r) => setTimeout(r, 1500));
      }, rowIdx);

      await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});
      await page.locator('[name="dialog-Markush-Enumerator"] [name="button-CANCEL"]').click();
      await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 10_000});
    }

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after multi-cell context-menu dispatch').toBe(0);
  });

  await softStep('Scenario 3 steps 1-7: Zip mode produces fewer results than Cartesian mode', async () => {

    await page.evaluate(async () => {

      const views = Array.from((window as any).grok.shell.tableViews);
      for (const tv of views) {
        if ((tv as any).dataFrame.name === 'Chem Enumeration') {
          await (window as any).grok.shell.closeTable((tv as any).dataFrame);
        }
      }
      await new Promise((r) => setTimeout(r, 500));

      const coreDf = (window as any).DG.DataFrame.fromCsv('ZipCore\nC([*:1])N([*:2])');
      (window as any).grok.shell.addTableView(coreDf);
      await (window as any).grok.data.detectSemanticTypes(coreDf);
      await new Promise((r) => setTimeout(r, 1000));

      const r1Df = (window as any).DG.DataFrame.fromCsv('ZipR1\nOCC[*:1]\nNCC[*:1]');
      (window as any).grok.shell.addTableView(r1Df);
      await (window as any).grok.data.detectSemanticTypes(r1Df);
      await new Promise((r) => setTimeout(r, 1000));

      const r2Df = (window as any).DG.DataFrame.fromCsv('ZipR2\nOCC[*:2]\nNCC[*:2]');
      (window as any).grok.shell.addTableView(r2Df);
      await (window as any).grok.data.detectSemanticTypes(r2Df);
      await new Promise((r) => setTimeout(r, 1000));
    });
    await page.waitForTimeout(1000);

    const tableNames = await page.evaluate(() => {
      return Array.from((window as any).grok.shell.tableViews).map((tv: any) => ({
        name: tv.dataFrame.name,
        cols: tv.dataFrame.columns.names(),
        rows: tv.dataFrame.rowCount,
      }));
    });

    const zipCoreTableName = tableNames.find((t: any) => t.cols.includes('ZipCore'))?.name;
    const zipR1TableName = tableNames.find((t: any) => t.cols.includes('ZipR1'))?.name;
    const zipR2TableName = tableNames.find((t: any) => t.cols.includes('ZipR2'))?.name;

    expect(zipCoreTableName, 'ZipCore table must be open').toBeTruthy();
    expect(zipR1TableName, 'ZipR1 table must be open').toBeTruthy();
    expect(zipR2TableName, 'ZipR2 table must be open').toBeTruthy();

    await page.evaluate(async () => {
      await (window as any).grok.functions.call('SequenceTranslator:chemEnumerateMarkushTopMenu', {});
      await new Promise((r) => setTimeout(r, 1500));
    });
    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    await importCoresFromTable(page, zipCoreTableName as string, 'ZipCore');
    await importRGroupsFromTable(page, zipR1TableName as string, 'ZipR1', 1);
    await importRGroupsFromTable(page, zipR2TableName as string, 'ZipR2', 2);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]').selectOption('Zip');

    const zipMode = await page.evaluate(() => {
      const sel = document.querySelector('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]') as HTMLSelectElement | null;
      return sel ? sel.value : null;
    });
    expect(zipMode, 'Enumerator type must be set to Zip').toBe('Zip');

    const zipOkEnabled = await page.evaluate(() => {
      const btn = document.querySelector('[name="dialog-Markush-Enumerator"] [name="button-OK"]') as HTMLButtonElement | null;
      return btn ? !btn.disabled : false;
    });
    expect(zipOkEnabled, 'OK button must be enabled for Zip mode (equal-length R-groups)').toBe(true);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 30_000});
    await page.waitForTimeout(3000);

    const zipResult = await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const enumViews = views.filter((tv: any) => tv.dataFrame.name === 'Chem Enumeration');
      if (enumViews.length === 0) return {found: false, rowCount: -1};
      const df = (enumViews[enumViews.length - 1] as any).dataFrame;
      const enumCol = df.col('Enumerated');
      let residualR = 0;
      if (enumCol) {
        for (let i = 0; i < df.rowCount; i++) {
          const v = enumCol.get(i);
          if (v && String(v).includes('[*:')) residualR++;
        }
      }
      return {found: true, rowCount: df.rowCount, enumColPresent: !!enumCol, residualRLabels: residualR};
    });

    expect(zipResult.found, 'Chem Enumeration result must be present after Zip enumeration').toBe(true);
    expect(zipResult.rowCount, 'Zip mode must produce at least 1 result').toBeGreaterThan(0);
    expect(zipResult.residualRLabels, 'no residual R-labels in Zip mode results').toBe(0);

    await page.evaluate(async () => {

      const views = Array.from((window as any).grok.shell.tableViews);
      for (const tv of views) {
        if ((tv as any).dataFrame.name === 'Chem Enumeration') {
          await (window as any).grok.shell.closeTable((tv as any).dataFrame);
        }
      }
      await new Promise((r) => setTimeout(r, 300));
    });

    await page.evaluate(async () => {
      await (window as any).grok.functions.call('SequenceTranslator:chemEnumerateMarkushTopMenu', {});
      await new Promise((r) => setTimeout(r, 1500));
    });
    await dropStaleMarkush(page);
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'visible', timeout: 15_000});

    await importCoresFromTable(page, zipCoreTableName as string, 'ZipCore');
    await importRGroupsFromTable(page, zipR1TableName as string, 'ZipR1', 1);
    await importRGroupsFromTable(page, zipR2TableName as string, 'ZipR2', 2);

    await page.locator('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]').selectOption('Cartesian');

    const cartesianMode = await page.evaluate(() => {
      const sel = document.querySelector('[name="dialog-Markush-Enumerator"] [name="input-Enumerator-type"]') as HTMLSelectElement | null;
      return sel ? sel.value : null;
    });
    expect(cartesianMode, 'Enumerator type must be set to Cartesian').toBe('Cartesian');

    await page.locator('[name="dialog-Markush-Enumerator"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Markush-Enumerator"]').waitFor({state: 'hidden', timeout: 30_000});
    await page.waitForTimeout(3000);

    const cartesianResult = await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.tableViews);
      const enumViews = views.filter((tv: any) => tv.dataFrame.name === 'Chem Enumeration');
      if (enumViews.length === 0) return {found: false, rowCount: -1};
      const df = (enumViews[enumViews.length - 1] as any).dataFrame;
      const enumCol = df.col('Enumerated');
      let residualR = 0;
      if (enumCol) {
        for (let i = 0; i < df.rowCount; i++) {
          const v = enumCol.get(i);
          if (v && String(v).includes('[*:')) residualR++;
        }
      }
      return {found: true, rowCount: df.rowCount, enumColPresent: !!enumCol, residualRLabels: residualR};
    });

    expect(cartesianResult.found, 'Chem Enumeration result must be present after Cartesian enumeration').toBe(true);
    expect(cartesianResult.rowCount, 'Cartesian mode must produce at least 1 result').toBeGreaterThan(0);
    expect(cartesianResult.residualRLabels, 'no residual R-labels in Cartesian mode results').toBe(0);

    expect(cartesianResult.rowCount, 'Cartesian mode must produce more results than Zip mode for same inputs').toBeGreaterThan(zipResult.rowCount);

    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Zip vs Cartesian mode comparison').toBe(0);
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
