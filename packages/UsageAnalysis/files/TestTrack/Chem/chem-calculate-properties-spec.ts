/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;

// Selector recon-notes (class-2: derived from the registered function signatures in
// public/packages/Chem/src/package.g.ts, read 2026-08-22; added to chem.md in the same change):
//   [name="dialog-Chemical-Properties"] — title "Chemical Properties", the top-menu editor dialog of
//     `Chem | Calculate | Chemical Properties...` (package.g.ts:1057-1075, //top-menu: line).
//   Its checkbox inputs follow the caption→`[name="input-<Caption>"]` convention, captions taken
//     verbatim from the //input: lines at package.g.ts:1060-1068 — input-MW, input-HBA, input-HBD,
//     input-Log-P, input-Log-S, input-PSA, input-Rotatable-bonds, input-Stereo-centers,
//     input-Molecule-charge. Registered defaults: MW = true, the other eight false.
//   [name="dialog-Toxicity-Risks"] — title "Toxicity Risks", top-menu editor dialog of
//     `Chem | Calculate | Toxicity Risks...` (package.g.ts:1088-1100).
//   Its checkbox inputs, captions from package.g.ts:1093-1096 — input-Mutagenicity,
//     input-Tumorigenicity, input-Irritating-effects, input-Reproductive-effects. Registered
//     defaults: mutagenicity = true, the other three false.
//   [name="input-host-Molecules"] .d4-column-selector-column (both dialogs) — the molecule-column
//     host is named after the `molecules` input's caption, the same shape chem.md documents for the
//     Descriptors dialog ("Chem | Calculate | Descriptors dialog").
//   [name="dialog-To-InchI"] / [name="dialog-To-InchI-Keys"] with [name="button-OK"] — the
//     function-call dialogs of `Chem | Calculate | To InchI...` / `To InchI Keys...`; there is no
//     immediate-execution branch.

test.use(specTestOptions);

const MOL_COL = 'canonical_smiles';
// technical: closed OCL risk-level vocabulary, verbatim from Chem/src/open-chem/ocl-service/consts.ts:16-21
const RISK_LEVELS = ['Unknown', 'None', 'Low', 'High'];

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.columns.names());
}

async function setBool(page: Page, dialogName: string, inputName: string, value: boolean): Promise<void> {
  await page.locator(`[name="${dialogName}"] [name="${inputName}"]`).setChecked(value);
}

async function waitForAddedColumns(page: Page, before: string[], minNew: number, capMs: number): Promise<string[]> {
  await page.waitForFunction(({before, minNew}) =>
    grok.shell.t.columns.names().filter((n: string) => !before.includes(n)).length >= minNew,
  {before, minNew}, {timeout: capMs}).catch(() => {});
  const after = await columnNames(page);
  return after.filter((n) => !before.includes(n));
}

test('Chem: Calculate — Chemical Properties, Toxicity Risks, To InChI, To InChI Keys', async ({page}) => {
  test.setTimeout(600_000);

  const consoleErrors: string[] = [];
  const AMBIENT = [/favicon/i, /ResizeObserver loop/i, /Permissions policy violation/i,
    /Unable to find element in cloned iframe/i];
  const isAmbient = (text: string) => AMBIENT.some((re) => re.test(text));
  page.on('console', (msg) => { if (msg.type() === 'error' && !isAmbient(msg.text())) consoleErrors.push(msg.text()); });
  page.on('pageerror', (e) => { if (!isAmbient(String(e))) consoleErrors.push(String(e)); });
  const expectNoConsoleErrors = (where: string) => {
    const seen = consoleErrors.splice(0, consoleErrors.length);
    expect(seen, `JavaScript console errors fired during ${where}`).toEqual([]);
  };

  await loginToDatagrok(page);

  let baseline = 0;

  await softStep('Setup: open smiles.csv, wait for Molecule semType + Chem menu; record baseline row count', async () => {
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      const detected = new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      grok.shell.addTableView(df);
      await detected;
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30000, state: 'attached'});
    await waitForChemMenu(page);
    await waitForMolecule(page);
    baseline = await page.evaluate(() => grok.shell.t.rowCount);
    expect(baseline, 'smiles.csv baseline row count must be established and positive').toBeGreaterThan(0);
    const molTyped = await page.evaluate((c) => grok.shell.t.col(c)?.semType === 'Molecule', MOL_COL);
    expect(molTyped, `the ${MOL_COL} column must be typed as Molecule before any calculator runs`).toBe(true);
  });

  await softStep('S1.1-3: Chem → Calculate → Chemical Properties → dialog opens, MW only, OK', async () => {
    await openChemMenuItem(page, 'Chemical Properties...', {delayMs: 700});
    await page.locator('[name="dialog-Chemical-Properties"]').waitFor({timeout: 15000});
    const molHost = await page.locator('[name="dialog-Chemical-Properties"] [name="input-host-Molecules"] .d4-column-selector-column').textContent();
    expect(molHost?.trim(), 'the Chemical Properties molecule selector must auto-detect the Molecule column').toBe(MOL_COL);
    for (const b of ['input-HBA', 'input-HBD', 'input-Log-P', 'input-Log-S', 'input-PSA',
      'input-Rotatable-bonds', 'input-Stereo-centers', 'input-Molecule-charge'])
      await expect(page.locator(`[name="dialog-Chemical-Properties"] [name="${b}"]`),
        `${b} must open unchecked, so that this run is the minimal MW-only case`).not.toBeChecked();
    await expect(page.locator('[name="dialog-Chemical-Properties"] [name="input-MW"]'),
      'MW must open checked — it is the only property whose registered default is true').toBeChecked();
    await page.locator('[name="dialog-Chemical-Properties"] [name="button-OK"]').click();
  });

  await softStep('S1.4: an MW column is appended with numeric values > 0; baseline row count unchanged; no console errors', async () => {
    await page.waitForFunction(() => !!grok.shell.t.col('MW'), null, {timeout: 60000}).catch(() => {});
    const probe = await page.evaluate((baseline) => {
      const t = grok.shell.t;
      const mw = t.col('MW');
      const vals = mw ? Array.from({length: t.rowCount}, (_: any, i: number) => mw.get(i)) : [];
      const numericPositive = vals.filter((v: any) => typeof v === 'number' && !Number.isNaN(v) && v > 0).length;
      return {hasMW: !!mw, mwType: mw?.type ?? null, numericPositive, rowCount: t.rowCount, baseline};
    }, baseline);
    console.log(`S1.4 MW: type=${probe.mwType} rows=${probe.rowCount} numericPositive=${probe.numericPositive} ` +
      `baseline=${probe.baseline}`);
    expect(probe.hasMW, 'an MW column must be appended after the MW-only Chemical Properties run').toBe(true);
    expect(probe.mwType, 'the MW column must hold numeric (double) values').toBe('double');
    expect(probe.numericPositive,
      'every row of the MW column must carry a numeric value greater than zero — smiles.csv holds no unparseable molecule')
      .toBe(probe.rowCount);
    expect(probe.rowCount, 'the grid row count must be unchanged after column-append').toBe(probe.baseline);
    expectNoConsoleErrors('the Chemical Properties dialog open, OK click, and MW column-append');
  });

  await softStep('S1.6: re-open Chemical Properties, select all nine properties, OK; all nine columns appended, numeric, row count == baseline; no console errors', async () => {
    await openChemMenuItem(page, 'Chemical Properties...', {delayMs: 700});
    await page.locator('[name="dialog-Chemical-Properties"]').waitFor({timeout: 15000});
    await expect(page.locator('[name="dialog-Chemical-Properties"] [name="input-MW"]'),
      'MW must re-open already checked — its registered default is true, so the full-property run writes nothing to it').toBeChecked();
    for (const b of ['input-HBA', 'input-HBD', 'input-Log-P', 'input-Log-S', 'input-PSA',
      'input-Rotatable-bonds', 'input-Stereo-centers', 'input-Molecule-charge'])
      await setBool(page, 'dialog-Chemical-Properties', b, true);
    const before = await columnNames(page);
    await page.locator('[name="dialog-Chemical-Properties"] [name="button-OK"]').click();
    const added = await waitForAddedColumns(page, before, 9, 90000);
    const probe = await page.evaluate(({added, baseline}) => {
      const t = grok.shell.t;
      const numericCols = added.filter((n: string) => {
        const c = t.col(n);
        if (!c) return false;
        const vals = Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i));
        return vals.filter((v: any) => typeof v === 'number' && !Number.isNaN(v)).length === t.rowCount;
      });
      const stable = ['HBA', 'HBD', 'LogP', 'LogS', 'PSA', 'Rotatable bonds', 'Stereo centers', 'Molecule charge']
        .filter((n) => added.includes(n));
      const dedupMW = added.filter((n: string) => /^MW\b/.test(n) && n !== 'MW');
      return {addedCount: added.length, numericColCount: numericCols.length, stablePresent: stable,
        dedupMW, addedNames: added, rowCount: t.rowCount, baseline};
    }, {added, baseline});
    console.log(`S1.6 appended ${probe.addedCount} columns ${JSON.stringify(probe.addedNames)} over ` +
      `${probe.rowCount} rows (baseline=${probe.baseline}); numericColCount=${probe.numericColCount} ` +
      `stablePresent=${JSON.stringify(probe.stablePresent)} dedupMW=${JSON.stringify(probe.dedupMW)}`);
    expect(probe.addedCount, 'all nine property columns must be appended when the full set is selected').toBe(9);
    expect(probe.stablePresent.length,
      'the eight fixed-name property columns (HBA, HBD, LogP, LogS, PSA, Rotatable bonds, Stereo centers, Molecule charge) must all be present').toBe(8);
    expect(probe.dedupMW,
      'exactly one appended column must be MW re-calculated under a de-duplicated name (e.g. "MW (2)") — that is the ' +
      'ninth column, and it must not be the literal "MW" already appended by S1.4').toHaveLength(1);
    expect(probe.numericColCount, 'every appended property column must contain numeric values for valid molecules').toBe(9);
    expect(probe.rowCount, 'the grid row count must equal the baseline after the full property set').toBe(probe.baseline);
    expectNoConsoleErrors('the Chemical Properties re-open, full-property OK click, and nine-column append');
  });

  await softStep('S2.1-3: Chem → Calculate → Toxicity Risks → Mutagenicity only (default), OK', async () => {
    await openChemMenuItem(page, 'Toxicity Risks...', {delayMs: 700});
    await page.locator('[name="dialog-Toxicity-Risks"]').waitFor({timeout: 15000});
    const molHost = await page.locator('[name="dialog-Toxicity-Risks"] [name="input-host-Molecules"] .d4-column-selector-column').textContent();
    expect(molHost?.trim(), 'the Toxicity Risks molecule selector must auto-detect the Molecule column').toBe(MOL_COL);
    for (const b of ['input-Tumorigenicity', 'input-Irritating-effects', 'input-Reproductive-effects'])
      await expect(page.locator(`[name="dialog-Toxicity-Risks"] [name="${b}"]`),
        `${b} must open unchecked, so that this run is Mutagenicity-only`).not.toBeChecked();
    await expect(page.locator('[name="dialog-Toxicity-Risks"] [name="input-Mutagenicity"]'),
      'Mutagenicity must open checked — it is the only risk whose registered default is true').toBeChecked();
    await page.locator('[name="dialog-Toxicity-Risks"] [name="button-OK"]').click();
  });

  await softStep('S2.4: a Mutagenicity column is appended; non-empty; at least one Low/High risk level; row count unchanged', async () => {
    await page.waitForFunction(() => !!grok.shell.t.col('Mutagenicity'), null, {timeout: 60000}).catch(() => {});
    const probe = await page.evaluate(({baseline, levels}) => {
      const t = grok.shell.t;
      const col = t.col('Mutagenicity');
      const vals = col ? Array.from({length: t.rowCount}, (_: any, i: number) => col.get(i)) : [];
      const nonEmpty = vals.filter((v: any) => v !== '' && v != null).length;
      const inVocabulary = vals.filter((v: any) => levels.includes(v)).length;
      const outOfVocabulary = Array.from(new Set(vals.filter((v: any) => !levels.includes(v))));
      const assessedRisk = vals.filter((v: any) => v === 'Low' || v === 'High').length;
      return {has: !!col, colType: col?.type ?? null, nonEmpty, inVocabulary, outOfVocabulary,
        assessedRisk, rowCount: t.rowCount, baseline};
    }, {baseline, levels: RISK_LEVELS});
    console.log(`S2.4 Mutagenicity: type=${probe.colType} rows=${probe.rowCount} nonEmpty=${probe.nonEmpty} ` +
      `inVocabulary=${probe.inVocabulary} outOfVocabulary=${JSON.stringify(probe.outOfVocabulary)} ` +
      `assessedRisk=${probe.assessedRisk}`);
    expect(probe.has, 'a Mutagenicity column must be appended').toBe(true);
    expect(probe.colType, 'the Mutagenicity column must be a string column, never numeric or boolean').toBe('string');
    expect(probe.nonEmpty, 'the Mutagenicity column must be non-empty for every row').toBe(probe.rowCount);
    expect(probe.outOfVocabulary,
      `no Mutagenicity value may fall outside the closed OCL risk vocabulary ${RISK_LEVELS.join('/')}`).toEqual([]);
    expect(probe.inVocabulary,
      `every row of the Mutagenicity column must read one of ${RISK_LEVELS.join('/')}`).toBe(probe.rowCount);
    expect(probe.assessedRisk,
      'at least one row must read "Low" or "High" — "Unknown" is the failed-to-assess level and an all-"None"/"Unknown" column means the risk engine produced no result').toBeGreaterThan(0);
    expect(probe.rowCount, 'the grid row count must be unchanged').toBe(probe.baseline);
  });

  await softStep('S2.6: re-open Toxicity Risks, select all four risks, OK; all four columns appended, non-empty, Low/High present, row count == baseline', async () => {
    await openChemMenuItem(page, 'Toxicity Risks...', {delayMs: 700});
    await page.locator('[name="dialog-Toxicity-Risks"]').waitFor({timeout: 15000});
    await expect(page.locator('[name="dialog-Toxicity-Risks"] [name="input-Mutagenicity"]'),
      'Mutagenicity must re-open already checked — its registered default is true, so the full-risk run writes nothing to it').toBeChecked();
    for (const b of ['input-Tumorigenicity', 'input-Irritating-effects', 'input-Reproductive-effects'])
      await setBool(page, 'dialog-Toxicity-Risks', b, true);
    const before = await columnNames(page);
    await page.locator('[name="dialog-Toxicity-Risks"] [name="button-OK"]').click();
    const added = await waitForAddedColumns(page, before, 4, 120000);
    const probe = await page.evaluate(({added, baseline, levels}) => {
      const t = grok.shell.t;
      const detail = added.map((n: string) => {
        const c = t.col(n);
        const vals = c ? Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i)) : [];
        return {name: n, has: !!c, colType: c?.type ?? null,
          nonEmpty: vals.filter((v: any) => v !== '' && v != null).length,
          inVocabulary: vals.filter((v: any) => levels.includes(v)).length,
          outOfVocabulary: Array.from(new Set(vals.filter((v: any) => !levels.includes(v)))),
          assessedRisk: vals.filter((v: any) => v === 'Low' || v === 'High').length};
      });
      const fixed = ['Tumorigenicity', 'Irritating effects', 'Reproductive effects'].filter((n) => added.includes(n));
      const dedupMutagenicity = added.filter((n: string) => /^Mutagenicity\b/.test(n) && n !== 'Mutagenicity');
      return {addedCount: added.length, addedNames: added, detail, fixedPresent: fixed, dedupMutagenicity,
        rowCount: t.rowCount, baseline};
    }, {added, baseline, levels: RISK_LEVELS});
    console.log(`S2.6 appended ${probe.addedCount} columns ${JSON.stringify(probe.addedNames)} over ` +
      `${probe.rowCount} rows; fixedPresent=${JSON.stringify(probe.fixedPresent)} ` +
      `dedupMutagenicity=${JSON.stringify(probe.dedupMutagenicity)}; per column ${JSON.stringify(probe.detail)}`);
    expect(probe.addedCount, 'all four risk columns must be appended when the full risk set is selected').toBe(4);
    expect(probe.fixedPresent.length,
      'the three fixed-name risk columns (Tumorigenicity, Irritating effects, Reproductive effects) must all be present').toBe(3);
    expect(probe.dedupMutagenicity,
      'exactly one appended column must be Mutagenicity re-calculated under a de-duplicated name (e.g. "Mutagenicity (2)") — ' +
      'that is the fourth column, and it must not be the literal "Mutagenicity" already appended by S2.4').toHaveLength(1);
    for (const d of probe.detail) {
      expect(d.has, `the ${d.name} risk column must be appended`).toBe(true);
      expect(d.colType, `the ${d.name} column must be a string column, never numeric or boolean`).toBe('string');
      expect(d.nonEmpty, `the ${d.name} column must be non-empty for every row`).toBe(probe.rowCount);
      expect(d.outOfVocabulary,
        `no ${d.name} value may fall outside the closed OCL risk vocabulary ${RISK_LEVELS.join('/')}`).toEqual([]);
      expect(d.inVocabulary,
        `every row of the ${d.name} column must read one of ${RISK_LEVELS.join('/')}`).toBe(probe.rowCount);
      expect(d.assessedRisk,
        `the ${d.name} engine must produce at least one "Low" or "High" value — "Unknown" means it failed to assess`).toBeGreaterThan(0);
    }
    expect(probe.rowCount, 'the grid row count must equal the baseline after the full risk set').toBe(probe.baseline);
  });

  await softStep('S3.3: Chem → Calculate → To InChI → OK; an InChI column of "InChI=" strings is appended, no blank on valid, row count == baseline', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'To InchI...', {delayMs: 700});
    await page.locator('[name="dialog-To-InchI"] [name="button-OK"]').waitFor({timeout: 15000});
    await page.locator('[name="dialog-To-InchI"] [name="button-OK"]').click();
    await waitForAddedColumns(page, before, 1, 60000);
    const probe = await page.evaluate((baseline) => {
      const t = grok.shell.t;
      const col = t.col('inchi');
      const vals = col ? Array.from({length: t.rowCount}, (_: any, i: number) => col.get(i)) : [];
      const validMol = t.col('canonical_smiles');
      const validRows = Array.from({length: t.rowCount}, (_: any, i: number) => validMol.get(i))
        .map((v: any, i: number) => ({i, valid: v != null && v !== ''}));
      const startsInchi = vals.filter((v: any) => typeof v === 'string' && v.startsWith('InChI=')).length;
      const blankOnValid = validRows.filter((r: any) => r.valid && (vals[r.i] === '' || vals[r.i] == null)).length;
      return {has: !!col, startsInchi, blankOnValid, rowCount: t.rowCount, baseline};
    }, baseline);
    expect(probe.has, 'an InChI column must be appended').toBe(true);
    expect(probe.startsInchi,
      'every valid molecule row must carry a string starting with "InChI="').toBe(probe.rowCount);
    expect(probe.blankOnValid,
      'no valid molecule row may have a blank or null InChI value').toBe(0);
    expect(probe.rowCount, 'the grid row count must equal the baseline').toBe(probe.baseline);
  });

  await softStep('S4.3: Chem → Calculate → To InChI Keys → OK; a 27-char InChIKey column is appended, no blank on valid, differs from InChI, row count == baseline', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'To InchI Keys...', {delayMs: 700});
    await page.locator('[name="dialog-To-InchI-Keys"] [name="button-OK"]').waitFor({timeout: 15000});
    await page.locator('[name="dialog-To-InchI-Keys"] [name="button-OK"]').click();
    await waitForAddedColumns(page, before, 1, 60000);
    const probe = await page.evaluate((baseline) => {
      const t = grok.shell.t;
      const col = t.col('inchi_key');
      const inchi = t.col('inchi');
      const re = /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/;
      const vals = col ? Array.from({length: t.rowCount}, (_: any, i: number) => col.get(i)) : [];
      const inchiVals = inchi ? Array.from({length: t.rowCount}, (_: any, i: number) => inchi.get(i)) : [];
      const validMol = t.col('canonical_smiles');
      const validRows = Array.from({length: t.rowCount}, (_: any, i: number) => validMol.get(i));
      const matchKey = vals.filter((v: any) => typeof v === 'string' && re.test(v)).length;
      const blankOnValid = validRows.filter((v: any, i: number) => v != null && v !== '' && (vals[i] === '' || vals[i] == null)).length;
      const differsFromInchi = vals.filter((v: any, i: number) => v !== inchiVals[i]).length;
      return {has: !!col, hasInchi: !!inchi, matchKey, blankOnValid, differsFromInchi, rowCount: t.rowCount, baseline};
    }, baseline);
    console.log(`S4.3 inchi_key: rows=${probe.rowCount} baseline=${probe.baseline} matchKey=${probe.matchKey} ` +
      `blankOnValid=${probe.blankOnValid} hasInchi=${probe.hasInchi} differsFromInchi=${probe.differsFromInchi}`);
    expect(probe.has, 'an InChIKey column must be appended').toBe(true);
    expect(probe.matchKey,
      'every valid molecule row must carry a 27-char InChIKey matching [A-Z]{14}-[A-Z]{10}-[A-Z]').toBe(probe.rowCount);
    expect(probe.blankOnValid, 'no valid molecule row may have a blank InChIKey').toBe(0);
    expect(probe.hasInchi,
      'the inchi column from S3.3 must still be present — without it the InChIKey/InChI comparison below would ' +
      'be satisfied by an absent column').toBe(true);
    expect(probe.differsFromInchi,
      'InChIKey values must differ from the InChI string column values').toBe(probe.rowCount);
    expect(probe.rowCount, 'the grid row count must equal the baseline').toBe(probe.baseline);
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
