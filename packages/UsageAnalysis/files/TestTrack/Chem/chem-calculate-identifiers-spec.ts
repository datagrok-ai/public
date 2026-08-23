/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;

// Selector recon-notes (class-2: live-MCP-observed 2026-08-20 on dev via chrome-devtools MCP; not in chem.md):
//   [name="dialog-Map-Identifiers"] with [name="input-From-Source"] / [name="input-To-Source"] SELECTs
//     (source ids lowercase: smiles/chembl/pubchem/...) and [name="input-host-Ids"] molecule column —
//     Map Identifiers dialog (Chem | Calculate | Map Identifiers...). The appended column is named after
//     the toSource value verbatim (chembl, pubchem). A bare `select.value = …` assignment leaves the
//     ui.input.choice model at its default (→ wrong-named column); Playwright's selectOption drives the
//     model correctly and is what this spec uses (re-probed on dev 2026-08-22: chembl, 993/1000 resolved).
//   [name="dialog-Biochemical-Properties"] with .biochem-calc-nav-item rows (each an input[type=checkbox]
//     + span label: Chemical Properties, logD, logP, pI, pKa) and [name="input-host-Molecules"] —
//     Biochemical Properties dialog (Chem | Calculate | Biochemical Properties, no "..." suffix).
//     Chemical Properties appends MW; logP appends clogP (both 1000/1000 finite on dev 2026-08-22).
//   [name="dialog-Generate-Conformers"] with a .d4-input-molecule-canvas molecule input that opens a
//     sketcher dialog (input[placeholder*="SMILES"]) on click — Generate Conformers (Chem | Calculate |
//     Generate Conformers...). Real fill + Enter + a real click on the sketcher's OK closes the sketcher
//     and carries the molecule into the script dialog (dev 2026-08-22, probed on the script's default
//     butane: 3 conformers, finite energies; the spec drives pentane so the run cannot be the default).

test.use(specTestOptions);

const MOL_COL = 'canonical_smiles';

// technical: pentane, deliberately NOT the script's default "CCCC" (Chem/scripts/generate_conformers.py:7);
// it is canonical under RDKit, so the output `smiles` column reads back verbatim.
const MOLECULE = 'CCCCC';

// technical: CHEM_PROP_MAP names (Chem/src/open-chem/ocl-service/calculations.ts:4); getUnusedName may add " (2)".
const OCL_PROPERTY_NAMES = ['MW', 'HBA', 'HBD', 'LogP', 'LogS', 'PSA',
  'Rotatable bonds', 'Stereo centers', 'Molecule charge'];

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.columns.names());
}

async function pollForColumn(page: Page, before: string[], capMs: number): Promise<string[]> {
  await page.waitForFunction((before) =>
    grok.shell.t.columns.names().some((n: string) => !before.includes(n)),
  before, {timeout: capMs}).catch(() => {});
  const after = await columnNames(page);
  return after.filter((n) => !before.includes(n));
}

async function setMapIdentifiersSources(page: Page, from: string, to: string): Promise<void> {
  const dialog = '[name="dialog-Map-Identifiers"]';
  await page.selectOption(`${dialog} [name="input-From-Source"]`, {label: from});
  await page.selectOption(`${dialog} [name="input-To-Source"]`, {label: to});
}

async function tablesWithName(page: Page, substr: string): Promise<string[]> {
  return page.evaluate((s) => grok.shell.tables
    .filter((t: any) => t.name.toLowerCase().includes(s.toLowerCase()))
    .map((t: any) => t.name), substr);
}

test('Chem: Calculate — Map Identifiers, Biochemical Properties, Generate Conformers', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

  let baseline = 0;
  let firstPassColumns: string[] = [];

  await softStep('Setup: open smiles.csv, wait for Molecule semType + Chem menu; record baseline row count', async () => {
    await page.evaluate(async () => {
      // technical: body.selenium keeps hover-only icons visible (core/client/xamgle/web/dock-manager.css:33).
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
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

  await softStep('S1.5: Chem → Calculate → Map Identifiers, smiles → chembl, OK; a chembl column resolves > 0 rows, no blank resolved cell, row count == baseline', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Map Identifiers...', {delayMs: 700});
    await page.locator('[name="dialog-Map-Identifiers"]').waitFor({timeout: 15000});
    const molInDialog = await page.evaluate((c) => document
      .querySelector('[name="dialog-Map-Identifiers"] [name="input-host-Ids"]')
      ?.textContent?.includes(c) ?? false, MOL_COL);
    expect(molInDialog, `the Map Identifiers dialog must show the detected ${MOL_COL} Molecule column`).toBe(true);
    await setMapIdentifiersSources(page, 'smiles', 'chembl');
    await page.locator('[name="dialog-Map-Identifiers"] [name="button-OK"]').click();
    const added = await pollForColumn(page, before, 90000);
    const probe = await page.evaluate(({added, baseline}) => {
      const t = grok.shell.t;
      const name = added.find((n: string) => n.toLowerCase().includes('chembl')) ?? added[0];
      const c = name ? t.col(name) : null;
      const vals = c ? Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i)) : [];
      const resolved = vals.filter((v: any) => v != null && v !== '').length;
      return {name, hasCol: !!c, resolved, rowCount: t.rowCount, baseline};
    }, {added, baseline});
    expect(probe.hasCol, 'an identifier column must be appended after Map Identifiers').toBe(true);
    expect(probe.name, 'the appended column name must reflect the chosen toSource (chembl)').toBe('chembl');
    expect(probe.resolved,
      'the chembl column must resolve at least one identifier — an all-blank column means resolution did not run').toBeGreaterThan(0);
    expect(probe.rowCount, 'the grid row count must be unchanged after column-append').toBe(probe.baseline);
  });

  await softStep('S1.6: re-open Map Identifiers, smiles → pubchem, OK; a second column is appended and differs from chembl in ≥ 1 row; row count == baseline', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Map Identifiers...', {delayMs: 700});
    await page.locator('[name="dialog-Map-Identifiers"]').waitFor({timeout: 15000});
    await setMapIdentifiersSources(page, 'smiles', 'pubchem');
    await page.locator('[name="dialog-Map-Identifiers"] [name="button-OK"]').click();
    const added = await pollForColumn(page, before, 90000);
    const probe = await page.evaluate(({added, baseline}) => {
      const t = grok.shell.t;
      const name = added.find((n: string) => n.toLowerCase().includes('pubchem')) ?? added[0];
      const pub = name ? t.col(name) : null;
      const chembl = t.col('chembl');
      let differ = 0;
      if (pub && chembl)
        for (let i = 0; i < t.rowCount; i++) if (pub.get(i) !== chembl.get(i)) differ++;
      const resolved = pub ? Array.from({length: t.rowCount}, (_: any, i: number) => pub.get(i))
        .filter((v: any) => v != null && v !== '').length : 0;
      return {name, hasCol: !!pub, hasChembl: !!chembl, differ, resolved, rowCount: t.rowCount, baseline};
    }, {added, baseline});
    expect(probe.hasCol, 'a second identifier column (pubchem) must be appended').toBe(true);
    expect(probe.name, 'the second column name must reflect the second toSource (pubchem)').toBe('pubchem');
    expect(probe.resolved, 'the pubchem column must resolve at least one identifier').toBeGreaterThan(0);
    expect(probe.differ,
      'the two identifier columns must differ in at least one row — identical columns mean the two passes did not target distinct sources').toBeGreaterThan(0);
    expect(probe.rowCount, 'the grid row count must be unchanged after the second pass').toBe(probe.baseline);
  });

  await softStep('S2.4: Chem → Calculate → Biochemical Properties, select Chemical Properties, OK; MW is appended and every one of its rows is numeric, none echoing the SMILES, row count == baseline', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Biochemical Properties', {delayMs: 700});
    await page.locator('[name="dialog-Biochemical-Properties"]').waitFor({timeout: 15000});
    const molInDialog = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Biochemical-Properties"]');
      return dlg?.querySelector('[name="input-host-Molecules"]')?.textContent?.includes('canonical_smiles') ?? false;
    });
    expect(molInDialog, 'the Biochemical Properties dialog must show the canonical_smiles Molecule column').toBe(true);
    const navItem = page.locator('[name="dialog-Biochemical-Properties"] .biochem-calc-nav-item')
      .filter({hasText: 'Chemical Properties'});
    await navItem.locator('input[type="checkbox"]').check({force: true});
    await page.locator('[name="dialog-Biochemical-Properties"] [name="button-OK"]').click();
    const added = await pollForColumn(page, before, 90000);
    const probe = await page.evaluate(({added, baseline}) => {
      const t = grok.shell.t;
      const numericRows: Record<string, number> = {};
      const numericCols = added.filter((n: string) => {
        const c = t.col(n);
        if (!c) return false;
        let hits = 0;
        for (let i = 0; i < t.rowCount; i++) {
          const v = c.get(i);
          if (typeof v === 'number' && !Number.isNaN(v)) hits++;
        }
        numericRows[n] = hits;
        return hits === t.rowCount;
      });
      const passthroughRows: Record<string, number> = {};
      const passthrough = added.filter((n: string) => {
        const c = t.col(n);
        if (!c) return false;
        const smi = t.col('canonical_smiles');
        let same = 0;
        for (let i = 0; i < t.rowCount; i++) if (c.get(i) === smi.get(i)) same++;
        passthroughRows[n] = same;
        return same > 0;
      });
      return {addedCount: added.length, numericCols, numericRows,
        passthroughCols: passthrough, passthroughRows, rowCount: t.rowCount, baseline};
    }, {added, baseline});
    console.log(`S2.4 appended ${JSON.stringify(added)} over ${probe.rowCount} rows; numeric rows per column ` +
      `${JSON.stringify(probe.numericRows)}; rows echoing ${MOL_COL} per column ${JSON.stringify(probe.passthroughRows)}`);
    expect(probe.addedCount, 'at least one biochemical property column must be appended').toBeGreaterThan(0);
    expect(probe.numericCols,
      `every row of every appended property column must carry a numeric value — not an error string, not a raw ` +
      `stack trace; numeric rows out of ${probe.rowCount}: ${JSON.stringify(probe.numericRows)}`)
      .toEqual(added);
    expect(probe.passthroughCols,
      `no appended column may echo the original SMILES on any row; matching rows: ${JSON.stringify(probe.passthroughRows)}`)
      .toEqual([]);
    expect(probe.rowCount, 'the grid row count must be unchanged after the calculator run').toBe(probe.baseline);
    // technical: pinning the name set rules out a column produced by some other calculator.
    const bare = added.map((n) => n.replace(/ \(\d+\)$/, ''));
    expect(bare.filter((n) => !OCL_PROPERTY_NAMES.includes(n)),
      `every column Chemical Properties appends must be one of its own properties; appended: ${JSON.stringify(added)}`)
      .toEqual([]);
    expect(bare,
      `MW is the one property enabled by default, so the default run must append it; appended: ${JSON.stringify(added)}`)
      .toContain('MW');
    firstPassColumns = added;
  });

  await softStep('S2.5: re-open Biochemical Properties, select a second (logP) calculator, OK; the second pass must append clogP and nothing else, finite on every row', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Biochemical Properties', {delayMs: 700});
    await page.locator('[name="dialog-Biochemical-Properties"]').waitFor({timeout: 15000});
    const navItem = page.locator('[name="dialog-Biochemical-Properties"] .biochem-calc-nav-item')
      .filter({hasText: /^logP$/});
    await navItem.locator('input[type="checkbox"]').check({force: true});
    await page.locator('[name="dialog-Biochemical-Properties"] [name="button-OK"]').click();
    const added = await pollForColumn(page, before, 120000);
    const probe = await page.evaluate(({added, baseline, firstPass}) => {
      const t = grok.shell.t;
      const finite = added.filter((n: string) => {
        const c = t.col(n);
        if (!c) return false;
        for (let i = 0; i < t.rowCount; i++) {
          const v = c.get(i);
          if (typeof v !== 'number' || !Number.isFinite(v)) return false;
        }
        return true;
      });
      return {finiteCols: finite, firstPassResolved: firstPass.filter((n: string) => !!t.col(n)).length,
        rowCount: t.rowCount, baseline};
    }, {added, baseline, firstPass: firstPassColumns});
    expect(probe.firstPassResolved,
      `the first pass's columns must still be on the table to contrast against; expected ${JSON.stringify(firstPassColumns)}`)
      .toBe(firstPassColumns.length);
    // technical: logP is OCL CrippenClogP — exactly one column, clogP, a different kind from MW.
    expect(added.map((n) => n.replace(/ \(\d+\)$/, '')),
      `the logP pass must append its own clogP column and nothing else; appended: ${JSON.stringify(added)}`)
      .toEqual(['clogP']);
    expect(probe.finiteCols,
      `every row of the appended clogP column must hold a finite number; appended: ${JSON.stringify(added)}`)
      .toEqual(added);
    expect(probe.rowCount, 'the grid row count must be unchanged after the second calculator').toBe(probe.baseline);
  });

  await softStep('S3.4: Chem → Calculate → Generate Conformers, sketch pentane (CCCCC), Run; the conformers table must land, hold the documented structure, and be built from the sketched molecule', async () => {
    await openChemMenuItem(page, 'Generate Conformers...', {delayMs: 700});
    await page.locator('[name="dialog-Generate-Conformers"]').waitFor({timeout: 15000});
    const defaults = await page.evaluate(() => {
      const input = (n: string) => document.querySelector(
        `[name="dialog-Generate-Conformers"] [name="input-${n}"]`) as HTMLInputElement | null;
      const value = (n: string) => { const e = input(n); return e ? e.value : null; };
      return {num: value('Num-conformers'), optimize: input('Optimize')?.checked ?? null,
        rms: value('RMS-threshold'), attempts: value('Max-attempts'), seed: value('Random-seed')};
    });
    // technical: Optimize=true is what makes every energy finite below, so a changed default must fail here.
    expect(defaults, 'the Generate Conformers dialog must open on its documented defaults')
      .toEqual({num: '50', optimize: true, rms: '0.1', attempts: '5000', seed: '42'});
    await page.locator('[name="dialog-Generate-Conformers"] .d4-input-molecule-canvas').click();
    const smiles = page.locator('.d4-dialog input[placeholder*="SMILES" i]');
    await smiles.waitFor({timeout: 10000});
    // technical: the script's own default is "CCCC" (Chem/scripts/generate_conformers.py:7), so typing
    // butane would pass every assertion below even if the sketcher never received the gesture.
    await smiles.fill(MOLECULE);
    await smiles.press('Enter');
    const sketched = await smiles.inputValue();
    console.log(`S3.4 sketcher molecule input reads "${sketched}" after fill+Enter (script default is "CCCC")`);
    expect(sketched,
      'the sketched molecule must be the non-default molecule this step drives — reading back the script default means the fill did not land')
      .toBe(MOLECULE);
    const sketcher = page.locator('.d4-dialog').filter({has: page.locator('input[placeholder*="SMILES" i]')});
    await sketcher.locator('[name="button-OK"]').first().click();
    await sketcher.first().waitFor({state: 'detached', timeout: 15000});
    await page.locator('[name="dialog-Generate-Conformers"] [name="button-OK"]').click();
    await page.waitForFunction(() =>
      grok.shell.tables.some((t: any) => t.name.toLowerCase().includes('conformer')),
    null, {timeout: 180000}).catch(() => {});
    const names = await tablesWithName(page, 'conformer');
    expect(names.length,
      'Generate Conformers must produce a conformers table within 180 s — no such table means the run did not complete')
      .toBeGreaterThan(0);
    {
      const probe = await page.evaluate(({name, molecule}) => {
        const t = grok.shell.tables.find((x: any) => x.name === name);
        const get = (n: string) => {
          const c = t.col(n);
          return c ? Array.from({length: t.rowCount}, (_: any, i: number) => c.get(i)) : null;
        };
        const conformer = get('conformer');
        const molblock = get('molblock');
        const energy = get('energy');
        const rmsd = get('rmsd');
        const smi = get('smiles');
        return {
          rowCount: t.rowCount,
          smilesDistinct: smi ? Array.from(new Set(smi)) : null,
          smilesRowsMatching: smi ? smi.filter((v: any) => v === molecule).length : 0,
          hasConformer: !!conformer, firstConformer: conformer ? conformer[0] : null,
          conformerValues: conformer ? conformer.slice(0, 20) : null,
          conformerIsSequential: conformer !== null && conformer.every((v: any, i: number) => v === i + 1),
          hasMolblock: !!molblock,
          molblockNonEmptyWithEnd: molblock ? molblock.filter((m: any) =>
            typeof m === 'string' && m.length > 0 && m.includes('M  END')).length : 0,
          hasEnergy: !!energy,
          energyFinite: energy ? energy.filter((v: any) =>
            typeof v === 'number' && Number.isFinite(v)).length : 0,
          energyFirst: energy ? energy.slice(0, 5) : null,
          hasRmsd: !!rmsd, rmsdFirst: rmsd ? rmsd[0] : null,
          rmsdHasNonZero: rmsd ? rmsd.some((v: any) => typeof v === 'number' && v > 0) : false,
        };
      }, {name: names[0], molecule: MOLECULE});
      console.log(`S3.4 conformers table "${names[0]}": ${probe.rowCount} rows, smiles column holds ` +
        `${JSON.stringify(probe.smilesDistinct)} (${probe.smilesRowsMatching} rows == "${MOLECULE}"), ` +
        `first energies ${JSON.stringify(probe.energyFirst)}`);
      expect(probe.rowCount, 'the conformers table must contain at least two rows for a flexible molecule').toBeGreaterThanOrEqual(2);
      expect(probe.smilesRowsMatching,
        `every row's smiles must be the sketched molecule "${MOLECULE}" — the script defaults to "CCCC", so a table built ` +
        `from the default proves the sketcher gesture never reached the run; observed: ${JSON.stringify(probe.smilesDistinct)}`)
        .toBe(probe.rowCount);
      expect(probe.firstConformer, 'the conformer column must start its numbering at 1').toBe(1);
      expect(probe.conformerIsSequential,
        `the conformer column must increment by one per row — a repeated index means enumeration broke; first values: ${JSON.stringify(probe.conformerValues)}`)
        .toBe(true);
      expect(probe.molblockNonEmptyWithEnd,
        'every molblock cell must be a non-empty molblock ending in "M  END"').toBe(probe.rowCount);
      expect(probe.hasEnergy, 'an energy column must be present').toBe(true);
      expect(probe.energyFinite,
        `Optimize is left at its default of true, so every MMFF94 energy must be finite — an all-NaN energy column means optimization silently failed; first values: ${JSON.stringify(probe.energyFirst)}`)
        .toBe(probe.rowCount);
      expect(probe.rmsdFirst, 'the reference conformer (first row) must have RMSD 0.0').toBe(0);
      expect(probe.rmsdHasNonZero,
        'at least one non-reference conformer must have a non-zero RMSD — geometric diversity').toBe(true);
    }
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
