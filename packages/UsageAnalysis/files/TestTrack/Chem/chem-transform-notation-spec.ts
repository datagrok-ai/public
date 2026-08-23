/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const MOL_COL = 'canonical_smiles';
const SMILES_FILE = 'System:DemoFiles/chem/smiles.csv';
const NAMES_FILE = 'System:AppData/Chem/tests/names_to_smiles.csv';
// A plainly non-chemical string, so no chemical registry can ever start resolving it.
const UNRESOLVABLE_NAME = 'ZZZZ NOT A COMPOUND ZZZZ';
// Each canonical comparison costs a round-trip through Chem:canonicalize per value, so the
// row-for-row claims are checked over the first 30 valid rows rather than all 1000.
const CANONICAL_SAMPLE = 30;
// The Chem converters signal a parse failure by RETURNING a placeholder, not by throwing:
// _convertMolNotation pre-seeds its result with one of these and returns it from the finally
// block (Chem/src/utils/convert-notation-utils.ts:28-32, :60-62; the SMILES-side string is also
// Chem/src/constants.ts:30 MESSAGE_MALFORMED). Both placeholders are non-empty, and the molblock
// one even carries the "V2000" and "M  END" markers — so every shape test below needs a guard
// through a different platform function.
const MALFORMED_MOL_V2000 = '\nMalformed\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END';
const MALFORMED_INPUT_VALUE = 'MALFORMED_INPUT_VALUE';

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.columns.names());
}

/** Substitutes for picking a column / choosing a list item in the dialog's control. */
async function setDialogInputByCaption(
  page: Page, caption: string, value: {column: string} | {literal: string},
): Promise<string> {
  return page.evaluate(({caption, value}) => {
    const dlgs = DG.Dialog.getOpenDialogs();
    const d = dlgs[dlgs.length - 1];
    const input = d.inputs.find((i: any) => i.caption === caption);
    if (!input) throw new Error(`dialog input '${caption}' not found`);
    input.value = 'column' in value ? grok.shell.t.col(value.column) : value.literal;
    const rb = input.value;
    return rb?.name ?? String(rb);
  }, {caption, value} as any);
}

async function dialogBoolByCaption(page: Page, caption: string): Promise<boolean | null> {
  return page.evaluate(({caption}) => {
    const dlgs = DG.Dialog.getOpenDialogs();
    const d = dlgs[dlgs.length - 1];
    const input = d.inputs.find((i: any) => i.caption === caption);
    return input ? !!input.value : null;
  }, {caption});
}

async function waitForAddedColumn(page: Page, before: string[], capMs: number): Promise<string | null> {
  await page.waitForFunction((before) =>
    grok.shell.t.columns.names().filter((n: string) => !before.includes(n)).length >= 1,
  before, {timeout: capMs}).catch(() => {});
  const after = await columnNames(page);
  const added = after.filter((n) => !before.includes(n));
  return added.length ? added[added.length - 1] : null;
}

async function gridCellRenderer(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({column}) => {
    const gridViewer = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Grid') as any;
    try { return gridViewer?.columns?.byName(column)?.cellType ?? null; } catch (e) { return null; }
  }, {column});
}

async function canonicalMatchCount(page: Page, colA: string, colB: string, sampleN: number): Promise<{
  validRows: number; total: number; matches: number; degenerate: number;
  mismatchSamples: string[]; degenerateSamples: string[];
}> {
  return page.evaluate(async ({colA, colB, sampleN, malformed}) => {
    const t = grok.shell.t;
    const a = t.col(colA); const b = t.col(colB);
    const n = Math.min(sampleN, t.rowCount);
    // Chem:canonicalize is _convertMolNotation(mol, Unknown, Smiles), which returns
    // MALFORMED_INPUT_VALUE instead of throwing — so both sides of the equality below collapse to
    // the same sentinel on a build where RDKit parses nothing, and the row-for-row claim would
    // hold without a single molecule being compared. Chem:isSmiles reaches RDKit's own parser
    // through a different platform function, so one value cannot satisfy both.
    const isStructure = (s: any) => {
      if (typeof s !== 'string' || s === '' || s === malformed) return false;
      try { return !!grok.chem.checkSmiles(s); } catch (e) { return false; }
    };
    let matches = 0; let total = 0; let degenerate = 0; let validRows = 0;
    const mismatchSamples: string[] = []; const degenerateSamples: string[] = [];
    for (let i = 0; i < n; i++) {
      const va = a.get(i); const vb = b.get(i);
      if (va == null || va === '') continue;
      validRows++;
      // A converted value that came back empty is a lost round-trip, not a row to skip: skipping it
      // would drop it from both counters and let one surviving row satisfy the row-for-row claim.
      if (vb == null || vb === '') {
        degenerate++;
        if (degenerateSamples.length < 3)
          degenerateSamples.push(`row ${i}: source ${JSON.stringify(va)} -> converted ${JSON.stringify(vb)}`);
        continue;
      }
      const ca = await grok.functions.call('Chem:canonicalize', {molecule: va});
      const cb = await grok.functions.call('Chem:canonicalize', {molecule: vb});
      if (!isStructure(ca) || !isStructure(cb)) {
        degenerate++;
        if (degenerateSamples.length < 3)
          degenerateSamples.push(`row ${i}: ${JSON.stringify(ca)} / ${JSON.stringify(cb)}`);
        continue;
      }
      total++;
      if (ca === cb) matches++;
      else if (mismatchSamples.length < 3) mismatchSamples.push(`row ${i}: ${ca} != ${cb}`);
    }
    return {validRows, total, matches, degenerate, mismatchSamples, degenerateSamples};
  }, {colA, colB, sampleN, malformed: MALFORMED_INPUT_VALUE});
}

test('Chem: Transform notation roundtrip — Convert Notation, Recalculate Coordinates, Names To Smiles', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await softStep('Setup: open smiles.csv, wait for Molecule semType + Chem menu', async () => {
    await page.evaluate(async ({file}) => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(file);
      const detected = new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      grok.shell.addTableView(df);
      await detected;
    }, {file: SMILES_FILE});
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30000, state: 'attached'});
    await waitForChemMenu(page);
    await waitForMolecule(page);
    const molTyped = await page.evaluate((c) => grok.shell.t.col(c)?.semType === 'Molecule', MOL_COL);
    expect(molTyped, `the ${MOL_COL} column must be typed as Molecule before Convert Notation runs`).toBe(true);
  });

  let molblockCol = '';
  let roundtripCol = '';
  let coordCol = '';
  let smilesCol = '';
  let namesTotal = 0;

  await softStep('Scenario 1 Step 1-3: Chem → Transform → Convert Notation → target molblock, Overwrite off, OK', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Convert Notation...', {delayMs: 700});
    await page.locator('[name="dialog-Convert-Notation"]').waitFor({timeout: 15000});
    const molSet = await setDialogInputByCaption(page, 'Molecules', {column: MOL_COL});
    expect(molSet, 'Convert Notation must operate on the smiles column').toBe(MOL_COL);
    const tgtSet = await setDialogInputByCaption(page, 'Target Notation', {literal: 'molblock'});
    expect(tgtSet, 'Target Notation must be set to molblock').toBe('molblock');
    const overwrite = await dialogBoolByCaption(page, 'Overwrite');
    expect(overwrite, 'Overwrite must be off so the conversion adds a column instead of replacing one').toBe(false);
    await page.locator('[name="dialog-Convert-Notation"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 60000);
    expect(added, 'a new molblock column must be appended (Overwrite off)').not.toBeNull();
    molblockCol = added!;
  });

  await softStep('Scenario 1 Step 4: the new column holds non-empty V2000 molblock strings, not raw SMILES; text differs from smiles', async () => {
    const probe = await page.evaluate(({molblockCol, origCol, malformedMolblock}) => {
      const t = grok.shell.t;
      const mb = t.col(molblockCol); const orig = t.col(origCol);
      const n = t.rowCount;
      let validRows = 0; let v2000 = 0; let textDiffers = 0; let looksSmiles = 0;
      let structures = 0; let sentinels = 0; const badSamples: string[] = [];
      // MALFORMED_MOL_V2000 carries "V2000" and "M  END", holds newlines and differs from the
      // source SMILES, so it satisfies every text-shape counter above. A real molblock declares
      // at least one atom in its counts line (the placeholder declares zero) and parses through
      // Chem:isSmiles, which reaches RDKit's own parser rather than Chem:convertMolNotation.
      for (let i = 0; i < n; i++) {
        const v = mb.get(i); const o = orig.get(i);
        if (o == null || o === '') continue;
        validRows++;
        if (typeof v === 'string' && v.includes('V2000') && v.includes('M  END')) v2000++;
        if (v !== o) textDiffers++;
        if (typeof v === 'string' && v !== '' && !v.includes('\n') && !v.includes('V2000')) looksSmiles++;
        const isSentinel = typeof v === 'string' && v.trim() === malformedMolblock.trim();
        if (isSentinel) sentinels++;
        const countsLine = typeof v === 'string' ? (v.split('\n')[3] ?? '') : '';
        const atoms = parseInt(countsLine.slice(0, 3), 10);
        let rdkitParsed = false;
        try { rdkitParsed = typeof v === 'string' && v !== '' && !!grok.chem.checkSmiles(v); } catch (e) {}
        if (!isSentinel && rdkitParsed && isFinite(atoms) && atoms > 0) structures++;
        else if (badSamples.length < 3)
          badSamples.push(`row ${i}: sentinel=${isSentinel} atoms=${atoms} parsed=${rdkitParsed}`);
      }
      return {n, validRows, v2000, textDiffers, looksSmiles, structures, sentinels, badSamples, semType: mb.semType};
    }, {molblockCol, origCol: MOL_COL, malformedMolblock: MALFORMED_MOL_V2000});
    console.log(`[convert-notation] validRows=${probe.validRows} v2000=${probe.v2000} ` +
      `structures=${probe.structures} sentinels=${probe.sentinels} textDiffers=${probe.textDiffers} ` +
      `looksSmiles=${probe.looksSmiles} badSamples=${JSON.stringify(probe.badSamples)}`);
    expect(probe.validRows, 'the input must contain valid molecule rows to convert').toBeGreaterThan(0);
    expect(probe.semType, 'the molblock column carries semType Molecule').toBe('Molecule');
    expect(probe.v2000,
      'every valid row must carry a V2000 molblock string (RDKit "V2000" + "M  END" markers)').toBe(probe.validRows);
    expect(probe.sentinels,
      `no row may hold the converter's MALFORMED_MOL_V2000 placeholder — it carries the same ` +
      `markers as a real molblock; samples: ${JSON.stringify(probe.badSamples)}`).toBe(0);
    expect(probe.structures,
      `every valid row's molblock must declare at least one atom and parse through Chem:isSmiles, ` +
      `not merely carry the V2000 markers; samples: ${JSON.stringify(probe.badSamples)}`).toBe(probe.validRows);
    expect(probe.looksSmiles,
      'no valid row may hold a raw SMILES string in the molblock column').toBe(0);
    expect(probe.textDiffers,
      'the stored molblock text must visibly differ from the original SMILES text on every valid row').toBe(probe.validRows);
  });

  await softStep('Scenario 1 Step 5a: Convert Notation on the molblock column back to smiles, OK', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Convert Notation...', {delayMs: 700});
    await page.locator('[name="dialog-Convert-Notation"]').waitFor({timeout: 15000});
    const molSet = await setDialogInputByCaption(page, 'Molecules', {column: molblockCol});
    expect(molSet, 'the second Convert Notation must read the molblock column').toBe(molblockCol);
    const tgtSet = await setDialogInputByCaption(page, 'Target Notation', {literal: 'smiles'});
    expect(tgtSet, 'Target Notation must be set back to smiles').toBe('smiles');
    await page.locator('[name="dialog-Convert-Notation"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 60000);
    expect(added, 'a roundtrip smiles column must be appended').not.toBeNull();
    // 'smiles' is the registered default of Target Notation (Chem/src/package.ts:1870), so reading
    // the input back cannot tell a working setter from a dead one. The product names the output
    // `${molecules.name}_${targetNotation}` (package.ts:1881), so the added column's name is the
    // one signal that carries BOTH chosen inputs.
    expect(added!.startsWith(`${molblockCol}_smiles`),
      `the appended column must be named after the chosen inputs — expected ` +
      `'${molblockCol}_smiles', got '${added}'`).toBe(true);
    roundtripCol = added!;
  });

  await softStep('Scenario 1 Step 5b: canonical SMILES of the roundtrip column match the original row-for-row (identity preserved)', async () => {
    const cmp = await canonicalMatchCount(page, MOL_COL, roundtripCol, CANONICAL_SAMPLE);
    console.log(`[roundtrip] canonical compare validRows=${cmp.validRows} total=${cmp.total} ` +
      `matches=${cmp.matches} degenerate=${cmp.degenerate} ` +
      `degenerateSamples=${JSON.stringify(cmp.degenerateSamples)}`);
    expect(cmp.validRows, 'the sample must contain valid source molecules to round-trip').toBeGreaterThan(0);
    expect(cmp.degenerate,
      `every sampled value on both sides must canonicalize to a parseable structure, or the ` +
      `equality below compares failure sentinels; got ${JSON.stringify(cmp.degenerateSamples)}`).toBe(0);
    expect(cmp.total,
      `every valid source row must have produced a comparable round-trip value — a converter that ` +
      `returns an empty string for most of the sample would otherwise pass on the few survivors`)
      .toBe(cmp.validRows);
    expect(cmp.matches,
      `canonical SMILES must match the original row-for-row after the molblock roundtrip; ` +
      `mismatches=${JSON.stringify(cmp.mismatchSamples)}`).toBe(cmp.total);
  });

  await softStep('Scenario 2 Step 1-2: Chem → Transform → Recalculate Coordinates → Method CoordGen, Join on, OK', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Recalculate Coordinates...', {delayMs: 700});
    await page.locator('[name="dialog-Recalculate-Coordinates"]').waitFor({timeout: 15000});
    const molSet = await setDialogInputByCaption(page, 'Molecules', {column: MOL_COL});
    expect(molSet, 'Recalculate Coordinates must operate on the smiles column').toBe(MOL_COL);
    const methodSet = await setDialogInputByCaption(page, 'Method', {literal: 'CoordGen'});
    expect(methodSet, 'the coordinate method must be CoordGen').toBe('CoordGen');
    const join = await dialogBoolByCaption(page, 'Join');
    expect(join, 'Join must be on so the recalculated molblocks land beside the original column').toBe(true);
    await page.locator('[name="dialog-Recalculate-Coordinates"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 90000);
    expect(added, 'a coordinate-recalculated column must be appended (Join on)').not.toBeNull();
    coordCol = added!;
  });

  await softStep('Scenario 2 Step 3a: the recalculated molblocks carry different coordinates from the default-coordinate molblocks while the atom/bond counts stay identical', async () => {
    const probe = await page.evaluate(({coordCol, defaultCol, origCol}) => {
      const t = grok.shell.t;
      const c = t.col(coordCol); const d = t.col(defaultCol); const orig = t.col(origCol);
      // V2000 layout: title, program, comment, counts line, then one line per atom carrying x/y/z.
      const parse = (v: any) => {
        if (typeof v !== 'string' || !v.includes('V2000')) return null;
        const lines = v.split('\n');
        const counts = lines[3] || '';
        const atoms = parseInt(counts.slice(0, 3), 10);
        const bonds = parseInt(counts.slice(3, 6), 10);
        if (!isFinite(atoms) || !isFinite(bonds) || atoms <= 0) return null;
        return {atoms, bonds, coords: lines.slice(4, 4 + atoms).map((l) => l.slice(0, 31)).join('|')};
      };
      let validRows = 0; let compared = 0; let sameGraph = 0; let coordsDiffer = 0;
      const unparsedSamples: string[] = [];
      for (let i = 0; i < t.rowCount; i++) {
        const o = orig.get(i);
        if (o == null || o === '') continue;
        validRows++;
        const a = parse(c.get(i)); const b = parse(d.get(i));
        if (a == null || b == null) {
          if (unparsedSamples.length < 3)
            unparsedSamples.push(`row ${i}: recalculated=${a == null ? 'unparsed' : 'ok'}, ` +
              `default=${b == null ? 'unparsed' : 'ok'}`);
          continue;
        }
        compared++;
        if (a.atoms === b.atoms && a.bonds === b.bonds) sameGraph++;
        if (a.coords !== b.coords) coordsDiffer++;
      }
      return {validRows, compared, sameGraph, coordsDiffer, unparsedSamples};
    }, {coordCol, defaultCol: molblockCol, origCol: MOL_COL});
    console.log(`[recalc-coords] validRows=${probe.validRows} compared=${probe.compared} ` +
      `sameGraph=${probe.sameGraph} coordsDiffer=${probe.coordsDiffer} ` +
      `unparsedSamples=${JSON.stringify(probe.unparsedSamples)}`);
    expect(probe.validRows,
      'the coordinate contrast must run over a non-empty set of valid input rows').toBeGreaterThan(0);
    expect(probe.compared,
      `every valid row must yield an atom-bearing V2000 molblock in BOTH the recalculated and the ` +
      `default column — parse() rejects the zero-atom failure placeholder, and skipping such rows ` +
      `would let a column of placeholders pass on one good row; samples: ` +
      `${JSON.stringify(probe.unparsedSamples)}`).toBe(probe.validRows);
    expect(probe.sameGraph,
      'recalculating coordinates must not change the atom or bond count on any row').toBe(probe.compared);
    expect(probe.coordsDiffer,
      'the recalculated coordinates must differ from the default ones on at least one row, otherwise CoordGen did nothing').toBeGreaterThan(0);
  });

  await softStep('Scenario 2 Step 3b: Convert the recalculated column back to smiles; canonical SMILES match the original row-for-row (atom-bond graph unchanged)', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Convert Notation...', {delayMs: 700});
    await page.locator('[name="dialog-Convert-Notation"]').waitFor({timeout: 15000});
    await setDialogInputByCaption(page, 'Molecules', {column: coordCol});
    await setDialogInputByCaption(page, 'Target Notation', {literal: 'smiles'});
    await page.locator('[name="dialog-Convert-Notation"] [name="button-OK"]').click();
    const backCol = await waitForAddedColumn(page, before, 60000);
    expect(backCol, 'a smiles column derived from the recalculated column must be appended').not.toBeNull();
    const cmp = await canonicalMatchCount(page, MOL_COL, backCol!, CANONICAL_SAMPLE);
    console.log(`[recalc-roundtrip] canonical compare validRows=${cmp.validRows} total=${cmp.total} ` +
      `matches=${cmp.matches} degenerate=${cmp.degenerate} ` +
      `degenerateSamples=${JSON.stringify(cmp.degenerateSamples)}`);
    expect(cmp.validRows, 'the sample must contain valid source molecules to round-trip').toBeGreaterThan(0);
    expect(cmp.degenerate,
      `every sampled value on both sides must canonicalize to a parseable structure, or the ` +
      `equality below compares failure sentinels; got ${JSON.stringify(cmp.degenerateSamples)}`).toBe(0);
    expect(cmp.total,
      `every valid source row must have produced a comparable round-trip value — a converter that ` +
      `returns an empty string for most of the sample would otherwise pass on the few survivors`)
      .toBe(cmp.validRows);
    expect(cmp.matches,
      `coordinate recalculation must not alter the atom-bond graph: canonical SMILES must match ` +
      `the original row-for-row; mismatches=${JSON.stringify(cmp.mismatchSamples)}`).toBe(cmp.total);
  });

  await softStep('Scenario 3 Step 1: open the names_to_smiles fixture, record the blank-name rows, the known-bad name row and the total', async () => {
    await page.evaluate(async ({file}) => {
      grok.shell.closeAll();
      const ndf = await grok.dapi.files.readCsv(file);
      grok.shell.addTableView(ndf);
    }, {file: NAMES_FILE});
    await page.waitForFunction(() => grok.shell.t?.col('Name') != null, null, {timeout: 30000});
    await waitForChemMenu(page);
    const probe = await page.evaluate(({badName}) => {
      const t = grok.shell.t;
      const nameCol = t.col('Name');
      let blankNames = 0; let badNameRows = 0;
      for (let i = 0; i < t.rowCount; i++) {
        const v = nameCol.get(i);
        if (v == null || String(v).trim() === '') blankNames++;
        else if (String(v).trim() === badName) badNameRows++;
      }
      return {rowCount: t.rowCount, blankNames, badNameRows};
    }, {badName: UNRESOLVABLE_NAME});
    expect(probe.rowCount, 'the names fixture must contain a non-degenerate set of rows').toBeGreaterThan(1);
    expect(probe.blankNames,
      'the fixture must contain at least one blank name, otherwise resolved<total cannot discriminate').toBeGreaterThan(0);
    expect(probe.blankNames,
      'the fixture must contain resolvable names as well as blank ones').toBeLessThan(probe.rowCount);
    expect(probe.badNameRows,
      `the fixture must carry exactly one row named '${UNRESOLVABLE_NAME}', otherwise the ` +
      `present-but-unresolvable check is vacuous`).toBe(1);
    namesTotal = probe.rowCount;
  });

  await softStep('Scenario 3 Step 3 (via Step 2 Names To Smiles dialog): Chem → Transform → Names To Smiles → Names column, OK; resolved-row count is >0 and <total, blank names stay blank and the known-bad name yields no structure', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Names To Smiles...', {delayMs: 700});
    await page.locator('[name="dialog-Names-To-Smiles"]').waitFor({timeout: 15000});
    const namesSet = await setDialogInputByCaption(page, 'Names', {column: 'Name'});
    expect(namesSet, 'Names To Smiles must operate on the Name column').toBe('Name');
    await page.locator('[name="dialog-Names-To-Smiles"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 180000);
    expect(added, 'a SMILES column must be appended by Names To Smiles').not.toBeNull();
    smilesCol = added!;
    const probe = await page.evaluate(async ({smilesCol, badName, malformed}) => {
      const t = grok.shell.t;
      const c = t.col(smilesCol); const names = t.col('Name');
      const n = t.rowCount;
      let resolved = 0; let parseable = 0; let blankNames = 0; let blankNamesResolved = 0;
      let badNameRows = 0; let badNameResolved = 0; const badNameValues: string[] = [];
      const unparseableSamples: string[] = [];
      for (let i = 0; i < n; i++) {
        const v = c.get(i);
        const nm = names.get(i);
        const blank = nm == null || String(nm).trim() === '';
        if (blank) blankNames++;
        else if (String(nm).trim() === badName) {
          badNameRows++;
          if (v != null && v !== '') {
            badNameResolved++;
            badNameValues.push(String(v));
          }
        }
        if (v != null && v !== '') {
          resolved++;
          if (blank) blankNamesResolved++;
          // Chem:canonicalize returns MALFORMED_INPUT_VALUE instead of throwing when RDKit cannot
          // parse the value, and that return is non-empty — so a non-empty test counts any garbage
          // as parseable. Chem:isSmiles goes through RDKit's own parser instead.
          let canon: any = null;
          try { canon = await grok.functions.call('Chem:canonicalize', {molecule: v}); } catch (e) {}
          let isStructure = false;
          if (typeof canon === 'string' && canon !== '' && canon !== malformed) {
            try { isStructure = !!grok.chem.checkSmiles(canon); } catch (e) {}
          }
          if (isStructure) parseable++;
          else if (unparseableSamples.length < 3)
            unparseableSamples.push(`row ${i}: ${JSON.stringify(v)} -> ${JSON.stringify(canon)}`);
        }
      }
      return {n, resolved, parseable, blankNames, blankNamesResolved, badNameRows, badNameResolved,
        badNameValues, unparseableSamples, semType: c.semType};
    }, {smilesCol: added, badName: UNRESOLVABLE_NAME, malformed: MALFORMED_INPUT_VALUE});
    console.log(`[names-to-smiles] rows=${probe.n} resolved=${probe.resolved} ` +
      `parseable=${probe.parseable} blankNames=${probe.blankNames} badNameRows=${probe.badNameRows} ` +
      `unparseableSamples=${JSON.stringify(probe.unparseableSamples)}`);
    expect(probe.semType, 'the resolved column carries semType Molecule').toBe('Molecule');
    expect(probe.resolved,
      'at least one name must resolve to a non-empty SMILES').toBeGreaterThan(0);
    expect(probe.resolved,
      'the resolved-row count must be strictly less than the total row count').toBeLessThan(namesTotal);
    expect(probe.parseable,
      `every resolved SMILES value must canonicalize to a structure RDKit can parse (Chem:isSmiles), ` +
      `not merely to a non-empty string; samples: ${JSON.stringify(probe.unparseableSamples)}`)
      .toBe(probe.resolved);
    expect(probe.blankNamesResolved,
      'a row whose Name is blank must yield a blank SMILES, not an arbitrary molecule').toBe(0);
    expect(probe.badNameRows,
      `the row named '${UNRESOLVABLE_NAME}' must still be present after the conversion`).toBe(1);
    expect(probe.badNameResolved,
      `a present but unresolvable name must yield an empty SMILES, not a structure; ` +
      `got ${JSON.stringify(probe.badNameValues)}`).toBe(0);
  });

  await softStep('Scenario 3 Step 4: the resolved column renders molecule cells (semType Molecule + Molecule cell renderer); rows with no SMILES stay empty', async () => {
    const probe = await page.evaluate(({smilesCol}) => {
      const t = grok.shell.t;
      const c = t.col(smilesCol);
      let empty = 0; let nonEmpty = 0;
      for (let i = 0; i < t.rowCount; i++) {
        const v = c.get(i);
        if (v == null || v === '') empty++;
        else nonEmpty++;
      }
      return {semType: c.semType, empty, nonEmpty};
    }, {smilesCol});
    const cellRenderer = await gridCellRenderer(page, smilesCol);
    expect(probe.semType, 'the resolved column must carry semType Molecule so cells render as structures').toBe('Molecule');
    expect(cellRenderer, 'the grid must use the Molecule cell renderer for the resolved column').toBe('Molecule');
    expect(probe.empty,
      'rows that produced no SMILES must hold an empty string').toBeGreaterThan(0);
    expect(probe.nonEmpty,
      'resolved rows must carry non-empty molecule text').toBeGreaterThan(0);
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
