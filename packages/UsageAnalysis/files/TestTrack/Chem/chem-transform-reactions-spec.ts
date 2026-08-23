/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;

test.use(specTestOptions);

const MOL_COL = 'canonical_smiles';
const SMILES_FILE = 'System:DemoFiles/chem/smiles.csv';
// Purpose-built two-component fixture: 10 carboxylic-acid x primary-amine pairs plus two
// negative-control rows (row 11 has no acid, row 12 has no amine), so Amide Coupling yields
// 10 products over 12 rows and the "fewer products than rows" assertion stays falsifiable.
const TWO_COMPONENT_FILE = 'System:AppData/Chem/tests/amide_coupling_2_columns.csv';
const ACID_COL = 'smiles1';
const AMINE_COL = 'smiles2';
// 12-row fixture: rows 1-10 are acid/amine pairs, rows 11-12 are negative controls
// (toluene with no acid, anisole with no amine).
const PAIR_ROWS = 10;
// Purpose-built deprotection fixture: rows 1-5 are Fmoc-protected amino acids (Gly, Ala, Val,
// Leu, Phe), rows 6-10 are the same acids and benzoic acid with no protecting group, so half
// the table must change and half must pass through.
const FMOC_FILE = 'System:AppData/Chem/tests/deprotect_fmoc.csv';
const FMOC_COL = 'smiles';
const FMOC_PROTECTED_ROWS = 5;
// The Deprotect dialog's shipped default fragment, mirrored from `defaultFragment` in
// Chem/src/analysis/deprotect.ts — Fmoc (fluorenylmethyl carbamate).
const FMOC_FRAGMENT = 'O=C(N[*:1])OCC1c2ccccc2-c2ccccc21';
// RDKit matches a dummy atom only against another dummy atom, so the attachment point has to be
// dropped before the fragment can be used as a substructure query (with [*:1] left in place the
// query matches nothing). Being derived from the mirrored constant, the query cannot see a change
// to the product's SHIPPED default — the match count would stay at FMOC_PROTECTED_ROWS and the
// matching-row guard would stay silent. What fails on a changed shipped default is the changed-rows
// assertion in Scenario 3 Step 3-4: a dialog defaulting to some other fragment leaves the Fmoc rows
// untouched, so the changed set no longer equals the matched set. This guard's own job is narrower —
// it catches a fixture that drifted away from the mirrored constant.
const FMOC_QUERY = FMOC_FRAGMENT.replace('[*:1]', '');
const ESTER_SMARTS = '[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3]';
// Water, hydrogen halides, halides and sulfuric acid written as standalone canonical fragments —
// what "Remove Water and Salts" must strip. Kept short and explicit so the assertion states an
// expectation rather than re-deriving the feature's own fragment table.
const SALT_FRAGMENTS = ['O', 'Cl', 'Br', 'I', 'F', 'O=S(=O)(O)O'];
// _convertMolNotation swallows a parse failure and RETURNS a sentinel instead of throwing
// (Chem/src/utils/convert-notation-utils.ts#L62 picks 'MALFORMED_INPUT_VALUE' for a SMILES target,
// #L28-36 the blank MALFORMED_MOL_V2000/V3000 molblocks whose title line reads `Malformed`). Both
// are non-empty and differ from every reactant, so with RDKit degraded a whole column of sentinels
// satisfies every "non-empty and different" comparison in this spec. Every place a converted or
// produced molecule is compared or accepted must reject these two forms first.
const MALFORMED_SMILES_SENTINEL = 'MALFORMED_INPUT_VALUE';
const MALFORMED_MOLBLOCK_TITLE = 'Malformed';
const SENTINELS = {smiles: MALFORMED_SMILES_SENTINEL, molblockTitle: MALFORMED_MOLBLOCK_TITLE};
// Two negative-control rows accompany the acid/amine pairs (row 11 has no acid, row 12 no amine).
const CONTROL_ROWS = 2;
const FIXTURE_ROWS = PAIR_ROWS + CONTROL_ROWS;
// Reactant halves of the shipped Amide Coupling SMARTS
// [C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3], atom maps dropped so each can be used as a
// substructure query. Which rows must react is then READ off the fixture — a row reacts iff its
// acid cell carries a carboxylic acid and its amine cell a primary amine — instead of being
// assumed from a count.
const ACID_QUERY = 'C(=O)[OH]';
// Written with the `;` low-precedence AND rather than as `[NH2]`: _isSmarts
// (Chem/src/utils/mol-creation_rdkit.ts#L21-26) only calls a string SMARTS when it carries one of
// `[.?#\d`, `$`, `&`, `;`, `,`, `!`, and `[NH2]` carries none of them — getQueryMolSafe then takes
// the plain-SMILES branch and the query becomes a molecule that matches no amine at all.
const AMINE_QUERY = '[NX3;H2]';
// Single-molecule positive controls for the two substructure channels below: a carboxylic acid the
// acid query must find and only it, a primary amine the amine query must find and only it.
const ACID_PROBE = 'CC(=O)O';
const AMINE_PROBE = 'NCC';
// Matches neither probe — pads the control column so no cell is empty.
const INERT_PROBE = 'Cc1ccccc1';

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.columns.names());
}

async function waitForAddedColumn(page: Page, before: string[], capMs: number): Promise<string | null> {
  await page.waitForFunction((before) =>
    grok.shell.t.columns.names().filter((n: string) => !before.includes(n)).length >= 1,
  before, {timeout: capMs}).catch(() => {});
  const after = await columnNames(page);
  const added = after.filter((n) => !before.includes(n));
  // Resolves the LAST added name: if an operation ever appends more than one column, the one under
  // test may not be the one returned, so the full set is logged before the choice is made.
  console.log(`[transform-reactions] columns added: ${JSON.stringify(added)}`);
  return added.length ? added[added.length - 1] : null;
}

async function openTable(page: Page, file: string): Promise<void> {
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
  }, {file});
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30000, state: 'attached'});
  await waitForChemMenu(page);
  await waitForMolecule(page);
}

async function selectReactionCard(page: Page, dialogSel: string, name: string): Promise<string | null> {
  return page.evaluate(({dialogSel, name}) => {
    const dlg = document.querySelector(dialogSel);
    if (!dlg) return null;
    const cards = Array.from(dlg.querySelectorAll('div.d4-flex-col.ui-div')).filter((d: any) =>
      d.firstElementChild && d.firstElementChild.tagName === 'CANVAS' && d.querySelector('label.ui-label'));
    const target = cards.find((c: any) => c.querySelector('label.ui-label').textContent.trim() === name);
    if (!target) return null;
    (target as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true}));
    const sel = dlg.querySelector('.reaction-card-selected');
    return sel ? sel.querySelector('label.ui-label')!.textContent!.trim() : null;
  }, {dialogSel, name});
}

async function reactionProductCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const texts = Array.from(document.querySelectorAll('.d4-balloon')).map((b) => b.textContent || '');
    const withProducts = texts.reverse().find((t) => /(\d+)\s+products/.test(t));
    if (!withProducts) return -1;
    const m = withProducts.match(/(\d+)\s+products/);
    return m ? parseInt(m[1], 10) : -1;
  });
}

async function gridCellRenderer(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({column}) => {
    const gridViewer = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Grid') as any;
    try { return gridViewer?.columns?.byName(column)?.cellType ?? null; } catch (e) { return null; }
  }, {column});
}

test('Chem: Transform Reactions — Remove Water and Salts, Transformation, Deprotect, Two-Component Reaction', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  let desaltCol = '';
  let deprotectCol = '';
  let reactedCol = '';
  let productCol = '';
  let rowCount = 0;
  let reactionProducts = -1;
  let reactionChangedCells = -1;
  let fmocMatchedRows: number[] = [];

  await softStep('Setup: open smiles.csv, wait for Molecule semType + Chem menu', async () => {
    await openTable(page, SMILES_FILE);
    const molTyped = await page.evaluate((c) => grok.shell.t.col(c)?.semType === 'Molecule', MOL_COL);
    expect(molTyped, `the ${MOL_COL} column must be typed as Molecule before the reactions run`).toBe(true);
  });

  await softStep('Scenario 1 Step 1-3: Chem → Transform → Reactions → Remove Water and Salts → OK', async () => {
    const before = await columnNames(page);
    await openChemMenuItem(page, 'Remove Water and Salts...', {delayMs: 700});
    await page.locator('[name="dialog-Remove-Water-and-Salts"]').waitFor({timeout: 15000});
    const molSet = await page.evaluate(() =>
      document.querySelector('[name="dialog-Remove-Water-and-Salts"] [name="input-host-Molecules"] .d4-column-selector-column')?.textContent);
    expect(molSet, 'Remove Water and Salts must operate on the smiles column').toBe(MOL_COL);
    await page.locator('[name="dialog-Remove-Water-and-Salts"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 90000);
    expect(added, 'a Desalted column must be appended').not.toBeNull();
    desaltCol = added!;
  });

  await softStep('Scenario 1 Step 4: salt and water fragments are stripped, no row gains a fragment, single-fragment inputs keep their molecule', async () => {
    const probe = await page.evaluate(({desaltCol, origCol, saltFragments}) => {
      const t = grok.shell.t;
      const src = t.col(origCol); const out = t.col(desaltCol);
      const frags = (v: any) => String(v).split('.');
      let validRows = 0; let emptyOnValid = 0; let multiFragInput = 0;
      let saltBearingInput = 0; let saltBearingOutput = 0; let gainedFragments = 0; let changedRows = 0;
      const singleFragRows: number[] = [];
      for (let i = 0; i < t.rowCount; i++) {
        const inp = src.get(i); const o = out.get(i);
        if (inp == null || inp === '') continue;
        validRows++;
        if (o == null || String(o).trim() === '') { emptyOnValid++; continue; }
        const inFrags = frags(inp); const outFrags = frags(o);
        if (inFrags.length > 1) multiFragInput++;
        else if (singleFragRows.length < 30) singleFragRows.push(i);
        if (inFrags.some((f) => saltFragments.includes(f))) saltBearingInput++;
        if (outFrags.some((f) => saltFragments.includes(f))) saltBearingOutput++;
        if (outFrags.length > inFrags.length) gainedFragments++;
        if (String(o) !== String(inp)) changedRows++;
      }
      return {validRows, emptyOnValid, multiFragInput, saltBearingInput, saltBearingOutput,
        gainedFragments, changedRows, singleFragRows, semType: out.semType};
    }, {desaltCol, origCol: MOL_COL, saltFragments: SALT_FRAGMENTS});
    expect(probe.semType, 'the Desalted column carries semType Molecule').toBe('Molecule');
    expect(probe.validRows, 'the input must contain valid molecule rows').toBeGreaterThan(0);
    expect(probe.multiFragInput,
      'smiles.csv must contain multi-fragment rows for the desalt assertion to be falsifiable').toBeGreaterThan(0);
    expect(probe.saltBearingInput,
      'smiles.csv must contain rows carrying a water/halide/sulfate fragment, otherwise nothing is there to strip').toBeGreaterThan(0);
    expect(probe.saltBearingOutput,
      'no Desalted value may still carry a water/halide/sulfate fragment').toBe(0);
    expect(probe.gainedFragments,
      'no Desalted value may hold more disconnected fragments than its input').toBe(0);
    expect(probe.changedRows,
      'the Desalted column must not be a row-for-row copy of its input').toBeGreaterThan(0);
    expect(probe.emptyOnValid,
      'no valid-input row may produce an empty Desalted cell').toBe(0);

    // Both ends of this comparison come out of the same RDKit conversion path, so a degraded RDKit
    // collapses them onto the same failure sentinel and the equality holds no matter what Remove
    // Water and Salts wrote. Reject the sentinel on each end before comparing.
    const passthrough = await page.evaluate(async ({desaltCol, origCol, rows, sentinels}) => {
      const t = grok.shell.t;
      const src = t.col(origCol); const out = t.col(desaltCol);
      const isSentinel = (v: any) => {
        const s = String(v ?? '');
        return s === sentinels.smiles || s.split('\n').some((l: string) => l.trim() === sentinels.molblockTitle);
      };
      let compared = 0; let sameMolecule = 0; let sentinelEnds = 0; let emptyEnds = 0;
      const sentinelSamples: string[] = [];
      for (const i of rows) {
        const canon = await grok.functions.call('Chem:canonicalize', {molecule: src.get(i)});
        const o = out.get(i);
        if (canon == null || String(canon) === '' || o == null || String(o).trim() === '') { emptyEnds++; continue; }
        if (isSentinel(canon) || isSentinel(o)) {
          sentinelEnds++;
          if (sentinelSamples.length < 3)
            sentinelSamples.push(`row ${i}: canonical input=${JSON.stringify(canon)}, output=${JSON.stringify(o)}`);
          continue;
        }
        compared++;
        if (String(canon) === String(o)) sameMolecule++;
      }
      return {compared, sameMolecule, sentinelEnds, emptyEnds, sentinelSamples};
    }, {desaltCol, origCol: MOL_COL, rows: probe.singleFragRows, sentinels: SENTINELS});
    console.log(`[transform-reactions] passthrough: singleFragRows=${probe.singleFragRows.length} ` +
      `compared=${passthrough.compared} sameMolecule=${passthrough.sameMolecule} ` +
      `sentinelEnds=${passthrough.sentinelEnds} emptyEnds=${passthrough.emptyEnds} ` +
      `sentinelSamples=${JSON.stringify(passthrough.sentinelSamples)}`);
    expect(passthrough.sentinelEnds,
      'neither the canonicalized input nor the Desalted output may be an RDKit failure sentinel — ' +
      `two sentinels compare equal, so the passthrough claim below would hold on a broken RDKit; ` +
      `samples=${JSON.stringify(passthrough.sentinelSamples)}`).toBe(0);
    expect(passthrough.compared,
      'the passthrough check must cover a non-degenerate set of single-fragment rows').toBeGreaterThan(0);
    expect(passthrough.sameMolecule,
      'a single-fragment input must pass through as the same molecule (canonical form unchanged)').toBe(passthrough.compared);
  });

  await softStep('Scenario 1 Step 5: the Desalted column renders molecule cells (semType Molecule + Molecule cell renderer)', async () => {
    const semType = await page.evaluate((c) => grok.shell.t.col(c).semType, desaltCol);
    const cellRenderer = await gridCellRenderer(page, desaltCol);
    expect(semType, 'the Desalted column must carry semType Molecule so cells render as structures').toBe('Molecule');
    expect(cellRenderer, 'the grid must use the Molecule cell renderer for the Desalted column').toBe('Molecule');
  });

  await softStep('Scenario 2 Step 1-2: Chem → Transform → Reactions → Transformation; create the ester-hydrolysis reaction via New Reaction', async () => {
    await openChemMenuItem(page, 'Transformation...', {delayMs: 700});
    await page.locator('[name="dialog-Run-Reaction"]').waitFor({timeout: 15000});
    const molSet = await page.evaluate(() =>
      document.querySelector('[name="dialog-Run-Reaction"] [name="input-host-Molecules"] .d4-column-selector-column')?.textContent);
    expect(molSet, 'Transformation must operate on the smiles column').toBe(MOL_COL);
    await page.locator('[name="dialog-Run-Reaction"] [name="button-+-New-Reaction"]').click();
    await page.locator('[name="dialog-New-Reaction"]').waitFor({timeout: 15000});
    // Substitutes for typing the reaction name and SMARTS into the two text inputs.
    await page.evaluate((smarts) => {
      const ed = document.querySelector('[name="dialog-New-Reaction"]')!;
      const set = (sel: string, v: string) => {
        const e = ed.querySelector(sel) as HTMLInputElement;
        e.value = v;
        e.dispatchEvent(new Event('input', {bubbles: true}));
        e.dispatchEvent(new Event('change', {bubbles: true}));
      };
      set('[name="input-Reaction-Name"]', 'EsterHydrolysisTest');
      set('[name="input-Reaction-SMARTS"]', smarts);
    }, ESTER_SMARTS);
    await page.locator('[name="dialog-New-Reaction"] [name="button-OK"]').click();
    await page.locator('[name="dialog-New-Reaction"]').waitFor({state: 'detached', timeout: 15000});
    const selected = await selectReactionCard(page, '[name="dialog-Run-Reaction"]', 'EsterHydrolysisTest');
    expect(selected,
      'the created ester-hydrolysis reaction card must be selected (reaction-card-selected armed)').toBe('EsterHydrolysisTest');
  });

  await softStep('Scenario 2 Step 4: run the transformation; the product column holds a product on every valid row and a changed product on some of them', async () => {
    const before = await columnNames(page);
    rowCount = await page.evaluate(() => grok.shell.t.rowCount);
    await page.locator('[name="dialog-Run-Reaction"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 120000);
    expect(added, 'a product column must be appended by the transformation').not.toBeNull();
    reactedCol = added!;
    await page.locator('.d4-balloon:has-text("products")').first().waitFor({timeout: 30000}).catch(() => {});
    reactionProducts = await reactionProductCount(page);
    // The transformation desalts and canonicalizes each reactant before running the SMARTS, so the
    // Desalted column from Scenario 1 — not the raw input — is what an unreacted row is copied from.
    const probe = await page.evaluate(({reactedCol, desaltCol, origCol}) => {
      const t = grok.shell.t;
      const src = t.col(origCol); const base = t.col(desaltCol); const p = t.col(reactedCol);
      let validRows = 0; let emptyOnValid = 0; let changed = 0;
      const changedSamples: string[] = [];
      for (let i = 0; i < t.rowCount; i++) {
        const inp = src.get(i);
        if (inp == null || inp === '') continue;
        validRows++;
        const v = p.get(i);
        if (v == null || String(v).trim() === '') { emptyOnValid++; continue; }
        if (String(v) !== String(base.get(i))) {
          changed++;
          if (changedSamples.length < 3) changedSamples.push(`row ${i}: ${base.get(i)} -> ${v}`);
        }
      }
      return {validRows, emptyOnValid, changed, changedSamples, semType: p.semType};
    }, {reactedCol: added, desaltCol, origCol: MOL_COL});
    reactionChangedCells = probe.changed;
    expect(probe.semType, 'the transformation product column carries semType Molecule').toBe('Molecule');
    expect(probe.validRows, 'the input must contain valid molecule rows to transform').toBeGreaterThan(0);
    expect(probe.emptyOnValid,
      'a row the SMARTS did not match must carry the desalted input, not an empty cell').toBe(0);
    expect(probe.changed,
      `at least one product cell must differ from its desalted input (the ester was hydrolysed); samples=${JSON.stringify(probe.changedSamples)}`).toBeGreaterThan(0);
    expect(probe.changed,
      'not every row may carry a changed product — smiles.csv is not all esters').toBeLessThan(probe.validRows);
  });

  await softStep('Scenario 2 Step 5: the reported product count is between 1 and the row count, and equals the number of changed product cells', async () => {
    expect(reactionProducts,
      'the balloon must report at least one product (some rows matched the SMARTS)').toBeGreaterThan(0);
    expect(reactionProducts,
      'the balloon must report fewer products than rows (some rows did not match the SMARTS)').toBeLessThan(rowCount);
    expect(reactionChangedCells,
      'the reported product count must equal the number of product cells that differ from their desalted input').toBe(reactionProducts);
  });

  await softStep('Scenario 3 Step 1-2: open the Fmoc fixture, count the rows the default fragment matches, then Chem → Transform → Reactions → Deprotect (default Fmoc fragment) → OK', async () => {
    await openTable(page, FMOC_FILE);
    const molTyped = await page.evaluate((c) => grok.shell.t.col(c)?.semType === 'Molecule', FMOC_COL);
    expect(molTyped, `the ${FMOC_COL} column must be typed as Molecule before Deprotect runs`).toBe(true);
    fmocMatchedRows = await page.evaluate(async ({col, query}) => {
      const t = grok.shell.t;
      const bs = await grok.chem.searchSubstructure(t.col(col), query);
      return Array.from({length: t.rowCount}, (_, i) => i).filter((i) => bs.get(i));
    }, {col: FMOC_COL, query: FMOC_QUERY});
    expect(fmocMatchedRows.length,
      `the fixture promises ${FMOC_PROTECTED_ROWS} rows carrying the dialog's default protecting group; ` +
      'a different count means the fixture no longer matches the shipped fragment and the scenario would ' +
      'otherwise pass on pass-through alone').toBe(FMOC_PROTECTED_ROWS);

    const before = await columnNames(page);
    await openChemMenuItem(page, 'Deprotect...', {delayMs: 700});
    await page.locator('[name="dialog-Deprotect"]').waitFor({timeout: 15000});
    const molSet = await page.evaluate(() =>
      document.querySelector('[name="dialog-Deprotect"] [name="input-host-Molecules"] .d4-column-selector-column')?.textContent);
    expect(molSet, 'Deprotect must operate on the fixture molecule column').toBe(FMOC_COL);
    await page.locator('[name="dialog-Deprotect"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 90000);
    expect(added, 'a Deprotected column must be appended').not.toBeNull();
    deprotectCol = added!;
  });

  await softStep('Scenario 3 Step 3-4: exactly the protected rows change and the clean rows pass through unchanged', async () => {
    // A matched row comes back as a molblock and an unmatched one as the input SMILES, so both
    // sides are canonicalized before they are compared — a raw string compare against a SMILES
    // fragment can never see a molblock product.
    const probe = await page.evaluate(async ({deprotectCol, origCol, query}) => {
      const t = grok.shell.t;
      const src = t.col(origCol); const out = t.col(deprotectCol);
      const canon = async (v: string) => grok.functions.call('Chem:canonicalize', {molecule: v});
      const changed: number[] = []; const unchanged: number[] = [];
      let validRows = 0; let emptyOnValid = 0; let uncanonicalizable = 0;
      for (let i = 0; i < t.rowCount; i++) {
        const inp = src.get(i);
        if (inp == null || inp === '') continue;
        validRows++;
        const o = out.get(i);
        if (o == null || String(o).trim() === '') { emptyOnValid++; continue; }
        const ci = await canon(String(inp)); const co = await canon(String(o));
        if (!ci || !co) { uncanonicalizable++; continue; }
        (ci === co ? unchanged : changed).push(i);
      }
      // The fragment check runs over the canonical SMILES of the output rather than over the raw
      // cells. A search over a hand-built canonical column is its own channel, so it gets its own
      // positive control: the same construction applied to the INPUT column must still find the
      // fragment on the protected rows. Without it, a channel that silently matches nothing would
      // make the empty result below unfalsifiable.
      const buildCanonCol = async (name: string, source: any) => {
        const c = t.columns.addNewString(t.columns.getUnusedName(name));
        for (let i = 0; i < t.rowCount; i++)
          c.set(i, await canon(String(source.get(i) ?? '')));
        c.semType = 'Molecule';
        return c;
      };
      const dropCol = (c: any) => {
        for (const n of t.columns.names().filter((n: string) => n === c.name || n.startsWith(`~${c.name}.`)))
          t.columns.remove(n);
      };
      const rowIdx = Array.from({length: t.rowCount}, (_, i) => i);
      const controlCol = await buildCanonCol('input_canonical', src);
      const controlBs = await grok.chem.searchSubstructure(controlCol, query);
      const fragmentInControl = rowIdx.filter((i) => controlBs.get(i));
      dropCol(controlCol);
      const canonCol = await buildCanonCol('deprotected_canonical', out);
      const bs = await grok.chem.searchSubstructure(canonCol, query);
      const fragmentInOutput = rowIdx.filter((i) => bs.get(i));
      dropCol(canonCol);
      return {validRows, emptyOnValid, uncanonicalizable, changed, unchanged,
        fragmentInControl, fragmentInOutput, semType: out.semType};
    }, {deprotectCol, origCol: FMOC_COL, query: FMOC_QUERY});
    expect(probe.semType, 'the Deprotected column carries semType Molecule').toBe('Molecule');
    expect(probe.validRows, 'the fixture must contain valid molecule rows to deprotect').toBe(10);
    expect(probe.emptyOnValid,
      'the Deprotect operation must run to completion: no valid-input row may yield an empty Deprotected cell').toBe(0);
    expect(probe.uncanonicalizable,
      'every input and output cell must canonicalize, otherwise the change comparison is meaningless').toBe(0);
    expect(probe.changed,
      'exactly the rows carrying the protecting group must change — a pass-through leaves this set empty')
      .toEqual(fmocMatchedRows);
    expect(probe.unchanged.length,
      'every row without the protecting group must pass through as the same molecule')
      .toBe(probe.validRows - FMOC_PROTECTED_ROWS);
    console.log(`[transform-reactions] deprotect fragment channel: control rows=` +
      `${JSON.stringify(probe.fragmentInControl)} expected=${JSON.stringify(fmocMatchedRows)} ` +
      `output rows=${JSON.stringify(probe.fragmentInOutput)}`);
    expect(probe.fragmentInControl,
      'positive control on the same channel: searching a hand-built canonical column made from the ' +
      'INPUT must still find the protecting group on the protected rows — otherwise an empty result ' +
      'below proves nothing about the output').toEqual(fmocMatchedRows);
    expect(probe.fragmentInOutput,
      'no deprotected output may still contain the protecting group').toEqual([]);
  });

  await softStep('Scenario 4 Step 1-3: open the acid/amine fixture, then Chem → Transform → Reactions → Two-Component Reaction → Amide Coupling, pairwise', async () => {
    await openTable(page, TWO_COMPONENT_FILE);
    await page.waitForFunction(({acid, amine}) =>
      [acid, amine].every((c) => grok.shell.t.col(c)?.semType === 'Molecule'),
    {acid: ACID_COL, amine: AMINE_COL}, {timeout: 45000}).catch(() => {});
    const semTypes = await page.evaluate(({acid, amine}) =>
      [grok.shell.t.col(acid)?.semType, grok.shell.t.col(amine)?.semType], {acid: ACID_COL, amine: AMINE_COL});
    expect(semTypes, 'both fixture columns must be typed as Molecule so the two reactant selectors can bind them')
      .toEqual(['Molecule', 'Molecule']);
    await openChemMenuItem(page, 'Two-Component Reaction...', {delayMs: 700});
    await page.locator('[name="dialog-Two-Component-Reaction"]').waitFor({timeout: 15000});
    const cols = await page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="dialog-Two-Component-Reaction"] .d4-column-selector-column')).map((c) => c.textContent));
    expect(cols[0], 'Reactant 1 must be the carboxylic-acid column').toBe(ACID_COL);
    expect(cols[1], 'Reactant 2 must be the primary-amine column').toBe(AMINE_COL);
    const mode = await page.evaluate(() =>
      (document.querySelector('[name="dialog-Two-Component-Reaction"] [name="input-Combination-Mode"]') as HTMLSelectElement).value);
    expect(mode, 'Combination Mode must be pairwise — each row is one acid/amine pair').toBe('pairwise');
    const selected = await selectReactionCard(page, '[name="dialog-Two-Component-Reaction"]', 'Amide Coupling');
    expect(selected, 'the Amide Coupling reaction card must be selected').toBe('Amide Coupling');
  });

  await softStep('Scenario 4 Step 4: Run pairwise; a product column is added whose non-empty rows are exactly the rows carrying both reactants, each product parsing as a molecule and differing from both the acid and the amine reactant', async () => {
    const before = await columnNames(page);
    // Which rows must react is read off the fixture before the reaction runs, by searching each
    // reactant column for the half of the Amide Coupling SMARTS it has to supply. A count alone is
    // satisfied by any ten-of-twelve distribution — including one where both negative controls react
    // and two genuine pairs do not — and by a twelve-row table that lost its controls entirely.
    const fixture = await page.evaluate(async ({colA, colB, acidQuery, amineQuery, probes}) => {
      const t = grok.shell.t;
      const rowIdx = Array.from({length: t.rowCount}, (_, i) => i);
      const acidBs = await grok.chem.searchSubstructure(t.col(colA), acidQuery);
      const amineBs = await grok.chem.searchSubstructure(t.col(colB), amineQuery);
      const acidRows = rowIdx.filter((i) => acidBs.get(i));
      const amineRows = rowIdx.filter((i) => amineBs.get(i));
      // Each half of the derivation gets a positive control on its own channel before the pair-row
      // set below is trusted: a query that silently stops being a query returns nothing, and an
      // empty result is otherwise indistinguishable from a fixture that lost its rows.
      const probeCol = t.columns.addNewString(t.columns.getUnusedName('reactant_query_control'));
      for (const i of rowIdx)
        probeCol.set(i, i === 0 ? probes.acid : i === 1 ? probes.amine : probes.inert);
      probeCol.semType = 'Molecule';
      const acidProbeBs = await grok.chem.searchSubstructure(probeCol, acidQuery);
      const amineProbeBs = await grok.chem.searchSubstructure(probeCol, amineQuery);
      const acidProbeRows = rowIdx.filter((i) => acidProbeBs.get(i));
      const amineProbeRows = rowIdx.filter((i) => amineProbeBs.get(i));
      for (const n of t.columns.names().filter((n: string) => n === probeCol.name || n.startsWith(`~${probeCol.name}.`)))
        t.columns.remove(n);
      return {
        rowCount: t.rowCount, acidRows, amineRows, acidProbeRows, amineProbeRows,
        pairRows: rowIdx.filter((i) => acidBs.get(i) && amineBs.get(i)),
        controlRows: rowIdx.filter((i) => !(acidBs.get(i) && amineBs.get(i))),
      };
    }, {colA: ACID_COL, colB: AMINE_COL, acidQuery: ACID_QUERY, amineQuery: AMINE_QUERY,
      probes: {acid: ACID_PROBE, amine: AMINE_PROBE, inert: INERT_PROBE}});
    console.log(`[transform-reactions] fixture shape: rows=${fixture.rowCount} ` +
      `acidRows=${JSON.stringify(fixture.acidRows)} amineRows=${JSON.stringify(fixture.amineRows)} ` +
      `pairRows=${JSON.stringify(fixture.pairRows)} controlRows=${JSON.stringify(fixture.controlRows)}`);
    console.log(`[transform-reactions] reactant query controls: acidQuery=${ACID_QUERY} ` +
      `matched=${JSON.stringify(fixture.acidProbeRows)} amineQuery=${AMINE_QUERY} ` +
      `matched=${JSON.stringify(fixture.amineProbeRows)} (row 0=${ACID_PROBE}, row 1=${AMINE_PROBE}, ` +
      `rest=${INERT_PROBE})`);
    expect(fixture.acidProbeRows,
      `positive control on the acid channel: ${ACID_QUERY} must find the known carboxylic acid ` +
      `${ACID_PROBE} at row 0 and nothing else — a query that matches nothing makes every ` +
      'acid-derived row set below unfalsifiable').toEqual([0]);
    expect(fixture.amineProbeRows,
      `positive control on the amine channel: ${AMINE_QUERY} must find the known primary amine ` +
      `${AMINE_PROBE} at row 1 and nothing else — a query that matches nothing makes every ` +
      'amine-derived row set below unfalsifiable').toEqual([1]);
    expect(fixture.rowCount,
      `the fixture promises ${PAIR_ROWS} acid/amine pairs plus ${CONTROL_ROWS} negative-control rows`)
      .toBe(FIXTURE_ROWS);
    expect(fixture.pairRows.length,
      `exactly ${PAIR_ROWS} rows may carry both a carboxylic acid and a primary amine; a different ` +
      'count means the fixture no longer holds the pairs and controls the scenario names, and the ' +
      'product-row assertion below would be asserting against a shape nobody checked').toBe(PAIR_ROWS);
    // The two control rows are asserted through the reactant columns, not through their count:
    // controlRows is the complement of pairRows, so with the row count and the pair count already
    // pinned its length is arithmetic. A fixture whose BOTH controls lacked an amine would satisfy
    // every count above while the scenario's "one missing each reactant" expectation lapsed.
    expect(fixture.acidRows.length,
      `exactly one row may lack a carboxylic acid: ${ACID_QUERY} must match ${FIXTURE_ROWS - 1} of the ` +
      `${FIXTURE_ROWS} fixture rows`).toBe(FIXTURE_ROWS - 1);
    expect(fixture.amineRows.length,
      `exactly one row may lack a primary amine: ${AMINE_QUERY} must match ${FIXTURE_ROWS - 1} of the ` +
      `${FIXTURE_ROWS} fixture rows — together with the acid count above this pins one control row ` +
      'per missing reactant class').toBe(FIXTURE_ROWS - 1);
    await page.locator('[name="dialog-Two-Component-Reaction"] [name="button-OK"]').click();
    const added = await waitForAddedColumn(page, before, 120000);
    expect(added, 'a two-component product column must be appended').not.toBeNull();
    productCol = added!;
    // Nothing below reads a product AS A MOLECULE unless it is canonicalized: the failure sentinels
    // are non-empty and differ from every reactant, so a column full of them satisfies the
    // non-empty and differs-from-both claims on their own.
    const probe = await page.evaluate(async ({productCol, colA, colB, sentinels}) => {
      const t = grok.shell.t;
      const p = t.col(productCol); const a = t.col(colA); const b = t.col(colB);
      const isSentinel = (v: any) => {
        const s = String(v ?? '');
        return s === sentinels.smiles || s.split('\n').some((l: string) => l.trim() === sentinels.molblockTitle);
      };
      const productRows: number[] = []; const sentinelSamples: string[] = [];
      let differsBoth = 0; let sentinelProducts = 0; let unparsable = 0;
      for (let i = 0; i < t.rowCount; i++) {
        const v = p.get(i);
        if (v == null || String(v).trim() === '') continue;
        productRows.push(i);
        if (isSentinel(v)) {
          sentinelProducts++;
          if (sentinelSamples.length < 3) sentinelSamples.push(`row ${i}: ${JSON.stringify(v)}`);
        }
        else {
          const c = await grok.functions.call('Chem:canonicalize', {molecule: String(v)});
          if (c == null || String(c) === '' || isSentinel(c)) {
            unparsable++;
            if (sentinelSamples.length < 3)
              sentinelSamples.push(`row ${i}: ${JSON.stringify(v)} -> ${JSON.stringify(c)}`);
          }
        }
        if (v !== a.get(i) && v !== b.get(i)) differsBoth++;
      }
      return {productRows, differsBoth, sentinelProducts, unparsable, sentinelSamples, semType: p.semType};
    }, {productCol: added, colA: ACID_COL, colB: AMINE_COL, sentinels: SENTINELS});
    console.log(`[transform-reactions] products: rows=${JSON.stringify(probe.productRows)} ` +
      `expected=${JSON.stringify(fixture.pairRows)} differsBoth=${probe.differsBoth} ` +
      `sentinelProducts=${probe.sentinelProducts} unparsable=${probe.unparsable} ` +
      `samples=${JSON.stringify(probe.sentinelSamples)}`);
    expect(probe.semType, 'the two-component product column carries semType Molecule').toBe('Molecule');
    expect(probe.sentinelProducts,
      'no product cell may hold an RDKit failure sentinel — semType is stamped on the column ' +
      `regardless of content, so a sentinel would pass every other check in this step; ` +
      `samples=${JSON.stringify(probe.sentinelSamples)}`).toBe(0);
    expect(probe.unparsable,
      'every non-empty product must canonicalize back to a molecule, otherwise the cell holds text ' +
      `that only looks like a product; samples=${JSON.stringify(probe.sentinelSamples)}`).toBe(0);
    // Identity, not count: the rows that produced a product must be exactly the rows that carry both
    // reactants. Ten of twelve is also satisfied by both controls reacting and two pairs failing.
    expect(probe.productRows,
      `exactly the ${PAIR_ROWS} rows carrying both a carboxylic acid and a primary amine must produce ` +
      `an amide product, and the ${CONTROL_ROWS} negative-control rows none`).toEqual(fixture.pairRows);
    expect(probe.differsBoth,
      'every non-empty product must differ from BOTH input reactants (a genuine coupling of two distinct molecules, not a passthrough)').toBe(probe.productRows.length);
  });

  await softStep('Scenario 4 Step 5: the product column renders molecule cells (semType Molecule + Molecule cell renderer)', async () => {
    const semType = await page.evaluate((c) => grok.shell.t.col(c).semType, productCol);
    const cellRenderer = await gridCellRenderer(page, productCol);
    expect(semType, 'the product column must carry semType Molecule so cells render as structures').toBe('Molecule');
    expect(cellRenderer, 'the grid must use the Molecule cell renderer for the product column').toBe('Molecule');
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
