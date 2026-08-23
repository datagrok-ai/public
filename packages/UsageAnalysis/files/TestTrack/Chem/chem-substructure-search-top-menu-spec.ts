/* ---
realizes: [chem.cp.substructure-search-top-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {openChemMenuItem} from '../helpers/chem';
import {finishSpec} from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

const MOL_COL = 'canonical_smiles';
const SUBSTRUCTURE_LEAF = 'Substructure Search...';
const SMILES_INPUT = '.d4-dialog input[placeholder*="SMILES" i]';
const PICKER_DIALOG = '[name="dialog-Substructure-search"]';
const FILTER_CARD = '[name="viewer-Filters"] .d4-filter';

// The home-page PowerPack widgets fail on dev with a broken datagrok_reader DB password and
// log both the failure and its translated stack trace. Excluded by name so that everything
// else the run logs still fails the teardown check.
const AMBIENT_CONSOLE = [
  /PowerPack:MostRecentEntities/,
  /PowerPack:RecentlySharedWithMe/,
  /password authentication failed for user "datagrok_reader"/,
];

interface CardState {
  cards: number;
  columns: (string | null)[];
  type: string | null;
  molBlock: string | null;
  empty: boolean;
  isFiltering: boolean | null;
}

interface Partition {
  patternValid: boolean;
  kept: number;
  keptMisses: number;
  excluded: number;
  excludedHits: number;
  unparseable: number;
}

async function readTrueCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.t.filter.trueCount as number);
}

async function readCardState(page: Page, molCol: string): Promise<CardState> {
  return page.evaluate((col) => {
    const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
    const filter = grok.shell.tv.getFiltersGroup({createDefaultFilters: false}).filters
      .find((f: any) => (f.columnName ?? f.column) === col);
    const state = filter ? filter.saveState() : null;
    return {
      cards: cards.length,
      columns: cards.map((c) => {
        const label = c.querySelector('.d4-filter-column-name');
        return label ? label.textContent : null;
      }),
      type: state ? state.type : null,
      molBlock: state ? state.molBlock : null,
      empty: !!state && state.molBlock === DG.WHITE_MOLBLOCK,
      isFiltering: filter ? filter.isFiltering : null,
    };
  }, molCol);
}

// Re-derives the substructure predicate row by row with RDKit and reports it against the
// filter's own partition, so a row count can be checked for being the RIGHT rows.
async function substructurePartition(page: Page, molCol: string, query: string): Promise<Partition> {
  return page.evaluate(async ({col, q}) => {
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const df = grok.shell.t;
    const molecules = df.col(col);
    const bitset = df.filter;
    const pattern = rdkit.get_qmol(q);
    const p = {patternValid: false, kept: 0, keptMisses: 0, excluded: 0, excludedHits: 0, unparseable: 0};
    try {
      p.patternValid = !!pattern && !!pattern.is_valid();
      for (let i = 0; i < df.rowCount; i++) {
        let mol = null;
        try { mol = rdkit.get_mol(molecules.get(i)); }
        catch (e) { p.unparseable++; continue; }
        if (!mol) { p.unparseable++; continue; }
        let hit = false;
        try {
          const match = mol.get_substruct_match(pattern);
          hit = !!match && match !== '{}';
        }
        finally { mol.delete(); }
        if (bitset.get(i)) { p.kept++; if (!hit) p.keptMisses++; }
        else { p.excluded++; if (hit) p.excludedHits++; }
      }
    }
    finally { pattern.delete(); }
    return p;
  }, {col: molCol, q: query});
}

// Re-parses the query the PRODUCT holds — the card's own saved molBlock, which per
// references/chem.md:239 may come back as raw SMILES rather than a molblock, both of which
// get_mol accepts — and reports whether the given fragment is present in it. Lets a query typed
// into the sketcher be told apart from one appended to a field that was never cleared.
async function cardQueryContains(page: Page, molBlock: string, smarts: string): Promise<boolean> {
  return page.evaluate(async ({mb, q}) => {
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const mol = rdkit.get_mol(mb);
    const pattern = rdkit.get_qmol(q);
    try {
      const match = mol.get_substruct_match(pattern);
      return !!match && match !== '{}';
    }
    finally { mol.delete(); pattern.delete(); }
  }, {mb: molBlock, q: smarts});
}

async function typeQueryIntoOpenSketcher(
  page: Page, smiles: string, opts: {expectZero?: boolean} = {},
): Promise<number> {
  const input = page.locator(SMILES_INPUT).first();
  await input.waitFor({state: 'visible', timeout: 20_000});
  const baseline = await readTrueCount(page);
  await input.click();
  await input.fill('');
  await input.click();
  await page.keyboard.type(smiles, {delay: 30});
  await page.keyboard.press('Enter');
  await page.evaluate(() => { (window as any).__ssPrev = undefined; });
  await page.waitForFunction(
    ({base, zero}) => {
      const c = grok.shell.t.filter.trueCount as number;
      const w = window as any;
      const settled = w.__ssPrev === c;
      w.__ssPrev = c;
      return zero ? (c === 0 && settled) : (c !== base && settled);
    },
    {base: baseline, zero: opts.expectZero === true},
    {timeout: 20_000, polling: 400},
  );
  return readTrueCount(page);
}

async function okSketcher(page: Page): Promise<void> {
  await page.locator('.d4-dialog [name="button-OK"]').first().click();
  await page.locator(SMILES_INPUT).first().waitFor({state: 'detached', timeout: 15_000});
}

test.use(specTestOptions);

test('Chem: Substructure Search via top menu — match / no-match / re-invoke replaces', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);

  let total = 0;
  let benzeneCount = 0;
  let cardsBefore = 0;

  await softStep('Setup: open smiles.csv, canonical_smiles Molecule column ready', async () => {
    await page.evaluate(async () => {
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      grok.shell.addTableView(df);
      (window as any).__consoleErrors = [];
      const orig = console.error;
      const collector = function (...args: any[]) {
        (window as any).__consoleErrors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
      // Tag carried by the wrapper itself, not by the array it fills: anything that replaces
      // console.error afterwards leaves the array in place and empty, which reads as a clean run.
      (collector as any).__ssCollector = true;
      console.error = collector;
    });
    await waitForChemMenu(page);
    await waitForMolecule(page);
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30_000});
    const info = await page.evaluate((molCol) => {
      const df = grok.shell.t;
      return {
        total: df.rowCount as number,
        hasMolCol: df.columns.toList().some((c: any) => c.name === molCol && c.semType === 'Molecule'),
        molCols: df.columns.toList()
          .filter((c: any) => c.semType === 'Molecule').map((c: any) => c.name as string),
      };
    }, MOL_COL);
    total = info.total;
    console.log(`[probe] setup: rowCount=${info.total} hasMolCol=${info.hasMolCol} ` +
      `molCols=${info.molCols.length} ${JSON.stringify(info.molCols)}`);
    expect(info.hasMolCol, `${MOL_COL} must be a Molecule column`).toBe(true);
    // The branch the editor takes is chosen on this count (Chem/src/package.ts:744-765): one
    // Molecule column calls the function directly, two or more raise the picker. Pinning it
    // here and again in Scenario 2 is what makes the two invocations' different behaviour a
    // measured state change rather than two contradictory claims about the same product.
    expect(info.molCols, `the fixture must start with exactly one Molecule column, got ${JSON.stringify(info.molCols)}`)
      .toEqual([MOL_COL]);
  });

  await softStep('Scenario 1, Step 2: top-menu Chem | Search | Substructure Search seeds an empty filter card and opens its sketcher', async () => {

    await openChemMenuItem(page, SUBSTRUCTURE_LEAF);
    await page.locator(SMILES_INPUT).first().waitFor({state: 'visible', timeout: 20_000});

    expect(await page.locator(SMILES_INPUT).count(),
      'the leaf must open exactly one sketcher — a second SMILES input means a stacked dialog').toBe(1);
    expect(await page.locator(PICKER_DIALOG).count(),
      'with one Molecule column there is nothing to choose, so the first invocation must skip the column picker').toBe(0);

    await page.locator(FILTER_CARD).first().waitFor({state: 'visible', timeout: 20_000});
    const card = await readCardState(page, MOL_COL);
    console.log(`[probe] S1S2 first invocation: sketchers=${await page.locator(SMILES_INPUT).count()} ` +
      `pickers=${await page.locator(PICKER_DIALOG).count()} cards=${card.cards} ` +
      `columns=${JSON.stringify(card.columns)} type=${card.type} empty=${card.empty} ` +
      `isFiltering=${card.isFiltering} trueCount=${await readTrueCount(page)}/${total}`);
    expect(card.cards, `the command must add exactly one filter card; rendered columns: ${JSON.stringify(card.columns)}`).toBe(1);
    expect(card.columns, 'the card must belong to the molecule column').toEqual([MOL_COL]);
    expect(card.type, 'the card must be the Chem substructure filter').toBe('Chem:substructureFilter');
    expect(card.empty, `the card must be seeded with the empty molblock, got ${JSON.stringify(card.molBlock)}`).toBe(true);
    expect(card.isFiltering, 'a card seeded empty must not be filtering yet').toBe(false);
    expect(await readTrueCount(page), 'the command prepares the search — it must not filter on invocation').toBe(total);
  });

  await softStep('Scenario 1, Step 3: benzene query filters to a proper non-empty subset', async () => {
    benzeneCount = await typeQueryIntoOpenSketcher(page, 'c1ccccc1');
    expect(benzeneCount, 'benzene must match at least one row').toBeGreaterThan(0);
    expect(benzeneCount, 'benzene must exclude at least one row (proper subset)').toBeLessThan(total);
    const rowCount = await page.evaluate(() => grok.shell.t.rowCount as number);
    expect(rowCount, 'no rows deleted — only filtered').toBe(total);

    const p = await substructurePartition(page, MOL_COL, 'c1ccccc1');
    console.log(`[probe] S1S3 benzene: trueCount=${benzeneCount}/${total} ` +
      `kept=${p.kept} keptMisses=${p.keptMisses} excluded=${p.excluded} ` +
      `excludedHits=${p.excludedHits} unparseable=${p.unparseable}`);
    expect(p.unparseable, 'every molecule must parse, or the row-level check below is measuring gaps').toBe(0);
    expect(p.kept, 'the row-level check must cover the whole surviving set').toBe(benzeneCount);
    expect(p.keptMisses, 'every surviving row must genuinely contain the benzene ring').toBe(0);
    expect(p.excludedHits, 'no excluded row may contain it — the surviving set must be the substructure predicate, not an arbitrary subset').toBe(0);
  });

  await softStep('Scenario 1, Step 4: a non-matching query filters to exactly zero, rows intact', async () => {

    const zeroCount = await typeQueryIntoOpenSketcher(page, '[Au]', {expectZero: true});
    expect(zeroCount, 'a query absent from every molecule filters to zero rows').toBe(0);
    const rowCount = await page.evaluate(() => grok.shell.t.rowCount as number);
    expect(rowCount, 'the empty filter deletes no rows — the dataset is intact').toBe(total);

    // The clear before typing is a single fill(''); if it no-ops the field holds c1ccccc1[Au],
    // which also matches zero rows, and the row-level check below would re-derive against the
    // literal '[Au]' instead of the query the product actually holds.
    const cardGold = await readCardState(page, MOL_COL);
    expect(cardGold.molBlock,
      'the card must hold a saved query to check the typed query against').toBeTruthy();
    const holdsGold = await cardQueryContains(page, cardGold.molBlock!, '[Au]');
    const holdsBenzene = await cardQueryContains(page, cardGold.molBlock!, 'c1ccccc1');
    console.log(`[probe] S1S4 card query: holdsGold=${holdsGold} holdsBenzene=${holdsBenzene} ` +
      `query=${JSON.stringify((cardGold.molBlock ?? '').slice(0, 120))}`);
    expect(holdsGold, 'the card must hold the gold query that was just typed').toBe(true);
    expect(holdsBenzene,
      'the sketcher must have been cleared before [Au] was typed — a leftover benzene ring means the row-level check below re-derives a query the product does not hold').toBe(false);

    // Zero is the right answer only if no molecule actually contains gold; without this the
    // count cannot be told apart from a search that silently failed and left the set empty.
    const p = await substructurePartition(page, MOL_COL, '[Au]');
    console.log(`[probe] S1S4 gold: trueCount=${zeroCount} rowCount=${rowCount}/${total} ` +
      `patternValid=${p.patternValid} excluded=${p.excluded} excludedHits=${p.excludedHits} ` +
      `unparseable=${p.unparseable}`);
    // A query RDKit failed to parse matches nothing, so excludedHits would be zero for a reason
    // that has nothing to do with the data.
    expect(p.patternValid, 'the [Au] query must parse into a valid RDKit query molecule, or the zero hit count below holds vacuously').toBe(true);
    expect(p.unparseable, 'every molecule must parse, or the row-level check below is measuring gaps').toBe(0);
    expect(p.excluded, 'the whole dataset must be excluded').toBe(total);
    expect(p.excludedHits, 'no molecule contains gold, so an empty result is the correct answer and not a silent failure').toBe(0);
  });

  await softStep('Scenario 2, Step 1: re-invoking the top menu replaces (does not stack) the filter', async () => {

    await okSketcher(page);
    const entryCount = await readTrueCount(page);
    console.log(`[probe] S2S1 entry state: trueCount=${entryCount}/${total}`);
    expect(entryCount,
      'Scenario 2 must enter from the zero-row state Scenario 1 left behind — otherwise the reset below is measured against an unknown starting point').toBe(0);
    cardsBefore = await page.locator(FILTER_CARD).count();
    // Positive control: with no Filter Panel rendered both counts would be 0 and the
    // no-second-card check below would pass on there being no cards at all.
    expect(cardsBefore, 'the Filter Panel must be rendering the card for the card-count comparison to mean anything').toBeGreaterThan(0);
    // Why this invocation behaves differently from Scenario 1's: running a search caches
    // helper columns on the dataframe (saveColumn, Chem/src/chem-searches.ts:176-187), one of
    // which holds canonical SMILES and is itself typed Molecule. The editor's branch is chosen
    // on the Molecule-column count, so it now raises the picker instead of calling directly.
    // Asserting the count here is what keeps that from reading as a contradiction with S1S2.
    const molColsNow = await page.evaluate(() => grok.shell.t.columns.toList()
      .filter((c: any) => c.semType === 'Molecule').map((c: any) => c.name as string));
    console.log(`[probe] S2S1 before re-invoke: molCols=${molColsNow.length} ` +
      `${JSON.stringify(molColsNow)} cardsBefore=${cardsBefore}`);
    expect(molColsNow.length,
      `the searches already run must have cached a second Molecule-typed column, or the picker ` +
      `below cannot appear; columns seen: ${JSON.stringify(molColsNow)}`).toBeGreaterThan(1);

    await openChemMenuItem(page, SUBSTRUCTURE_LEAF);
    const picker = page.locator(`${PICKER_DIALOG} [name="button-OK"]`).first();
    await picker.waitFor({state: 'visible', timeout: 20_000});

    expect(await page.locator(PICKER_DIALOG).count(), 're-invoking the top-menu leaf must open the column-picker dialog').toBe(1);
    await picker.click();
    await page.locator(SMILES_INPUT).first().waitFor({state: 'visible', timeout: 20_000});
    await page.waitForFunction((t) => grok.shell.t.filter.trueCount === t, total, {timeout: 20_000});
    const cardsAfter = await page.locator(FILTER_CARD).count();
    const cardAfter = await readCardState(page, MOL_COL);
    const resetCount = await readTrueCount(page);
    console.log(`[probe] S2S1 after OK: cardsAfter=${cardsAfter} (before=${cardsBefore}) ` +
      `columns=${JSON.stringify(cardAfter.columns)} type=${cardAfter.type} empty=${cardAfter.empty} ` +
      `isFiltering=${cardAfter.isFiltering} molBlock=${JSON.stringify((cardAfter.molBlock ?? '').slice(0, 120))} ` +
      `trueCount=${resetCount}/${total}`);
    expect(cardsAfter, 're-invoking must not add a second substructure filter card').toBe(cardsBefore);
    // The picker offers the user-visible column, not the cached ~-prefixed helper; confirming
    // it must therefore leave the filter on canonical_smiles.
    expect(cardAfter.columns,
      `confirming the picker must prepare the filter on the user-visible molecule column, got ${JSON.stringify(cardAfter.columns)}`)
      .toEqual([MOL_COL]);
    expect(cardAfter.empty,
      `confirming the picker must re-prepare the card empty, got ${JSON.stringify(cardAfter.molBlock)}`).toBe(true);
    expect(resetCount, 're-invoke resets the filter to the full table before the new query').toBe(total);
  });

  await softStep('Scenario 2, Step 2: carboxylic-acid query yields a different non-empty subset', async () => {
    const acidCount = await typeQueryIntoOpenSketcher(page, 'C(=O)O');
    await okSketcher(page);
    expect(acidCount, 'carboxylic acid must match at least one row').toBeGreaterThan(0);
    expect(acidCount, 'carboxylic acid must be a proper subset').toBeLessThan(total);
    expect(acidCount, 'the new query replaced the benzene filter — counts must differ').not.toBe(benzeneCount);

    // Read AFTER the commit, which is the event the reuse claim is about: a second card added on
    // commit and left empty ANDs to exactly the acid row set, so every count below still holds.
    const cardCommit = await readCardState(page, MOL_COL);
    console.log(`[probe] S2S2 card after commit: cards=${cardCommit.cards} (before=${cardsBefore}) ` +
      `columns=${JSON.stringify(cardCommit.columns)} type=${cardCommit.type} empty=${cardCommit.empty} ` +
      `isFiltering=${cardCommit.isFiltering}`);
    expect(cardCommit.cards,
      `committing the new query must update the existing card, not add a second one; rendered columns: ${JSON.stringify(cardCommit.columns)}`)
      .toBe(cardsBefore);
    expect(cardCommit.isFiltering,
      'the card on the molecule column must be the one filtering — a second, empty card ahead of it would filter nothing').toBe(true);

    const p = await substructurePartition(page, MOL_COL, 'C(=O)O');
    console.log(`[probe] S2S2 acid: trueCount=${acidCount}/${total} benzeneWas=${benzeneCount} ` +
      `kept=${p.kept} keptMisses=${p.keptMisses} excludedHits=${p.excludedHits} unparseable=${p.unparseable}`);
    expect(p.unparseable, 'every molecule must parse, or the row-level check below is measuring gaps').toBe(0);
    expect(p.kept, 'the row-level check must cover the whole surviving set').toBe(acidCount);
    expect(p.keptMisses, 'every surviving row must genuinely contain the carboxylic acid fragment').toBe(0);
    expect(p.excludedHits, 'no excluded row may contain it — the benzene filter must be gone, not intersected').toBe(0);
  });

  await softStep('Teardown: the run logs no console error of its own', async () => {
    const collectorInstalled = await page.evaluate(() => Array.isArray((window as any).__consoleErrors) &&
      (console.error as any).__ssCollector === true);
    console.log(`[probe] teardown: collectorInstalled=${collectorInstalled}`);
    expect(collectorInstalled,
      'the collector installed in Setup must still BE console.error — an absent array, or a console.error replaced after Setup, reads exactly like a clean run').toBe(true);
    const errs = await page.evaluate(() => ((window as any).__consoleErrors ?? []) as string[]);
    // A translated stack trace is logged separately from its message, keyed by the id the
    // message carries — drop those alongside the ambient messages they belong to.
    const ambientIds = errs
      .filter((e) => AMBIENT_CONSOLE.some((r) => r.test(e)))
      .map((e) => (e.match(/ID = (\w+)/) ?? [])[1])
      .filter((id) => !!id);
    const relevant = errs.filter((e) => !AMBIENT_CONSOLE.some((r) => r.test(e)) &&
      !ambientIds.some((id) => e.startsWith(`Stack trace ${id}`)));
    expect(relevant.length, `console errors during the run: ${JSON.stringify(relevant.slice(0, 5))}`).toBe(0);
    await page.evaluate(() => grok.shell.closeAll());
  });

  finishSpec();
});
