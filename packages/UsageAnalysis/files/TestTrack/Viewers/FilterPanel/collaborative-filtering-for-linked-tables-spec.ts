/* ---
realizes: [ filters.cp.linked-tables-collaborative, GROK-19137 ]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const LINK_DIALOG = '.d4-dialog[name="dialog-Link-Tables"]';
const LINK_TYPE_EDITOR = `${LINK_DIALOG} select[name="input-Link-Type"]`;
const BYSTANDER_SETTLE_MS = 1500;

interface FrameState {
  rowCount: number;
  filtered: number;
  selected: number;
}

interface Fixture {
  linkedRowCount: number;
  linkAndFilterCount: number;
  bothCriteriaCount: number;
  noSelectionBaseline: number;
  linked1RowCount: number;
}

async function readFrame(page: Page, name: string): Promise<FrameState | null> {
  return page.evaluate((n) => {
    for (const view of (window as any).grok.shell.tableViews) {
      const df = view.dataFrame;
      if (df.name === n)
        return {rowCount: df.rowCount, filtered: df.filter.trueCount, selected: df.selection.trueCount};
    }
    return null;
  }, name);
}

async function frameState(page: Page, name: string): Promise<FrameState> {
  const state = await readFrame(page, name);
  expect(state, `${name} is no longer open as a table view, so its counts cannot be read`).not.toBeNull();
  return state!;
}

async function expectUnmoved(page: Page, name: string, expected: FrameState, why: string): Promise<void> {
  expect(await readFrame(page, name), why).toEqual(expected);
  await page.waitForTimeout(BYSTANDER_SETTLE_MS);
  expect(await readFrame(page, name), `${why} (second reading, ${BYSTANDER_SETTLE_MS} ms later)`)
    .toEqual(expected);
}

async function filteredCountOf(page: Page, name: string): Promise<number> {
  return (await readFrame(page, name))?.filtered ?? -1;
}

async function survivingCategories(
  page: Page, frameName: string, column: string,
): Promise<string[] | null> {
  return page.evaluate(({n, col}) => {
    for (const view of (window as any).grok.shell.tableViews) {
      const df = view.dataFrame;
      if (df.name !== n) continue;
      const c = df.col(col);
      if (!c) return null;
      const seen = new Set<string>();
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i)) seen.add(String(c.get(i)));
      return Array.from(seen);
    }
    return null;
  }, {n: frameName, col: column});
}

async function deriveFixture(page: Page): Promise<Fixture> {
  return page.evaluate(() => {
    const byName = (n: string) => {
      for (const view of (window as any).grok.shell.tableViews)
        if (view.dataFrame.name === n) return view.dataFrame;
      throw new Error(`no open table view holds ${n}, so the fixture cannot be derived from its columns`);
    };
    const master = byName('spgi-100');
    const linked1 = byName('SPGI-linked1');
    const linked2 = byName('SPGI-linked2');
    const KEYS = ['Sample Name', 'link column 1', 'link column 2', 'link column 3'];
    const keyCols1 = KEYS.map((k) => linked1.col(k));
    const keyCols2 = KEYS.map((k) => linked2.col(k));
    const keyOf = (cols: any[], i: number) => cols.map((c) => String(c.get(i))).join('');

    const idCol = master.col('Id');
    const top5 = new Set<string>();
    for (let i = 0; i < 5; i++) top5.add(String(idCol.get(i)));

    const vIIKeys = new Set<string>();
    const l2c3 = linked2.col('link column 3');
    for (let i = 0; i < linked2.rowCount; i++)
      if (String(l2c3.get(i)) === 'v ii') vIIKeys.add(keyOf(keyCols2, i));

    const concept = linked1.col('Concept Id');
    const pampa = linked1.col('PAMPA Classification');
    let linkedRowCount = 0;
    let linkAndFilterCount = 0;
    let bothCriteriaCount = 0;
    let noSelectionBaseline = 0;
    for (let i = 0; i < linked1.rowCount; i++) {
      const bySelection = top5.has(String(concept.get(i)));
      const byLink = vIIKeys.has(keyOf(keyCols1, i));
      const byPanel = String(pampa.get(i)) === 'inconclusive';
      if (bySelection) linkedRowCount++;
      if (bySelection && byLink) linkAndFilterCount++;
      if (bySelection && byLink && byPanel) bothCriteriaCount++;
      if (byLink && byPanel) noSelectionBaseline++;
    }
    return {
      linkedRowCount, linkAndFilterCount, bothCriteriaCount, noSelectionBaseline,
      linked1RowCount: linked1.rowCount,
    };
  });
}

async function switchToView(page: Page, name: string): Promise<void> {
  const found = await page.evaluate((n) => {
    for (const view of (window as any).grok.shell.tableViews)
      if (view.dataFrame.name === n) { (window as any).grok.shell.v = view; return true; }
    return false;
  }, name);
  expect(found,
    `no open table view holds ${name}, so switching to it was a silent no-op and every gesture in ` +
    'this step would land on whatever view happened to be active')
    .toBe(true);
  await expect.poll(() => page.evaluate(() => (window as any).grok.shell.tv?.dataFrame?.name ?? ''), {
    timeout: 15_000,
    message: `grok.shell.tv never became ${name} after the view switch`,
  }).toBe(name);
}

async function openFilterPanelViaRibbon(page: Page): Promise<void> {
  await page.locator('.d4-ribbon-panel [name="icon-filter"]').first().click();
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15_000});
}

async function selectTopFiveRows(page: Page): Promise<FrameState & {name: string}> {
  return page.evaluate(() => {
    const df = (window as any).grok.shell.tv.dataFrame;
    df.selection.setAll(false);
    for (let i = 0; i < 5; i++) df.selection.set(i, true);
    return {
      name: df.name, rowCount: df.rowCount,
      filtered: df.filter.trueCount, selected: df.selection.trueCount,
    };
  });
}

async function linkTypeEditorsWithValue(page: Page, value: string): Promise<number[]> {
  return page.evaluate(({sel, val}) => {
    const editors = Array.from(document.querySelectorAll(sel)) as HTMLSelectElement[];
    return editors.map((e, i) => ({i, v: e.value})).filter((x) => x.v === val).map((x) => x.i);
  }, {sel: LINK_TYPE_EDITOR, val: value});
}

test('Collaborative Filtering for Linked Tables', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  let fx: Fixture = {
    linkedRowCount: -1, linkAndFilterCount: -1, bothCriteriaCount: -1,
    noSelectionBaseline: -1, linked1RowCount: -1,
  };
  const masterSelected: FrameState = {rowCount: 100, filtered: 100, selected: 5};
  const masterCleared: FrameState = {rowCount: 100, filtered: 100, selected: 0};
  let linked2Rest: FrameState = {rowCount: 224, filtered: 224, selected: 0};

  try {
    await softStep('Setup: open 3 tables and link them', async () => {
      const result = await page.evaluate(async () => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        document.body.classList.add('selenium');
        grok.shell.settings.showFiltersIconsConstantly = true;
        grok.shell.windows.simpleMode = true;
        grok.shell.closeAll();

        const df1 = await grok.data.files.openTable('System:AppData/Chem/tests/spgi-100.csv');
        const df2 = await grok.data.files.openTable('System:AppData/ApiTests/datasets/SPGI-linked1.csv');
        const df3 = await grok.data.files.openTable('System:AppData/ApiTests/datasets/SPGI-linked2.csv');

        grok.data.linkTables(df3, df2,
          ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
          ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
          [DG.SYNC_TYPE.FILTER_TO_FILTER]);
        grok.data.linkTables(df1, df2, ['Id'], ['Concept Id'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);

        grok.shell.addTableView(df1);
        grok.shell.addTableView(df2);
        grok.shell.addTableView(df3);

        for (const df of [df1, df2, df3]) {
          await new Promise((resolve) => {
            const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
            setTimeout(resolve, 5000);
          });
        }
        const hasBioChem = [df1, df2, df3].some((df) =>
          Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
            .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule'));
        if (hasBioChem) {
          for (let i = 0; i < 50; i++) {
            if (document.querySelector('[name="viewer-Grid"] canvas')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 5000));
        }
        const views = Array.from((window as any).grok.shell.tableViews).map((tv: any) => tv.dataFrame.name);
        return {df1: df1.rowCount, df2: df2.rowCount, df3: df3.rowCount,
          names: [df1.name, df2.name, df3.name], views};
      });
      expect(result.df1, 'spgi-100 did not load its full 100 rows — the fixture changed').toBe(100);
      expect(result.df2, 'SPGI-linked1 did not load its full 3624 rows — the fixture changed').toBe(3624);
      expect(result.df3, 'SPGI-linked2 did not load its full 224 rows — the fixture changed').toBe(224);
      expect(result.names, 'the three frames are not named as the later steps address them').toEqual(
        ['spgi-100', 'SPGI-linked1', 'SPGI-linked2']);
      for (const name of ['spgi-100', 'SPGI-linked1', 'SPGI-linked2'])
        expect(result.views, `no table view holds ${name}, so its counts cannot be read`).toContain(name);

      fx = await deriveFixture(page);
      console.log(`derived fixture: ${JSON.stringify(fx)}`);
      expect(fx.linkedRowCount,
        `no SPGI-linked1 row carries a Concept Id among the first five spgi-100 Id values, so the ` +
        'selection-to-filter link has nothing to narrow to and Step 2 would assert an empty set')
        .toBeGreaterThan(0);
      expect(fx.linkedRowCount,
        `the first five spgi-100 Id values match all ${fx.linked1RowCount} SPGI-linked1 rows, so the ` +
        'selection link would not be a narrowing at all')
        .toBeLessThan(fx.linked1RowCount);
      expect(fx.linkAndFilterCount,
        'no SPGI-linked1 row is both key-matched to the selected spgi-100 rows and carried by a "v ii" ' +
        'SPGI-linked2 row, so Step 6 would assert an empty set')
        .toBeGreaterThan(0);
      expect(fx.linkAndFilterCount,
        `the "v ii" criterion drops none of the ${fx.linkedRowCount} key-matched rows, so Step 6 could ` +
        'not tell the filter-to-filter link apart from the selection link alone')
        .toBeLessThan(fx.linkedRowCount);
      expect(fx.bothCriteriaCount,
        'no SPGI-linked1 row satisfies selection + "v ii" + "inconclusive" at once, so Step 8 would ' +
        'assert an empty set')
        .toBeGreaterThan(0);
      expect(fx.bothCriteriaCount,
        `the "inconclusive" criterion drops none of the ${fx.linkAndFilterCount} rows Step 6 leaves, so ` +
        'Step 8 would not show the panel criterion composing with the link-driven narrowing')
        .toBeLessThan(fx.linkAndFilterCount);
      expect(fx.noSelectionBaseline,
        `the derived no-selection-link count ${fx.noSelectionBaseline} does not exceed the ` +
        `${fx.bothCriteriaCount} rows of Step 8, so releasing the selection would show no movement`)
        .toBeGreaterThan(fx.bothCriteriaCount);
      expect(fx.noSelectionBaseline,
        `the derived no-selection-link count covers all ${fx.linked1RowCount} SPGI-linked1 rows, so ` +
        'Step 9 could not tell "the link released" apart from "every filter was dropped"')
        .toBeLessThan(fx.linked1RowCount);
    });

    await softStep('Step 1: Select the top 5 rows in spgi-100', async () => {
      await switchToView(page, 'spgi-100');
      const result = await selectTopFiveRows(page);
      expect(result.name, 'the active view is not spgi-100 — the selection landed on the wrong frame')
        .toBe('spgi-100');
      expect(result.rowCount, 'the active view is not spgi-100 — the selection landed on the wrong frame')
        .toBe(100);
      expect(result.selected, 'the top 5 rows are not selected, so nothing can propagate over the link')
        .toBe(5);
      expect(result.filtered, 'spgi-100 narrowed itself on a selection — selection leaked into its own filter')
        .toBe(result.rowCount);
    });

    await softStep('Step 2: SPGI-linked1 shows 9 filtered rows', async () => {
      await switchToView(page, 'SPGI-linked1');
      await expect.poll(() => filteredCountOf(page, 'SPGI-linked1'), {
        timeout: 15_000,
        message: 'the selection-to-filter link did not narrow SPGI-linked1 to the key-matched rows',
      }).toBe(fx.linkedRowCount);
      const view = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.name);
      expect(view, 'the active view is not SPGI-linked1 — the later panel steps would target the wrong frame')
        .toBe('SPGI-linked1');
      await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
        'SPGI-linked2 moved when the spgi-100 selection propagated to SPGI-linked1 — it is on the other ' +
        'end of a filter-to-filter link and nothing was asked of it');
      await expectUnmoved(page, 'spgi-100', masterSelected,
        'spgi-100 moved while its own selection propagated outward — the selection-to-filter link fed ' +
        'back into its master');
    });

    await softStep('Step 5: Filter link column 3 to "v ii" on SPGI-linked2', async () => {
      await switchToView(page, 'SPGI-linked2');
      await openFilterPanelViaRibbon(page);
      const {filteredCount} = await v.applyCategoricalFilter(page, 'link column 3', ['v ii']);
      expect(filteredCount, 'the "v ii" criterion filtered SPGI-linked2 down to nothing')
        .toBeGreaterThan(0);
      expect(filteredCount, 'the "v ii" criterion did not narrow SPGI-linked2 at all').toBeLessThan(224);
      const categories = await survivingCategories(page, 'SPGI-linked2', 'link column 3');
      expect(categories, 'SPGI-linked2 or its "link column 3" column is missing').not.toBeNull();
      expect(categories, 'rows outside the "v ii" category survived the criterion').toEqual(['v ii']);
      linked2Rest = await frameState(page, 'SPGI-linked2');
      await expectUnmoved(page, 'spgi-100', masterSelected,
        'spgi-100 moved when a criterion was applied two links away on SPGI-linked2 — no link runs in ' +
        'that direction');
    });

    await softStep('Step 6: SPGI-linked1 shows 5 filtered rows', async () => {
      await switchToView(page, 'SPGI-linked1');
      await expect.poll(() => filteredCountOf(page, 'SPGI-linked1'), {
        timeout: 15_000,
        message: 'the filter-to-filter link did not propagate the SPGI-linked2 criterion down to SPGI-linked1',
      }).toBe(fx.linkAndFilterCount);
      await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
        'SPGI-linked2 changed while its own criterion propagated to SPGI-linked1 — the filter-to-filter ' +
        'link fed back into its source');
      await expectUnmoved(page, 'spgi-100', masterSelected,
        'spgi-100 changed while SPGI-linked1 was narrowed by the link — nothing propagates back to the ' +
        'selection master');
    });

    await softStep('Step 8: PAMPA Classification = "inconclusive" — 2 rows', async () => {
      await openFilterPanelViaRibbon(page);
      const {filteredCount} = await v.applyCategoricalFilter(page, 'PAMPA Classification', ['inconclusive']);
      expect(filteredCount,
        'the panel criterion did not compose with the link-driven narrowing on SPGI-linked1')
        .toBe(fx.bothCriteriaCount);
      const categories = await survivingCategories(page, 'SPGI-linked1', 'PAMPA Classification');
      expect(categories, 'SPGI-linked1 or its "PAMPA Classification" column is missing').not.toBeNull();
      expect(categories, 'rows outside the "inconclusive" category survived the panel criterion')
        .toEqual(['inconclusive']);
      await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
        'SPGI-linked2 narrowed when a panel criterion was applied on SPGI-linked1 — the filter-to-filter ' +
        'link ran backwards, which is the cross-contamination the link direction exists to prevent');
      await expectUnmoved(page, 'spgi-100', masterSelected,
        'spgi-100 moved when a panel criterion was applied on SPGI-linked1 — the selection-to-filter ' +
        'link ran backwards into its master');
    });

    await softStep('Step 9: Clearing the spgi-100 selection releases only the link-driven narrowing',
      async () => {
        await switchToView(page, 'spgi-100');
        await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.selection.setAll(false));
        await expect.poll(() => filteredCountOf(page, 'SPGI-linked1'), {
          timeout: 15_000,
          message: 'clearing the spgi-100 selection did not widen SPGI-linked1 — either the link never ' +
            'narrowed it or the release did not propagate, and Step 11 could not tell those apart',
        }).toBeGreaterThan(fx.bothCriteriaCount);
        const linked1 = await frameState(page, 'SPGI-linked1');
        expect(linked1.filtered,
          `SPGI-linked1 settled at ${linked1.filtered} rows once the selection was released, not at the ` +
          `${fx.noSelectionBaseline} rows its "v ii" and "inconclusive" criteria alone select — the ` +
          'release either dropped criteria it should have kept or kept narrowing it should have dropped')
          .toBe(fx.noSelectionBaseline);
        await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
          'SPGI-linked2 moved when the spgi-100 selection was cleared — it is on neither end of that link');
        await expectUnmoved(page, 'spgi-100', masterCleared,
          'spgi-100 did not settle at zero selected rows over its full 100 rows after its selection was ' +
          'cleared');
      });

    await softStep('Step 10: Re-selecting the same 5 rows narrows SPGI-linked1 back to 2', async () => {
      await switchToView(page, 'spgi-100');
      const result = await selectTopFiveRows(page);
      expect(result.name,
        'the re-selection landed on a frame other than spgi-100, so "5 rows are selected" says nothing ' +
        'about the link under test')
        .toBe('spgi-100');
      expect(result.rowCount, 'the active view is not spgi-100 — the re-selection landed on the wrong frame')
        .toBe(100);
      expect(result.selected, 'the top 5 rows of spgi-100 are not selected again').toBe(5);
      await expect.poll(() => filteredCountOf(page, 'SPGI-linked1'), {
        timeout: 15_000,
        message: 'the selection link no longer narrows SPGI-linked1, so "the narrowing is gone" in Step 11 ' +
          'would prove nothing',
      }).toBe(fx.bothCriteriaCount);
      await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
        'SPGI-linked2 moved when the spgi-100 selection was restored — it is on neither end of that link');
      await expectUnmoved(page, 'spgi-100', masterSelected,
        'spgi-100 narrowed itself on its own selection — selection leaked into its own filter');
    });

    await softStep('Step 11: Changing the link type clears the previous type\'s filtering (GROK-19137)',
      async () => {
        const opened = await v.driveTopMenuLeaf(page, ['Data', 'Link Tables...']);
        expect(opened, 'the Data > Link Tables... menu leaf could not be actuated').toBe(true);
        await page.locator(LINK_DIALOG).waitFor({timeout: 15_000});

        const tab = page.locator(`${LINK_DIALOG} .d4-tab-header`)
          .filter({hasText: 'spgi-100 -> SPGI-linked1'});
        expect(await tab.count(),
          'the Link Tables dialog offers no tab for the spgi-100 -> SPGI-linked1 link')
          .toBeGreaterThan(0);
        await tab.first().click();
        await page.locator(LINK_TYPE_EDITOR).first().waitFor({timeout: 10_000});

        const before = await linkTypeEditorsWithValue(page, 'selection to filter');
        expect(before,
          'expected exactly one Link Type editor holding "selection to filter" — the spgi-100 link tab ' +
          'is not the one on screen, or another editor holds the same value and would be driven instead')
          .toHaveLength(1);
        await page.locator(LINK_TYPE_EDITOR).nth(before[0]).selectOption('selection to selection');

        await expect.poll(() => linkTypeEditorsWithValue(page, 'selection to selection'), {
          timeout: 10_000,
          message: 'no Link Type editor reads "selection to selection" — the change did not commit',
        }).not.toHaveLength(0);
        await expect.poll(() => linkTypeEditorsWithValue(page, 'selection to filter'), {
          timeout: 10_000,
          message: 'a Link Type editor still reads "selection to filter" — the old type was not replaced',
        }).toHaveLength(0);

        await expect.poll(async () => (await readFrame(page, 'SPGI-linked1'))?.selected ?? -1, {
          timeout: 15_000,
          message: 'SPGI-linked1 selection does not hold the key-matched rows — the new ' +
            'selection-to-selection link was not applied',
        }).toBe(fx.linkedRowCount);

        await expect.poll(() => filteredCountOf(page, 'SPGI-linked1'), {
          timeout: 15_000,
          message: 'SPGI-linked1 did not return to the no-selection-link count derived from its own ' +
            'columns after the link type changed — the previous link type\'s filtering was left in ' +
            'place (GROK-19137)',
        }).toBe(fx.noSelectionBaseline);
        const linked1 = await frameState(page, 'SPGI-linked1');
        expect(linked1.filtered,
          'SPGI-linked1 is still carrying the selection-to-filter narrowing after the link type changed ' +
          '(GROK-19137)')
          .not.toBe(fx.bothCriteriaCount);
        await expectUnmoved(page, 'SPGI-linked2', linked2Rest,
          'SPGI-linked2 moved when the spgi-100 -> SPGI-linked1 link type changed — it is on neither end ' +
          'of that link');
        await expectUnmoved(page, 'spgi-100', masterSelected,
          'spgi-100 lost its own 5-row selection or narrowed itself when the link type changed — the new ' +
          'selection-to-selection link fed back into its master');
      });
  } finally {
    await softStep('Teardown: close the Link Tables dialog and the probe tables', async () => {
      const closeButton = page.locator(`${LINK_DIALOG} [name="button-CLOSE"]`);
      if (await closeButton.count() > 0)
        await closeButton.first().click();
      await expect.poll(async () => page.locator(LINK_DIALOG).count(), {
        timeout: 10_000,
        message: 'the Link Tables dialog stayed open and would bleed into the next spec',
      }).toBe(0);
      await page.evaluate(() => (window as any).grok.shell.closeAll());
      await expect.poll(async () => page.evaluate(() =>
        Array.from((window as any).grok.shell.tableViews).length), {
        timeout: 10_000,
        message: 'the probe table views did not close, so the probe links stay attached',
      }).toBe(0);
    });
  }

  v.finishSpec();
});
