import { test, expect } from './helpers';
import {
  clickDialogPrimary, clickTopMenuLeaf, openDemoCsv, resetShell, waitForDialog,
} from './helpers';
import type { Page } from '@playwright/test';

// Control comparisons (Dunnett / Holm-Welch): ML > Analyze > Group Comparison > Control Comparisons...
//
// The second test is a characterisation test for GROK-20795 — the dialog computes on the whole
// column and ignores `df.filter`. It asserts the CURRENT (wrong) behaviour on purpose, so that
// fixing the bug makes it fail loudly. When GROK-20795 lands, invert it: the group sizes must
// then shrink with the filter and sum to the visible row count.

const RESULT_TABLE = 'Control comparisons result';

async function resultTableCount(page: Page): Promise<number> {
  return page.evaluate((name: string) => Array.from((window as any).grok.shell.tables as any[])
    .filter((t) => String(t.name).startsWith(name)).length, RESULT_TABLE);
}

async function runControlComparisons(page: Page): Promise<void> {
  await clickTopMenuLeaf(page, 'div-ML---Analyze---Group-Comparison---Control-Comparisons...');
  await waitForDialog(page, 'Control comparisons');
  await clickDialogPrimary(page, ['Run', 'RUN', 'OK']);
  await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /^Control comparisons$/i }))
    .toHaveCount(0, { timeout: 15_000 });
}

/** Per-group sizes from the results table produced by the run that took the count past `before`. */
async function readGroupSizes(page: Page, before: number): Promise<number[]> {
  await page.waitForFunction((args: { name: string, before: number }) =>
    Array.from((window as any).grok.shell.tables as any[])
      .filter((t) => String(t.name).startsWith(args.name)).length > args.before,
  { name: RESULT_TABLE, before }, { timeout: 30_000 });

  return page.evaluate((name: string) => {
    const tables = Array.from((window as any).grok.shell.tables as any[])
      .filter((t) => String(t.name).startsWith(name));
    return Array.from(tables[tables.length - 1].col('n').toList() as number[]);
  }, RESULT_TABLE);
}

test.describe.serial('EDA / Control comparisons', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('Control comparisons on demog.csv produces a Box plot and the results table', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'demog.csv');

    const before = await resultTableCount(page);
    await runControlComparisons(page);

    await page.waitForFunction(() => {
      const tv = (window as any).grok?.shell?.tv;
      return !!tv && Array.from(tv.viewers).some((v: any) => /box ?plot/i.test(String(v?.type ?? '')));
    }, undefined, { timeout: 30_000 });

    const columns = await page.evaluate((name: string) => {
      const tables = Array.from((window as any).grok.shell.tables as any[])
        .filter((t) => String(t.name).startsWith(name));
      return tables[tables.length - 1].columns.names() as string[];
    }, RESULT_TABLE);

    expect(await readGroupSizes(page, before)).not.toHaveLength(0);
    expect(columns).toEqual(expect.arrayContaining(['Group', 'n', 't', 'df', 'p (raw)', 'p (adj)']));
  });

  test('GROK-20795: computes on every row, ignoring the active filter', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'demog.csv');

    const unfilteredRows = await page.evaluate(() => (window as any).grok.shell.t.filter.trueCount);

    let before = await resultTableCount(page);
    await runControlComparisons(page);
    const sizesUnfiltered = await readGroupSizes(page, before);
    expect(sizesUnfiltered.length).toBeGreaterThan(0);

    // Same `df.filter` the filter panel writes.
    await page.evaluate(() => {
      const t = (window as any).grok.shell.t;
      const sex = t.col('SEX');
      const first = sex.get(0);
      t.filter.init((i: number) => sex.get(i) === first);
    });
    const filteredRows = await page.evaluate(() => (window as any).grok.shell.t.filter.trueCount);
    expect(filteredRows).toBeLessThan(unfilteredRows);

    before = await resultTableCount(page);
    await runControlComparisons(page);
    const sizesFiltered = await readGroupSizes(page, before);

    // GROK-20795: group sizes do not shrink with the filter, and still cover more rows than are
    // visible. Invert both assertions when the dialog becomes filter-aware.
    expect(sizesFiltered).toEqual(sizesUnfiltered);
    expect(sizesFiltered.reduce((a, b) => a + b, 0)).toBeGreaterThan(filteredRows);
  });
});
