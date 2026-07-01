/* ---
sub_features_covered: [general.toolbox-search, general.toolbox-search-compound, general.toolbox-search-date, general.toolbox-search-whitespace, general.toolbox-search-quoted, general.toolbox-search-null-date-crash]
--- */

import {test, expect, Page, Locator} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

// Expands the Toolbox > Search pane (if collapsed) and returns the search <input>.
async function openSearch(page: Page): Promise<Locator> {
  const input = page.locator('input[placeholder="Search"]');
  if (!(await input.isVisible().catch(() => false)))
    await page.locator('[name="div-section--Search"]').first().click();
  await expect(input).toBeVisible({timeout: 10_000});
  return input;
}

// Resets the filter, types `query` into the Search box via the UI, presses Enter, returns trueCount.
async function runSearch(page: Page, input: Locator, query: string): Promise<number> {
  await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.setAll(true));
  await input.click();
  await input.fill(query);
  await input.press('Enter');
  await page.waitForTimeout(300); // let the filter recompute
  return page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.trueCount as number);
}

test('Toolbox Search — compound and date-value filtering (GROK-20229)', async ({page}) => {
  test.setTimeout(300_000);

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await page.evaluate(() => {
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  let input!: Locator;

  await softStep('Setup: open demog.csv and expand Toolbox > Search', async () => {
    const rows = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      return df.rowCount as number;
    }, demogPath);
    expect(rows).toBe(5850);
    input = await openSearch(page);
  });

  await softStep('Single conditions match rows (sanity)', async () => {
    expect(await runSearch(page, input, 'AGE > 50'), 'AGE > 50').toBe(2176);
    expect(await runSearch(page, input, 'SEX = M'), 'SEX = M').toBe(2607);
    expect(await runSearch(page, input, 'CONTROL = true'), 'CONTROL = true').toBe(39);
  });

  await softStep('BUG 1a: compound AND filters to intersection, not 0', async () => {
    // Each condition alone matches; AND should be the intersection (856). Current dev returns 0.
    expect(await runSearch(page, input, 'AGE > 50 and SEX = M'), 'AGE > 50 and SEX = M').toBe(856);
  });

  await softStep('BUG 1b: compound OR filters to union, not 0', async () => {
    expect(await runSearch(page, input, 'AGE > 50 or SEX = M'), 'AGE > 50 or SEX = M').toBe(3927);
  });

  await softStep('BUG 1c: tautology OR keeps all rows, not 0', async () => {
    expect(await runSearch(page, input, 'SEX = M or SEX = F'), 'SEX = M or SEX = F').toBe(5850);
  });

  await softStep('BUG 2: year-only date value matches same rows as explicit date', async () => {
    const explicit = await runSearch(page, input, 'STARTED > 1/1/1991');
    expect(explicit, 'STARTED > 1/1/1991 (explicit date)').toBeGreaterThan(2000);
    const yearOnly = await runSearch(page, input, 'STARTED > 1991');
    expect(yearOnly, 'STARTED > 1991 (year only) should match like the explicit date').toBeGreaterThan(2000);
  });

  await softStep('BUG 4: leading/trailing whitespace is trimmed, not treated as a non-match', async () => {
    const trimmed = await runSearch(page, input, 'Asian');
    expect(trimmed, 'Asian (substring, incl. Caucasian)').toBe(5339);
    // A user who types a stray space around the term should get the same matches, not 0.
    expect(await runSearch(page, input, ' Asian'), 'leading space " Asian"').toBe(5339);
    expect(await runSearch(page, input, 'Asian '), 'trailing space "Asian "').toBe(5339);
    expect(await runSearch(page, input, '  Asian  '), 'padded "  Asian  "').toBe(5339);
  });

  await softStep('BUG 5: quoted equality value matches the same as the unquoted form', async () => {
    const unquoted = await runSearch(page, input, 'RACE = Asian');
    expect(unquoted, 'RACE = Asian (exact category)').toBe(72);
    // Quoting a string value must not change the match — quotes should not be taken literally.
    expect(await runSearch(page, input, 'RACE = "Asian"'), 'RACE = "Asian" (quoted)').toBe(72);
  });

  await softStep('BUG 3 (GROK-20229): date comparison on a datetime column with nulls must not crash', async () => {
    // Build a synthetic table whose datetime column contains null values. The year is irrelevant —
    // the nulls are the trigger. A valid date comparison currently throws an uncaught TypeError
    // ("Cannot read properties of undefined (reading 'a')") and filters to 0.
    const errsBefore = consoleErrors.length;
    const total = await page.evaluate(() => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      grok.shell.closeAll();
      const n = 200;
      const dates: (Date | null)[] = [];
      for (let i = 0; i < n; i++)
        dates.push(i % 7 === 0 ? null : new Date(2019, i % 12, (i % 27) + 1));
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList('datetime', 'EVENT_DATE', dates),
        DG.Column.fromList('int', 'NUM', dates.map((_: unknown, i: number) => i)),
      ]);
      df.name = 'synth2019';
      grok.shell.addTableView(df);
      return df.rowCount as number;
    });
    expect(total).toBe(200);
    const synthInput = await openSearch(page);
    // Sanity: a non-date search on the same table works and filters. NUM is 0..199, so NUM > 50 = 149.
    expect(await runSearch(page, synthInput, 'NUM > 50'), 'NUM > 50 (sanity)').toBe(149);
    // The date comparison must return a sane count (dates are all in 2019) and, above all, not crash.
    const dated = await runSearch(page, synthInput, 'EVENT_DATE > 1/1/2019');
    expect(dated, 'EVENT_DATE > 1/1/2019 on a column with nulls').toBeGreaterThan(0);
    const newErrors = consoleErrors.slice(errsBefore);
    expect(newErrors, `date search over a null-containing datetime column crashed:\n${newErrors.join('\n')}`).toEqual([]);
  });

  await softStep('No fatal console errors during search', async () => {
    expect(consoleErrors, consoleErrors.join('\n')).toEqual([]);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
