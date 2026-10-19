/* Applying a model and what it does to the user's tables: the predictions against the stored training metrics, the
   cells a scenario empties, the API's application, and claims about a named open table while another view is current.
   In-page lookups go through the library's runtime (`__bdd.table()`, `tableNamed()`, `col()`), which throw with the
   tables or columns there are. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, gestures, locate, viewers} from '@datagrok-libraries/bdd/runtime';
import {guardIdentity, listOf, modelNumber} from './grok-s.js';
import {contextPanelTitle} from './results-state.js';

declare const grok: any;
declare const DG: any;

/** A statistic of the actual column [a] and the predicted column [b] of the open table [table], computed in the page
 * over the rows where [a] has a value: the share of equal values, or the mean absolute difference. An empty [b] cell
 * in such a row is the answer instead: a skipped prediction would otherwise leave the statistic of the rest. */
function pairStatistic(page: Page, table: string, a: string, b: string, kind: 'agree' | 'mae'): Promise<number | string> {
  return viewers.evaluate(page, ([t, x, y, k]: [string, string, string, string]) => {
    const bdd = (window as any).__bdd;
    const df = bdd.tableNamed(t);
    const [ca, cb] = [bdd.col(x, df), bdd.col(y, df)];
    let sum = 0;
    let n = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (ca.isNone(i))
        continue;
      if (cb.isNone(i))
        return `"${y}" is empty in row ${i + 1} of "${t}", where "${x}" has a value`;
      sum += k === 'agree' ? (String(ca.get(i)) === String(cb.get(i)) ? 1 : 0) : Math.abs(ca.get(i) - cb.get(i));
      n++;
    }
    if (n === 0)
      throw new Error(`no row of "${t}" has a value of "${x}"`);
    return sum / n;
  }, [table, a, b, kind]);
}

async function expectStatistic(page: Page, table: string, a: string, b: string, kind: 'agree' | 'mae', model: string,
  metric: string, within: number): Promise<void> {
  await guardIdentity(page);
  const [stored, value] = await Promise.all([modelNumber(model, `metrics.train.${metric}`),
    pairStatistic(page, table, a, b, kind)]);
  if (typeof value === 'string')
    throw new Error(value);
  const what = kind === 'agree' ? 'the share of equal rows' : 'the mean absolute difference';
  expect(Math.abs(value - stored) <= within,
    `${what} of "${a}" and "${b}" (${value}) against the train ${metric} of "${model}" (${stored}), within ${within}`)
    .toBe(true);
}

/* The manual's formula column and its average in the context panel, computed in the page instead. */
export const shareAgrees = Then('the share of rows where {string} and {string} columns of the table {string} agree should ' +
  'equal the Train Accuracy of the Forge model {string} within {float}',
(page: Page, a: string, b: string, table: string, model: string, within: number) =>
  expectStatistic(page, table, a, b, 'agree', model, 'accuracy', within),
{tier: 'api', description: 'the share of the rows with an actual value (the first column) where the prediction (the ' +
  'second) agrees, against metrics.train.accuracy of the model row (grok s); an empty prediction in such a row fails'});

export const maeMatches = Then('the mean absolute difference of {string} and {string} columns of the table {string} should ' +
  'equal the Train MAE of the Forge model {string} within {float}',
(page: Page, a: string, b: string, table: string, model: string, within: number) =>
  expectStatistic(page, table, a, b, 'mae', model, 'mae', within),
{tier: 'api', description: 'over the rows with an actual value (the first column), against metrics.train.mae of the ' +
  'model row (grok s); an empty prediction in such a row fails'});

/* LIB-PROMOTE: the library sets a cell to a text; an empty cell is what a user leaves after deleting the text. */
export const clearCells = When('user clears {string} column in rows {string}', async (page: Page, column: string, rows: string) => {
  const left = await viewers.evaluate(page, ([c, r]: [string, number[]]) => {
    const col = (window as any).__bdd.col(c);
    for (const row of r)
      col.set(row - 1, null);
    return r.filter((row) => !col.isNone(row - 1)).join(', ');
  }, [column, listOf(rows).map(Number)]);
  expect(left, `rows of "${column}" still holding a value`).toBe('');
}, {tier: 'api', description: 'the cells of the current table in the rows (1-based, comma-separated) set empty, read back'});

/* LIB-PROMOTE: the library claims a cell's exact text; a prediction is one of several values. */
export const valueMatches = Then('the value of {string} column in row {int} should match {string}', async (page: Page,
  column: string, row: number, pattern: string) => {
  const text = await viewers.evaluate(page, ([c, r]: [string, number]) => {
    const col = (window as any).__bdd.col(c);
    return col.isNone(r - 1) ? '' : String(col.get(r - 1));
  }, [column, row]);
  expect(text, `"${column}" in row ${row}`).toMatch(new RegExp(pattern));
}, {description: 'rows count from 1; the value as text (empty for a missing value) against the regular expression'});

/** For each of [columns] of the current table, its empty rows (1-based). */
function emptyRowsOf(page: Page, columns: string[]): Promise<number[][]> {
  return viewers.evaluate(page, (names: string[]) => {
    const bdd = (window as any).__bdd;
    const t = bdd.table();
    return names.map((n) => {
      const col = bdd.col(n, t);
      const empty: number[] = [];
      for (let r = 0; r < t.rowCount; r++) {
        if (col.isNone(r))
          empty.push(r + 1);
      }
      return empty;
    });
  }, columns);
}

/* LIB-PROMOTE: the library claims a column complete or not; the user reads the count in the column's statistics. */
export const missingValueCount = Then('{string} column should have {int} missing value(s)', async (page: Page,
  column: string, n: number) => {
  expect((await emptyRowsOf(page, [column]))[0].length, `missing values in "${column}"`).toBe(n);
}, {description: 'the empty cells of the column of the current table; read once, after the end signal of what could ' +
  'change it'});

/* LIB-PROMOTE: a gap count over several columns (the rows the Skip rows option leaves out). */
export const rowsWithGaps = Then('the rows with a missing value in any of {string} columns should be {int}', async (page: Page,
  list: string, n: number) => {
  const gaps = new Set((await emptyRowsOf(page, listOf(list))).flat());
  expect(gaps.size, `rows of the current table with a missing value in any of ${list}`).toBe(n);
}, {description: 'the rows of the current table where at least one of the columns (comma-separated) is empty'});

/* LIB-PROMOTE: a derived column against the gaps of its inputs, row by row, both ways. */
export const emptyExactlyInGaps = Then('the {string} column should be empty exactly in the rows where any of {string} ' +
  'columns is empty', async (page: Page, column: string, list: string) => {
  const [empty, ...inputs] = await emptyRowsOf(page, [column, ...listOf(list)]);
  const gaps = new Set(inputs.flat());
  const emptySet = new Set(empty);
  const emptyWithoutGap = empty.filter((r) => !gaps.has(r));
  const filledWithGap = [...gaps].filter((r) => !emptySet.has(r)).sort((a, b) => a - b);
  expect({emptyWithoutGap: emptyWithoutGap.slice(0, 10), filledWithGap: filledWithGap.slice(0, 10)},
    `"${column}" against the gaps of ${list}: ${emptyWithoutGap.length} rows empty without a gap, ` +
    `${filledWithGap.length} rows filled with a gap (the first 10 of each, 1-based)`)
    .toEqual({emptyWithoutGap: [], filledWithGap: []});
}, {description: 'rows count from 1; fails on an empty cell in a row where every listed column has a value and on a ' +
  'filled one where some listed column is empty, naming the first rows of each'});

/* LIB-PROMOTE: the library's column claims read the current table; this names the open table. Read once, after the
   end signal of what could have added the column. */
export const namedTableNoColumnMatching = Then('the table {string} should have no column matching {string}', async (page: Page,
  table: string, pattern: string) => {
  const names = await viewers.evaluate(page, (t: string) => (window as any).__bdd.tableNamed(t).columns.names() as string[],
    table);
  expect(names.filter((n) => new RegExp(pattern).test(n)), `the columns of "${table}" matching /${pattern}/`).toEqual([]);
}, {description: 'the open table of that name, read once'});

/* LIB-PROMOTE: the missing scenario-end hook. A scenario stops at its first failing step, so a table it opened for
   itself stays open when it fails: the next scenario closes it first, or every later claim about the open tables fails
   with it. */
export const tableClosed = Given('the table {string} is closed', async (page: Page, table: string) => {
  await page.evaluate((t) => {
    for (const df of Array.from(grok.shell.tables as any[]).filter((d: any) => d.name === t))
      grok.shell.closeTable(df);
  }, table);
  await expect.poll(() => page.evaluate((t) => Array.from(grok.shell.tables as any[]).filter((d: any) => d.name === t).length,
    table), {message: `the open tables named "${table}"`}).toBe(0);
}, {tier: 'api', description: 'every open table of that name closed with its views (none: nothing to do), read back'});

/* The console line of the manual (3.8), without the console: the call a script makes. */
export const applyThroughApi = When('user applies the Forge model {string} to the current table through the API', async (page: Page,
  name: string) => {
  await page.evaluate((model) => grok.functions.call('Forge:applyModel',
    {model, table: grok.shell.t, columnNamesMap: {}, showProgress: true}), name);
}, {tier: 'api', description: 'Forge:applyModel {model: the name, table: the current table, columnNamesMap: {}, ' +
  'showProgress: true}, awaited'});

export const applyThroughApiFails = Then('applying the Forge model {string} to the table {string} through the API should ' +
  'fail with {string}', async (page: Page, name: string, table: string, message: string) => {
  const outcome = await viewers.evaluate(page, async ([model, t]: [string, string]) => {
    const df = (window as any).__bdd.tableNamed(t);
    try {
      await grok.functions.call('Forge:applyModel', {model, table: df, columnNamesMap: {}, showProgress: true});
      return 'it applied';
    }
    catch (e) {
      return String((e as {message?: unknown} | null)?.message ?? e);
    }
  }, [name, table]);
  expect(outcome, `applying "${name}" to "${table}"`).toBe(message);
}, {tier: 'api', description: 'Forge:applyModel on the open table of that name is refused with exactly that message'});

/* LIB-PROMOTE: the library claims a tag's value; this claims its absence. */
export const columnLacksTag = Then('the {string} column should not have the tag {string}', async (page: Page, column: string,
  tag: string) => {
  const value = await viewers.evaluate(page, ([c, t]: [string, string]) =>
    (window as any).__bdd.col(c).getTag(t) as string | null, [column, tag]);
  expect(value ?? null, `the ${tag} tag of "${column}"`).toBeNull();
}, {description: 'the column of the current table has no tag of that name'});

/* LIB-PROMOTE: a Dart radio input marks the option chosen only on its own <input type="radio"> (no value attribute:
   every option reads "on"), so "should have value" cannot tell the options apart. */
export const radioChosen = Then('{element} should have {string} chosen', async (page: Page, target: ElementRef, option: string) => {
  const radios = (await locate(page, target)).first().locator('input[type="radio"]');
  await expect.poll(() => radios.evaluateAll((all) => all.filter((r) => (r as HTMLInputElement).checked)
    .map((r) => r.getAttribute('data-value') ?? '').join(', ')), {message: `the option chosen in ${target.phrase}`})
    .toBe(option);
}, {description: 'the checked radio of a Dart radio input, by its option (data-value)'});

const PICKER = '.d4-column-grid:has(.d4-column-selector-backdrop)';

/* LIB-PROMOTE: the library picks a column and reads what the selector took; what the picker offers is the rows of its
   grid (`__name`) that pass its filter, which no library step reads. */
export const columnPickerOffers = Then('the column picker of {element} should offer exactly {string}', async (page: Page,
  target: ElementRef, list: string) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  const selector = await loc.evaluate((e) => e.classList.contains('d4-column-selector')) ? loc :
    loc.locator('.d4-column-selector').first();
  await gestures.openColumnSelector(page, selector);
  const grid = page.locator(PICKER).filter({visible: true}).last().locator('[name="viewer-Grid"]').first();
  try {
    await expect.poll(() => grid.evaluate((el) => {
      const df = DG.Widget.find(el).dataFrame;
      const names = df.col('__name');
      const offered: string[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if (df.filter.get(i))
          offered.push(String(names.get(i)));
      }
      return offered.join(', ');
    }), {message: `the columns the picker of ${target.phrase} offers`}).toBe(listOf(list).join(', '));
  }
  finally {
    // An Escape that reaches the document closes the picker; one inside the dialog would close the dialog.
    await page.evaluate(() => (document.activeElement as HTMLElement | null)?.blur?.());
    await page.keyboard.press('Escape');
    await expect(page.locator(PICKER).filter({visible: true}), 'the column picker after Escape').toHaveCount(0);
  }
}, {tier: 'ui', description: 'opens the column picker of the input, reads the column names it lists (in its order), ' +
  'and closes it with Escape, leaving the selector as it was'});

interface ImputedCheck {gapRows: number; differ: string[]; differCount: number}

/* Records one application of the model (table "forge-test-impute-check", source api): place it after the counts. */
export const imputedPredictions = Then('the {string} column should hold the predictions of the Forge model {string} on ' +
  '{string} imputed with {int} neighbors and {string} distance', async (page: Page, column: string, model: string,
  features: string, neighbors: number, distance: string) => {
  const result = await viewers.evaluate(page, async ([c, m, f, k, d]: [string, string, string[], number, string]):
    Promise<ImputedCheck | string> => {
    const bdd = (window as any).__bdd;
    const t = bdd.table();
    const actual = bdd.col(c, t);
    const cols = f.map((n) => bdd.col(n, t));
    let gapRows = 0;
    for (let r = 0; r < t.rowCount; r++) {
      if (cols.some((col: any) => col.isNone(r)))
        gapRows++;
    }
    const copy = DG.DataFrame.fromColumns(cols.map((col: any) => col.clone()));
    copy.name = 'forge-test-impute-check';
    const imputable = ['int', 'double', 'string', 'datetime', 'qnum'];
    await grok.functions.call('Eda:knnImpute', {table: copy,
      columns: cols.filter((col: any) => col.stats.missingValueCount > 0).map((col: any) => col.name),
      features: cols.filter((col: any) => imputable.includes(col.type)).map((col: any) => col.name),
      neighbors: k, distance: d});
    const before = copy.columns.length;
    await grok.functions.call('Forge:applyModel', {model: m, table: copy, columnNamesMap: {}, showProgress: false});
    if (copy.columns.length !== before + 1)
      return 'Forge:applyModel added no column to the imputed copy';
    const expected = copy.columns.byIndex(before);
    const differ: string[] = [];
    for (let r = 0; r < t.rowCount; r++) {
      const a = actual.isNone(r) ? null : actual.get(r);
      const e = expected.isNone(r) ? null : expected.get(r);
      const same = a === null || e === null ? a === e : typeof a === 'number' && typeof e === 'number' ?
        Math.abs(a - e) <= 1e-9 * Math.max(1, Math.abs(e)) : String(a) === String(e);
      if (!same)
        differ.push(`row ${r + 1}: ${String(a)} against ${String(e)}`);
    }
    return {gapRows, differ: differ.slice(0, 10), differCount: differ.length};
  }, [column, model, listOf(features), neighbors, distance]);
  if (typeof result === 'string')
    throw new Error(result);
  expect(result.gapRows > 0, `rows with a gap in ${features}`).toBe(true);
  expect(result.differ, `"${column}" against the predictions on the imputed copy (${result.gapRows} rows with a gap; ` +
    `${result.differCount} rows differ, the first 10)`).toEqual([]);
}, {tier: 'api', description: 'Impute that fills nothing still predicts every row (the method reads an empty float as a ' +
  'tiny number): the features of the current table cloned, imputed with Eda:knnImpute as Forge does (the gapped ' +
  'columns filled, every imputable feature used to find the neighbours), predicted with Forge:applyModel (one more ' +
  'application of the model); every row of the column equals that prediction, and some row has a gap'});

/* LIB-PROMOTE: "the context panel should show" finds the name in the panel's text, and the panel of the column
   "Species (predicted 2)" holds "Species" too: a claim about the new panel waits for its title. */
export const panelTitled = Then('the context panel title should be {string}', async (page: Page, name: string) => {
  await expect.poll(async () => {
    const [object, title] = await Promise.all([page.evaluate(() => String(grok.shell.o?.name ?? grok.shell.o)),
      contextPanelTitle(page)]);
    return `${object} | ${title}`;
  }, {message: 'the current object and the title of the context panel'}).toBe(`${name} | ${name}`);
}, {description: 'the current object\'s name and the context panel\'s title are both exactly the text'});
