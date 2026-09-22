/* The columns of the current table, as facts: what the detectors set (semantic type, units and
   the other tags), the storage type, and the cells — every filled value against a pattern, a
   range, a length; one row's value; the maximum's row; distinct counts — and the current row.
   A step that changes the current row baselines every viewer first (see data.ts). */
import {Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Then, When} from '../../src/registry.js';
import {changeAll} from '../../src/runtime/viewers.js';

declare const grok: any;

interface ColumnFacts {
  type: string;
  semType: string;
  tags: Record<string, string>;
  rows: number;
  missing: number;
  distinct: number;
  values: string[];
  numbers: (number | null)[];
}

/** Everything a column claim reads, in one evaluate; values are strings, numbers those of a
 * numeric column (null where missing). */
function columnFacts(page: Page, column: string): Promise<ColumnFacts> {
  return page.evaluate((c) => {
    const df = grok.shell.t;
    if (!df)
      throw new Error('no table is open');
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const tags: Record<string, string> = {};
    for (const [k, v] of Array.from(col.tags as Iterable<[string, string]>))
      tags[String(k)] = String(v);
    const values: string[] = [];
    const numbers: (number | null)[] = [];
    let missing = 0;
    const numeric = col.matches('numerical');
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i)) {
        missing++;
        values.push('');
        numbers.push(null);
        continue;
      }
      const v = col.get(i);
      values.push(String(v));
      numbers.push(numeric ? Number(v) : null);
    }
    return {type: String(col.type), semType: String(col.semType ?? ''), tags, rows: df.rowCount, missing, distinct: new Set(values.filter((v) => v !== '')).size, values, numbers};
  }, column);
}

const filled = (f: ColumnFacts): number[] => f.values.map((_, i) => i).filter((i) => f.values[i] !== '');

export const columnSemType = Then('{string} column should have semantic type {string}', async (page: Page, column: string, semType: string) => {
  await expect.poll(async () => (await columnFacts(page, column)).semType, {message: `semantic type of "${column}"`}).toBe(semType);
}, {description: 'what the detectors set (Macromolecule, Molecule, Monomer); polled, since detection runs after the column appears'});

export const columnUnits = Then('{string} column should have units {string}', async (page: Page, column: string, units: string) => {
  await expect.poll(async () => (await columnFacts(page, column)).tags['units'] ?? '', {message: `units of "${column}"`}).toBe(units);
}, {description: 'the `units` tag — a sequence column\'s notation (fasta, separator, helm), a molecule column\'s molblock'});

export const columnTag = Then('{string} column should have tag {string} equal to {string}', async (page: Page, column: string, tag: string, value: string) => {
  await expect.poll(async () => (await columnFacts(page, column)).tags[tag] ?? '', {message: `tag "${tag}" of "${column}"`}).toBe(value);
});

export const columnTagLists = Then('the {string} tag of {string} column should list at least {int} values', async (page: Page, tag: string, column: string, count: number) => {
  const values = ((await columnFacts(page, column)).tags[tag] ?? '').split(',').map((s) => s.trim()).filter((s) => s !== '');
  expect(values.length, `values in the "${tag}" tag of "${column}": ${values.slice(0, 5).join(', ')}${values.length > 5 ? ', …' : ''}`).toBeGreaterThanOrEqual(count);
}, {description: 'a comma-separated tag (the .positionNames a numbering run writes on the aligned column)'});

export const columnType = Then('{string} column should have type {string}', async (page: Page, column: string, type: string) => {
  expect((await columnFacts(page, column)).type, `type of "${column}"`).toBe(type);
}, {description: 'the storage type: string, int, double, ...'});

export const columnComplete = Then('{string} column should have no missing values', async (page: Page, column: string) => {
  expect((await columnFacts(page, column)).missing, `missing values in "${column}"`).toBe(0);
}, {description: 'of the current table'});

export const columnIncomplete = Then('{string} column should have missing values', async (page: Page, column: string) => {
  expect((await columnFacts(page, column)).missing, `missing values in "${column}"`).toBeGreaterThan(0);
});

export const everyValueMatches = Then('every value of {string} column should match {string}', async (page: Page, column: string, pattern: string) => {
  const f = await columnFacts(page, column);
  const re = new RegExp(pattern);
  const bad = filled(f).filter((i) => !re.test(f.values[i]));
  expect(bad.map((i) => `row ${i + 1}: ${f.values[i].slice(0, 60)}`), `values of "${column}" not matching /${pattern}/ (missing values skipped)`).toEqual([]);
  expect(filled(f).length, `filled values of "${column}"`).toBeGreaterThan(0);
}, {description: 'a regular expression over every filled cell; a column with no filled cell fails'});

export const fewerDistinctThanRows = Then('{string} column should have fewer distinct values than the table has rows', async (page: Page, column: string) => {
  const f = await columnFacts(page, column);
  const distinct = new Set(f.values.filter((_, i) => filled(f).includes(i))).size;
  expect(distinct, `distinct values of "${column}" against the ${f.values.length} rows`).toBeLessThan(f.values.length);
}, {description: 'a column that groups the rows rather than naming each of them'});

export const someValueMatches = Then('some value of {string} column should match {string}', async (page: Page, column: string, pattern: string) => {
  const f = await columnFacts(page, column);
  const re = new RegExp(pattern);
  expect(filled(f).filter((i) => re.test(f.values[i])).length, `values of "${column}" matching /${pattern}/`).toBeGreaterThan(0);
}, {description: 'a regular expression that at least one filled cell matches'});

export const someValueDiffers = Then('some value of {string} column should differ from {string} column in the same row', async (page: Page, x: string, y: string) => {
  const [a, b] = [await columnFacts(page, x), await columnFacts(page, y)];
  expect(a.values.filter((v, i) => v !== b.values[i]).length, `rows where "${x}" and "${y}" hold different text`).toBeGreaterThan(0);
}, {description: 'the cells\' text, row by row'});

export const everyValueContains = Then('every value of {string} column should contain {string}', async (page: Page, column: string, text: string) => {
  const f = await columnFacts(page, column);
  const bad = filled(f).filter((i) => !f.values[i].includes(text));
  expect(bad.map((i) => `row ${i + 1}: ${f.values[i].slice(0, 60)}`), `values of "${column}" without "${text}" (missing values skipped)`).toEqual([]);
  expect(filled(f).length, `filled values of "${column}"`).toBeGreaterThan(0);
});

export const everyValueBetween = Then('every value of {string} column should lie between {float} and {float}', async (page: Page, column: string, lo: number, hi: number) => {
  const f = await columnFacts(page, column);
  const nums = f.numbers.filter((n): n is number => n !== null);
  expect(nums.length, `numbers in "${column}"`).toBeGreaterThan(0);
  expect(nums.filter((n) => n < lo || n > hi), `values of "${column}" outside ${lo}..${hi}`).toEqual([]);
});

export const everyValueSameLength = Then('every value of {string} column should have the same length', async (page: Page, column: string) => {
  const f = await columnFacts(page, column);
  expect([...new Set(filled(f).map((i) => f.values[i].length))], `distinct lengths in "${column}"`).toHaveLength(1);
}, {description: 'the filled cells — an aligned column'});

export const sameLengthPerGroup = Then('every value of {string} column should have the same length within each {string} value', async (page: Page, column: string, group: string) => {
  const f = await columnFacts(page, column);
  const g = await columnFacts(page, group);
  const lengths = new Map<string, Set<number>>();
  for (const i of filled(f))
    lengths.set(g.values[i], (lengths.get(g.values[i]) ?? new Set()).add(f.values[i].length));
  const uneven = [...lengths].filter(([, ls]) => ls.size > 1).map(([k, ls]) => `${group} ${k}: ${[...ls].join(', ')}`);
  expect(uneven, `groups of "${group}" whose "${column}" values differ in length`).toEqual([]);
}, {description: 'an alignment run per cluster: one width per cluster, not one width overall'});

export const distinctLengths = Then('the values of {string} column should have at least {int} distinct lengths', async (page: Page, column: string, count: number) => {
  const f = await columnFacts(page, column);
  const lengths = [...new Set(filled(f).map((i) => f.values[i].length))].sort((a, b) => a - b);
  expect(lengths.length, `distinct lengths of the filled values of "${column}": ${lengths.join(', ') || 'none'}`).toBeGreaterThanOrEqual(count);
}, {description: 'the filled cells — an alignment run per cluster pads each cluster to its own width, one global alignment to one'});

export const someValueContains = Then('some value of {string} column should contain {string}', async (page: Page, column: string, text: string) => {
  const f = await columnFacts(page, column);
  expect(filled(f).some((i) => f.values[i].includes(text)), `a value of "${column}" containing "${text}"`).toBe(true);
});

export const columnsEqual = Then('{string} column should hold the same values as {string} column', async (page: Page, a: string, b: string) => {
  const fa = await columnFacts(page, a);
  const fb = await columnFacts(page, b);
  expect(fa.rows, 'rows of the current table').toBeGreaterThan(0);
  const bad = fa.values.map((x, i) => [i, x, fb.values[i]] as const).filter(([, x, y]) => x !== y)
    .map(([i, x, y]) => `row ${i + 1}: ${x.slice(0, 30)} ≠ ${y.slice(0, 30)}`);
  expect(bad, `rows where "${a}" and "${b}" differ`).toEqual([]);
}, {description: 'row by row, as text; a missing value equals only a missing value'});

export const displayedInRow = Then('the {string} cell of row {int} should be displayed as {string}',
  async (page: Page, column: string, row: number, text: string) => {
    const shown = await page.evaluate(([c, r]) => {
      const t = grok.shell.t;
      const grid = grok.shell.tv?.grid;
      const i = (r as number) - 1;
      if (i < 0 || i >= t.rowCount)
        throw new Error(`row ${r} is outside the table's ${t.rowCount} rows`);
      return grid != null && grid.col(c) != null
        ? String(grid.cell(c as string, i).cell.valueString)
        : String(t.col(c as string).getString(i));
    }, [column, row] as [string, number]);
    expect(shown, `the grid's text for "${column}" in row ${row}`).toBe(text);
  }, {description: 'what the grid draws in the cell — the column\'s format applied, unlike "the value of … column in row …", which reads the raw value'});

export const valueInRow = Then('the value of {string} column in row {int} should be {string}', async (page: Page, column: string, row: number, value: string) => {
  const f = await columnFacts(page, column);
  if (row < 1 || row > f.rows)
    throw new Error(`row ${row} is outside the table's ${f.rows} rows`);
  expect(f.values[row - 1], `"${column}" in row ${row}`).toBe(value);
}, {description: 'rows count from 1; a number reads as the platform prints it'});

export const maxInRow = Then('{string} column should have its maximum in row {int}', async (page: Page, column: string, row: number) => {
  const f = await columnFacts(page, column);
  const max = Math.max(...f.numbers.filter((n): n is number => n !== null));
  expect(f.numbers[row - 1], `"${column}" in row ${row} (its maximum is ${max})`).toBe(max);
});

export const distinctValues = Then('{string} column should have at least {int} distinct values', async (page: Page, column: string, count: number) => {
  expect((await columnFacts(page, column)).distinct, `distinct values of "${column}"`).toBeGreaterThanOrEqual(count);
});

const columnNames = (page: Page): Promise<string[]> => page.evaluate(() => grok.shell.t?.columns.names() ?? []);

export const hasColumn = Then('the table should have a column {string}', async (page: Page, column: string) => {
  await expect.poll(() => columnNames(page), {message: 'columns of the current table', timeout: pollMs(60000)}).toContain(column);
}, {description: 'the current table, by exact name; a column a computation produces arrives when the computation ends, so the claim carries that budget'});

export const hasNoColumn = Then('the table should not have a column {string}', async (page: Page, column: string) => {
  expect(await columnNames(page), 'columns of the current table').not.toContain(column);
});

export const columnCount = Then('the table should have {int} column(s)', async (page: Page, count: number) => {
  await expect.poll(async () => (await columnNames(page)).length, {message: 'columns of the current table'}).toBe(count);
});

// --- the current row ---------------------------------------------------------------------------------

function makeCurrent(page: Page, row: number | 'last'): Promise<void> {
  return changeAll(page, (r) => {
    const df = grok.shell.t;
    const idx = r === 'last' ? df.rowCount - 1 : r - 1;
    if (idx < 0 || idx >= df.rowCount)
      throw new Error(`row ${r} is outside the table's ${df.rowCount} rows`);
    df.currentRowIdx = idx;
  }, row);
}

export const makeRowCurrent = When('user makes row {int} current', (page: Page, row: number) => makeCurrent(page, row),
  {tier: 'api', description: 'rows count from 1; the UI path is a click on the row in the grid'});

export const makeLastRowCurrent = When('user makes the last row current', (page: Page) => makeCurrent(page, 'last'), {tier: 'api'});

export const currentRowIs = Then('row {int} should be current', async (page: Page, row: number) => {
  await expect.poll(() => page.evaluate(() => grok.shell.t.currentRowIdx + 1), {message: 'the current row'}).toBe(row);
});

export const mouseOverRowIs = Then('row {int} should be under the mouse', async (page: Page, row: number) => {
  await expect.poll(() => page.evaluate(() => grok.shell.t.mouseOverRowIdx + 1), {message: 'the row under the mouse'}).toBe(row);
}, {description: 'rows count from 1 — the table\'s mouse-over row, which every viewer of the table highlights'});

export const joinedValues = Then('every value of {string} column should be {string} and {string} of the same row joined by {string}', async (page: Page, column: string, a: string, b: string, sep: string) => {
  const bad: string[] = await page.evaluate(([c, x, y, s]) => {
    const df = grok.shell.t;
    for (const n of [c, x, y])
      if (!df.col(n))
        throw new Error(`no "${n}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const out: string[] = [];
    for (let i = 0; i < df.rowCount; i++) {
      const want = `${df.col(x).get(i) ?? ''}${s}${df.col(y).get(i) ?? ''}`;
      if (String(df.col(c).get(i) ?? '') !== want)
        out.push(`row ${i + 1}`);
    }
    return out;
  }, [column, a, b, sep] as [string, string, string, string]);
  expect(bad, `rows of "${column}" that are not "${a}${sep}${b}"`).toEqual([]);
}, {description: 'a pairing column: each cell is the two source cells of its row with the separator between'});

export const setColumnSemType = When('user sets the semantic type of {string} column to {string}', async (page: Page, column: string, semType: string) => {
  await page.evaluate(([c, s]) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}; it has: ${grok.shell.t.columns.names().join(', ')}`);
    col.semType = s;
  }, [column, semType] as [string, string]);
}, {tier: 'api', description: 'the column\'s semantic type written directly — a type the detectors would not give the column, for a claim that something keeps it'});
