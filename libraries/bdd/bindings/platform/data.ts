/* The current table's rows and columns through the JS API: selection, filter, the current row,
   deleting rows, calculated columns, column colors, and the tables the workspace holds. A step
   that changes rows or colors takes every open viewer's baseline first, so the viewer checks that
   follow ("should have repainted", "a narrower value range than before") read against the state
   before the change, and returns once every viewer has drawn the change, so the next step's
   baseline is the state after it. */
import {expect, Page} from '@playwright/test';
import {Then, When} from '../../src/registry.js';
import {baselineAll, settleAll} from '../../src/runtime/viewers.js';

declare const grok: any;
declare const DG: any;

interface Match {
  matching: number;
  selected: number;
  selectedTotal: number;
}

/** Rows where the column reads one of the values, and how many of them (and of all rows) are
 * selected. */
function matchSelection(page: Page, column: string, values: string[]): Promise<Match> {
  return page.evaluate(([c, vs]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let matching = 0;
    let selected = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (!vs.includes(String(col.get(i) ?? '')))
        continue;
      matching++;
      if (df.selection.get(i))
        selected++;
    }
    return {matching, selected, selectedTotal: df.selection.trueCount};
  }, [column, values] as [string, string[]]);
}

/** The match, or a failure when no row reads the value at all: a category that does not exist
 * (a typo) must not pass a check about its rows. */
async function matchSome(page: Page, column: string, values: string[]): Promise<Match> {
  const m = await matchSelection(page, column, values);
  if (m.matching === 0)
    throw new Error(`no row of ${await tableName(page)} has ${column} ${values.length > 1 ? `in ${values.join(', ')}` : `= ${values[0]}`}`);
  return m;
}

const tableName = (page: Page): Promise<string> => page.evaluate(() => String(grok.shell.t.name));
const selectedCount = (page: Page): Promise<number> => page.evaluate(() => grok.shell.t.selection.trueCount as number);
const filteredCount = (page: Page): Promise<number> => page.evaluate(() => grok.shell.t.filter.trueCount as number);
const list = (s: string): string[] => s.split(/\s*,\s*/).filter((x) => x.length > 0);

/** A change to the table every viewer answers: baselines first, the change, then every viewer
 * has drawn it. */
async function changeTable(page: Page, body: (arg: any) => void, arg: unknown): Promise<void> {
  await baselineAll(page);
  await page.evaluate(body, arg);
  await settleAll(page);
}

// --- selection ---------------------------------------------------------------------------------------

export const clearSelection = When('user clears the row selection', (page: Page) =>
  changeTable(page, () => { grok.shell.t.selection.setAll(false); }, null), {tier: 'api'});

/** The rows where the column reads one of the values become the selection, nothing else; a
 * category no row has fails (a typo must not pass as an empty selection). */
function selectWhere(page: Page, column: string, values: string[]): Promise<void> {
  return changeTable(page, ([c, vs]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let hits = 0;
    df.selection.init((i: number) => {
      const hit = vs.includes(String(col.get(i) ?? ''));
      if (hit)
        hits++;
      return hit;
    });
    if (hits === 0)
      throw new Error(`no row of ${df.name} has ${c} ${vs.length > 1 ? `in ${vs.join(', ')}` : `= ${vs[0]}`}`);
  }, [column, values] as [string, string[]]);
}

export const selectWhereIs = When('user selects rows where {string} is {string}', (page: Page, column: string, value: string) =>
  selectWhere(page, column, [value]), {tier: 'api', description: 'the table\'s selection bitset — every row of the category and nothing else; the UI path is a click on the category in a viewer'});

export const selectWhereOneOf = When('user selects rows where {string} is one of {string}', (page: Page, column: string, values: string) =>
  selectWhere(page, column, list(values)), {tier: 'api', description: 'comma-separated categories'});

export const noneSelected = Then('no rows should be selected', (page: Page) =>
  expect.poll(() => selectedCount(page), {message: 'rows selected'}).toBe(0));

export const someSelected = Then('some rows should be selected', (page: Page) =>
  expect.poll(() => selectedCount(page), {message: 'no row is selected'}).toBeGreaterThan(0));

export const allOfSelected = Then('all rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => {
    const m = await matchSome(page, column, [value]);
    return `${m.selected} of ${m.matching}`;
  }, {message: `rows where ${column} is ${value} selected`}).toMatch(/^(\d+) of \1$/);
}, {description: 'every row of the category, whatever else is selected'});

export const onlyOfSelected = Then('only rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => {
    const m = await matchSome(page, column, [value]);
    return m.selected === m.matching && m.selectedTotal === m.matching ? 'exactly' : `${m.selected} of ${m.matching}, ${m.selectedTotal} in all`;
  }, {message: `rows where ${column} is ${value} selected`}).toBe('exactly');
}, {description: 'every row of the category and nothing else'});

export const onlyOfAnySelected = Then('only rows where {string} is one of {string} should be selected', async (page: Page, column: string, values: string) => {
  await expect.poll(async () => {
    const m = await matchSome(page, column, list(values));
    return m.selected === m.matching && m.selectedTotal === m.matching ? 'exactly' : `${m.selected} of ${m.matching}, ${m.selectedTotal} in all`;
  }, {message: `rows where ${column} is one of ${values} selected`}).toBe('exactly');
}, {description: 'every row of these categories (comma-separated) and nothing else — a union built with Control clicks'});

export const noneOfSelected = Then('no rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await matchSome(page, column, [value])).selected, {message: `rows where ${column} is ${value} selected`}).toBe(0);
});

export const someOfSelected = Then('some rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await matchSome(page, column, [value])).selected, {message: `rows where ${column} is ${value} selected`}).toBeGreaterThan(0);
});

export const hasCurrentRow = Then('the table should have a current row', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.currentRowIdx as number), {message: 'the current row index'}).toBeGreaterThanOrEqual(0));

export const currentRowValue = Then('{string} of the current row should be {string}', async (page: Page, column: string, value: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    return df.currentRowIdx < 0 ? '(no current row)' : String(col.get(df.currentRowIdx) ?? '');
  }, column), {message: `"${column}" of the current row`}).toBe(value);
}, {description: 'the column\'s value in the current row, as text'});

// --- filter ------------------------------------------------------------------------------------------

export const filterBetween = When('user filters rows where {string} is between {float} and {float}', (page: Page, column: string, lo: number, hi: number) =>
  changeTable(page, ([c, a, b]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.filter.init((i: number) => { const x = col.get(i); return x != null && x >= a && x <= b; });
  }, [column, lo, hi] as [string, number, number]),
{tier: 'api', description: 'the table\'s filter bitset, as a filter viewer would set it; the filter panel\'s own handles are not driven'});

export const filterTo = When('user filters rows where {string} is {string}', (page: Page, column: string, value: string) =>
  changeTable(page, ([c, v]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.filter.init((i: number) => String(col.get(i) ?? '') === v);
    if (df.filter.trueCount === 0)
      throw new Error(`no row of ${df.name} has ${c} = ${v}`);
  }, [column, value] as [string, string]), {tier: 'api', description: 'keeps the category\'s rows only — the table\'s filter bitset'});

export const filterOut = When('user filters out rows where {string} is {string}', (page: Page, column: string, value: string) =>
  changeTable(page, ([c, v]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.filter.init((i: number) => String(col.get(i) ?? '') !== v);
  }, [column, value] as [string, string]), {tier: 'api', description: 'the table\'s filter bitset'});

export const resetFilter = When('user resets the filter', (page: Page) =>
  changeTable(page, () => { grok.shell.t.filter.setAll(true); }, null), {tier: 'api'});

export const filterPasses = Then('{int} row(s) should pass the filter', (page: Page, count: number) =>
  expect.poll(() => filteredCount(page), {message: 'rows passing the filter'}).toBe(count));

export const filterPassesFewer = Then('fewer than {int} rows should pass the filter', (page: Page, count: number) =>
  expect.poll(() => filteredCount(page), {message: 'rows passing the filter'}).toBeLessThan(count));

export const filterPassesAll = Then('all rows should pass the filter', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.filter.trueCount === grok.shell.t.rowCount), {message: 'every row passes the filter'}).toBe(true));

export const filterIsExactly = Then('the filter should pass exactly the rows where {string} is between {float} and {float}',
  async (page: Page, column: string, lo: number, hi: number) => {
    const off = await page.evaluate(([c, a, b]) => {
      const df = grok.shell.t;
      const col = df.col(c);
      let wrong = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const x = col.get(i);
        if (df.filter.get(i) !== (x != null && x >= a && x <= b))
          wrong++;
      }
      return wrong;
    }, [column, lo, hi] as [string, number, number]);
    expect(off, `rows whose filter bit disagrees with ${column} in [${lo}, ${hi}]`).toBe(0);
  }, {description: 'the table filter bit by bit — a viewer\'s own filter must leave it alone'});

export const filterIsExactlyCategory = Then('the filter should pass exactly the rows where {string} is {string}', async (page: Page, column: string, value: string) => {
  const off = await page.evaluate(([c, v]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let wrong = 0;
    let matching = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const hit = String(col.get(i) ?? '') === v;
      if (hit)
        matching++;
      if (df.filter.get(i) !== hit)
        wrong++;
    }
    if (matching === 0)
      throw new Error(`no row of ${df.name} has ${c} = ${v}`);
    return wrong;
  }, [column, value] as [string, string]);
  expect(off, `rows whose filter bit disagrees with ${column} = ${value}`).toBe(0);
}, {description: 'the table filter bit by bit: the category\'s rows pass, no other row does'});

// --- rows --------------------------------------------------------------------------------------------

export const deleteSelected = When('user deletes the selected rows', (page: Page) =>
  changeTable(page, () => {
    const df = grok.shell.t;
    const set = new Set<number>(Array.from(df.selection.getSelectedIndexes() as Iterable<number>));
    df.rows.removeWhereIdx((i: number) => set.has(i));
  }, null), {tier: 'api', description: 'df.rows.removeWhereIdx — the UI path is the grid\'s Delete Rows command'});

export const rowCount = Then('the table should have {int} row(s)', (page: Page, count: number) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.rowCount as number), {message: 'rows in the table'}).toBe(count));

export const noRowsWhere = Then('the table should have no rows where {string} is {string}', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await matchSelection(page, column, [value])).matching, {message: `rows where ${column} is ${value}`}).toBe(0);
});

// --- columns -----------------------------------------------------------------------------------------

export const addCalculated = When('user adds a calculated column {string} with formula {string}', async (page: Page, name: string, formula: string) => {
  await page.evaluate(async ([n, f]) => { await grok.shell.t.columns.addNewCalculated(n, f); }, [name, formula] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.t.columns.names().includes(n), name), {message: `"${name}" in the table's columns`}).toBe(true);
}, {tier: 'api', description: 'a formula in the platform\'s syntax: ${HEIGHT} * 2'});

export const removeColumn = When('user removes {string} column', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const df = grok.shell.t;
    if (!df.columns.names().includes(n))
      throw new Error(`no "${n}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.columns.remove(n);
  }, name);
}, {tier: 'api'});

function color(page: Page, column: string, apply: string, arg: unknown): Promise<void> {
  return changeTable(page, ([c, how, a]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const colors = col.meta.colors;
    if (how === 'linear')
      colors.setLinear((a as any).scheme, (a as any).range);
    else if (how === 'conditional')
      colors.setConditional(a);
    else if (how === 'categorical')
      colors.setCategorical(a);
    else
      colors.setDisabled();
  }, [column, apply, arg] as [string, string, unknown]);
}

export const colorLinear = When('user colors {string} column linearly from {string} to {string}', (page: Page, column: string, from: string, to: string) =>
  color(page, column, 'linear', {scheme: [DG_COLOR(from), DG_COLOR(to)]}),
{tier: 'api', description: 'the column\'s color coding, as the grid\'s Color Coding menu sets it; #rrggbb colors'});

export const colorLinearOver = When('user colors {string} column linearly from {string} to {string} over {float} to {float}',
  (page: Page, column: string, from: string, to: string, min: number, max: number) =>
    color(page, column, 'linear', {scheme: [DG_COLOR(from), DG_COLOR(to)], range: {min, max}}), {tier: 'api'});

export const colorConditional = When('user colors {string} column conditionally:', (page: Page, column: string, table: string[][]) =>
  color(page, column, 'conditional', Object.fromEntries(table)), {tier: 'api', description: '| range | color | rows: 50-90 | #00ff00'});

export const colorCategorical = When('user colors {string} column categorically:', (page: Page, column: string, table: string[][]) =>
  color(page, column, 'categorical', Object.fromEntries(table)), {tier: 'api', description: '| category | color | rows'});

export const colorOff = When('user removes the coloring of {string} column', (page: Page, column: string) => color(page, column, 'off', null), {tier: 'api'});

/** The column's color-coding type tag: '' when none. */
function colorCodingType(page: Page, column: string): Promise<string> {
  return page.evaluate((c) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const type = col.getTag('.color-coding-type');
    return type == null || type === 'Off' ? '' : String(type);
  }, column);
}

export const noColorCoding = Then('{string} column should have no color coding', async (page: Page, column: string) => {
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe('');
}, {description: 'the column\'s color-coding tag is unset or Off'});

export const colorCodedCategorically = Then('{string} column should be color-coded categorically', async (page: Page, column: string) => {
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe('Categorical');
});

/** `#rrggbb` → the ARGB int the color API takes. */
function DG_COLOR(hex: string): number {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex);
  if (!m)
    throw new Error(`"${hex}" is not a #rrggbb color`);
  return (0xFF000000 | parseInt(m[1], 16)) >>> 0;
}

// --- tables ------------------------------------------------------------------------------------------

function tableInfo(page: Page, name: string): Promise<{columns: string[]; rows: number} | string[]> {
  return page.evaluate((n) => {
    const tables: any[] = grok.shell.tables;
    const t = tables.find((x) => x.name === n);
    return t ? {columns: t.columns.names(), rows: t.rowCount} : tables.map((x) => x.name);
  }, name);
}

export const tableOpen = Then('table {string} should be open', async (page: Page, name: string) => {
  await expect.poll(async () => {
    const info = await tableInfo(page, name);
    return Array.isArray(info) ? `not open; open: ${info.join(' | ')}` : 'open';
  }, {message: `table "${name}"`}).toBe('open');
}, {description: 'in grok.shell.tables, by its exact name'});

export const tableColumns = Then('table {string} should have columns {string}', async (page: Page, name: string, list: string) => {
  const info = await tableInfo(page, name);
  if (Array.isArray(info))
    throw new Error(`table "${name}" is not open; open: ${info.join(' | ')}`);
  expect(info.columns, `columns of "${name}"`).toEqual(list.split(/\s*,\s*/));
}, {description: 'exactly these, in this order (comma-separated)'});

export const tableRows = Then('table {string} should have {int} row(s)', async (page: Page, name: string, count: number) => {
  const info = await tableInfo(page, name);
  if (Array.isArray(info))
    throw new Error(`table "${name}" is not open; open: ${info.join(' | ')}`);
  expect(info.rows, `rows of "${name}"`).toBe(count);
});

function missingCount(page: Page, name: string, column: string): Promise<number> {
  return page.evaluate(([n, c]) => {
    const t = grok.shell.tables.find((x: any) => x.name === n);
    if (!t)
      throw new Error(`table "${n}" is not open; open: ${grok.shell.tables.map((x: any) => x.name).join(' | ')}`);
    const col = t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in "${n}"; it has: ${t.columns.names().join(', ')}`);
    let missing = 0;
    for (let i = 0; i < t.rowCount; i++) {
      if (col.isNone(i))
        missing++;
    }
    return missing;
  }, [name, column] as [string, string]);
}

export const tableColumnComplete = Then('table {string} should have no missing values in {string} column', async (page: Page, name: string, column: string) => {
  expect(await missingCount(page, name, column), `missing values in "${column}" of "${name}"`).toBe(0);
}, {description: 'a computed column is filled in, not blank — a result table whose fit failed has the right shape and no numbers'});

export const tableColumnIncomplete = Then('table {string} should have missing values in {string} column', async (page: Page, name: string, column: string) => {
  expect(await missingCount(page, name, column), `missing values in "${column}" of "${name}"`).toBeGreaterThan(0);
}, {description: 'at least one blank — the precondition of a scenario about empty categories'});
