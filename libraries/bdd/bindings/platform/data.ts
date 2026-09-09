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

export const selectWhereBetween = When('user selects rows where {string} is between {float} and {float}', (page: Page, column: string, min: number, max: number) =>
  changeTable(page, ([c, lo, hi]) => {
    const df = grok.shell.t;
    const col = df.col(c as string);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let hits = 0;
    df.selection.init((i: number) => {
      const v = col.get(i);
      const hit = v != null && !isNaN(v) && v >= (lo as number) && v <= (hi as number);
      if (hit)
        hits++;
      return hit;
    });
    if (hits === 0)
      throw new Error(`no row of ${df.name} has ${c} between ${lo} and ${hi}`);
  }, [column, min, max] as [string, number, number]),
{tier: 'api', description: 'a numeric range, both ends included — the selection becomes exactly those rows; a range no row matches fails'});

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

/** How many rows of a category the filter keeps — throws when the category is not in the column, so
 * a typo fails instead of passing as "none". */
async function matchFilter(page: Page, column: string, value: string): Promise<{matching: number; passing: number}> {
  const m = await page.evaluate(([c, v]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let matching = 0;
    let passing = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (String(col.get(i) ?? '') !== v)
        continue;
      matching++;
      if (df.filter.get(i))
        passing++;
    }
    return {matching, passing};
  }, [column, value] as [string, string]);
  if (m.matching === 0)
    throw new Error(`no row of ${await tableName(page)} has ${column} = ${value}`);
  return m;
}

export const noneOfFiltered = Then('no rows where {string} is {string} should pass the filter', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await matchFilter(page, column, value)).passing,
    {message: `rows where ${column} is ${value} passing the filter`}).toBe(0);
}, {description: 'the category is filtered out — what the filter actually keeps, not what a card says it keeps'});

export const allOfFiltered = Then('all rows where {string} is {string} should pass the filter', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => {
    const m = await matchFilter(page, column, value);
    return m.passing === m.matching ? 'all' : `${m.passing} of ${m.matching}`;
  }, {message: `rows where ${column} is ${value} passing the filter`}).toBe('all');
});

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

export const filterNotNull = When('user filters rows where {string} is not null', (page: Page, column: string) =>
  changeTable(page, ([c]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.filter.init((i: number) => { return !col.isNone(i) });
  }, [column] as [string]),
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

export const filterToAnyOf = When('user filters rows where {string} is one of {string}', (page: Page, column: string, values: string) =>
  changeTable(page, ([c, vs]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    df.filter.init((i: number) => vs.includes(String(col.get(i) ?? '')));
    if (df.filter.trueCount === 0)
      throw new Error(`no row of ${df.name} has ${c} in ${vs.join(', ')}`);
  }, [column, list(values)] as [string, string[]]),
{tier: 'api', description: 'keeps the rows of these categories (comma-separated) and no other'});

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

export const setCell = When('user sets {string} column in row {int} to {string}', (page: Page, column: string, row: number, value: string) =>
  changeTable(page, ([c, r, v]) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}; it has: ${grok.shell.t.columns.names().join(', ')}`);
    const i = (r as number) - 1;
    if (i < 0 || i >= grok.shell.t.rowCount)
      throw new Error(`row ${r} is outside ${grok.shell.t.name} (${grok.shell.t.rowCount} rows)`);
    const text = v as string;
    if (text === '')
      col.set(i, null);
    else if (col.type === 'double' || col.type === 'int' || col.type === 'float' || col.type === 'bigint')
      col.set(i, text === 'NaN' ? Number.NaN : text === 'Infinity' ? Number.POSITIVE_INFINITY
        : text === '-Infinity' ? Number.NEGATIVE_INFINITY : Number(text));
    else
      col.set(i, text);
  }, [column, row, value] as [string, number, string]),
{tier: 'api', description: 'one cell, 1-based row; on a numeric column "NaN", "Infinity" and "-Infinity" write those values and "" writes a blank — what a viewer must survive'});

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

export const renameColumn = When('user renames {string} column to {string}', async (page: Page, from: string, to: string) => {
  await page.evaluate(([a, b]) => {
    const df = grok.shell.t;
    const col = df.col(a);
    if (!col)
      throw new Error(`no "${a}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    col.name = b;
  }, [from, to] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.t.columns.names().includes(n), to), {message: `"${to}" in the table's columns`}).toBe(true);
}, {tier: 'api', description: 'the column\'s name through the API — the UI path is the header\'s Column Properties dialog'});

export const currentColumnIs = Then('the current column should be {string}', (page: Page, name: string) =>
  expect.poll(() => page.evaluate(() => String(grok.shell.t.currentCol?.name ?? '')), {message: 'the current column'}).toBe(name),
{description: 'the table\'s current column by name ("" when none)'});

function color(page: Page, column: string, apply: string, arg: unknown): Promise<void> {
  return changeTable(page, ([c, how, a]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const colors = col.meta.colors;
    // a null argument switches the type on and names no colors, so the stored ones stay
    if (how === 'linear')
      colors.setLinear(a == null ? null : (a as any).scheme, a == null ? null : (a as any).range);
    else if (how === 'conditional')
      colors.setConditional(a as any);
    else if (how === 'categorical')
      colors.setCategorical(a as any);
    else if (how === 'linked') {
      if (!df.col(a as string))
        throw new Error(`no "${a}" column to link the coloring of "${c}" to`);
      col.setTag('.color-coding-type', 'Linked');
      col.setTag('.%color-coding-linked-column-name', a as string);
    }
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

export const colorAgain = When('user colors {string} column {word} again', (page: Page, column: string, how: string) => {
  const apply = {linearly: 'linear', conditionally: 'conditional', categorically: 'categorical'}[how];
  if (!apply)
    throw new Error(`a column is colored linearly, conditionally or categorically again — not "${how}"`);
  return color(page, column, apply, null);
}, {tier: 'api', description: 'switches the type back on without naming colors, so what the column already stored has to survive'});

export const colorLinked = When('user colors {string} column linked to {string} column', (page: Page, column: string, source: string) =>
  color(page, column, 'linked', source), {tier: 'api', description: 'the column takes the source column\'s colors, as the Color Coding menu\'s Linked type does'});

export const colorLinkedText = When('user colors the text of {string} column linked to {string} column', async (page: Page, column: string, source: string) => {
  await color(page, column, 'linked', source);
  await page.evaluate((c) => {
    const gc = grok.shell.tv?.grid?.col(c);
    if (!gc)
      throw new Error(`no "${c}" column in the grid`);
    gc.isTextColorCoded = true;
  }, column);
}, {tier: 'api', description: 'the Linked type with "Apply to: Text" — the grid column paints the letters, not the cell'});

export const colorPickUp = When('user applies the coloring of {string} column to {string} column', (page: Page, from: string, to: string) =>
  changeTable(page, ([a, b]) => {
    const df = grok.shell.t;
    for (const name of [a, b])
      if (!df.col(name))
        throw new Error(`no "${name}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    // Pick Up / Apply copies every color-coding tag the source carries, and nothing else
    for (const tag of Array.from(df.col(a).tags.keys()).filter((k: any) => String(k).includes('color-coding')))
      df.col(b).setTag(tag as string, df.col(a).getTag(tag as string));
  }, [from, to] as [string, string]),
{tier: 'api', description: 'the grid header menu\'s Color Coding > Pick Up / Apply pair — every color-coding tag of the source lands on the target'});

export const colorInverted = When('user inverts the color scheme of {string} column', (page: Page, column: string) =>
  changeTable(page, (c) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}`);
    const scheme = col.getTag('.color-coding-linear');
    if (!scheme)
      throw new Error(`"${c}" has no linear color scheme to invert`);
    col.setTag('.color-coding-linear', JSON.stringify((JSON.parse(scheme) as unknown[]).reverse()));
  }, column), {tier: 'api', description: 'the arrows icon next to the scheme — the stops in the opposite order'});

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

const CODING = {linearly: 'Linear', conditionally: 'Conditional', categorically: 'Categorical', linked: 'Linked'};

export const colorCodedAs = Then('{string} column should be color-coded {word}', async (page: Page, column: string, how: string) => {
  const want = CODING[how as keyof typeof CODING];
  if (!want)
    throw new Error(`a column is color-coded ${Object.keys(CODING).join(', ')} — not "${how}"`);
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe(want);
}, {description: 'linearly, conditionally, categorically or linked — the column\'s color-coding type tag'});

export const colorLinkedTo = Then('the coloring of {string} column should be linked to {string} column', async (page: Page, column: string, source: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}`);
    return `${col.getTag('.color-coding-type') ?? 'Off'}/${col.getTag('.%color-coding-linked-column-name') ?? ''}`;
  }, column), {message: `the coloring of "${column}"`}).toBe(`Linked/${source}`);
}, {description: 'the type and the source column a Linked coloring names'});

export const textColorCoded = Then('the text of {string} column should be color-coded', async (page: Page, column: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const gc = grok.shell.tv?.grid?.col(c);
    if (!gc)
      throw new Error(`no "${c}" column in the grid`);
    return gc.isTextColorCoded === true;
  }, column), {message: `"${column}" paints its text`}).toBe(true);
}, {description: 'the grid column applies its coloring to the letters, not the cell background'});

export const colorSchemeIs = Then('the color scheme of {string} column should be {string}', async (page: Page, column: string, colors: string) => {
  const want = list(colors).map((c) => c.replace(/^#/, '').toUpperCase());
  await expect.poll(() => page.evaluate((c) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}`);
    const tag = col.getTag('.color-coding-linear');
    if (!tag)
      return ['(no linear scheme)'];
    // the tag keeps ARGB ints from the API and #rrggbb strings from the picker
    return (JSON.parse(tag) as unknown[]).map((v) => typeof v === 'number'
      ? (v >>> 0).toString(16).toUpperCase().padStart(8, '0').slice(2) : String(v).replace(/^#/, '').toUpperCase());
  }, column), {message: `the linear color scheme of "${column}"`}).toEqual(want);
}, {description: 'the stops of a linear scheme in order, comma-separated #rrggbb — an inverted scheme reads back reversed'});

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

// --- the filter panel -------------------------------------------------------------------------------

/** The filters of the current view's filter panel, by column name. */
const filterColumns = (page: Page): Promise<string[]> => page.evaluate(() => {
  // getFiltersGroup creates one when the view has none, which would make "0 filters" resurrect a panel
  const open = Array.from(grok.shell.tv?.viewers ?? []).some((v: any) => String(v.type) === 'Filters');
  const group = open ? grok.shell.tv.getFiltersGroup({createDefaultFilters: false}) : null;
  return group ? (group.filters as any[]).map((f) => String(f.columnName ?? f.column?.name ?? '')) : [];
});

export const filterPanelCount = Then('the filter panel should have {int} filter(s)', async (page: Page, count: number) => {
  await expect.poll(() => filterColumns(page), {message: 'filters of the filter panel'}).toHaveLength(count);
}, {description: 'the filters the current view\'s filter panel holds (none is created for the check)'});

export const filterPanelHas = Then('the filter panel should have a filter on {string} column', async (page: Page, column: string) => {
  await expect.poll(() => filterColumns(page), {message: 'filters of the filter panel'}).toContain(column);
});

export const filterIsExactlyContains = Then('the filter should pass exactly the rows where {string} contains {string}', async (page: Page, column: string, text: string) => {
  await expect.poll(() => page.evaluate(([c, t]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let wrong = 0;
    let matching = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const hit = String(col.get(i) ?? '').includes(t);
      if (hit)
        matching++;
      if (hit !== df.filter.get(i))
        wrong++;
    }
    return matching === 0 ? `no row contains "${t}"` : wrong === 0 ? 'exactly' : `${wrong} rows off`;
  }, [column, text] as [string, string]), {message: `the filter against rows where ${column} contains "${text}"`}).toBe('exactly');
}, {description: 'every row whose value contains the text passes and no other; a text no row contains fails'});

// --- selection facts ---------------------------------------------------------------------------------

export const selectedRowCount = Then('{int} row(s) should be selected', (page: Page, count: number) =>
  expect.poll(() => selectedCount(page), {message: 'rows selected'}).toBe(count), {description: 'exactly that many, whichever rows'});

export const rowsRangeSelected = Then('rows {int} to {int} should be selected', async (page: Page, from: number, to: number) => {
  await expect.poll(() => page.evaluate(([a, b]) => {
    const df = grok.shell.t;
    let inRange = 0;
    for (let i = a - 1; i < b && i < df.rowCount; i++) {
      if (df.selection.get(i))
        inRange++;
    }
    return `${inRange} of ${b - a + 1} in the range, ${df.selection.trueCount} in all`;
  }, [from, to] as [number, number]), {message: `rows ${from} to ${to} selected`}).toBe(`${to - from + 1} of ${to - from + 1} in the range, ${to - from + 1} in all`);
}, {description: 'rows counted from 1 as the grid shows them, every row of the range and nothing else'});

export const allRowsSelected = Then('all rows should be selected', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.selection.trueCount === grok.shell.t.rowCount), {message: 'every row is selected'}).toBe(true));

/** The names of the selected grid columns of the current table view, in grid order. */
const selectedColumns = (page: Page): Promise<string[]> => page.evaluate(() => {
  const cols = grok.shell.tv.grid.columns;
  const out: string[] = [];
  for (let i = 0; i < cols.length; i++) {
    const c = cols.byIndex(i);
    if (c?.selected && c.name)
      out.push(String(c.name));
  }
  return out;
});

export const columnsSelected = Then('columns {string} should be selected', async (page: Page, names: string) => {
  await expect.poll(async () => (await selectedColumns(page)).slice().sort(), {message: 'the selected grid columns'})
    .toEqual(list(names).slice().sort());
}, {description: 'exactly these grid columns are selected (comma-separated, in any order — a selection is a set) — a Shift or Control click on a header'});

export const noColumnsSelected = Then('no columns should be selected', async (page: Page) => {
  await expect.poll(() => selectedColumns(page), {message: 'the selected grid columns'}).toEqual([]);
});

export const selectedPassFilter = Then('every selected row should pass the filter', async (page: Page) => {
  const off = await page.evaluate(() => {
    const df = grok.shell.t;
    let wrong = 0;
    for (const i of df.selection.getSelectedIndexes() as Iterable<number>) {
      if (!df.filter.get(i))
        wrong++;
    }
    return wrong;
  });
  expect(off, 'selected rows the filter drops').toBe(0);
});

export const selectFirstRows = When('user selects the first {int} rows', (page: Page, count: number) =>
  changeTable(page, (n: number) => { grok.shell.t.selection.init((i: number) => i < n); }, count),
{tier: 'api', description: 'the first N rows in table order become the selection, nothing else'});

// --- other tables and links ------------------------------------------------------------------------

function tableCount(page: Page, name: string, what: 'filter' | 'selection'): Promise<number> {
  return page.evaluate(([n, w]) => {
    const t = grok.shell.tables.find((x: any) => x.name === n);
    if (!t)
      throw new Error(`table "${n}" is not open; open: ${grok.shell.tables.map((x: any) => x.name).join(' | ')}`);
    return (w === 'filter' ? t.filter : t.selection).trueCount as number;
  }, [name, what] as [string, 'filter' | 'selection']);
}

export const tableFilterCount = Then('{int} rows of table {string} should pass the filter', (page: Page, count: number, name: string) =>
  expect.poll(() => tableCount(page, name, 'filter'), {message: `rows of "${name}" passing its filter`}).toBe(count),
{description: 'another open table\'s filter, without switching to its view — what a link between tables carries'});

export const tableSelectedCount = Then('{int} rows of table {string} should be selected', (page: Page, count: number, name: string) =>
  expect.poll(() => tableCount(page, name, 'selection'), {message: `rows of "${name}" selected`}).toBe(count));

export const linkTables = When('user links table {string} to table {string} as {string}:', async (page: Page, from: string, to: string, type: string, table: string[][]) => {
  await page.evaluate(([a, b, t, keys]) => {
    const find = (n: string) => {
      const df = grok.shell.tables.find((x: any) => x.name === n);
      if (!df)
        throw new Error(`table "${n}" is not open; open: ${grok.shell.tables.map((x: any) => x.name).join(' | ')}`);
      return df;
    };
    const norm = (s: string) => s.toLowerCase().replace(/[^a-z]/g, '');
    const sync = (Object.values(DG.SYNC_TYPE) as string[]).find((s) => norm(s) === norm(t));
    if (!sync)
      throw new Error(`no "${t}" link type; the platform has: ${(Object.values(DG.SYNC_TYPE) as string[]).join(', ')}`);
    grok.data.linkTables(find(a), find(b), keys.map((k) => k[0]), keys.map((k) => k[1]), [sync]);
  }, [from, to, type, table] as [string, string, string, string[][]]);
}, {tier: 'api', description: '| key column in the first | key column in the second | rows; the type as the platform names it: "filter to filter", "selection to filter", "selection to selection", "row to row"…; the link lives as long as the tables'});

export const columnTrueWhereFiltered = Then('{string} column should be true exactly where the filter passes', async (page: Page, column: string) => {
  const off = await page.evaluate((c) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let wrong = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if ((col.get(i) === true) !== df.filter.get(i))
        wrong++;
    }
    return wrong;
  }, column);
  expect(off, `rows where "${column}" disagrees with the filter`).toBe(0);
}, {description: 'a boolean column written from the filter (Filter to Column) row by row against the filter bitset'});

export const categoricalColorIs = Then('the categorical color of {string} in {string} column should be {string}', async (page: Page, category: string, column: string, color: string) => {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  await expect.poll(() => page.evaluate(([c, cat]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const tag = col.getTag('.color-coding-categorical');
    const map: Record<string, unknown> = tag ? JSON.parse(tag) : {};
    const v = map[cat];
    if (v === undefined)
      return `(no color for ${cat}; colored: ${Object.keys(map).join(', ') || 'nothing'})`;
    // the tag keeps what was written: an ARGB int from the API, a #rrggbb string from the picker
    const n = typeof v === 'number' ? v : parseInt(String(v).replace(/^#/, ''), 16);
    return '#' + (n >>> 0 & 0xFFFFFF).toString(16).padStart(6, '0').toUpperCase();
  }, [column, category] as [string, string]), {message: `the categorical color of ${category} in "${column}"`}).toBe(want);
}, {description: 'the column\'s categorical color-coding map (what a legend picker writes), as #rrggbb'});

// --- the filter panel through the API ------------------------------------------------------------------

/** The current view's filter group, creating the panel with or without its default cards. */
async function filterGroup(page: Page, defaults: boolean): Promise<void> {
  await page.evaluate((d) => { grok.shell.tv.getFiltersGroup({createDefaultFilters: d}); }, defaults);
  await page.locator('[name="viewer-Filters"]').filter({visible: true}).first().waitFor();
  await settleAll(page);
}

export const openFilterPanel = When('user opens the filter panel', (page: Page) => filterGroup(page, true),
  {tier: 'api', description: 'the view\'s Filters viewer with a default card per column, as the ribbon\'s filter icon opens it; done when every viewer has drawn'});

export const openEmptyFilterPanel = When('user opens an empty filter panel', (page: Page) => filterGroup(page, false),
  {tier: 'api', description: 'the Filters viewer with no cards — the entry point for scenarios that add every card themselves'});

/** A card criterion through the group's own API — the state a category click or a range drag
 * would leave; the gesture itself is the filter panel features' subject. */
function filterState(page: Page, state: Record<string, unknown>): Promise<void> {
  return changeTable(page, (s: Record<string, unknown>) => {
    const fg = grok.shell.tv.getFiltersGroup({createDefaultFilters: false});
    const col = grok.shell.t.col(String(s.column ?? ''));
    if (s.column !== undefined && !col)
      throw new Error(`no "${s.column}" column in ${grok.shell.t.name}; it has: ${grok.shell.t.columns.names().join(', ')}`);
    fg.updateOrAdd(s, true);
  }, state);
}

export const addCategoricalFilter = When('user adds a categorical filter on {string} keeping {string}', (page: Page, column: string, values: string) =>
  filterState(page, {type: 'categorical', column, selected: list(values)}),
{tier: 'api', description: 'a categorical card with these categories (comma-separated) checked and the rest unchecked, as the category clicks leave it; updates the card when it exists'});

export const addRangeFilter = When('user adds a range filter on {string} from {float} to {float}', (page: Page, column: string, min: number, max: number) =>
  filterState(page, {type: 'histogram', column, min, max}), {tier: 'api', description: 'a histogram card narrowed to the range, as the handles leave it'});

export const configureHierarchical = When('user configures the hierarchical filter with columns {string}', (page: Page, columns: string) =>
  filterState(page, {type: 'hierarchical', colNames: list(columns), allEnabled: true}),
{tier: 'api', description: 'the hierarchical card\'s levels in this order (comma-separated), every node checked'});

/** The filters of the current view's panel with their state, by column name. */
const filterFacts = (page: Page): Promise<{column: string; filtering: boolean; selected?: string[]}[]> => page.evaluate(() => {
  const group = grok.shell.tv?.getFiltersGroup?.({createDefaultFilters: false});
  return group ? (group.filters as any[]).map((f) => ({
    column: String(f.columnName ?? f.column?.name ?? ''),
    filtering: f.isFiltering === true,
    selected: (() => {
      try {
        const s = f.saveState?.();
        return Array.isArray(s?.selected) ? s.selected.map(String) : undefined;
      }
      catch {
        return undefined;
      }
    })(),
  })) : [];
});

export const filterKeepsOnly = Then('the filter on {string} column should keep only {string}', async (page: Page, column: string, values: string) => {
  await expect.poll(async () => {
    const f = (await filterFacts(page)).find((x) => x.column === column);
    return f ? (f.selected ? [...f.selected].sort().join(', ') : '(no categories)') : `(no filter on ${column})`;
  }, {message: `the categories the filter on "${column}" keeps`}).toBe(list(values).sort().join(', '));
}, {description: 'the categorical card\'s checked categories (comma-separated, any order) — its own state, not the table\'s rows'});

export const filterIsFiltering = Then('the filter on {string} column should be filtering', async (page: Page, column: string) => {
  await expect.poll(async () => {
    const f = (await filterFacts(page)).find((x) => x.column === column);
    return f ? (f.filtering ? 'filtering' : 'not filtering') : `(no filter on ${column})`;
  }, {message: `the filter on "${column}"`}).toBe('filtering');
}, {description: 'the card restricts rows (its criterion is not "everything") — the header counter counts these'});

export const filterNotFiltering = Then('the filter on {string} column should not be filtering', async (page: Page, column: string) => {
  await expect.poll(async () => {
    const f = (await filterFacts(page)).find((x) => x.column === column);
    return f ? (f.filtering ? 'filtering' : 'not filtering') : `(no filter on ${column})`;
  }, {message: `the filter on "${column}"`}).toBe('not filtering');
});

export const onlyStartingWithSelected = Then('only rows where {string} starts with {string} should be selected', async (page: Page, column: string, prefix: string) => {
  await expect.poll(() => page.evaluate(([c, p]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    let matching = 0;
    let wrong = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const hit = String(col.get(i) ?? '').startsWith(p);
      if (hit)
        matching++;
      if (hit !== df.selection.get(i))
        wrong++;
    }
    return matching === 0 ? `no row starts with "${p}"` : wrong === 0 ? 'exactly' : `${wrong} rows off`;
  }, [column, prefix] as [string, string]), {message: `the selection against rows where ${column} starts with "${prefix}"`}).toBe('exactly');
}, {description: 'every row whose value starts with the text is selected and no other'});
