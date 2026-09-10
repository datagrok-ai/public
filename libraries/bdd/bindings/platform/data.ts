/* The current table's rows and columns through the JS API: selection, filter, the current row,
   deleting rows, calculated columns, column colors, and the tables the workspace holds. A step
   that changes rows or colors takes every open viewer's baseline first, so the viewer checks that
   follow ("should have repainted", "a narrower value range than before") read against the state
   before the change, and returns once every viewer has drawn the change, so the next step's
   baseline is the state after it. The in-page runtime's `table`, `col`, `rowFacts` and `setRows`
   do the lookups and the row scans, and say which column or category is missing. */
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
import {Then, When} from '../../src/registry.js';
import {baselineAll, evaluate, RowFacts, RowTest, settleAll} from '../../src/runtime/viewers.js';

declare const grok: any;
declare const DG: any;

const list = (s: string): string[] => s.split(/\s*,\s*/).filter((x) => x.length > 0);
const facts = (page: Page, column: string, test: RowTest, mustMatch = true): Promise<RowFacts> =>
  evaluate(page, ([c, t, m]) => (window as any).__bdd.rowFacts(c, t, m), [column, test, mustMatch] as [string, RowTest, boolean]);
const selectedCount = (page: Page): Promise<number> => page.evaluate(() => grok.shell.t.selection.trueCount as number);
const filteredCount = (page: Page): Promise<number> => page.evaluate(() => grok.shell.t.filter.trueCount as number);

/** A change to the table every viewer answers: baselines first, the change, then every viewer
 * has drawn it. */
async function changeTable(page: Page, body: (arg: any) => void, arg: unknown): Promise<void> {
  await baselineAll(page);
  await page.evaluate(body as (arg: unknown) => void, arg);
  await settleAll(page);
}

/** The selection or the filter becomes exactly the rows a test names; a category that matches no
 * row fails (a typo), unless the step allows an empty result. */
const setRows = (page: Page, what: 'selection' | 'filter', column: string, test: RowTest, negate = false, mustMatch = true): Promise<void> =>
  changeTable(page, ([w, c, t, n, m]) => { (window as any).__bdd.setRows(w, c, t, n, m); },
    [what, column, test, negate, mustMatch] as [typeof what, string, RowTest, boolean, boolean]);

// --- selection ---------------------------------------------------------------------------------------

export const clearSelection = When('user clears the row selection', (page: Page) =>
  changeTable(page, () => { grok.shell.t.selection.setAll(false); }, null), {tier: 'api'});

export const selectWhereIs = When('user selects rows where {string} is {string}', (page: Page, column: string, value: string) =>
  setRows(page, 'selection', column, {eq: value}), {tier: 'api', description: 'the table\'s selection bitset — every row of the category and nothing else; the UI path is a click on the category in a viewer'});

export const selectWhereOneOf = When('user selects rows where {string} is one of {string}', (page: Page, column: string, values: string) =>
  setRows(page, 'selection', column, {in: list(values)}), {tier: 'api', description: 'comma-separated categories'});

export const selectWhereBetween = When('user selects rows where {string} is between {float} and {float}', (page: Page, column: string, min: number, max: number) =>
  setRows(page, 'selection', column, {between: [min, max]}),
{tier: 'api', description: 'a numeric range, both ends included — the selection becomes exactly those rows; a range no row matches fails'});

export const selectFirstRows = When('user selects the first {int} rows', (page: Page, count: number) =>
  changeTable(page, (n: number) => { grok.shell.t.selection.init((i: number) => i < n); }, count),
{tier: 'api', description: 'the first N rows in table order become the selection, nothing else'});

export const selectAllRows = When('user selects all rows', (page: Page) =>
  page.evaluate(() => { grok.shell.t.selection.setAll(true); }),
{tier: 'api', description: 'the whole table selected — "all rows should be selected" had no counterpart to reach it'});

export const selectNoRows = When('user selects no rows', (page: Page) =>
  page.evaluate(() => { grok.shell.t.selection.setAll(false); }), {tier: 'api', description: 'the selection cleared'});

export const noneSelected = Then('no rows should be selected', (page: Page) =>
  expect.poll(() => selectedCount(page), {message: 'rows selected'}).toBe(0));

export const someSelected = Then('some rows should be selected', (page: Page) =>
  expect.poll(() => selectedCount(page), {message: 'no row is selected'}).toBeGreaterThan(0));

export const selectedRowCount = Then('{int} row(s) should be selected', (page: Page, count: number) =>
  expect.poll(() => selectedCount(page), {message: 'rows selected'}).toBe(count), {description: 'exactly that many, whichever rows'});

export const allRowsSelected = Then('all rows should be selected', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.selection.trueCount === grok.shell.t.rowCount), {message: 'every row is selected'}).toBe(true));

export const allOfSelected = Then('all rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => {
    const m = await facts(page, column, {eq: value});
    return `${m.selected} of ${m.matching}`;
  }, {message: `rows where ${column} is ${value} selected`}).toMatch(/^(\d+) of \1$/);
}, {description: 'every row of the category, whatever else is selected'});

/** Every row the test names and nothing else. */
async function expectOnlySelected(page: Page, column: string, test: RowTest, what: string): Promise<void> {
  await expect.poll(async () => {
    const m = await facts(page, column, test);
    return m.wrongSelection === 0 ? 'exactly' : `${m.selected} of ${m.matching}, ${m.selectedTotal} in all`;
  }, {message: `rows where ${column} ${what} selected`}).toBe('exactly');
}

export const onlyOfSelected = Then('only rows where {string} is {string} should be selected', (page: Page, column: string, value: string) =>
  expectOnlySelected(page, column, {eq: value}, `is ${value}`), {description: 'every row of the category and nothing else'});

export const onlyOfAnySelected = Then('only rows where {string} is one of {string} should be selected', (page: Page, column: string, values: string) =>
  expectOnlySelected(page, column, {in: list(values)}, `is one of ${values}`),
{description: 'every row of these categories (comma-separated) and nothing else — a union built with Control clicks'});

export const onlyStartingWithSelected = Then('only rows where {string} starts with {string} should be selected', (page: Page, column: string, prefix: string) =>
  expectOnlySelected(page, column, {startsWith: prefix}, `starts with "${prefix}"`), {description: 'every row whose value starts with the text is selected and no other'});

export const noneOfSelected = Then('no rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await facts(page, column, {eq: value})).selected, {message: `rows where ${column} is ${value} selected`}).toBe(0);
});

export const someOfSelected = Then('some rows where {string} is {string} should be selected', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await facts(page, column, {eq: value})).selected, {message: `rows where ${column} is ${value} selected`}).toBeGreaterThan(0);
});

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

// --- the current row and column -------------------------------------------------------------------------

export const hasCurrentRow = Then('the table should have a current row', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.currentRowIdx as number), {message: 'the current row index'}).toBeGreaterThanOrEqual(0));

export const currentRowValue = Then('{string} of the current row should be {string}', async (page: Page, column: string, value: string) => {
  await expect.poll(() => evaluate(page, (c) => {
    const df = (window as any).__bdd.table();
    const col = (window as any).__bdd.col(c);
    return df.currentRowIdx < 0 ? '(no current row)' : String(col.get(df.currentRowIdx) ?? '');
  }, column), {message: `"${column}" of the current row`}).toBe(value);
}, {description: 'the column\'s value in the current row, as text'});

export const currentColumnIs = Then('the current column should be {string}', (page: Page, name: string) =>
  expect.poll(() => page.evaluate(() => String(grok.shell.t.currentCol?.name ?? '')), {message: 'the current column'}).toBe(name),
{description: 'the table\'s current column by name ("" when none)'});

// --- filter ------------------------------------------------------------------------------------------

export const filterBetween = When('user filters rows where {string} is between {float} and {float}', (page: Page, column: string, lo: number, hi: number) =>
  setRows(page, 'filter', column, {between: [lo, hi]}, false, false),
{tier: 'api', description: 'the table\'s filter bitset, as a filter viewer would set it; a range no row falls in empties the table on purpose'});

export const filterNotNull = When('user filters rows where {string} is not null', (page: Page, column: string) =>
  setRows(page, 'filter', column, {notNull: true}, false, false), {tier: 'api', description: 'the table\'s filter bitset, as a filter viewer would set it'});

export const filterTo = When('user filters rows where {string} is {string}', (page: Page, column: string, value: string) =>
  setRows(page, 'filter', column, {eq: value}), {tier: 'api', description: 'keeps the category\'s rows only — the table\'s filter bitset'});

export const filterToAnyOf = When('user filters rows where {string} is one of {string}', (page: Page, column: string, values: string) =>
  setRows(page, 'filter', column, {in: list(values)}), {tier: 'api', description: 'keeps the rows of these categories (comma-separated) and no other'});

export const filterOut = When('user filters out rows where {string} is {string}', (page: Page, column: string, value: string) =>
  setRows(page, 'filter', column, {eq: value}, true), {tier: 'api', description: 'the table\'s filter bitset'});

export const resetFilter = When('user resets the filter', (page: Page) =>
  changeTable(page, () => { grok.shell.t.filter.setAll(true); }, null), {tier: 'api'});

export const filterPasses = Then('{int} row(s) should pass the filter', (page: Page, count: number) =>
  expect.poll(() => filteredCount(page), {message: 'rows passing the filter'}).toBe(count));

export const filterPassesFewer = Then('fewer than {int} rows should pass the filter', (page: Page, count: number) =>
  expect.poll(() => filteredCount(page), {message: 'rows passing the filter'}).toBeLessThan(count));

export const filterPassesAll = Then('all rows should pass the filter', (page: Page) =>
  expect.poll(() => page.evaluate(() => grok.shell.t.filter.trueCount === grok.shell.t.rowCount), {message: 'every row passes the filter'}).toBe(true));

/** The table filter bit by bit: the rows the test names pass, no other row does. */
async function expectFilterExactly(page: Page, column: string, test: RowTest, what: string): Promise<void> {
  await expect.poll(async () => {
    const m = await facts(page, column, test);
    return m.wrongFilter === 0 ? 'exactly' : `${m.wrongFilter} rows off`;
  }, {message: `the filter against rows where ${column} ${what}`}).toBe('exactly');
}

export const filterIsExactly = Then('the filter should pass exactly the rows where {string} is between {float} and {float}',
  (page: Page, column: string, lo: number, hi: number) => expectFilterExactly(page, column, {between: [lo, hi]}, `is in [${lo}, ${hi}]`),
  {description: 'the table filter bit by bit — a viewer\'s own filter must leave it alone'});

export const filterIsExactlyCategory = Then('the filter should pass exactly the rows where {string} is {string}', (page: Page, column: string, value: string) =>
  expectFilterExactly(page, column, {eq: value}, `is ${value}`), {description: 'the table filter bit by bit: the category\'s rows pass, no other row does'});

export const filterIsExactlyContains = Then('the filter should pass exactly the rows where {string} contains {string}', (page: Page, column: string, text: string) =>
  expectFilterExactly(page, column, {contains: text}, `contains "${text}"`), {description: 'every row whose value contains the text passes and no other; a text no row contains fails'});

export const noneOfFiltered = Then('no rows where {string} is {string} should pass the filter', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => (await facts(page, column, {eq: value})).passing, {message: `rows where ${column} is ${value} passing the filter`}).toBe(0);
}, {description: 'the category is filtered out — what the filter actually keeps, not what a card says it keeps'});

export const allOfFiltered = Then('all rows where {string} is {string} should pass the filter', async (page: Page, column: string, value: string) => {
  await expect.poll(async () => {
    const m = await facts(page, column, {eq: value});
    return m.passing === m.matching ? 'all' : `${m.passing} of ${m.matching}`;
  }, {message: `rows where ${column} is ${value} passing the filter`}).toBe('all');
});

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
  await expect.poll(async () => (await facts(page, column, {eq: value}, false)).matching, {message: `rows where ${column} is ${value}`}).toBe(0);
});

// --- columns -----------------------------------------------------------------------------------------

export const setCell = When('user sets {string} column in row {int} to {string}', (page: Page, column: string, row: number, value: string) =>
  changeTable(page, ([c, r, v]) => {
    const b = (window as any).__bdd;
    const df = b.table();
    const col = b.col(c);
    const i = r - 1;
    if (i < 0 || i >= df.rowCount)
      throw new Error(`row ${r} is outside ${df.name} (${df.rowCount} rows)`);
    if (v === '')
      col.set(i, null);
    else if (col.type === 'double' || col.type === 'int' || col.type === 'float' || col.type === 'bigint')
      col.set(i, v === 'NaN' ? Number.NaN : v === 'Infinity' ? Number.POSITIVE_INFINITY : v === '-Infinity' ? Number.NEGATIVE_INFINITY : Number(v));
    else
      col.set(i, v);
  }, [column, row, value] as [string, number, string]),
{tier: 'api', description: 'one cell, 1-based row; on a numeric column "NaN", "Infinity" and "-Infinity" write those values and "" writes a blank — what a viewer must survive'});

export const addCalculated = When('user adds a calculated column {string} with formula {string}', async (page: Page, name: string, formula: string) => {
  await page.evaluate(async ([n, f]) => { await grok.shell.t.columns.addNewCalculated(n, f); }, [name, formula] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.t.columns.names().includes(n), name), {message: `"${name}" in the table's columns`}).toBe(true);
}, {tier: 'api', description: 'a formula in the platform\'s syntax: ${HEIGHT} * 2'});

export const removeColumn = When('user removes {string} column', (page: Page, name: string) =>
  evaluate(page, (n) => { const b = (window as any).__bdd; b.col(n); b.table().columns.remove(n); }, name), {tier: 'api'});

export const renameColumn = When('user renames {string} column to {string}', async (page: Page, from: string, to: string) => {
  await evaluate(page, ([a, b]) => { (window as any).__bdd.col(a).name = b; }, [from, to] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.t.columns.names().includes(n), to), {message: `"${to}" in the table's columns`}).toBe(true);
}, {tier: 'api', description: 'the column\'s name through the API — the UI path is the header\'s Column Properties dialog'});

// --- column colors -------------------------------------------------------------------------------------

/** `#rrggbb` → the ARGB int the color API takes. */
function argb(hex: string): number {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex);
  if (!m)
    throw new Error(`"${hex}" is not a #rrggbb color`);
  return (0xFF000000 | parseInt(m[1], 16)) >>> 0;
}

function color(page: Page, column: string, apply: string, arg: unknown): Promise<void> {
  return changeTable(page, ([c, how, a]) => {
    const b = (window as any).__bdd;
    const col = b.col(c);
    const colors = col.meta.colors;
    // a null argument switches the type on and names no colors, so the stored ones stay
    if (how === 'linear')
      colors.setLinear(a == null ? null : (a as any).scheme, a == null ? null : (a as any).range);
    else if (how === 'conditional')
      colors.setConditional(a as any);
    else if (how === 'categorical')
      colors.setCategorical(a as any);
    else if (how === 'linked') {
      b.col(a as string);
      col.setTag('.color-coding-type', 'Linked');
      col.setTag('.%color-coding-linked-column-name', a as string);
    }
    else
      colors.setDisabled();
  }, [column, apply, arg] as [string, string, unknown]);
}

export const colorLinear = When('user colors {string} column linearly from {string} to {string}', (page: Page, column: string, from: string, to: string) =>
  color(page, column, 'linear', {scheme: [argb(from), argb(to)]}),
{tier: 'api', description: 'the column\'s color coding, as the grid\'s Color Coding menu sets it; #rrggbb colors'});

export const colorLinearOver = When('user colors {string} column linearly from {string} to {string} over {float} to {float}',
  (page: Page, column: string, from: string, to: string, min: number, max: number) =>
    color(page, column, 'linear', {scheme: [argb(from), argb(to)], range: {min, max}}), {tier: 'api'});

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

/** Whether the grid column paints its text, set or read; a column the grid does not show fails. */
const textColorCoding = (page: Page, column: string, set?: boolean): Promise<boolean> => page.evaluate(([c, s]) => {
  const gc = grok.shell.tv?.grid?.col(c);
  if (!gc)
    throw new Error(`no "${c}" column in the grid`);
  if (s !== null)
    gc.isTextColorCoded = s;
  return gc.isTextColorCoded === true;
}, [column, set ?? null] as [string, boolean | null]);

export const colorLinkedText = When('user colors the text of {string} column linked to {string} column', async (page: Page, column: string, source: string) => {
  await color(page, column, 'linked', source);
  await textColorCoding(page, column, true);
}, {tier: 'api', description: 'the Linked type with "Apply to: Text" — the grid column paints the letters, not the cell'});

export const colorPickUp = When('user applies the coloring of {string} column to {string} column', (page: Page, from: string, to: string) =>
  changeTable(page, ([a, b]) => {
    const bdd = (window as any).__bdd;
    const src = bdd.col(a);
    const dst = bdd.col(b);
    // Pick Up / Apply copies every color-coding tag the source carries, and nothing else
    for (const tag of Array.from(src.tags.keys()).filter((k: any) => String(k).includes('color-coding')))
      dst.setTag(tag as string, src.getTag(tag as string));
  }, [from, to] as [string, string]),
{tier: 'api', description: 'the grid header menu\'s Color Coding > Pick Up / Apply pair — every color-coding tag of the source lands on the target'});

export const colorInverted = When('user inverts the color scheme of {string} column', (page: Page, column: string) =>
  changeTable(page, (c) => {
    const col = (window as any).__bdd.col(c);
    const scheme = col.getTag('.color-coding-linear');
    if (!scheme)
      throw new Error(`"${c}" has no linear color scheme to invert`);
    col.setTag('.color-coding-linear', JSON.stringify((JSON.parse(scheme) as unknown[]).reverse()));
  }, column), {tier: 'api', description: 'the arrows icon next to the scheme — the stops in the opposite order'});

/** The column's color-coding type tag: '' when none. */
const colorCodingType = (page: Page, column: string): Promise<string> => evaluate(page, (c) => {
  const type = (window as any).__bdd.col(c).getTag('.color-coding-type');
  return type == null || type === 'Off' ? '' : String(type);
}, column);

const CODING = {linearly: 'Linear', conditionally: 'Conditional', categorically: 'Categorical', linked: 'Linked'};

export const noColorCoding = Then('{string} column should have no color coding', async (page: Page, column: string) => {
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe('');
}, {description: 'the column\'s color-coding tag is unset or Off'});

export const colorCodedCategorically = Then('{string} column should be color-coded categorically', async (page: Page, column: string) => {
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe('Categorical');
});

export const colorCodedAs = Then('{string} column should be color-coded {word}', async (page: Page, column: string, how: string) => {
  const want = CODING[how as keyof typeof CODING];
  if (!want)
    throw new Error(`a column is color-coded ${Object.keys(CODING).join(', ')} — not "${how}"`);
  await expect.poll(() => colorCodingType(page, column), {message: `the color coding of "${column}"`}).toBe(want);
}, {description: 'linearly, conditionally, categorically or linked — the column\'s color-coding type tag'});

export const colorLinkedTo = Then('the coloring of {string} column should be linked to {string} column', async (page: Page, column: string, source: string) => {
  await expect.poll(() => evaluate(page, (c) => {
    const col = (window as any).__bdd.col(c);
    return `${col.getTag('.color-coding-type') ?? 'Off'}/${col.getTag('.%color-coding-linked-column-name') ?? ''}`;
  }, column), {message: `the coloring of "${column}"`}).toBe(`Linked/${source}`);
}, {description: 'the type and the source column a Linked coloring names'});

export const textColorCoded = Then('the text of {string} column should be color-coded', async (page: Page, column: string) => {
  await expect.poll(() => textColorCoding(page, column), {message: `"${column}" paints its text`}).toBe(true);
}, {description: 'the grid column applies its coloring to the letters, not the cell background'});

export const colorSchemeIs = Then('the color scheme of {string} column should be {string}', async (page: Page, column: string, colors: string) => {
  const want = list(colors).map((c) => c.replace(/^#/, '').toUpperCase());
  await expect.poll(() => evaluate(page, (c) => {
    const tag = (window as any).__bdd.col(c).getTag('.color-coding-linear');
    if (!tag)
      return ['(no linear scheme)'];
    // the tag keeps ARGB ints from the API and #rrggbb strings from the picker
    return (JSON.parse(tag) as unknown[]).map((v) => typeof v === 'number'
      ? (v >>> 0).toString(16).toUpperCase().padStart(8, '0').slice(2) : String(v).replace(/^#/, '').toUpperCase());
  }, column), {message: `the linear color scheme of "${column}"`}).toEqual(want);
}, {description: 'the stops of a linear scheme in order, comma-separated #rrggbb — an inverted scheme reads back reversed'});

export const categoricalColorIs = Then('the categorical color of {string} in {string} column should be {string}', async (page: Page, category: string, column: string, color: string) => {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  await expect.poll(() => evaluate(page, ([c, cat]) => {
    const tag = (window as any).__bdd.col(c).getTag('.color-coding-categorical');
    const map: Record<string, unknown> = tag ? JSON.parse(tag) : {};
    const v = map[cat];
    if (v === undefined)
      return `(no color for ${cat}; colored: ${Object.keys(map).join(', ') || 'nothing'})`;
    // the tag keeps what was written: an ARGB int from the API, a #rrggbb string from the picker
    const n = typeof v === 'number' ? v : parseInt(String(v).replace(/^#/, ''), 16);
    return '#' + (n >>> 0 & 0xFFFFFF).toString(16).padStart(6, '0').toUpperCase();
  }, [column, category] as [string, string]), {message: `the categorical color of ${category} in "${column}"`}).toBe(want);
}, {description: 'the column\'s categorical color-coding map (what a legend picker writes), as #rrggbb'});

// --- other tables ------------------------------------------------------------------------------------

const tableInfo = (page: Page, name: string): Promise<{columns: string[]; rows: number}> => evaluate(page, (n) => {
  const t = (window as any).__bdd.tableNamed(n);
  return {columns: t.columns.names() as string[], rows: t.rowCount as number};
}, name);

export const tableOpen = Then('table {string} should be open', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate((n) => (grok.shell.tables as any[]).some((x) => x.name === n) ? 'open' :
    `not open; open: ${(grok.shell.tables as any[]).map((x) => x.name).join(' | ')}`, name), {message: `table "${name}"`}).toBe('open');
}, {description: 'in grok.shell.tables, by its exact name'});

export const tableColumns = Then('table {string} should have columns {string}', async (page: Page, name: string, columns: string) => {
  expect((await tableInfo(page, name)).columns, `columns of "${name}"`).toEqual(list(columns));
}, {description: 'exactly these, in this order (comma-separated)'});

export const tableRows = Then('table {string} should have {int} row(s)', async (page: Page, name: string, count: number) => {
  expect((await tableInfo(page, name)).rows, `rows of "${name}"`).toBe(count);
});

const missingCount = (page: Page, name: string, column: string): Promise<number> => evaluate(page, ([n, c]) => {
  const b = (window as any).__bdd;
  const t = b.tableNamed(n);
  const col = b.col(c, t);
  let missing = 0;
  for (let i = 0; i < t.rowCount; i++) {
    if (col.isNone(i))
      missing++;
  }
  return missing;
}, [name, column] as [string, string]);

export const tableColumnComplete = Then('table {string} should have no missing values in {string} column', async (page: Page, name: string, column: string) => {
  expect(await missingCount(page, name, column), `missing values in "${column}" of "${name}"`).toBe(0);
}, {description: 'a computed column is filled in, not blank — a result table whose fit failed has the right shape and no numbers'});

export const tableColumnIncomplete = Then('table {string} should have missing values in {string} column', async (page: Page, name: string, column: string) => {
  expect(await missingCount(page, name, column), `missing values in "${column}" of "${name}"`).toBeGreaterThan(0);
}, {description: 'at least one blank — the precondition of a scenario about empty categories'});

export const tableFilterCount = Then('{int} rows of table {string} should pass the filter', (page: Page, count: number, name: string) =>
  expect.poll(() => evaluate(page, (n) => (window as any).__bdd.tableNamed(n).filter.trueCount as number, name), {message: `rows of "${name}" passing its filter`}).toBe(count),
{description: 'another open table\'s filter, without switching to its view — what a link between tables carries'});

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
    if (s.column !== undefined)
      (window as any).__bdd.col(String(s.column));
    grok.shell.tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(s, true);
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
