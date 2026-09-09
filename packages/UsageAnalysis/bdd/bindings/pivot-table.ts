/* The steps only the pivot table needs. Everything else its features use is the library's `viewers`
   tier and the platform's data steps, read from the regions and readings the pivot reports
   (`group by chip <name>`, `add aggregate`, `counts`, `grid cell <r> of <col>`, `group by`,
   `aggregated rows`, `text of grid cell <r> of <col>`, … — see
   `core/client/d4/lib/src/viewers/pivot_viewer/CLAUDE.md`).

   What is left here:
   - the column picker the `+` of a tag row opens (a `ColumnComboBox` popup whose search box is
     created by the first keystroke on the focused icon) — a thin wrapper, to be deleted once the
     picker announces its readiness in the core;
   - the arithmetic: the aggregation the pivot shows, compared against an independent `groupBy`
     over the same table, which is the only check that reads the numbers rather than their shape;
   - the two popups that are not menus and not dialogs: the history menu of the command bar and
     the viewer-column combo of the Aggregate row;
   - the saved-parameters store, which lives in `localStorage` and outlives a feature.

   Should be promoted to the library, none of them pivot business:
   - `user drags the {string} area of {widget} onto the {string} area of {widget}` — the existing
     drag step works within one widget; this is the cross-widget one (a grid column header onto
     another viewer);
   The features here also use `the {string} reading of {widget} should (not )contain {string}`,
   which the tile viewer's bindings register in this package and which belongs in the library:
   membership in a reading that is a list (the chips of a row, the offered history), where pinning
   the whole ordered list would claim more than the scenario means. */
import {expect, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {el, ElementRef, exactText, viewers} from '@datagrok-libraries/bdd/runtime';

const HISTORY_KEY = 'grok-aggregation-history';

/** The pivot re-renders through an asynchronous aggregation; the settle returns when it landed. */
async function settle(page: Page, target: ElementRef): Promise<void> {
  const loc = await viewers.viewerLocator(page, target);
  await loc.evaluate((e) => (window as any).__bdd.settle(e, 3000));
}

// --- the tag rows' column picker ------------------------------------------------------------------

/** The `+` of a tag row opens a `ColumnComboBox` popup (`.d4-column-grid`): the icon keeps the
 * focus, the first keystroke creates the popup's search box, and Enter takes the name TYPED (not
 * the row under the pointer — hence the pointer leaves the popup first). */
export const addToRow = When('user adds {string} to the {string} row of pivot table viewer',
  async (page: Page, column: string, row: string) => {
    const target = el('pivot table viewer');
    const c = viewers.centerOf(await viewers.hitArea(page, target, `add ${row}`, true));
    await page.mouse.click(c.x, c.y);
    const popup = page.locator('.d4-column-grid').last();
    await popup.waitFor({state: 'visible', timeout: 10000});
    await page.mouse.move(2, 2);
    await page.keyboard.press(column[0]);
    const search = page.locator('input.d4-column-selector-search-input');
    await expect(search, `the column picker's search box of the "${row}" row`).toBeVisible({timeout: 10000});
    if (column.length > 1)
      await page.keyboard.type(column.slice(1));
    await expect(search).toHaveValue(column, {timeout: 5000});
    await page.keyboard.press('Enter');
    await popup.waitFor({state: 'detached', timeout: 5000}).catch(() => {});
    await settle(page, target);
  }, {tier: 'ui', description: 'the + of a tag row, the column typed and committed — as a user picks it'});

// --- a drag between two widgets -------------------------------------------------------------------

/** A drag from one widget's hit area to another's — a grid column header onto the pivot's Group by
 * row. The platform's drag starts on the first move after the button goes down, so the pointer
 * travels in steps; a single jump drops nothing. */
export const dragAreaOntoWidget = When('user drags the {string} area of {widget} onto the {string} area of {widget}',
  async (page: Page, from: string, source: ElementRef, to: string, target: ElementRef) => {
    const a = viewers.centerOf(await viewers.hitArea(page, source, from, true));
    const b = viewers.centerOf(await viewers.hitArea(page, target, to, true));
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    for (let i = 1; i <= 10; i++)
      await page.mouse.move(a.x + (b.x - a.x) * i / 10, a.y + (b.y - a.y) * i / 10);
    await page.mouse.up();
    await settle(page, target);
  }, {tier: 'ui', description: 'cross-widget drag and drop — should be promoted to the library'});

// --- the arithmetic --------------------------------------------------------------------------------

/** Every value the pivot draws, against an independent `groupBy` over the source table: the key of
 * each row is the text of its key cell, and every other cell is compared with the aggregation of
 * the group that key names (the shown text is rounded, so a cell matches within half of its last
 * digit). A cell the pivot leaves blank must have no value in the reference, and every reference
 * value must be somewhere on the grid. */
async function expectAggregation(page: Page, target: ElementRef, aggregation: string, keyColumn: string,
    pivotColumn: string, filteredOnly: boolean): Promise<void> {
  await viewers.installViewerRuntime(page);
  const loc = await viewers.viewerLocator(page, target);
  const parsed = /^\s*([A-Za-z0-9#]+)\s*\(\s*(.+?)\s*\)\s*$/.exec(aggregation);
  if (parsed == null)
    throw new Error(`"${aggregation}" is not an aggregation — write it as "avg(AGE)"`);
  const problems: string[] = await loc.evaluate((element, [aggType, column, keyCol, pivotCol, onlyFiltered]) => {
    const bdd = (window as any).__bdd;
    const viewer = bdd.viewerOf(element);
    let source = viewer.dataFrame;
    if (onlyFiltered === 'yes')
      source = source.clone(source.filter);
    const builder = source.groupBy([keyCol]);
    if (pivotCol !== '')
      builder.pivot(pivotCol);
    builder.add(aggType, column);
    const reference = builder.aggregate();
    const expected: Record<string, Record<string, number | null>> = {};
    const refKey = reference.col(keyCol);
    for (let r = 0; r < reference.rowCount; r++) {
      const row: Record<string, number | null> = {};
      for (const name of reference.columns.names()) {
        if (name === keyCol)
          continue;
        const col = reference.col(name);
        row[name] = col.isNone(r) ? null : col.get(r);
      }
      expected[String(refKey.get(r))] = row;
    }
    const values = viewer.getWidgetStatus().values;
    const shown: Record<number, Record<string, string>> = {};
    for (const name of Object.keys(values)) {
      const m = /^text of grid cell (\d+) of (.+)$/.exec(name);
      if (m == null)
        continue;
      (shown[Number(m[1])] ??= {})[m[2]] = String(values[name] ?? '');
    }
    const out: string[] = [];
    const seen = new Set<string>();
    for (const r of Object.keys(shown).map(Number).sort((x, y) => x - y)) {
      const cells = shown[r];
      const key = cells[keyCol];
      if (key === undefined) {
        out.push(`row ${r} shows no "${keyCol}" cell`);
        continue;
      }
      const group = expected[key];
      if (group === undefined) {
        out.push(`row ${r} is "${key}", which the reference aggregation does not have`);
        continue;
      }
      seen.add(key);
      for (const name of Object.keys(cells)) {
        if (name === keyCol)
          continue;
        if (!(name in group)) {
          out.push(`the reference aggregation has no "${name}" column`);
          continue;
        }
        const want = group[name];
        const text = cells[name];
        if (text === '') {
          if (want != null)
            out.push(`${key} / ${name}: the pivot shows nothing, the reference has ${want}`);
          continue;
        }
        const value = Number(text);
        if (want == null) {
          out.push(`${key} / ${name}: the pivot shows ${text}, the reference has nothing`);
          continue;
        }
        const digits = (/\.(\d+)$/.exec(text) ?? ['', ''])[1].length;
        const tolerance = Math.pow(10, -digits) / 2 + Math.abs(want) * 1e-9;
        if (Math.abs(value - want) > tolerance)
          out.push(`${key} / ${name}: the pivot shows ${text}, the reference has ${want}`);
      }
    }
    for (const key of Object.keys(expected))
      if (!seen.has(key))
        out.push(`the group "${key}" of the reference aggregation is not on the grid`);
    return out;
  }, [parsed[1], parsed[2], keyColumn, pivotColumn, filteredOnly ? 'yes' : 'no'] as [string, string, string, string, string]);
  expect(problems, `${target.phrase} does not show ${aggregation} grouped by ${keyColumn}`
    + (pivotColumn === '' ? '' : ` pivoted on ${pivotColumn}`) + (filteredOnly ? ' over the filtered rows' : '')).toEqual([]);
}

export const aggregationMatches = Then('the aggregated values of {widget} should match {string} grouped by {string}',
  (page: Page, target: ElementRef, aggregation: string, key: string) => expectAggregation(page, target, aggregation, key, '', false),
  {description: 'every cell the pivot draws equals an independent groupBy of the source table'});

export const pivotedAggregationMatches = Then('the aggregated values of {widget} should match {string} grouped by {string} pivoted on {string}',
  (page: Page, target: ElementRef, aggregation: string, key: string, pivot: string) => expectAggregation(page, target, aggregation, key, pivot, false),
  {description: 'the cross tab against an independent groupBy().pivot() — column by column, group by group'});

export const filteredAggregationMatches = Then('the aggregated values of {widget} should match {string} grouped by {string} over the filtered rows',
  (page: Page, target: ElementRef, aggregation: string, key: string) => expectAggregation(page, target, aggregation, key, '', true),
  {description: 'the same, over the rows the source filter passes — what Row Source = Filtered aggregates'});

// --- the command bar's history ---------------------------------------------------------------------

/** The history icon opens a menu of its own (`pivot-history`): "Save parameters" and one entry per
 * saved configuration whose columns the table still has. */
export const pickFromHistory = When('user picks {string} from the history menu of pivot table viewer',
  async (page: Page, item: string) => {
    const target = el('pivot table viewer');
    const c = viewers.centerOf(await viewers.hitArea(page, target, 'history', true));
    await page.mouse.click(c.x, c.y);
    const menu = page.locator('.d4-menu-popup[name="pivot-history"]');
    await menu.waitFor({state: 'visible', timeout: 5000});
    const entry = menu.locator('[role="menuitem"]').filter({hasText: exactText(item)}).first();
    if (await entry.count() === 0)
      throw new Error(`no "${item}" in the history menu; it offers: ${(await menu.locator('[role="menuitem"]').allTextContents()).join(', ')}`);
    await entry.click();
    await settle(page, target);
  }, {tier: 'ui', description: 'the command bar\'s history icon and one of its entries'});

export const clearSavedParameters = Given('user clears the saved pivot table parameters',
  (page: Page) => page.evaluate((key) => {
    window.localStorage.removeItem(key);
  }, HISTORY_KEY), {tier: 'api', description: 'the saved configurations live in localStorage and outlive a feature'});

// --- a checkbox no kind covers ----------------------------------------------------------------------

/** The select-all checkbox of the grid's "Order or Hide Columns" dialog is a bare
 * `input[type="checkbox"]` with no class and no label, which the library's `checkbox` kind (u2 bool
 * inputs and `.ui-input-bool`) does not see — the kind should learn about it, and then this step
 * goes away. */
export const clickPlainCheckbox = When('user clicks the plain checkbox in the {string} dialog',
  async (page: Page, title: string) => {
    const dialog = page.locator('.d4-dialog').filter({has: page.locator('.d4-dialog-title', {hasText: exactText(title)})}).last();
    await dialog.waitFor({state: 'visible', timeout: 5000});
    const box = dialog.locator('input[type="checkbox"]').filter({visible: true}).first();
    await expect(box, `a checkbox in the "${title}" dialog`).toBeVisible({timeout: 5000});
    await box.click();
  }, {tier: 'ui', description: 'a bare checkbox of a Dart dialog — should be promoted to the library as part of the checkbox kind'});

// --- the in-cell viewer columns ----------------------------------------------------------------------

/** The combo on the Aggregate row's title adds an in-cell viewer column per pivot category. */
export const addViewerColumn = When('user adds a {string} viewer column to pivot table viewer',
  async (page: Page, viewerName: string) => {
    const target = el('pivot table viewer');
    const c = viewers.centerOf(await viewers.hitArea(page, target, 'viewer selector', true));
    await page.mouse.click(c.x, c.y);
    const popup = page.locator('.d4-combo-popup-expanded').last();
    await popup.waitFor({state: 'visible', timeout: 5000});
    const item = popup.locator('.d4-list-item').filter({hasText: exactText(viewerName)}).first();
    if (await item.count() === 0)
      throw new Error(`no "${viewerName}" in the viewer selector; it offers: ${(await popup.locator('.d4-list-item').allTextContents()).join(', ')}`);
    await item.click();
    await settle(page, target);
  }, {tier: 'ui', description: 'the viewer picker of the Aggregate row — one viewer column per pivot category'});
