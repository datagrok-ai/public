/* What only the Queries features need: the walk over every column of a schema, and the visual
   query builder, which publishes no status of its own.
   Everything generic these features use (queries on the server, the editor's code, the current
   view's type) is the library's. */
import {Page} from '@playwright/test';
import {Given, Then, When, kind} from '@datagrok-libraries/bdd';
import {exactText, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* --- columns of a schema ----------------------------------------------------------------------
   The case clicks every column of every table and watches the context panel follow. The tree is
   walked with real clicks; each column click is done when the current object (grok.shell.o) is
   that column and the panel's accordion has been constructed — the panes a feature adds on that
   event (Database meta) come a task after the column is current — and the panes the panel then
   shows are recorded for the claim that follows. */

interface InspectedColumn {table: string; column: string; shown: string; panes: string[]}
const inspected = new WeakMap<Page, InspectedColumn[]>();

const norm = (s: string) => s.toLowerCase().replace(/[^a-z0-9]/g, '');

const directChildren = (page: Page, p: string) => page.evaluate((pre) => Array.from(document.querySelectorAll(`[name^="${pre}---"]`))
  .map((e) => e.getAttribute('name')!).filter((n) => !n.slice(pre.length + 3).includes('---')), p);

/** Expands each table node in turn, clicks each of its columns and records what the context panel shows. */
async function clickColumnsOf(page: Page, tables: string[], prefix: string): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    w.__bddPanelAccordions = 0;
    w.__bddPanelAccordionSub?.unsubscribe();
    w.__bddPanelAccordionSub = grok.events.onAccordionConstructed.subscribe((a: any) => {
      if (a?.root?.closest?.('.grok-prop-panel'))
        w.__bddPanelAccordions++;
    });
  });
  const constructed = () => page.evaluate(() => (window as any).__bddPanelAccordions as number);
  const results: InspectedColumn[] = [];
  for (const table of tables) {
    const node = page.locator(`[name="${table}"]`).first();
    await node.scrollIntoViewIfNeeded();
    const tri = node.locator('.d4-tree-view-tri').first();
    if (!await tri.evaluate((e) => e.classList.contains('d4-tree-view-tri-expanded')))
      await tri.click();
    await expect.poll(async () => (await directChildren(page, table)).length, {message: `the columns listed under ${table}`, timeout: pollMs(15000)}).toBeGreaterThan(0);
    for (const col of [...new Set(await directChildren(page, table))]) {
      const colNode = page.locator(`[name="${col}"]`).first();
      await colNode.scrollIntoViewIfNeeded();
      const want = col.slice(table.length + 3);
      let shown = '';
      const before = await constructed();
      // the panel drops a change within 2 s of a property edit and while the object is frozen
      // (property_panel.dart), so a click that lands in that window makes nothing current: it is
      // made again before the column counts as one the panel did not follow
      for (let attempt = 0; attempt < 2 && norm(shown) !== norm(want); attempt++) {
        await colNode.click({position: {x: 60, y: 8}});
        await expect.poll(async () => (shown = await page.evaluate(() => String(grok.shell.o?.name ?? ''))) && norm(shown) === norm(want),
          {message: `the current object after clicking ${col}`, timeout: pollMs(attempt === 0 ? 5000 : 10000)}).toBe(true).catch(() => undefined);
      }
      if (norm(shown) === norm(want))
        await expect.poll(constructed, {message: `the context panel's accordion for ${col}`, timeout: pollMs(10000)}).toBeGreaterThan(before);
      const panes: string[] = await page.evaluate(() => Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
        .map((e) => (e.firstChild?.textContent ?? '').trim()));
      results.push({table: table.slice(prefix.length + 3), column: want, shown, panes});
    }
    await tri.click();
  }
  await page.evaluate(() => (window as any).__bddPanelAccordionSub?.unsubscribe());
  inspected.set(page, results);
}

// package: the columns-inspect loop (the survey's step 11)
export const clickEveryColumn = When('user clicks every column of every table under {element}', async (page: Page, target: ElementRef) => {
  const prefix = String(await (await locate(page, target)).first().getAttribute('name'));
  const tables = [...new Set(await directChildren(page, prefix))];
  expect(tables.length, `tables under ${target.phrase}`).toBeGreaterThan(0);
  await clickColumnsOf(page, tables, prefix);
}, {tier: 'ui', description: 'expands each table of the schema node in turn, clicks each of its columns and records what the context panel then shows'});

// package
export const clickEveryColumnOfTable = When('user clicks every column of {element}', async (page: Page, target: ElementRef) => {
  const table = String(await (await locate(page, target)).first().getAttribute('name'));
  await clickColumnsOf(page, [table], table.slice(0, table.lastIndexOf('---')));
}, {tier: 'ui', description: 'expands the table node, clicks each of its columns and records what the context panel then shows'});

// package
export const everyColumnShown = Then('every clicked column should have been shown in the context panel with {string}', async (page: Page, list: string) => {
  const results = inspected.get(page) ?? [];
  expect(results.length, 'columns clicked by "user clicks every column of …"').toBeGreaterThan(0);
  const want = list.split(',').map((s) => s.trim());
  const bad = results.filter((r) => norm(r.shown) !== norm(r.column) || want.some((p) => !r.panes.includes(p)))
    .map((r) => `${r.table}.${r.column}: panel shows "${r.shown}" with ${r.panes.join(', ') || 'no panes'}`);
  expect(bad, `columns of ${results.length} clicked whose context panel did not follow`).toEqual([]);
}, {description: 'each clicked column became the current object and its panel listed every named pane'});

// package: the actions browser renders a function as its markup, which carries no name of its own
export const openAction = When('user opens the {string} action of the transformations browser', async (page: Page, name: string) => {
  await page.locator('.grok-actions-browser').filter({visible: true}).getByText(name, {exact: true}).first().click();
}, {tier: 'ui', description: 'clicks the function in the Transformations tab\'s actions browser, which opens its dialog'});

/* --- the visual query builder ------------------------------------------------------------------
   Each row of the builder (Data, Where, Group by, Aggregate, Pivot, Having, Order by) has a `+`
   named `div-add-<Row>`; the picker it opens is the platform's column grid. The tags it adds carry
   the column name (`span-tag-rows-<column>`). Package until the builder publishes a status of its
   own (the pivot viewer's `getWidgetStatus` is the model). */

const plusName = (row: string) => `div-add-${row.trim().replace(/\s+/g, '-')}`;

export const addToBuilderRow = When('user adds {string} to the {string} row of the visual query', async (page: Page, column: string, row: string) => {
  const plus = page.locator(`[name="${plusName(row)}"]`).filter({visible: true}).first();
  await expect(plus, `the + of the "${row}" row of the visual query`).toBeVisible({timeout: pollMs(15000)});
  await plus.click();
  await page.mouse.move(2, 2);
  await gestures.pickInColumnGrid(page, column, `the "${row}" row of the visual query`);
}, {tier: 'ui', description: 'clicks the + of the row and picks the column in the grid the platform opens'});

// the builder of the view shown: a closed editor's can stay in the page
const builderPanel = (page: Page, row: string) => page.locator('.grok-pivot-root').filter({visible: true}).first()
  .locator('.grok-pivot-column-panel').filter({has: page.locator(`[name="${plusName(row)}"]`)}).first();

/** The captions of a row's tags: a Where or Having tag also holds its condition and parameter inputs,
 * whose own text is no caption. */
const builderTags = (page: Page, row: string) => builderPanel(page, row).locator('.d4-tag').evaluateAll((tags) => tags.map((t) => {
  const copy = t.cloneNode(true) as Element;
  copy.querySelectorAll('.ui-input-root').forEach((i) => i.remove());
  return (copy.textContent ?? '').replace(/×|✕/g, '').trim();
}).filter(Boolean));

export const builderRowHolds = Then('the {string} row of the visual query should hold {string}', async (page: Page, row: string, list: string) => {
  const want = list.split(',').map((s) => s.trim()).filter(Boolean);
  await expect.poll(() => builderTags(page, row), {message: `the tags of the "${row}" row of the visual query`, timeout: pollMs(15000)}).toEqual(want);
}, {description: 'the tags the row shows, in order'});

export const removeFromBuilderRow = When('user removes {string} from the {string} row of the visual query', async (page: Page, tag: string, row: string) => {
  const tags = await builderTags(page, row);
  expect(tags, `the tags of the "${row}" row of the visual query`).toContain(tag);
  const chip = builderPanel(page, row).locator('.d4-tag').nth(tags.indexOf(tag));
  await chip.hover();
  await chip.locator('.d4-tag-delete').click();
}, {tier: 'ui', description: 'clicks the close icon of the tag'});

/* The Where and Having conditions of the builder are named by row and column (visual_query.dart:
   `input-where-condition-name`, `input-where-param-name`, a Having one by its aggregation as well,
   `input-having-condition-values-id`), and annotate() dashes the column's underscores. */
kind('where condition', {
  selector: 'input[name^="input-where-condition-"]',
  match: ['dart'],
  dartNames: ['input-where-condition-{q}'],
  description: 'the condition typed beside a column of the Where row of the visual query ("name" where condition)',
});
kind('having condition', {
  selector: 'input[name^="input-having-condition-"]',
  match: ['dart'],
  dartNames: ['input-having-condition-{q}'],
  description: 'the condition beside an aggregation of the Having row, by "<aggregation>-<column>" ("values-id" having condition)',
});
kind('where parameter checkbox', {
  selector: 'input[name^="input-where-param-"]',
  match: ['dart'],
  dartNames: ['input-where-param-{q}'],
  description: 'the "Expose as function parameter" box of a Where condition, by its column',
});

export const builderRowShown = Then('the {string} row of the visual query should be shown', async (page: Page, row: string) => {
  await expect(builderPanel(page, row), `the "${row}" row of the visual query`).toBeVisible({timeout: pollMs(15000)});
});

export const builderRowHidden = Then('the {string} row of the visual query should be hidden', async (page: Page, row: string) => {
  await expect(builderPanel(page, row), `the "${row}" row of the visual query`).toBeHidden({timeout: pollMs(15000)});
}, {description: 'the Pivot and Having rows wait for an aggregation'});

export const builderRowOffers = Then('the {string} row of the visual query should offer the columns {string}', async (page: Page, row: string, list: string) => {
  const want = list.split(',').map((s) => s.trim()).filter(Boolean);
  const plus = page.locator(`[name="${plusName(row)}"]`).filter({visible: true}).first();
  await plus.click();
  await page.mouse.move(2, 2);
  const popup = page.locator('.d4-column-grid:has(.d4-column-selector-backdrop)').filter({visible: true}).last();
  await popup.waitFor({state: 'visible', timeout: pollMs(10000)});
  const offered = await gestures.columnGridNames(popup.locator('[name="viewer-Grid"]').first());
  await page.keyboard.press('Escape');
  await popup.waitFor({state: 'detached', timeout: pollMs(5000)});
  expect(offered, `the columns the + of the "${row}" row of the visual query offers`).toEqual(want);
}, {tier: 'ui', description: 'opens the picker of the row\'s +, reads every column it offers in its order and dismisses it'});

export const clickBuilderPlus = When('user clicks on the + of the {string} row of the visual query', async (page: Page, row: string) => {
  const plus = page.locator(`[name="${plusName(row)}"] i`).filter({visible: true}).first();
  await expect(plus, `the + of the "${row}" row of the visual query`).toBeVisible({timeout: pollMs(15000)});
  await plus.click();
  await page.mouse.move(2, 2);
}, {tier: 'ui', description: 'the + of the row, or the Data row\'s icon that adds a table to join; what it opens is the open menu'});

/* With a table joined, a row's + opens a menu of the tables, each a group whose submenu is a column
   grid the keyboard does not reach (the focus stays on the +): the column's row is clicked where the
   grid reports it drew the name, its "text of cell N of __name" reading and "cell N of __name" area —
   in a submenu the grid's `cell().bounds` put its rows 16 px high against the 28 px it draws. */
export const pickJoinedColumn = When('user picks the {string} column of the {string} table from the open menu', async (page: Page, column: string, table: string) => {
  const group = page.locator('.d4-menu-item').filter({visible: true})
    .filter({has: page.locator(':scope > .d4-menu-item-label', {hasText: exactText(table)})}).last();
  await expect(group, `the "${table}" group of the open menu`).toBeVisible({timeout: pollMs(10000)});
  await group.hover();
  const grid = page.locator('.d4-menu-popup [name="viewer-Grid"]').filter({visible: true}).last();
  await expect(grid, `the columns of "${table}" in the open menu`).toBeVisible({timeout: pollMs(10000)});
  await viewers.settleAll(page);
  let box: {x: number; y: number} | string = '';
  await expect.poll(async () => {
    box = await grid.evaluate((el, n) => {
      const w = window as any;
      const values = w.__bdd.viewerOf(el).getWidgetStatus()?.values ?? {};
      const key = Object.keys(values).find((k) => /^text of cell \d+ of __name$/.test(k) && values[k] === n);
      if (!key)
        return Object.keys(values).filter((k) => /^text of cell \d+ of __name$/.test(k)).map((k) => values[k]).join(', ');
      const area = w.__bdd.areas(el)[key.replace(/^text of /, '')];
      return area ? {x: area.x + area.width / 2, y: area.y + area.height / 2} : 'no area for ' + key;
    }, column);
    return typeof box !== 'string';
  }, {message: `the "${column}" row in the columns of "${table}"`, timeout: pollMs(5000)}).toBe(true).catch(() => {
    throw new Error(`the "${table}" table offers no "${column}" row; the grid shows: ${box}`);
  });
  const at = box as unknown as {x: number; y: number};
  // into the submenu level with the group first, so the pointer crosses no other item of the menu
  const from = (await group.boundingBox())!;
  await page.mouse.move(at.x, from.y + from.height / 2);
  await page.mouse.move(at.x, at.y);
  await page.mouse.click(at.x, at.y);
  await expect(grid, `the menu after "${column}" was picked`).toBeHidden({timeout: pollMs(5000)});
}, {tier: 'ui', description: 'hovers the table\'s group of the open menu and clicks the column\'s row in the column grid it opens; the menu closes on the pick'});

export const builderPreviewRows = Then('the visual query should preview {int} row(s)', async (page: Page, rows: number) => {
  const counts = page.locator('.grok-pivot-counts').filter({visible: true}).first();
  await expect.poll(async () => /(\d+) rows/.exec(await counts.innerText().catch(() => ''))?.[1] ?? 'nothing',
    {message: 'the row count over the visual query\'s preview', timeout: pollMs(15000)}).toBe(String(rows));
}, {description: 'the "N rows" the builder shows over its preview grid'});
