/* What only the Queries features need: the walk over every column of a schema, and the visual
   query builder, which publishes no status of its own.
   Everything generic these features use (queries on the server, the editor's code, the current
   view's type) is the library's. */
import {expect, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* --- columns of a schema ----------------------------------------------------------------------
   The case clicks every column of every table and watches the context panel follow. The tree is
   walked with real clicks; each column click is done when the current object (grok.shell.o) is
   that column, and the panes the panel then shows are recorded for the claim that follows. */

interface InspectedColumn {table: string; column: string; shown: string; panes: string[]}
const inspected = new WeakMap<Page, InspectedColumn[]>();

const norm = (s: string) => s.toLowerCase().replace(/[^a-z0-9]/g, '');

// package: the columns-inspect loop (the survey's step 11)
export const clickEveryColumn = When('user clicks every column of every table under {element}', async (page: Page, target: ElementRef) => {
  const schema = (await locate(page, target)).first();
  const prefix = String(await schema.getAttribute('name'));
  const directChildren = (p: string) => page.evaluate((pre) => Array.from(document.querySelectorAll(`[name^="${pre}---"]`))
    .map((e) => e.getAttribute('name')!).filter((n) => !n.slice(pre.length + 3).includes('---')), p);
  const tables = [...new Set(await directChildren(prefix))];
  expect(tables.length, `tables under ${target.phrase}`).toBeGreaterThan(0);
  const results: InspectedColumn[] = [];
  for (const table of tables) {
    const node = page.locator(`[name="${table}"]`).first();
    await node.scrollIntoViewIfNeeded();
    const tri = node.locator('.d4-tree-view-tri').first();
    if (!await tri.evaluate((e) => e.classList.contains('d4-tree-view-tri-expanded')))
      await tri.click();
    await expect.poll(async () => (await directChildren(table)).length, {message: `the columns listed under ${table}`, timeout: pollMs(15000)}).toBeGreaterThan(0);
    for (const col of [...new Set(await directChildren(table))]) {
      const colNode = page.locator(`[name="${col}"]`).first();
      await colNode.scrollIntoViewIfNeeded();
      await colNode.click({position: {x: 60, y: 8}});
      const want = col.slice(table.length + 3);
      let shown = '';
      await expect.poll(async () => (shown = await page.evaluate(() => String(grok.shell.o?.name ?? ''))) && norm(shown) === norm(want),
        {message: `the current object after clicking ${col}`, timeout: pollMs(10000)}).toBe(true).catch(() => undefined);
      const panes: string[] = await page.evaluate(() => Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
        .map((e) => (e.firstChild?.textContent ?? '').trim()));
      results.push({table: table.slice(prefix.length + 3), column: want, shown, panes});
    }
    await tri.click();
  }
  inspected.set(page, results);
}, {tier: 'ui', description: 'expands each table of the schema node in turn, clicks each of its columns and records what the context panel then shows'});

// package
export const everyColumnShown = Then('every clicked column should have been shown in the context panel with {string}', async (page: Page, list: string) => {
  const results = inspected.get(page) ?? [];
  expect(results.length, 'columns clicked by "user clicks every column of every table under …"').toBeGreaterThan(0);
  const want = list.split(',').map((s) => s.trim());
  const bad = results.filter((r) => norm(r.shown) !== norm(r.column) || want.some((p) => !r.panes.includes(p)))
    .map((r) => `${r.table}.${r.column}: panel shows "${r.shown}" with ${r.panes.join(', ') || 'no panes'}`);
  expect(bad, `columns of ${results.length} clicked whose context panel did not follow`).toEqual([]);
}, {description: 'each clicked column became the current object and its panel listed every named pane'});

// package, temporary (CORE-SIGNAL): the actions browser of the Transformations tab takes about a
// second after the tab opens before a click on a function runs it — the first click only makes
// the row current — and nothing in the DOM says when it is ready. The click is repeated until the
// function's dialog is up; replace with a plain click once the browser publishes a ready state.
export const openAction = When('user opens the {string} action of the transformations browser', async (page: Page, name: string) => {
  const link = page.locator('.grok-actions-browser').filter({visible: true}).getByText(name, {exact: true}).first();
  const dialog = page.locator(`[name="dialog-${name.trim().replace(/\s+/g, '-')}"]`).filter({visible: true});
  await expect.poll(async () => {
    if (await dialog.count() === 0)
      await link.click();
    return dialog.count();
  }, {message: `the "${name}" dialog after clicking its action`, timeout: pollMs(10000), intervals: [300]}).toBeGreaterThan(0);
}, {tier: 'ui', description: 'clicks the function in the Transformations tab\'s actions browser until its dialog opens (the browser has no ready signal)'});

/* --- the visual query builder ------------------------------------------------------------------
   Each row of the builder (Data, Where, Group by, Aggregate, Pivot, Having, Order by) has a `+`
   named `div-add-<Row>`; the picker it opens is the platform's column grid. The tags it adds carry
   the column name (`span-tag-rows-<column>`). Package until the builder publishes a status of its
   own (the pivot viewer's `getWidgetStatus` is the model). */

// package, temporary (CORE-SIGNAL: a status for the visual query editor)
export const addToBuilderRow = When('user adds {string} to the {string} row of the visual query', async (page: Page, column: string, row: string) => {
  const plus = page.locator(`[name="div-add-${row.trim().replace(/\s+/g, '-')}"]`).filter({visible: true}).first();
  await expect(plus, `the + of the "${row}" row of the visual query`).toBeVisible({timeout: pollMs(15000)});
  await plus.click();
  await page.mouse.move(2, 2);
  await gestures.pickInColumnGrid(page, column, `the "${row}" row of the visual query`);
}, {tier: 'ui', description: 'clicks the + of the row and picks the column in the grid the platform opens'});

// package, temporary (same CORE-SIGNAL)
export const builderRowHolds = Then('the {string} row of the visual query should hold {string}', async (page: Page, row: string, list: string) => {
  const want = list.split(',').map((s) => s.trim()).filter(Boolean);
  const panel = page.locator('.grok-pivot-column-panel').filter({visible: true})
    .filter({has: page.locator(`[name="div-add-${row.trim().replace(/\s+/g, '-')}"]`)}).first();
  await expect.poll(async () => (await panel.locator('.d4-tag').allTextContents()).map((t) => t.replace(/×|✕/g, '').trim()).filter(Boolean),
    {message: `the tags of the "${row}" row of the visual query`, timeout: pollMs(15000)}).toEqual(want);
}, {description: 'the tags the row shows, in order'});
