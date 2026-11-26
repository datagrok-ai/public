/* Steps the "Databases and queries" how-to clips (features/howto/12-databases) need and the library's vocabulary
   lacks: the visual query builder, which publishes no status of its own. The builder's steps follow UsageAnalysis'
   bdd/bindings/queries.ts (same names: each row's `+` is `div-add-<Row>`, its picker the platform's column grid),
   filmed as a person works it: the + lit before the press, the pointer resting beside it, clear of the picker. */
import type {Page} from '@playwright/test';
import {kind, Then, When} from '@datagrok-libraries/bdd';
import {expect, gestures, guideHop, pollMs} from '@datagrok-libraries/bdd/runtime';

/* The schema view names each table box's host (`div-table-<table>`, db_views.dart), which no library kind reaches:
   the same noun as UsageAnalysis' bdd/bindings/connections.ts (a library candidate). */
kind('schema table', {
  selector: '.d4-host',
  match: ['dart'],
  dartNames: ['div-table-{q}'],
  description: 'a table box of the schema view (Schemas > <schema> > Browse), by its table name',
});

const plusName = (row: string) => `div-add-${row.trim().replace(/\s+/g, '-')}`;

// the builder of the view shown: a closed editor's can stay in the page
const builderPanel = (page: Page, row: string) => page.locator('.grok-pivot-root').filter({visible: true}).first()
  .locator('.grok-pivot-column-panel').filter({has: page.locator(`[name="${plusName(row)}"]`)}).first();

/** The captions of a row's tags: a Where or Having tag also holds its condition inputs, whose text is no caption. */
const builderTags = (page: Page, row: string) => builderPanel(page, row).locator('.d4-tag').evaluateAll((tags) => tags.map((t) => {
  const copy = t.cloneNode(true) as Element;
  copy.querySelectorAll('.ui-input-root').forEach((i) => i.remove());
  return (copy.textContent ?? '').replace(/×|✕/g, '').trim();
}).filter(Boolean));

export const addToVisualQueryRow = When('user adds {string} to the {string} row of the visual query', async (page: Page, column: string, row: string) => {
  const plus = page.locator(`[name="${plusName(row)}"]`).filter({visible: true}).first();
  await expect(plus, `the + of the "${row}" row of the visual query`).toBeVisible({timeout: pollMs(15000)});
  await guideHop(page, plus);
  const box = (await plus.boundingBox())!;
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  // beside the +, to its left: the picker opens below and to the right of it
  await page.mouse.move(Math.max(2, box.x - 24), box.y + box.height / 2);
  await gestures.pickInColumnGrid(page, column, `the "${row}" row of the visual query`);
}, {tier: 'ui', description: 'clicks the + of the row and picks the column in the grid the platform opens (typed into its search)'});

/** The query editor's play icon runs the query, and the result gives the view a table: the ribbon is rebuilt with
 * the table's menus, and "Select" lands under the pointer still resting on the play icon — Chrome's move under a
 * still pointer then opens it over the editor. So the pointer leaves the icon at once, down into the editor's empty
 * part, before the result comes back. */
export const runWithPlayIcon = When('user runs the query with the play icon', async (page: Page) => {
  const play = page.locator('[name="icon-play"]').filter({visible: true}).first();
  await expect(play, 'the play icon of the query editor').toBeVisible({timeout: pollMs(15000)});
  await guideHop(page, play);
  const box = (await play.boundingBox())!;
  // the Layout tab shows no editor: the pointer then goes straight down from the icon
  const editor = await page.locator('.CodeMirror').filter({visible: true}).first().boundingBox({timeout: 500}).catch(() => null);
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  const rest = editor ? {x: editor.x + editor.width / 3, y: editor.y + Math.min(160, editor.height / 2)} :
    {x: box.x + box.width / 2, y: box.y + 160};
  await page.mouse.move(rest.x, rest.y, {steps: 4});
}, {tier: 'ui', description: 'clicks the play icon of the query editor, then rests the pointer in the editor below, off the ribbon the result rebuilds'});

export const visualQueryRowHolds =Then('the {string} row of the visual query should hold {string}', async (page: Page, row: string, list: string) => {
  const want = list.split(',').map((s) => s.trim()).filter(Boolean);
  await expect.poll(() => builderTags(page, row), {message: `the tags of the "${row}" row of the visual query`, timeout: pollMs(15000)}).toEqual(want);
}, {description: 'the tags the row shows, in order'});

export const visualQueryPreviewRows = Then('the visual query should preview {int} row(s)', async (page: Page, rows: number) => {
  const counts = page.locator('.grok-pivot-counts').filter({visible: true}).first();
  await expect.poll(async () => /(\d+) rows/.exec(await counts.innerText().catch(() => ''))?.[1] ?? 'nothing',
    {message: 'the row count over the visual query\'s preview', timeout: pollMs(15000)}).toBe(String(rows));
}, {description: 'the "N rows" the builder shows over its preview grid'});
