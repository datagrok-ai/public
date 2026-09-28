/* What only the derived-table, integration, augment and complex Projects features need. */
import type {Page} from '@playwright/test';
import {element, When} from '@datagrok-libraries/bdd';
import {expect, gestures, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The Join Tables dialog (Data > Join Tables...) names its two table choices and the key column
   selectors of its first key row, but not through an input host a phrase reaches: the tables share
   one "Tables" host, the keys one "Key Columns" host. */
const JOIN = '[name="dialog-Join-Tables"]';
element('join left table selector', {selector: `${JOIN} select[name="input-selectTableLeft"]`});
element('join right table selector', {selector: `${JOIN} select[name="input-selectTableRight"]`});
element('join left key selector', {selector: `${JOIN} [name="div-selectKeyCol0Row1"]`,
  description: 'the key column of the left table, first key row: a Dart column selector'});
element('join right key selector', {selector: `${JOIN} [name="div-selectKeyCol1Row1"]`,
  description: 'the key column of the right table, first key row'});

/* The table name at the left of the status bar: a click makes the current table the current object,
   so the context panel shows its panes (Actions: Clone, Rename...). Its name carries the table's
   size, which differs between stands. */
element('status bar table name', {selector: '.layout-status-bar .d4-status-bar-panel.d4-toolbar-items-text'});

/* A click in the grid of the current table view makes a cell current, which releases the object a
   rename left current in the context panel. The view's own grid is found through the view: with a
   pivot table on a view, the first grid in the page can be the pivot's. */
export const clickCurrentGrid = When('user clicks in the grid of the {string} table view', async (page: Page, table: string) => {
  const at = await page.evaluate((t) => {
    const tv = grok.shell.tv;
    if (tv == null || tv.dataFrame.name !== t)
      return `the current table view shows "${tv?.dataFrame.name ?? 'nothing'}"`;
    const r = (tv.grid.root as HTMLElement).getBoundingClientRect();
    return r.width > 0 && r.height > 0 ? {x: r.x + Math.min(r.width / 2, 200), y: r.y + Math.min(r.height / 2, 60)} : `the grid of "${t}" is not shown`;
  }, table);
  if (typeof at === 'string')
    throw new Error(at);
  await page.mouse.click(at.x, at.y);
}, {tier: 'ui', description: 'a click in the upper rows of that view\'s own grid; fails when the current table view shows another table'});

/* Data > Aggregate Rows... puts a second pivot table next to the one the toolbox added, so the tag
   rows of either are reached by an ordinal ("second pivot table viewer"); the pivot's own binding
   takes only the first. The same gesture: the + of the row, the pointer off it, the column picked
   (Aggregate Rows docks the Columns pane, a column grid too: `pickColumnCounted`). */
export const addToPivotRow = When('user adds {string} to the {string} row of {widget}',
  async (page: Page, column: string, row: string, target: ElementRef) => {
    const c = viewers.centerOf(await viewers.hitArea(page, target, `add ${row}`, true));
    // the picker opens its search on a keydown at the + icon it was opened from (pivot_grid.dart
    // passes the icon as its origin); a key sent to the page can land before the icon has focus
    const plus = (await viewers.viewerLocator(page, target))
      .locator(`[name="div-add-${row.replace(/ /g, '-')}" i] [tabindex]`).first();
    await gestures.pickColumnCounted(page, async () => {
      await page.mouse.click(c.x, c.y);
      await page.mouse.move(2, 2);
    }, column, `the "${row}" row`, plus);
    await expect.poll(async () => String(await viewers.readingOf(page, target, row)),
      {message: `the "${row}" reading of ${target.phrase} after "${column}" was picked`, timeout: pollMs(10000)}).toContain(column);
  }, {tier: 'ui', description: 'the + of a tag row of that pivot table, the column typed and committed; done when the row\'s reading names the column'});
