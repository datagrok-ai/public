/* Claims about the user's table while another view is current: the Train view is not a table view, so the library's
   "the table" (grok.shell.t) is not the table the view trains on; its mouse-over row against a viewer's hovered row,
   and a scatter plot's point hovered by its row. */
import type {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, viewers} from '@datagrok-libraries/bdd/runtime';
import {listOf} from './grok-s.js';
import {slidePointer} from './pointer.js';

declare const DG: any;

/** The named open table's mouse-over row (0-based, -1 for none), or why there is none to read. */
function mouseOverRow(page: Page, name: string): Promise<number> {
  return viewers.evaluate(page, (n) => (window as any).__bdd.tableNamed(n).mouseOverRowIdx as number, name);
}

/* LIB-PROMOTE: the library's mouse-over claims read grok.shell.t. A baseline for the hover claim below: a row left
   standing by an earlier gesture would pass it. */
export const noMouseOverRow = Then('the table {string} should have no mouse-over row', async (page: Page, name: string) => {
  expect(await mouseOverRow(page, name), `the mouse-over row of "${name}"`).toBe(-1);
}, {description: 'the named open table\'s mouseOverRowIdx is -1, read once'});

/* LIB-PROMOTE: a hovered point of a dense plot can lie under another one, so which row the gesture aimed at is no
   claim; that the user's table has the row the plot says is hovered is (the plot is over the user's columns). */
export const mouseOverRowIsHovered = Then('the mouse-over row of the table {string} should be the {string} reading of {widget}',
  async (page: Page, name: string, reading: string, target: ElementRef) => {
    let hovered: unknown;
    await expect.poll(async () => {
      hovered = await viewers.readingOf(page, target, reading);
      const row = await mouseOverRow(page, name);
      return typeof hovered === 'number' && hovered > 0 && row + 1 === hovered;
    }, {message: `the mouse-over row of "${name}" (1-based) against "${reading}" of ${target.phrase} (${String(hovered)})`})
      .toBe(true);
  }, {description: 'the table\'s mouseOverRowIdx + 1 equals the viewer\'s 1-based reading, which is not 0'});

/* A viewer of the results draws the rows a training kept, a copy without the skipped rows: its hovered row is not
   the table's row number once a row was skipped before it, but it holds the same values. */
export const mouseOverRowHoldsHovered = Then('the mouse-over row of the table {string} should hold the {string} values ' +
  'of the hovered row of {widget}', async (page: Page, name: string, list: string, target: ElementRef) => {
  await expect.poll(() => viewers.onViewer(page, target, (el, [n, columns]: [string, string[]]) => {
    const b = (window as any).__bdd;
    const frame = DG.Widget.find(el).dataFrame;
    const table = b.tableNamed(n);
    const [hovered, row] = [frame.mouseOverRowIdx as number, table.mouseOverRowIdx as number];
    if (hovered < 0 || row < 0)
      return `the viewer's hovered row ${hovered + 1}, the table's mouse-over row ${row + 1}`;
    const differ = columns.filter((c) => String(b.col(c, frame).get(hovered)) !== String(b.col(c, table).get(row)));
    return differ.length === 0 ? 'same' : `row ${row + 1} of "${n}" and the hovered row ${hovered + 1} differ in ` +
      differ.join(', ');
  }, [name, listOf(list)]), {message: `the mouse-over row of "${name}" against the hovered row ` +
    `of ${target.phrase}`}).toBe('same');
}, {description: 'both have a hovered row, and the table\'s row holds the viewer\'s values in each listed column ' +
  '(comma-separated), compared as text'});

/* LIB-PROMOTE: a scatter plot reports `marker of row <n>` areas for its first rows only; a row further down (past the
   first skipped row of a training) is reached through the plot's own pointToScreen. With the library's `marker of row`
   areas falling back to pointToScreen for the rows past the listed ones, this step goes. */
export const hoverPointOfRow = When('user hovers over the point of row {int} of {widget}', async (page: Page, row: number,
  target: ElementRef) => {
  const at = await viewers.onViewer(page, target, (el, r: number) => {
    const point = DG.Widget.find(el).pointToScreen(r - 1);
    const canvas = el.querySelector('canvas');
    if (canvas === null)
      throw new Error('the viewer has no canvas');
    const box = canvas.getBoundingClientRect();
    return {x: box.left + point.x, y: box.top + point.y};
  }, row);
  await slidePointer(page, at, {x: at.x - 20, y: at.y});
}, {tier: 'ui', description: 'the scatter plot\'s point of the row (1-based, of the frame it draws), approached from ' +
  'the left in five moves; a point drawn over it may be the one hovered'});
