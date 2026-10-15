/* Claims about the user's table while another view is current: the Train view is not a table view, so the library's
   "the table" (grok.shell.t) is not the table the view trains on. */
import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The named open table's mouse-over row (0-based, -1 for none), or why there is none to read. */
function mouseOverRow(page: Page, name: string): Promise<number | string> {
  return page.evaluate((n) => {
    const table = Array.from(grok.shell.tables as any[]).find((t: any) => t.name === n);
    return table === undefined ? `no open table "${n}"` : table.mouseOverRowIdx as number;
  }, name);
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
      return typeof row === 'number' && typeof hovered === 'number' && hovered > 0 && row + 1 === hovered;
    }, {message: `the mouse-over row of "${name}" (1-based) against "${reading}" of ${target.phrase} (${String(hovered)})`})
      .toBe(true);
  }, {description: 'the table\'s mouseOverRowIdx + 1 equals the viewer\'s 1-based reading, which is not 0'});
