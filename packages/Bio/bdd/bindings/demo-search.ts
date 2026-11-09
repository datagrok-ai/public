/* Steps of the search and analysis demo clips (features/demos/15-19): the grid a search viewer keeps inside
   it, sized off camera, the main grid scrolled to the rows a clip is about, and what a diversity subset spans. */
import type {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {expect, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** The similarity and diversity viewers list their sequences in a grid of their own; the similarity viewer
 * asks for a 450 px column on its first render, which the grid's first layout narrows again, so the last
 * residues are cut until the next render. The width is set here, off camera, on the viewer of that type in
 * the current view. */
export const viewerGridColumnWidth = When('the {string} column of the {string} viewer is {int} px wide',
  async (page: Page, column: string, type: string, width: number) => {
    silent(page);
    await page.evaluate(([c, t, w]) => {
      const viewer = Array.from(grok.shell.tv?.viewers ?? []).find((v: any) => v.type === t) as any;
      if (!viewer)
        throw new Error(`the current view has no "${t}" viewer`);
      const host = viewer.root.querySelector('[name="viewer-Grid"]');
      const grid = host ? DG.Widget.find(host) : null;
      const col = grid?.col(c);
      if (!col)
        throw new Error(`the "${t}" viewer shows no grid with a "${c}" column`);
      col.width = w;
      grid.invalidate();
    }, [column, type, width] as [string, string, number]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a filmed demo leaves out: the width of a column in the grid a viewer shows inside it'});

/** The rows a clip selects or compares may lie below the first screen of the grid: the grid of the current view
 * starts at that table row, off camera, so they are in view when the clip acts on them. */
export const gridScrolledToRow = When('the grid is scrolled down to row {int}', async (page: Page, row: number) => {
  silent(page);
  await page.evaluate((r) => {
    const grid = grok.shell.tv?.grid;
    if (!grid || r < 1 || r > grid.dataFrame.rowCount)
      throw new Error(`the grid of the current view has no row ${r}`);
    grid.vertScroll.scrollTo(grid.tableRowToGrid(r - 1));
  }, row);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the grid\'s first visible row is that table row (from 1)'});

/** The rows a diversity search picked (its `subset` reading, table rows from 0) over the values a column of the
 * table holds in them: a subset that spans the table's families reads every one of them. */
export const subsetSpans = Then('the subset of {widget} should span {int} values of the {string} column',
  async (page: Page, target: ElementRef, count: number, column: string) => {
    await expect.poll(async () => {
      const subset = await viewers.readingOf(page, target, 'subset');
      if (typeof subset !== 'string' || subset === '')
        return String(subset);
      return page.evaluate(([ids, c]) => {
        const col = grok.shell.t?.col(c);
        if (!col)
          return `the current table has no "${c}" column`;
        return new Set(ids.split(',').map((i) => String(col.get(Number(i))))).size;
      }, [subset, column] as [string, string]);
    }, {message: `the "${column}" values of the subset of ${target.phrase}`}).toBe(count);
  }, {description: 'the distinct values of the column over the rows of the viewer\'s "subset" reading'});
