// Steps of the "every sequence is a molecule" demo clips (features/demos/05-07): the set-up a filmed demo leaves out
// (the context panel widened, sequences made molecules or cut to a region by Bio's own commands, the grid sorted), a
// click on one marker of a scatter plot of the molecules' properties and a nested context panel pane closed.
import type {Page} from '@playwright/test';
import {Given, When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {atFeatureEnd, el, expect, locate, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The context panel docks at a quarter of the window each time it is shown (`PropertyPanel.buildPanel`), and Bio's
 * 3D pane is as tall as the panel is wide: widened by its splitter for the clip, the panel is docked afresh at feature
 * end (hidden, then shown again), so a later feature meets the default width. */
export const contextPanelWider = Given('the context panel is {int} px wider', async (page: Page, px: number) => {
  silent(page);
  const panel = page.locator('.grok-prop-panel').filter({visible: true}).first();
  await expect(panel, 'the context panel').toBeVisible();
  const before = (await panel.boundingBox())!.width;
  const splitter = (await locate(page, el('context panel splitter'))).filter({visible: true}).first();
  await expect(splitter, 'the splitter of the context panel').toBeVisible();
  const box = (await splitter.boundingBox())!;
  await viewers.dragFrom(page, {x: box.x + box.width / 2, y: box.y + box.height / 2}, {dx: -px, dy: 0});
  atFeatureEnd(page, async () => {
    await page.evaluate(() => {
      if (!grok.shell.windows.showContextPanel)
        return;
      grok.shell.windows.showContextPanel = false;
      grok.shell.windows.showContextPanel = true;
    }).catch(() => undefined);
  });
  await expect.poll(async () => Math.round(((await panel.boundingBox())?.width ?? before) - before),
    {message: 'how much wider the context panel is'}).toBeGreaterThan(px / 2);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'set-up a filmed demo leaves out: the splitter left of the context panel dragged left by that many px; the panel is docked afresh at its default width at feature end'});

/** A clip about what Chem does with the molecules starts from them: Bio | Transform | To Atomic Level with the
 * dialog's defaults (Non-linear on: rings and bridges drawn as rings; no monomer highlight), run off camera. */
export const atomicLevelDefaults = Given('the {string} column is converted to atomic level', async (page: Page, column: string) => {
  silent(page);
  await page.evaluate(async (name) => {
    const t = grok.shell.t;
    if (!t?.col(name))
      throw new Error(`the current table has no "${name}" column`);
    await grok.functions.call('Bio:toAtomicLevel', {table: t, seqCol: t.col(name), nonlinear: true, highlight: false});
  }, column);
  await expect.poll(() => page.evaluate((name) => grok.shell.t.col(`molfile(${name})`) != null, column),
    {message: `the molfile(${column}) column`}).toBe(true);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'Bio:toAtomicLevel on the current table\'s column with Non-linear on and Highlight monomers off, the dialog\'s defaults'});

/** Bio | Calculate | Extract Region (`Bio:getRegionTopMenu`), run off camera: positions start..end of every sequence,
 * named as the sequence counts them (from 1), in a new sequence column of that name. */
export const regionExtracted = Given('the region {int} to {int} of the {string} column is extracted as {string}',
  async (page: Page, start: number, end: number, column: string, name: string) => {
    silent(page);
    await page.evaluate(async ([s, e, c, n]) => {
      const t = grok.shell.t;
      if (!t?.col(c))
        throw new Error(`the current table has no "${c}" column`);
      await grok.functions.call('Bio:getRegionTopMenu', {table: t, sequence: t.col(c), start: String(s), end: String(e), name: n});
    }, [start, end, column, name] as [number, number, string, string]);
    await expect.poll(() => page.evaluate((n) => grok.shell.t.col(n)?.semType ?? 'no such column', name),
      {message: `the "${name}" column`}).toBe('Macromolecule');
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a filmed demo leaves out: Bio\'s Extract Region on the current table, the region in a new Macromolecule column'});

/** The grid of the current view sorted on a column, as a click on its header sorts it, off camera: a clip that filters
 * a table shows a mix of its kinds of rows in the first screen, not the first family of the file. */
export const gridSortedBy = Given('the grid is sorted by the {string} column in {word} order', async (page: Page, column: string, order: string) => {
  silent(page);
  if (order !== 'ascending' && order !== 'descending')
    throw new Error(`a grid is sorted in ascending or descending order, not "${order}"`);
  await page.evaluate(([c, asc]) => {
    const grid = grok.shell.tv?.grid;
    if (!grid?.dataFrame.col(c))
      throw new Error(`the grid of the current view has no "${c}" column`);
    grid.sort([c], [asc]);
  }, [column, order === 'ascending'] as [string, boolean]);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the current view\'s grid sorted on the column, ascending or descending'});

/** A scatter plot reports a hit area for the markers of its first twenty rows only (`marker of row <n>`, in table
 * order): a marker further down the table is found through the plot's own coordinates (`worldToScreen` of the row's x
 * and y values, relative to its canvas). The pointer rests there until the plot says that row is under it and its
 * tooltip shows, then clicks: the row becomes the current one. */
export const clickMarker = When('user clicks on the marker of row {int} of {widget}', async (page: Page, row: number, target: ElementRef) => {
  await viewers.onViewer(page, target, (e) => (window as any).__bdd.settle(e, 300), undefined).catch(() => undefined);
  const point = await viewers.onViewer(page, target, (e, r) => {
    const sp = (window as any).__bdd.viewerOf(e);
    const df = sp.dataFrame;
    const xCol = df.col(sp.props.xColumnName);
    const yCol = df.col(sp.props.yColumnName);
    if (!xCol || !yCol)
      return 'it plots no x and y columns';
    if (!df.filter.get(r - 1))
      return `row ${r} does not pass the filter`;
    const p = sp.worldToScreen(xCol.get(r - 1), yCol.get(r - 1));
    const canvas = Array.from(sp.root.querySelectorAll('canvas') as NodeListOf<HTMLCanvasElement>)
      .find((c) => c.getBoundingClientRect().width > 0);
    if (!canvas)
      return 'it shows no canvas';
    const box = canvas.getBoundingClientRect();
    return {x: box.x + p.x, y: box.y + p.y};
  }, row) as {x: number, y: number} | string;
  if (typeof point === 'string')
    throw new Error(`${target.phrase} draws no marker of row ${row}: ${point}`);
  await page.mouse.move(point.x - 3, point.y - 3);
  await page.mouse.move(point.x, point.y);
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx + 1),
    {message: `the row under the pointer on ${target.phrase}`}).toBe(row);
  await expect(page.locator('.d4-tooltip').filter({visible: true}), `the tooltip of the marker of row ${row}`).toHaveCount(1);
  await page.mouse.down();
  await page.mouse.up();
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx + 1),
    {message: `the current row after a click on the marker of row ${row}`}).toBe(row);
}, {tier: 'ui', description: 'the marker of that table row (from 1), where the plot draws it, hovered until its tooltip shows, then clicked; done when the row is current'});
