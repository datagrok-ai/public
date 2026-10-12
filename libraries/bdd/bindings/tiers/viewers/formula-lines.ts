/* The `viewers` tier, formula lines: the Formula Lines dialog a viewer's Tools > Formula Lines... opens
   (PowerPack's dialogs/formula-lines.ts) — its list of lines, a Line pane that edits the selected one's
   formula, a Format pane with its range — and what a viewer draws of the lines its look or its table
   holds: a "formula line <title>" / "formula band <title>" hit area for each item it drew on screen. */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from '../../../src/runtime/patience.js';
import {Then, When} from '../../../src/registry.js';
import type {ElementRef} from '../../../src/runtime/args.js';
import {hasFocus} from '../../../src/runtime/gestures.js';
import * as v from '../../../src/runtime/viewers.js';

const formulaDialog = (page: Page): Locator => page.locator('[name="dialog-Formula-Lines"]').filter({visible: true}).last();

/** The formulas the dialog's list holds on its current tab, read from the list's own grid. */
async function listedFormulas(page: Page): Promise<string[]> {
  return formulaDialog(page).locator('[name="viewer-Grid"]').filter({visible: true}).first().evaluate((el) => {
    const values = (window as any).__bdd.viewerOf(el).getWidgetStatus()?.values ?? {};
    return Object.keys(values).filter((k) => /^text of cell \d+ of formula$/.test(k))
      .sort((a, b) => Number(/\d+/.exec(a)![0]) - Number(/\d+/.exec(b)![0])).map((k) => String(values[k]));
  });
}

/** "Add new" > Line puts a line with a default formula in the list and selects it; its formula is then
 * typed into the editor's Line pane, which the list takes on blur. */
export const addFormulaLine = When('user adds the formula line {string} in the Formula Lines dialog', async (page: Page, formula: string) => {
  await v.settleAll(page);
  const dialog = formulaDialog(page);
  // the list may hold the formula already (a line over another range): the new one is counted
  const times = async (): Promise<number> => (await listedFormulas(page)).filter((f) => f === formula).length;
  const before = await times();
  await dialog.locator('[name="button-Add-new"]').click();
  await v.pickMenuPath(page, 'Line');
  const editor = dialog.locator('[name="pane-Line"] textarea').first();
  await expect(editor, 'the formula editor of the new line').toBeVisible({timeout: pollMs(5000)});
  await editor.fill(formula);
  await editor.press('Tab');
  await expect.poll(times, {message: `the lines of the Formula Lines dialog with the formula (${before} before)`}).toBe(before + 1);
}, {tier: 'ui', description: 'Add new > Line, the formula typed into the Line pane of the new line; done when the list shows one more line with it'});

export const setFormulaLineRange = When('user sets the range of the selected formula line to {string} .. {string}', async (page: Page, min: string, max: string) => {
  const host = formulaDialog(page).locator('[name="input-host-Range"]').first();
  const inputs = [host.locator('input').first(), host.locator('xpath=following-sibling::div[1]//input').first()];
  for (const [input, value] of [[inputs[0], min], [inputs[1], max]] as [Locator, string][]) {
    await input.fill(value);
    await input.press('Tab');
    await expect(input).toHaveValue(value, {timeout: pollMs(5000)});
  }
}, {tier: 'ui', description: 'the min and the max field of the Range row of the Format pane'});

/** A line of the list is found by its formula: the list's first row is clicked and the Down arrow walks
 * the list, each row showing its item in the editor, until the Line pane holds that formula; the new
 * formula is then typed there. The list shows only a few rows, so a line further down is reached the
 * way a keyboard user reaches it. */
export const editFormulaLine = When('user changes the formula line {string} in the Formula Lines dialog to {string}', async (page: Page, from: string, to: string) => {
  const dialog = formulaDialog(page);
  const grid = dialog.locator('[name="viewer-Grid"]').filter({visible: true}).first();
  const status = () => grid.evaluate((el) => {
    const s = (window as any).__bdd.viewerOf(el).getWidgetStatus();
    const c = (s?.parts?.canvas ?? el.querySelector('canvas')).getBoundingClientRect();
    const a = s?.hitAreas?.['cell 1 of formula'];
    return {rows: Number(s?.values?.['rows'] ?? 0), current: Number(s?.values?.['current row'] ?? 0),
      first: a ? {x: c.x + a.x + a.width / 2, y: c.y + a.y + a.height / 2} : null};
  });
  // the list can still be in its first layout as the dialog opens (its overlay 300×150): a click aimed
  // before it lands on the bare canvas, which takes the focus the Down arrow needs
  await grid.evaluate((el) => (window as any).__bdd.settle(el, 300));
  const st = await status();
  if (!st.first)
    throw new Error('the Formula Lines dialog lists no line');
  // the dialog makes the line of the viewer's axes current on a timer it starts as it opens (DG.delay(1),
  // formula-lines.ts): a timer started now fires after it, so a line picked from here on stays picked
  await page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 1)));
  await page.mouse.click(st.first.x, st.first.y);
  await expect.poll(async () => (await status()).current, {message: 'the current line of the list after the click'}).toBe(1);
  await expect.poll(() => hasFocus(grid), {message: 'the list holding the focus the Down arrow goes to'}).toBe(true);
  const editor = dialog.locator('[name="pane-Line"] textarea').first();
  const seen: string[] = [];
  for (let i = 0; i < st.rows; i++) {
    // a constant line ("${Y} = 100") is edited in a Constant line pane, not the Line pane
    const text = await editor.isVisible() ? await editor.inputValue() : '(not a formula line)';
    seen.push(text);
    if (text === from || i === st.rows - 1)
      break;
    await page.keyboard.press('ArrowDown');
    await expect.poll(async () => (await status()).current, {message: 'the current line of the list after the Down arrow'}).toBe(i + 2);
  }
  if (seen[seen.length - 1] !== from)
    throw new Error(`the Formula Lines dialog has no line "${from}"; it walked: ${seen.join(' | ')}`);
  await editor.fill(to);
  await editor.press('Tab');
  await expect.poll(() => editor.inputValue(), {message: 'the formula of the edited line'}).toBe(to);
}, {tier: 'ui', description: 'the first line of the list clicked, the Down arrow pressed until the editor shows that formula, the new one typed into the Line pane'});

export const drawsFormulaLines = Then('{widget} should draw {int} formula line(s)', async (page: Page, target: ElementRef, count: number) => {
  let names: string[] = [];
  await expect.poll(async () => (names = Object.keys(await v.hitAreas(page, target)).filter((k) => /^formula (line|band) /.test(k))).length,
    {message: `formula line areas of ${target.phrase}`}).toBe(count).catch(() => {
    throw new Error(`${target.phrase} draws ${names.length} formula line(s), not ${count}: ${names.join(', ') || 'none'}`);
  });
}, {description: 'the "formula line <title>" / "formula band <title>" hit areas, which a viewer reports only for an item it drew on screen'});

export const formulaLineRanges = Then('the formula lines of {widget} should hold {string} over the ranges {string}', async (page: Page, target: ElementRef, formula: string, ranges: string) => {
  const want = ranges.split(/\s*;\s*/).filter(Boolean).sort();
  await expect.poll(async () => {
    const raw = await v.readProperty(page, target, 'formulaLines');
    return (JSON.parse(raw || '[]') as any[]).filter((l) => l.formula === formula).map((l) => `${l.min ?? ''}..${l.max ?? ''}`).sort();
  }, {message: `the items of "formulaLines" of ${target.phrase} with that formula, as min..max`}).toEqual(want);
}, {description: 'the items of the viewer\'s formulaLines look with exactly this formula, their min..max pairs (";"-separated, any order)'});

/** A formula line keeps the title it was given when it was made; what it computes is its formula. */
const formulasOf = (raw: string): string[] => (JSON.parse(raw || '[]') as any[]).map((l) => String(l.formula ?? ''));

/** The formulas that name the text, or why there is nothing to read: some lines must be there. */
function namingFormulas(raw: string, text: string): string[] {
  const formulas = formulasOf(raw);
  return formulas.length === 0 ? ['(no formula lines)'] : formulas.filter((f) => f.includes(text));
}

export const viewerFormulasLack = Then('no formula line of {widget} should name {string}', async (page: Page, target: ElementRef, text: string) => {
  await expect.poll(async () => namingFormulas(await v.readProperty(page, target, 'formulaLines'), text),
    {message: `formulas of ${target.phrase} that name "${text}"`}).toEqual([]);
}, {description: 'the formulas of the viewer\'s formula lines (not their titles), of which there must be some'});

export const tableFormulasLack = Then('no formula line of the table should name {string}', async (page: Page, text: string) => {
  await expect.poll(async () => namingFormulas(await page.evaluate(() => String((window as any).grok.shell.tv?.dataFrame?.getTag('.formula-lines') ?? '')), text),
    {message: `formulas of the current table that name "${text}"`}).toEqual([]);
}, {description: 'the formulas of the table\'s ".formula-lines" (not their titles), of which there must be some'});

/** The constant a `${column} = <number>` line of the viewer's look sits at, and the X axis range the
 * viewer reports, in the axis's own units: a date written in other units than the axis (milliseconds
 * against microseconds) lands far outside it. */
export const lineWithinXAxis = Then('the {string} line of {widget} should lie within its x axis', async (page: Page, column: string, target: ElementRef) => {
  let shown = '';
  const holds = async (): Promise<boolean> => {
    try {
      const look: string = await v.onViewer(page, target, (el) => (window as any).__bdd.viewerOf(el).props.formulaLines ?? '');
      const prefix = `\${${column}} = `;
      const line = (JSON.parse(look || '[]') as {formula?: string}[]).find((l) => l.formula?.startsWith(prefix));
      const min = Number(await v.readValue(page, target, 'x axis min'));
      const max = Number(await v.readValue(page, target, 'x axis max'));
      if (!line) {
        shown = `no "${prefix}…" line in its formulaLines: ${look}`;
        return false;
      }
      const value = Number(line.formula!.substring(prefix.length));
      shown = `the line is at ${value}, the x axis runs ${min} to ${max}`;
      return Number.isFinite(value) && value >= min && value <= max;
    }
    catch (e) {
      shown = (e as Error).message;
      return false;
    }
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${column}" line of ${target.phrase} is not within its x axis: ${shown}`);
  }
}, {description: 'the value of a constant line on the X column against the "x axis min" and "x axis max" readings'});
