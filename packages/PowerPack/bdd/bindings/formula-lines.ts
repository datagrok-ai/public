/* The steps the Formula Lines regressions need beyond the library's `viewers` tier: the full SPGI
   demo table (its date columns), and a reading of one viewer compared with the same reading of
   another — the preview inside PowerPack's Formula Lines dialog is a viewer of its own, and the
   claim is that it shows the axes of the viewer it was opened from. */
import {Page} from '@playwright/test';
import {dataset, Then} from '@datagrok-libraries/bdd';
import {ElementRef, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

dataset('spgi-full', {path: 'System:DemoFiles/chem/SPGI.csv',
  description: 'the full SPGI demo table: 3624 compounds; Competition assay Date is a datetime column, Average Mass, TPSA and Chemical Space X numeric ones'});

async function expectAcross(page: Page, name: string, a: ElementRef, b: ElementRef, same: boolean): Promise<void> {
  let shown = '';
  const holds = async (): Promise<boolean> => {
    try {
      const x = String(await viewers.readValue(page, a, name));
      const y = String(await viewers.readValue(page, b, name));
      shown = `${a.phrase} reads ${x}, ${b.phrase} reads ${y}`;
      return (x === y) === same;
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
    throw new Error(`the "${name}" readings should ${same ? 'be the same' : 'differ'}: ${shown}`);
  }
}

export const readingSameAcross = Then('the {string} reading of {widget} should be the same as on {widget}',
  (page: Page, name: string, a: ElementRef, b: ElementRef) => expectAcross(page, name, a, b, true),
  {description: 'one reading of two viewers, as text — an axis range the preview must share with the viewer it was opened from'});

export const readingDiffersAcross = Then('the {string} reading of {widget} should differ from the one on {widget}',
  (page: Page, name: string, a: ElementRef, b: ElementRef) => expectAcross(page, name, a, b, false));

/** The constant a `${column} = <number>` line of the viewer's look sits at, and the X axis range the
 * viewer reports, in the axis's own units: a date written in other units than the axis (milliseconds
 * against microseconds) lands far outside it. */
export const lineWithinXAxis = Then('the {string} line of {widget} should lie within its x axis',
  async (page: Page, column: string, target: ElementRef) => {
    let shown = '';
    const holds = async (): Promise<boolean> => {
      try {
        const look: string = await viewers.onViewer(page, target, (el) => (window as any).__bdd.viewerOf(el).props.formulaLines ?? '');
        const prefix = `\${${column}} = `;
        const line = (JSON.parse(look || '[]') as {formula?: string}[]).find((l) => l.formula?.startsWith(prefix));
        const min = Number(await viewers.readValue(page, target, 'x axis min'));
        const max = Number(await viewers.readValue(page, target, 'x axis max'));
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
