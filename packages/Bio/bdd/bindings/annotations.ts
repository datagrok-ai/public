/* The annotations a macromolecule column carries: the column's `.annotations` tag holds a JSON
   list, the liability scanner's per-row hits live in a `~<column>_annotations` companion column
   (a JSON list of hits with the monomers matched and where). The steps read both the way the
   package writes them, so a claim about what a scan found is a claim about the data. */
import {expect, Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

declare const grok: any;

interface AnnotationFacts {
  count: number;
  hits: number;
  misplaced: string[];
}

function annotationFacts(page: Page, column: string): Promise<AnnotationFacts> {
  return page.evaluate((c) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const tag = col.getTag('.annotations');
    const count = tag ? JSON.parse(tag).length : 0;
    const rows = df.col(`~${c}_annotations`);
    let hits = 0;
    const misplaced: string[] = [];
    if (rows) {
      for (let i = 0; i < df.rowCount; i++) {
        const raw = rows.get(i);
        if (!raw)
          continue;
        const seq = String(col.get(i) ?? '');
        for (const hit of JSON.parse(raw)) {
          hits++;
          const m = hit.matchedMonomers;
          if (typeof m === 'string' && typeof hit.positionIndex === 'number' && seq.substr(hit.positionIndex, m.length) !== m)
            misplaced.push(`row ${i + 1} at ${hit.positionIndex}: ${m}`);
        }
      }
    }
    return {count, hits, misplaced};
  }, column);
}

const remembered = new WeakMap<Page, number>();

export const carriesAtLeast = Then('{string} column should carry at least {int} annotation(s)', async (page: Page, column: string, count: number) => {
  const f = await annotationFacts(page, column);
  remembered.set(page, f.count);
  expect(f.count, `annotations on "${column}"`).toBeGreaterThanOrEqual(count);
}, {description: 'the column\'s .annotations tag; remembered for "one annotation fewer than before"'});

export const carriesNone = Then('{string} column should carry no annotations', async (page: Page, column: string) => {
  expect((await annotationFacts(page, column)).count, `annotations on "${column}"`).toBe(0);
});

export const oneFewer = Then('{string} column should carry one annotation fewer than before', async (page: Page, column: string) => {
  const before = remembered.get(page);
  if (before === undefined)
    throw new Error('no annotation count was read before ("should carry at least N annotations")');
  const now = (await annotationFacts(page, column)).count;
  remembered.set(page, now);
  expect(now, `annotations on "${column}" (before: ${before})`).toBe(before - 1);
});

const rememberedHits = new WeakMap<Page, number>();

export const hitsMatch = Then('every liability hit on {string} column should match its motif at the position it reports', async (page: Page, column: string) => {
  const f = await annotationFacts(page, column);
  rememberedHits.set(page, f.hits);
  expect(f.hits, `liability hits on "${column}"`).toBeGreaterThan(0);
  expect(f.misplaced, 'hits whose matched monomers are not at the position reported').toEqual([]);
}, {description: 'every hit of the ~<column>_annotations column reads its matched monomers at its position index'});

export const countBelowHits = Then('the total of {string} column should be fewer than the liability hits found before', async (page: Page, counts: string) => {
  const before = rememberedHits.get(page);
  if (before === undefined)
    throw new Error('no liability hits were read before ("every liability hit … should match its motif")');
  const total: number = await page.evaluate((c) => {
    const col = grok.shell.t.col(c);
    let sum = 0;
    for (let i = 0; i < col.length; i++)
      sum += col.isNone(i) ? 0 : Number(col.get(i));
    return sum;
  }, counts);
  expect(total, `the total of "${counts}"`).toBeGreaterThan(0);
  expect(total, `the total of "${counts}" against the ${before} hits the full rule set found`).toBeLessThan(before);
}, {description: 'a narrower rule set finds fewer liabilities than the full one did (the per-row column is rewritten by every scan, so the earlier count is the one remembered)'});

export const listMatchesAnnotations = Then('annotations list in {string} dialog should have as many items as {string} column has annotations', async (page: Page, dialog: string, column: string) => {
  const f = await annotationFacts(page, column);
  await expect(page.locator(`[name="dialog-${dialog.replace(/\s+/g, '-')}"] [data-u2-name="annotations"] [data-u2="item"]`), `rows of the annotations list against ${f.count} annotations`).toHaveCount(f.count);
});
