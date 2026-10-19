/* The annotations a macromolecule column carries: the column's `.annotations` tag holds a JSON
   list (the liability scanner's per-row hits live in a `~<column>_annotations` companion column,
   tested with the scanner in src/tests/liability-scanner-tests.ts). */
import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

function annotationNames(page: Page, column: string): Promise<string[]> {
  return page.evaluate((c) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const tag = col.getTag('.annotations');
    return tag ? JSON.parse(tag).map((a: any) => String(a.name)) : [];
  }, column);
}

const remembered = new WeakMap<Page, number>();

export const carriesAtLeast = Then('{string} column should carry at least {int} annotation(s)', async (page: Page, column: string, count: number) => {
  let seen = 0;
  await expect.poll(async () => (seen = (await annotationNames(page, column)).length), {message: `annotations on "${column}"`})
    .toBeGreaterThanOrEqual(count);
  remembered.set(page, seen);
}, {description: 'the column\'s .annotations tag, polled: a scan writes it after the column it adds; remembered for "one annotation fewer than before"'});

export const carriesNone = Then('{string} column should carry no annotations', async (page: Page, column: string) => {
  expect((await annotationNames(page, column)).length, `annotations on "${column}"`).toBe(0);
});

export const oneFewer = Then('{string} column should carry one annotation fewer than before', async (page: Page, column: string) => {
  const before = remembered.get(page);
  if (before === undefined)
    throw new Error('no annotation count was read before ("should carry at least N annotations")');
  const now = (await annotationNames(page, column)).length;
  remembered.set(page, now);
  expect(now, `annotations on "${column}" (before: ${before})`).toBe(before - 1);
});

export const allNamedLike = Then('every annotation on {string} column should be named like {string}', async (page: Page, column: string, pattern: string) => {
  let names: string[] = [];
  await expect.poll(async () => (names = await annotationNames(page, column)).length, {message: `annotations on "${column}"`}).toBeGreaterThan(0);
  expect(names.filter((n) => !new RegExp(pattern).test(n)), `annotations on "${column}" not named like /${pattern}/`).toEqual([]);
}, {description: 'the names in the column\'s .annotations tag, polled for the first: a scan with some rules unchecked lists only the rules left'});

export const listMatchesAnnotations = Then('annotations list in {string} dialog should have as many items as {string} column has annotations', async (page: Page, dialog: string, column: string) => {
  const count = (await annotationNames(page, column)).length;
  await expect(page.locator(`[name="dialog-${dialog.replace(/\s+/g, '-')}"] [data-u2-name="annotations"] [role="listitem"]`), `rows of the annotations list against ${count} annotations`).toHaveCount(count);
});
