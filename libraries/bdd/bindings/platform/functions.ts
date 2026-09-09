/* A package function called through the platform (`grok.functions.call`) and its result — the
   API surface a package promises other packages, checked the way they use it. The last result
   is kept for the checks that follow; a rejected call fails the step with the platform's
   message. A package's own readings of the result build on `readResult` from the runtime. */
import {expect, Page} from '@playwright/test';
import {Then, When} from '../../src/registry.js';
import {callFunction, readResult} from '../../src/runtime/functions.js';

export const callWith = When('user calls {string} function with:', (page: Page, name: string, rows: string[][]) => callFunction(page, name, rows),
  {tier: 'api', description: '"Package:function"; a | name | value | table of arguments — `column:X` and `table` name the current table\'s column and the table itself'});

export const call = When('user calls {string} function', (page: Page, name: string) => callFunction(page, name, []), {tier: 'api'});

export const resultEmpty = Then('the result should be empty', async (page: Page) => {
  expect(await readResult(page, `value == null || value === '' ? 'empty' : typeof value + ': ' + String(value).slice(0, 80)`), 'the last result').toBe('empty');
}, {description: 'null, undefined or an empty string'});

export const resultContains = Then('the result should contain text {string}', async (page: Page, text: string) => {
  expect(String(await readResult(page, `String(value ?? '').slice(0, 4000)`)), 'the last result').toContain(text);
});

export const resultMatches = Then('the result should match {string}', async (page: Page, pattern: string) => {
  expect(String(await readResult(page, `String(value ?? '').slice(0, 4000)`)), 'the last result').toMatch(new RegExp(pattern));
});

export const resultHasMethods = Then('the result should have methods {string}', async (page: Page, list: string) => {
  const names = list.split(/\s*,\s*/).filter(Boolean);
  const missing = await readResult(page, `arg.filter((m) => typeof value?.[m] !== 'function')`, names);
  expect(missing, 'methods missing on the last result').toEqual([]);
}, {description: 'comma-separated names, each a function of the result — a service object another package would call'});

export const resultIsList = Then('the result should be a list of {int} or more items', async (page: Page, count: number) => {
  expect(await readResult(page, `Array.isArray(value) ? value.length : (value?.length ?? -1)`), 'items in the last result').toBeGreaterThanOrEqual(count);
});

export const resultTableColumns = Then('the result should be a table with columns {string}', async (page: Page, list: string) => {
  expect(await readResult(page, `value?.columns?.names?.() ?? ('not a table: ' + typeof value)`), 'columns of the result table').toEqual(list.split(/\s*,\s*/));
});

export const resultTableFilled = Then('every column of the result table should be filled in row {int}', async (page: Page, row: number) => {
  const blanks = await readResult(page, `value.columns.toList().filter((c) => c.isNone(arg - 1)).map((c) => c.name)`, row);
  expect(blanks, `blank columns in row ${row} of the result table`).toEqual([]);
});

export const resultProperty = Then('the result should have a {string} of {string}', async (page: Page, name: string, value: string) => {
  expect(String(await readResult(page, `String(value?.[arg] ?? '')`, name)), `"${name}" of the last result`).toBe(value);
}, {description: 'a property of the result object by name, read as text'});
