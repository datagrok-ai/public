/* A package function called through the platform (`grok.functions.call`) and its result — the
   API surface a package promises other packages, checked the way they use it. The last result
   is kept for the checks that follow; a rejected call fails the step with the platform's
   message. A package's own readings of the result build on `readResult` from the runtime. */
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
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

export const resultIsNumber = Then('the result should be the number {float}', async (page: Page, expected: number) => {
  const value = await readResult(page, 'value');
  expect(typeof value, `the result ${JSON.stringify(value)}`).toBe('number');
  expect(value as number, 'the last result').toBeCloseTo(expected, 6);
});

export const resultNumberBetween = Then('the result should be a number between {float} and {float}', async (page: Page, min: number, max: number) => {
  const value = await readResult(page, 'value');
  expect(typeof value, `the result ${JSON.stringify(value)}`).toBe('number');
  expect(value as number, 'the last result').toBeGreaterThanOrEqual(min);
  expect(value as number, 'the last result').toBeLessThanOrEqual(max);
});

export const resultColumnLength = Then('the result should be a column of {int} values', async (page: Page, count: number) => {
  expect(await readResult(page, `value?.length != null && typeof value?.get === 'function' ? value.length : 'not a column: ' + typeof value`),
    'values in the result column').toBe(count);
}, {description: 'a column the function returned, not added to any table'});

export const resultColumnValue = Then('row {int} of the result column should be {string}', async (page: Page, row: number, expected: string) => {
  const actual = await readResult(page, `value?.get ? String(value.get(arg - 1) ?? '') : 'not a column: ' + typeof value`, row);
  expect(actual, `row ${row} of the returned column`).toBe(expected);
}, {description: 'rows count from 1'});

export const resultEveryStarts = Then('every value of the result should start with {string}', async (page: Page, prefix: string) => {
  const facts: {n: number; bad: string[]} = await readResult(page,
    `value?.get ? {n: value.length, bad: Array.from({length: value.length}, (_, i) => String(value.get(i) ?? '')).map((s, i) => [i + 1, s]).filter(([, s]) => !s.startsWith(arg)).slice(0, 5).map(([i, s]) => i + ': ' + s.slice(0, 40))} : {n: 0, bad: ['not a column: ' + typeof value]}`, prefix);
  expect(facts.bad, `rows of the result not starting with "${prefix}"`).toEqual([]);
  expect(facts.n, 'values in the result column').toBeGreaterThan(0);
}, {description: 'a returned column; an empty one fails'});
