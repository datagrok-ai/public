/* What a Forge feature leaves on the server, and the claims about it. Every read goes through `grok s` (grok-s.ts) and
   names the feature's own data: a model by its exact name, the runs by the table they were trained on; a global count is
   never a claim, since the user trains on the same stand and four workers train at once. A model is deleted the Forge
   way, in the page (`Forge:deleteModel`: the row, its model file folder, its uploaded training copy, its applications);
   `grok s functions run` would run on the server, where a JS package function does not exist. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect} from '@datagrok-libraries/bdd/runtime';
import {countOf, expectModelCount, expectTables, folderFiles, grokS, guardIdentity, MODEL_FOLDER, modelRow, quoted,
  requireModelRow, rowsOf, RUN, SERVER_POLL, sweepModels, sweepRuns, tablesNamed} from './grok-s.js';
import {settled, shownRun, shownSeed} from './results-state.js';

declare const grok: any;

const namesOf = (list: string): string[] => list.split(',').map((s) => s.trim()).filter((s) => s !== '');

export const noModelsNamed = Given('no Forge models named {string} are on the server', async (page: Page, list: string) => {
  const names = namesOf(list);
  await sweepModels(page, names);
  atFeatureEnd(page, () => sweepModels(page, names));
}, {tier: 'api', description: 'grok s domains query forge.model by the names (comma-separated) and their families older ' +
  'than an hour; each deleted in the page with Forge:deleteModel, now and at feature end; the rows, the model file ' +
  'folders and the training copies read back gone'});

export const noRunsOf = Given('no Forge training runs of the table {string} are on the server', async (page: Page, table: string) => {
  await sweepRuns(page, table);
  atFeatureEnd(page, () => sweepRuns(page, table));
}, {tier: 'api', description: 'grok s domains delete forge.training_run --filter dataset_name = the table, and the tables ' +
  'of its family older than an hour, now and at feature end; the count of each read back 0'});

export const deleteThroughApi = When('the Forge model {string} is deleted through the API', async (page: Page, name: string) => {
  await guardIdentity(page);
  await page.evaluate((modelName) => grok.functions.call('Forge:deleteModel', {model: modelName}), name);
  await expectModelCount(name, 0, `Forge models named "${name}"`);
}, {tier: 'api', description: 'Forge:deleteModel by the name in the page; the row read back gone'});

export const modelCount = Then('{int} Forge model(s) named {string} should be on the server', async (page: Page, n: number, name: string) => {
  await guardIdentity(page);
  await expectModelCount(name, n, `Forge models named "${name}"`);
}, {tier: 'api', description: 'grok s domains count forge.model --filter name = …'});

/** The value at the dotted [path] of [row] ("dataset_fingerprint.rowCount"). */
function valueAt(row: Record<string, unknown>, path: string): unknown {
  let value: unknown = row;
  for (const key of path.split('.'))
    value = typeof value === 'object' && value !== null ? (value as Record<string, unknown>)[key] : undefined;
  return value;
}

const isEmpty = (v: unknown) => v === null || v === undefined || v === '';

/** `| field | value |` rows against [row]: `empty`, `present`, or the value as text. Returns the mismatches. */
function mismatches(row: Record<string, unknown>, table: string[][]): string[] {
  const wrong: string[] = [];
  for (const [field, expected] of table) {
    const actual = valueAt(row, field);
    const holds = expected === 'empty' ? isEmpty(actual) : expected === 'present' ? !isEmpty(actual) :
      !isEmpty(actual) && (typeof actual === 'object' ? JSON.stringify(actual) : String(actual)) === expected;
    if (!holds)
      wrong.push(`${field} is ${JSON.stringify(actual)}, not ${expected}`);
  }
  return wrong;
}

export const modelShouldHave = Then('the Forge model {string} on the server should have:', async (page: Page, name: string, table: string[][]) => {
  await guardIdentity(page);
  await expect.poll(async () => {
    const row = await modelRow(name);
    return typeof row === 'string' ? [row] : mismatches(row, table);
  }, {...SERVER_POLL, message: `the Forge model "${name}" on the server`}).toEqual([]);
}, {tier: 'api', description: '| field | value | rows (a dotted path into a json column; empty, present or the value as ' +
  'text) against the one row grok s domains query forge.model reads by the name'});

export const modelFolderHolds = Then('the model file folder of {string} should hold exactly {string}', async (page: Page, name: string, files: string) => {
  await guardIdentity(page);
  const row = await requireModelRow(name);
  const folder = MODEL_FOLDER.exec(String(row.blob ?? ''))?.[1];
  if (folder === undefined)
    throw new Error(`the Forge model "${name}" has no model file in Forge's folders: ${JSON.stringify(row.blob)}`);
  await expect.poll(async () => (await folderFiles(folder)).sort().join(', '),
    {...SERVER_POLL, message: `the files of ${folder}`}).toBe(namesOf(files).sort().join(', '));
}, {tier: 'api', description: 'the folder of the model row\'s blob, listed by grok s files list'});

/* The copy is the model's own table: its id is the model's dataset_table_id; grok s lists its size, and the page reads
   its column names (the CLI lists no columns). */
export const trainingCopyIs = Then('the training copy of the Forge model {string} should be the table {string} of {int} rows ' +
  'and the columns {string}', async (page: Page, name: string, table: string, rows: number, columns: string) => {
  await guardIdentity(page);
  const [listed, model] = await Promise.all([tablesNamed(table), requireModelRow(name)]);
  if (listed.length !== 1)
    throw new Error(`${listed.length} server tables named "${table}"`);
  const [copy] = listed;
  const names = await page.evaluate(async (id) => (await grok.dapi.tables.getTable(id)).columns.names().join(', '),
    String(copy.id));
  expect({id: copy.id, rowCount: copy.rowCount, colCount: copy.colCount, columns: names},
    `the server table "${table}" against the model's dataset_table_id`)
    .toEqual({id: model.dataset_table_id, rowCount: rows, colCount: namesOf(columns).length, columns: namesOf(columns).join(', ')});
}, {tier: 'api', description: 'grok s tables list by the name (one table): its id is the model\'s dataset_table_id, its ' +
  'rowCount and colCount; the column names from grok.dapi.tables.getTable in the page'});

const FINGERPRINT_KEYS = ['columnCount', 'columns', 'hash', 'rowCount'];
const COLUMN_KEYS = new Set(['name', 'type', 'missingCount', 'min', 'max', 'mean', 'categories']);

/** What is wrong with [fp] as a fingerprint of [columns]: its keys, the column entries, an array of values. */
function fingerprintProblems(fp: unknown, columns: string[]): string[] {
  if (typeof fp !== 'object' || fp === null)
    return [`no dataset_fingerprint: ${JSON.stringify(fp)}`];
  const f = fp as {rowCount?: unknown; hash?: unknown; columns?: unknown};
  const problems: string[] = [];
  if (Object.keys(f).sort().join(',') !== FINGERPRINT_KEYS.join(','))
    problems.push(`keys ${Object.keys(f).sort().join(', ')}`);
  if (typeof f.hash !== 'string' || !/^[0-9a-f]{8}$/.test(f.hash))
    problems.push(`hash ${JSON.stringify(f.hash)}`);
  const entries = Array.isArray(f.columns) ? f.columns as Record<string, unknown>[] : [];
  if (entries.map((c) => c.name).join(', ') !== columns.join(', '))
    problems.push(`columns ${entries.map((c) => String(c.name)).join(', ')}`);
  for (const c of entries) {
    const extra = Object.keys(c).filter((k) => !COLUMN_KEYS.has(k));
    if (extra.length > 0)
      problems.push(`${String(c.name)} has ${extra.join(', ')}`);
    for (const [key, value] of Object.entries(c)) {
      const bad = key === 'categories' ?
        !Array.isArray(value) || new Set(value).size !== value.length || value.length >= Number(f.rowCount) :
        typeof value === 'object' && value !== null;
      if (bad)
        problems.push(`${String(c.name)}.${key} is ${JSON.stringify(value).slice(0, 80)}`);
    }
  }
  return problems;
}

export const fingerprintDescribes = Then('the fingerprint of the Forge model {string} should describe the columns {string} ' +
  'and hold no row values', async (page: Page, name: string, columns: string) => {
  await guardIdentity(page);
  const row = await requireModelRow(name);
  expect(fingerprintProblems(row.dataset_fingerprint, namesOf(columns)), `the dataset_fingerprint of "${name}"`).toEqual([]);
}, {tier: 'api', description: 'exactly rowCount, columnCount, an 8-digit hex hash and one entry per column in order; an ' +
  'entry holds its name, type, missing count, min/max/mean or distinct categories, nothing else, so no row values'});

export const noTableNamed = Then('no table named {string} should be on the server', async (page: Page, name: string) => {
  await guardIdentity(page);
  await expectTables(name, 0);
}, {tier: 'api', description: 'grok s tables list --filter, the names compared exactly'});

export const tableNamed = Then('a table named {string} should be on the server', async (page: Page, name: string) => {
  await guardIdentity(page);
  await expectTables(name, 1);
}, {tier: 'api', description: 'grok s tables list --filter, the names compared exactly; exactly one'});

async function newestRun(table: string): Promise<Record<string, unknown> | string> {
  const rows = rowsOf(await grokS(['domains', 'query', RUN, '--filter', `dataset_name = ${quoted(table)}`,
    '--sort', '!created_on', '--limit', '1']));
  return rows[0] ?? `no training run of the table "${table}"`;
}

export const newestRunShouldHave = Then('the newest training run of the table {string} should have:', async (page: Page, table: string, fields: string[][]) => {
  await guardIdentity(page);
  await expect.poll(async () => {
    const row = await newestRun(table);
    return typeof row === 'string' ? [row] : mismatches(row, fields);
  }, {...SERVER_POLL, message: `the newest training run of the table "${table}"`}).toEqual([]);
}, {tier: 'api', description: '| field | value | rows against the newest forge.training_run row of the table (grok s ' +
  'domains query --sort !created_on --limit 1)'});

export const shownIsNewestRun = Then('the shown training should be the newest training run of the table {string}, with:',
  async (page: Page, table: string, fields: string[][]) => {
    await guardIdentity(page);
    const shown = await shownRun(page);
    if (shown === null)
      throw new Error('the Train view shows no training\'s result (no data-forge-run-id)');
    await expect.poll(async () => {
      const row = await newestRun(table);
      if (typeof row === 'string')
        return [row];
      return row.id === shown ? mismatches(row, fields) : [`the newest run is ${String(row.id)}, the view shows ${shown}`];
    }, {...SERVER_POLL, message: `the shown training and the newest training run of the table "${table}"`}).toEqual([]);
  }, {tier: 'api', description: 'the run id the results show (data-forge-run-id) is the newest forge.training_run row of ' +
  'the table, and | field | value | rows hold of it'});

export const newestRunLinked = Then('the newest training run of the table {string} should be linked to the model {string}',
  async (page: Page, table: string, name: string) => {
    await guardIdentity(page);
    await expect.poll(async () => {
      const [run, model] = await Promise.all([newestRun(table), modelRow(name)]);
      if (typeof run === 'string' || typeof model === 'string')
        return typeof run === 'string' ? run : model;
      return run.model_id === model.id ? 'linked' : `model_id ${JSON.stringify(run.model_id)}, the model is ${model.id}`;
    }, {...SERVER_POLL, message: `the newest training run of "${table}" and the model "${name}"`}).toBe('linked');
  }, {tier: 'api', description: 'the run\'s model_id is the model row\'s id'});

const rememberedRuns = new WeakMap<Page, Map<string, number>>();

export const rememberRunCount = When('user remembers the training run count of the table {string}', async (page: Page, table: string) => {
  await guardIdentity(page);
  await settled(page);
  const map = rememberedRuns.get(page) ?? new Map<string, number>();
  rememberedRuns.set(page, map);
  map.set(table, await countOf(RUN, `dataset_name = ${quoted(table)}`));
}, {tier: 'api', description: 'once the Train view has settled (aria-busy not true), grok s domains count ' +
  'forge.training_run --filter dataset_name = the table'});

async function expectRunCount(page: Page, table: string, more: number): Promise<void> {
  await guardIdentity(page);
  await settled(page);
  const before = rememberedRuns.get(page)?.get(table);
  if (before === undefined)
    throw new Error(`the training run count of "${table}" was not remembered`);
  await expect.poll(() => countOf(RUN, `dataset_name = ${quoted(table)}`),
    {...SERVER_POLL, message: `training runs of "${table}", ${before} remembered`}).toBe(before + more);
}

export const runCountAsRemembered = Then('the training run count of the table {string} should be as remembered',
  (page: Page, table: string) => expectRunCount(page, table, 0),
  {tier: 'api', description: 'read once the Train view has settled (aria-busy not true, valid or not): a run it ' +
  'records is written before its result is shown'});

export const runCountOneMore = Then('the training run count of the table {string} should be one more than remembered',
  (page: Page, table: string) => expectRunCount(page, table, 1),
  {tier: 'api', description: 'read once the Train view has settled, as above'});

export const shownSeedIsModels = Then('the shown seed should be the seed of the Forge model {string} on the server',
  async (page: Page, name: string) => {
    await guardIdentity(page);
    const shown = await shownSeed(page);
    await expect.poll(async () => {
      const row = await modelRow(name);
      return typeof row === 'string' ? row : String(row.seed);
    }, {...SERVER_POLL, message: `the seed of the Forge model "${name}", the Performance panel shows ${shown}`}).toBe(shown);
  }, {tier: 'api', description: 'the seed in the Performance panel against the seed column of the model row'});
