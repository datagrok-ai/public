/* What a Forge feature leaves on the server, and the claims about it. Every read goes through `grok s` (grok-s.ts) and
   names the feature's own data: a model by its exact name, the runs by the table they were trained on; a global count is
   never a claim, since the user trains on the same stand and four workers train at once. A model is deleted the Forge
   way, in the page (`Forge:deleteModel`: the row, its model file folder, its uploaded training copy, its applications);
   `grok s functions run` would run on the server, where a JS package function does not exist. */
import {readdirSync, statSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, type ElementRef, expect, locate, viewers} from '@datagrok-libraries/bdd/runtime';
import {MINUTE_FORMAT, PREDICTION_TAG} from '../../src/constants.js';
import {applicationCount, countOf, currentCliUser, expectCreatedModel, expectModelCount, expectTables, folderFiles,
  grokS, guardIdentity, knownModel, listOf, MODEL, modelNumber, modelRow, newestApplication, newestRun, quoted,
  requireModelFolder, requireModelRow, rowsOf, RUN, SERVER_POLL, sharesOf, sharingAccount, sweepModels, sweepRuns,
  tablesNamed, valueAt} from './grok-s.js';
import {PANE_SEED, readingText, RESULTS_SEED, settled, shownRun, shownSeed} from './results-state.js';

declare const grok: any;

export const noModelsNamed = Given('no Forge models named {string} are on the server', async (page: Page, list: string) => {
  const names = listOf(list);
  await sweepModels(page, names);
  atFeatureEnd(page, () => sweepModels(page, names));
}, {tier: 'api', description: 'grok s domains query forge.model by the names (comma-separated) and their families older ' +
  'than an hour; each deleted in the page with Forge:deleteModel, now and at feature end; the rows, the model file ' +
  'folders and the training copies read back gone'});

export const noRunsOf = Given('no Forge training runs of the table(s) {string} are on the server', async (page: Page,
  list: string) => {
  const tables = listOf(list);
  await sweepRuns(page, tables);
  atFeatureEnd(page, () => sweepRuns(page, tables));
}, {tier: 'api', description: 'grok s domains delete forge.training_run --filter dataset_name = each table (comma-' +
  'separated), and the tables of their families older than an hour, now and at feature end; the count read back 0'});

export const deleteThroughApi = When('the Forge model {string} is deleted through the API', async (page: Page, name: string) => {
  await guardIdentity(page);
  await page.evaluate((modelName) => grok.functions.call('Forge:deleteModel', {model: modelName}), name);
  await expectModelCount(name, 0, `Forge models named "${name}"`);
}, {tier: 'api', description: 'Forge:deleteModel by the name in the page; the row read back gone'});

export const modelCount = Then('{int} Forge model(s) named {string} should be on the server', async (page: Page, n: number, name: string) => {
  await guardIdentity(page);
  await expectModelCount(name, n, `Forge models named "${name}"`);
}, {tier: 'api', description: 'grok s domains count forge.model --filter name = …'});

const isEmpty = (v: unknown) => v === null || v === undefined || v === '';

/** `| field | value |` rows against [row]: `empty`, `present`, `number` (a finite number, not below 0), or the value
 * as text. Returns the mismatches. */
function mismatches(row: Record<string, unknown>, table: string[][]): string[] {
  const wrong: string[] = [];
  for (const [field, expected] of table) {
    const actual = valueAt(row, field);
    const holds = expected === 'empty' ? isEmpty(actual) : expected === 'present' ? !isEmpty(actual) :
      expected === 'number' ? typeof actual === 'number' && Number.isFinite(actual) && actual >= 0 :
        !isEmpty(actual) && (typeof actual === 'object' ? JSON.stringify(actual) : String(actual)) === expected;
    if (!holds)
      wrong.push(`${field} is ${JSON.stringify(actual)}, not ${expected}`);
  }
  return wrong;
}

/** Polls the row [read] gives (or why there is none) until the | field | value | rows of [table] hold of it. */
async function expectFields(page: Page, read: () => Promise<Record<string, unknown> | string>, table: string[][],
  message: string): Promise<void> {
  await guardIdentity(page);
  await expect.poll(async () => {
    const row = await read();
    return typeof row === 'string' ? [row] : mismatches(row, table);
  }, {...SERVER_POLL, message}).toEqual([]);
}

export const modelShouldHave = Then('the Forge model {string} on the server should have:', (page: Page, name: string,
  table: string[][]) => expectFields(page, () => modelRow(name), table, `the Forge model "${name}" on the server`),
{tier: 'api', description: '| field | value | rows (a dotted path into a json column; empty, present, number or the ' +
  'value as text) against the one row grok s domains query forge.model reads by the name'});

export const readingIsStoredMetric = Then('the {string} reading of {widget} should be the {string} of the Forge model ' +
  '{string} on the server to three decimals', async (page: Page, name: string, target: ElementRef, path: string,
  model: string) => {
  await guardIdentity(page);
  const value = await modelNumber(model, path);
  await expect.poll(() => readingText(page, target, name),
    {message: `"${name}" of ${target.phrase} against ${path} of "${model}" (${value})`}).toBe(value.toFixed(3));
}, {tier: 'api', description: 'a grid cell\'s text against the dotted path of the model row (grok s), formatted with ' +
  'three decimals'});

export const modelFolderHolds = Then('the model file folder of {string} should hold exactly {string}', async (page: Page, name: string, files: string) => {
  await guardIdentity(page);
  const row = await requireModelRow(name);
  const folder = requireModelFolder(name, String(row.blob ?? ''));
  await expect.poll(async () => (await folderFiles(folder)).sort().join(', '),
    {...SERVER_POLL, message: `the files of ${folder}`}).toBe(listOf(files).sort().join(', '));
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
    .toEqual({id: model.dataset_table_id, rowCount: rows, colCount: listOf(columns).length, columns: listOf(columns).join(', ')});
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
  expect(fingerprintProblems(row.dataset_fingerprint, listOf(columns)), `the dataset_fingerprint of "${name}"`).toEqual([]);
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

export const newestRunShouldHave = Then('the newest training run of the table {string} should have:', (page: Page,
  table: string, fields: string[][]) =>
  expectFields(page, () => newestRun(table), fields, `the newest training run of the table "${table}"`),
{tier: 'api', description: '| field | value | rows against the newest forge.training_run row of the table (grok s ' +
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

/* A row a user makes by hand in the Domains view (its ribbon's New model...): the required columns, no model file. */
export const handMadeRow = Given('a Forge model row {string} without a model file is on the server', async (page: Page, name: string) => {
  await guardIdentity(page);
  atFeatureEnd(page, () => sweepModels(page, [name]));
  await grokS(['domains', 'insert', MODEL, `name=${name}`, 'engine_name=XGBoost', 'engine_namespace=Eda',
    'engine_kind=function', 'task=classification', 'target_name=Species', 'storage_mode=none']);
  await expectCreatedModel(name, `the hand-made Forge model row "${name}"`);
}, {tier: 'api', description: 'grok s domains insert forge.model with the name, XGBoost, Eda, function, classification, ' +
  'Species and storage none (no blob, no features, no metrics); swept at feature end'});

export const applicationsCount = Then('the model {string} should have {int} application(s) on the server',
  async (page: Page, name: string, n: number) => {
    await guardIdentity(page);
    const {id} = await knownModel(name);
    await expect.poll(() => applicationCount(id), {...SERVER_POLL, message: `applications of the model "${name}" (${id})`})
      .toBe(n);
  }, {tier: 'api', description: 'grok s domains count forge.application --filter model_id = the model\'s id (the id ' +
  'last read by its name once the model is deleted)'});

export const newestApplicationShouldHave = Then('the newest application of the model {string} should have:',
  (page: Page, name: string, fields: string[][]) => expectFields(page,
    async () => newestApplication((await knownModel(name)).id), fields, `the newest application of the model "${name}"`),
  {tier: 'api', description: '| field | value | rows (empty, present, number or the value as text) against the ' +
  'newest forge.application row of the model (grok s domains query --sort !created_on --limit 1)'});

export const startsWithNewestApplication = Then('{element} should start with the time of the newest application of the ' +
  'model {string}', async (page: Page, target: ElementRef, name: string) => {
  await guardIdentity(page);
  const {id} = await knownModel(name);
  const row = await newestApplication(id);
  if (typeof row === 'string')
    throw new Error(row);
  const when = await page.evaluate(([iso, format]) => (window as any).dayjs(iso).format(format) as string,
    [String(row.created_on), MINUTE_FORMAT]);
  await expect.poll(async () => ((await (await locate(page, target)).first().textContent()) ?? '').trim(),
    {message: `${target.phrase} against the newest application of "${name}" (${String(row.created_on)})`})
    .toMatch(new RegExp(`^${when}`));
}, {tier: 'api', description: 'the element\'s text starts with the created_on of the model\'s newest forge.application ' +
  'row (grok s), in the product\'s minute format (MINUTE_FORMAT, YYYY-MM-DD HH:mm) and the page\'s time zone'});

/* LIB-PROMOTE: the library keeps the files a page downloads in the temp directory as `bdd-<pid>-<time>-<name>` and
   claims their text; a binary file is claimed here by its size against the server's file. */
export const downloadedModelFile = Then('the downloaded file {string} should have the size of the model file of the Forge ' +
  'model {string}', async (page: Page, file: string, name: string) => {
  await guardIdentity(page);
  const {blob} = await knownModel(name);
  const folder = requireModelFolder(name, blob);
  const fileName = blob.split('/').pop();
  const stored = rowsOf(await grokS(['files', 'list', folder])).find((f) => String(f.name) === fileName);
  const size = Number(stored?.size);
  const prefix = `bdd-${process.pid}-`;
  const local = readdirSync(tmpdir()).filter((f) => f.startsWith(prefix) && f.endsWith(`-${file}`))
    .map((f) => join(tmpdir(), f)).sort((a, b) => statSync(a).mtimeMs - statSync(b).mtimeMs).pop();
  if (local === undefined)
    throw new Error(`no downloaded "${file}" of this worker in ${tmpdir()}`);
  expect({downloaded: statSync(local).size, positive: size > 0}, `${local} against ${folder}/${fileName}`)
    .toEqual({downloaded: size, positive: true});
}, {tier: 'api', description: 'the newest file of that name this worker downloaded has as many bytes as the model ' +
  'row\'s blob on the server (grok s files list), which is not empty'});

export const deletedModelFolderGone = Then('the model file folder of the deleted Forge model {string} should be gone',
  async (page: Page, name: string) => {
    await guardIdentity(page);
    const {blob} = await knownModel(name);
    const folder = requireModelFolder(name, blob);
    await expectModelCount(name, 0, `Forge models named "${name}"`);
    await expect.poll(() => folderFiles(folder), {...SERVER_POLL, message: `the files of ${folder}`}).toEqual([]);
  }, {tier: 'api', description: 'no row of the name, and the folder of the blob it had lists no file (grok s files list)'});

const AUTHOR_GRANTS = 'Delete, Edit, Share, View';

export const modelIsEntity = Then('the Forge model {string} should be an entity on the server with the author\'s grants',
  async (page: Page, name: string) => {
    await guardIdentity(page);
    const [{id}, me] = await Promise.all([knownModel(name), currentCliUser()]);
    const [entities, shares] = await Promise.all([grokS(['raw', 'GET', `entities/${id}`]), sharesOf(id)]);
    const [entity] = rowsOf(entities);
    expect({name: entity?.name, friendlyName: entity?.friendlyName}, `the entity ${id}`)
      .toEqual({name: `${MODEL}:${id}`, friendlyName: name});
    const own = shares.filter((s) => s.groupId === me.groupId && s.personal === true);
    expect(own.map((s) => String(s.access)).sort().join(', '), `the grants of ${me.login}'s personal group on ${id}`)
      .toBe(AUTHOR_GRANTS);
  }, {tier: 'api', description: 'grok s raw GET entities/<id> is forge.model:<id> named after the model, and grok s ' +
  'shares list <id> gives the CLI user\'s personal group View, Edit, Delete and Share'});

async function expectSharedWithSharingUser(page: Page, name: string, shared: boolean): Promise<void> {
  await guardIdentity(page);
  const [{id}, user] = await Promise.all([knownModel(name), sharingAccount()]);
  await expect.poll(async () => (await sharesOf(id)).some((s) => s.groupId === user.groupId),
    {...SERVER_POLL, message: `a grant of ${user.name} on the Forge model "${name}"`}).toBe(shared);
}

export const sharedWithSharingUser = Then('the Forge model {string} should be shared with the sharing user on the server',
  (page: Page, name: string) => expectSharedWithSharingUser(page, name, true),
  {tier: 'api', description: 'grok s shares list <id> has a row of the personal group of DATAGROK_SHARING_LOGIN'});

export const notSharedWithSharingUser = Then('the Forge model {string} should not be shared with the sharing user on the server',
  (page: Page, name: string) => expectSharedWithSharingUser(page, name, false),
  {tier: 'api', description: 'grok s shares list <id> has no row of the personal group of DATAGROK_SHARING_LOGIN'});

export const columnTaggedWithModel = Then('the {string} column should be tagged with the id of the Forge model {string}',
  async (page: Page, column: string, name: string) => {
    await guardIdentity(page);
    const {id} = await knownModel(name);
    const tag = await viewers.evaluate(page, ([c, t]) => String((window as any).__bdd.col(c).getTag(t) ?? ''),
      [column, PREDICTION_TAG]);
    expect(tag, `the ${PREDICTION_TAG} tag of "${column}"`).toBe(id);
  }, {tier: 'api', description: 'the column\'s prediction tag (PREDICTION_TAG) is the id of the model row read by its name'});

async function expectSeedOfModel(page: Page, name: string, where: string): Promise<void> {
  await guardIdentity(page);
  const [shown, model] = await Promise.all([shownSeed(page, where), requireModelRow(name)]);
  expect(shown, `the seed ${where} shows, against the seed of the Forge model "${name}"`).toBe(String(model.seed));
}

export const shownSeedIsModels = Then('the shown seed should be the seed of the Forge model {string} on the server',
  (page: Page, name: string) => expectSeedOfModel(page, name, RESULTS_SEED),
  {tier: 'api', description: 'the seed in the Performance panel of the results against the seed column of the model row'});

export const paneSeedIsModels = Then('the seed in the context panel should be the seed of the Forge model {string} on ' +
  'the server', (page: Page, name: string) => expectSeedOfModel(page, name, PANE_SEED),
{tier: 'api', description: 'the seed in the model\'s Performance pane against the seed column of the model row'});
