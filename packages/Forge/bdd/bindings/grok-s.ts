/* The server claims of the Forge features read the server through the `grok s` CLI, run from Node as a child process
   (no shell, so nothing rewrites a path). The CLI is the workspace's own datagrok-tools; it signs in with the keypair of
   the alias FORGE_BDD_GROK_HOST (`local` by default, ~/.grok/config.yaml), and that user must be the page's: the guard
   below fails the first claim otherwise, since the rows a user sees are the rows granted to them. The reads and the
   sweeps shared by the step modules live here: a step module imported by another would be loaded twice. */
import {appendFileSync} from 'node:fs';
import {execFile} from 'node:child_process';
import {createRequire} from 'node:module';
import type {Page} from '@playwright/test';
import {atFeatureEnd, expect, fixtureFamilies, isStaleFixture, pollMs, serverNow} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const GROK_HOST = process.env.FORGE_BDD_GROK_HOST ?? 'local';
/** A poll of the server: every call signs in, so the reads are spaced, not the page's 100 ms. */
export const SERVER_POLL = {intervals: [1000, 2000, 2000, 3000], timeout: pollMs(30000)};
const CALL_TIMEOUT_MS = 60000;

export const MODEL = 'forge.model';
export const RUN = 'forge.training_run';
const APPLICATION = 'forge.application';
const MODEL_FOLDER = /^file:\/\/(System:DomainFiles\/forge\/model\/[0-9a-f-]{36})\/[^/]+$/i;
const COPY_SUFFIX = ' (training data)';
const DELETE_PAGES = 20;

let bin: string | undefined;

function grokBin(): string {
  if (bin === undefined) {
    try {
      bin = createRequire(import.meta.url).resolve('datagrok-tools/bin/grok.js');
    }
    catch {
      throw new Error('the grok CLI (datagrok-tools) does not resolve from the Forge package: run `pnpm install` in public/');
    }
  }
  return bin;
}

/** `grok s <args> --host <alias> --output json`, parsed; a failed call is an error with the CLI's own message.
 * FORGE_BDD_GROK_LOG=<file> appends one line per call (the arguments and the milliseconds it took). */
export function grokS(args: string[]): Promise<unknown> {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    execFile(process.execPath, [grokBin(), 's', ...args, '--host', GROK_HOST, '--output', 'json'],
      {timeout: CALL_TIMEOUT_MS, maxBuffer: 64 * 1024 * 1024, windowsHide: true}, (error, stdout, stderr) => {
        const log = process.env.FORGE_BDD_GROK_LOG;
        if (log)
          appendFileSync(log, JSON.stringify({args, ms: Date.now() - start, ok: error === null}) + '\n');
        if (error !== null) {
          // Node kills the child on a full buffer too, so the code is read before `killed`.
          const code = (error as {code?: unknown}).code;
          reject(new Error(code === 'ERR_CHILD_PROCESS_STDIO_MAXBUFFER' ? `grok s printed more than 64 MB: ${args.join(' ')}` :
            error.killed ? `grok s timed out: ${args.join(' ')}` :
              `grok s ${args.join(' ')}: ${(stderr || stdout || error.message).trim()}`));
          return;
        }
        try {
          resolve(stdout.trim() === '' ? null : JSON.parse(stdout));
        }
        catch {
          reject(new Error(`grok s ${args.join(' ')} printed no JSON: ${stdout.slice(0, 300)}`));
        }
      });
  });
}

/** A JSON value the CLI printed as a list of rows. */
export function rowsOf(value: unknown): Record<string, unknown>[] {
  if (!Array.isArray(value))
    throw new Error(`grok s printed ${JSON.stringify(value)?.slice(0, 200)} where a list was expected`);
  return value.filter((r): r is Record<string, unknown> => typeof r === 'object' && r !== null);
}

/** `domains count` prints `{"count": N}`. */
export async function countOf(table: string, filter: string): Promise<number> {
  const value = await grokS(['domains', 'count', table, '--filter', filter]);
  const count = (value as {count?: unknown} | null)?.count;
  if (typeof count !== 'number')
    throw new Error(`grok s domains count ${table} printed ${JSON.stringify(value)}`);
  return count;
}

/** The items of a comma-separated list in a step, trimmed, the empty ones dropped. */
export const listOf = (list: string): string[] => list.split(',').map((s) => s.trim()).filter((s) => s !== '');

export const quoted = (s: string): string => `"${s.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}"`;

interface CliUser {login: string; groupId: string}
let cliUser: Promise<CliUser> | undefined;

/** The id of a user's personal group, as `grok s` prints the user. */
function groupIdOf(user: {login?: unknown; group?: unknown}): string {
  const id = (user.group as {id?: unknown} | undefined)?.id;
  if (typeof id !== 'string' || id === '')
    throw new Error(`grok s printed the user ${String(user.login)} without a personal group`);
  return id;
}

/** The CLI's user: its login and its personal group, read once. */
export function currentCliUser(): Promise<CliUser> {
  cliUser ??= grokS(['raw', 'GET', 'users/current']).then((u) => {
    const user = (u ?? {}) as {login?: unknown; group?: unknown};
    return {login: String(user.login ?? ''), groupId: groupIdOf(user)};
  });
  return cliUser;
}

/** Fails when the CLI's user is not the page's: the claims would read another user's rows. The page's login is read
 * every time (a feature can sign in as another account); the CLI's once. */
export async function guardIdentity(page: Page): Promise<void> {
  const [{login: cli}, shown] = await Promise.all([currentCliUser(), page.evaluate(() => String(grok.shell.user.login))]);
  if (cli !== shown) {
    throw new Error(`the grok s alias "${GROK_HOST}" is "${cli}", the page is "${shown}": the server claims would read ` +
      'another user\'s rows (set FORGE_BDD_GROK_HOST to an alias of the page\'s user)');
  }
}

/* EMS returns `created_on` with its zone (`2026-10-06T16:12:19.015442Z` in `--output json`), so Date.parse reads the
   server's instant. */
const createdOn = (row: Record<string, unknown>): number => Date.parse(String(row.created_on ?? '')) || 0;

/** The `forge-test-…` rows of [table] (a name column [column]) that are a stale member of the family of [names]: the
 * ones a killed run of the feature left an hour or more ago, by the server's clock. Read oldest first, so the limit
 * cuts the newest. */
export async function fixtureRows(page: Page, table: string, column: string, names: string[], columns: string[]):
  Promise<Record<string, unknown>[]> {
  const [listed, now] = await Promise.all([grokS(['domains', 'query', table, '--filter', `${column} like "forge-test-%"`,
    '--columns', [column, ...columns].join(','), '--sort', 'created_on', '--limit', '1000']), serverNow(page)]);
  const families = fixtureFamilies(names);
  return rowsOf(listed).filter((r) => {
    const name = String(r[column] ?? '');
    return isStaleFixture({name, friendlyName: name, createdOn: createdOn(r)}, families, now);
  });
}

/** The [terms] (`column = "value"`) in filters of 20 terms each. */
function filtersOf(terms: string[]): string[] {
  const filters: string[] = [];
  for (let i = 0; i < terms.length; i += 20)
    filters.push(terms.slice(i, i + 20).join(' or '));
  return filters;
}

/** The count of [table] rows matching any of the [terms]. */
async function countWhere(table: string, terms: string[]): Promise<number> {
  const counts = await Promise.all(filtersOf(terms).map((filter) => countOf(table, filter)));
  return counts.reduce((a, b) => a + b, 0);
}

const termsOf = (column: string, values: string[]): string[] => values.map((v) => `${column} = ${quoted(v)}`);

/** Every model row read by its name in this worker, as last read: the claims about a deleted model (its applications,
 * its folder) need it once the row is gone, and the fields a model keeps (id, blob, metrics, seed, fingerprint, table)
 * are read from it. A sweep forgets its names. */
const knownModels = new Map<string, Record<string, unknown>>();

/** The one model row named [name], read now, or why there is not exactly one; a row found is remembered by its name.
 * The read for a field that can change (tags, description). */
export async function modelRow(name: string): Promise<Record<string, unknown> | string> {
  const rows = rowsOf(await grokS(['domains', 'query', MODEL, '--filter', `name = ${quoted(name)}`, '--limit', '2']));
  if (rows.length !== 1)
    return `${rows.length} Forge models named "${name}"`;
  knownModels.set(name, rows[0]);
  return rows[0];
}

/** The row of the model [name] as last read, read now if it never was: the read for the fields a model keeps. */
export async function requireModelRow(name: string): Promise<Record<string, unknown>> {
  const known = knownModels.get(name);
  if (known !== undefined)
    return known;
  const row = await modelRow(name);
  if (typeof row === 'string')
    throw new Error(row);
  return row;
}

/** The id and the model file of the model [name], also once it is deleted. */
export async function knownModel(name: string): Promise<{id: string; blob: string}> {
  const row = await requireModelRow(name);
  return {id: String(row.id), blob: String(row.blob ?? '')};
}

/** Waits for the one row of the model [name] a step has just created, and remembers it. */
export async function expectCreatedModel(name: string, message: string): Promise<void> {
  await expect.poll(async () => {
    const row = await modelRow(name);
    return typeof row === 'string' ? row : 'one row';
  }, {...SERVER_POLL, message}).toBe('one row');
}

/** The number at the dotted [path] of the model [name]'s row (a stored metric); fails without one. */
export async function modelNumber(name: string, path: string): Promise<number> {
  const value = valueAt(await requireModelRow(name), path);
  if (typeof value !== 'number')
    throw new Error(`the Forge model "${name}" has no number at ${path}: ${JSON.stringify(value)}`);
  return value;
}

/** The newest row of [table] whose [column] is [value], or what was found instead. */
async function newestRow(table: string, column: string, value: string): Promise<Record<string, unknown> | string> {
  const rows = rowsOf(await grokS(['domains', 'query', table, '--filter', `${column} = ${quoted(value)}`,
    '--sort', '!created_on', '--limit', '1']));
  return rows[0] ?? `no ${table} row with ${column} = ${value}`;
}

/** The applications of the model [id] on the server: their count, or the newest. */
export const applicationCount = (id: string): Promise<number> => countOf(APPLICATION, `model_id = ${quoted(id)}`);

export const newestApplication = (id: string): Promise<Record<string, unknown> | string> =>
  newestRow(APPLICATION, 'model_id', id);

/** The newest training run of the table [table]. */
export const newestRun = (table: string): Promise<Record<string, unknown> | string> =>
  newestRow(RUN, 'dataset_name', table);

/** The training run [id], or why there is none. */
export const runOf = (id: string): Promise<Record<string, unknown> | string> => newestRow(RUN, 'id', id);

/** The value at the dotted [path] of [row] ("dataset_fingerprint.rowCount"). */
export function valueAt(row: Record<string, unknown>, path: string): unknown {
  let value: unknown = row;
  for (const key of path.split('.'))
    value = typeof value === 'object' && value !== null ? (value as Record<string, unknown>)[key] : undefined;
  return value;
}

/** The folder of the model file [blob] (a model row's `blob`), or null when it is not in Forge's folders. */
export const modelFolderOf = (blob: string): string | null => MODEL_FOLDER.exec(blob)?.[1] ?? null;

/** The folder of the model file of the model [name], whose blob was [blob]; fails when it is not in Forge's folders. */
export function requireModelFolder(name: string, blob: string): string {
  const folder = modelFolderOf(blob);
  if (folder === null)
    throw new Error(`the Forge model "${name}" has no model file in Forge's folders: ${JSON.stringify(blob)}`);
  return folder;
}

/** `shares list <id>`: one row per group and permission ({group, groupId, access, personal, inherited}). */
export async function sharesOf(id: string): Promise<Record<string, unknown>[]> {
  return rowsOf(await grokS(['shares', 'list', id]));
}

interface SharingAccount {login: string; groupId: string; name: string}
let sharingUser: Promise<SharingAccount> | undefined;

/** The second account of a sharing feature (DATAGROK_SHARING_LOGIN): its personal group and the name the platform
 * shows for it, read once. */
export function sharingAccount(): Promise<SharingAccount> {
  const login = process.env.DATAGROK_SHARING_LOGIN;
  if (!login)
    throw new Error('no second account: set DATAGROK_SHARING_LOGIN');
  sharingUser ??= grokS(['users', 'list', '--filter', login]).then((listed) => {
    const user = rowsOf(listed).find((u) => u.login === login);
    if (user === undefined)
      throw new Error(`no user "${login}" on the server`);
    return {login, groupId: groupIdOf(user), name: String(user.friendlyName ?? login)};
  });
  return sharingUser;
}

/** The grants of the group [groupId] on the model [id], as the page reads them: `grok.dapi.permissions.get` reads the
 * project the platform keeps the row's shares in (made by the first share), as `revoke` writes it; an empty list
 * comes back undefined. Null when the page gets no entity of [id] (a deleted model). */
function pageGrantsOf(page: Page, id: string, groupId: string): Promise<number | null> {
  return page.evaluate(async ([entityId, group]) => {
    const [entity] = await grok.dapi.getEntities([entityId]);
    if (entity == null)
      return null;
    const {view, edit} = await grok.dapi.permissions.get(entity);
    return [...view ?? [], ...edit ?? []].filter((g: {id: string}) => g.id === group).length;
  }, [id, groupId]);
}

/** Takes back every grant of the second account on the model [id] named [name] (a run that failed between the share
 * and its removal): `grok.dapi.permissions.revoke` in the page, read back gone the same way. A deleted model's entity
 * keeps its grants but no longer reaches the page: `grok s shares list` then fails on a grant left. */
export async function revokeSharingGrant(page: Page, id: string, name: string): Promise<void> {
  const {groupId, name: shown} = await sharingAccount();
  const grants = await pageGrantsOf(page, id, groupId);
  if (grants === null && (await sharesOf(id)).some((s) => s.groupId === groupId))
    throw new Error(`the deleted Forge model "${name}" (${id}) keeps a grant of ${shown}, which the page can no longer revoke`);
  if (!grants)
    return;
  await page.evaluate(async ([entityId, group]) => {
    const [entity] = await grok.dapi.getEntities([entityId]);
    await grok.dapi.permissions.revoke(entity, await grok.dapi.groups.find(group));
  }, [id, groupId]);
  await expect.poll(() => pageGrantsOf(page, id, groupId),
    {message: `the grants of ${shown} on the Forge model "${name}" after the revoke`}).toBe(0);
}

export async function expectModelCount(name: string, n: number, message: string): Promise<void> {
  await expect.poll(() => countOf(MODEL, `name = ${quoted(name)}`), {...SERVER_POLL, message}).toBe(n);
}

/** The file names in the server folder [folder]. */
export async function folderFiles(folder: string): Promise<string[]> {
  return rowsOf(await grokS(['files', 'list', folder])).map((f) => String(f.name));
}

/** The server tables whose name or friendly name is exactly [name]. */
export async function tablesNamed(name: string): Promise<Record<string, unknown>[]> {
  return rowsOf(await grokS(['tables', 'list', '--filter', name])).filter((t) => t.friendlyName === name || t.name === name);
}

export async function expectTables(name: string, n: number): Promise<void> {
  await expect.poll(async () => (await tablesNamed(name)).length,
    {...SERVER_POLL, message: `server tables named "${name}"`}).toBe(n);
}

/** Deletes the feature's models (read by their exact names) and the stale ones of their families in the page, then,
 * when it deleted any, reads the server back by the exact names and ids: no such row, no model file folder of theirs,
 * no training copy of theirs. */
export async function sweepModels(page: Page, names: string[]): Promise<void> {
  await guardIdentity(page);
  for (const name of names)
    knownModels.delete(name);
  const columns = ['name', 'blob', 'dataset_table_id'];
  const [stale, ...exact] = await Promise.all([fixtureRows(page, MODEL, 'name', names, columns.slice(1)),
    ...filtersOf(termsOf('name', names)).map((filter) => grokS(['domains', 'query', MODEL, '--filter', filter,
      '--columns', columns.join(','), '--limit', '1000']))]);
  const rows = [...new Map([...exact.flatMap(rowsOf), ...stale].map((r) => [String(r.id), r])).values()];
  if (rows.length === 0)
    return;
  const ids = rows.map((r) => String(r.id));
  await page.evaluate((list) => Promise.all(list.map((id) => grok.functions.call('Forge:deleteModel', {model: id}))), ids);
  await expect.poll(() => countWhere(MODEL, [...termsOf('name', names), ...termsOf('id', ids)]),
    {...SERVER_POLL, message: `Forge models named ${names.join(', ')} or deleted by id (${ids.length}) left on the server`})
    .toBe(0);
  await Promise.all(rows.map(async (row) => {
    const folder = modelFolderOf(String(row.blob ?? ''));
    if (folder !== null) {
      await expect.poll(() => folderFiles(folder),
        {...SERVER_POLL, message: `the model file folder ${folder} of "${row.name}"`}).toEqual([]);
    }
    if (row.dataset_table_id)
      await expectTables(`${row.name}${COPY_SUFFIX}`, 0);
  }));
}

/** Deletes the runs of [names] (tables, and the stale tables of their families) and reads the count of all back. */
export async function sweepRuns(page: Page, names: string[]): Promise<void> {
  await guardIdentity(page);
  const stale = (await fixtureRows(page, RUN, 'dataset_name', names, [])).map((r) => String(r.dataset_name));
  const tables = [...new Set([...names, ...stale])];
  await Promise.all(tables.map(async (name) => {
    let pages = 0;
    for (let more = true; more; pages++) {
      if (pages === DELETE_PAGES)
        throw new Error(`the runs of "${name}" still have more after ${DELETE_PAGES} deletes of 1000`);
      const result = await grokS(['domains', 'delete', RUN, '--filter', `dataset_name = ${quoted(name)}`, '--limit', '1000']);
      more = (result as {hasMore?: unknown} | null)?.hasMore === true;
    }
  }));
  await expect.poll(() => countWhere(RUN, termsOf('dataset_name', tables)),
    {...SERVER_POLL, message: `training runs of the tables ${tables.join(', ')} left on the server`}).toBe(0);
}

const endSweeps = new WeakMap<Page, {models: Set<string>; tables: Set<string>}>();

/** Adds the models [models] and the runs of the tables [tables] to the feature-end sweep of the page: one feature-end
 * cleanup sweeps all the feature named, the models and the runs at once. */
export function sweepAtFeatureEnd(page: Page, models: string[], tables: string[]): void {
  let sweep = endSweeps.get(page);
  if (sweep === undefined) {
    const s = {models: new Set<string>(), tables: new Set<string>()};
    endSweeps.set(page, sweep = s);
    atFeatureEnd(page, async () => {
      endSweeps.delete(page);
      await Promise.all([s.models.size > 0 ? sweepModels(page, [...s.models]) : null,
        s.tables.size > 0 ? sweepRuns(page, [...s.tables]) : null]);
    });
  }
  for (const name of models)
    sweep.models.add(name);
  for (const table of tables)
    sweep.tables.add(table);
}
