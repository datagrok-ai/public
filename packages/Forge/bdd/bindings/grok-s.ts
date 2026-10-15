/* The server claims of the Forge features read the server through the `grok s` CLI, run from Node as a child process
   (no shell, so nothing rewrites a path). The CLI is the workspace's own datagrok-tools; it signs in with the keypair of
   the alias FORGE_BDD_GROK_HOST (`local` by default, ~/.grok/config.yaml), and that user must be the page's: the guard
   below fails the first claim otherwise, since the rows a user sees are the rows granted to them. The reads and the
   sweeps shared by the step modules live here: a step module imported by another would be loaded twice. */
import {appendFileSync} from 'node:fs';
import {execFile} from 'node:child_process';
import {createRequire} from 'node:module';
import type {Page} from '@playwright/test';
import {expect, fixtureFamilies, isStaleFixture, pollMs, serverNow} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const GROK_HOST = process.env.FORGE_BDD_GROK_HOST ?? 'local';
/** A poll of the server: every call signs in, so the reads are spaced, not the page's 100 ms. */
export const SERVER_POLL = {intervals: [1000, 2000, 2000, 3000], timeout: pollMs(30000)};
const CALL_TIMEOUT_MS = 60000;

export const MODEL = 'forge.model';
export const RUN = 'forge.training_run';
export const MODEL_FOLDER = /^file:\/\/(System:DomainFiles\/forge\/model\/[0-9a-f-]{36})\/[^/]+$/i;
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

export const quoted = (s: string): string => `"${s.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}"`;

let cliLogin: Promise<string> | undefined;

/** Fails when the CLI's user is not the page's: the claims would read another user's rows. The page's login is read
 * every time (a feature can sign in as another account); the CLI's once. */
export async function guardIdentity(page: Page): Promise<void> {
  cliLogin ??= grokS(['raw', 'GET', 'users/current']).then((u) => String((u as {login?: unknown} | null)?.login ?? ''));
  const [cli, shown] = await Promise.all([cliLogin, page.evaluate(() => String(grok.shell.user.login))]);
  if (cli !== shown) {
    throw new Error(`the grok s alias "${GROK_HOST}" is "${cli}", the page is "${shown}": the server claims would read ` +
      'another user\'s rows (set FORGE_BDD_GROK_HOST to an alias of the page\'s user)');
  }
}

/* EMS returns `created_on` with its zone (`2026-10-06T16:12:19.015442Z` in `--output json`), so Date.parse reads the
   server's instant. */
const createdOn = (row: Record<string, unknown>): number => Date.parse(String(row.created_on ?? '')) || 0;

/** The `forge-test-…` rows of [table] (a name column [column]) that are [names] or a stale member of their family:
 * the ones a killed run of the feature left an hour or more ago, by the server's clock. */
export async function fixtureRows(page: Page, table: string, column: string, names: string[], columns: string[]):
  Promise<Record<string, unknown>[]> {
  const [listed, now] = await Promise.all([grokS(['domains', 'query', table, '--filter', `${column} like "forge-test-%"`,
    '--columns', [column, ...columns].join(','), '--limit', '1000']), serverNow(page)]);
  const families = fixtureFamilies(names);
  return rowsOf(listed).filter((r) => {
    const name = String(r[column] ?? '');
    return names.includes(name) || isStaleFixture({name, friendlyName: name, createdOn: createdOn(r)}, families, now);
  });
}

/** The count of [table] rows matching any of the [terms] (`column = "value"`), in filters of 20 terms. */
async function countWhere(table: string, terms: string[]): Promise<number> {
  const chunks: string[][] = [];
  for (let i = 0; i < terms.length; i += 20)
    chunks.push(terms.slice(i, i + 20));
  const counts = await Promise.all(chunks.map((chunk) => countOf(table, chunk.join(' or '))));
  return counts.reduce((a, b) => a + b, 0);
}

const termsOf = (column: string, values: string[]): string[] => values.map((v) => `${column} = ${quoted(v)}`);

/** The one model row named [name], or why there is not exactly one. */
export async function modelRow(name: string): Promise<Record<string, unknown> | string> {
  const rows = rowsOf(await grokS(['domains', 'query', MODEL, '--filter', `name = ${quoted(name)}`, '--limit', '2']));
  return rows.length === 1 ? rows[0] : `${rows.length} Forge models named "${name}"`;
}

export async function requireModelRow(name: string): Promise<Record<string, unknown>> {
  const row = await modelRow(name);
  if (typeof row === 'string')
    throw new Error(row);
  return row;
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

/** Deletes the feature's models (and the stale ones of their families) in the page, then reads the server back by the
 * exact names and ids: no such row, no model file folder of theirs, no training copy of theirs. */
export async function sweepModels(page: Page, names: string[]): Promise<void> {
  await guardIdentity(page);
  const rows = await fixtureRows(page, MODEL, 'name', names, ['blob', 'dataset_table_id']);
  const ids = rows.map((r) => String(r.id));
  await page.evaluate((list) => Promise.all(list.map((id) => grok.functions.call('Forge:deleteModel', {model: id}))), ids);
  await expect.poll(() => countWhere(MODEL, [...termsOf('name', names), ...termsOf('id', ids)]),
    {...SERVER_POLL, message: `Forge models named ${names.join(', ')} or deleted by id (${ids.length}) left on the server`})
    .toBe(0);
  await Promise.all(rows.map(async (row) => {
    const folder = MODEL_FOLDER.exec(String(row.blob ?? ''))?.[1];
    if (folder !== undefined) {
      await expect.poll(() => folderFiles(folder),
        {...SERVER_POLL, message: `the model file folder ${folder} of "${row.name}"`}).toEqual([]);
    }
    if (row.dataset_table_id)
      await expectTables(`${row.name}${COPY_SUFFIX}`, 0);
  }));
}

/** Deletes the runs of [table] (and of the stale tables of its family) and reads the count of every one back. */
export async function sweepRuns(page: Page, table: string): Promise<void> {
  await guardIdentity(page);
  const stale = (await fixtureRows(page, RUN, 'dataset_name', [table], [])).map((r) => String(r.dataset_name));
  const tables = [...new Set([table, ...stale])];
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
