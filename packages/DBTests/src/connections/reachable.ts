import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {delay, tests} from '@datagrok-libraries/test/src/test';

// A provider that refuses the login or does not answer costs every test on it the whole test timeout:
// grok_connect retries the connection for 30 s before it reports the failure. So a connection is tested
// once per run, and every test on it starts from that answer: it fails at once, with the provider's reason.

/** How long a category waits for the answers about its connections before its first test. */
export const CATEGORY_WAIT = 60000;
/** How long a test waits for an answer its category did not get. */
const TEST_WAIT = 20000;
/** How many connections are tested at a time. */
const CONCURRENCY = 4;

const probes = new Map<string, Promise<string>>();
const waiting: (() => void)[] = [];
let running = 0;

/** The result of the connection's test, 'ok' or the reason; the connection is tested once per run. */
function probe(connectionId: string): Promise<string> {
  let result = probes.get(connectionId);
  if (!result) {
    result = (async () => {
      if (running >= CONCURRENCY)
        await new Promise<void>((resolve) => waiting.push(resolve));
      running++;
      try {
        return await (await grok.dapi.connections.find(connectionId)).test();
      } catch (e: any) {
        return String(e?.message ?? e);
      } finally {
        running--;
        waiting.shift()?.();
      }
    })();
    probes.set(connectionId, result);
  }
  return result;
}

/** 'ok' when the connection answers its test, otherwise the reason. Waits for the answer for at most
 * [wait] ms; the test of the connection itself goes on, and a later call gets its answer. */
export async function reachable(connection: DG.DataConnection, wait: number = CATEGORY_WAIT): Promise<string> {
  return await Promise.race([probe(connection.id),
    delay(wait).then(() => `the connection test has not answered in ${Math.round(wait / 1000)} s`)]);
}

/** Fails with the provider's reason when the connection does not answer its test. */
export async function requireConnection(connection: DG.DataConnection | string, wait?: number): Promise<void> {
  const c: DG.DataConnection = typeof connection === 'string' ? await grok.functions.eval(connection) : connection;
  const result = await reachable(c, wait);
  if (result !== 'ok')
    throw new Error(`${c.friendlyName || c.name} is not reachable: ${result}`);
}

/** Starts the tests of the connections without waiting for them. */
export function prewarm(connections: DG.DataConnection[]): void {
  for (const c of connections)
    probe(c.id);
}

/** Waits until every connection has answered its test, for at most [wait] ms. */
export async function answered(connections: DG.DataConnection[], wait: number = CATEGORY_WAIT): Promise<void> {
  await Promise.all(connections.map((c) => reachable(c, wait)));
}

const guarded = new WeakSet<object>();

/** The query tests of the package (`--test:` annotations) with the connection each one runs on. Every such
 * test first asks whether its connection is reachable, and its category gets the answers before it starts.
 * Returns the connections the tests need, by category. */
export async function guardQueryTests(packageId: string): Promise<{[category: string]: DG.DataConnection[]}> {
  const queries = await grok.dapi.queries.filter(`package.id = "${packageId}"`).include('connection').list();
  const connectionOf = new Map<string, DG.DataConnection>(queries.map((q) => [q.name, q.connection]));
  const used: {[category: string]: DG.DataConnection[]} = {};
  for (const [name, cat] of Object.entries(tests)) {
    for (const t of cat.tests ?? []) {
      // the library names a query's tests after the query, numbered when it has several, and marks them `node`
      const connection = t.options?.node ?
        connectionOf.get(t.name) ?? connectionOf.get(t.name.replace(/ \d+$/, '')) : null;
      if (!connection)
        continue;
      used[name] ??= [];
      if (!used[name].some((c) => c.id === connection.id))
        used[name].push(connection);
      if (guarded.has(t))
        continue;
      guarded.add(t);
      const run = t.test;
      t.test = async () => {
        await requireConnection(connection, TEST_WAIT);
        return await run();
      };
    }
    if (used[name] && !guarded.has(cat)) {
      guarded.add(cat);
      const before = cat.before;
      cat.before = async () => {
        await answered(used[name]);
        await before?.();
      };
    }
  }
  return used;
}
