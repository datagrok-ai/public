/* Fills the sources' platform seam (Q5) at dg-index import: functions through `DG.Func`, server
   collections and entities through `grok.dapi.*`, the open tables through `grok.shell`. The
   sources themselves stay platform-free and headless-testable — this is the only file that knows
   what a Datagrok function or a dapi collection actually is. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {backends, FuncDescriptorLike} from '../sources/backends.js';
import type {DataFrameLike} from '../sources/df-bindings.js';
import type {DapiPagerSourceLike} from './dapi-source.js';

type Collection = () => DapiPagerSourceLike<unknown>;

/** Curated on purpose: everything here is a paged `HttpDataSource`, which `dapiPager` drives.
 * `grok.dapi.files` is deliberately absent — it is a file browser, not a listing source. */
const DAPI: Record<string, Collection> = {
  users: () => grok.dapi.users,
  groups: () => grok.dapi.groups,
  projects: () => grok.dapi.projects,
  queries: () => grok.dapi.queries,
  connections: () => grok.dapi.connections,
  scripts: () => grok.dapi.scripts,
  spaces: () => grok.dapi.spaces,
} as unknown as Record<string, Collection>;

/** `Func.byName` matches the bare name; a package function is addressed `Package:Name`. */
function findFunc(name: string): DG.Func | null {
  const colon = name.indexOf(':');
  const found = colon < 0 ? DG.Func.find({name}) :
    DG.Func.find({package: name.slice(0, colon), name: name.slice(colon + 1)});
  return found.length > 0 ? found[0] : null;
}

function kindOf(func: DG.Func): FuncDescriptorLike['kind'] {
  return func instanceof DG.TableQuery ? 'table-query' : func instanceof DG.DataQuery ? 'query' :
    func instanceof DG.Script ? 'script' : 'function';
}

/** What one query entity is capped by right now: the limit it had before the first design-time run
 * touched it, and how many runs are still holding the cap. Keyed by the DART handle — `Func.find`
 * answers a fresh JS wrapper every lookup while `limit` lives on the entity behind it. */
const capping = new Map<unknown, {previous: number | undefined, runs: number}>();

/** OR-1: the only row cap the platform offers at this level is `TableQuery.limit`, which lives on
 * the shared query entity and is read while the call EXECUTES — so it must stay set across the
 * await, not only across `prepare`. That makes it shared state: the platform cannot cancel a run,
 * so two overlapping design-time runs interleave, and a naive save/restore would leave the query
 * capped for the rest of the session — the user's own execution of it included. Only the first
 * run captures the previous limit and only the last one puts it back. A plain-SQL query or an
 * arbitrary function runs uncapped: there is no FuncCall-level option to honor yet. */
async function capped<T>(func: DG.Func, maxRows: number | undefined,
  body: () => Promise<T>): Promise<T> {
  if (maxRows === undefined || !(func instanceof DG.TableQuery))
    return body();
  const query = func as DG.TableQuery;
  const key = query.dart;
  let held = capping.get(key);
  if (held === undefined) {
    held = {previous: query.limit, runs: 0};
    capping.set(key, held);
  }
  held.runs++;
  query.limit = maxRows;
  try {
    return await body();
  } finally {
    if (--held.runs === 0) {
      capping.delete(key);
      query.limit = held.previous;
    }
  }
}

backends.funcRunner = {
  find: (name) => {
    const func = findFunc(name);
    return func === null ? null :
      {name, kind: kindOf(func), inputs: func.inputs, outputs: func.outputs};
  },
  run: async (name, params, _abort, options) => {
    const func = findFunc(name);
    if (func === null)
      throw new Error(`"${name}" is not a known function`);
    const call = await capped(func, options?.maxRows, () => func.prepare(params).call());
    const outputs: Record<string, unknown> = {};
    for (const output of func.outputs)
      outputs[output.name] = call.outputs[output.name];
    return outputs;
  },
};

backends.dapi = (collection) => {
  const source = DAPI[collection];
  return source === undefined ? null : source();
};

backends.dapiFind = async (collection, id) => {
  const source = DAPI[collection];
  return source === undefined ? null :
    (source() as unknown as {find(id: string): Promise<unknown>}).find(id);
};

backends.workspace = {
  table: (name) => (grok.shell.tableByName(name) as unknown as DataFrameLike) ?? null,
  tableNames: () => grok.shell.tableNames,
  onTablesChanged: {
    subscribe: (next) => {
      const subs = [grok.events.onTableAdded, grok.events.onTableRemoved]
        .map((event) => event.subscribe(next));
      return {unsubscribe: () => {
        for (const sub of subs)
          sub.unsubscribe();
      }};
    },
  },
};

backends.tableFromRows = (rows) =>
  (DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.create(0)) as unknown as DataFrameLike;
