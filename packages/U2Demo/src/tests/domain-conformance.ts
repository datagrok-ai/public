/* The backend conformance scenarios (`u2/tests/conformance/scenarios.mjs`) run against the
   PLATFORM backend — the same file the headless suite runs over `MemoryTable`, so a memory
   backend that drifts from the server fails on one side or the other. Fixture: a throwaway
   user-managed schema with the contract's `folder` table (hierarchy, soft delete, `code` as the
   business key); the category skips cleanly without the CreateDomainSchema privilege. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectArray, expectExceptionAsync, expectObject, test}
  from '@datagrok-libraries/test/src/test';
import {backends} from '@datagrok-libraries/u2';
import type {DomainTableLike, DomainTransactionOpLike} from '@datagrok-libraries/u2';
import '@datagrok-libraries/u2/src/dg/index.js';
// extensionless: tsc reads `scenarios.d.ts` next to it, webpack resolves the `.mjs` itself
import {scenarios} from '@datagrok-libraries/u2/tests/conformance/scenarios';
import type {ConformanceAssert, ConformanceSeedRow} from '@datagrok-libraries/u2/tests/conformance/scenarios';

const TABLE = 'folder';
const REF = '$ref:';

category('U2: domain conformance', () => {
  const name = `zzu${`${Date.now()}`.slice(-8)}`;
  const address = `${name}.${TABLE}`;
  let table: DomainTableLike;
  let skip: string | null = null;

  before(async () => {
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Conformance probe'});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403)) {
        skip = 'no CreateDomainSchema privilege';
        return;
      }
      throw e;
    }
    await grok.dapi.domains.schema(name).apply({tables: {
      folder: {
        hierarchy: true,
        softDelete: true,
        businessKey: ['code'],
        columns: {
          code: {type: 'string', required: true, unique: true},
          name: {type: 'string', isName: true},
          parent_id: {type: 'ref', ref: 'folder'},
        },
      },
    }});
    table = await backends.domain!.table(address);
  });

  after(async () => {
    if (skip == null)
      await grok.dapi.domains.schema(name).delete();
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  /** What the fixture answers: the two flags its manifest declares, and for anything else the
   * optional seam member the table carries. A scenario asking for more is skipped, not failed. */
  const declares = (required: string): boolean =>
    required === 'hierarchy' || required === 'softDelete' ? true : (table as any)[required] !== undefined;

  /** The rows before this scenario, newest first: a child is created after its parent, so this
   * order never deletes a parent a live row still refers to. The business key's unique index
   * skips deleted rows, so the next scenario seeds the same codes. */
  async function purge(): Promise<void> {
    const rows = await table.query({sort: '!created_on', limit: 1000});
    const ops: DomainTransactionOpLike[] = rows.map((row) => ({op: 'delete', table: TABLE, id: String(row.id)}));
    if (ops.length > 0)
      await table.transaction(ops);
  }

  /** The scenario's seed as ONE transaction: `$ref:<key>` is the row that key named, which the
   * platform resolves through the op's own `ref`. */
  async function seed(rows: ConformanceSeedRow[]): Promise<{[key: string]: string}> {
    if (rows.length === 0)
      return {};
    const ops: DomainTransactionOpLike[] = rows.map(({key, ...columns}) => {
      const values: {[column: string]: any} = {};
      for (const [column, value] of Object.entries(columns))
        values[column] = typeof value === 'string' && value.startsWith(REF) ? `$${value.slice(REF.length)}` : value;
      return {op: 'insert', table: TABLE, ref: key, values};
    });
    const results = await table.transaction(ops);
    const ids: {[key: string]: string} = {};
    rows.forEach((row, i) => ids[row.key] = results[i].id!);
    return ids;
  }

  /** `node:assert/strict` as the scenarios call it, over the package test library. */
  const asserts: ConformanceAssert = {
    ok: (value, message) => expect(!!value, true, message),
    equal: (actual, expected, message) => expect(actual, expected, message),
    deepEqual: (actual, expected, message) => {
      if (Array.isArray(expected))
        expectArray(actual as any[], expected);
      else if (expected !== null && typeof expected === 'object')
        expectObject(actual as any, expected as any);
      else
        expect(actual, expected, message);
    },
    rejects: (fn, code) => expectExceptionAsync(
      async () => { await fn(); }, code === undefined ? undefined : (e: any) => e?.code === code),
  };

  for (const scenario of scenarios) {
    test(scenario.name, async () => {
      if (skipped())
        return;
      const missing = scenario.requires.filter((required: string) => !declares(required));
      if (missing.length > 0) {
        console.log(`skipped: the fixture answers no ${missing.join(', ')}`);
        return;
      }
      await purge();
      await scenario.run(table, await seed(scenario.seed), asserts);
    });
  }
}, {owner: 'askalkin@datagrok.ai'});
