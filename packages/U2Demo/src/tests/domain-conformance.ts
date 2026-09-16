/* The backend conformance scenarios (`u2/tests/conformance/scenarios.mjs`) run against the
   PLATFORM backend — the same file the headless suite runs over `MemoryTable`, so a memory
   backend that drifts from the server fails on one side or the other. Fixture: a throwaway
   user-managed schema carrying ONE copy of the contract's `folder` table per scenario
   (hierarchy, soft delete, `code` as the business key), all created by a single `apply`; the
   category skips cleanly without the CreateDomainSchema privilege. */
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

const REF = '$ref:';

/** The fixture table a scenario gets to itself. A scenario must start from an EMPTY table, its
 * trash included — the headless harness builds a fresh backend per scenario, while a soft delete
 * here would leave rows the next scenario's `deleted: 'only'` read still sees. One table per
 * scenario is the cheapest way to say that: a single `apply` creates them all. */
const TABLE = (index: number): string => `folder${index}`;

category('U2: domain conformance', () => {
  const name = `zzu${`${Date.now()}`.slice(-8)}`;
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
    const tables: {[table: string]: object} = {};
    for (let index = 0; index < scenarios.length; index++)
      tables[TABLE(index)] = {
        hierarchy: true,
        softDelete: true,
        businessKey: ['code'],
        columns: {
          code: {type: 'string', required: true, unique: true},
          name: {type: 'string', isName: true},
          parent_id: {type: 'ref', ref: TABLE(index)},
        },
      };
    await grok.dapi.domains.schema(name).apply({tables});
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

  /** What the fixture answers: the two flags its manifest declares, then the table's own
   * `support`, then the optional seam member it carries — the same derivation the headless
   * harness makes, so a scenario cannot run on one side and be skipped on the other. */
  const declares = (table: DomainTableLike, required: string): boolean =>
    required === 'hierarchy' || required === 'softDelete' ? true :
      (table.support as unknown as Record<string, unknown>)[required] === true ||
      (table as any)[required] !== undefined;

  /** The scenario's seed as ONE transaction: `$ref:<key>` is the row that key named, which the
   * platform resolves through the op's own `ref`. */
  async function seed(table: DomainTableLike, rows: ConformanceSeedRow[]): Promise<{[key: string]: string}> {
    if (rows.length === 0)
      return {};
    const ops: DomainTransactionOpLike[] = rows.map(({key, ...columns}) => {
      const values: {[column: string]: any} = {};
      for (const [column, value] of Object.entries(columns))
        values[column] = typeof value === 'string' && value.startsWith(REF) ? `$${value.slice(REF.length)}` : value;
      return {op: 'insert', table: table.address, ref: key, values};
    });
    const results = await table.transaction(ops);
    const ids: {[key: string]: string} = {};
    rows.forEach((row, i) => ids[row.key] = results[i].id!);
    return ids;
  }

  /** `assert.equal` as the scenarios call it. NOT `expect`: its `expected` parameter DEFAULTS to
   * `true`, so `t.equal(value, undefined)` would compare the value against `true` instead. */
  const equal = (actual: unknown, expected: unknown, message?: string): void => {
    if (actual !== expected)
      throw new Error(`${message ? `${message}, ` : ''}Expected "${expected}", got "${actual}"`);
  };

  /** `node:assert/strict` as the scenarios call it, over the package test library. */
  const asserts: ConformanceAssert = {
    ok: (value, message) => expect(!!value, true, message),
    equal,
    deepEqual: (actual, expected, message) => {
      if (Array.isArray(expected))
        expectArray(actual as any[], expected);
      else if (expected !== null && typeof expected === 'object')
        expectObject(actual as any, expected as any);
      else
        equal(actual, expected, message);
    },
    rejects: (fn, code) => expectExceptionAsync(
      async () => { await fn(); }, code === undefined ? undefined : (e: any) => e?.code === code),
  };

  for (let index = 0; index < scenarios.length; index++) {
    const scenario = scenarios[index];
    const address = `${name}.${TABLE(index)}`;
    test(scenario.name, async () => {
      if (skipped())
        return;
      const table = await backends.domain!.table(address);
      const missing = scenario.requires.filter((required: string) => !declares(table, required));
      if (missing.length > 0) {
        console.log(`skipped: the fixture answers no ${missing.join(', ')}`);
        return;
      }
      await scenario.run(table, await seed(table, scenario.seed), asserts);
    });
  }
}, {owner: 'askalkin@datagrok.ai'});
