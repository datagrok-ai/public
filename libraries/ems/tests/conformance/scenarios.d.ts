import type {DomainTableLike} from '../../src/sources/domain-backend.js';

/** The assertions a scenario has — `node:assert/strict` in the headless harness, `expect` on the
 * stand, so a scenario never names a test framework. */
export interface ConformanceAssert {
  ok(value: unknown, message?: string): void;
  equal(actual: unknown, expected: unknown, message?: string): void;
  deepEqual(actual: unknown, expected: unknown, message?: string): void;
  /** The call refused — with that `code` on the refusal where one is named. */
  rejects(fn: () => Promise<unknown>, code?: string): Promise<void>;
}

/** One row of a scenario's fixture: `key` names it for `$ref:` and for the ids the run receives,
 * the rest are column values; a value of `'$ref:<key>'` is the id an earlier row was given. */
export interface ConformanceSeedRow {
  key: string;
  [column: string]: unknown;
}

/** One case both backends must answer the same way. `table` is the fixture table as the seam sees
 * it (`DomainTableLike`), `ids` maps each seed row's `key` to the id it was assigned. */
export interface ConformanceScenario {
  name: string;
  /** What the fixture table must declare for the scenario to mean anything — `'hierarchy'`,
   * `'softDelete'`, `'updateWhere'`; a harness whose fixture declares less skips it. */
  requires: string[];
  seed: ConformanceSeedRow[];
  run(table: DomainTableLike, ids: Record<string, string>, t: ConformanceAssert): Promise<void>;
}

export declare const scenarios: ConformanceScenario[];
