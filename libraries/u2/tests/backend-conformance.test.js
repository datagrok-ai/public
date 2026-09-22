/* The backend conformance suite run headless over the memory backend (holistic review 4, A4): the
   scenarios in conformance/scenarios.mjs are what BOTH backends must answer, and U2Demo's
   `U2: domain conformance` runs the same module against the platform's over the stand. Server
   semantics asserted here belong in the scenarios, never in memory-domain.test.js — a memory-only
   assertion of a server rule is how the backends drifted four times in one day. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {scenarios} from './conformance/scenarios.mjs';

/** The fixture table, column for column as the stand harness declares it in its throwaway schema. */
const SCHEMA = {name: 'conformance', tables: {folder: {
  hierarchy: true, singularName: 'Folder', businessKey: ['code'],
  columns: {
    code: {type: 'string', required: true, searchable: true},
    name: {type: 'string', isName: true, searchable: true},
    parent_id: {type: 'ref', ref: 'folder'},
  },
}}};

/** What this fixture answers, derived the way the stand harness derives it — the two flags the
 * manifest declares, then the table's own `support`, then the optional seam member it carries. A
 * constant list here would silently skip a scenario the stand runs (and did, for `updateWhere`). */
const declares = (table, required) => required === 'hierarchy' || required === 'softDelete' ? true :
  table.support[required] === true || table[required] !== undefined;

const t = {
  ok: (value, message) => assert.ok(value, message),
  equal: (actual, expected, message) => assert.equal(actual, expected, message),
  deepEqual: (actual, expected, message) => assert.deepEqual(actual, expected, message),
  rejects: (fn, code) => assert.rejects(fn, (e) => code === undefined || e.code === code),
};

/** The rows inserted in order, `'$ref:<key>'` replaced by the id the row it names was given — the
 * resolution both harnesses do, so a scenario never writes an id of its own. */
async function seed(table, rows) {
  const ids = {};
  for (const {key, ...values} of rows) {
    const resolved = Object.fromEntries(Object.entries(values).map(([column, value]) =>
      [column, typeof value === 'string' && value.startsWith('$ref:') ? ids[value.slice(5)] : value]));
    const [result] = await table.transaction([{op: 'insert', table: table.address, values: resolved}]);
    ids[key] = result.id;
  }
  return ids;
}

for (const scenario of scenarios) {
  test(`conformance: ${scenario.name}`, async () => {
    const table = new MemoryDomainBackend(SCHEMA).tableSync('conformance.folder');
    assert.deepEqual(scenario.requires.filter((name) => !declares(table, name)), [],
      'the memory fixture declares everything the scenarios ask for — it never skips');
    await scenario.run(table, await seed(table, scenario.seed), t);
  });
}
