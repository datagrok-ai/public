// the evaluation half of u2's filter corpus, run over the memory backend

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';

const corpus = JSON.parse(readFileSync(fileURLToPath(
  new URL(import.meta.resolve('@datagrok-libraries/u2/tests/filter-grammar.corpus.json'))), 'utf8'));
/* The corpus' second half: every entry carrying `matches` EVALUATED against the fixture, so the
   grammar's meaning is pinned and not only its shape. The SQL twin is
   `core/server/datlas/test/services/domain_filter_corpus_test.dart`, which builds the same table
   from the same block and runs the same entries against Postgres. */
const fixture = corpus.fixture;
/** The corpus' column type → the property type a reader must have built it as. */
const PROPERTY_TYPE = {string: 'string', int: 'int', float: 'double', bool: 'bool',
  datetime: 'datetime', ref: 'string'};
const be = new MemoryDomainBackend({name: 'corpus', tables: {[fixture.table]: {
  hierarchy: fixture.hierarchy === true, businessKey: fixture.businessKey, columns: fixture.columns}}});
const table = be.tableSync(`corpus.${fixture.table}`);
/** `code` → the id the row was given; an `under` value naming a code is that row's id. */
const ids = {};
for (const {parent, ...values} of fixture.rows) {
  const [row] = await table.transaction([{op: 'insert', table: fixture.table,
    values: {...values, ...(parent === undefined ? {} : {parent_id: ids[parent]})}}]);
  ids[values.code] = row.id;
}

test('corpus fixture: the table the evaluators build is exactly the corpus\' columns and types', () => {
  const declared = Object.entries(fixture.columns).map(([name, c]) => `${name}:${PROPERTY_TYPE[c.type]}`);
  const built = table.properties.filter((p) => fixture.columns[p.name] !== undefined)
    .map((p) => `${p.name}:${p.propertyType ?? p.type}`);
  assert.deepEqual(built, declared, 'a fixture change is a red test on both sides, never a silent drift');
  assert.equal(table.info.hierarchy, fixture.hierarchy === true);
  assert.deepEqual(table.info.businessKey, fixture.businessKey);
  assert.ok(corpus.entries.filter((e) => e.matches !== undefined).length >= 25);
});

for (const entry of corpus.entries.filter((e) => e.matches !== undefined)) {
  test(`corpus ${entry.id}: selects ${JSON.stringify(entry.matches)}`, async () => {
    const input = entry.input.replace(/(\bunder\s+)"([^"]*)"/g,
      (all, head, code) => ids[code] === undefined ? all : `${head}"${ids[code]}"`);
    const rows = await table.query({filter: input, limit: 1000});
    assert.deepEqual(rows.map((r) => r.code).sort(), entry.matches);
  });
}
