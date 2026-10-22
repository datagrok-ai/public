/* What a source sends the server, over the memory backend: a host operator's condition answered by the app as
   the ids it matches (under and/or/not, from a tree or from text), a search spelled like the caption column's
   formatted number, and a builder's incomplete rows kept out of the query. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Filters} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    const unregister = [];
    try {
      await body(unregister);
    } finally {
      for (const off of unregister)
        off();
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const CHEM = {
  name: 'chem',
  tables: {
    compound: {
      businessKey: ['code'], singularName: 'Compound', pluralName: 'Compounds',
      columns: {
        code: {type: 'int', format: 'CP-0000'},
        structure: {type: 'string', semType: 'Molecule', required: true},
        label: {type: 'string', searchable: true},
        weight: {type: 'float'},
      },
    },
  },
};
const COMPOUNDS = [
  {id: 'c1', code: 1, structure: 'CCO', label: 'ethanol', weight: 46.07},
  {id: 'c2', code: 2, structure: 'CC(=O)O', label: 'acetic acid'},
  {id: 'c12', code: 12, structure: 'c1ccccc1O', label: 'phenol', weight: 94.11},
];

/** `has` answers by substring — chemistry is the app's; the calls are counted. */
async function compounds() {
  emsBackends.domain = new MemoryDomainBackend(CHEM, {rows: {compound: COMPOUNDS.map((r) => ({...r}))}});
  const table = await domains.table('chem.compound');
  const calls = [];
  const off = table.operators.add({id: 'has', label: 'has', semType: 'Molecule',
    parameter: {name: 'threshold', label: 'At least', default: 0.5},
    rows: async (c) => {
      calls.push(c);
      return COMPOUNDS.filter((r) => r.structure.includes(String(c.value))).map((r) => r.id);
    }});
  return {table, calls, off};
}

const ids = (source) => source.rows.items.value.map((r) => r.id).sort();

scoped('a host condition is sent as the ids it matches, under and, or and not', async (unregister) => {
  const {table, calls, off} = await compounds();
  unregister.push(off);
  const source = table.source({query: Filters.group('and', [Filters.cond('structure', 'has', 'CC')])});
  await flush();
  assert.deepEqual(ids(source), ['c1', 'c2']);
  assert.equal(source.total.value, 2);
  assert.equal(calls.length, 1, 'answered once per load: the count and the page share it');
  assert.deepEqual(source.readScope().filter, Filters.toDomainTree(
    Filters.group('and', [Filters.cond('id', 'in', ['c1', 'c2'])])), 'the server never sees the host operator');
  source.query.value = Filters.group('or', [Filters.cond('structure', 'has', 'O)'), Filters.cond('label', '=', 'phenol')]);
  await flush();
  assert.deepEqual(ids(source), ['c12', 'c2']);
  source.query.value = {...Filters.group('and', [Filters.cond('structure', 'has', 'CC')]), not: true};
  await flush();
  assert.deepEqual(ids(source), ['c12']);
  source.query.value = Filters.group('and', [Filters.cond('structure', 'has', 'N')]);
  await flush();
  assert.deepEqual(ids(source), [], 'no hits is no rows');
  source.dispose();
});

scoped('text naming a host operator parses, threshold and all, and formats back the same', async (unregister) => {
  const {table, calls, off} = await compounds();
  unregister.push(off);
  const source = table.source({query: 'structure has(0.8) "CC" and weight > 40'});
  await flush();
  assert.deepEqual(ids(source), ['c1']);
  assert.equal(calls[0].options?.threshold, 0.8);
  const {root, problems} = Filters.parse('structure has(0.8) "CC"', source.schema, 'domain');
  assert.deepEqual(problems, []);
  assert.equal(Filters.format(root), 'structure has(0.8) "CC"');
  source.dispose();
});

scoped('a search spelled like the formatted key finds that row by its number', async (unregister) => {
  const {table, off} = await compounds();
  unregister.push(off);
  const source = table.source({search: 'CP-0012'});
  await flush();
  assert.deepEqual(ids(source), ['c12']);
  assert.equal(source.readScope().search, undefined, 'the server is not asked to search the text');
  source.search.value = 'cp-12';
  await flush();
  assert.deepEqual(ids(source), ['c12']);
  source.search.value = 'phen';
  await flush();
  assert.deepEqual(ids(source), ['c12'], 'other text is still a text search');
  source.search.value = '12';
  await flush();
  assert.deepEqual(ids(source), ['c12'], 'a bare number is the key too');
  assert.equal(source.readScope().search, undefined);
  source.search.value = '12a';
  await flush();
  assert.equal(source.readScope().search, '12a', 'text with digits in it is a text search');
  source.dispose();
});

scoped('the builder applies only the rows that are complete', async (unregister) => {
  const {table, off} = await compounds();
  unregister.push(off);
  const source = table.source();
  const filters = domains.filters(source, {mode: 'builder'});
  document.body.append(filters.root);
  await flush();
  const builder = filters.input.value;
  builder.value.value = Filters.group('and', [Filters.cond('label', '=')]);
  await flush();
  assert.equal(source.query.value, '', 'a condition without a value is not applied');
  assert.equal(builder.value.value.nodes.length, 1, 'and it stays in the builder');
  assert.equal(filters.problem.value, null, 'the rows answer what is applied: nothing is stale');
  builder.value.value = Filters.group('and', [Filters.cond('label', '=', 'phenol')]);
  await flush();
  assert.equal(source.query.value, 'label = "phenol"');
  await flush();
  assert.deepEqual(ids(source), ['c12']);
  filters.dispose();
  source.dispose();
});

scoped('a column of a table that answers no host operator keeps the core operators, whichever opened first', async (unregister) => {
  const LAB = {name: 'chem', tables: {
    compound: CHEM.tables.compound,
    addition: {singularName: 'Addition', pluralName: 'Additions',
      columns: {name: {type: 'string', isName: true}, structure: {type: 'string', semType: 'Molecule'}}},
  }};
  emsBackends.domain = new MemoryDomainBackend(LAB, {rows: {compound: COMPOUNDS.map((r) => ({...r})),
    addition: [{id: 'a1', name: 'HCl', structure: 'Cl'}]}});
  const compounds = await domains.table('chem.compound');
  unregister.push(compounds.operators.add({id: 'has', label: 'has', semType: 'Molecule', rows: async () => []}));
  const additions = await domains.table('chem.addition');
  const offer = async (table) => {
    const source = table.source();
    const filters = domains.filters(source, {mode: 'builder'});
    document.body.append(filters.root);
    await flush();
    const schema = filters.input.value.options.schema;
    const prop = Filters.property(schema, 'structure');
    const ids = schema.operators(prop, Filters.operators.for(prop)).map((o) => o.id);
    const valid = Filters.validate(Filters.group('and', [Filters.cond('structure', 'like', 'Cl')]), schema, 'domain');
    filters.dispose();
    source.dispose();
    return {ids, valid};
  };
  const compound = await offer(compounds);
  assert.deepEqual(compound.ids, ['has', 'is null', 'is not null'], 'the table that answers it offers that set alone');
  const addition = await offer(additions);
  assert.equal(addition.ids.includes('has'), false, 'another table is not offered what it cannot answer');
  for (const id of ['=', 'like', 'in', 'is null'])
    assert.equal(addition.ids.includes(id), true, `and keeps "${id}"`);
  assert.deepEqual(addition.valid, [], 'a core condition on its column validates');
});

scoped('a builder row whose parameter is wrong holds the filter applied before it, and says why', async (unregister) => {
  const {table, off} = await compounds();
  unregister.push(off);
  const source = table.source({query: 'structure has(0.8) "CC"'});
  const filters = domains.filters(source, {mode: 'builder'});
  document.body.append(filters.root);
  await flush();
  await flush();
  const builder = filters.input.value;
  const row = builder.value.value.nodes[0];
  const threshold = async (value) => {
    builder.value.value = {...builder.value.value, nodes: [{...row, options: {threshold: value}}]};
    await flush();
  };
  await threshold(null);
  assert.equal(source.query.value, 'structure has(0.8) "CC"', 'emptied to be retyped: the applied filter stands');
  assert.equal(filters.problem.value, 'At least: a number', 'and the list is told it does not answer the builder');
  await threshold(0.9);
  assert.equal(source.query.value, 'structure has(0.9) "CC"');
  assert.equal(filters.problem.value, null);
  builder.value.value = Filters.group('and', [...builder.value.value.nodes, Filters.cond('label', '=')]);
  await flush();
  assert.equal(source.query.value, 'structure has(0.9) "CC"', 'a row still being written is left out, as before');
  assert.equal(filters.problem.value, null);
  filters.dispose();
  source.dispose();
});
