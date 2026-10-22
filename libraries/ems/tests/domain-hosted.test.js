/* What a source sends the server, over the memory backend: a host operator's condition answered by the app as
   the ids it matches (under and/or/not and nested groups, from a tree, the builder or text) and reused until the
   table moves, a search spelled like the caption column's formatted number, and a builder's incomplete rows kept
   out of the query. */

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
  const source = table.source();
  const searched = async (text) => {
    source.search.value = text;
    await flush();
    return [ids(source), source.readScope().search];
  };
  assert.deepEqual(await searched('CP-0012'), [['c12'], undefined], 'the server is not asked to search the text');
  assert.deepEqual(await searched('cp-12'), [['c12'], undefined]);
  assert.deepEqual(await searched('12'), [['c12'], undefined], 'a bare number is the key too');
  assert.deepEqual(await searched('phen'), [['c12'], 'phen'], 'other text is still a text search');
  assert.deepEqual((await searched('12a'))[1], '12a', 'and text with digits in it');
  source.dispose();
});

scoped('nested groups through the builder: an incomplete row left out, each group\'s op and not kept, as text too',
  async (unregister) => {
    const {table, off} = await compounds();
    unregister.push(off);
    const source = table.source();
    const filters = domains.filters(source, {mode: 'builder'});
    document.body.append(filters.root);
    await flush();
    const builder = filters.input.value;
    const {cond, group} = Filters;
    const eth = () => group('and', [cond('weight', '>', 40), cond('label', 'like', 'eth')]);
    const apply = async (tree) => {
      builder.value.value = tree;
      await flush();
      await flush();
      const again = table.source({query: source.query.value});
      await flush();
      assert.deepEqual(ids(again), ids(source), `"${source.query.value}" read back as text answers the same rows`);
      again.dispose();
      return ids(source);
    };
    assert.deepEqual(await apply(group('or', [cond('structure', 'has', 'c1'), eth(), cond('label', '=')])),
      ['c1', 'c12']);
    assert.equal(source.query.value, 'structure has "c1" or (weight > 40 and label like "eth")',
      'the row without a value is not applied');
    assert.equal(builder.value.value.nodes.length, 3, 'and it stays in the builder');
    assert.equal(filters.problem.value, null, 'the rows answer what is applied: nothing is stale');
    const sent = JSON.stringify(source.readScope().filter);
    assert.match(sent, /"c12"/);
    assert.doesNotMatch(sent, /has/, 'the server never sees the host operator');
    assert.deepEqual(await apply({...builder.value.value, op: 'and'}), []);
    assert.deepEqual(await apply(group('and', [cond('structure', 'has', 'C'), {...eth(), not: true}])), ['c2'],
      'a null weight passes the not as SQL reads it: NULL and false is false');
    filters.dispose();
    source.dispose();
  });

scoped('a column of a table that answers no host operator keeps the core operators', async (unregister) => {
  const LAB = {name: 'chem', tables: {compound: CHEM.tables.compound,
    addition: {columns: {name: {type: 'string', isName: true}, structure: {type: 'string', semType: 'Molecule'}}}}};
  emsBackends.domain = new MemoryDomainBackend(LAB, {rows: {compound: [], addition: []}});
  const [compounds, additions] = await Promise.all([domains.table('chem.compound'), domains.table('chem.addition')]);
  unregister.push(compounds.operators.add({id: 'has', label: 'has', semType: 'Molecule', rows: async () => []}));
  const offered = async (table) => {
    const source = table.source();
    const filters = domains.filters(source, {mode: 'builder'});
    await flush();
    const schema = filters.input.value.options.schema;
    const prop = Filters.property(schema, 'structure');
    filters.dispose();
    source.dispose();
    return schema.operators(prop, Filters.operators.for(prop)).map((o) => o.id);
  };
  assert.deepEqual(await offered(compounds), ['has', 'is null', 'is not null'],
    'the table that answers it offers that set alone');
  const other = await offered(additions);
  assert.equal(other.includes('has'), false, 'another table is not offered what it cannot answer');
  assert.equal(other.includes('like'), true);
});

scoped('a host condition is answered again only once the table moves: a new order or search reuses its ids',
  async (unregister) => {
    const {table, calls, off} = await compounds();
    unregister.push(off);
    const source = table.source({query: 'structure has "CC"'});
    await flush();
    assert.deepEqual([ids(source), source.total.value, calls.length], [['c1', 'c2'], 2, 1],
      'answered once per load: the count and the page share it');
    source.sort.value = '!code';
    await flush();
    source.search.value = 'acid';
    await flush();
    assert.equal(calls.length, 1, 'the order and the search changed, the table did not');
    assert.deepEqual(ids(source), ['c2']);
    const other = table.source();
    await flush();
    other.newRow({code: 4, structure: 'CCC(=O)O', label: 'propionic acid'});
    assert.equal(await other.save(), true);
    await flush();
    assert.equal(calls.length, 2, 'a write moves the table: the condition is answered again');
    source.search.value = '';
    source.query.value = 'structure has "N"';
    await flush();
    assert.deepEqual(ids(source), [], 'no hits is no rows');
    other.dispose();
    source.dispose();
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
