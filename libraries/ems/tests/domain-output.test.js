import {emsBackends} from '../src/core/domain-backend.js';
import {mock} from 'node:test';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope} from '@datagrok-libraries/u2';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');

// the pickers' frames: the stub's DataFrame and Column stand in for the platform's typed ones
DG.Column.fromType ??= (type, name) => new DG.Column(name, type);
DG.DataFrame.fromColumns ??= (columns) => new DG.DataFrame(columns.map((c) => ({name: c.name, type: c.dart.type})));

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    const savedTable = grok.dapi.domains.table;
    mock.timers.enable({apis: ['setInterval']});
    try {
      await body();
    } finally {
      mock.timers.reset();
      grok.dapi.domains.table = savedTable;
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
    assert.equal(DomainApp.live.size, 0, 'no app left behind');
  });
}

/** The memory backend with every probe recorded by table, and `grok.dapi.domains.table(address).queryDf` — what
 * the aggregations read — answered from its rows: the rows whose `path` names one of the ids, grouped by it. */
function served() {
  const inner = backend();
  const probes = [];
  const be = {
    table: async (address) => {
      const t = await inner.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {
        probe: (spec) => {
          probes.push(address);
          return t.probe(spec);
        },
      });
    },
    saveAll: (edits) => inner.saveAll(edits),
  };
  let reads = 0;
  grok.dapi.domains.table = (address) => ({queryDf: async ({filter}) => {
    reads++;
    const rows = (await (await inner.table(address)).query({})).filter((r) => filter.value.includes(r[filter.property]));
    return {rowCount: rows.length, groupBy: ([path]) => {
      const measures = [];
      const by = new Map();
      for (const r of rows)
        by.set(r[path], [...by.get(r[path]) ?? [], r]);
      const keys = [...by.keys()];
      return {
        add(aggregation, column, name) {
          measures.push({aggregation, column, name});
          return this;
        },
        aggregate: () => ({rowCount: keys.length, get: (name, i) => {
          if (name === path)
            return keys[i];
          const m = measures.find((m) => m.name === name);
          const group = by.get(keys[i]);
          return m.aggregation === 'count' ? group.length : group.reduce((s, r) => s + (r[m.column] ?? 0), 0);
        }}),
      };
    }};
  }});
  return {be, probes, reads: () => reads};
}

async function settle() {
  for (let j = 0; j < 4; j++)
    await flush();
}

async function tick() {
  mock.timers.tick(1000);
  await settle();
}

/** The Projects app as a table with count(Issues.Title) and sum(Issues.Number) in the Output pane. */
async function projects(stand) {
  emsBackends.domain = stand.be;
  const table = await domains.table('grit.project');
  const app = testOwner().build(() => domains.app({table, base: '/apps/T/Projects', pageSize: 10, live: true,
    liveMs: 1000}));
  await flush();
  app.view.value = 'table';
  app.output.show();
  await settle();
  app.output.aggregations.value = [{column: 'Issues.Title', aggregation: 'count'},
    {column: 'Issues.Number', aggregation: 'sum'}];
  await settle();
  return app;
}

const shown = (app, id) => {
  const row = app.output.values.value.get(id) ?? {};
  return [row['count(Issues.Title)'], row['sum(Issues.Number)']];
};

// EQ-063: the aggregation columns live in the child tables, so a write there has to reach them
scoped('an aggregation column follows a child row saved in this session (EQ-063)', async () => {
  const app = await projects(served());
  assert.deepEqual(shown(app, 'p1'), [2, 3]);
  const issues = await domains.table('grit.issue');
  const child = testOwner().build(() => issues.source({query: 'project_id = "p1"', pageSize: 10}));
  await flush();
  child.newRow({project_id: 'p1', number: 7, title: 'Added from the project page'});
  assert.equal(await child.save(), true);
  await settle();
  assert.deepEqual(shown(app, 'p1'), [3, 10], 'a child row saved in this session re-fetches the aggregations');
  child.dispose();
  app.dispose();
});

scoped('an aggregation column follows a peer\'s child write on the live probe, one probe per child table (EQ-063)', async () => {
  const stand = served();
  const app = await projects(stand);
  await tick();
  const before = stand.probes.length;
  const reads = stand.reads();
  await tick();
  assert.equal(stand.reads(), reads, 'nothing moved, nothing is read again');
  const issue = await stand.be.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p2', number: 4, title: 'Theirs'}}]);
  await tick();
  assert.deepEqual(shown(app, 'p2'), [2, 5], 'a peer\'s child insert seen by the live probe re-fetches the aggregations');
  assert.deepEqual(stand.probes.slice(before).sort(), ['grit.issue', 'grit.issue', 'grit.project', 'grit.project'],
    'one probe per table per tick, not one per aggregation');
  app.dispose();
});
