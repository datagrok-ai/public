// domainFilterSchema — moved from u2's FilterSchemas.forDomainTable

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {Scope, Filters, KIND} from '@datagrok-libraries/u2';
import {TYPE} from 'datagrok-api/u2core';
import {Property} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {domainFilterSchema} = await import('../src/dg/filter-schema.js');

const abort = new AbortController().signal;

const propOf = (schema, name) => schema.properties.find((p) => p.name === name);

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      delete grok.dapi.domains;
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

// --- forDomainTable ---

function domainsDouble(filters = 'full') {
  const calls = {rowProperties: [], facets: []};
  const tables = {
    'plates.plate': [
      new Property('name', TYPE.STRING, {friendlyName: 'Name', nullable: false,
        get: (r) => r.name, set: (r, v) => r.name = v}),
      new Property('project_id', TYPE.STRING, {friendlyName: 'Project', semType: 'Core.projects'}),
      new Property('volume', TYPE.FLOAT),
    ],
    'Core.projects': [new Property('name', TYPE.STRING), new Property('created_on', TYPE.DATE_TIME)],
  };
  grok.dapi.domains = {
    registry: {
      rowProperties: async (address) => {
        calls.rowProperties.push(address);
        return tables[address];
      },
    },
    table: (address) => ({
      access: async () => ({support: {filters}}),
      facets: async (spec) => {
        calls.facets.push([address, spec]);
        const column = spec.facets[0].column;
        // an external ('basic') facet answers no unfiltered total
        return {facets: {v: {categories: column === 'project_id' ?
          [{value: 'p1', display: 'Alpha', total: 3, filtered: 1}] :
          [{value: 'Ann', total: filters === 'basic' ? null : 2, filtered: 2}]}}};
      },
    }),
  };
  return calls;
}

smoke('forDomainTable: registry properties are copied field by field, a dotted semType is the ref', async () => {
  const calls = domainsDouble();
  const schema = await domainFilterSchema('plates.plate');
  assert.deepEqual(calls.rowProperties, ['plates.plate']);

  const source = (await grok.dapi.domains.registry.rowProperties('plates.plate'))[0];
  assert.deepEqual(Object.keys({...source}), ['dart'], 'a spread of the platform property yields only its handle');
  const name = propOf(schema, 'name');
  assert.deepEqual([name.propertyType, name.friendlyName, name.nullable], [TYPE.STRING, 'Name', false]);
  assert.equal(typeof name.get, 'function');
  assert.equal(name.ref, undefined);

  const project = propOf(schema, 'project_id');
  assert.equal(project.ref, 'Core.projects');
  assert.equal(Filters.kindOf(project), KIND.REF);
  assert.equal(schema.renderer(project).caption({type: 'Core.projects', id: 'p1', name: 'Alpha'}), 'Alpha');
  assert.equal(schema.renderer(project).caption({type: 'Core.projects', id: 'p1'}), 'p1');
});

smoke('forDomainTable: values go through the categories facet, refs come back as FilterRef', async () => {
  const calls = domainsDouble();
  const schema = await domainFilterSchema('plates.plate');

  const refs = await schema.values(propOf(schema, 'project_id'), 'al', abort);
  assert.deepEqual(calls.facets, [['plates.plate', {facets: [
    {id: 'v', kind: 'categories', column: 'project_id', search: 'al', limit: 50}]}]]);
  assert.deepEqual(refs, [{value: {type: 'Core.projects', id: 'p1', name: 'Alpha'}, label: 'Alpha', count: 3}]);

  const names = await schema.values(propOf(schema, 'name'), '', abort);
  assert.deepEqual(names, [{value: 'Ann', label: 'Ann', count: 2}]);
});

smoke('forDomainTable: a basic backend is asked for values on strings and refs only; a null total carries no count',
  async () => {
    const calls = domainsDouble('basic');
    const schema = await domainFilterSchema('plates.plate');

    assert.deepEqual(await schema.values(propOf(schema, 'volume'), '1', abort), []);
    assert.deepEqual(calls.facets, [], 'a float column issues no facets call');

    const names = await schema.values(propOf(schema, 'name'), 'an', abort);
    assert.equal(calls.facets.length, 1);
    assert.deepEqual(names, [{value: 'Ann', label: 'Ann'}], 'no count key at all, never a null one');
    const refs = await schema.values(propOf(schema, 'project_id'), '', abort);
    assert.equal(calls.facets.length, 2, 'a ref column is still asked');
    assert.equal(refs[0].count, 3);
  });

smoke('forDomainTable: resolveRef loads the target table', async () => {
  const calls = domainsDouble();
  const schema = await domainFilterSchema('plates.plate');
  const project = propOf(schema, 'project_id');

  const target = await schema.resolveRef(project);
  assert.deepEqual(calls.rowProperties, ['plates.plate', 'Core.projects']);
  assert.deepEqual(target.properties.map((p) => p.name), ['name', 'created_on']);
  assert.equal(Filters.kindOf(propOf(target, 'created_on')), KIND.DATE_TIME);
});
