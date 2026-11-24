import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {ImportLinks} = await import('../src/dg/import-links.js');

// a table keyed by a whole number the app shows as DG-000000
const KEY = {name: 'molregno', propertyType: 'int', type: 'int', format: 'DG-000000', friendlyName: 'Grok ID'};
const link = {
  prop: {name: 'compound_id', semType: 'moltrack.compound'},
  target: {address: 'moltrack.compound', info: {singularName: 'Compound', pluralName: 'Compounds'}, properties: [KEY]},
  fields: [{prop: KEY, path: 'molregno'}],
};
const STORED = [{id: 'c72', molregno: 72}, {id: 'c73', molregno: 73}];

test('an import links a key typed as its format to the record the bare number names, and to no other', async () => {
  grok.dapi.domains.table = () => ({query: async ({filter}) => STORED.filter((r) => filter[0].value.includes(r.molregno))});
  const found = await ImportLinks.resolve(link, link.fields, [['DG-000072'], ['dg-000073'], ['72'], ['nonesuch']]);
  assert.deepEqual(found.map((m) => m.id), ['c72', 'c73', 'c72', null]);
  assert.equal(found[3].error, 'No compound with Grok ID nonesuch');
  const {DomainTable} = await import('../src/dg/facade.js');
  assert.equal(DomainTable.parseKey({...KEY, format: 'C10-0000'}, 'C10-0042'), 42, 'a zero in the literal prefix is not a digit');
});
